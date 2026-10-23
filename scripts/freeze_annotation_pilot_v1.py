"""Freeze private pilot row selectors and measured features before modeling.

No expression matrix is copied. Clinical outcomes are never read. Unknown
material/case qualifications remain explicit rather than being guessed.
"""
import argparse
import collections
import datetime
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np
import pandas as pd

from annotation_pilot_sampling_v1 import balanced_cell_selection, stable_seed
from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_panel_aggregate_v1 import STATE_PANELS, IDENTITY_CANDIDATE_PANELS


MAJOR_POOL = {
    'T cell': 'T_NK', 'T/NK lymphoid': 'T_NK', 'T': 'T_NK', 'NK': 'T_NK',
    'T_NK_boundary': 'T_NK', 'Myeloid': 'Myeloid', 'Monocyte': 'Myeloid',
    'Macrophage': 'Myeloid', 'Epithelial cell': 'Epithelial', 'Epithelial': 'Epithelial',
}
ORTHOGONAL_MARKER_READOUTS = {
    'T_NK_related': ['CD2', 'CD7', 'TRBC1'],
    'B_related': ['CD19', 'CD22', 'CD79B'],
    'Plasma_related': ['FKBP11', 'DERL3', 'SDC1'],
    'Myeloid_related': ['CSF1R', 'LILRB1', 'FCGR2A'],
    'Epithelial_related': ['KRT7', 'KRT18', 'KRT19'],
    'Fibroblast_related': ['COL1A2', 'COL3A1', 'LUM'],
    'Endothelial_related': ['ESAM', 'EMCN', 'KDR'],
    'Mural_related': ['MCAM', 'PDGFRB'],
}


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def truth(values):
    return np.asarray([str(v).lower() in {'true', '1'} for v in values])


def valid_strings(values):
    return np.asarray([str(v) not in {'', 'None', 'nan', '__SOURCE_VALUE_MISSING__'} for v in values])


def eligibility_and_units(reader, entry, legacy, breast):
    dataset = entry['dataset']
    count = len(reader.rows(dataset))
    available = reader.field_names(dataset)
    if reader.is_legacy(dataset):
        sample = reader.frozen_expression_units(dataset).astype(str)
        metadata = legacy.set_index('sampleID')
        if not set(sample) <= set(metadata.index):
            raise ValueError('Native/frozen sampling units cannot be joined')
        per_cell = metadata.loc[sample]
        if entry['material_rule'] == 'official_breast_tumor_Total':
            mapped = breast.set_index('sampleID').loc[sample]
            eligible = (mapped.official_material_role.eq('tumor_mammary_unresolved_repeat_specimens') &
                        mapped.cell_population.eq('Total')).to_numpy()
            units = ('author_token_not_patient_claim::' + mapped.author_patient_token).to_numpy()
            basis = 'exact_author_token_blocks_global_patient_identity_not_promoted'
        elif entry['material_rule'] == 'frozen_primary_or_unspecified':
            eligible = per_cell.sample_class.eq('primary_or_unspecified_tumor').to_numpy()
            units = np.asarray(['sample_fallback::' + v for v in sample], dtype=object)
            basis = 'physical_expression_unit_fallback_not_patients'
        elif entry['material_rule'] == 'diagnostic_material_unknown_stress':
            eligible = np.ones(count, dtype=bool)
            units = np.asarray(['sample_fallback::' + v for v in sample], dtype=object)
            basis = 'sample_fallback_material_not_qualified_for_primary_tumor_claim'
        else:
            raise ValueError('Unsupported frozen legacy material rule')
    else:
        sample = reader.frozen_expression_units(dataset).astype(str)
        eligible = np.ones(count, dtype=bool)
        if 'tissue_equals' in entry:
            eligible &= reader.field(dataset, 'tissue').astype(str) == entry['tissue_equals']
        for name in ['ex_vivo_stimulation', 'annotation_exclude_from_primary', 'residual_doublet_flag']:
            if name in available:
                eligible &= ~truth(reader.field(dataset, name))
        units = np.asarray(['library_fallback::' + v for v in sample], dtype=object)
        verified = truth(reader.field(dataset, 'patient_identity_verified'))
        if verified.any():
            block = reader.field(dataset, 'patient_block').astype(str)
            if not valid_strings(block[verified]).all():
                raise ValueError('Verified cells lack their declared patient block')
            units[verified] = np.asarray(['reviewed_block::' + v for v in block[verified]], dtype=object)
        basis = 'verified_patient_blocks_where_provided_else_original_library_fallback'
    if not valid_strings(sample[eligible]).all():
        raise ValueError('Eligible cell lacks original expression unit')
    return eligible, units, sample, basis


def model_gene_plan(reader, entries):
    lookups = {entry['dataset']: reader.feature_lookup(entry['dataset']) for entry in entries}
    common = None
    for lookup in lookups.values():
        names = {name for name, cols in lookup.items()
                 if len(cols) == 1 and re.fullmatch(r'[A-Za-z][A-Za-z0-9_.-]*', name)
                 and not name.startswith(('ENSG', 'AMBIGUOUS'))}
        common = names if common is None else common & names
    collisions = set()
    for lookup in lookups.values():
        by_physical = collections.defaultdict(list)
        for name in common:
            by_physical[next(iter(lookup[name]))].append(name)
        collisions.update(name for names in by_physical.values() if len(names) > 1 for name in names)
    holdout = {gene for panels in [STATE_PANELS, IDENTITY_CANDIDATE_PANELS, ORTHOGONAL_MARKER_READOUTS]
               for genes in panels.values() for gene in genes}
    def control(name):
        return (name.startswith(('MT-', 'RPL', 'RPS', 'ERCC')) or
                re.match(r'^(IGH[VDJ]|IGK[VJ]|IGL[VJ]|TRA[VJ]|TRB[VDJ]|TRG[VJ]|TRD[VDJ])', name))
    candidates = sorted(name for name in common - collisions - holdout if not control(name))
    if len(candidates) < 5000:
        raise ValueError('Insufficient unambiguous measured common model features')
    return {'common_model_feature_candidates': candidates,
            'n_common_unique_names': len(common), 'n_collision_alias_names_excluded': len(collisions),
            'heldout_RNA_readouts': {'state': STATE_PANELS, 'identity_candidate': IDENTITY_CANDIDATE_PANELS,
                                     'marker': ORTHOGONAL_MARKER_READOUTS},
            'heldout_readouts_are_not_independent_ground_truth': True,
            'old_labels_used_for_inclusion_are_not_test_truth': True,
            'source_measurement_coverage': {
                source: {gene: len(lookup.get(gene, set())) == 1 for gene in sorted(holdout)}
                for source, lookup in lookups.items()},
            'HVG_and_scaling': 'not fitted here; training-only in source-holdout mapping; all-source discovery separately declared',
            'filters_do_not_modify_archive_gene_universe': True}


def main(args):
    start = time.monotonic()
    output = Path(args.output)
    spec_path = Path(args.spec)
    spec = json.loads(spec_path.read_text())
    source = Path(args.source_run)
    manifest = source / 'config/full_context_collection_v1.json'
    audit = json.loads((source / 'provenance/independent_stage1_audit_v1.json').read_text())
    if audit['status'] != 'INDEPENDENT_STAGE1_STRUCTURE_AND_FROZEN_LABEL_COUNT_AUDIT_PASS':
        raise ValueError('Native collection not independently accepted')
    for suite, expected, receipt in [('test_full_context_atlas_v1.py', 6, 'namespace_reader_tests_v2.json'),
                                      ('test_annotation_pilot_sampling_v1.py', 4, 'sampling_tests_v1.json')]:
        test = subprocess.run([sys.executable, '-B', str(output / 'tests' / suite), '-v'],
                              capture_output=True, text=True, timeout=120)
        save(output / 'provenance' / receipt,
             {'actual_exit_code': test.returncode, 'stdout': test.stdout, 'stderr': test.stderr})
        if test.returncode or ('Ran %d tests' % expected) not in test.stderr or 'skipped' in test.stderr.lower():
            raise ValueError('Native namespace/sampling counterexamples not accepted')
    for role in ['official_breast_map', 'official_breast_map_audit']:
        entry = spec['bindings'][role]
        if sha256_file(entry['path']) != entry['sha256']:
            raise ValueError('Changed accepted material bridge')
    legacy = pd.read_csv(args.legacy_metadata, sep='\t', dtype=str, keep_default_na=False)
    breast = pd.read_csv(spec['bindings']['official_breast_map']['path'], sep='\t', dtype=str, keep_default_na=False)
    reader = FullContextAtlas(manifest)
    save(output / 'config/frozen_pre_model_spec_v1.json', spec)
    (output / 'selectors').mkdir()
    summaries = []
    for entry in spec['sources']:
        dataset = entry['dataset']
        eligible, units, samples, basis = eligibility_and_units(reader, entry, legacy, breast)
        labels = reader.frozen_labels(dataset)
        pools = np.asarray([MAJOR_POOL.get(str(v), 'Other_context_or_unresolved') for v in labels['major']], dtype=object)
        chosen = []
        pool_counts = {}
        for pool, cap in spec['per_source_caps'].items():
            rows = balanced_cell_selection(units, labels['fine'], eligible & (pools == pool),
                                           cap, stable_seed(spec['sampling_seed'], dataset, pool),
                                           rare_reserve=spec['rare_fine_reserve_per_unit'])
            chosen.extend(rows.tolist())
            pool_counts[pool] = {'eligible': int((eligible & (pools == pool)).sum()), 'selected': len(rows)}
        positions = np.sort(np.asarray(chosen, dtype=np.int64))
        if len(np.unique(positions)) != len(positions) or not eligible[positions].all():
            raise ValueError('Pilot repeated or crossed eligible cells')
        name = output / 'selectors' / (dataset + '_context_positions.npy')
        with name.open('xb') as handle:
            np.save(handle, positions, allow_pickle=False)
        summaries.append({'dataset': dataset, 'family': entry['family'], 'role': entry['role'],
                          'material_rule': entry['material_rule'], 'fold_group': entry['fold_group'],
                          'sampling_unit_basis': basis, 'n_eligible_cells': int(eligible.sum()),
                          'n_selected_cells': len(positions), 'pool_counts': pool_counts,
                          'n_selected_sampling_blocks_not_patient_count': len(set(units[positions])),
                          'n_selected_expression_units': len(set(samples[positions])),
                          'selector_file': str(name.relative_to(output)), 'selector_sha256': sha256_file(name),
                          'source_native_cell_row_hash': hashlib.sha256(reader.rows(dataset)[positions].tobytes()).hexdigest(),
                          'selected_fine_strata_not_evaluation_truth': True})
    if sum(x['n_selected_cells'] for x in summaries) > spec['total_cell_cap']:
        raise ValueError('Whole pilot exceeds frozen total cap')
    gene_plan = model_gene_plan(reader, spec['sources'])
    save(output / 'config/model_features_and_heldouts_v1.json', gene_plan)
    save(output / 'outputs/pilot_input_inventory_v1.json', summaries)
    total_bytes = sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
    if total_bytes > args.output_cap_bytes:
        raise ValueError('Pilot metadata/selector byte cap exceeded')
    save(output / 'provenance/terminal_v1.json',
         {'status': 'PILOT_ROWS_AND_COMMON_MEASURED_FEATURES_FROZEN_NO_MODELS_FITTED',
          'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'n_sources': len(summaries), 'n_selected_cells': sum(x['n_selected_cells'] for x in summaries),
          'n_selected_T_NK': sum(x['pool_counts']['T_NK']['selected'] for x in summaries),
          'n_model_feature_candidates': len(gene_plan['common_model_feature_candidates']),
          'n_fold_groups': len({x['fold_group'] for x in spec['sources']}),
          'spec_sha256': sha256_file(spec_path), 'native_manifest_sha256': sha256_file(manifest),
          'feature_plan_sha256': sha256_file(output / 'config/model_features_and_heldouts_v1.json'),
          'inventory_sha256': sha256_file(output / 'outputs/pilot_input_inventory_v1.json'),
          'metadata_and_selector_bytes_before_terminal': total_bytes,
          'expression_arrays_copied_or_written': False, 'A_B_C_fitted': False,
          'alias_and_patient_independence_promoted': False,
          'wall_seconds': time.monotonic() - start})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for flag in ['spec', 'source-run', 'legacy-metadata', 'output']:
        parser.add_argument('--' + flag, required=True)
    parser.add_argument('--output-cap-bytes', type=int, default=16 * 1024 ** 2)
    args = parser.parse_args()
    try:
        main(args)
    except Exception as error:
        save(Path(args.output) / 'provenance/failure_v1.json',
             {'status': 'FAILED_NOT_ACCEPTED', 'error_class': type(error).__name__, 'error': str(error)})
        raise
