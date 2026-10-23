"""Independently audit pilot selectors against native unit keys and material."""
import argparse
import collections
import datetime
import json
from pathlib import Path

import numpy as np
import pandas as pd

from full_context_atlas_v1 import FullContextAtlas, sha256_file


def fair_allocation_check(sizes, observed, target):
    """Independent integer water-level characterization; no RNG replay."""
    sizes = np.asarray(sizes, dtype=np.int64)
    observed = np.asarray(observed, dtype=np.int64)
    target = min(int(target), int(sizes.sum()))
    if observed.shape != sizes.shape or (observed < 0).any() or (observed > sizes).any() or observed.sum() != target:
        raise ValueError('Sampling allocation changed population/budget')
    low, high = 0, int(sizes.max()) + 1 if len(sizes) else 1
    while low + 1 < high:
        middle = (low + high) // 2
        if np.minimum(sizes, middle).sum() <= target:
            low = middle
        else:
            high = middle
    base = np.minimum(sizes, low)
    if np.any(observed < base) or np.any(observed > np.minimum(base + 1, sizes)):
        raise ValueError('Declared unit balance not met')


def main(args):
    output = Path(args.output)
    spec = json.loads((output / 'config/pilot_spec_v1.json').read_text())
    original = json.loads(Path(args.original_spec).read_text())
    if spec != original:
        raise ValueError('Namespace repair changed the frozen scientific specification')
    supervisor = json.loads((output / 'provenance/supervisor_terminal_v1.json').read_text())
    if supervisor['actual_child_exit_code'] != 0 or supervisor['wall_timeout']:
        raise ValueError('Pilot freezing process did not finish successfully')
    terminal = json.loads((output / 'provenance/terminal_v1.json').read_text())
    inventory_path = output / 'outputs/pilot_input_inventory_v1.json'
    feature_path = output / 'config/model_features_and_heldouts_v1.json'
    if sha256_file(inventory_path) != terminal['inventory_sha256'] or sha256_file(feature_path) != terminal['feature_plan_sha256']:
        raise ValueError('Changed completed pilot inventory')
    request = json.loads((output / 'config/supervisor_request_v1.json').read_text())
    for name, digest in request['script_SHA256'].items():
        if sha256_file(output / name) != digest:
            raise ValueError('Changed executed source script')
    reader = FullContextAtlas(Path(args.source_run) / 'config/full_context_collection_v1.json')
    legacy = pd.read_csv(args.legacy_metadata, sep='\t', dtype=str, keep_default_na=False).set_index('sampleID')
    breast = pd.read_csv(spec['bindings']['official_breast_map']['path'], sep='\t', dtype=str, keep_default_na=False).set_index('sampleID')
    inventory = json.loads(inventory_path.read_text())
    source_entries = {entry['dataset']: entry for entry in spec['sources']}
    if len(inventory) != len(source_entries):
        raise ValueError('Pilot source axis differs')
    checked = []
    native_rows_by_store = collections.defaultdict(set)
    for row in inventory:
        dataset = row['dataset']
        entry = source_entries[dataset]
        path = output / row['selector_file']
        if sha256_file(path) != row['selector_sha256']:
            raise ValueError('Changed row selector')
        selected = np.load(path, allow_pickle=False)
        all_rows = reader.rows(dataset)
        if (selected.dtype.kind not in 'iu' or len(selected) != len(np.unique(selected)) or
                len(selected) != row['n_selected_cells'] or
                (len(selected) and (selected.min() < 0 or selected.max() >= len(all_rows)))):
            raise ValueError('Invalid actual pilot rows')
        count = len(all_rows)
        if reader.is_legacy(dataset):
            # Independent axis: the producer used NPZ category codes, not this field.
            sample = reader.field(dataset, 'sample_uid').astype(str)
            if not np.array_equal(sample, reader.frozen_expression_units(dataset).astype(str)):
                raise ValueError('Native sample_uid/frozen-roster correspondence failed')
            meta = legacy.loc[sample]
            if entry['material_rule'] == 'official_breast_tumor_Total':
                raw = breast.loc[sample]
                eligible = (raw.official_material_role.eq('tumor_mammary_unresolved_repeat_specimens') &
                            raw.cell_population.eq('Total')).to_numpy()
                blocks = ('author_token_not_patient_claim::' + raw.author_patient_token).to_numpy()
            else:
                eligible = (meta.sample_class.eq('primary_or_unspecified_tumor').to_numpy()
                            if entry['material_rule'] == 'frozen_primary_or_unspecified' else np.ones(count, dtype=bool))
                blocks = np.asarray(['sample_fallback::' + value for value in sample], dtype=object)
        else:
            sample = reader.field(dataset, 'library_id').astype(str)
            eligible = np.ones(count, dtype=bool)
            if 'tissue_equals' in entry:
                eligible &= reader.field(dataset, 'tissue').astype(str) == entry['tissue_equals']
            fields = reader.field_names(dataset)
            for name in ['ex_vivo_stimulation', 'annotation_exclude_from_primary', 'residual_doublet_flag']:
                if name in fields:
                    eligible &= ~np.isin(np.char.lower(reader.field(dataset, name).astype(str)), ['true', '1'])
            verified = np.isin(np.char.lower(reader.field(dataset, 'patient_identity_verified').astype(str)), ['true', '1'])
            blocks = np.asarray(['library_fallback::' + value for value in sample], dtype=object)
            if verified.any():
                patient_blocks = reader.field(dataset, 'patient_block').astype(str)
                blocks[verified] = np.asarray(['reviewed_block::' + value for value in patient_blocks[verified]], dtype=object)
        if not eligible[selected].all() or int(eligible.sum()) != row['n_eligible_cells']:
            raise ValueError('Actual cells violate material/QC gates')
        labels = reader.frozen_labels(dataset)
        major = np.asarray(labels['major']).astype(str)
        masks = {'T_NK': np.isin(major, ['T cell', 'T/NK lymphoid', 'T', 'NK', 'T_NK_boundary']),
                 'Myeloid': np.isin(major, ['Myeloid', 'Monocyte', 'Macrophage']),
                 'Epithelial': np.isin(major, ['Epithelial cell', 'Epithelial'])}
        masks['Other_context_or_unresolved'] = ~(masks['T_NK'] | masks['Myeloid'] | masks['Epithelial'])
        for pool, mask in masks.items():
            universe = eligible & mask
            chosen = selected[mask[selected]]
            if len(chosen) != min(int(universe.sum()), spec['per_source_caps'][pool]):
                raise ValueError('Pool budget not preserved')
            names, codes = np.unique(blocks[universe].astype(str), return_inverse=True)
            sizes = np.bincount(codes, minlength=len(names))
            values = pd.Series(blocks[chosen]).value_counts().reindex(names, fill_value=0).to_numpy(dtype=np.int64)
            fair_allocation_check(sizes, values, spec['per_source_caps'][pool])
        physical = str(reader.store(dataset).store.path)
        rows = set(map(int, all_rows[selected]))
        if native_rows_by_store[physical] & rows:
            raise ValueError('Pilot contains repeated physical cells across sources')
        native_rows_by_store[physical].update(rows)
        checked.append({'dataset': dataset, 'material_QC_unit_and_pool_gates_independently_pass': True,
                        'n_selected_cells': len(selected), 'fold_group': entry['fold_group']})
    features = json.loads(feature_path.read_text())
    genes = features['common_model_feature_candidates']
    heldout = {gene for panels in features['heldout_RNA_readouts'].values() for values in panels.values() for gene in values}
    if len(genes) != len(set(genes)) or set(genes) & heldout:
        raise ValueError('Common model genes contain duplicate or held-out marker')
    for dataset in source_entries:
        lookup = reader.feature_lookup(dataset)
        if any(len(lookup.get(gene, set())) != 1 for gene in genes):
            raise ValueError('Model feature is unmeasured or ambiguous')
        if len({next(iter(lookup[gene])) for gene in genes}) != len(genes):
            raise ValueError('Two model names map to one physical feature')
    if any(p.suffix in {'.h5ad', '.npz', '.h5'} for p in output.rglob('*') if p.is_file()) or list(output.rglob('*.zarr')):
        raise ValueError('Unexpected expression object inside pilot index')
    receipt = {'status': 'INDEPENDENT_PILOT_INPUT_MATERIAL_NAMESPACE_AND_FEATURE_AUDIT_PASS',
               'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'auditor_sha256': sha256_file(__file__),
               'actual_freeze_child_exit_code': supervisor['actual_child_exit_code'],
               'same_scientific_specification_as_failed_namespace_attempt': True,
               'n_sources': len(checked), 'n_selected_cells': sum(x['n_selected_cells'] for x in checked),
               'n_unique_common_measured_model_features': len(genes),
               'selected_material_and_unit_allocations_independently_checked': True,
               'all_model_features_measured_unique_and_disjoint_from_heldout_readouts': True,
               'expression_arrays_copied': False, 'A_B_C_model_results_accepted': False,
               'source_checks': checked, 'terminal_sha256': sha256_file(output / 'provenance/terminal_v1.json')}
    with (output / 'provenance/independent_input_audit_v1.json').open('x') as stream:
        json.dump(receipt, stream, indent=2)
        stream.write('\n')
    print(json.dumps(receipt))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['output', 'source-run', 'legacy-metadata', 'original-spec']:
        parser.add_argument('--' + name, required=True)
    main(parser.parse_args())
