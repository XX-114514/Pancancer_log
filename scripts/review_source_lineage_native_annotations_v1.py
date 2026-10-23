"""Native review of source-by-lineage partitions; retains branch eligibility.

No counts/probability matrices are written. Reference proposals and descriptive
markers do not constitute independent validation or inferential DE. The pilot
is stratified, so its fractions are not source population abundance estimates.
"""
import argparse
import collections
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
from celltypist.models import Model

from annotation_representation_pilot_v1 import normalized_log_counts, save
from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_panel_aggregate_v1 import gather_csr_rows, STATE_PANELS
from native_annotation_review_v1 import (
    REVIEW_GENES, LINEAGE_PANELS, TYPE_PANELS, native_measured_genes,
    native_group_evidence, marker_calls, classify_group, descriptive_markers,
    predict_reference, model_fingerprint,
)


from review_lineage_anchor_sensitivity_v1 import EXTRA_GENES, anchored_calls

_base_marker_calls = marker_calls
REVIEW_GENES = sorted(set(REVIEW_GENES + EXTRA_GENES))


def marker_calls(detected, genes, measured):
    calls = _base_marker_calls(detected, genes, measured)
    calls.update(anchored_calls(detected, genes, measured, calls))
    return calls


def counts(values):
    return {str(k): int(v) for k, v in sorted(collections.Counter(values).items())}


def compact_reference(result):
    return {k: result[k] for k in ['best_code', 'unique_code', 'n_positive_labels', 'maximum_probability']}


def main(root):
    start = time.monotonic()
    cfg = json.loads((root / 'config/review_spec_v1.json').read_text())
    for item in cfg['bindings'].values():
        if sha256_file(item['path']) != item['sha256']:raise ValueError('Frozen review input changed')
    test = subprocess.run([sys.executable, '-B', str(root / 'tests/test_native_annotation_review_v1.py')],
                          capture_output=True, text=True, timeout=180)
    (root / 'provenance/test_stdout.txt').write_text(test.stdout)
    (root / 'provenance/test_stderr.txt').write_text(test.stderr)
    good = test.returncode == 0 and 'Ran 5 tests' in test.stderr and 'skipped' not in test.stderr.lower()
    save(root / 'provenance/numerical_tests_v1.json', {'actual_exit_code': test.returncode, 'n_tests': 5, 'no_skips': good})
    if not good:raise ValueError('Actual reference/numerical tests did not pass')
    anchor_test = subprocess.run([sys.executable, '-B', str(root / 'tests/test_lineage_anchor_sensitivity_v1.py')],
                                 capture_output=True, text=True, timeout=180)
    (root / 'provenance/anchor_test_stdout.txt').write_text(anchor_test.stdout)
    (root / 'provenance/anchor_test_stderr.txt').write_text(anchor_test.stderr)
    if anchor_test.returncode != 0 or 'Ran 3 tests' not in anchor_test.stderr or 'skipped' in anchor_test.stderr.lower():
        raise ValueError('Anchored lineage tests did not pass')
    reader = FullContextAtlas(cfg['manifest'])
    pilot, discovery = Path(cfg['pilot']), Path(cfg['discovery'])
    inventory = json.loads((pilot / 'outputs/pilot_input_inventory_v1.json').read_text())
    with np.load(discovery / 'outputs/derived_axis_and_model_mask.npz', allow_pickle=False) as f:
        sources = f['source'];modeled = f['modeled'];routes = f['inherited_route']
        native_routes = f['provisional_native_route']
        fine_fit = {'A': f['A_fine_discovery_eligible'], 'B': modeled, 'C': f['C_fine_discovery_eligible']}
        if any(len(v) != len(sources) for v in [native_routes, *fine_fit.values()]):
            raise ValueError('Lineage branch eligibility axes differ')
    with np.load(discovery / ('outputs/candidate_partitions_' + str(cfg['seed']) + '.npz'), allow_pickle=False) as f:
        partitions = {m: f[m] for m in ['A', 'B', 'C']}
    if sum(row['n_selected_cells'] for row in inventory) != len(sources):
        raise ValueError('Frozen inventory/partition count differs')
    models = {};fingerprints = {};model_metadata = {}
    for name, item in cfg['models'].items():
        if sha256_file(item['path']) != item['sha256']:raise ValueError('Reference model changed')
        model = Model.load(item['path']);models[name] = model;fingerprints[name] = model_fingerprint(model)
        model_metadata[name] = {'file_sha256': item['sha256'], 'n_features': len(model.features),
                                'labels': list(model.cell_types), 'description': model.description,
                                'reference_training_overlap_with_pilot': 'not_ruled_out',
                                'independent_annotation_truth': False}
    save(root / 'outputs/reference_models_v1.json', model_metadata)
    index = [];offset = 0
    for item in inventory:
        ds = item['dataset'];selector_path = pilot / item['selector_file']
        if sha256_file(selector_path) != item['selector_sha256']:raise ValueError('Frozen selector changed')
        positions = np.load(selector_path, allow_pickle=False);n = len(positions)
        sl = slice(offset, offset + n)
        if not np.all(sources[sl] == ds):raise ValueError('Frozen source order differs')
        inherited = reader.frozen_labels(ds)
        oldfine = inherited['fine'][positions].astype(str);oldmajor = inherited['major'][positions].astype(str)
        units = reader.frozen_expression_units(ds)[positions].astype(str)
        names = reader.field_names(ds)
        author_field = 'ev_author_label' if reader.is_legacy(ds) else (reader.primary.partitions[ds].get('annotation', {}).get('field') or 'source_author_label')
        author = reader.field(ds, author_field)[positions].astype(str) if author_field in names else np.full(n, '', dtype='U1')
        author_origin = author_field if author_field in names else 'missing_no_author_field'
        author_present = ~np.isin(np.char.lower(author), ['', 'nan', 'none', 'na', 'unknown', 'unassigned'])
        genes, columns = native_measured_genes(reader, ds)
        gene_pos = {g: j for j, g in enumerate(genes)}
        lookup = reader.feature_lookup(ds)
        measured = np.array([len(lookup.get(g, set())) == 1 for g in REVIEW_GENES])
        review_columns = np.array([next(iter(lookup[g])) if measured[j] else 0 for j, g in enumerate(REVIEW_GENES)])
        if len(set(review_columns[measured])) != measured.sum():raise ValueError('Review marker physical aliases')
        local = {}
        for method in ['A', 'B', 'C']:
            labels, codes = np.unique(partitions[method][sl], return_inverse=True)
            k = len(labels)
            local[method] = {'labels': labels, 'codes': codes, 'n': np.zeros(k, np.int64),
                             'cp_sum': np.zeros((k, len(genes))), 'log_sum': np.zeros((k, len(genes))),
                             'detected': np.zeros((k, len(genes))),
                             'review_log_sum': np.zeros((k, len(REVIEW_GENES))),
                             'review_detected': np.zeros((k, len(REVIEW_GENES))),
                             'panel_positive': {p: np.zeros(k, np.int64) for p in {**LINEAGE_PANELS, **TYPE_PANELS}}}
        source_models = ['immune_low', 'immune_high', cfg['source_organ_model'][ds]]
        predictions = {name: {'best_code': np.full(n, -1, np.int16), 'unique_code': np.full(n, -1, np.int16),
                               'n_positive_labels': np.full(n, -1, np.int16), 'maximum_probability': np.full(n, np.nan, np.float32)}
                       for name in source_models}
        reference_coverage = {};immune_eligible = np.zeros(n, bool);positive_RNA = np.zeros(n, bool)
        z = reader.store(ds)
        key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        nnz = 0
        for low in range(0, n, cfg['read_block_rows']):
            high = min(n, low + cfg['read_block_rows']);chunk = slice(low, high)
            raw = gather_csr_rows(z[key], reader.rows(ds)[positions[chunk]])
            nnz += raw.nnz
            totals = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
            x = raw[:, columns]
            log, positive = normalized_log_counts(x, totals);positive_RNA[chunk] = positive
            marker_raw = raw[:, review_columns].toarray();marker_raw[:, ~measured] = 0
            detected = marker_raw > 0
            scale = np.zeros(len(totals));scale[positive] = 10000 / totals[positive]
            review_log = np.log1p(marker_raw * scale[:, None])
            calls = marker_calls(detected, REVIEW_GENES, measured)
            native_immune = np.logical_or.reduce([calls[p] for p in ['T', 'NK_related', 'B', 'Plasma', 'Myeloid', 'DC_related', 'Neutrophil_related', 'Mast']])
            prior_immune = np.isin(routes[sl][chunk], ['T_NK', 'B', 'Plasma', 'Myeloid', 'Neutrophil', 'Mast'])
            eligible = positive & (native_immune | prior_immune);immune_eligible[chunk] = eligible
            for name in source_models:
                scope = eligible if name.startswith('immune_') else positive
                pred = predict_reference(log[scope], genes, models[name], cfg['minimum_model_feature_fraction'])
                reference_coverage[name] = {k: pred[k] for k in ['feature_fraction', 'n_model_features_measured', 'status']}
                for field, value in compact_reference(pred).items():predictions[name][field][low + np.flatnonzero(scope)] = value
            for method, agg in local.items():
                codes = agg['codes'][chunk];k = len(agg['labels'])
                stats = native_group_evidence(x, genes, totals, codes, k)
                for field in ['n', 'cp_sum', 'log_sum', 'detected']:agg[field] += stats[field]
                np.add.at(agg['review_log_sum'], codes[positive], review_log[positive])
                np.add.at(agg['review_detected'], codes[positive], detected[positive])
                for panel, value in calls.items():np.add.at(agg['panel_positive'][panel], codes[positive], value[positive])
        if not np.array_equal(modeled[sl], positive_RNA):raise ValueError('Original native RNA model mask changed')
        side = {'positive_RNA': positive_RNA, 'immune_reference_scope': immune_eligible}
        side.update({name + '__' + field: value for name, arrays in predictions.items() for field, value in arrays.items()})
        sidepath = root / 'outputs' / (ds + '_reference_proposals_v1.npz')
        with sidepath.open('xb') as f:np.savez_compressed(f, **side)
        records = []
        for method, agg in local.items():
            reviewed = []
            for j, label in enumerate(agg['labels']):
                member = agg['codes'] == j;n_group = int(member.sum());nrna = int(agg['n'][j])
                fractions = {p: float(v[j] / nrna) if nrna else 0. for p, v in agg['panel_positive'].items()}
                detection = {g: float(agg['review_detected'][j, k] / nrna) for k, g in enumerate(REVIEW_GENES) if measured[k] and nrna}
                references = {}
                for name in source_models:
                    codes = predictions[name]['unique_code'][member]
                    reference_counts = counts([str(models[name].cell_types[c]) for c in codes if c >= 0])
                    pos = predictions[name]['n_positive_labels'][member]
                    references[name] = {'unique_label_counts': reference_counts, 'unassigned': int((pos == 0).sum()),
                                        'multiple_labels': int((pos > 1).sum()), 'not_evaluable_or_out_of_scope': int((pos < 0).sum()),
                                        'not_independent_truth': True}
                review = classify_group(fractions, detection, references['immune_low']['unique_label_counts'], nrna)
                if label < 0:review.update(local_major='unresolved', local_celltype='unresolved', confidence='unresolved', annotation_evidence_status='no_frozen_candidate_partition')
                stable = hashlib.sha256((ds + '\0' + method + '\0' + str(int(label)) + '\0' + str(cfg['seed']) + '\0' + item['selector_sha256'] + '\0' + cfg['bindings']['partitions']['sha256']).encode()).hexdigest()
                unit_counts = counts(units[member])
                record = {'source': ds, 'method': method, 'seed': cfg['seed'], 'partition_code': int(label),
                          'source_local_group_id': ds + ':' + method + ':' + str(int(label)), 'group_sha256': stable,
                          'n_pilot_cells': n_group, 'n_native_RNA_evaluable': nrna,
                          'n_expression_units_not_patients': len(unit_counts),
                          'n_units_at_least20_group_cells': sum(v >= 20 for v in unit_counts.values()),
                          'pilot_sampling_not_population_abundance': True, 'annotation_stage': 'source_lineage_native_review_proposal_v1',
                          'provisional_native_route_counts_not_truth': counts(native_routes[sl][member]),
                          'n_cells_with_fine_discovery_fit': int(fine_fit[method][sl][member].sum()),
                          'n_cells_in_coarse_only_branch': int((~fine_fit[method][sl][member] & positive_RNA[member]).sum()),
                          'n_author_label_available': int(author_present[member].sum()),
                          'inherited_major_counts_not_truth': counts(oldmajor[member]),
                          'inherited_fine_crosswalk_not_truth': counts(oldfine[member]),
                          'author_evidence_counts_not_truth': counts(author[member]), 'author_field_origin': author_origin,
                          'native_multi_marker_fraction': fractions, 'native_review_gene_detection': detection,
                          'reference_proposals': references, **review}
                record['state_evidence_separate_from_celltype'] = {
                    name: {'n_measured_genes': int(sum(measured[REVIEW_GENES.index(g)] for g in panel)),
                           'mean_native_logCP10K': float(np.mean([agg['review_log_sum'][j, REVIEW_GENES.index(g)] / nrna for g in panel if measured[REVIEW_GENES.index(g)]]))
                           if nrna and any(measured[REVIEW_GENES.index(g)] for g in panel) else None,
                           'functional_diagnosis': False} for name, panel in STATE_PANELS.items()}
                reviewed.append(record)
            for j, record in enumerate(reviewed):
                parent = record['local_major']
                others = [k for k, v in enumerate(reviewed) if k != j and parent != 'unresolved' and v['local_major'] == parent and agg['n'][k] > 0]
                contrast = 'other_source_pilot_groups_in_same_native_supported_parent'
                if not others:
                    others = [k for k in range(len(reviewed)) if k != j and agg['n'][k] > 0]
                    contrast = 'other_source_pilot_groups_all_parents_no_same_parent_comparator'
                nrna = agg['n'][j];denominator = agg['n'][others].sum()
                if nrna and denominator:
                    markers = descriptive_markers(genes, agg['cp_sum'][j] / nrna, agg['detected'][j] / nrna,
                                                  agg['cp_sum'][others].sum(axis=0) / denominator,
                                                  agg['detected'][others].sum(axis=0) / denominator)
                else:markers = []
                record.update(marker_candidates=markers, marker_origin='own_source_pilot_group',
                              marker_comparator=contrast if denominator else 'not_available',
                              marker_evidence_type='own_group_expression_specificity_not_inferential_DEG',
                              valid_replicate_aware_DEG_available_in_this_run=False,
                              display_number=j + 1, display_number_is_biological_type=False,
                              display_marker=markers[0]['gene'] if markers else 'marker_unresolved')
                record['display_label'] = record['local_celltype'] + '_' + record['display_marker'] + '_' + str(j + 1)
                records.append(record)
        output = root / 'outputs' / (ds + '_native_annotation_review_v1.json.gz')
        payload = {'source': ds, 'n_pilot_cells': n, 'n_uniquely_measured_genes': len(genes),
                   'review_gene_measurement': dict(zip(REVIEW_GENES, measured.tolist())),
                   'reference_coverage': reference_coverage, 'groups': records}
        with gzip.open(output, 'xt') as f:json.dump(payload, f, ensure_ascii=False, allow_nan=False)
        summary = {'source': ds, 'n_pilot_cells': n, 'n_native_nnz_read': nnz,
                   'n_immune_reference_scope': int(immune_eligible.sum()), 'output': output.name,
                   'output_sha256': sha256_file(output), 'reference_proposal_sha256': sha256_file(sidepath),
                   'methods': {m: {'n_source_groups': sum(r['method'] == m for r in records),
                                   'celltype_proposal_cell_counts': dict(collections.Counter({}))}
                               for m in ['A', 'B', 'C']}}
        for m in summary['methods']:
            acc = collections.Counter()
            for r in records:
                if r['method'] == m:acc[r['local_celltype']] += r['n_pilot_cells']
            summary['methods'][m]['celltype_proposal_cell_counts'] = dict(acc)
        index.append(summary);offset += n
        save(root / 'checkpoints' / (ds + '_complete.json'), {'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                                                          'source': ds, 'n_cells': n, 'native_counts_copied': False})
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) > cfg['output_cap_bytes']:
            raise ValueError('Review sidecar cap exceeded')
    if offset != len(sources):raise ValueError('Not all frozen pilot sources processed')
    if any(model_fingerprint(models[k]) != v for k, v in fingerprints.items()):raise ValueError('Reference model mutated')
    save(root / 'outputs/source_annotation_review_index_v1.json', index)
    save(root / 'provenance/terminal_v1.json', {'status': 'SOURCE_LINEAGE_NATIVE_REVIEW_PROPOSALS_COMPLETE_NOT_NEW_AUTHORITY',
         'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_sources': len(index),
         'n_pilot_cells': offset, 'seed': cfg['seed'], 'n_tests_actual_no_skips': 8,
         'index_sha256': sha256_file(root / 'outputs/source_annotation_review_index_v1.json'),
         'wall_seconds': time.monotonic() - start, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
         'reference_models_immutable': True, 'reference_predictions_independent_truth': False,
         'native_counts_copied': False, 'all_full_archive_rows_annotated': False, 'clinical_or_CM_inputs_used': False,
         'method_winner_selected': False, 'replicate_aware_DEG_computed': False})


if __name__ == '__main__':
    p = argparse.ArgumentParser();p.add_argument('--output', required=True)
    root = Path(p.parse_args().output)
    try:main(root)
    except Exception as error:
        save(root / 'provenance/failure_v1.json', {'error_class': type(error).__name__, 'error': str(error)})
        raise
