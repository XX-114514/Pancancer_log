"""Stream native evidence and scope-qualified cached reference proposals.

All selected native archive rows are represented. This stage is evidence
extraction, not final annotation, malignant status, model eligibility or release.
"""
import argparse
import collections
import datetime
import gc
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np
from celltypist.models import Model

from annotation_representation_pilot_v1 import ROUTES, normalized_log_counts, save
from full_context_atlas_v1 import sha256_file
from native_archive_collection_v1 import NativeArchiveCollection, row_axis_sha256, text_axis_sha256
from native_annotation_review_v1 import (REVIEW_GENES, native_measured_genes, marker_calls,
    predict_reference, model_fingerprint)
from native_lineage_routing_v1 import FAMILIES, native_family_support
from native_panel_aggregate_v1 import STATE_PANELS, gather_csr_rows
from review_lineage_anchor_sensitivity_v1 import EXTRA_GENES, anchored_calls


GENES = sorted(set(REVIEW_GENES + EXTRA_GENES))
IMMUNE_ROUTES = ['T_NK', 'B', 'Plasma', 'Myeloid', 'Neutrophil', 'Mast']
STATUS_FIELDS = ['residual_doublet_flag', 'annotation_exclude_from_primary', 'preintegration_cell_eligible',
                 'ambient_high_loss_review', 'analysis_correction_applied', 'ex_vivo_stimulation']


def reference_scopes(calls, inherited_route, positive, organ_available):
    immune_native = np.logical_or.reduce([calls[k] for k in
        ['T', 'NK_related', 'B', 'Plasma', 'Myeloid', 'DC_related', 'Neutrophil_related', 'Mast']])
    immune = positive & (immune_native | np.isin(inherited_route, IMMUNE_ROUTES))
    nonimmune_native = np.logical_or.reduce([calls[k] for k in
        ['Endothelial', 'Lymphatic_related', 'Fibroblast', 'Mural', 'Epithelial',
         'Glial_Oligodendrocyte_related', 'Astrocyte_related', 'Neuron_related']])
    # Organ references are additional evidence for native nonimmune candidates.
    # Mixed-lineage rows remain flagged; no reference assigns final identity.
    organ = positive & nonimmune_native if organ_available else np.zeros(len(positive), bool)
    return immune, organ


def categorical(values):
    names, codes = np.unique(np.asarray(values).astype(str), return_inverse=True)
    return names.tolist(), codes.astype(np.int32)


def pack_support(calls, names):
    if len(names) > 32:
        raise ValueError('Panel bitmask capacity exceeded')
    bits = np.zeros(len(next(iter(calls.values()))), np.uint32)
    for j, name in enumerate(names):
        bits |= np.asarray(calls[name], np.uint32) << j
    return bits


def main(root):
    started = time.monotonic()
    cfg = json.loads((root / 'config/full_archive_reference_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']:
            raise ValueError('Frozen full-archive evidence input changed')
    def checkpoint(stage, **extra):
        save(root / 'checkpoints' / (str(time.time_ns()) + '.json'), {
            'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'stage': stage,
            'wall_seconds': time.monotonic() - started, **extra})
    test_records = []
    for filename, number in [('test_native_annotation_review_v1.py', 5),
                             ('test_lineage_anchor_sensitivity_v1.py', 3),
                             ('test_full_archive_reference_scopes_v1.py', 3)]:
        result = subprocess.run([sys.executable, '-B', str(root / 'tests' / filename)],
                                capture_output=True, text=True, timeout=180)
        (root / 'provenance' / (filename + '.stdout.txt')).write_text(result.stdout)
        (root / 'provenance' / (filename + '.stderr.txt')).write_text(result.stderr)
        if result.returncode != 0 or ('Ran %d tests' % number) not in result.stderr or 'skipped' in result.stderr.lower():
            raise ValueError('Full archive reference counterexamples failed')
        test_records.append({'test': filename, 'n_tests': number, 'actual_exit_code': 0, 'no_skips': True})
    save(root / 'provenance/numerical_tests_v1.json', {'tests': test_records, 'n_tests': 11, 'actual_exit_code': 0})
    reader = NativeArchiveCollection(cfg['archive_manifest'])
    models, fingerprints, model_meta = {}, {}, {}
    for name, item in cfg['models'].items():
        if sha256_file(item['path']) != item['sha256']:
            raise ValueError('Cached reference model changed')
        model = Model.load(item['path']); models[name] = model; fingerprints[name] = model_fingerprint(model)
        model_meta[name] = {'file_sha256': item['sha256'], 'fingerprint': fingerprints[name],
            'labels': list(model.cell_types), 'n_features': len(model.features), 'description': model.description,
            'reference_training_overlap_unresolved': True, 'probability_is_not_calibrated_identity_confidence': True}
    save(root / 'outputs/reference_models_v1.json', model_meta)
    pilot = Path(cfg['pilot']); prior = Path(cfg['accepted_pilot_review'])
    inventory = {e['dataset']: e for e in json.loads((pilot / 'outputs/pilot_input_inventory_v1.json').read_text())}
    previous_index = {e['source']: e for e in json.loads((prior / 'outputs/source_annotation_review_index_v1.json').read_text())}
    index, total_rows, total_nnz, compared_values = [], 0, 0, 0
    for ds in cfg['source_order']:
        began = time.monotonic(); rows = reader.rows(ds); n = len(rows)
        part = reader.archive_partitions[ds]
        identity_sha = text_axis_sha256(reader.native_ids(ds))
        if identity_sha != part['native_ids_sha256']:
            raise ValueError('Native identity axis changed')
        inherited = reader.frozen_labels(ds)
        inherited_route = np.array([ROUTES.get(str(v), 'unresolved_context') for v in inherited['major']])
        names = set(reader.field_names(ds)); author_field = 'ev_author_label' if reader.is_legacy(ds) else (
            reader.primary.partitions[ds].get('annotation', {}).get('field') or 'source_author_label')
        author = reader.field(ds, author_field).astype(str) if author_field in names else np.full(n, '', dtype='<U1')
        author_nonmissing = ~np.isin(np.char.lower(author), ['', 'nan', 'none', 'na', 'unknown', 'unassigned'])
        metadata_fields, arrays = {}, {'native_row': rows, 'native_total_RNA': np.zeros(n, np.float64),
            'positive_RNA': np.zeros(n, bool), 'immune_reference_scope': np.zeros(n, bool),
            'organ_reference_scope': np.zeros(n, bool), 'native_family_bitmask': np.zeros(n, np.uint32),
            'native_panel_bitmask': np.zeros(n, np.uint32),
            'native_gene_detection_packbits': np.zeros((n, (len(GENES) + 7) // 8), np.uint8),
            'state_native_logCP10K': np.full((n, len(STATE_PANELS)), np.nan, np.float32)}
        for field in STATUS_FIELDS:
            if field not in names:
                metadata_fields[field] = {'status': 'field_absent_not_false'}
                continue
            dictionary, codes = categorical(reader.field(ds, field))
            arrays['carried_status__' + field] = codes
            metadata_fields[field] = {'status': 'carried_not_automatically_model_exclusion', 'dictionary': dictionary,
                                     'counts': np.bincount(codes, minlength=len(dictionary)).tolist()}
        membership = reader.membership(ds)
        dictionary, codes = categorical(membership['membership_status'])
        arrays['archive_membership_code'] = codes
        membership_meta = {'dictionary': dictionary, 'counts': np.bincount(codes, minlength=len(dictionary)).tolist()}
        organ = cfg['source_organ_model'].get(ds)
        source_models = ['immune_low', 'immune_high'] + ([organ] if organ else [])
        for model in source_models:
            for field in ['best_code', 'unique_code', 'n_positive_labels']:
                arrays[model + '__' + field] = np.full(n, -1, np.int16)
            arrays[model + '__maximum_probability'] = np.full(n, np.nan, np.float32)
        genes, columns = native_measured_genes(reader, ds)
        lookup = reader.feature_lookup(ds)
        measured = np.array([len(lookup.get(g, set())) == 1 for g in GENES])
        marker_columns = np.array([next(iter(lookup[g])) if measured[j] else 0 for j, g in enumerate(GENES)])
        if len(set(marker_columns[measured])) != int(measured.sum()):
            raise ValueError('Native review marker physical aliases')
        state_columns = [[GENES.index(g) for g in panel if measured[GENES.index(g)]] for panel in STATE_PANELS.values()]
        state_coverage = {name: {'measured_genes': [g for g in panel if measured[GENES.index(g)]],
                               'minimum_measured_genes': 3} for name, panel in STATE_PANELS.items()}
        z = reader.store(ds); key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        nnz, coverage, panel_order = 0, {}, None
        for low in range(0, n, cfg['read_block_rows']):
            high = min(n, low + cfg['read_block_rows']); sl = slice(low, high)
            raw = gather_csr_rows(z[key], rows[sl]); nnz += raw.nnz
            if not np.isfinite(raw.data).all() or (raw.data < 0).any() or not np.equal(raw.data, np.floor(raw.data)).all():
                raise ValueError('Native count domain changed')
            total = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
            log, positive = normalized_log_counts(raw[:, columns], total)
            marker_raw = raw[:, marker_columns].toarray(); marker_raw[:, ~measured] = 0
            detected = marker_raw > 0
            calls = marker_calls(detected, GENES, measured); calls.update(anchored_calls(detected, GENES, measured, calls))
            if panel_order is None:
                panel_order = sorted(calls)
            immune_scope, organ_scope = reference_scopes(calls, inherited_route[sl], positive, organ is not None)
            arrays['native_total_RNA'][sl] = total; arrays['positive_RNA'][sl] = positive
            arrays['immune_reference_scope'][sl] = immune_scope; arrays['organ_reference_scope'][sl] = organ_scope
            arrays['native_family_bitmask'][sl] = pack_support(native_family_support(calls), FAMILIES)
            arrays['native_panel_bitmask'][sl] = pack_support(calls, panel_order)
            arrays['native_gene_detection_packbits'][sl] = np.packbits(detected, axis=1, bitorder='little')
            scale = np.zeros(len(total)); scale[positive] = 10000 / total[positive]
            marker_log = np.log1p(marker_raw * scale[:, None])
            for j, positions in enumerate(state_columns):
                if len(positions) >= 3:
                    values = marker_log[:, positions].mean(axis=1); values[~positive] = np.nan
                    arrays['state_native_logCP10K'][sl, j] = values
            for model in source_models:
                scope = immune_scope if model.startswith('immune_') else organ_scope
                pred = predict_reference(log[scope], genes, models[model], cfg['minimum_model_feature_fraction'])
                record = {k: pred[k] for k in ['status', 'feature_fraction', 'n_model_features_measured']}
                if model in coverage and coverage[model] != record:
                    raise ValueError('Model feature coverage changed within source')
                coverage[model] = record
                selected = low + np.flatnonzero(scope)
                for field in ['best_code', 'unique_code', 'n_positive_labels', 'maximum_probability']:
                    arrays[model + '__' + field][selected] = pred[field]
            if low == 0 or high == n or (low // cfg['read_block_rows']) % 64 == 0:
                checkpoint('full_archive_native_reference_block', source=ds, n_source_rows_processed=high,
                           n_source_rows=n, completed_sources=len(index), completed_archive_rows=total_rows + high)
        replay = {'status': 'not_in_frozen_pilot', 'prediction_values_replayed': 0}
        if ds in inventory:
            inv = inventory[ds]; old = previous_index[ds]
            if sha256_file(pilot / inv['selector_file']) != inv['selector_sha256']:
                raise ValueError('Pilot selector changed')
            native = reader.context_rows(ds)[np.load(pilot / inv['selector_file'], allow_pickle=False)]
            pos = np.searchsorted(rows, native); np.testing.assert_array_equal(rows[pos], native)
            prior_path = prior / 'outputs' / (ds + '_reference_proposals_v1.npz')
            if sha256_file(prior_path) != old['reference_proposal_sha256']:
                raise ValueError('Accepted pilot reference changed')
            with np.load(prior_path, allow_pickle=False) as f:
                np.testing.assert_array_equal(arrays['positive_RNA'][pos], f['positive_RNA'])
                np.testing.assert_array_equal(arrays['immune_reference_scope'][pos], f['immune_reference_scope'])
                values = 0
                for model in source_models:
                    scope = arrays['immune_reference_scope'][pos] if model.startswith('immune_') else arrays['organ_reference_scope'][pos]
                    for field in ['best_code', 'unique_code', 'n_positive_labels', 'maximum_probability']:
                        first, second = arrays[model + '__' + field][pos][scope], f[model + '__' + field][scope]
                        if field == 'maximum_probability':
                            np.testing.assert_allclose(first, second, atol=2e-5, rtol=2e-5, equal_nan=True)
                        else:
                            np.testing.assert_array_equal(first, second)
                        values += len(first)
                replay = {'status': 'accepted_pilot_predictions_reproduced_on_common_scope', 'prediction_values_replayed': values,
                          'organ_scope_restricted_to_native_nonimmune_candidates_in_this_stage': True}
                compared_values += values
        for name, model in models.items():
            if model_fingerprint(model) != fingerprints[name]:
                raise ValueError('Shared reference mutated during full archive scoring')
        output = root / 'outputs' / (ds + '_native_reference_evidence.npz')
        with output.open('xb') as f:
            np.savez_compressed(f, **arrays)
        with np.load(output, allow_pickle=False) as f:
            for name, value in arrays.items():
                np.testing.assert_array_equal(f[name], value)
        meta = {'source': ds, 'n_archive_rows': n, 'native_rows_sha256': row_axis_sha256(rows), 'native_ids_sha256': identity_sha,
            'archive_membership': membership_meta, 'author_field': author_field if author_field in names else None,
            'n_nonmissing_author_labels': int(author_nonmissing.sum()), 'author_labels_remain_in_native_source': True,
            'carried_status_fields': metadata_fields, 'review_genes': GENES, 'gene_measured': measured.tolist(),
            'detection_packbits_order': 'little', 'native_panel_bit_order': panel_order, 'native_family_bit_order': FAMILIES,
            'state_order': list(STATE_PANELS), 'state_measurement_coverage': state_coverage,
            'state_scores_descriptive_not_functional_diagnosis': True, 'source_reference_models': source_models,
            'reference_feature_coverage': coverage, 'reference_scope_policy': cfg['reference_scope_policy'],
            'organ_reference_unavailable_or_not_yet_qualified': organ is None, 'pilot_replay': replay,
            'native_count_origin': key, 'native_nnz_read': nnz, 'n_positive_RNA': int(arrays['positive_RNA'].sum()),
            'n_immune_reference_scope': int(arrays['immune_reference_scope'].sum()),
            'n_organ_reference_scope': int(arrays['organ_reference_scope'].sum()),
            'n_cells_without_reference_scope': int((arrays['positive_RNA'] & ~arrays['immune_reference_scope'] & ~arrays['organ_reference_scope']).sum()),
            'no_reference_scope_is_not_a_final_unknown_identity_call': True,
            'model_and_composition_eligibility_not_assigned': True, 'malignancy_not_inferred': True,
            'new_authoritative_labels': False, 'wall_seconds': time.monotonic() - began}
        meta_path = output.with_suffix('.json'); save(meta_path, meta)
        record = {'source': ds, 'n_archive_rows': n, 'n_native_nnz': nnz, 'arrays': output.name,
            'arrays_sha256': sha256_file(output), 'metadata': meta_path.name, 'metadata_sha256': sha256_file(meta_path),
            'pilot_prediction_values_replayed': replay['prediction_values_replayed']}
        save(root / 'outputs' / (ds + '_completion.json'), record); index.append(record)
        total_rows += n; total_nnz += nnz
        written = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
        if written > cfg['output_cap_bytes']:
            raise ValueError('Full archive reference sidecar cap exceeded')
        checkpoint('full_archive_native_reference_source_complete', source=ds, n_archive_rows=total_rows,
                   completed_sources=len(index), written_bytes=written)
        del arrays, inherited, inherited_route, author, membership, log, raw, marker_raw, detected, pred
        reader._fields.clear(); gc.collect()
    if set(cfg['source_order']) != set(reader.archive_partitions) or total_rows != sum(p['n_archive_rows'] for p in reader.archive_partitions.values()):
        raise ValueError('Full selected native archive not covered')
    save(root / 'outputs/full_archive_reference_index_v1.json', index)
    save(root / 'provenance/terminal_v1.json', {'status': 'FULL_NATIVE_ARCHIVE_REFERENCE_EVIDENCE_COMPLETE_NOT_FINAL_ANNOTATION',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_sources': len(index),
        'n_archive_rows': total_rows, 'native_nnz_read': total_nnz, 'pilot_prediction_values_replayed': compared_values,
        'index_sha256': sha256_file(root / 'outputs/full_archive_reference_index_v1.json'),
        'wall_seconds': time.monotonic() - started, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        'native_counts_copied': False, 'source_qualified_reference_expansion_and_local_groups_still_required': True,
        'new_authoritative_labels': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    args = parser.parse_args()
    try:
        main(Path(args.output))
    except Exception as error:
        save(Path(args.output) / 'provenance/failure_v1.json', {'status': 'FAILED_NOT_ACCEPTED',
             'error_class': type(error).__name__, 'error': str(error)})
        raise
