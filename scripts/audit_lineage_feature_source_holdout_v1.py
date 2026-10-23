"""Replay frozen native inputs, train-only features and heldout predictions.

No model is refitted. Full summaries are replayed; a size-selected query source
also receives direct-distance neighbor checks with boundary ties reported.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import resource
import time

import numpy as np
from scipy.spatial.distance import cdist

from annotation_representation_pilot_v1 import save, select_features
from full_context_atlas_v1 import sha256_file
from load_annotation_pilot_native_v1 import load_pilot
from lineage_holdout_diagnostics_v1 import VARIANTS, fold_members, training_feature_sets, paired_profiles


def same_record(first, second):
    if isinstance(first, dict):
        if set(first) != set(second):
            raise ValueError('Replayed record keys differ')
        for key in first:
            same_record(first[key], second[key])
    elif isinstance(first, list):
        if len(first) != len(second):
            raise ValueError('Replayed record length differs')
        for a, b in zip(first, second):
            same_record(a, b)
    elif isinstance(first, float):
        np.testing.assert_allclose(first, second, atol=1e-11, rtol=1e-10)
    elif first != second:
        raise ValueError('Replayed record value differs')


def saved_reference_fingerprint(values):
    digest = hashlib.sha256()
    digest.update('\0'.join(values['genes']).encode())
    for key in ['mean', 'std', 'loadings', 'C', 'Nr', 'sigma']:
        a = np.asarray(values[key])
        digest.update(str(a.dtype).encode()); digest.update(str(a.shape).encode()); digest.update(a.tobytes())
    return digest.hexdigest()


def main(root, producer):
    started = time.monotonic()
    request = json.loads((root / 'config/audit_spec_v1.json').read_text())
    for value in request['bindings'].values():
        if sha256_file(value['path']) != value['sha256']:
            raise ValueError('Frozen holdout audit dependency changed')
    terminal = json.loads((producer / 'provenance/supervisor_terminal_v1.json').read_text())
    if terminal['actual_child_exit_code'] != 0 or terminal['wall_timeout']:
        raise ValueError('Holdout producer not complete')
    cfg = json.loads((producer / 'config/lineage_source_holdout_spec_v1.json').read_text())
    for value in cfg['bindings'].values():
        if sha256_file(value['path']) != value['sha256']:
            raise ValueError('Producer input changed before audit')
    def checkpoint(stage, **extra):
        save(root / 'checkpoints' / (str(time.time_ns()) + '.json'), {
            'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'stage': stage, **extra})
    data = load_pilot(cfg['pilot'], cfg['manifest'], checkpoint, cfg['read_block_rows'])
    with np.load(Path(cfg['source_lineage']) / 'outputs/derived_axis_and_model_mask.npz', allow_pickle=False) as f:
        source, route, positive = f['source'], f['provisional_native_route'], f['modeled']
    np.testing.assert_array_equal(source, data['source']); np.testing.assert_array_equal(positive, data['native_total'] > 0)
    direct = np.full(len(source), 'unresolved', dtype='<U80')
    graph_index = json.loads((Path(cfg['graphs']) / 'outputs/graph_comparison_index_v1.json').read_text())
    for entry in graph_index:
        path = Path(cfg['graphs']) / 'outputs' / entry['axes']
        if sha256_file(path) != entry['axes_sha256']:
            raise ValueError('Direct-candidate input changed')
        with np.load(path, allow_pickle=False) as f:
            np.testing.assert_array_equal(f['pilot_row'], np.flatnonzero(positive & (route == entry['route'])))
            direct[f['pilot_row']] = f['direct_candidate']
    index = json.loads((producer / 'outputs/lineage_source_holdout_index_v1.json').read_text())
    programs = sorted(data['state']); targets = np.column_stack([data['state'][p] for p in programs])
    names, counts = np.unique(source, return_counts=True)
    numerical_source = names[np.argmin(counts)]
    all_profiles, feature_cache, cases = [], {}, []
    complete_route_folds = 0; feature_values = 0; max_error = 0.; finite_values = 0; boundary_ties = 0
    expected = {(fold, wanted) for fold in cfg['expected_fold_groups'] for wanted in cfg['routes']}
    if {(e['fold'], e['route']) for e in index} != expected or len(index) != len(expected):
        raise ValueError('Holdout fold/route index incomplete or duplicated')
    for entry in index:
        train, query = fold_members(source, data['fold'], positive, entry['fold'])
        tr, qu = train[route[train] == entry['route']], query[route[query] == entry['route']]
        if entry['status'] != 'complete':
            if len(tr) >= cfg['minimum_training_group'] and len(qu):
                raise ValueError('Eligible route silently omitted')
            continue
        folder = producer / entry['relative_folder']
        for name, digest in entry['files_SHA256'].items():
            if sha256_file(folder / name) != digest:
                raise ValueError('Holdout artifact changed')
        with np.load(folder / 'axes_and_readouts.npz', allow_pickle=False) as f:
            axes = {key: f[key] for key in f.files}
        np.testing.assert_array_equal(axes['training_pilot_rows'], tr)
        np.testing.assert_array_equal(axes['query_pilot_rows'], qu)
        anchors = axes['anchor_pilot_rows']
        if len(np.unique(anchors)) != len(anchors) or not np.isin(anchors, tr).all():
            raise ValueError('Anchor repeats or leaks query')
        for key, value in {'query_source': source[qu], 'anchor_source': source[anchors],
                'query_direct_candidate': direct[qu], 'anchor_direct_candidate': direct[anchors],
                'query_native_totals': data['native_total'][qu], 'query_RNA_readouts': targets[qu],
                'anchor_RNA_readouts': targets[anchors]}.items():
            np.testing.assert_array_equal(axes[key], value)
        if entry['fold'] not in feature_cache:
            feature_cache[entry['fold']] = select_features(data['x'][train], source[train], cfg['n_features'])
        features = training_feature_sets(data['x'], source, route, train, entry['route'], cfg['n_features'], feature_cache[entry['fold']])
        metadata = json.loads((folder / 'models_and_coverage.json').read_text())
        predictions, support, coordinates = {}, {}, {}
        for model in metadata['models']:
            variant = model['variant']; selected = features[variant]
            np.testing.assert_array_equal(model['selected_feature_indices'], selected); feature_values += len(selected)
            if model['selected_feature_names'] != [data['genes'][j] for j in selected]:
                raise ValueError('Feature name axis differs')
            if model['training_sources'] != sorted(set(source[tr])) or model['query_sources'] != sorted(set(source[qu])):
                raise ValueError('Model source qualification differs')
            with np.load(folder / (variant + '_reference.npz'), allow_pickle=False) as f:
                fingerprint = saved_reference_fingerprint(f)
                np.testing.assert_array_equal(f['feature_candidate_indices'], selected)
            if fingerprint != model['reference_fingerprint_before'] or fingerprint != model['reference_fingerprint_after']:
                raise ValueError('Stored reference differs from frozen mapping fingerprint')
            with np.load(folder / (variant + '_coordinates.npz'), allow_pickle=False) as f:
                coordinates[variant] = (f['anchor_embedding'], f['query_embedding'])
                support[variant] = f['query_supported']
            ref, mapped = coordinates[variant]
            if len(ref) != len(anchors) or len(mapped) != len(qu) or not np.isfinite(ref).all():
                raise ValueError('Stored coordinate axes differ')
            np.testing.assert_array_equal(np.isfinite(mapped).all(axis=1), support[variant])
        for scope_record in metadata['scopes']:
            scope = scope_record['scope']
            qm = np.ones(len(qu), bool) if scope == 'all_route_cells_not_pure_type' else direct[qu] == scope
            am = np.ones(len(anchors), bool) if scope == 'all_route_cells_not_pure_type' else direct[anchors] == scope
            if scope_record['n_query_cells'] != qm.sum() or scope_record['n_anchor_cells'] != am.sum():
                raise ValueError('Scope denominator differs')
            with np.load(folder / (scope + '_predictions.npz'), allow_pickle=False) as f:
                predictions[scope] = {v: f[v] for v in VARIANTS}
            for variant in VARIANTS:
                pred = predictions[scope][variant]
                outside = ~qm | ~support[variant]
                if np.isfinite(pred[outside]).any():
                    raise ValueError('Prediction crosses query scope')
                if am.sum() < cfg['minimum_anchor_cells']:
                    if np.isfinite(pred).any():
                        raise ValueError('Prediction below anchor threshold')
                    continue
                selected_query = np.flatnonzero(qm & support[variant] & (source[qu] == numerical_source))[:32]
                if not len(selected_query):
                    continue
                ref, mapped = coordinates[variant]
                distances = cdist(mapped[selected_query].astype(float), ref[am].astype(float))
                k = min(cfg['state_neighbors'], int(am.sum()))
                order = np.argsort(distances, axis=1, kind='stable')
                tied = (np.isclose(np.take_along_axis(distances, order[:, k-1:k], axis=1)[:, 0],
                                   np.take_along_axis(distances, order[:, k:k+1], axis=1)[:, 0], atol=1e-9, rtol=1e-7)
                        if k < distances.shape[1] else np.zeros(len(selected_query), bool))
                boundary_ties += int(tied.sum())
                chosen = order[~tied, :k]; d = np.take_along_axis(distances[~tied], chosen, axis=1)
                weight = 1. / np.maximum(d, 1e-6); weight /= weight.sum(axis=1, keepdims=True)
                replay = np.sum(targets[anchors[am]][chosen] * weight[:, :, None], axis=1)
                recorded = pred[selected_query[~tied]]
                np.testing.assert_array_equal(np.isfinite(replay), np.isfinite(recorded))
                valid = np.isfinite(replay)
                if valid.any():
                    error = float(np.max(np.abs(replay[valid] - recorded[valid]))); max_error = max(max_error, error)
                    np.testing.assert_allclose(replay[valid], recorded[valid], atol=3e-5, rtol=3e-5)
                    finite_values += int(valid.sum())
                cases.append({'route': entry['route'], 'scope': scope, 'variant': variant,
                              'n_query_rows': len(selected_query), 'n_boundary_ties_not_numerically_adjudicated': int(tied.sum())})
        replay = paired_profiles(source[qu], direct[qu], targets[qu], data['native_total'][qu], predictions,
            support, programs, cfg['minimum_query_evaluation'])
        same_record(replay, json.loads((folder / 'paired_source_profiles.json').read_text()))
        all_profiles.extend(dict(fold=entry['fold'], route=entry['route'], **r) for r in replay)
        complete_route_folds += 1
        checkpoint('new_route_holdout_audit_route_complete', fold=entry['fold'], route=entry['route'])
    same_record(all_profiles, json.loads((producer / 'outputs/paired_source_holdout_profiles_v1.json').read_text()))
    if finite_values == 0:
        raise ValueError('No independent numerical predictions verified')
    report = {'completed_route_folds': complete_route_folds, 'profile_records_replayed': len(all_profiles),
        'selected_feature_indices_replayed': feature_values, 'independent_query_source_selection': 'smallest_frozen_pilot_source',
        'independent_prediction_cases': cases, 'independent_finite_prediction_values': finite_values,
        'independent_prediction_max_abs_error': max_error, 'boundary_tied_query_rows_not_adjudicated': boundary_ties,
        'native_readouts_axes_and_all_training_feature_scopes_replayed': True,
        'all_minimal_reference_fingerprints_preserved': True, 'reference_models_refitted': False,
        'identity_accuracy_or_independent_patients_claimed': False}
    save(root / 'outputs/lineage_feature_holdout_audit_v1.json', report)
    save(root / 'provenance/terminal_v1.json', {
        'status': 'NEW_ROUTE_FEATURE_HOLDOUT_AXES_FEATURES_PROFILES_AND_DIRECT_PREDICTIONS_ACCEPTED',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'output_sha256': sha256_file(root / 'outputs/lineage_feature_holdout_audit_v1.json'),
        'completed_route_folds': complete_route_folds, 'profile_records_replayed': len(all_profiles),
        'independent_finite_prediction_values': finite_values, 'independent_prediction_max_abs_error': max_error,
        'boundary_tied_query_rows_not_adjudicated': boundary_ties,
        'wall_seconds': time.monotonic() - started, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        'new_authoritative_labels': False, 'native_counts_copied': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True); parser.add_argument('--producer', required=True)
    args = parser.parse_args()
    try:
        main(Path(args.output), Path(args.producer))
    except Exception as error:
        save(Path(args.output) / 'provenance/failure_v1.json', {'status': 'FAILED_NOT_ACCEPTED',
             'error_class': type(error).__name__, 'error': str(error)})
        raise
