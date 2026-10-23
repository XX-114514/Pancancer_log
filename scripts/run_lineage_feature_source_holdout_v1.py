"""Frozen new-route, training-only shared versus lineage-feature references.

Private bounded sidecars contain model coordinates and RNA program readouts,
never native count copies. This is conditional reference transfer, not a new
identity gold standard or independent clinical validation.
"""
import argparse
import datetime
import gc
import json
from pathlib import Path
import resource
import subprocess
import sys
import time

import numpy as np

from annotation_pilot_sampling_v1 import balanced_cell_selection, stable_seed
from annotation_representation_pilot_v1 import save, select_features
from full_context_atlas_v1 import sha256_file
from frozen_symphony_reference_v1 import load_symphonypy, build_reference, map_query, save_minimal_reference, reference_fingerprint
from load_annotation_pilot_native_v1 import load_pilot
from run_annotation_source_holdout_v1 import inverse_neighbor_prediction
from lineage_holdout_diagnostics_v1 import VARIANTS, fold_members, training_feature_sets, paired_profiles


def write_npz(path, **values):
    with path.open('xb') as f:
        np.savez_compressed(f, **values)
    with np.load(path, allow_pickle=False) as f:
        for key, value in values.items():
            np.testing.assert_array_equal(f[key], value)


def main(root):
    started = time.monotonic()
    cfg = json.loads((root / 'config/lineage_source_holdout_spec_v1.json').read_text())
    for binding in cfg['bindings'].values():
        if sha256_file(binding['path']) != binding['sha256']:
            raise ValueError('Frozen lineage holdout input changed')
    for path in cfg['accepted_dependencies']:
        terminal = json.loads((Path(path) / 'provenance/supervisor_terminal_v1.json').read_text())
        if terminal['actual_child_exit_code'] != 0 or terminal['wall_timeout']:
            raise ValueError('Lineage holdout dependency incomplete')
    def checkpoint(stage, **extra):
        save(root / 'checkpoints' / (str(time.time_ns()) + '.json'), {
            'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'stage': stage, 'wall_seconds': time.monotonic() - started, **extra})
    result = subprocess.run([sys.executable, '-B', str(root / 'tests/test_lineage_source_holdout_v1.py')],
                            capture_output=True, text=True, timeout=180)
    (root / 'provenance/test_stdout.txt').write_text(result.stdout)
    (root / 'provenance/test_stderr.txt').write_text(result.stderr)
    if result.returncode != 0 or 'Ran 6 tests' not in result.stderr or 'skipped' in result.stderr.lower():
        raise ValueError('Lineage holdout counterexamples failed')
    n_tests = 6
    if cfg.get('feature_universe_binding') is not None:
        contract_test = subprocess.run([sys.executable, '-B', str(root / 'tests/test_common_measured_feature_contract_v1.py')],
                                      capture_output=True, text=True, timeout=180)
        (root / 'provenance/measurement_test_stdout.txt').write_text(contract_test.stdout)
        (root / 'provenance/measurement_test_stderr.txt').write_text(contract_test.stderr)
        if contract_test.returncode != 0 or 'Ran 3 tests' not in contract_test.stderr or 'skipped' in contract_test.stderr.lower():
            raise ValueError('Measurement feature contract counterexamples failed')
        n_tests += 3
    save(root / 'provenance/numerical_tests_v1.json', {'actual_exit_code': 0, 'n_tests': n_tests, 'no_skips': True})
    symphony = load_symphonypy(cfg['wheel_path'], cfg['wheel_SHA256'])
    data = load_pilot(cfg['pilot'], cfg['manifest'], checkpoint, cfg['read_block_rows'], cfg.get('feature_universe_binding'))
    save(root / 'outputs/model_feature_universe_v1.json', {
        'n_candidate_genes': len(data['genes']), 'candidate_genes': data['genes'],
        'measurement_metadata_binding': cfg.get('feature_universe_binding'),
        'query_expression_used_for_feature_ranking': False, 'missing_genes_imputed': False})
    with np.load(Path(cfg['source_lineage']) / 'outputs/derived_axis_and_model_mask.npz', allow_pickle=False) as f:
        source, route, positive = f['source'], f['provisional_native_route'], f['modeled']
        np.testing.assert_array_equal(source, data['source'])
        np.testing.assert_array_equal(positive, data['native_total'] > 0)
    direct = np.full(len(source), 'unresolved', dtype='<U80')
    graph_index = json.loads((Path(cfg['graphs']) / 'outputs/graph_comparison_index_v1.json').read_text())
    for wanted in cfg['routes']:
        entry = next(e for e in graph_index if e['route'] == wanted)
        path = Path(cfg['graphs']) / 'outputs' / entry['axes']
        if sha256_file(path) != entry['axes_sha256']:
            raise ValueError('Frozen direct-candidate axes changed')
        with np.load(path, allow_pickle=False) as f:
            np.testing.assert_array_equal(f['pilot_row'], np.flatnonzero(positive & (route == wanted)))
            np.testing.assert_array_equal(f['source'], source[f['pilot_row']])
            direct[f['pilot_row']] = f['direct_candidate']
    folds = sorted(set(data['fold']))
    if folds != cfg['expected_fold_groups']:
        raise ValueError('Conservative holdout fold axis changed')
    programs = sorted(data['state'])
    target = np.column_stack([data['state'][p] for p in programs])
    save(root / 'outputs/readout_definition_v1.json', {'programs': programs, 'genes': data['state_genes'],
        'native_RNA_unchanged': True, 'candidate_labels_and_routes_are_conditional_not_test_truth': True,
        'frozen_candidate_feature_universe_uses_measurement_availability_not_query_expression': True})
    all_profiles, index = [], []
    for number, heldout in enumerate(folds):
        train, query = fold_members(source, data['fold'], positive, heldout)
        shared = select_features(data['x'][train], source[train], cfg['n_features'])
        for wanted in cfg['routes']:
            tr, qu = train[route[train] == wanted], query[route[query] == wanted]
            if len(tr) < cfg['minimum_training_group'] or not len(qu):
                index.append({'fold': heldout, 'route': wanted, 'status': 'training_or_query_scope_unavailable',
                              'n_training_cells': len(tr), 'n_query_cells': len(qu)})
                continue
            folder = root / 'models' / ('fold_' + str(number)) / wanted
            folder.mkdir(parents=True)
            features = training_feature_sets(data['x'], source, route, train, wanted, cfg['n_features'], shared)
            anchors = balanced_cell_selection(source[tr], direct[tr], np.ones(len(tr), bool),
                cfg['anchor_cap_per_route'], stable_seed(cfg['seed'], heldout, wanted), rare_reserve=2)
            anchor_rows = tr[anchors]
            write_npz(folder / 'axes_and_readouts.npz', training_pilot_rows=tr, query_pilot_rows=qu,
                anchor_pilot_rows=anchor_rows, query_source=source[qu], query_direct_candidate=direct[qu],
                anchor_source=source[anchor_rows], anchor_direct_candidate=direct[anchor_rows],
                query_native_totals=data['native_total'][qu], query_RNA_readouts=target[qu],
                anchor_RNA_readouts=target[anchor_rows])
            coords, supported, models = {}, {}, []
            for variant in VARIANTS:
                began = time.monotonic(); selected = features[variant]
                genes = [data['genes'][j] for j in selected]
                reference, raw, corrected, pca = build_reference(data['x'][tr][:, selected], source[tr],
                    genes, cfg['seed'], cfg, symphony)
                before = reference_fingerprint(reference)
                raw_query, mapped, ok = map_query(reference, data['x'][qu][:, selected], genes, source[qu], symphony)
                after = reference_fingerprint(reference)
                if before != after:
                    raise ValueError('Heldout mapping changed reference')
                coords[variant] = (corrected[anchors].copy(), mapped.copy())
                supported[variant] = ok.copy()
                meta = {'fold': heldout, 'route': wanted, 'variant': variant,
                    'training_sources': sorted(set(source[tr])), 'query_sources': sorted(set(source[qu])),
                    'n_training_cells': len(tr), 'n_query_cells': len(qu), 'n_anchors': len(anchors),
                    'feature_selection_scope': 'all_training_lineages' if variant == VARIANTS[0] else 'training_route_only',
                    'selected_feature_indices': selected.tolist(), 'selected_feature_names': genes,
                    'reference_fingerprint_before': before, 'reference_fingerprint_after': after,
                    'harmony_converged': None if reference.uns['harmony']['converged'] is None else bool(reference.uns['harmony']['converged']),
                    'n_numerically_supported_queries': int(ok.sum()), 'wall_seconds': time.monotonic() - began,
                    'query_expression_excluded_from_feature_PCA_Harmony_fit': True,
                    'mapping_adapts_query_batch_reference_frozen': True}
                save_minimal_reference(folder / (variant + '_reference.npz'), reference, selected, meta)
                write_npz(folder / (variant + '_coordinates.npz'), anchor_embedding=coords[variant][0],
                    query_embedding=coords[variant][1], query_supported=ok)
                models.append(meta)
                checkpoint('lineage_holdout_reference_complete', fold=heldout, route=wanted, variant=variant)
                del reference, raw, corrected, pca, raw_query, mapped; gc.collect()
            predictions, scope_coverage = {}, []
            scopes = ['all_route_cells_not_pure_type'] + sorted(set(direct[qu]) - {'unresolved'})
            for scope in scopes:
                am = np.ones(len(anchors), bool) if scope == scopes[0] else direct[anchor_rows] == scope
                qm = np.ones(len(qu), bool) if scope == scopes[0] else direct[qu] == scope
                scope_coverage.append({'scope': scope, 'n_query_cells': int(qm.sum()), 'n_anchor_cells': int(am.sum()),
                    'status': 'evaluated' if am.sum() >= cfg['minimum_anchor_cells'] else 'insufficient_anchor_coverage'})
                predictions[scope] = {}
                for variant in VARIANTS:
                    pred = np.full((len(qu), len(programs)), np.nan)
                    use = qm & supported[variant]
                    if am.sum() >= cfg['minimum_anchor_cells'] and use.any():
                        zref, zquery = coords[variant]
                        pred[use], _ = inverse_neighbor_prediction(zref[am], target[anchor_rows[am]], zquery[use],
                                                                   cfg['state_neighbors'], cfg['threads'])
                    predictions[scope][variant] = pred
                write_npz(folder / (scope + '_predictions.npz'), **predictions[scope])
            profiles = paired_profiles(source[qu], direct[qu], target[qu], data['native_total'][qu],
                predictions, supported, programs, cfg['minimum_query_evaluation'])
            save(folder / 'models_and_coverage.json', {'models': models, 'scopes': scope_coverage})
            save(folder / 'paired_source_profiles.json', profiles)
            all_profiles.extend(dict(fold=heldout, route=wanted, **r) for r in profiles)
            files = {p.name: sha256_file(p) for p in sorted(folder.iterdir()) if p.is_file()}
            index.append({'fold': heldout, 'route': wanted, 'status': 'complete', 'relative_folder': str(folder.relative_to(root)),
                'n_training_cells': len(tr), 'n_query_cells': len(qu), 'n_anchor_cells': len(anchors), 'files_SHA256': files})
            checkpoint('lineage_holdout_route_complete', fold=heldout, route=wanted)
            written = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
            if written > cfg['output_cap_bytes']:
                raise ValueError('Lineage holdout output cap exceeded')
            del coords, predictions, supported; gc.collect()
        checkpoint('lineage_holdout_fold_complete', fold=heldout)
    save(root / 'outputs/lineage_source_holdout_index_v1.json', index)
    save(root / 'outputs/paired_source_holdout_profiles_v1.json', all_profiles)
    save(root / 'provenance/terminal_v1.json', {
        'status': 'NEW_ROUTE_TRAIN_ONLY_FEATURE_HOLDOUT_COMPLETE_INDEPENDENT_AUDIT_PENDING',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_folds': len(folds),
        'n_completed_route_folds': sum(r['status'] == 'complete' for r in index), 'n_source_profile_records': len(all_profiles),
        'index_sha256': sha256_file(root / 'outputs/lineage_source_holdout_index_v1.json'),
        'wall_seconds': time.monotonic() - started, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        'native_counts_copied': False, 'final_identity_accuracy_claimed': False, 'new_authoritative_labels': False,
        'candidate_routes_from_transductive_review_condition_evaluation': True,
        'known_or_possible_source_overlap_retained_in_same_fold': True, 'production_winner_selected': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--output', required=True)
    args = parser.parse_args()
    try:
        main(Path(args.output))
    except Exception as error:
        save(Path(args.output) / 'provenance/failure_v1.json', {
            'status': 'FAILED_NOT_ACCEPTED', 'error_class': type(error).__name__, 'error': str(error)})
        raise
