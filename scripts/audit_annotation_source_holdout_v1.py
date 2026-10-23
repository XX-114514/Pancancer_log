"""Audit heldout axes and replay native prototype scores on a size-selected fold.

All-fold feature selection is recomputed from training counts. Numeric query
replay uses whole query batches and saved immutable models, not new refits.
Reference-cell neighbor scores are checked for axes/bounds, not numerically
reproduced by this prototype replay. Individual data/models remain private.
"""
import argparse
import collections
import datetime
import json
from pathlib import Path
import resource
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial.distance import cdist
from scipy.stats import rankdata

from annotation_pilot_sampling_v1 import balanced_cell_selection, stable_seed
from annotation_representation_pilot_v1 import select_features
from full_context_atlas_v1 import FullContextAtlas, sha256_file
from frozen_symphony_reference_v1 import load_symphonypy, map_query, reference_fingerprint
from load_annotation_pilot_native_v1 import load_pilot


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:json.dump(value, stream, indent=2, allow_nan=False);stream.write('\n')


def restore_reference(path, expected_genes, expected_indices):
    import anndata as ad
    with np.load(path, allow_pickle=False) as f:
        values = {k: f[k].copy() for k in f.files}
    if (values['genes'].tolist() != list(expected_genes) or
            not np.array_equal(values['feature_candidate_indices'], expected_indices)):
        raise ValueError('Frozen reference feature axis changed')
    genes, PCs, centers, nr = values['genes'], values['loadings'], values['C'], values['Nr']
    # symphonypy compresses centers as soft clusters × PCA coordinates.
    # This is not Harmony's alternate transposed representation.
    if PCs.ndim != 2 or PCs.shape[0] != len(genes) or centers.shape != (len(nr), PCs.shape[1]):
        raise ValueError('Reference PCA/Harmony axes differ')
    if values['sigma'].shape != nr.shape or (nr < 0).any() or not np.all(values['std'] == 1):
        raise ValueError('Invalid compressed reference parameters')
    for k in ['mean', 'std', 'loadings', 'C', 'Nr', 'sigma']:
        if not np.isfinite(values[k]).all():raise ValueError('Nonfinite reference parameter')
    if not np.allclose(PCs.T @ PCs, np.eye(PCs.shape[1]), atol=2e-4, rtol=2e-4):
        raise ValueError('Reference PCA loadings not orthonormal')
    reference = ad.AnnData(sp.csr_matrix((0, len(genes)), dtype=np.float32),
                          var=pd.DataFrame({'mean': values['mean'], 'std': values['std']}, index=genes))
    reference.varm['PCs'] = PCs
    reference.uns['harmony'] = {'C': centers, 'Nr': nr, 'K': len(nr), 'sigma': values['sigma'],
                                'ref_basis_loadings': 'PCs', 'ref_basis_adjusted': 'X_pca_harmony', 'converged': None}
    return reference, json.loads(str(values['metadata_json']))


def prototype_prediction(centers, outcomes, query, neighbors):
    """Independent Euclidean distances and inverse-distance reconstruction."""
    if len(centers) != len(outcomes) or centers.shape[1] != query.shape[1] or not len(centers):
        raise ValueError('Prototype/target/query axes differ')
    chunks = []
    for low in range(0, len(query), 1024):
        distances = cdist(query[low:low + 1024], centers, metric='euclidean')
        indices = np.argsort(distances, axis=1, kind='stable')[:, :min(neighbors, len(centers))]
        near = np.take_along_axis(distances, indices, axis=1)
        weights = 1 / np.maximum(near, 1e-6)
        weights /= weights.sum(axis=1, keepdims=True)
        chunks.append((outcomes[indices] * weights[..., None]).sum(axis=1))
    return np.vstack(chunks) if chunks else np.empty((0, outcomes.shape[1]))


def independent_scores(observed, predicted, totals, minimum):
    keep = np.isfinite(observed) & np.isfinite(predicted) & np.isfinite(totals)
    n = int(keep.sum())
    if n < minimum:return {'status': 'insufficient_evaluable_query_cells', 'n_cells': n}
    y, pred, depth = observed[keep].astype(float), predicted[keep].astype(float), np.log1p(totals[keep])
    depth -= depth.mean()
    def residual(v):
        centered = v - v.mean()
        norm = depth @ depth
        return centered - depth * (depth @ centered / norm) if norm > 0 else centered
    ry, rp = residual(y), residual(pred)
    out = {'status': 'estimated', 'n_cells': n,
           'MAE_native_logCP10K': float(np.abs(y - pred).mean()),
           'RMSE_native_logCP10K': float(np.sqrt(np.square(y - pred).mean()))}
    if ry.std() < 1e-8 or rp.std() < 1e-8:
        out['depth_adjusted_structure_status'] = 'no_variable_observed_or_predicted_signal'
    else:
        out.update(depth_adjusted_structure_status='estimated',
                   depth_adjusted_Pearson=float(np.corrcoef(ry, rp)[0, 1]),
                   depth_adjusted_Spearman=float(np.corrcoef(rankdata(ry), rankdata(rp))[0, 1]))
    return out


def check_score_bounds(row):
    if row['n_cells'] < 0 or row['n_cells'] > row['n_query_route_cells']:
        raise ValueError('Query score denominator outside scope')
    if row['status'] not in {'estimated', 'insufficient_evaluable_query_cells'}:
        raise ValueError('Unexpected score status')
    for field in ['depth_adjusted_Pearson', 'depth_adjusted_Spearman']:
        if field in row and (not np.isfinite(row[field]) or not -1 <= row[field] <= 1):
            raise ValueError('Invalid structure correlation')
    for field in ['MAE_native_logCP10K', 'RMSE_native_logCP10K']:
        if field in row and (not np.isfinite(row[field]) or row[field] < 0):
            raise ValueError('Invalid reconstruction error')


def main(args):
    start = time.monotonic();run, out, pilot = Path(args.run), Path(args.output), Path(args.pilot)
    cfg = json.loads((run / 'config/source_holdout_spec_v1.json').read_text())
    launch = json.loads((run / 'provenance/launch_receipt_v1.json').read_text())
    supervisor = json.loads((run / 'provenance/supervisor_terminal_v1.json').read_text())
    if supervisor['actual_child_exit_code'] != 0 or supervisor['wall_timeout']:raise ValueError('Heldout producer not successful')
    if sha256_file(run / 'config/source_holdout_spec_v1.json') != launch['config_SHA256']:raise ValueError('Model configuration changed')
    for file, digest in launch['script_SHA256'].items():
        if sha256_file(run / file) != digest:raise ValueError('Executed model code changed')
    for role, binding in cfg['bindings'].items():
        if sha256_file(binding['path']) != binding['sha256']:raise ValueError('Bound input changed: ' + role)
    total_bytes = sum(p.stat().st_size for p in run.rglob('*') if p.is_file())
    if total_bytes > cfg['output_cap_bytes']:raise ValueError('Final output cap exceeded including final summary')
    symphony = load_symphonypy(cfg['wheel_path'], cfg['wheel_SHA256'])
    def checkpoint(stage, **extra):
        save(out / 'checkpoints' / (str(time.time_ns()) + '.json'),
             {'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'stage': stage, **extra})
    data = load_pilot(pilot, args.manifest, checkpoint, cfg['read_block_rows'])
    reader = FullContextAtlas(args.manifest)
    x, sources, routes, eligible = data['x'], data['source'], data['route'], data['native_total'] > 0
    programs = sorted(data['state'])
    definitions = json.loads((run / 'config/common_measured_state_readouts_v1.json').read_text())
    if definitions['genes'] != data['state_genes']:raise ValueError('Common RNA readout definitions changed')
    folds = sorted(set(data['fold']))
    if folds != cfg['expected_fold_groups']:raise ValueError('Source fold grouping changed')
    summaries = json.loads((run / 'outputs/source_heldout_RNA_structure_v1.json').read_text())
    coverage = json.loads((run / 'outputs/source_heldout_qualification_v1.json').read_text())
    if len(coverage) != len(folds):raise ValueError('Missing heldout fold coverage')
    keys = ['heldout_fold', 'query_source', 'inherited_route', 'probe', 'program']
    summary_map = {tuple(row[k] for k in keys): row for row in summaries}
    if len(summary_map) != len(summaries):raise ValueError('Repeated heldout scientific result key')
    for row in summaries:check_score_bounds(row)
    expected_keys = set();fold_checks = []
    # Size-selected before inspecting any effect/diagnostic score.
    replay_fold = min(folds, key=lambda f: int((eligible & (data['fold'] == f)).sum()))
    replay = None
    for number, heldout in enumerate(folds):
        train = np.flatnonzero(eligible & (data['fold'] != heldout))
        query = np.flatnonzero(eligible & (data['fold'] == heldout))
        train_sources, query_sources = sorted(set(sources[train])), sorted(set(sources[query]))
        if set(train_sources) & set(query_sources):raise ValueError('Study overlap across train/query')
        c = coverage[number]
        if c != json.loads((run / 'outputs' / ('fold_' + str(number) + '_qualification.json')).read_text()):
            raise ValueError('Per-fold and combined qualifications disagree')
        if (c['fold'] != heldout or c['training_sources'] != train_sources or c['query_sources'] != query_sources or
                c['n_training_cells'] != len(train) or c['n_query_cells'] != len(query)):
            raise ValueError('Recorded fold cells or source axes disagree')
        selected = select_features(x[train], sources[train], cfg['n_features'])
        genes = [data['genes'][j] for j in selected]
        # Confirm measurement availability and physical disjointness, not names alone.
        for ds in sorted(set(sources)):
            lookup = reader.feature_lookup(ds)
            if any(len(lookup.get(g, set())) != 1 for g in genes):raise ValueError('Model feature not uniquely measured')
            model_cols = {next(iter(lookup[g])) for g in genes}
            heldout_cols = {col for g in data['readout_genes'] for col in lookup.get(g, set())}
            if model_cols & heldout_cols:raise ValueError('RNA validation gene alias used for training')
        folder = run / 'models' / ('fold_' + str(number))
        axes = json.loads((folder / 'reference_axes_and_parameters_v1.json').read_text())
        if (axes['training_sources'] != train_sources or axes['query_sources'] != query_sources or
                axes['selected_model_feature_indices'] != selected.tolist()):
            raise ValueError('Training-only selected features or sources differ')
        ref, metadata = restore_reference(folder / 'A_B_shared_minimal_reference.npz', genes, selected)
        expected_fp = reference_fingerprint(ref)
        if (expected_fp != axes['reference_fingerprint_before_query'] or
                expected_fp != axes['reference_fingerprint_after_query']):raise ValueError('Frozen reference changed')
        if metadata['heldout_fold'] != heldout or metadata['training_sources'] != train_sources or metadata['query_sources'] != query_sources:
            raise ValueError('Minimal reference metadata differs')
        expected_mean = np.asarray(x[train][:, selected].mean(axis=0)).ravel()
        if not np.allclose(ref.var['mean'], expected_mean, atol=3e-5, rtol=3e-5):
            raise ValueError('Reference centering differs from training-only expression')
        if not np.isclose(np.sum(ref.uns['harmony']['Nr']), len(train), rtol=1e-5):raise ValueError('Compressed reference count differs')
        if heldout == replay_fold:replay = (number, train, query, selected, genes, ref)
        route_checks = []
        for wanted in cfg['reference_routes']:
            rt = train[routes[train] == wanted];rq = query[routes[query] == wanted]
            entries = [q for q in c['routes'] if q['route'] == wanted]
            if len(entries) != 1:raise ValueError('Missing/repeated route qualification')
            cr = entries[0]
            if cr['n_training_cells'] != len(rt) or cr['n_query_cells'] != len(rq):raise ValueError('Route cell axis differs')
            if cr['status'] == 'reference_or_query_scope_unavailable':continue
            cref, cm = restore_reference(folder / (wanted + '_minimal_reference.npz'), genes, selected)
            if cm['route'] != wanted or cm['heldout_fold'] != heldout:raise ValueError('Route model identity differs')
            if not np.allclose(cref.var['mean'], np.asarray(x[rt][:, selected].mean(axis=0)).ravel(), atol=3e-5, rtol=3e-5):
                raise ValueError('Route centering differs from training-only cells')
            if not np.isclose(np.sum(cref.uns['harmony']['Nr']), len(rt), rtol=1e-5):raise ValueError('Route reference mass differs')
            with np.load(folder / (wanted + '_training_anchor_positions.npz'), allow_pickle=False) as f:anchors = f['frozen_pilot_positions']
            expected = rt[balanced_cell_selection(sources[rt], data['fine'][rt], np.ones(len(rt), bool), cfg['anchor_cap_per_route'],
                          stable_seed(cfg['seed'], heldout, wanted), rare_reserve=2)]
            if not np.array_equal(anchors, expected):raise ValueError('Anchors include query/wrong route or changed sampling')
            probes = ['A_B_shared_reference_neighbors', 'C_reference_neighbors']
            for method, prefix in [('A', 'A_conditional_partition_prototypes'), ('B', 'B_conditional_partition_prototypes'),
                                   ('C', 'C_partition_prototypes')]:
                meta = json.loads((folder / (prefix + '_' + wanted + '_metadata.json')).read_text())
                with np.load(folder / (prefix + '_' + wanted + '.npz'), allow_pickle=False) as f:
                    centers, targets = f['centers'], f['native_RNA_program_means']
                if (centers.shape != (len(meta), ref.varm['PCs'].shape[1]) or targets.shape != (len(meta), len(programs)) or
                        sum(t['n_training_cells'] for t in meta) != len(rt) or
                        any(t['conditional_on_inherited_route'] != wanted for t in meta) or
                        len({t['partition_code'] for t in meta}) != len(meta)):
                    raise ValueError('Prototype count, dimension or route differs')
                if not np.isfinite(centers).all() or not np.isfinite(targets).all():raise ValueError('Nonfinite prototype')
                probes.append(method + '_conditional_partition_prototypes')
            for ds in sorted(set(sources[rq])):
                n = int((sources[rq] == ds).sum())
                denominators = set()
                for probe in probes:
                    for program in programs:
                        key = (heldout, ds, wanted, probe, program);expected_keys.add(key)
                        row = summary_map[key]
                        if row['n_query_route_cells'] != n:raise ValueError('Query source/route denominator differs')
                        denominators.add(row['n_cells'])
                if len(denominators) != 1:raise ValueError('Compared probes have different evaluable query denominators')
            route_checks.append({'route': wanted, 'axes_and_training_anchors_pass': True,
                                 'producer_converged': cr['C_harmony_converged']})
        outside = dict(collections.Counter(routes[query][~np.isin(routes[query], cfg['reference_routes'])]))
        if c['outside_initial_lineage_reference_scope'] != outside:raise ValueError('Unmodeled routes omitted')
        fold_checks.append({'fold': heldout, 'training_only_features_recomputed': True,
                            'reference_centers_source_axes_and_anchors_pass': True,
                            'producer_AB_converged': c['A_B_harmony_converged'], 'routes': route_checks})
        checkpoint('heldout_fold_axes_and_train_only_features_audited', fold=heldout)
    if expected_keys != set(summary_map):raise ValueError('Missing or extra scientific result records')
    # Replay all three route/prototype methods in the smallest complete query
    # fold. Mapping uses exactly its original whole query batch composition.
    number, train, query, selected, genes, ref = replay
    folder = run / 'models' / ('fold_' + str(number))
    _, mapped, supported = map_query(ref, x[query][:, selected], genes, sources[query], symphony)
    replay_count, max_difference = 0, 0.
    for wanted in cfg['reference_routes']:
        qlocal = np.flatnonzero(routes[query] == wanted)
        if not len(qlocal):continue
        cref, _ = restore_reference(folder / (wanted + '_minimal_reference.npz'), genes, selected)
        _, c_mapped, c_supported = map_query(cref, x[query[qlocal]][:, selected], genes, sources[query[qlocal]], symphony)
        common = supported[qlocal] & c_supported
        for method, prefix, coordinates in [
                ('A', 'A_conditional_partition_prototypes', mapped[qlocal]),
                ('B', 'B_conditional_partition_prototypes', mapped[qlocal]),
                ('C', 'C_partition_prototypes', c_mapped)]:
            with np.load(folder / (prefix + '_' + wanted + '.npz'), allow_pickle=False) as f:
                prediction = prototype_prediction(f['centers'], f['native_RNA_program_means'], coordinates[common], cfg['prototype_neighbors'])
            chosen = query[qlocal[common]]
            for ds in sorted(set(sources[chosen])):
                keep = sources[chosen] == ds
                for j, program in enumerate(programs):
                    score = independent_scores(data['state'][program][chosen[keep]], prediction[keep, j], data['native_total'][chosen[keep]], cfg['minimum_query_evaluation'])
                    original = summary_map[(replay_fold, ds, wanted, method + '_conditional_partition_prototypes', program)]
                    for key, value in score.items():
                        if isinstance(value, float):
                            error = abs(value - original[key]);max_difference = max(max_difference, error)
                            if not np.isclose(value, original[key], atol=2e-5, rtol=2e-5):raise ValueError('Prototype native RNA numeric replay differs')
                        elif original[key] != value:raise ValueError('Prototype numeric replay status/denominator differs')
                    replay_count += 1
        checkpoint('smallest_fold_whole_query_prototype_replay_complete', route=wanted)
    output = out / 'outputs/private_source_holdout_independent_audit_v1.json'
    save(output, {'fold_checks': fold_checks, 'n_scientific_result_keys': len(summaries),
                  'size_selected_numeric_replay_fold': replay_fold, 'prototype_replay_records': replay_count,
                  'largest_replayed_numeric_error': max_difference,
                  'training_cell_neighbor_scores_numerically_replayed': False})
    save(out / 'provenance/terminal_v1.json', {
        'status': 'HELDOUT_AXES_TRAINING_FEATURES_AND_REPRESENTATIVE_PROTOTYPE_REPLAY_PASS_NOT_IDENTITY',
        'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_folds_audited': len(folds),
        'all_fold_source_gene_model_anchor_and_scoring_axes_pass': True,
        'all_fold_training_only_feature_selection_recomputed': True,
        'size_selected_whole_query_numeric_replay_fold': replay_fold, 'prototype_replay_records': replay_count,
        'largest_replayed_numeric_error': max_difference, 'native_query_loader_shared_with_qualified_producer': True,
        'training_cell_neighbor_scores_numerically_replayed': False, 'identity_or_method_winner_accepted': False,
        'production_final_total_bytes': total_bytes, 'production_cap_bytes': cfg['output_cap_bytes'],
        'output_sha256': sha256_file(output), 'wall_seconds': time.monotonic() - start,
        'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024})


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    for name in ['run', 'pilot', 'manifest', 'output']:parser.add_argument('--' + name, required=True)
    main(parser.parse_args())
