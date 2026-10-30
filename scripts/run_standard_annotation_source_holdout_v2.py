"""Study-held-out native RNA structure transfer, not identity accuracy.

HVGs/PCA/Harmony/prototypes use training folds only. Symphony query adaptation
keeps the reference frozen. Common within-route anchors control the probe's
coarse-label information and reference sampling across representations.
"""
import argparse
import collections
import datetime
import json
from pathlib import Path
import resource
import time
import sys,subprocess

import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.neighbors import NearestNeighbors

from annotation_pilot_sampling_v1 import balanced_cell_selection, stable_seed
from standard_annotation_order_primitives_v1 import select_features,pca_fit,coarse_routes_from_local
from annotation_order_standard_HVG_pilot_v2 import cluster
from full_context_atlas_v1 import sha256_file
from frozen_symphony_reference_standard_v2 import load_symphonypy, build_reference, map_query, save_minimal_reference, reference_fingerprint
from load_standard_annotation_pilot_native_v2 import load_pilot


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False);f.write('\n')


def inverse_neighbor_prediction(z_ref, targets, z_query, neighbors, threads=4):
    if (z_ref.ndim != 2 or z_query.ndim != 2 or targets.ndim != 2 or
            z_ref.shape[0] != targets.shape[0] or z_ref.shape[1] != z_query.shape[1]):
        raise ValueError('Reference/query/target axes differ')
    if not len(z_ref) or not len(z_query):return np.full((len(z_query), targets.shape[1]), np.nan), np.full(len(z_query), np.nan)
    number = min(int(neighbors), len(z_ref))
    search = NearestNeighbors(n_neighbors=number, n_jobs=threads).fit(z_ref)
    answer, distance = [], []
    for low in range(0, len(z_query), 2048):
        d, indices = search.kneighbors(z_query[low:low + 2048])
        weights = 1. / np.maximum(d, 1e-6)
        weights /= weights.sum(axis=1, keepdims=True)
        answer.append(np.sum(targets[indices] * weights[:, :, None], axis=1))
        distance.append(d[:, 0])
    return np.vstack(answer), np.concatenate(distance)


def conditional_scores(observed, predicted, totals, minimum):
    good = np.isfinite(observed) & np.isfinite(predicted) & np.isfinite(totals)
    if good.sum() < minimum:return {'status': 'insufficient_evaluable_query_cells', 'n_cells': int(good.sum())}
    y, pred, depth = observed[good], predicted[good], np.log1p(totals[good])
    design = np.column_stack([np.ones(len(depth)), depth - depth.mean()])
    residual_y = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    residual_p = pred - design @ np.linalg.lstsq(design, pred, rcond=None)[0]
    out = {'status': 'estimated', 'n_cells': len(y), 'MAE_native_logCP10K': float(np.mean(np.abs(y - pred))),
           'RMSE_native_logCP10K': float(np.sqrt(np.mean((y - pred) ** 2)))}
    if residual_y.std() < 1e-8 or residual_p.std() < 1e-8:
        out['depth_adjusted_structure_status'] = 'no_variable_observed_or_predicted_signal'
    else:
        out.update(depth_adjusted_structure_status='estimated',
                   depth_adjusted_Pearson=float(pearsonr(residual_y, residual_p).statistic),
                   depth_adjusted_Spearman=float(spearmanr(residual_y, residual_p).statistic))
    return out


def cluster_prototypes(z, partition, target, route, wanted):
    """Conditional diagnostic, not reannotation of mixed clusters.

    A/B clusters may cross coarse boundaries; subset by the same native coarse RNA route
    for this probe and report that conditioning. Do not hide cluster mixing.
    """
    centers, outcomes, metadata = [], [], []
    for code in np.unique(partition[route == wanted]):
        mask = (partition == code) & (route == wanted)
        if code < 0 or not mask.any():continue
        centers.append(z[mask].mean(axis=0));outcomes.append(target[mask].mean(axis=0))
        metadata.append({'partition_code': int(code), 'n_training_cells': int(mask.sum()),
                         'conditional_on_native_coarse_RNA_route': wanted})
    return np.asarray(centers), np.asarray(outcomes), metadata


def main(args):
    start = time.monotonic();root = Path(args.output)
    cfg = json.loads((root / 'config/source_holdout_spec_v1.json').read_text())
    for role, item in cfg['bindings'].items():
        if sha256_file(item['path']) != item['sha256']:raise ValueError('Frozen holdout input changed: ' + role)
    symphony = load_symphonypy(cfg['wheel_path'], cfg['wheel_SHA256'])
    def checkpoint(stage, **extra):
        save(root / 'checkpoints' / (str(time.time_ns()) + '.json'),
             {'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'stage': stage, 'wall_seconds': time.monotonic() - start, **extra})
    sys.path.insert(0,str(Path(cfg['joined'])/'scripts'))
    from extract_full_source_local_marker_evidence_v1 import storage_gate
    storage_gate(root,cfg)
    test=subprocess.run([sys.executable,'-B',str(root/'tests/test_standard_annotation_reference_v2.py')],capture_output=True,text=True,timeout=180)
    (root/'provenance/reference_test_stdout.txt').write_text(test.stdout);(root/'provenance/reference_test_stderr.txt').write_text(test.stderr)
    assert test.returncode==0 and 'Ran 2 tests' in test.stderr and 'skipped' not in test.stderr
    save(root/'provenance/reference_test_receipt_v1.json',dict(actual_exit_code=0,n_tests=2,no_skips=True))
    data = load_pilot(args.pilot, args.manifest, checkpoint, cfg['read_block_rows'],discovery=cfg['discovery'],native_auditor=cfg['native_auditor'])
    x, source, route, fine = data['x'], data['source'], data['route'], data['fine']
    eligible = data['native_total'] > 0
    programs = sorted(data['state'])
    target = np.column_stack([data['state'][name] for name in programs])
    save(root / 'config/common_measured_state_readouts_v1.json',
         {'genes': data['state_genes'], 'measurement_definitions_fixed_before_fits': True,
          'RNA_readouts_not_orthogonal_identity_truth': True})
    folds = sorted(set(data['fold']))
    if folds != cfg['expected_fold_groups']:raise ValueError('Conservative source fold axis changed')
    summaries, qualifications = [], []
    for fold_number, heldout in enumerate(folds):
        train = np.flatnonzero(eligible & (data['fold'] != heldout))
        query = np.flatnonzero(eligible & (data['fold'] == heldout))
        if set(source[train]) & set(source[query]):raise ValueError('Study leaked between training/query')
        selected,hvg_table = select_features(data['counts'][train],source[train],cfg['n_features'],data['genes'],cfg['candidate_HVGs'])
        genes = [data['genes'][j] for j in selected]
        xtrain, xquery = x[train][:, selected], x[query][:, selected]
        # This global model is one reference shared by A/B; no query expression
        # is used to fit genes, centering, PCA, Harmony or reference partitions.
        reference, raw, corrected, pca = build_reference(xtrain, source[train], genes, cfg['seed'], cfg, symphony)
        fingerprint = reference_fingerprint(reference)
        raw_query, mapped, supported = map_query(reference, xquery, genes, source[query], symphony)
        labels_B, n_B = cluster(corrected, cfg['seed'], cfg)
        labels_A = np.full(len(train), -1, dtype=np.int32);offset = 0
        for ds in np.unique(source[train]):
            member = np.flatnonzero(source[train] == ds)
            if len(member) < cfg['minimum_training_group']:continue
            local, _ = pca_fit(xtrain[member], cfg['seed'], cfg['components'], cfg['pca_iterations'])
            labels, number = cluster(local, cfg['seed'], cfg)
            labels_A[member] = labels + offset;offset += number
        train_route, train_route_records = coarse_routes_from_local(source[train],labels_A,data['native_flags'][train],data['source_assessed'])
        query_A = np.full(len(query),-1,np.int32);query_offset=0
        for ds in np.unique(source[query]):
            local=np.flatnonzero(source[query]==ds)
            if len(local)<cfg['minimum_training_group']:continue
            zq,_=pca_fit(xquery[local],cfg['seed'],cfg['components'],cfg['pca_iterations'])
            codes,nq=cluster(zq,cfg['seed'],cfg);query_A[local]=codes+query_offset;query_offset+=nq
        query_route,query_route_records=coarse_routes_from_local(source[query],query_A,data['native_flags'][query],data['source_assessed'])
        fold_dir = root / 'models' / ('fold_' + str(fold_number))
        fold_dir.mkdir(parents=True)
        hvg_table.to_csv(fold_dir/'training_only_standard_HVG_v2.tsv',sep='\t')
        save(fold_dir/'source_local_native_coarse_routes_v2.json',dict(training=train_route_records,query=query_route_records,query_source_local_clustering_uses_training_selected_genes=True,not_identity_truth=True))
        with (fold_dir/'private_fold_axes_and_partitions_v2.npz').open('xb') as f:np.savez_compressed(f,train_rows=train,query_rows=query,training_coarse_RNA_route=train_route,query_coarse_RNA_route=query_route,A_train=labels_A,B_train=labels_B,A_query=query_A)
        metadata = {'heldout_fold': heldout, 'training_sources': sorted(set(source[train])),
                    'query_sources': sorted(set(source[query])), 'fit_excludes_query_expression': True,
                    'native_manifest_SHA256': sha256_file(args.manifest),
                    'method': 'A_B_shared_reference', 'source_holdout_not_independent_patient_proof': True}
        save_minimal_reference(fold_dir / 'A_B_shared_minimal_reference.npz', reference, selected, metadata)
        save(fold_dir / 'reference_axes_and_parameters_v1.json',
             {**metadata, 'n_training_cells': len(train), 'n_query_cells': len(query),
              'reference_fingerprint_before_query': fingerprint,
              'reference_fingerprint_after_query': reference_fingerprint(reference),
              'converged': None if reference.uns['harmony']['converged'] is None else bool(reference.uns['harmony']['converged']),
              'selected_model_feature_indices': selected.tolist(), 'n_A_partitions': offset, 'n_B_partitions': n_B})
        checkpoint('train_only_A_B_reference_and_partitions_complete', fold=heldout)
        fold_coverage = {'fold': heldout, 'training_sources': sorted(set(source[train])),
                         'query_sources': sorted(set(source[query])), 'n_training_cells': len(train),
                         'n_query_cells': len(query), 'routes': [],
                         'query_to_reference_cell_ratio': len(query) / len(train),
                         'Symphony_not_assumed_equivalent_to_joint_reintegration': True,
                         'A_B_harmony_converged': None if reference.uns['harmony']['converged'] is None else bool(reference.uns['harmony']['converged'])}
        for wanted in cfg['reference_routes']:
            local_train = np.flatnonzero(train_route == wanted)
            local_query = np.flatnonzero(query_route == wanted)
            ntrain, nquery = len(local_train), len(local_query)
            if ntrain < cfg['minimum_training_group'] or nquery == 0:
                fold_coverage['routes'].append({'route': wanted, 'status': 'reference_or_query_scope_unavailable',
                                                'n_training_cells': ntrain, 'n_query_cells': nquery})
                continue
            c_ref, c_raw, c_corrected, c_pca = build_reference(xtrain[local_train], source[train[local_train]],
                                                            genes, cfg['seed'], cfg, symphony)
            c_before = reference_fingerprint(c_ref)
            c_raw_query, c_mapped, c_supported = map_query(c_ref, xquery[local_query], genes, source[query[local_query]], symphony)
            labels_C, n_C = cluster(c_corrected, cfg['seed'], cfg)
            save_minimal_reference(fold_dir / (wanted + '_minimal_reference.npz'), c_ref, selected,
                                   {**metadata, 'method': 'C_within_source_local_native_coarse_RNA', 'route': wanted})
            # The same training cells anchor the within-route A/B and C probes.
            anchors = balanced_cell_selection(source[train[local_train]], fine[train[local_train]],
                np.ones(ntrain, bool), cfg['anchor_cap_per_route'], stable_seed(cfg['seed'], heldout, wanted), rare_reserve=2)
            anchor_target = target[train[local_train[anchors]]]
            truth = target[query[local_query]]
            common = supported[local_query] & c_supported
            predictions = {}
            for name, zref, zquery in [
                ('A_B_shared_reference_neighbors', corrected[local_train[anchors]], mapped[local_query]),
                ('C_reference_neighbors', c_corrected[anchors], c_mapped)]:
                pred = np.full(truth.shape, np.nan)
                pred[common], _ = inverse_neighbor_prediction(zref, anchor_target, zquery[common], cfg['state_neighbors'], cfg['threads'])
                predictions[name] = pred
            for name, labels in [('A_conditional_partition_prototypes', labels_A), ('B_conditional_partition_prototypes', labels_B)]:
                centers, outcome, proto = cluster_prototypes(corrected, labels, target[train], train_route, wanted)
                pred = np.full(truth.shape, np.nan)
                pred[common], _ = inverse_neighbor_prediction(centers, outcome, mapped[local_query][common], cfg['prototype_neighbors'], cfg['threads'])
                predictions[name] = pred
                with (fold_dir / (name + '_' + wanted + '.npz')).open('xb') as f:
                    np.savez_compressed(f, centers=centers, native_RNA_program_means=outcome)
                save(fold_dir / (name + '_' + wanted + '_metadata.json'), proto)
            centers, outcome, proto = cluster_prototypes(c_corrected, labels_C, target[train[local_train]],
                                                        np.repeat(wanted, ntrain), wanted)
            pred = np.full(truth.shape, np.nan)
            pred[common], distance = inverse_neighbor_prediction(centers, outcome, c_mapped[common], cfg['prototype_neighbors'], cfg['threads'])
            predictions['C_conditional_partition_prototypes'] = pred
            with (fold_dir / ('C_partition_prototypes_' + wanted + '.npz')).open('xb') as f:
                np.savez_compressed(f, centers=centers, native_RNA_program_means=outcome)
            save(fold_dir / ('C_partition_prototypes_' + wanted + '_metadata.json'), proto)
            with (fold_dir / (wanted + '_training_anchor_positions.npz')).open('xb') as f:
                np.savez_compressed(f, frozen_pilot_positions=train[local_train[anchors]])
            for ds in np.unique(source[query[local_query]]):
                member = (source[query[local_query]] == ds) & common
                for name, pred in predictions.items():
                    for j, program in enumerate(programs):
                        summaries.append({'heldout_fold': heldout, 'query_source': str(ds), 'source_local_native_coarse_route': wanted,
                                          'probe': name, 'program': program,
                                          'n_query_route_cells': int((source[query[local_query]] == ds).sum()),
                                          'same_evaluable_query_subset_across_compared_probes': True,
                                          **conditional_scores(truth[member, j], pred[member, j],
                                                               data['native_total'][query[local_query[member]]], cfg['minimum_query_evaluation'])})
            with (fold_dir/(wanted+'_private_scoring_readback_v2.npz')).open('xb') as f:
                np.savez_compressed(f,query_pilot_rows=query[local_query],anchor_pilot_rows=train[local_train[anchors]],query_source=source[query[local_query]],query_native_total=data['native_total'][query[local_query]],query_RNA=truth,anchor_RNA=anchor_target,common_supported=common,B_anchor=corrected[local_train[anchors]],B_query=mapped[local_query],C_anchor=c_corrected[anchors],C_query=c_mapped,**predictions)
            fold_coverage['routes'].append({'route': wanted, 'status': 'train_only_reference_mapped_not_identity_accuracy',
                'n_training_cells': ntrain, 'n_query_cells': nquery, 'n_common_numerically_supported_query_cells': int(common.sum()),
                'n_training_anchor_cells_identical_A_B_C': len(anchors), 'n_C_partitions': n_C,
                'C_query_to_reference_cell_ratio': nquery / ntrain,
                'C_harmony_converged': None if c_ref.uns['harmony']['converged'] is None else bool(c_ref.uns['harmony']['converged']),
                'C_reference_unchanged': c_before == reference_fingerprint(c_ref),
                'query_source_local_coarse_RNA_is_conditioning_input_not_validation_truth': True})
            checkpoint('source_holdout_route_complete', fold=heldout, route=wanted)
        outside = ~np.isin(query_route, cfg['reference_routes'])
        fold_coverage['outside_initial_lineage_reference_scope'] = dict(collections.Counter(query_route[outside]))
        fold_coverage['outside_scope_not_forced_to_fine_labels'] = True
        qualifications.append(fold_coverage)
        save(root / 'outputs' / ('fold_' + str(fold_number) + '_qualification.json'), fold_coverage)
        checkpoint('source_holdout_fold_complete', fold=heldout)
        storage_gate(root,cfg)
        size = sum(p.stat().st_size for p in root.rglob('*') if p.is_file())
        if size > cfg['output_cap_bytes']:raise ValueError('Bounded reference-model output cap exceeded')
    save(root / 'outputs/source_heldout_RNA_structure_v1.json', summaries)
    save(root / 'outputs/source_heldout_qualification_v1.json', qualifications)
    save(root / 'provenance/terminal_v1.json',
        {'status': 'STANDARD_TRAIN_ONLY_NATIVE_COARSE_STUDY_HELDOUT_RNA_STRUCTURE_COMPLETE_REQUIRES_INDEPENDENT_AUDIT',
         'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'n_folds': len(folds),
         'query_expression_excluded_from_HVG_PCA_Harmony_reference_partition_fits': True,
         'query_native_RNA_readouts_used_only_for_scoring': True,
         'known_or_possible_overlapping_sources_kept_in_same_fold': True,
         'all_sources_not_claimed_independent_patients': True, 'C_native_coarse_RNA_routes_conditional_not_new_truth': True,
         'A_B_shared_neighbor_probe_one_observation_only': True, 'prototype_complexity_not_controlled_in_this_primary_probe': True,
         'method_winner_or_identity_accuracy_accepted': False,
         'wall_seconds': time.monotonic() - start, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
         'native_expression_copied_to_disk': False})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['output', 'pilot', 'manifest']:p.add_argument('--' + name, required=True)
    a = p.parse_args()
    try:main(a)
    except Exception as error:
        save(Path(a.output) / 'provenance/failure_v1.json',
             {'status': 'FAILED_NOT_ACCEPTED', 'error_class': type(error).__name__, 'error': str(error)})
        raise
