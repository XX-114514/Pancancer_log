"""Controlled discovery partitions; these are not accepted biological labels.

Native expression stays in RAM. Same measured feature candidates are used for
all routes. Discovery uses all frozen sources; source-held-out reference mapping
is a separate required stage and is never claimed by this program.
"""
import argparse
import collections
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata, pearsonr, spearmanr
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from full_context_atlas_v1 import FullContextAtlas, sha256_file
from native_panel_aggregate_v1 import gather_csr_rows, LINEAGE_PANELS


ROUTES = {
    'T cell': 'T_NK', 'T/NK lymphoid': 'T_NK', 'T': 'T_NK', 'NK': 'T_NK',
    'T_NK_boundary': 'T_NK', 'B cell': 'B', 'B': 'B',
    'Plasma cell': 'Plasma', 'Plasma': 'Plasma',
    'Myeloid': 'Myeloid', 'Monocyte': 'Myeloid', 'Macrophage': 'Myeloid',
    'DC': 'Myeloid', 'Dendritic cell': 'Myeloid', 'Mast cell': 'Mast', 'Mast': 'Mast',
    'Neutrophil': 'Neutrophil', 'Erythroid': 'Erythroid',
    'Epithelial cell': 'Epithelial', 'Epithelial': 'Epithelial',
    'Endothelial cell': 'Endothelial', 'Endothelial': 'Endothelial',
    'Fibroblast': 'Fibroblast', 'Fibroblast cell': 'Fibroblast',
    'Pericyte': 'Mural', 'Smooth muscle cell': 'Mural', 'Mural': 'Mural',
    'Mural cell': 'Mural', 'Pericyte_SMC': 'Mural',
    'Oligodendrocyte': 'Glial',
}


def save(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write('\n')


def normalized_log_counts(counts, native_totals):
    x = counts.tocsr().astype(np.float32)
    totals = np.asarray(native_totals, dtype=np.float64)
    if (x.shape[0] != len(totals) or not np.isfinite(totals).all() or
            (totals < 0).any() or not np.isfinite(x.data).all() or (x.data < 0).any() or
            not np.equal(x.data, np.floor(x.data)).all()):
        raise ValueError('Invalid native count domain or denominator')
    if np.any(np.asarray(x.sum(axis=1, dtype=np.float64)).ravel() > totals + 1e-5):
        raise ValueError('Feature sum exceeds native RNA denominator')
    positive = totals > 0
    scale = np.zeros(len(totals), dtype=np.float32)
    scale[positive] = 10000.0 / totals[positive]
    x = x.multiply(scale[:, None]).tocsr()
    x.data = np.log1p(x.data)
    return x, positive


def feature_columns(lookup, model, readout):
    columns = [next(iter(lookup[g])) if len(lookup.get(g, set())) == 1 else None for g in model]
    if any(c is None for c in columns) or len(set(columns)) != len(columns):
        raise ValueError('Model genes are missing, ambiguous or physical aliases')
    held = {next(iter(lookup[g])) for g in readout if len(lookup.get(g, set())) == 1}
    if set(columns) & held:
        raise ValueError('A held-out RNA readout aliases a physical model feature')
    return np.asarray(columns, dtype=np.int64)


def select_features(training, sources, number):
    """Mean within-source variance/mean rank; input must be training rows only.

    This is a declared sparse log-expression dispersion heuristic, not a claim
    to implement Seurat's HVG estimator. No source/query/clinical labels used as
    identity truth. All-source discovery and held-out fitting must call apart.
    """
    if training.shape[0] != len(sources):
        raise ValueError('Training/source axis differs')
    ranks = []
    for source in np.unique(sources):
        x = training[np.asarray(sources) == source]
        mean = np.asarray(x.mean(axis=0)).ravel()
        second = np.asarray(x.multiply(x).mean(axis=0)).ravel()
        variance = np.maximum(second - mean * mean, 0)
        score = np.divide(variance, mean + 1e-3)
        score[np.asarray((x > 0).sum(axis=0)).ravel() < 5] = -np.inf
        ranks.append(rankdata(-score, method='average') / training.shape[1])
    average = np.mean(ranks, axis=0)
    return np.lexsort((np.arange(len(average)), average))[:min(number, len(average))]


def pca_fit(x, seed, components, iterations):
    dense = x.toarray().astype(np.float32, copy=False)
    n = min(components, dense.shape[0] - 1, dense.shape[1] - 1)
    if n < 2:
        raise ValueError('Too few observations/features for PCA')
    model = PCA(n_components=n, svd_solver='randomized', random_state=seed,
                iterated_power=iterations, copy=False)
    embedding = model.fit_transform(dense).astype(np.float32, copy=False)
    if not np.isfinite(embedding).all():
        raise ValueError('Non-finite PCA')
    return embedding, model


def harmonize(embedding, sources, seed, parameters):
    if len(np.unique(sources)) == 1:
        return embedding.copy(), {'status': 'one_source_no_correction'}
    import harmonypy
    fit = harmonypy.run_harmony(embedding, pd.DataFrame({'source': sources}), ['source'],
                               theta=parameters['theta'], lamb=parameters['lambda'],
                               max_iter_harmony=parameters['max_iter_harmony'],
                               ncores=parameters['threads'], random_state=seed, verbose=False)
    corrected = np.asarray(fit.Z_corr)
    if corrected.shape == embedding.T.shape:
        corrected = corrected.T
    if corrected.shape != embedding.shape or not np.isfinite(corrected).all():
        raise ValueError('Unexpected Harmony orientation or non-finite output')
    return corrected.astype(np.float32), {'status': 'discovery_all_sources_jointly_fit',
                                        'implementation_version': harmonypy.__version__}


def cluster(embedding, seed, parameters):
    import anndata as ad
    import scanpy as sc
    import numba
    sc.settings.n_jobs = parameters['threads']
    numba.set_num_threads(parameters['threads'])
    a = ad.AnnData(sp.csr_matrix((len(embedding), 0)))
    a.obsm['X_pilot'] = embedding
    sc.pp.neighbors(a, use_rep='X_pilot', n_neighbors=min(parameters['neighbors'], len(embedding) - 1),
                    random_state=seed, transformer='pynndescent', metric='euclidean')
    sc.tl.leiden(a, resolution=parameters['resolution'], random_state=seed,
                 flavor='igraph', n_iterations=2, directed=False)
    labels = a.obs['leiden'].to_numpy().astype(str)
    unique, codes = np.unique(labels, return_inverse=True)
    return codes.astype(np.int32), len(unique)


def conditional_readout_metrics(embedding, states, native_totals, candidate, indices, neighbors):
    """State recovery by same-source/same-prior-route neighbors, depth adjusted.

    Fixed subsamples are identical for every representation. Readouts are not
    lineage ground truth, causal function or an external validation endpoint.
    """
    indices = np.asarray(indices, dtype=np.int64)
    if len(indices) < max(30, neighbors + 2):
        return {'status': 'insufficient_conditional_cells'}
    z = embedding[indices]
    _, near = NearestNeighbors(n_neighbors=neighbors + 1, n_jobs=4).fit(z).kneighbors(z)
    # Tied distances may omit the query itself; remove exact row ID rather than
    # blindly dropping column 0 and accidentally retaining a self-neighbor.
    near = np.asarray([row[row != i][:neighbors] for i, row in enumerate(near)])
    depth = np.log1p(np.asarray(native_totals)[indices])
    design = np.column_stack([np.ones(len(depth)), depth - depth.mean()])
    out = {'status': 'conditional_RNA_diagnostic_not_independent_truth', 'n_cells': len(indices), 'states': {}}
    for name, values in states.items():
        y = np.asarray(values)[indices]
        if not np.isfinite(y).all():
            out['states'][name] = {'status': 'measurement_incomplete'}
            continue
        residual = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
        pred = residual[near].mean(axis=1)
        if residual.std() < 1e-8 or pred.std() < 1e-8:
            out['states'][name] = {'status': 'no_variable_signal'}
        else:
            out['states'][name] = {'status': 'estimated',
                                   'depth_adjusted_neighbor_Pearson': float(pearsonr(residual, pred).statistic),
                                   'depth_adjusted_neighbor_Spearman': float(spearmanr(residual, pred).statistic)}
    rare = np.asarray(candidate)[indices]
    n = int(rare.sum())
    out['gamma_delta_RNA_candidate'] = {'n_candidates_not_verified_identity': n,
        'conditional_neighbor_candidate_fraction': float(rare[near[rare]].mean()) if n else None,
        'conditional_candidate_prevalence': float(rare.mean()),
        'definition': 'detected_TRDC_and_TRGC1_or_TRGC2_and_CD3D_not_TCR_validation'}
    return out


def marker_summary(labels, panel_log, detected, genes, measured):
    result = []
    for code in np.unique(labels):
        member = labels == code
        values = {}
        for name, panel in LINEAGE_PANELS.items():
            columns = [genes.index(g) for g in panel if g in genes and measured[genes.index(g)]]
            values[name] = {'n_measured': len(columns), 'n_declared': len(panel),
                            'mean_native_logCP10K': float(panel_log[member][:, columns].mean()) if columns else None,
                            'mean_detection_fraction': float(detected[member][:, columns].mean()) if columns else None}
        result.append({'partition_code': int(code), 'n_cells': int(member.sum()),
                       'lineage_panel_evidence': values, 'biological_identity_assigned': False})
    return result


def main(args):
    start = time.monotonic()
    root = Path(args.output)
    cfg = json.loads((root / 'config/model_spec_v1.json').read_text())
    pilot = Path(args.pilot)
    for role, item in cfg['bindings'].items():
        if sha256_file(item['path']) != item['sha256']:
            raise ValueError('Changed frozen input: ' + role)
    audit = json.loads((pilot / 'provenance/input_auditor_exit_v1.json').read_text())
    if audit['actual_exit_code'] != 0:
        raise ValueError('Pilot input audit not accepted')
    reader = FullContextAtlas(args.manifest)
    inventory = json.loads((pilot / 'outputs/pilot_input_inventory_v1.json').read_text())
    plan = json.loads((pilot / 'config/model_features_and_heldouts_v1.json').read_text())
    model_genes = plan['common_model_feature_candidates']
    heldout = {g for panels in plan['heldout_RNA_readouts'].values() for panel in panels.values() for g in panel}
    evidence_genes = sorted(heldout | {g for p in LINEAGE_PANELS.values() for g in p})
    pieces, evidence, detected, totals, sources, routes, fine, measured_per_source = [], [], [], [], [], [], [], {}
    def checkpoint(stage, **extra):
        save(root / 'checkpoints' / (str(time.time_ns()) + '.json'),
             {'observed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'stage': stage, 'wall_seconds': time.monotonic() - start, **extra})
    for row in inventory:
        ds = row['dataset']
        positions = np.load(pilot / row['selector_file'], allow_pickle=False)
        if sha256_file(pilot / row['selector_file']) != row['selector_sha256']:
            raise ValueError('Changed selected cells')
        lookup = reader.feature_lookup(ds)
        model_columns = feature_columns(lookup, model_genes, sorted(heldout))
        available = np.asarray([len(lookup.get(g, set())) == 1 for g in evidence_genes])
        columns = np.asarray([next(iter(lookup[g])) if available[i] else 0 for i, g in enumerate(evidence_genes)])
        measured_per_source[ds] = available
        local, ev, det, native = [], [], [], []
        z = reader.store(ds)
        key = reader.legacy_config['counts_key'] if reader.is_legacy(ds) else reader.primary.partitions[ds]['counts_key']
        for low in range(0, len(positions), cfg['read_block_rows']):
            chosen = positions[low:low + cfg['read_block_rows']]
            raw = gather_csr_rows(z[key], reader.rows(ds)[chosen])
            if (not np.isfinite(raw.data).all() or (raw.data < 0).any() or
                    not np.equal(raw.data, np.floor(raw.data)).all()):
                raise ValueError('Native raw counts fail integer assumption')
            total = np.asarray(raw.sum(axis=1, dtype=np.float64)).ravel()
            x, positive = normalized_log_counts(raw[:, model_columns], total)
            y = raw[:, columns].toarray().astype(np.float32)
            y[:, ~available] = 0
            fraction = np.zeros(len(total), dtype=np.float32)
            fraction[positive] = 10000 / total[positive]
            ev.append(np.log1p(y * fraction[:, None]))
            det.append(y > 0)
            native.append(total)
            local.append(x)
        labels = reader.frozen_labels(ds)
        major = labels['major'][positions].astype(str)
        pieces.append(sp.vstack(local, format='csr'))
        evidence.append(np.vstack(ev)); detected.append(np.vstack(det)); totals.append(np.concatenate(native))
        sources.extend([ds] * len(positions))
        routes.extend([ROUTES.get(v, 'unresolved_context') for v in major])
        fine.extend(labels['fine'][positions].astype(str))
        checkpoint('native_pilot_loaded', dataset=ds, selected_cells=len(positions))
    x = sp.vstack(pieces, format='csr'); del pieces
    ev = np.vstack(evidence); det = np.vstack(detected); del evidence, detected
    total = np.concatenate(totals)
    source = np.asarray(sources); route = np.asarray(routes); oldfine = np.asarray(fine)
    fit_rows = total > 0
    if fit_rows.sum() < len(total):
        checkpoint('zero_RNA_rows_preserved_but_not_modeled', n_excluded=int((~fit_rows).sum()))
    # These small position/label arrays are private, not exported to Git/chat.
    with (root / 'outputs/derived_axis_and_model_mask.npz').open('xb') as f:
        np.savez_compressed(f, modeled=fit_rows, source=source, inherited_route=route)
    hvg = select_features(x[fit_rows], source[fit_rows], cfg['n_features'])
    save(root / 'outputs/discovery_feature_selection_v1.json',
         {'all_source_discovery_not_holdout': True, 'method': 'mean_source_log_variance_over_mean_rank',
          'selected_feature_indices': hvg.tolist(), 'selected_feature_names': [model_genes[i] for i in hvg]})
    x = x[:, hvg]
    states = {}
    for name, genes in plan['heldout_RNA_readouts']['state'].items():
        value = np.full(len(total), np.nan, dtype=np.float32)
        for ds, available in measured_per_source.items():
            mask = source == ds
            cols = [evidence_genes.index(g) for g in genes if available[evidence_genes.index(g)]]
            if len(cols) >= cfg['minimum_state_genes_measured']:
                value[mask] = ev[mask][:, cols].mean(axis=1)
        states[name] = value
    dg = {g: evidence_genes.index(g) for g in ['TRDC', 'TRGC1', 'TRGC2', 'CD3D']}
    candidate = det[:, dg['TRDC']] & (det[:, dg['TRGC1']] | det[:, dg['TRGC2']]) & det[:, dg['CD3D']]
    anchors = {}
    rng = np.random.default_rng(cfg['diagnostic_sampling_seed'])
    for ds in np.unique(source):
        for r in np.unique(route):
            positions = np.flatnonzero((source == ds) & (route == r) & fit_rows)
            if len(positions) >= cfg['minimum_diagnostic_stratum']:
                anchors[(ds, r)] = np.sort(rng.choice(positions, min(len(positions), cfg['diagnostic_cells_per_stratum']), replace=False))
    results = []
    for seed in cfg['model_seeds']:
        partitions = {}
        global_raw, model = pca_fit(x[fit_rows], seed, cfg['components'], cfg['pca_iterations'])
        checkpoint('global_PCA_fit', seed=seed)
        global_corrected, harmony_info = harmonize(global_raw, source[fit_rows], seed, cfg)
        # Aligned A/B use precisely this same joint geometry; A annotations come
        # from local graphs and B from the joint graph. Geometry scores are one
        # observation, never two independent pieces of evidence.
        embeddings = {'raw_global_control': np.full((len(total), global_raw.shape[1]), np.nan, np.float32),
                      'A_B_shared_aligned': np.full((len(total), global_raw.shape[1]), np.nan, np.float32),
                      'A_local': np.full((len(total), cfg['components']), np.nan, np.float32),
                      'C_within_prior_route': np.full((len(total), cfg['components']), np.nan, np.float32)}
        embeddings['raw_global_control'][fit_rows] = global_raw
        embeddings['A_B_shared_aligned'][fit_rows] = global_corrected
        b, nb = cluster(global_corrected, seed, cfg)
        full = np.full(len(total), -1, np.int32); full[fit_rows] = b; partitions['B'] = full
        checkpoint('B_joint_graph_partitioned', seed=seed, n_partitions=nb)
        model_dir = root / 'models'; model_dir.mkdir(exist_ok=True)
        with (model_dir / ('global_discovery_PCA_' + str(seed) + '.npz')).open('xb') as f:
            np.savez_compressed(f, components=model.components_, mean=model.mean_, explained_variance=model.explained_variance_)
        for method, grouping, correct in [('A', source, False), ('C', route, True)]:
            labels = np.full(len(total), -1, np.int32)
            offset = 0
            for group in np.unique(grouping):
                member = np.flatnonzero(fit_rows & (grouping == group))
                if len(member) < cfg['minimum_fit_group']:
                    checkpoint('small_group_not_force_clustered', method=method, group=str(group), n_cells=len(member), seed=seed)
                    continue
                local, fitted = pca_fit(x[member], seed, cfg['components'], cfg['pca_iterations'])
                if correct:
                    local, _ = harmonize(local, source[member], seed, cfg)
                target = 'A_local' if method == 'A' else 'C_within_prior_route'
                embeddings[target][member, :local.shape[1]] = local
                # Small groups can have fewer PCs: pad missing coordinates with
                # zero, not NaN, without inventing extra fitted components.
                embeddings[target][member, local.shape[1]:] = 0
                lab, n = cluster(local, seed, cfg)
                labels[member] = lab + offset; offset += n
                checkpoint('local_graph_partitioned', method=method, group=str(group), n_cells=len(member), n_partitions=n, seed=seed)
            partitions[method] = labels
        summary = {'seed': seed, 'harmony': harmony_info, 'conditional_RNA_diagnostics': [],
                   'inherited_label_overlap_is_diagnostic_not_accuracy': True,
                   'C_routing_uses_inherited_coarse_labels_not_newly_accepted_truth': True}
        for (ds, r), positions in anchors.items():
            for name, embedding in embeddings.items():
                if not np.isfinite(embedding[positions]).all():
                    metrics = {'status': 'representation_unavailable_small_group'}
                else:
                    metrics = conditional_readout_metrics(embedding, states, total, candidate, positions, cfg['diagnostic_neighbors'])
                summary['conditional_RNA_diagnostics'].append({'source': str(ds), 'inherited_route': str(r), 'representation': name, **metrics})
        with (root / 'outputs' / ('candidate_partitions_' + str(seed) + '.npz')).open('xb') as f:
            np.savez_compressed(f, **partitions)
        tables = {}
        for method, labels in partitions.items():
            rows = []
            for code in np.unique(labels):
                member = labels == code
                rows.append({'partition_code': int(code), 'n_cells': int(member.sum()),
                             'inherited_route_counts': dict(collections.Counter(route[member])),
                             'source_counts': dict(collections.Counter(source[member])),
                             'inherited_fine_counts_not_truth': dict(collections.Counter(oldfine[member]))})
            tables[method] = rows
            # Source-specific measurement masks remain explicit; aggregate only
            # within a source. Do not treat source-unmeasured evidence as zero.
            for ds in np.unique(source):
                mask = (source == ds) & (labels >= 0) & fit_rows
                table = marker_summary(labels[mask], ev[mask], det[mask], evidence_genes, measured_per_source[ds])
                path = root / 'outputs' / ('native_partition_panels_' + str(seed) + '_' + method + '_' + ds + '.json.gz')
                with gzip.open(path, 'xt') as f:json.dump(table, f, allow_nan=False)
        save(root / 'outputs' / ('partition_crosswalk_diagnostic_' + str(seed) + '.json'), tables)
        save(root / 'outputs' / ('representation_diagnostics_' + str(seed) + '.json'), summary)
        results.append({'seed': seed, 'partitions': {m: len(set(v[v >= 0])) for m, v in partitions.items()}})
        checkpoint('seed_discovery_complete_not_annotation_acceptance', seed=seed)
        if sum(p.stat().st_size for p in root.rglob('*') if p.is_file()) > cfg['output_cap_bytes']:
            raise ValueError('Bounded pilot output cap exceeded')
    save(root / 'provenance/terminal_v1.json',
         {'status': 'CONTROLLED_DISCOVERY_PARTITIONS_COMPLETE_NOT_FINAL_ANNOTATIONS_OR_HOLDOUT',
          'completed_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'results': results,
          'wall_seconds': time.monotonic() - start, 'peak_RSS_bytes_Linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
          'native_counts_copied_to_disk': False, 'source_holdout_complete': False,
          'method_winner_selected': False, 'full_ABC_identity_validation_complete': False,
          'existing_authoritative_labels_or_CM_changed': False})


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for name in ['output', 'pilot', 'manifest']:p.add_argument('--' + name, required=True)
    a = p.parse_args()
    try:main(a)
    except Exception as error:
        save(Path(a.output) / 'provenance/failure_v1.json',
             {'status': 'FAILED_NOT_ACCEPTED', 'error_class': type(error).__name__, 'error': str(error)})
        raise
