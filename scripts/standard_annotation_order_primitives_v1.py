"""Standard count-HVG/PCA and source-conditional RNA diagnostics.

Native typing rules remain provisional RNA evidence. No clinical endpoint,
inherited fine identity, or method result chooses the feature set or parameters.
"""
import re
from types import SimpleNamespace

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import pearsonr, spearmanr
from sklearn.neighbors import NearestNeighbors

from TME_exact_globin_feature_exclusion_v1 import is_exact_globin_symbol
from TME_all_lineage_cluster_policy_v1 import GATES, GATE_ORDER, BROAD, annotate_cluster
from TME_assay_aware_T_anchor_policy_v1 import ACCESSORY, SIGNAL, CHAIN


def technical(gene):
    return (gene.upper() == 'MALAT1' or is_exact_globin_symbol(gene) or
            bool(re.search(r'^(MT-|RP[SL][0-9]|HB[ABDEGQMZ][0-9]|ERCC-|IG[HKL][VDJC]|TR[ABGD][VJ])', gene, re.I)))


def select_features(counts, sources, number, genes, candidates):
    import anndata as ad
    import scanpy as sc
    if (not sp.issparse(counts) or counts.shape != (len(sources), len(genes)) or
            len(set(genes)) != len(genes) or not np.isfinite(counts.data).all() or
            (counts.data < 0).any() or not np.equal(counts.data, np.floor(counts.data)).all()):
        raise ValueError('Seurat-v3 must receive original sparse integer counts with explicit axes')
    a = ad.AnnData(counts.astype(np.float32),
                   obs=pd.DataFrame({'source': pd.Categorical(sources)}),
                   var=pd.DataFrame(index=genes))
    sc.pp.highly_variable_genes(a, flavor='seurat_v3', batch_key='source',
                               n_top_genes=candidates, check_values=True)
    table = a.var.copy()
    table['technical_gene'] = [technical(g) for g in genes]
    eligible = table.loc[~table.technical_gene & table.highly_variable_rank.notna()].sort_values(
        ['highly_variable_rank', 'highly_variable_nbatches'], ascending=[True, False], kind='mergesort')
    if len(eligible) < number:
        raise ValueError('Insufficient qualified official HVGs; no dispersion fallback')
    selected = eligible.index[:number].tolist()
    table['selected_HVG'] = table.index.isin(selected)
    return np.asarray([genes.index(g) for g in selected], np.int64), table


def pca_fit(x, seed, components, iterations):
    import anndata as ad
    import scanpy as sc
    n = min(components, x.shape[0] - 1, x.shape[1] - 1)
    if n < 2:
        raise ValueError('Insufficient rows/features for centered PCA')
    a = ad.AnnData(x.tocsr().astype(np.float32))
    sc.pp.pca(a, n_comps=n, zero_center=True, svd_solver='arpack',
              random_state=seed, dtype='float32', mask_var=None)
    z = np.asarray(a.obsm['X_pca'], np.float32)
    assert np.isfinite(z).all()
    model = SimpleNamespace(components_=np.asarray(a.varm['PCs']).T,
        mean_=np.asarray(a.X.mean(axis=0)).ravel(), explained_variance_=a.uns['pca']['variance'])
    return z, model


COARSE_ROUTING_GENES = frozenset(ACCESSORY + SIGNAL + CHAIN + [
    g for name in BROAD for module in GATES[name] for g in module])


def coarse_routes_from_local(source, partitions, flags, assessed):
    """Source-local fine groups -> RNA coarse routing, including unresolved."""
    routes = np.full(len(source), 'unresolved_context', dtype='<U40')
    records = []
    for ds in np.unique(source):
        for k in np.unique(partitions[source == ds]):
            if k < 0:
                continue
            member = (source == ds) & (partitions == k)
            counts = dict(zip(GATE_ORDER, flags[member].sum(axis=0).astype(int).tolist()))
            proposal = annotate_cluster(int(member.sum()), counts, assessed[ds])
            major = proposal['major']
            route = 'T_NK' if major in ('T', 'NK', 'T_NK') else major
            if route == 'unresolved':
                route = 'unresolved_context'
            routes[member] = route
            records.append(dict(source=str(ds),source_local_partition=int(k),n_cells=int(member.sum()),
                coarse_RNA_route=route,status=proposal['status'],not_independent_identity_truth=True))
    return routes, records


def conditional_readout_metrics(embedding, states, native_totals, candidate, indices, neighbors, *, space_ids=None):
    import scanpy as sc
    indices = np.asarray(indices, dtype=np.int64)
    if len(indices) < max(30, neighbors + 2):
        return {'status': 'insufficient_conditional_cells'}
    if space_ids is not None and len(np.unique(np.asarray(space_ids)[indices])) != 1:
        raise ValueError('Neighbor diagnostics cannot mix separately fitted coordinate spaces')
    z = embedding[indices]
    _, near = NearestNeighbors(n_neighbors=neighbors + 1, n_jobs=4).fit(z).kneighbors(z)
    near = np.asarray([row[row != i][:neighbors] for i, row in enumerate(near)])
    assert near.shape == (len(indices), neighbors)
    graph = sp.csr_matrix((np.ones(near.size),
        (np.repeat(np.arange(len(near)), neighbors), near.ravel())), shape=(len(near), len(near)))
    assert not graph.diagonal().any()
    depth = np.log1p(np.asarray(native_totals)[indices])
    design = np.column_stack([np.ones(len(depth)), depth - depth.mean()])
    out = dict(status='conditional_reserved_RNA_diagnostic_not_orthogonal_truth',
               n_cells=len(indices), states={})
    for name, values in states.items():
        y = np.asarray(values)[indices]
        if not np.isfinite(y).all():
            out['states'][name] = {'status': 'measurement_or_unused_gene_coverage_incomplete'}
            continue
        residual = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
        pred = residual[near].mean(axis=1)
        if residual.std() < 1e-8 or pred.std() < 1e-8:
            out['states'][name] = {'status': 'no_variable_signal'}
        else:
            out['states'][name] = dict(status='estimated',
                depth_adjusted_neighbor_Pearson=float(pearsonr(residual, pred).statistic),
                depth_adjusted_neighbor_Spearman=float(spearmanr(residual, pred).statistic),
                depth_adjusted_official_Scanpy_Morans_I=float(sc.metrics.morans_i(graph, residual)))
    rare = np.asarray(candidate)[indices]
    if not np.isfinite(rare).all():
        out['gamma_delta_RNA_candidate'] = dict(status='unmeasured_not_negative',
            independent_identity_accuracy_not_claimed=True)
    else:
        rare = rare.astype(bool)
        n = int(rare.sum())
        out['gamma_delta_RNA_candidate'] = dict(status='conditional_RNA_candidate_diagnostic',
            n_candidates_not_verified_identity=n,
            conditional_neighbor_candidate_fraction=float(rare[near[rare]].mean()) if n else None,
            conditional_candidate_prevalence=float(rare.mean()),
            candidate_RNA_genes_overlap_coarse_routing_not_independent_identity_validation=True)
    return out
