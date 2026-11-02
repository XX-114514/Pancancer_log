"""CoVarNet spatial CM scores shared by the Visium and single-cell neighbourhood analyses.

Definitions follow Shi et al., Nature 2025, Methods:
- CM activity: component-subset abundances weighted by the NMF basis W, then
  scaled so the 99th percentile over all CMs (pooled values) equals 1.
- Colocalization score: for a CM, S = Spearman correlations of all subset pairs
  with at least one subset in the CM; r = median correlation of pairs with both
  subsets in the CM; score = fraction of S that is <= r.
- Aggregation score: as above with global bivariate Moran's I (spdep moran_bv:
  scaled variables, row-standardised weights) instead of Spearman.
Additions here, not in the paper: a size-matched random-member null and
per-edge sign summaries across sections.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree
from scipy.stats import rankdata


def load_cm(fit_dir, seed=20261030):
    """W (nodes x CMs) and author-network members of the frozen CT75 fit."""
    d = Path(fit_dir) / 'outputs' / f'seed_{seed}'
    W = pd.read_csv(d / 'private_W_v1.tsv', sep='\t', index_col=0)
    net = pd.read_csv(d / 'author_compatibility_each_node_v1.tsv', sep='\t')
    edges = pd.read_csv(d / 'author_compatibility_each_edge_v1.tsv', sep='\t')
    members = {cm: sorted(g.subCluster) for cm, g in net.groupby('cm')}
    edge_list = {cm: [tuple(sorted((r.subCluster1, r.subCluster2))) for r in g.itertuples()] for cm, g in edges.groupby('cm')}
    return W, members, edge_list


def cm_activity(abund, W, members):
    """abund: DataFrame units x subsets. Returns unscaled W-weighted member sums."""
    out = {}
    for cm, nodes in members.items():
        use = [n for n in nodes if n in abund.columns]
        out[cm] = abund[use].to_numpy() @ W.loc[use, cm].to_numpy() if use else np.full(len(abund), np.nan)
    return pd.DataFrame(out, index=abund.index)


def scale_activity(act):
    q = np.nanpercentile(act.to_numpy().ravel(), 99)
    return act / q if q > 0 else act


def spearman_matrix(X):
    R = np.apply_along_axis(rankdata, 0, X)
    R = R - R.mean(axis=0)
    sd = np.sqrt((R ** 2).sum(axis=0))
    with np.errstate(invalid='ignore', divide='ignore'):
        C = (R.T @ R) / np.outer(sd, sd)
    return C


def knn_weights(coords, k=6):
    tree = cKDTree(coords)
    _, idx = tree.query(coords, k=k + 1)
    n = len(coords); rows = np.repeat(np.arange(n), k); cols = idx[:, 1:].ravel()
    return sp.csr_matrix((np.full(n * k, 1.0 / k), (rows, cols)), shape=(n, n))


def bivariate_moran_matrix(X, Wsp):
    """Symmetrised global bivariate Moran's I for all column pairs."""
    with np.errstate(invalid='ignore', divide='ignore'):
        Z = (X - X.mean(axis=0)) / X.std(axis=0)
    Z = np.nan_to_num(Z)
    M = (Z.T @ (Wsp @ Z)) / len(Z)
    return (M + M.T) / 2


def coloc_score(C, names, members):
    """Paper score from a pair-statistic matrix C over subsets `names`."""
    ix = [names.index(m) for m in members if m in names]
    if len(ix) < 2:
        return np.nan, np.nan
    inside = np.zeros(len(names), bool); inside[ix] = True
    iu = np.triu_indices(len(names), 1)
    keep = inside[iu[0]] | inside[iu[1]]
    S = C[iu][keep]; S = S[np.isfinite(S)]
    both = inside[iu[0]] & inside[iu[1]]
    within = C[iu][both]; within = within[np.isfinite(within)]
    if not len(within) or not len(S):
        return np.nan, np.nan
    r = np.median(within)
    return float((S <= r).mean()), float(r)


def random_member_null(C, names, members, pool, n_draw=1000, seed=0, strata=None):
    """Scores for size-matched random member sets drawn from `pool` (excluding real members).

    With `strata` (node -> parent compartment) each draw keeps the members'
    parent composition, so same-lineage co-occurrence is controlled.
    """
    rng = np.random.default_rng(seed)
    use = [m for m in members if m in names]
    cand = [p for p in pool if p in names and p not in members]
    if len(use) < 2 or len(cand) < len(use):
        return np.array([])
    if strata is None:
        return np.array([coloc_score(C, names, list(rng.choice(cand, len(use), replace=False)))[0] for _ in range(n_draw)])
    need = pd.Series([strata[m] for m in use]).value_counts()
    by = {k: [c for c in cand if strata.get(c) == k] for k in need.index}
    if any(len(by[k]) < n for k, n in need.items()):
        return np.array([])
    return np.array([coloc_score(C, names, [x for k, n in need.items() for x in rng.choice(by[k], n, replace=False)])[0]
                     for _ in range(n_draw)])


def edge_values(C, names, edges):
    return {f'{a}~{b}': (float(C[names.index(a), names.index(b)]) if a in names and b in names else np.nan) for a, b in edges}
