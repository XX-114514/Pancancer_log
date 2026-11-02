#!/usr/bin/env python3
"""Cellular-neighbourhood CM validation on single-cell spatial sections annotated with the CT75 scheme.

Per section, with labels = scheme nodes having >= --min-cells cells (rarer nodes
pooled as rare_<parent>):
1. Neighbourhood enrichment: radius graph (--radius um, default 50, about one
   Visium spot) between cell centroids; observed neighbour-pair counts per label
   pair vs --n-perm label permutations within the section -> z matrix.
2. Pseudo-spots: --bin um square bins (default 50) with >= 3 cells; per-bin label
   counts give the paper's Spearman colocalization and bivariate-Moran (8-NN bins)
   aggregation scores, as on Visium.
3. Per-cell neighbourhood CM activity: W-weighted member counts within the radius.
CM scores use only members the panel measures (held-out F1 >= --min-f1); a CM is
scored when >= 2 members are measurable and present. Random size-matched member
sets from measurable CT75 nodes give the null percentile.

Usage: cell_neighbourhood_cm_v1.py --annotation DIR --fit FIT_RUN --output DIR
"""
import argparse, glob, json, sys, zlib
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cm_spatial_scores_v1 import (load_cm, spearman_matrix, knn_weights, bivariate_moran_matrix,
                                  coloc_score, random_member_null, edge_values)


def nhood_z(pairs, lab, K, n_perm, rng):
    i, j = pairs[:, 0], pairs[:, 1]
    obs = np.bincount(lab[i] * K + lab[j], minlength=K * K) + np.bincount(lab[j] * K + lab[i], minlength=K * K)
    sims = np.empty((n_perm, K * K))
    for p in range(n_perm):
        l = rng.permutation(lab)
        sims[p] = np.bincount(l[i] * K + l[j], minlength=K * K) + np.bincount(l[j] * K + l[i], minlength=K * K)
    with np.errstate(invalid='ignore', divide='ignore'):
        z = (obs - sims.mean(0)) / sims.std(0)
    return z.reshape(K, K)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--annotation', required=True)
    ap.add_argument('--fit', required=True)
    ap.add_argument('--output', required=True)
    ap.add_argument('--radius', type=float, default=50.0)
    ap.add_argument('--bin', type=float, default=50.0)
    ap.add_argument('--min-cells', type=int, default=20)
    ap.add_argument('--min-f1', type=float, default=0.5)
    ap.add_argument('--n-perm', type=int, default=200)
    ap.add_argument('--n-null', type=int, default=1000)
    a = ap.parse_args()
    out = Path(a.output); out.mkdir(parents=True, exist_ok=True)
    W, members, edges = load_cm(a.fit)
    ct75 = W.index.tolist()
    rows, erows, mrows = [], [], []
    for f in sorted(glob.glob(str(Path(a.annotation) / '*.cells.tsv.gz'))):
        name = Path(f).name.replace('.cells.tsv.gz', '')
        c = pd.read_csv(f, sep='\t', index_col=0)
        key = c.panel.iloc[0]
        f1 = pd.read_csv(Path(a.annotation) / f'panel_{key}_heldout_node_F1.tsv', sep='\t').set_index('node').F1
        measurable = [n for n in ct75 if f1.get(n, 0) >= a.min_f1]
        vc = c.node.value_counts()
        lab_name = np.where(c.node.map(vc).to_numpy() >= a.min_cells, c.node.to_numpy(), 'rare_' + c.parent.astype(str).to_numpy())
        names = sorted(set(lab_name)); K = len(names); lab = pd.Series(range(K), index=names).loc[lab_name].to_numpy()
        xy = c[['x_centroid', 'y_centroid']].to_numpy()
        pairs = cKDTree(xy).query_pairs(a.radius, output_type='ndarray')
        rng = np.random.default_rng(zlib.crc32(name.encode()))
        Z = nhood_z(pairs, lab, K, a.n_perm, rng)
        # pseudo-spots
        b = np.floor(xy / a.bin).astype(np.int64); bid = pd.factorize(pd.Series(b[:, 0]).astype(str) + '_' + pd.Series(b[:, 1]).astype(str))[0]
        B = sp.csr_matrix((np.ones(len(lab)), (bid, lab)), shape=(bid.max() + 1, K)).toarray()
        ncell = B.sum(1); keep = ncell >= 3; B = B[keep]
        df = pd.DataFrame(xy, columns=['x', 'y']); df['b'] = bid; centers = df.groupby('b')[['x', 'y']].mean().to_numpy()[keep]
        C = spearman_matrix(B); M = bivariate_moran_matrix(B, knn_weights(centers, 8))
        # per-cell neighbourhood activity
        A = sp.coo_matrix((np.ones(len(pairs) * 2), (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])), shape=(len(c), len(c))).tocsr() + sp.eye(len(c))
        unit = dict(section=name, patient=str(c.patient.iloc[0]) if 'patient' in c else name, panel=key, n_cells=len(c),
                    n_pseudospots=int(keep.sum()), n_labels=K)
        for cm, mem in members.items():
            use = [m for m in mem if m in measurable and m in names]
            pool = [n for n in measurable if n in names]
            for stat, mat in [('neighbourhood_enrichment_z', Z), ('pseudospot_colocalization_Spearman', C), ('pseudospot_aggregation_Moran', M)]:
                score, r = coloc_score(mat, names, use) if len(use) >= 2 else (np.nan, np.nan)
                null = random_member_null(mat, names, use, pool, a.n_null, seed=zlib.crc32(f'{name}|{cm}|{stat}'.encode())) if len(use) >= 2 else np.array([])
                rows.append(dict(unit, cm=cm, statistic=stat, members_used=';'.join(use), n_members_used=len(use), score=score,
                                 within_CM_median=r, null_median=float(np.median(null)) if len(null) else np.nan,
                                 null_percentile=float((null < score).mean()) if len(null) else np.nan))
            onehot = np.isin(lab_name, use).astype(float) * pd.Series(lab_name).map(W[cm]).fillna(0).to_numpy()
            act = A @ onehot
            mrows.append(dict(unit, cm=cm, n_members_used=len(use), mean_cell_neighbourhood_activity=float(act.mean()),
                              frac_cells_any_member_within_radius=float((act > 0).mean())))
            for e in edges[cm]:
                ok = all(x in measurable and x in names for x in e)
                k = f'{e[0]}~{e[1]}'
                erows.append(dict(unit, cm=cm, edge=k, both_measurable=ok,
                                  nhood_z=edge_values(Z, names, [e])[k] if ok else np.nan,
                                  pseudospot_spearman=edge_values(C, names, [e])[k] if ok else np.nan,
                                  n_a=int(vc.get(e[0], 0)), n_b=int(vc.get(e[1], 0))))
        print(name, len(c), K, len(pairs), flush=True)
    pd.DataFrame(rows).to_csv(out / 'neighbourhood_CM_section_scores_v1.tsv', sep='\t', index=False)
    pd.DataFrame(erows).to_csv(out / 'neighbourhood_CM_edge_section_values_v1.tsv', sep='\t', index=False)
    pd.DataFrame(mrows).to_csv(out / 'neighbourhood_CM_activity_v1.tsv', sep='\t', index=False)
    json.dump(vars(a), open(out / 'neighbourhood_parameters_v1.json', 'w'), indent=2)


if __name__ == '__main__':
    main()
