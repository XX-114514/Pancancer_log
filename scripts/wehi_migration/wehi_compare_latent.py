#!/usr/bin/env python3
"""Compare scVI runs that cannot be bitwise equal (different GPU architecture or GPU run-to-run).

Usage: wehi_compare_latent.py OUT_TSV NAME=OUTPUTS_DIR NAME=OUTPUTS_DIR [NAME=OUTPUTS_DIR ...]
Each OUTPUTS_DIR holds X_scVI_v1.npy and history_*.tsv from the same recipe and cell axis.
For every pair: per-epoch history relative differences (max and last epoch), the fraction of each
cell's 15 nearest latent neighbours shared between the two runs (mean and 5th percentile), and the
Spearman correlation of pairwise latent distances among 5000 fixed-seed cells.

Fixed 2026-10-31 before looking at results: with a WEHI replicate pair (same GPU type) as the
run-to-run reference, old-vs-WEHI agreement is accepted when, for every old-vs-WEHI pair,
knn15_shared_mean >= replicate value - 0.05 and pairwise_dist_spearman >= replicate value - 0.02.
These margins are project choices for a migration check, not scVI benchmarks.
"""
import glob, itertools, os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors

K, N_DIST, SEED = 15, 5000, 20261031

def knn(x):
    return NearestNeighbors(n_neighbors=K + 1).fit(x).kneighbors(x, return_distance=False)[:, 1:]

def main():
    out, runs = sys.argv[1], dict(a.split('=', 1) for a in sys.argv[2:])
    lat = {k: np.load(os.path.join(v, 'X_scVI_v1.npy')) for k, v in runs.items()}
    shapes = {v.shape for v in lat.values()}
    assert len(shapes) == 1, f'latent shapes differ: {shapes}'
    n = next(iter(shapes))[0]
    nb = {k: knn(v) for k, v in lat.items()}
    pick = np.random.default_rng(SEED).choice(n, size=min(N_DIST, n), replace=False)
    iu = np.triu_indices(len(pick), 1)
    dist = {k: np.linalg.norm(v[pick][:, None] - v[pick][None], axis=2)[iu] for k, v in lat.items()}
    rows = []
    for a, b in itertools.combinations(runs, 2):
        shared = np.array([len(set(x) & set(y)) / K for x, y in zip(nb[a], nb[b])])
        rho = spearmanr(dist[a], dist[b]).statistic
        hist = []
        for f in sorted(glob.glob(os.path.join(runs[a], 'history_*.tsv'))):
            g = os.path.join(runs[b], os.path.basename(f))
            if not os.path.exists(g):
                continue
            x = pd.read_csv(f, sep='\t', index_col=0).iloc[:, 0].to_numpy(float)
            y = pd.read_csv(g, sep='\t', index_col=0).iloc[:, 0].to_numpy(float)
            if len(x) != len(y):
                hist.append((os.path.basename(f), np.nan, np.nan)); continue
            rel = np.abs(x - y) / np.maximum(np.abs(x), 1e-12)
            hist.append((os.path.basename(f), float(rel.max()), float(rel[-1])))
        rows.append(dict(pair=f'{a}~{b}', n_cells=n, knn15_shared_mean=shared.mean(),
                         knn15_shared_p05=np.quantile(shared, .05), pairwise_dist_spearman=rho,
                         max_abs_latent_diff=float(np.abs(lat[a] - lat[b]).max()),
                         history_max_rel_diff=max((h[1] for h in hist), default=np.nan),
                         history_last_epoch_max_rel_diff=max((h[2] for h in hist), default=np.nan),
                         history_detail=';'.join(f'{h[0]}:{h[1]:.3g}/{h[2]:.3g}' for h in hist)))
    df = pd.DataFrame(rows)
    df.to_csv(out, sep='\t', index=False)
    print(df.drop(columns='history_detail').to_string(index=False))

if __name__ == '__main__':
    main()
