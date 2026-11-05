"""Robustness of V7 primary K9 CMs: leave-one-source-out, source-equal weighting, seed, and nesting into the cophenetic-rule K4."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import linear_sum_assignment

run = Path(sys.argv[1]); o = run / 'outputs'; ref_seed = 'seed_20261103'
read = lambda d, m='W': pd.read_csv(d / f'private_{m}_v1.tsv', sep='\t', index_col=0)
W = read(o / 'covarnet_primary_strict' / ref_seed); H = read(o / 'covarnet_primary_strict' / ref_seed, 'H')

def match(Wb):
    n = W.index.intersection(Wb.index); C = np.corrcoef(W.loc[n].T.values, Wb.loc[n].T.values)[:W.shape[1], W.shape[1]:]
    r, c = linear_sum_assignment(-C); return pd.Series(C[r, c], index=W.columns[r]), pd.Series(Wb.columns[c], index=W.columns[r])

rows = {}
rows['seed_20261104'] = match(read(o / 'covarnet_primary_strict/seed_20261104'))[0]
se = o / 'covarnet_primary_strict_K9_source_equal' / ref_seed
if (se / 'private_W_v1.tsv').exists(): rows['source_equal_weights'] = match(read(se))[0]
loso = pd.read_csv(run / 'inputs/loso/loso_inputs.tsv', sep='\t')
lo = {}
for s in loso.source:
    d = o / 'loso' / s / ref_seed
    if (d / 'private_W_v1.tsv').exists(): lo[s] = match(read(d))[0]
lo = pd.DataFrame(lo).T
summary = pd.DataFrame(rows)
summary['LOSO_n'] = lo.notna().sum(); summary['LOSO_median_r'] = lo.median(); summary['LOSO_min_r'] = lo.min()
summary['LOSO_frac_r_ge_0.8'] = (lo >= 0.8).mean(); summary['LOSO_worst_source'] = lo.idxmin()
# K9 -> K4 nesting: each K9 CM's best K4 parent by W correlation and by sample activity
W4 = read(o / 'covarnet_primary_strict_K4' / ref_seed); H4 = read(o / 'covarnet_primary_strict_K4' / ref_seed, 'H')
Cw = pd.DataFrame(np.corrcoef(W.T.values, W4.loc[W.index].T.values)[:9, 9:], index=W.columns, columns=W4.columns)
Ch = pd.DataFrame(np.corrcoef(H.values, H4[H.columns].values)[:9, 9:], index=H.index, columns=H4.index)
summary['K4_parent_by_W'] = Cw.idxmax(1); summary['K4_parent_W_r'] = Cw.max(1); summary['K4_parent_by_H'] = Ch.idxmax(1); summary['K4_parent_H_r'] = Ch.max(1)
k4_seed = match_k4 = None
W4b = read(o / 'covarnet_primary_strict_K4/seed_20261104'); C4 = np.corrcoef(W4.T.values, W4b.loc[W4.index].T.values)[:4, 4:]
r, c = linear_sum_assignment(-C4)
out = run / 'results'; summary.round(3).to_csv(out / 'v7_primary_K9_robustness.tsv', sep='\t'); lo.round(3).to_csv(out / 'v7_primary_K9_LOSO_r.tsv', sep='\t')
Cw.round(3).to_csv(out / 'v7_K9_vs_K4_W_corr.tsv', sep='\t')
print(summary.round(2).to_string()); print('K4 seed stability r:', np.round(C4[r, c], 3))
