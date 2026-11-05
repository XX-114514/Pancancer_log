"""Align V7 CMs to frozen V5 K9 CMs by sample-activity (H) Pearson over shared samples and W cosine over shared node names."""
import sys, os
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import linear_sum_assignment
run, seed_dir = Path(sys.argv[1]), Path(sys.argv[2]); v5 = Path(os.environ['PC_ROOT']) / 'Project_v5/frozen/covarnet'
W7 = pd.read_csv(seed_dir / 'private_W_v1.tsv', sep='\t', index_col=0); H7 = pd.read_csv(seed_dir / 'private_H_v1.tsv', sep='\t', index_col=0)
W5 = pd.read_csv(v5 / 'nmf_W_subcluster_by_cm.tsv', sep='\t', index_col=0); H5 = pd.read_csv(v5 / 'nmf_H_cm_by_sample.tsv', sep='\t', index_col=0)
s = H7.columns.intersection(H5.columns); n = W7.index.intersection(W5.index)
Hc = pd.DataFrame(np.corrcoef(H7[s].values, H5[s].values)[:9, 9:], index=H7.index, columns=H5.index)
w7 = W7.loc[n] / np.linalg.norm(W7.loc[n], axis=0); w5 = W5.loc[n] / np.linalg.norm(W5.loc[n], axis=0)
Wc = pd.DataFrame(w7.T.values @ w5.values, index=W7.columns, columns=W5.columns)
comb = (Hc + Wc) / 2; r, c = linear_sum_assignment(-comb.values)
out = pd.DataFrame({'V7_CM': comb.index[r], 'V5_CM': comb.columns[c], 'H_pearson': Hc.values[r, c], 'W_cosine_shared_nodes': Wc.values[r, c], 'mean': comb.values[r, c]})
out['best_V5_by_H'] = Hc.idxmax(1).values[r]; out['best_V7_for_V5_by_H'] = Hc.idxmax(0).reindex(out.V5_CM).values
res = run / 'results'; res.mkdir(exist_ok=True); tag = f'{seed_dir.parent.name}_{seed_dir.name}'
out.to_csv(res / f'v7_v5_cm_alignment_{tag}.tsv', sep='\t', index=False)
Hc.to_csv(res / f'v7_v5_H_corr_{tag}.tsv', sep='\t'); Wc.to_csv(res / f'v7_v5_W_cosine_{tag}.tsv', sep='\t')
print(f'shared samples {len(s)}/{H7.shape[1]} V7, {H5.shape[1]} V5; shared nodes {len(n)}/{len(W7)} V7, {len(W5)} V5'); print(out.round(3).to_string(index=False))
