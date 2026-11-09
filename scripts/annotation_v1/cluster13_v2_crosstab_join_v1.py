import sys, h5py, pandas as pd
from anndata.io import read_elem
sub, lab, out = sys.argv[1:4]
with h5py.File(sub, 'r') as f:
    ob = read_elem(f['obs'])
print([c for c in ob.columns if 'sub' in c or 'tme' in c or 'covar' in c])
ob = ob.reset_index(drop=True)
print(ob.covarnet_v1_include.value_counts().to_dict() if 'covarnet_v1_include' in ob else 'no include')
fl = [c for c in ob.columns if c in ('tme_flag',)][0]; ob = ob[ob[fl].astype(str) == 'non_malignant'][['cell_id', 'sub_r1']].copy()
L = pd.read_csv(lab, sep='\t', low_memory=False, dtype=str)
print(L.columns.tolist(), len(L))
m = ob.merge(L, on='cell_id', how='left', validate='one_to_one')
assert len(m) == 68374; print('missing', int(m.major_label.isna().sum()))
m.to_pickle(out)
print(m.major_label.value_counts(dropna=False).head(30)); print(m.group.value_counts(dropna=False))
for c in ['r1_action','r2_action','origin_pool']: print(m[c].value_counts(dropna=False).head(12).to_dict())
