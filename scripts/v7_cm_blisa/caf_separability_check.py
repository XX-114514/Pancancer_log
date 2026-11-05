"""Fibroblast-only check: are V7 CAF labels separable with all reference genes? Held-out-sample confusion + marker means."""
import sys
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import f1_score, confusion_matrix
run = Path(sys.argv[1])
ref = ad.read_h5ad(run / 'reference/v7_spatial_reference_counts.h5ad')
f = ref[ref.obs.label.str.startswith('Fibroblast_') | (ref.obs.label == 'Other_Fibroblast')].copy()
X = f.X.astype(np.float64); s = np.asarray(X.sum(1)).ravel(); s[s == 0] = 1; X = np.log1p(sp.diags(1e4 / s) @ X).toarray()
tr = f.obs.split.eq('train').to_numpy(); y = f.obs.label.astype(str).to_numpy()
sc = StandardScaler().fit(X[tr]); clf = LogisticRegression(max_iter=2000, C=0.1).fit(sc.transform(X[tr]), y[tr])
p = clf.predict(sc.transform(X[~tr])); labs = clf.classes_
print('genes', X.shape[1]); print(pd.Series(f1_score(y[~tr], p, labels=labs, average=None), index=labs).round(2).to_string())
cm = pd.DataFrame(confusion_matrix(y[~tr], p, labels=labs), index=labs, columns=labs); print((cm.div(cm.sum(1), axis=0) * 100).round(0).to_string())
g = [x for x in ['ACTA2', 'TAGLN', 'MYH11', 'RGS5', 'COL1A1', 'COL3A1', 'POSTN', 'FAP', 'CXCL12', 'IL6', 'TOP2A', 'MKI67'] if x in f.var_names]
m = pd.DataFrame(X[:, [list(f.var_names).index(x) for x in g]], columns=g); m['label'] = y
print(m.groupby('label').mean().round(2).to_string())
