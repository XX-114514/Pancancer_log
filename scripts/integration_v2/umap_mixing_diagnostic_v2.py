# Diagnostic for integration v2: per compartment (marker-based, coarse), how source-specific are the
# Harmony neighbourhoods? Subsample N cells, 30-NN on X_pca_harmony among the subsample, fraction of
# neighbours from the same source vs the expectation under perfect mixing. Read-only.
# Usage: umap_mixing_diagnostic_v2.py RUN_DIR OUT_JSON [N]
import sys, glob, json
import numpy as np, pandas as pd, anndata as ad, zarr
from anndata.io import read_elem
from sklearn.neighbors import NearestNeighbors

run, outj = sys.argv[1:3]
N = int(sys.argv[3]) if len(sys.argv) > 3 else 300000
g = zarr.open(f'{run}/integration_v2.zarr', 'r')
src = np.asarray(read_elem(g['obs']['source'])).astype(str)
tech = np.asarray(read_elem(g['obs']['technology'])).astype(str)
H = g['obsm']['X_pca_harmony'][:]
M = {'immune': ['PTPRC'], 'epithelial': ['EPCAM', 'KRT19', 'KRT8', 'KRT18'], 'stromal': ['COL1A1', 'COL1A2', 'PECAM1', 'VWF', 'ACTA2'],
     'melanocytic': ['PMEL', 'MLANA'], 'glial': ['GFAP', 'OLIG2', 'SOX2']}
genes = sorted({x for v in M.values() for x in v})
vals = {k: [] for k in genes}
for f in sorted(glob.glob(f'{run}/shards/*.h5ad')):
    a = ad.read_h5ad(f, backed='r')
    idx = [a.var_names.get_loc(x) if x in a.var_names else -1 for x in genes]
    X = a.X[:, [i for i in idx if i >= 0]]; X = X.toarray() if hasattr(X, 'toarray') else np.asarray(X)
    k = 0
    for x, i in zip(genes, idx):
        vals[x].append(X[:, k] if i >= 0 else np.zeros(a.n_obs, np.float32)); k += i >= 0
    a.file.close()
E = pd.DataFrame({x: np.concatenate(v) for x, v in vals.items()})
score = pd.DataFrame({c: E[v].mean(1) for c, v in M.items()})
comp = np.where(score.max(1) < 0.5, 'unassigned', score.idxmax(1))
rng = np.random.default_rng(0); s = rng.choice(len(src), size=min(N, len(src)), replace=False)
nn = NearestNeighbors(n_neighbors=31).fit(H[s]); _, ind = nn.kneighbors(H[s])
same = (src[s][ind[:, 1:]] == src[s][:, None]).mean(1)
freq = pd.Series(src[s]).value_counts(normalize=True)
expected = freq.reindex(src[s]).values
d = pd.DataFrame({'comp': comp[s], 'same': same, 'expected': expected, 'tech': tech[s]})
res = {'n_sub': int(len(s)), 'by_compartment': d.groupby('comp').agg(cells=('same', 'size'), same_source_frac_median=('same', 'median'),
       same_source_frac_mean=('same', 'mean'), expected_mean=('expected', 'mean')).round(3).to_dict('index'),
       'compartment_fraction_all_cells': pd.Series(comp).value_counts(normalize=True).round(3).to_dict(),
       'by_technology_same_source_mean': d.groupby('tech')['same'].mean().round(3).to_dict()}
json.dump(res, open(outj, 'w'), indent=1); print(json.dumps(res, indent=1))
