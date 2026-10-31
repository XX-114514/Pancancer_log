# Read-only: classify non-RNA features in the core fullgene store against the source's published
# 10x feature_types, and measure how much of each cell's counts they carry.
# Usage: core_feature_type_audit.py <features_dir with *_features.tsv.gz> <harmonized_per_gse dir> <core zarr> <GSE>
import sys, glob, os, json
import numpy as np, pandas as pd, zarr, anndata as ad
from anndata.io import read_elem

feat_dir, per_gse, core, gse = sys.argv[1:5]
types, adt = {}, set()
for f in sorted(glob.glob(os.path.join(feat_dir, '*_features.tsv.gz'))):
    d = pd.read_csv(f, sep='\t', header=None)
    for t, n in d[2].value_counts().items():
        types[t] = types.get(t, 0) + int(n)
    adt |= {x.upper() for x in d.loc[d[2] == 'Antibody Capture', 1]}

cv = read_elem(zarr.open(core, 'r')['var'])
names = pd.Series(cv.index.astype(str))
is_adt = names.str.upper().isin(adt).values
a = ad.read_h5ad(os.path.join(per_gse, f'{gse}_symbol_counts.h5ad'))
m = np.array([v.upper() in adt for v in a.var_names])
tot = np.asarray(a.X.sum(1)).ravel()
frac = np.asarray(a.X[:, m].sum(1)).ravel() / np.maximum(tot, 1)
print(json.dumps({
    'gse': gse, 'n_feature_files': len(glob.glob(os.path.join(feat_dir, '*_features.tsv.gz'))),
    'published_feature_type_rows': types, 'distinct_antibody_capture_features': len(adt),
    'core_features_matching': int(is_adt.sum()),
    'core_measured_in_n_gse_of_matches': {str(k): int(v) for k, v in cv.loc[is_adt, 'measured_in_n_gse'].value_counts().items()},
    'source_cells': int(a.n_obs), 'cells_with_antibody_counts': int((frac > 0).sum()),
    'antibody_count_fraction_max': float(frac.max()),
}, indent=1))
