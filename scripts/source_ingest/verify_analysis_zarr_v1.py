# Independent read-back of pancancer_analysis_ready_v1.zarr against its input shards:
# unique cell ids, integer counts, per-source cell numbers and total-count sums, measured mask shape.
# Usage: verify_analysis_zarr_v1.py ZARR EXISTING_DIR SUPPLEMENT_QC_DIR OUT_JSON
import sys, os, glob, json
import numpy as np, pandas as pd, zarr, anndata as ad
from anndata.io import read_elem

z, ed, sd, outj = sys.argv[1:5]
g = zarr.open(z, 'r')
obs = read_elem(g['obs']); var = read_elem(g['var'])
X = g['X']; indptr = X['indptr'][:]
res = {'n_obs': len(obs), 'n_var': len(var), 'cell_id_unique': bool(obs['cell_id'].is_unique),
       'X_dtype': str(X['data'].dtype), 'X_shape': X.attrs['shape'],
       'measured_mask_shape': list(g['varm']['measured_in_source'].shape)}
checks, ok = [], True
src = obs['source'].astype(str).values
starts = np.r_[0, np.flatnonzero(src[1:] != src[:-1]) + 1, len(src)]
blocks = {src[a]: (a, b) for a, b in zip(starts[:-1], starts[1:])}
for f in sorted(glob.glob(f'{ed}/*.h5ad') + glob.glob(f'{sd}/*.h5ad')):
    s = os.path.basename(f)[:-5]
    a = ad.read_h5ad(f, backed='r')
    n = a.n_obs
    lo, hi = blocks[s]
    data = X['data'][indptr[lo]:indptr[hi]]
    tot_z = int(data.sum()); tot_s = int(a.X[:].sum())
    integer = bool(np.all(data == np.round(data)))
    c = {'source': s, 'cells_shard': n, 'cells_zarr': hi - lo, 'sum_shard': tot_s, 'sum_zarr': tot_z, 'integer': integer}
    c['ok'] = c['cells_shard'] == c['cells_zarr'] and tot_s == tot_z and integer
    ok &= c['ok']; checks.append(c)
    a.file.close()
res['per_source'] = checks
res['all_ok'] = bool(ok and res['cell_id_unique'])
json.dump(res, open(outj, 'w'), indent=1, default=int)
print(json.dumps({k: v for k, v in res.items() if k != 'per_source'}, default=int))
print([c['source'] for c in checks if not c['ok']])
