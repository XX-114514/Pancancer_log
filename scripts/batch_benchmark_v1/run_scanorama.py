# Batch-correction benchmark v1, step 3b (scvi env + scanorama 1.7.4): Scanorama on lognorm HVG, batch source,
# 50-dim integrated embedding. Usage: run_scanorama.py SUB.h5ad OUT_DIR
import sys, time, json
import numpy as np, anndata as ad, scanorama

subp, out = sys.argv[1], sys.argv[2]
A = ad.read_h5ad(subp)
t = time.time()
src = A.obs['source'].astype(str).values
order = [np.flatnonzero(src == s) for s in dict.fromkeys(src)]
mats = [A.layers['lognorm'][i] for i in order]
genes = [list(A.var_names)] * len(mats)
integ, _ = scanorama.integrate(mats, genes, dimred=50, seed=0)
Z = np.zeros((A.n_obs, 50), np.float32)
for i, z in zip(order, integ):
    Z[i] = z
np.save(f'{out}/emb_scanorama.npy', Z)
json.dump({'scanorama': round(time.time() - t, 1)}, open(f'{out}/times_scanorama.json', 'w'))
print('scanorama', round(time.time() - t, 1), 's')
