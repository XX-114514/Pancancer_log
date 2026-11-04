# Batch-correction benchmark v1, step 2 (CPU, scanpy env): embeddings / graphs on the shared subsample.
#   harmony_source_t2       Harmony on X_pca, batch source, theta 2 (integration v2 recipe)
#   harmony_source_unit     Harmony source theta 2 + unit theta 1 (TME v1 recipe, rerun on subsample)
#   combat                  ComBat (batch source) on lognorm HVG, then PCA 50
#   bbknn_source            BBKNN on X_pca, batch source (3 neighbours per batch, trimmed)
#   harmony_bbknn           BBKNN on harmony_source_t2
# Embeddings -> OUT/emb_<name>.npy; BBKNN graphs -> OUT/graph_<name>.npz (connectivities, csr)
# Usage: run_cpu_methods.py SUB.h5ad OUT_DIR SEED
import sys, os, time, json
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc, harmonypy as hm, bbknn

subp, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
os.makedirs(out, exist_ok=True)
A = ad.read_h5ad(subp)
pca = A.obsm['X_pca'].astype(np.float64)
times = {}


def harmony(vars_, theta):
    h = hm.run_harmony(pca, A.obs[vars_].astype(str), vars_use=vars_, theta=theta, lamb=1.0, max_iter_harmony=20, random_state=seed)
    Z = np.asarray(h.Z_corr)
    return (Z.T if Z.shape[0] == pca.shape[1] else Z).astype(np.float32)


def timed(name, f):
    t = time.time(); r = f(); times[name] = round(time.time() - t, 1); print(name, times[name], 's', flush=True); return r


np.save(f'{out}/emb_unintegrated.npy', A.obsm['X_pca']); np.save(f'{out}/emb_v1_fulldata_harmony.npy', A.obsm['X_v1_harmony'])
H = timed('harmony_source_t2', lambda: harmony(['source'], 2.0)); np.save(f'{out}/emb_harmony_source_t2.npy', H)
np.save(f'{out}/emb_harmony_source_unit.npy', timed('harmony_source_unit', lambda: harmony(['source', 'unit_id'], [2.0, 1.0])))


def combat():
    B = ad.AnnData(X=A.layers['lognorm'].toarray(), obs=A.obs[['source']].copy())
    sc.pp.combat(B, key='source')
    sc.pp.pca(B, n_comps=50, svd_solver='arpack', random_state=seed)
    return B.obsm['X_pca'].astype(np.float32)


np.save(f'{out}/emb_combat.npy', timed('combat', combat))


def run_bbknn(rep):
    B = ad.AnnData(X=sp.csr_matrix((A.n_obs, 1), dtype=np.float32), obs=A.obs[['source']].copy())
    B.obsm['rep'] = rep
    bbknn.bbknn(B, batch_key='source', use_rep='rep', neighbors_within_batch=3, n_pcs=50, computation='pynndescent', trim=30)
    return B.obsp['connectivities'].tocsr()


sp.save_npz(f'{out}/graph_bbknn_source.npz', timed('bbknn_source', lambda: run_bbknn(A.obsm['X_pca'])))
sp.save_npz(f'{out}/graph_harmony_bbknn.npz', timed('harmony_bbknn', lambda: run_bbknn(H)))
json.dump(times, open(f'{out}/times_cpu.json', 'w'), indent=1)
