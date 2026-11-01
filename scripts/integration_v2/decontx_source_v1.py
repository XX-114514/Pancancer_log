# Uniform ambient-RNA correction (decontX) for one source of the analysis-ready set, run per unit
# (sample / qualification unit) with the project's existing worker scripts/decontx_library_worker.R.
# Cluster labels for decontX come from a quick per-unit clustering (normalize_total 1e4, log1p,
# 2000 HVG, 30 PCs, 15-NN, Leiden 1.0); units with <2 clusters use a 2-way k-means on the PCs.
# Output: one h5ad per source with X = decontX counts (float32), layers['counts'] = raw counts,
# obs['decontX_contamination'] / ['decontX_status']. Cells are not removed here.
# Usage: decontx_source_v1.py SHARD.h5ad WORK_DIR OUT.h5ad RSCRIPT WORKER_R SEED
import sys, os, json, subprocess
import numpy as np, pandas as pd, scipy.sparse as sp, scipy.io, anndata as ad, scanpy as sc

shard, work, out, rscript, worker, seed = sys.argv[1:7]
seed = int(seed)
a = ad.read_h5ad(shard)
X = a.X.tocsr().astype(np.float32)
blocks, cont, status = [], np.full(a.n_obs, np.nan, np.float32), np.array(['not_run'] * a.n_obs, dtype=object)
order = []
for u, idx in a.obs.groupby('unit_id', observed=True).indices.items():
    idx = np.asarray(idx)
    sub = X[idx]
    keep_g = np.asarray((sub > 0).sum(0)).ravel() > 0
    sub = sub[:, keep_g]
    genes = a.var_names[keep_g].astype(str)
    d = os.path.join(work, f'u{len(order):04d}')
    os.makedirs(d, exist_ok=True)
    b = ad.AnnData(X=sub.copy())
    sc.pp.normalize_total(b, target_sum=1e4); sc.pp.log1p(b)
    sc.pp.highly_variable_genes(b, n_top_genes=min(2000, b.n_vars), flavor='seurat')
    sc.pp.pca(b, n_comps=min(30, b.n_obs - 1, int(b.var.highly_variable.sum()) - 1), mask_var='highly_variable', random_state=seed)
    sc.pp.neighbors(b, n_neighbors=15, random_state=seed)
    sc.tl.leiden(b, resolution=1.0, random_state=seed, flavor='igraph', n_iterations=2, directed=False)
    z = b.obs['leiden'].astype(str).values
    if len(set(z)) < 2:
        from sklearn.cluster import KMeans
        z = KMeans(2, random_state=seed, n_init=10).fit_predict(b.obsm['X_pca']).astype(str)
    scipy.io.mmwrite(os.path.join(d, 'input_counts_private.mtx'), sp.csr_matrix(sub).T.tocoo())
    pd.DataFrame({'stable_cell_id': a.obs_names[idx].astype(str), 'broad_group': z}).to_csv(os.path.join(d, 'cells_private.tsv'), sep='\t', index=False)
    open(os.path.join(d, 'genes_private.tsv'), 'w').write('\n'.join(genes) + '\n')
    r = subprocess.run([rscript, worker, d, str(seed)], capture_output=True, text=True)
    if r.returncode != 0:
        status[idx] = 'decontX_failed'
        print(json.dumps({'unit': u, 'status': 'failed', 'stderr': r.stderr[-300:]}), flush=True)
        blocks.append((idx, keep_g, None))
    else:
        cm = scipy.io.mmread(os.path.join(d, 'corrected_counts_private.mtx')).T.tocsr().astype(np.float32)
        info = pd.read_csv(os.path.join(d, 'decontX_cells_private.tsv'), sep='\t')
        cont[idx] = info['contamination'].values
        status[idx] = 'decontX_done'
        blocks.append((idx, keep_g, cm))
        s = json.load(open(os.path.join(d, 'decontX_summary.json')))
        print(json.dumps({'unit': u, 'cells': len(idx), 'median_contamination': s['median_contamination'], 'converged': s['converged']}), flush=True)
    order.append(u)
    for f in ('input_counts_private.mtx', 'corrected_counts_private.mtx'):
        if os.path.exists(os.path.join(d, f)):
            os.remove(os.path.join(d, f))
# assemble decontX matrix in the shard's cell and gene order; failed units keep raw counts
parts = []
rowsel = np.empty(a.n_obs, dtype=np.int64)
pos = 0
for idx, keep_g, cm in blocks:
    if cm is None:
        m = X[idx]
    else:
        cols = np.flatnonzero(keep_g)
        cm = cm.tocoo()
        m = sp.csr_matrix((cm.data, (cm.row, cols[cm.col])), shape=(len(idx), X.shape[1]), dtype=np.float32)
    parts.append(m); rowsel[idx] = np.arange(pos, pos + len(idx)); pos += len(idx)
D = sp.vstack(parts).tocsr()[rowsel]
D.eliminate_zeros()
a.layers['counts'] = a.X
a.X = D
a.obs['decontX_contamination'] = cont
a.obs['decontX_status'] = pd.Categorical(status)
a.write_h5ad(out, compression='gzip')
print(json.dumps({'source': os.path.basename(shard), 'cells': int(a.n_obs), 'units': len(order),
                  'failed_units': int(sum(c is None for _, _, c in blocks)),
                  'median_contamination': float(np.nanmedian(cont))}), flush=True)
