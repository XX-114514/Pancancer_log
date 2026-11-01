# Integration v2 of the cleaned analysis-ready set (decontX-corrected per source; see
# docs/ANALYSIS_READY_ZARR_V1.md and the integration v2 run record). Recipe follows the previous TME
# standard integration (20261029_TME_standard_HVG_PCA_Harmony_v1: seurat_v3 HVG 6500 -> 4000
# non-technical, normalize_total 1e4 + log1p, centred ARPACK PCA 50, Harmony by source theta 2,
# lambda 1, 20 iterations), then 15-NN on Harmony and UMAP.
# Cells: decontX finished and contamination < 0.5. Genes for HVG: measured in every source.
# seurat_v3 is computed per source on raw counts and combined exactly as scanpy's batch_key does
# (per-batch rank among the top n, nanmedian rank, number of batches), so sources never share memory.
# Outputs: OUT/integration_v2.zarr (X = log1p-normalised decontX counts, layers counts / decontx,
# obsm X_pca / X_pca_harmony / X_umap), OUT/figures/*.png, OUT/report/*.
# Usage: integrate_v2.py DECONTX_DIR OUT_DIR REPO SEED
import sys, os, glob, json, re, gc, time, datetime, subprocess
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc, harmonypy as hm, zarr
from anndata.io import read_elem, write_elem
from anndata.experimental import concat_on_disk

ddir, out, repo, seed = sys.argv[1:5]
seed = int(seed)
N_CAND, N_HVG, N_PC, MAX_CONT = 6500, 4000, 50, 0.5
os.makedirs(f'{out}/shards', exist_ok=True); os.makedirs(f'{out}/figures', exist_ok=True); os.makedirs(f'{out}/report', exist_ok=True)
log = lambda **k: print(json.dumps({'t': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'), **k}, default=str), flush=True)
technical = lambda g: g.upper() == 'MALAT1' or bool(re.search(r'^(MT-|RP[SL][0-9]|HB[ABDEGQMZ][0-9]|ERCC-)', g, re.I))
files = sorted(glob.glob(f'{ddir}/*.h5ad'))
sources = [os.path.basename(f)[:-5] for f in files]

# common genes
gsets = [set(ad.read_h5ad(f, backed='r').var_names.astype(str)) for f in files]
common = sorted(set.intersection(*gsets))
log(step='common_genes', n=len(common), sources=len(sources))


def load(f):
    a = ad.read_h5ad(f)
    keep = (a.obs['decontX_status'].astype(str) == 'decontX_done').values & (a.obs['decontX_contamination'].values < MAX_CONT)
    tot = np.asarray(a.X.sum(1)).ravel()
    keep &= tot > 0
    return a, keep


# pass 1: per-source seurat_v3 ranks on raw counts of common genes
ranks, filt = {}, {}
for f, s in zip(files, sources):
    a, keep = load(f)
    filt[s] = {'cells_in': int(a.n_obs), 'cells_kept': int(keep.sum()), 'removed_high_contamination': int(((a.obs['decontX_contamination'].values >= MAX_CONT)).sum())}
    b = ad.AnnData(X=a.layers['counts'][keep][:, a.var_names.get_indexer(common)].astype(np.float32), var=pd.DataFrame(index=common))
    sc.pp.highly_variable_genes(b, n_top_genes=N_CAND, flavor='seurat_v3')
    r = b.var['highly_variable_rank']
    ranks[s] = r.values.astype(float)
    log(step='hvg_source', source=s, **filt[s])
    del a, b; gc.collect()
R = pd.DataFrame(ranks, index=common)
hv = pd.DataFrame({'highly_variable_rank': np.nanmedian(R.values, axis=1), 'highly_variable_nbatches': R.notna().sum(1).values}, index=common)
hv['technical'] = [technical(g) for g in common]
cand = hv.sort_values(['highly_variable_nbatches', 'highly_variable_rank'], ascending=[False, True]).index[:N_CAND]
elig = hv.loc[cand][~hv.loc[cand, 'technical']].sort_values(['highly_variable_rank', 'highly_variable_nbatches'], ascending=[True, False], kind='mergesort')
selected = list(elig.index[:N_HVG])
hv['candidate'] = hv.index.isin(cand); hv['selected_HVG'] = hv.index.isin(selected)
hv.to_csv(f'{out}/report/hvg_v2.tsv', sep='\t')
log(step='hvg_selected', n=len(selected))

# pass 2: normalise decontX counts (all genes), write integration shards, collect HVG matrix
blocks, obs_all = [], []
for f, s in zip(files, sources):
    a, keep = load(f)
    a = a[keep].copy()
    a.layers['decontx'] = a.X.astype(np.float32)
    n = ad.AnnData(X=a.layers['decontx'].copy())
    sc.pp.normalize_total(n, target_sum=1e4); sc.pp.log1p(n)
    a.X = sp.csr_matrix(n.X, dtype=np.float32)
    a.obs['total_counts_decontx'] = np.asarray(a.layers['decontx'].sum(1)).ravel().astype(np.float32)
    blocks.append(a[:, selected].X.copy())
    a.write_h5ad(f'{out}/shards/{s}.h5ad')
    obs_all.append(a.obs[['source', 'data_batch', 'analysis_main', 'sampling_site', 'treatment_status', 'patient_id', 'unit_id']].copy())
    log(step='normalised', source=s, cells=int(a.n_obs))
    del a, n; gc.collect()
Y = sp.vstack(blocks).tocsr(); del blocks; gc.collect()
obs = pd.concat(obs_all); del obs_all
for c in obs.columns:
    obs[c] = obs[c].astype(str).astype('category')
log(step='matrix', shape=Y.shape, nnz=int(Y.nnz))

# PCA, Harmony, neighbours, UMAP
B = ad.AnnData(X=Y, obs=obs, var=pd.DataFrame(index=selected))
sc.pp.pca(B, n_comps=N_PC, zero_center=True, svd_solver='arpack', random_state=seed, dtype='float32')
np.save(f'{out}/report/X_pca.npy', B.obsm['X_pca']); log(step='pca_done')
ho = hm.run_harmony(data_mat=np.asarray(B.obsm['X_pca'], np.float64), meta_data=pd.DataFrame({'source': B.obs['source'].values}),
                    vars_use=['source'], theta=2.0, lamb=1.0, max_iter_harmony=20, random_state=seed)
Z = np.asarray(ho.Z_corr).T if ho.Z_corr.shape[0] == N_PC else np.asarray(ho.Z_corr)
B.obsm['X_pca_harmony'] = Z.astype(np.float32)
np.save(f'{out}/report/X_pca_harmony.npy', B.obsm['X_pca_harmony']); log(step='harmony_done')
B.X = sp.csr_matrix(B.shape, dtype=np.float32)  # free memory; graph uses the embedding only
sc.pp.neighbors(B, n_neighbors=15, use_rep='X_pca_harmony', random_state=seed)
log(step='neighbors_done')
sc.tl.umap(B, random_state=seed)
np.save(f'{out}/report/X_umap.npy', B.obsm['X_umap']); log(step='umap_done')

# figures
sc.settings.figdir = f'{out}/figures'
rng = np.random.default_rng(seed)
plot = B[rng.permutation(B.n_obs)]
for key in ['analysis_main', 'source', 'data_batch', 'sampling_site', 'treatment_status']:
    sc.pl.umap(plot, color=key, size=0.5, frameon=False, legend_fontsize=5, show=False, save=f'_{key}.png', title=f'UMAP ({B.n_obs:,} cells): {key}')
log(step='figures_done')

# final zarr
zpath = f'{out}/integration_v2.zarr'
concat_on_disk([f'{out}/shards/{s}.h5ad' for s in sources], zpath, join='outer', max_loaded_elems=200_000_000)
g = zarr.open(zpath, 'r+')
assert list(read_elem(g['obs']).index) == list(B.obs_names), 'cell order mismatch'
for k in ['X_pca', 'X_pca_harmony', 'X_umap']:
    write_elem(g['obsm'] if 'obsm' in g else g.create_group('obsm'), k, B.obsm[k])
var = read_elem(g['var'])
var['selected_HVG'] = var.index.isin(selected)
del g['var']; write_elem(g, 'var', var)
commit = subprocess.run(['git', '-C', repo, 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
write_elem(g, 'uns', {'schema': 'pancancer_integration_v2', 'code_commit': commit, 'seed': str(seed),
                      'X': 'log1p(normalize_total 1e4) of decontX-corrected counts', 'layers': 'counts = raw integer; decontx = decontX-corrected',
                      'hvg': f'seurat_v3 per source on raw counts of {len(common)} common genes; {N_CAND} candidates -> {N_HVG} non-technical',
                      'pca': f'{N_PC} PCs, centred, ARPACK, no scaling', 'harmony': 'batch=source, theta 2, lambda 1, max_iter 20',
                      'graph': '15-NN on X_pca_harmony; UMAP scanpy default', 'cell_filter': f'decontX done and contamination < {MAX_CONT}',
                      'source_filter': json.dumps(filt)})
zarr.consolidate_metadata(zpath)
json.dump({'cells': int(B.n_obs), 'sources': len(sources), 'common_genes': len(common), 'hvg': len(selected), 'filter': filt},
          open(f'{out}/report/integration_summary.json', 'w'), indent=1)
log(step='done', cells=int(B.n_obs))
