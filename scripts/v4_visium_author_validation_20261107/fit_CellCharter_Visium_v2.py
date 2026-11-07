"""Per-section expression/position niches, independent of V4 CM scores and labels."""
import sys, json, hashlib
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
import squidpy as sq
import scvi
import torch
import cellcharter as cc
from importlib.metadata import version

r = Path(sys.argv[1])
key = sys.argv[2]
start = datetime.now(timezone.utc).isoformat()
torch.set_num_threads(4)
seed = 2026110701
scvi.settings.seed = seed
np.random.seed(seed)
torch.manual_seed(seed)
job = next(j for j in json.loads((r / 'inputs/section_jobs_PRIVATE.json').read_text()) if j['section_code'] == key)
out = r / 'niches_retry2' / key
out.mkdir(parents=True, exist_ok=False)
source = ad.read_h5ad(job['input_path'], backed='r')
mask = source.obs[job['group_column']].astype(str).eq(job['group_value']).to_numpy() if job['group_column'] else np.ones(source.n_obs, bool)
if 'in_tissue' in source.obs:
    mask &= pd.to_numeric(source.obs.in_tissue, errors='coerce').fillna(0).eq(1).to_numpy()
a = source[np.flatnonzero(mask), :].to_memory()
source.file.close()
counts = sp.csr_matrix(a.layers[job['counts_layer']] if job['counts_layer'] != 'X' else a.X, dtype=np.float32)
assert np.isfinite(counts.data).all() and (counts.data >= 0).all() and np.allclose(counts.data, np.rint(counts.data))
if job['coordinates'] == 'obsm/spatial':
    xy = np.asarray(a.obsm['spatial'])[:, :2].astype(float)
elif job['coordinates'].endswith('pxl_row_in_fullres'):
    xy = a.obs[['pxl_col_in_fullres', 'pxl_row_in_fullres']].astype(float).to_numpy()
else:
    xy = a.obs[['pxl_col', 'pxl_row']].astype(float).to_numpy()
# Match the model's prespecified spot eligibility without using fitted abundance.
signature_genes = set(pd.read_csv(r / 'reference/model_genes.tsv', sep='\t').gene)
shared = a.var_names.isin(signature_genes)
keep = (np.asarray(counts[:, shared].sum(axis=1)).ravel() > 0) & np.isfinite(xy).all(axis=1)
counts = counts[keep]
xy = xy[keep]
assert len(xy) > 60 and len(np.unique(xy, axis=0)) == len(xy)
b = ad.AnnData(counts, obs=pd.DataFrame(index=a.obs_names[keep]), var=pd.DataFrame(index=a.var_names))
b.obsm['spatial'] = xy
genes = (~b.var_names.str.upper().str.startswith('MT-')) & (np.asarray(b.X.sum(axis=0)).ravel() >= 3)
b = b[:, genes].copy()
sc.pp.highly_variable_genes(b, flavor='seurat_v3', n_top_genes=min(2000, b.n_vars))
b = b[:, b.var.highly_variable].copy()
n_before_niche_information_filter = b.n_obs
niche_eligible = np.asarray(b.X.sum(axis=1)).ravel() > 0
pd.DataFrame({'niche_eligible_nonzero_HVG': niche_eligible}, index=b.obs_names).to_csv(out / 'spot_niche_eligibility_PRIVATE.tsv.gz', sep='\t')
b = b[niche_eligible].copy()
xy = xy[niche_eligible]
assert b.n_obs > 60

scvi.model.SCVI.setup_anndata(b)
model = scvi.model.SCVI(b, n_latent=10)
model.train(max_epochs=400, early_stopping=True, train_size=.9, validation_size=.1,
            batch_size=128, accelerator='cpu', enable_progress_bar=False)
latent = model.get_latent_representation(b).astype(np.float32)
assert latent.shape == (b.n_obs, 10) and np.isfinite(latent).all()
b.obsm['X_scVI'] = latent
# Six-neighbour Visium niche graph is distinct from the author's 60NN Moran graph.
sq.gr.spatial_neighbors(b, coord_type='generic', n_neighs=6, delaunay=False, set_diag=False)
cc.gr.aggregate_neighbors(b, n_layers=3, use_rep='X_scVI', out_key='X_cellcharter')
features = np.asarray(b.obsm['X_cellcharter'])
assert features.shape == (b.n_obs, 40) and np.isfinite(features).all()
autok = cc.tl.ClusterAutoK(n_clusters=(2, 10), max_runs=10, convergence_tol=.001,
    model_params={'random_state': seed, 'trainer_params': {'accelerator': 'cpu', 'enable_progress_bar': False}})
autok.fit(b, use_rep='X_cellcharter')
labels = np.asarray(autok.predict(b, use_rep='X_cellcharter')).astype(str)
assert len(labels) == b.n_obs and len(np.unique(labels)) >= 2
stability = [{'k': int(k), 'stability_mean': float(np.mean(v)), 'n_pair_scores': len(v)}
             for k, v in zip(autok.n_clusters[1:-1], autok.stability)]
best = int(autok.best_k)
assert 2 <= best <= 10
pd.DataFrame(stability).to_csv(out / 'K_stability_PRIVATE.tsv', sep='\t', index=False)
pd.DataFrame({'niche': labels, 'x': xy[:, 0], 'y': xy[:, 1]}, index=b.obs_names).to_csv(out / 'spot_niches_PRIVATE.tsv.gz', sep='\t')
np.save(out / 'scVI_latent_PRIVATE.npy', latent)
np.save(out / 'CellCharter_features_PRIVATE.npy', features)
sp.save_npz(out / 'spatial_6NN_PRIVATE.npz', b.obsp['spatial_connectivities'])
np.save(out / 'scVI_ELBO_history_PRIVATE.npy', np.asarray(model.history['elbo_train'], dtype=float))
pd.DataFrame({'gene': b.var_names}).to_csv(out / 'selected_genes.tsv', sep='\t', index=False)
model.save(str(out / 'scVI_model_PRIVATE'), overwrite=False)
autok.save(out / 'CellCharter_model_PRIVATE')
receipt = {
    'status': 'EXPRESSION_POSITION_ONLY_CELLCHARTER_NICHES_FITTED_PENDING_INDEPENDENT_AUDIT',
    'UTC_start': start, 'UTC_end': datetime.now(timezone.utc).isoformat(),
    'section_code': key, 'dataset': job['dataset'], 'n_spots': b.n_obs, 'n_HVG': b.n_vars,
    'n_spots_before_niche_information_filter': n_before_niche_information_filter,
    'n_spots_zero_in_selected_HVG_unassessable_for_niche': int((~niche_eligible).sum()),
    'scVI_max_epochs': 400, 'scVI_early_stopping': True,
    'scVI_epochs_recorded': len(model.history['elbo_train']), 'scVI_latent_dim': 10,
    'niche_graph_k': 6, 'n_layers': 3, 'feature_dim': 40,
    'eligible_K': [2, 10], 'helper_K': list(map(int, autok.n_clusters)),
    'max_runs': 10, 'convergence_tol': .001, 'best_k': best,
    'n_observed_niches': len(np.unique(labels)), 'seed': seed,
    'uses_CM_labels_abundance_activity_or_clinical_outcomes_to_fit': False,
    'CoVarNet_paper_exact_niche_parameters_claim': False,
    'spot_labels_SHA256': hashlib.sha256((out / 'spot_niches_PRIVATE.tsv.gz').read_bytes()).hexdigest(),
    'versions': {n: version(n) for n in ['cellcharter', 'scvi-tools', 'squidpy', 'anndata', 'torch']},
}
(out / 'complete.json').open('x').write(json.dumps(receipt, indent=2))
print(json.dumps(receipt))
