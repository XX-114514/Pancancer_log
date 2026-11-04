# Batch-correction benchmark v1, step 3 (GPU, scvi env): scVI on raw counts of the 4,000 TME HVGs.
#   scvi_source   batch_key = source
#   scvi_unit     batch_key = unit_id (sample)
# n_latent 30, n_layers 2, gene_likelihood nb, max 400 epochs with early stopping (scvi-tools defaults otherwise)
# Usage: run_scvi.py SUB.h5ad OUT_DIR SEED
import sys, os, time, json
import numpy as np, anndata as ad, scvi, torch

subp, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
scvi.settings.seed = seed
torch.set_float32_matmul_precision('high')
A = ad.read_h5ad(subp)
A.X = A.X.astype(np.float32)
times = {}
for name, key in (('scvi_source', 'source'), ('scvi_unit', 'unit_id')):
    t = time.time()
    B = A.copy()
    scvi.model.SCVI.setup_anndata(B, batch_key=key)
    m = scvi.model.SCVI(B, n_latent=30, n_layers=2, gene_likelihood='nb')
    m.train(max_epochs=400, early_stopping=True, batch_size=1024, accelerator='gpu')
    np.save(f'{out}/emb_{name}.npy', m.get_latent_representation().astype(np.float32))
    m.save(f'{out}/model_{name}', overwrite=True)
    times[name] = round(time.time() - t, 1); print(name, times[name], 's', flush=True)
json.dump(times, open(f'{out}/times_scvi.json', 'w'), indent=1)
