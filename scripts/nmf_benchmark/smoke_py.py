"""Smoke test: scikit-learn NMF and cNMF on one simulated Poisson count matrix (cells x genes, K=4)."""
import sys, time, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
from sklearn.decomposition import NMF
from cnmf import cNMF
out = sys.argv[1]
rng = np.random.default_rng(1); G, N, K = 500, 400, 4
Y = rng.poisson(rng.gamma(0.5, 1, (N, K)) @ rng.gamma(0.3, 1, (K, G)) * 3)
adata = ad.AnnData(sp.csr_matrix(Y.astype(np.float32)), obs=pd.DataFrame(index=[f"c{i}" for i in range(N)]),
                   var=pd.DataFrame(index=[f"g{i}" for i in range(G)]))
adata = adata[:, np.asarray(Y.sum(0)).ravel() > 0].copy()
t0 = time.time()
for solver, loss in [("cd", "frobenius"), ("mu", "kullback-leibler")]:
    m = NMF(K, init="nndsvda", solver=solver, beta_loss=loss, max_iter=500, random_state=1)
    W = m.fit_transform(np.log1p(adata.X.toarray()))
    print("sklearn", solver, loss, W.shape, m.reconstruction_err_)
counts = f"{out}/sim_counts.h5ad"; adata.write_h5ad(counts)
c = cNMF(output_dir=out, name="smoke")
c.prepare(counts_fn=counts, components=[3, 4], n_iter=10, seed=1, num_highvar_genes=300)
c.factorize(worker_i=0, total_workers=1)
c.combine()
c.k_selection_plot()
c.consensus(k=4, density_threshold=0.5)
usage, spectra_scores, spectra_tpm, top_genes = c.load_results(K=4, density_threshold=0.5)
print("cNMF usage", usage.shape, "spectra", spectra_tpm.shape, f"total {time.time() - t0:.1f}s")
