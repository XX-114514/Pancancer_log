# Annotation v2 (2026-11-07): checks on the lineage objects, no labels are changed.
#  1. LISI (Korsunsky 2019, perplexity 30, 90-NN on a 200,000-cell random subset) on source and sample, for X_pca and X_emb.
#  2. Harmony without early stopping (epsilon -inf, 20 iterations; methods file B5) against the run's embedding: ARI at resolution 1, LISI.
#  3. Detection of rare-subtype and identity genes at resolutions 1.5, 2, 3, with the crosswalk to the main resolution.
# Usage: lineage_batch_and_rare_checks_v1.py OBJ_ROOT OUT_DIR SEED POOL [POOL...]   (env HARMONY_POOLS="Endothelial,Myeloid")
import sys, os, json
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, harmonypy as hm, pynndescent
from sklearn.metrics import adjusted_rand_score

root, out, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]); pools = sys.argv[4:]
os.makedirs(out, exist_ok=True)
HP = set(x for x in os.environ.get('HARMONY_POOLS', '').split(',') if x)
GENES = {
    'T_NK': 'CD3D CD3E CD4 CD8A CD8B CCR7 SELL LEF1 TCF7 IL7R FOXP3 IL2RA CTLA4 GZMK GZMB GZMH GNLY PRF1 NKG7 FGFBP2 FCGR3A KLRF1 KLRD1 NCAM1 XCL1 CX3CR1 PDCD1 HAVCR2 LAG3 TOX CXCL13 CXCR5 BCL6 ICOS TRDC TRGC1 TRGC2 TRDV1 TRDV2 TRGV9 SLC4A10 KLRB1 ZBTB16 RORC CCR6 IL17A IL23R IFNG TBX21 IL4I1 KIT MKI67'.split(),
    'B_Plasma': 'MS4A1 CD79A CD19 CD27 IGHD IGHM TCL1A FCER2 TNFRSF13B ITGAX TBX21 FCRL4 FCRL5 RGS13 BCL6 AICDA MKI67 MZB1 SDC1 XBP1 PRDM1 IGHG1 IGHG3 IGHG4 IGHA1 IGHA2 JCHAIN IGHE'.split(),
    'Stromal': 'COL1A1 DCN LUM PDGFRA FAP POSTN COL10A1 CTHRC1 MMP11 ACTA2 TAGLN MYH11 ACTG2 DES RGS5 NOTCH3 PDGFRB KCNJ8 PI16 CD34 CFD APOD C3 CXCL12 CXCL14 CCL19 CCL21 MMP1 IL6 HLA-DRA CD74 SFRP2 SFRP4 COMP LRRC15 GREM1 MSLN KRT18 ESR1 MKI67'.split(),
}
def lisi(Z, labels, nn):
    y = pd.factorize(labels)[0]; n = y.max() + 1; D, I = nn
    D = D[:, 1:] ** 2; I = I[:, 1:]; out_ = np.zeros(len(Z))
    for i in range(len(Z)):
        lo, hi, beta = 0., np.inf, 1.; d = D[i]
        for _ in range(50):
            p = np.exp(-d * beta); s = p.sum()
            if s <= 0:
                p = np.ones_like(d); s = p.sum(); break
            p /= s; h = -(p[p > 0] * np.log(p[p > 0])).sum() - np.log(30)
            if abs(h) < 1e-5:
                break
            if h > 0:
                lo = beta; beta = beta * 2 if hi == np.inf else (beta + hi) / 2
            else:
                hi = beta; beta = (beta + lo) / 2
        p = p / p.sum(); w = np.bincount(y[I[i]], weights=p, minlength=n)
        out_[i] = 1 / (w ** 2).sum()
    return out_
def lisi_row(pool, rep, Z, ob, rng):
    idx = rng.choice(len(Z), min(200_000, len(Z)), replace=False); z = np.asarray(Z[idx], np.float32)
    nn = pynndescent.NNDescent(z, n_neighbors=91, random_state=seed, n_jobs=-1).neighbor_graph
    nn = (nn[1], nn[0])
    r = {'pool': pool, 'rep': rep}
    for k in ('source', 'sample_id'):
        v = lisi(z, ob[k].astype(str).values[idx], nn); r[f'{k}_lisi_median'] = round(float(np.median(v)), 2); r[f'{k}_lisi_mean'] = round(float(v.mean()), 2); r[f'{k}_n'] = int(ob[k].astype(str).iloc[idx].nunique())
    return r
rows = []; ev = []; cross = []
for pool in pools:
    A = ad.read_h5ad(f'{root}/{pool}/object.h5ad'); rng = np.random.default_rng(seed)
    print(json.dumps({'pool': pool, 'shape': A.shape}), flush=True)
    rows.append(lisi_row(pool, 'X_pca', A.obsm['X_pca'], A.obs, rng)); rows.append(lisi_row(pool, 'X_emb', A.obsm['X_emb'], A.obs, rng))
    if pool in HP:
        h = hm.run_harmony(np.asarray(A.obsm['X_pca'], np.float64), A.obs[['source']].astype(str), vars_use=['source'], theta=2.0, lamb=1.0, max_iter_harmony=20, epsilon_cluster=-np.inf, epsilon_harmony=-np.inf, random_state=seed, verbose=False)
        Z = np.asarray(h.Z_corr); Z = (Z.T if Z.shape[0] == A.obsm['X_pca'].shape[1] else Z).astype(np.float32)
        r = lisi_row(pool, 'X_emb_no_early_stop', Z, A.obs, rng)
        B = ad.AnnData(np.zeros((A.n_obs, 1), np.float32), obs=A.obs[['leiden_r1']].copy()); B.obsm['Z'] = Z
        sc.pp.neighbors(B, n_neighbors=15, use_rep='Z', random_state=seed); sc.tl.leiden(B, resolution=1.0, random_state=seed, flavor='igraph', n_iterations=2, directed=False)
        r['ARI_r1_vs_run'] = round(adjusted_rand_score(B.obs['leiden_r1'], B.obs['leiden']), 3); r['clusters_r1'] = int(B.obs['leiden'].nunique()); rows.append(r); del B, Z
    pd.DataFrame(rows).to_csv(f'{out}/lisi_harmony.tsv', sep='\t', index=False)
    if pool in GENES:
        g = [x for x in GENES[pool] if x in A.var_names]; X = (A[:, g].layers['counts'] > 0).astype(np.float32).toarray()
        for r_ in ('1', '1.5', '2', '3'):
            c = A.obs[f'leiden_r{r_}'].astype(str).values; t = pd.DataFrame(X, columns=g).groupby(c).mean().round(3)
            t.insert(0, 'cells', pd.Series(c).value_counts()); t.insert(0, 'resolution', r_); ev.append(t.rename_axis('cluster').reset_index().assign(pool=pool))
            if r_ != '1':
                x = pd.crosstab(c, A.obs['leiden_r1'].astype(str).values); m = x.idxmax(1); cross.append(pd.DataFrame({'pool': pool, 'resolution': r_, 'cluster': x.index, 'main_cluster': m.values, 'share_in_main': (x.max(1) / x.sum(1)).round(3).values}))
        pd.concat(ev).to_csv(f'{out}/rare_gene_detection.tsv', sep='\t', index=False); pd.concat(cross).to_csv(f'{out}/resolution_crosswalk.tsv', sep='\t', index=False)
    del A
print('done', flush=True)
