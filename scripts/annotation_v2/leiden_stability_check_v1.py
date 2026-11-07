# Annotation v2 (2026-11-07): is the low ARI between Harmony with and without early stopping (lisi_harmony.tsv) an effect of
# Harmony, or ordinary Leiden variability at resolution 1? Baseline: the run's own embedding, neighbours recomputed, Leiden with
# other seeds. Compared at cluster level (ARI) and at reviewed-label level (each new cluster takes the majority reviewed label).
# Usage: leiden_stability_check_v1.py OBJ_ROOT REVIEW_TSV OUT_TSV POOL [POOL...]
import sys
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, harmonypy as hm
from sklearn.metrics import adjusted_rand_score as ari, normalized_mutual_info_score as nmi

root, rev, out = sys.argv[1:4]; pools = sys.argv[4:]
R = pd.read_csv(rev, sep='\t', dtype={'cluster': str}); rows = []
def clus(Z, seed):
    B = ad.AnnData(np.zeros((len(Z), 1), np.float32)); B.obsm['Z'] = Z
    sc.pp.neighbors(B, n_neighbors=15, use_rep='Z', random_state=seed); sc.tl.leiden(B, resolution=1.0, random_state=seed, flavor='igraph', n_iterations=2, directed=False)
    return B.obs['leiden'].astype(str).values
for pool in pools:
    A = ad.read_h5ad(f'{root}/{pool}/object.h5ad', backed='r'); ob = A.obs; E = np.asarray(A.obsm['X_emb'], np.float32); P = np.asarray(A.obsm['X_pca'], np.float64)
    run = ob['leiden_r1'].astype(str).values; lab = pd.Series(run).map(R[R.pool == pool].set_index('cluster')['r2_label']).values
    h = hm.run_harmony(P, ob[['source']].astype(str), vars_use=['source'], theta=2.0, lamb=1.0, max_iter_harmony=20, epsilon_cluster=-np.inf, epsilon_harmony=-np.inf, random_state=2026110701, verbose=False)
    Z = np.asarray(h.Z_corr); Z = (Z.T if Z.shape[0] == P.shape[1] else Z).astype(np.float32)
    for name, emb, seed in [('same_embedding_same_seed', E, 2026110701), ('same_embedding_seed_1', E, 1), ('same_embedding_seed_2', E, 2), ('no_early_stop_same_seed', Z, 2026110701), ('no_early_stop_seed_1', Z, 1)]:
        c = clus(emb, seed); maj = pd.crosstab(c, lab).idxmax(1); ml = pd.Series(c).map(maj).values
        rows.append({'pool': pool, 'variant': name, 'clusters': len(set(c)), 'ARI_clusters': round(ari(run, c), 3), 'NMI_clusters': round(nmi(run, c), 3), 'ARI_labels': round(ari(lab, ml), 3), 'label_agreement': round(float((lab == ml).mean()), 3)})
        print(rows[-1], flush=True); pd.DataFrame(rows).to_csv(out, sep='\t', index=False)
