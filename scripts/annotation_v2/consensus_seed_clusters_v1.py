# Annotation v2 (2026-11-07, user decision): multi-seed consensus, stage 1. Leiden at resolution 1 changes about 30% of cell labels
# when only the seed changes (checks/leiden_stability.tsv), so the clustering is repeated with N seeds on the run's embedding
# (X_emb): each seed recomputes the 15-NN graph and Leiden (igraph, 2 iterations), as in major_by_cancer_v1.py. Only cluster ids
# are saved; labels are mapped in stage 2 (consensus_labels_v1.py) after the clusters are reviewed.
# Usage: consensus_seed_clusters_v1.py OBJ_H5AD OUT_TSV_GZ N_SEEDS BASE_SEED
import sys
import numpy as np, pandas as pd, anndata as ad, scanpy as sc, h5py
from anndata.io import read_elem

obj, out, n, base = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
with h5py.File(obj, 'r') as f:
    ob = read_elem(f['obs']); Z = f['obsm/X_emb'][:].astype(np.float32)
D = pd.DataFrame({'run_r1': ob['leiden_r1'].astype(str).values}, index=ob.index.astype(str))
B = ad.AnnData(np.zeros((len(Z), 1), np.float32)); B.obsm['Z'] = Z
for i in range(n):
    s = base + i
    sc.pp.neighbors(B, n_neighbors=15, use_rep='Z', random_state=s)
    sc.tl.leiden(B, resolution=1.0, key_added='l', random_state=s, flavor='igraph', n_iterations=2, directed=False)
    D[f'seed_{s}'] = B.obs['l'].astype(str).values
    print(f'seed {s}: {D[f"seed_{s}"].nunique()} clusters', flush=True)
    D.to_csv(out, sep='\t', compression='gzip', index_label='cell_id')
print('done', flush=True)
