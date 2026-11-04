# Batch-correction benchmark v1, step 4 (scvi env): one evaluation for every method.
#  A graph metrics (all methods; embeddings -> 15-NN UMAP connectivities as in scanpy, BBKNN -> its own graph):
#    weighted same-source / same-unit neighbour fraction, graph source iLISI (inverse Simpson of weighted
#    neighbour sources), label purity for confident coarse and fine marker labels, graph connectivity per
#    coarse label (share of a label's cells in its largest connected component; scib definition)
#  B scib-metrics Benchmarker (embedding methods only) on a 100k subset of cells with a confident coarse label:
#    batch = source, label = label_coarse
# Labels are per-cell marker labels computed from expression only (prep_subsample.py), so no method is favoured
# by how labels were made. Usage: evaluate.py SUB.h5ad BENCH_DIR SEED
import sys, os, glob, json
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc
from scipy.sparse.csgraph import connected_components

subp, bdir, seed = sys.argv[1], sys.argv[2], int(sys.argv[3])
A = ad.read_h5ad(subp, backed='r'); obs = A.obs.copy(); n = len(obs)
src = pd.Categorical(obs['source']).codes; unit = pd.Categorical(obs['unit_id']).codes
lc = obs['label_coarse'].astype(str).values; lf = obs['label_fine'].astype(str).values
conf_c = lc != 'ambiguous'; conf_f = ~pd.Series(lf).str.contains('ambiguous').values
graphs = {}
for f in sorted(glob.glob(f'{bdir}/emb_*.npy')):
    k = os.path.basename(f)[4:-4]
    B = ad.AnnData(X=sp.csr_matrix((n, 1), dtype=np.float32)); B.obsm['e'] = np.load(f)
    sc.pp.neighbors(B, n_neighbors=15, use_rep='e', random_state=seed)
    graphs[k] = B.obsp['connectivities'].tocsr(); print('graph', k, flush=True)
for f in sorted(glob.glob(f'{bdir}/graph_*.npz')):
    graphs[os.path.basename(f)[6:-4]] = sp.load_npz(f).tocsr()


def gmetrics(G):
    G = G.copy(); G.setdiag(0); G.eliminate_zeros()
    W = sp.diags(1 / np.maximum(np.asarray(G.sum(1)).ravel(), 1e-12)) @ G
    rows = np.repeat(np.arange(n), np.diff(W.indptr)); w = W.data; nb = W.indices
    same = lambda codes: np.bincount(rows, weights=w * (codes[nb] == codes[rows]), minlength=n)
    P = sp.csr_matrix((w, (rows, src[nb])), shape=(n, src.max() + 1))
    ilisi = 1 / np.asarray(P.multiply(P).sum(1)).ravel()
    lcc = pd.Categorical(lc).codes; lfc = pd.Categorical(lf).codes
    gc_ = {}
    for lab in sorted(set(lc[conf_c])):
        idx = np.flatnonzero(lc == lab)
        _, comp = connected_components(G[idx][:, idx], directed=False)
        gc_[lab] = round(float(np.bincount(comp).max() / len(idx)), 3)
    ss = same(src)
    return {'same_source_median': round(float(np.median(ss)), 3), 'same_unit_median': round(float(np.median(same(unit))), 3),
            'source_iLISI_median': round(float(np.median(ilisi)), 3),
            'purity_coarse': round(float(same(lcc)[conf_c].mean()), 3), 'purity_fine': round(float(same(lfc)[conf_f].mean()), 3),
            'graph_connectivity_mean': round(float(np.mean(list(gc_.values()))), 3),
            'same_source_by_label': pd.Series(ss[conf_c]).groupby(lc[conf_c]).median().round(3).to_dict(), 'graph_connectivity': gc_}


res = {}
for k, G in graphs.items():
    res[k] = gmetrics(G); print(k, json.dumps({a: b for a, b in res[k].items() if not isinstance(b, dict)}), flush=True)
json.dump(res, open(f'{bdir}/graph_metrics.json', 'w'), indent=1)
T = pd.DataFrame({k: {a: b for a, b in v.items() if not isinstance(b, dict)} for k, v in res.items()}).T
T.sort_values('purity_fine', ascending=False).to_csv(f'{bdir}/graph_metrics.tsv', sep='\t')
print(T.to_string(), flush=True)

# scib-metrics on embeddings
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
rng = np.random.default_rng(seed)
idx = np.sort(rng.choice(np.flatnonzero(conf_c), min(100_000, conf_c.sum()), replace=False))
C = ad.AnnData(X=sp.csr_matrix((len(idx), 1), dtype=np.float32), obs=obs.iloc[idx][['source', 'label_coarse']].copy())
C.obs['label_coarse'] = C.obs['label_coarse'].astype(str).astype('category')
keys = []
for f in sorted(glob.glob(f'{bdir}/emb_*.npy')):
    k = os.path.basename(f)[4:-4]; C.obsm[k] = np.load(f)[idx]; keys.append(k)
bm = Benchmarker(C, batch_key='source', label_key='label_coarse', embedding_obsm_keys=keys, pre_integrated_embedding_obsm_key='unintegrated',
                 bio_conservation_metrics=BioConservation(), batch_correction_metrics=BatchCorrection(), n_jobs=8)
bm.benchmark()
R = bm.get_results(min_max_scale=False)
R.to_csv(f'{bdir}/scib_metrics.tsv', sep='\t'); print(R.to_string(), flush=True)
