# TME v1 batch-correction diagnostic (2026-11-04). The v1 run (Harmony on source theta 2 + unit theta 1)
# lowered source mixing (source iLISI 4.26 -> 3.02, same-source neighbour fraction 0.35 -> 0.48, k=90,
# 200k cells) while improving unit mixing; facets show source-specific strips. This script:
#  1 per Leiden cluster: number of sources, top-source share, median depth/genes/MT/contamination, uncertain share
#  2 on a fixed 300k subsample, compares embeddings: PCA (none); Harmony source theta 2 (= v2 recipe);
#    Harmony unit theta 2; Harmony source 2 + unit 1 (v1 recipe, rerun on the subsample); and the full-data
#    v1 embedding subset. Metrics on a 30-NN graph: source iLISI, same-source / same-unit neighbour
#    fractions, and biology conservation against a per-cell marker label computed from expression
#    (independent of any embedding): label purity and cLISI.
# Usage: benchmark_tme_batch_v1.py RUN_DIR SEED
import sys, os, glob, json
import numpy as np, pandas as pd, anndata as ad, zarr, harmonypy as hm
from anndata.io import read_elem
from sklearn.neighbors import NearestNeighbors

run, seed = sys.argv[1], int(sys.argv[2])
out = f'{run}/benchmark'; os.makedirs(out, exist_ok=True)
g = zarr.open(f'{run}/tme_v1.zarr', 'r')
obs = read_elem(g['obs'])
keep = ['source', 'unit_id', 'leiden_tme', 'compartment', 'tme_flag', 'total_counts', 'n_genes', 'pct_mt', 'decontX_contamination']
obs = obs[keep]
# 1 cluster diagnostics
grp = obs.groupby('leiden_tme', observed=True)
top = grp['source'].agg(lambda x: x.value_counts(normalize=True).iloc[0])
nsrc = grp['source'].agg(lambda x: int((x.value_counts(normalize=True) > 0.01).sum()))
cd = pd.DataFrame({'cells': grp.size(), 'compartment': grp['compartment'].first(), 'top_source': grp['source'].agg(lambda x: x.value_counts().index[0]),
                   'top_source_share': top.round(3), 'sources_gt1pct': nsrc,
                   'median_umi': grp['total_counts'].median(), 'median_genes': grp['n_genes'].median(),
                   'median_pct_mt': grp['pct_mt'].median().round(2), 'median_contamination': grp['decontX_contamination'].median().round(3),
                   'uncertain_share': grp['tme_flag'].agg(lambda x: (x.astype(str) != 'non_malignant').mean()).round(3)})
cd.sort_values('top_source_share', ascending=False).to_csv(f'{out}/cluster_diagnostics.tsv', sep='\t')
print(cd.sort_values('top_source_share', ascending=False).head(25).to_string(), flush=True)

# 2 subsample + per-cell marker label from the TME shards (same order as the zarr)
rng = np.random.default_rng(seed)
N = 300_000
sub = np.sort(rng.choice(len(obs), N, replace=False))
MK = {'T_NK': ['CD3D', 'CD3E', 'NKG7', 'GNLY'], 'B': ['MS4A1', 'CD79A'], 'plasma': ['JCHAIN', 'MZB1'], 'myeloid': ['CD68', 'LYZ', 'CD14', 'C1QA'],
      'mast': ['TPSAB1', 'CPA3'], 'fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'], 'endothelial': ['PECAM1', 'VWF', 'CDH5'],
      'pericyte_SMC': ['RGS5', 'ACTA2', 'MYH11'], 'epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19']}
genes = sorted({x for v in MK.values() for x in v})
sources = list(dict.fromkeys(obs['source'].astype(str)))
offs = np.cumsum([0] + [int((obs['source'].astype(str) == s).sum()) for s in sources])
rows = []
for i, s in enumerate(sources):
    w = sub[(sub >= offs[i]) & (sub < offs[i + 1])] - offs[i]
    a = ad.read_h5ad(f'{run}/shards/{s}.h5ad', backed='r')
    assert a.n_obs == offs[i + 1] - offs[i]
    X = a[w][:, genes].to_memory().X
    rows.append(pd.DataFrame(np.asarray(X.todense()), columns=genes))
E = pd.concat(rows, ignore_index=True)
S = pd.DataFrame({k: E[v].mean(1) for k, v in MK.items()})
srt = np.sort(S.values, 1)
lab = np.where((srt[:, -1] >= 1.0) & (srt[:, -1] - srt[:, -2] >= 0.5), S.columns[S.values.argmax(1)], 'ambiguous')
o = obs.iloc[sub].reset_index(drop=True)
o['marker_cell'] = lab
print(json.dumps(pd.Series(lab).value_counts().to_dict()), flush=True)

pca = g['obsm']['X_pca'][:][sub].astype(np.float64)
emb = {'none_pca': pca, 'v1_full_data_subset': g['obsm']['X_pca_harmony'][:][sub]}


def harmony(vars_, theta):
    h = hm.run_harmony(pca, o[vars_].astype(str), vars_use=vars_, theta=theta, lamb=1.0, max_iter_harmony=20, random_state=seed)
    Z = np.asarray(h.Z_corr)
    return Z.T if Z.shape[0] == pca.shape[1] else Z


emb['harmony_source_t2'] = harmony(['source'], 2.0)
emb['harmony_unit_t2'] = harmony(['unit_id'], 2.0)
emb['harmony_source_t2_unit_t1'] = harmony(['source', 'unit_id'], [2.0, 1.0])
emb['harmony_source_t4'] = harmony(['source'], 4.0)

conf = lab != 'ambiguous'
src = pd.Categorical(o['source']).codes; unit = pd.Categorical(o['unit_id']).codes; lc = pd.Categorical(lab).codes
res = {}
for k, Z in emb.items():
    idx = NearestNeighbors(n_neighbors=31, n_jobs=-1).fit(Z).kneighbors(Z, return_distance=False)[:, 1:]
    Ls = src[idx]
    ilisi = 1 / np.array([np.sum((np.bincount(r) / 30) ** 2) for r in Ls])
    Ll = lc[idx]
    clisi = 1 / np.array([np.sum((np.bincount(r) / 30) ** 2) for r in Ll[conf]])
    res[k] = {'source_iLISI_median': round(float(np.median(ilisi)), 3), 'same_source_frac_median': round(float(np.median((Ls == src[:, None]).mean(1))), 3),
              'same_unit_frac_median': round(float(np.median((unit[idx] == unit[:, None]).mean(1))), 3),
              'label_purity_mean': round(float((Ll[conf] == lc[conf, None]).mean()), 3), 'label_cLISI_median': round(float(np.median(clisi)), 3),
              'same_source_by_label': pd.Series((Ls == src[:, None]).mean(1)).groupby(lab).median().round(3).to_dict()}
    print(k, json.dumps({kk: v for kk, v in res[k].items() if kk != 'same_source_by_label'}), flush=True)
json.dump(res, open(f'{out}/benchmark_metrics.json', 'w'), indent=1)
print('note: higher source iLISI / lower same-source = better mixing; higher label purity / cLISI near 1 = biology kept', flush=True)
