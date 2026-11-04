# Batch-correction benchmark v1 (2026-11-04), step 1: one shared 300k-cell subsample of TME v1 for all
# methods. Writes SUB.h5ad with
#   X = raw counts of the 4,000 TME v1 HVGs (scVI input); layers['lognorm'] = log1p(1e4) of decontX counts
#   obsm X_pca (full-data TME PCA, subset) and X_v1_harmony (full-data v1 Harmony, subset)
#   obs source, unit_id, analysis_main, leiden_tme, label_coarse / label_fine (per-cell marker labels from
#   expression alone, independent of any embedding; 'ambiguous' when no clear winner)
# Usage: prep_subsample.py TME_RUN_DIR OUT.h5ad N SEED
import sys, json
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, zarr
from anndata.io import read_elem

run, outp, N, seed = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
g = zarr.open(f'{run}/tme_v1.zarr', 'r')
obs = read_elem(g['obs'])[['source', 'unit_id', 'patient_id', 'analysis_main', 'leiden_tme', 'compartment', 'tme_flag']]
hv = pd.read_csv(f'{run}/report/hvg_tme_v1.tsv', sep='\t', index_col=0)
hvg = list(hv.index[hv['selected_HVG']])
COARSE = {'T_NK': ['CD3D', 'CD3E', 'NKG7', 'GNLY'], 'B': ['MS4A1', 'CD79A'], 'plasma': ['JCHAIN', 'MZB1'],
          'myeloid': ['CD68', 'LYZ', 'CD14', 'C1QA'], 'mast': ['TPSAB1', 'CPA3'], 'fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'],
          'endothelial': ['PECAM1', 'VWF', 'CDH5'], 'pericyte_SMC': ['RGS5', 'ACTA2', 'MYH11'], 'epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19']}
FINE = {'T_NK': {'CD8T': ['CD8A', 'CD8B'], 'CD4T': ['CD4', 'IL7R'], 'Treg': ['FOXP3', 'IL2RA'], 'NK': ['KLRF1', 'FCGR3A']},
        'myeloid': {'Mono_classical': ['CD14', 'S100A8', 'VCAN'], 'Macro': ['C1QA', 'C1QC', 'APOE'], 'cDC1': ['CLEC9A', 'XCR1'],
                    'cDC2': ['CD1C', 'FCER1A'], 'Neutrophil': ['CSF3R', 'FCGR3B']},
        'fibroblast': {'myCAF': ['ACTA2', 'TAGLN', 'POSTN'], 'iCAF': ['CFD', 'C3', 'PI16']}}
mk = sorted({x for v in COARSE.values() for x in v} | {x for d in FINE.values() for v in d.values() for x in v})
rng = np.random.default_rng(seed)
sub = np.sort(rng.choice(len(obs), N, replace=False))
sources = list(dict.fromkeys(obs['source'].astype(str)))
src_arr = obs['source'].astype(str).values
offs = {s: np.flatnonzero(src_arr == s) for s in sources}
C, L, M = [], [], []
for s in sources:
    rows = offs[s]; w = np.flatnonzero(np.isin(rows, sub))
    if not len(w):
        continue
    assert np.all(np.diff(rows) == 1)
    a = ad.read_h5ad(f'{run}/shards/{s}.h5ad', backed='r')
    b = a[w].to_memory()
    gi = b.var_names.get_indexer(hvg); assert (gi >= 0).all(), s
    C.append(b.layers['counts'][:, gi].astype(np.float32)); L.append(sp.csr_matrix(b.X[:, gi], dtype=np.float32))
    m = pd.DataFrame(0.0, index=range(b.n_obs), columns=mk)
    have = [x for x in mk if x in b.var_names]
    m[have] = np.asarray(b[:, have].X.todense())
    M.append(m)
    print(s, len(w), flush=True)
M = pd.concat(M, ignore_index=True)
S = pd.DataFrame({k: M[v].mean(1) for k, v in COARSE.items()})
srt = np.sort(S.values, 1)
coarse = np.where((srt[:, -1] >= 1.0) & (srt[:, -1] - srt[:, -2] >= 0.5), S.columns[S.values.argmax(1)], 'ambiguous')
fine = coarse.astype(object).copy()
for comp, subs in FINE.items():
    w = coarse == comp
    F = pd.DataFrame({k: M.loc[w, v].mean(1) for k, v in subs.items()})
    fs = np.sort(F.values, 1)
    fine[w] = np.where((fs[:, -1] >= 0.75) & (fs[:, -1] - fs[:, -2] >= 0.4), comp + ':' + F.columns[F.values.argmax(1)], comp + ':ambiguous')
o = obs.iloc[sub].copy()
o['label_coarse'] = pd.Categorical(coarse); o['label_fine'] = pd.Categorical(fine)
for c in o.columns:
    if o[c].dtype == object:
        o[c] = o[c].astype(str).astype('category')
A = ad.AnnData(X=sp.vstack(C).tocsr(), obs=o, var=pd.DataFrame(index=hvg))
A.layers['lognorm'] = sp.vstack(L).tocsr()
A.obsm['X_pca'] = g['obsm']['X_pca'][:][sub].astype(np.float32)
A.obsm['X_v1_harmony'] = g['obsm']['X_pca_harmony'][:][sub].astype(np.float32)
A.write_h5ad(outp)
print(json.dumps({'cells': A.n_obs, 'genes': A.n_vars, 'sources': int(o.source.nunique()), 'units': int(o.unit_id.nunique()),
                  'coarse': pd.Series(coarse).value_counts().to_dict(), 'fine': pd.Series(fine).value_counts().to_dict()}, default=int))
