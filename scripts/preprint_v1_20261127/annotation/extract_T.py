"""Read-only extraction of T cells from immune atlas v1 for fine re-annotation.
Cells: major=='T', not excluded/doublet. Discovery = 4 CM02 sources (CM-eligible);
replication = GSE169246 (CD45-enriched) and GSE176021 (lymphocyte-enriched, no ex vivo).
Per-unit (sampleID) random cap keeps within-unit proportions unbiased.
"""
import os
import zarr, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp, numcodecs, sys, json, time
numcodecs.blosc.set_nthreads(8)
Z = os.environ['PROJECT_V5_ROOT'] + '/runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'
OUT = sys.argv[1]; CAP = int(sys.argv[2])
g = zarr.open(Z, mode='r')

def col(name):
    a = g['obs'][name]
    if hasattr(a, 'keys') and 'categories' in a:
        cats = np.array(a['categories'][:], dtype=object); codes = a['codes'][:]
        out = cats[np.clip(codes, 0, None)]; out[codes < 0] = None; return out
    return a[:]

cols = ['dataset', 'major', 'subcluster', 'capture_scope', 'tissue', 'timepoint', 'sampleID', 'patient_id',
        'patient_block', 'patient_identity_verified', 'CM_composition_eligible', 'annotation_exclude_from_primary',
        'residual_doublet_flag', 'treatment', 'source_response', 'native_response_endpoint', 'ex_vivo_stimulation']
obs = pd.DataFrame({c: col(c) for c in cols})
obs['row'] = np.arange(len(obs))
s = lambda c: obs[c].astype(str)
base = (s('major') == 'T') & ~s('annotation_exclude_from_primary').str.lower().isin(['true', '1']) \
       & ~s('residual_doublet_flag').str.lower().isin(['true', '1'])
disc = base & s('dataset').isin(['GSE294109', 'GSE256326', 'GSE305078', 'GSE286827']) & (s('CM_composition_eligible') == 'True')
rep = base & ((s('dataset') == 'GSE169246') | ((s('dataset') == 'GSE176021') & (s('capture_scope') == 'lymphocyte_enriched')))
sel = obs[disc | rep].copy()
sel['cohort'] = np.where(disc[sel.index], 'discovery', 'replication')
rng = np.random.default_rng(20261125)
keep = []
for u, idx in sel.groupby('sampleID').indices.items():
    idx = np.asarray(idx)
    keep.append(idx if len(idx) <= CAP else rng.choice(idx, CAP, replace=False))
sel = sel.iloc[np.sort(np.concatenate(keep))]
# full per-unit T counts (pre-cap) and all-cell counts for denominators
sel['unit_T_total'] = sel.sampleID.map(obs[base].groupby('sampleID').size())
sel['unit_all_total'] = sel.sampleID.map(obs.groupby('sampleID').size())
rows = sel.row.to_numpy()
print('selected', len(rows), sel.groupby(['cohort', 'dataset']).size().to_dict(), flush=True)

C = g['layers']['counts']; ip = C['indptr'][:]
sym = g['var']['gene_symbol']
genes = (np.array(sym['categories'][:], dtype=object)[sym['codes'][:]] if hasattr(sym, 'keys') else sym[:]).astype(str)
data_l, ind_l, nnz = [], [], [0]
W = 20000; t0 = time.time(); i = 0
while i < len(rows):
    lo = rows[i]; hi = lo + W
    j = np.searchsorted(rows, hi)
    block = rows[i:j]
    d = C['data'][ip[lo]:ip[min(hi, len(ip) - 1)]]; x = C['indices'][ip[lo]:ip[min(hi, len(ip) - 1)]]
    for r in block:
        a, b = ip[r] - ip[lo], ip[r + 1] - ip[lo]
        data_l.append(d[a:b]); ind_l.append(x[a:b]); nnz.append(nnz[-1] + (b - a))
    i = j
    if len(data_l) % 50000 < len(block): print('rows', len(data_l), round(time.time() - t0), 's', flush=True)
X = sp.csr_matrix((np.concatenate(data_l).astype(np.float32), np.concatenate(ind_l), np.array(nnz)),
                  shape=(len(rows), len(genes)))
A = ad.AnnData(X=X, obs=sel.drop(columns=['row']).astype(str).set_index(pd.Index(sel.row.astype(str))),
               var=pd.DataFrame(index=pd.Index(genes, name='gene_symbol')))
A.obs['atlas_row'] = sel.row.values
A.var_names_make_unique()
A.layers['counts'] = A.X.copy()
A.write_h5ad(OUT, compression='gzip')
print('done', A.shape, round(time.time() - t0), 's')
