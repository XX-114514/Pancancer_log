"""Arm gene selection by the pre-specified rule (PREREG_bulk_v2.md §1). Single-cell only; no outcome data."""
import os
import zarr, numpy as np, pandas as pd, anndata as ad, numcodecs, json
numcodecs.blosc.set_nthreads(4)
W = os.environ['TFINE_ROOT']
Z = os.environ['PROJECT_V5_ROOT'] + '/runs/20261013_immunotherapy_atlas_v1/outputs/immune_validation_authority_v1.zarr'
cand = {'Treg': ['FOXP3','IL2RA','IKZF2','TNFRSF18','TNFRSF4','CCR8','IL1R2','CTLA4','BATF','LAYN'],
        'Tph': ['CXCL13','CXCR5','MAF','TOX2','PDCD1','ICOS','BTLA','IL21','CD200','SH2D1A'],
        'Tex': ['HAVCR2','LAG3','ENTPD1','LAYN','PDCD1','TOX','CXCL13','TNFRSF9','KIR2DL4','PHLDA1','CTLA4','ITGAE','CD8A','CD8B']}
genes = sorted(set(sum(cand.values(), [])))
S4 = ['GSE294109', 'GSE256326', 'GSE305078', 'GSE286827']
# T states: curated cells from the 4 discovery sources
A = ad.read_h5ad(W + '/T_cells_phaseA.h5ad', backed='r')
st = pd.read_csv(W + '/cells_phaseB_states.tsv.gz', sep='\t', index_col=0)
st.index = st.index.astype(str)
keep = (A.obs.cohort.astype(str) == 'discovery').values
idx = np.where(keep)[0]
gi = [list(A.var_names).index(g) for g in genes]
Xt = A.X[idx][:, gi]; Xt = Xt.toarray() if hasattr(Xt, 'toarray') else np.asarray(Xt)   # log1p CP10k
T = pd.DataFrame(Xt, columns=genes); T['grp'] = pd.Series(st.loc[A.obs_names[idx], 'state'].values).replace({'CD4_CXCL13_Tph': 'Tph'}).values
# non-T majors: sample from atlas counts
o = pd.read_parquet(W + '/obs_small.parquet')
o['row'] = np.arange(len(o))
nonT = o[o.dataset.isin(S4) & ~o.major.isin(['T', 'T_NK_boundary', 'Unresolved', None])]
rng = np.random.default_rng(1)
rows = np.sort(np.concatenate([rng.choice(g.row.values, min(len(g), 4000), replace=False) for _, g in nonT.groupby('major')]))
g = zarr.open(Z, mode='r'); C = g['layers']['counts']; ip = C['indptr'][:]
sym = g['var']['gene_symbol']; vs = (np.array(sym['categories'][:], dtype=object)[sym['codes'][:]] if hasattr(sym, 'keys') else sym[:]).astype(str)
col = {v: i for i, v in enumerate(vs) if v in genes}
out = np.zeros((len(rows), len(genes)), np.float32)
for k, r in enumerate(rows):
    d = C['data'][ip[r]:ip[r + 1]]; x = C['indices'][ip[r]:ip[r + 1]]; tot = d.sum()
    for j, gname in enumerate(genes):
        hit = d[x == col[gname]]
        if len(hit): out[k, j] = np.log1p(hit[0] / tot * 1e4)
N = pd.DataFrame(out, columns=genes); N['grp'] = o.major.values[rows]
allg = pd.concat([T, N])
mean = allg.groupby('grp').mean(); det = (allg.set_index('grp') > 0).groupby(level=0).mean()
nonTm = [m for m in ['Myeloid','B','Plasma','NK','Mast','Epithelial','Endothelial','Fibroblast','Pericyte_SMC'] if m in mean.index]  # exactly PREREG §1
res = []; panels = {}
for arm, gl in cand.items():
    others = [a for a in cand if a != arm]
    kept = []
    for gname in gl:
        a = det.loc[arm, gname] >= 0.25
        b = all(mean.loc[arm, gname] >= 2 * mean.loc[x, gname] for x in others)
        c = all(mean.loc[arm, gname] >= 3 * mean.loc[m, gname] for m in nonTm)
        res.append(dict(arm=arm, gene=gname, det=round(det.loc[arm, gname], 3), mean_arm=round(mean.loc[arm, gname], 3),
                        max_other_arm=round(max(mean.loc[x, gname] for x in others), 3),
                        max_nonT=round(max(mean.loc[m, gname] for m in nonTm), 3), worst_nonT=mean.loc[nonTm, gname].idxmax(),
                        pass_det=a, pass_arm=b, pass_nonT=c, kept=a and b and c))
        if a and b and c: kept.append(gname)
    panels[arm] = kept
R = pd.DataFrame(res); R.to_csv(W + '/bulk_v2/gene_specificity.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); print(R.to_string(index=False)); print(panels)
print('nonT cells per major', N.grp.value_counts().to_dict())
json.dump(panels, open(W + '/bulk_v2/frozen_arm_panels.json', 'x'), indent=1)
