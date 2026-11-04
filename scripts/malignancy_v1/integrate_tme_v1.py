# TME integration v1: integration v2 cells minus consensus-malignant cells (malignancy_v1; uncertain
# cells kept and flagged, per user decision 2026-11-02), re-integrated with the v2 recipe plus the
# changes from the v2 UMAP review:
#  - HVG recomputed on TME cells only (seurat_v3 per source on raw counts, 6500 -> 4000 non-technical)
#  - Harmony batch variables from env TME_HARMONY (v1 run: source theta 2 + unit theta 1; run v2 after the
#    2026-11-05 benchmark: source theta 2 only)
#  - Leiden clustering on the Harmony graph and a coarse marker-based compartment per cluster
#  - integration metrics before/after Harmony (iLISI on source and unit, cLISI on compartment, same-
#    source neighbour fraction) on a fixed subsample
#  - UMAP min_dist 0.3; figures with per-category facets, shuffled draw order, distinct palette
# Inputs: integration v2 shards (X log1p-normalised decontX, layers counts / decontx) and the per-source
# malignancy TSVs. Outputs: OUT/tme_v1.zarr, OUT/figures/*.png, OUT/report/*.
# Usage: integrate_tme_v1.py SHARD_DIR MALIG_DIR OUT_DIR REPO SEED
import sys, os, glob, json, re, gc, datetime, subprocess
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc, harmonypy as hm, zarr
from anndata.io import read_elem, write_elem
from anndata.experimental import concat_on_disk
from sklearn.neighbors import NearestNeighbors

sdir, mdir, out, repo, seed = sys.argv[1:6]
seed = int(seed)
N_CAND, N_HVG, N_PC, N_SUB = 6500, 4000, 50, 200_000
for d in ('shards', 'figures', 'report'):
    os.makedirs(f'{out}/{d}', exist_ok=True)
log = lambda **k: print(json.dumps({'t': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'), **k}, default=str), flush=True)
technical = lambda g: g.upper() == 'MALAT1' or bool(re.search(r'^(MT-|RP[SL][0-9]|HB[ABDEGQMZ][0-9]|ERCC-)', g, re.I))
MARKERS = {'T': ['CD3D', 'CD3E', 'CD2'], 'NK': ['NKG7', 'GNLY', 'KLRD1'], 'B': ['MS4A1', 'CD79A', 'BANK1'],
           'plasma': ['JCHAIN', 'MZB1', 'IGHG1'], 'myeloid': ['CD68', 'LYZ', 'CD14', 'C1QA'], 'mast': ['TPSAB1', 'CPA3', 'KIT'],
           'pDC': ['LILRA4', 'IL3RA'], 'fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'], 'endothelial': ['PECAM1', 'VWF', 'CDH5'],
           'pericyte_SMC': ['RGS5', 'ACTA2', 'MYH11'], 'epithelial': ['EPCAM', 'KRT8', 'KRT18', 'KRT19'],
           'melanocytic': ['PMEL', 'MLANA', 'TYR'], 'glial_neural': ['GFAP', 'OLIG2', 'SOX2'], 'hepatocyte': ['ALB', 'APOA1'],
           'cycling': ['MKI67', 'TOP2A']}
# sources left out of the TME: GSE179994 public data are only the authors' computational T-cell subset
# (GEO protocol has no sorting; non-T cells not deposited), so its composition cannot enter sample-level
# co-variation; kept in integration v2, not deleted (inventories/tme_v1/source_decisions_v1.tsv)
_ex = pd.read_csv(f'{repo}/inventories/tme_v1/scope_exclusions_v1.tsv', sep='\t', dtype=str).fillna('')  # source- and unit-level scope decisions
EXCLUDE = {r.source: r.reason.split(' ')[0] for r in _ex.itertuples() if r.level == 'source'}
EXCLUDE_UNITS = {r.unit_id: r.reason.split(' ')[0] for r in _ex.itertuples() if r.level == 'unit'}
files = sorted(f for f in glob.glob(f'{sdir}/*.h5ad') if os.path.basename(f)[:-5] not in EXCLUDE)
sources = [os.path.basename(f)[:-5] for f in files]
missing = [s for s in sources if not os.path.exists(f'{mdir}/{s}.tsv.gz')]
assert not missing, f'malignancy output missing for {missing}'
MCOLS = ['cell_id', 'quick_cluster', 'lineage', 'scATOMIC_pred', 'scATOMIC_tumour', 'cnv_score', 'cnv_corr', 'cnv_call', 'cnv_threshold',
         'malignancy_call', 'malignancy_rule']


def load(f, s):
    a = ad.read_h5ad(f)
    m = pd.read_csv(f'{mdir}/{s}.tsv.gz', sep='\t', usecols=MCOLS, dtype={'quick_cluster': str}).set_index('cell_id')
    assert m.index.is_unique and set(m.index) == set(a.obs_names), f'{s}: malignancy cells != shard cells'
    m = m.reindex(a.obs_names)
    for c in m.columns:
        a.obs['malig_' + c if not c.startswith('malignancy') else c] = m[c].values
    a.obs['tme_flag'] = np.where(m['malignancy_call'].values == 'uncertain', 'uncertain_malignancy_kept', 'non_malignant')
    return a, (m['malignancy_call'].values != 'malignant') & ~a.obs['unit_id'].astype(str).isin(list(EXCLUDE_UNITS)).values


common = sorted(set.intersection(*[set(ad.read_h5ad(f, backed='r').var_names.astype(str)) for f in files]))
log(step='common_genes', n=len(common), sources=len(sources))

# pass 1: per-source seurat_v3 on raw counts of TME cells; malignancy composition per source
ranks, comp = {}, []
for f, s in zip(files, sources):
    a, keep = load(f, s)
    comp.append(a.obs.groupby('malignancy_call', observed=True).size().rename(s))
    b = ad.AnnData(X=a.layers['counts'][keep][:, a.var_names.get_indexer(common)].astype(np.float32), var=pd.DataFrame(index=common))
    sc.pp.highly_variable_genes(b, n_top_genes=N_CAND, flavor='seurat_v3')
    ranks[s] = b.var['highly_variable_rank'].values.astype(float)
    log(step='hvg_source', source=s, cells=int(a.n_obs), tme=int(keep.sum()))
    del a, b; gc.collect()
pd.concat(comp, axis=1).T.fillna(0).astype(int).to_csv(f'{out}/report/malignancy_by_source.tsv', sep='\t')
R = pd.DataFrame(ranks, index=common)
hv = pd.DataFrame({'highly_variable_rank': np.nanmedian(R.values, axis=1), 'highly_variable_nbatches': R.notna().sum(1).values}, index=common)
hv['technical'] = [technical(g) for g in common]
cand = hv.sort_values(['highly_variable_nbatches', 'highly_variable_rank'], ascending=[False, True]).index[:N_CAND]
elig = hv.loc[cand][~hv.loc[cand, 'technical']].sort_values(['highly_variable_rank', 'highly_variable_nbatches'], ascending=[True, False], kind='mergesort')
selected = list(elig.index[:N_HVG])
hv['candidate'] = hv.index.isin(cand); hv['selected_HVG'] = hv.index.isin(selected)
hv.to_csv(f'{out}/report/hvg_tme_v1.tsv', sep='\t')
log(step='hvg_selected', n=len(selected))

# pass 2: TME shards (X already log1p-normalised decontX in v2 shards), HVG matrix, marker matrix
mk_genes = sorted({g for gs in MARKERS.values() for g in gs if g in common})
blocks, mblocks, obs_all = [], [], []
OBS = ['source', 'data_batch', 'analysis_main', 'oncotree_code', 'sampling_site', 'treatment_status', 'patient_id', 'unit_id',
       'technology', 'fresh_frozen', 'tme_flag', 'malig_lineage']
for f, s in zip(files, sources):
    a, keep = load(f, s)
    a = a[keep].copy()
    a.write_h5ad(f'{out}/shards/{s}.h5ad')
    blocks.append(a[:, selected].X.copy()); mblocks.append(a[:, mk_genes].X.copy())
    obs_all.append(a.obs[OBS].copy())
    log(step='tme_shard', source=s, cells=int(a.n_obs))
    del a; gc.collect()
Y = sp.vstack(blocks).tocsr(); M = sp.vstack(mblocks).tocsr(); del blocks, mblocks; gc.collect()
obs = pd.concat(obs_all); del obs_all
for c in obs.columns:
    obs[c] = obs[c].astype(str).astype('category')
log(step='matrix', shape=Y.shape)

B = ad.AnnData(X=Y, obs=obs, var=pd.DataFrame(index=selected))
sc.pp.pca(B, n_comps=N_PC, zero_center=True, svd_solver='arpack', random_state=seed, dtype='float32')
log(step='pca_done')
# Harmony batch variables, env TME_HARMONY="var:theta,..." (default = v1 run: source 2 + unit_id 1). The
# batch-correction benchmark v1 (2026-11-05) chose source:2 for the next run (inventories/tme_v1/batch_benchmark_v1/)
HV = [(k, float(t)) for k, t in (x.split(':') for x in os.environ.get('TME_HARMONY', 'source:2,unit_id:1').split(','))]
log(step='harmony_vars', vars=HV)
ho = hm.run_harmony(data_mat=np.asarray(B.obsm['X_pca'], np.float64), meta_data=B.obs[[k for k, _ in HV]].astype(str),
                    vars_use=[k for k, _ in HV], theta=[t for _, t in HV], lamb=1.0, max_iter_harmony=20, random_state=seed)
Z = np.asarray(ho.Z_corr).T if ho.Z_corr.shape[0] == N_PC else np.asarray(ho.Z_corr)
B.obsm['X_pca_harmony'] = Z.astype(np.float32); del ho; gc.collect()
log(step='harmony_done')
B.X = sp.csr_matrix(B.shape, dtype=np.float32)
sc.pp.neighbors(B, n_neighbors=15, use_rep='X_pca_harmony', random_state=seed)
sc.tl.leiden(B, resolution=1.0, random_state=seed, flavor='igraph', n_iterations=2, directed=False, key_added='leiden_tme')
log(step='leiden_done', clusters=int(B.obs['leiden_tme'].nunique()))
sc.tl.umap(B, min_dist=0.3, random_state=seed)
log(step='umap_done')

# coarse compartment per cluster: highest mean marker-set z-score (z across clusters)
cl = B.obs['leiden_tme'].values
Mc = pd.DataFrame(M.toarray(), columns=mk_genes).groupby(cl).mean()
S = pd.DataFrame({k: Mc[[g for g in gs if g in mk_genes]].mean(1) for k, gs in MARKERS.items() if k != 'cycling'})
Zs = (S - S.mean()) / S.std(ddof=0)
lab = Zs.idxmax(1)
cyc = Mc[[g for g in MARKERS['cycling'] if g in mk_genes]].mean(1)
lab = lab.where(Zs.max(1) >= 1.0, 'unresolved')
comp_name = {'T': 'T_NK', 'NK': 'T_NK', 'B': 'B_plasma', 'plasma': 'B_plasma', 'myeloid': 'myeloid', 'mast': 'myeloid', 'pDC': 'myeloid',
             'fibroblast': 'stromal', 'endothelial': 'endothelial', 'pericyte_SMC': 'stromal', 'epithelial': 'epithelial_nonmalignant',
             'melanocytic': 'melanocytic_nonmalignant', 'glial_neural': 'glial_neural', 'hepatocyte': 'epithelial_nonmalignant', 'unresolved': 'unresolved'}
ctab = pd.DataFrame({'cells': pd.Series(cl).value_counts(), 'marker_label': lab, 'compartment': lab.map(comp_name),
                     'cycling_score': cyc, 'max_z': Zs.max(1),
                     'uncertain_frac': B.obs.groupby('leiden_tme', observed=True)['tme_flag'].apply(lambda x: (x != 'non_malignant').mean())})
ctab = ctab.join(S.add_prefix('mean_')).sort_values('cells', ascending=False)
ctab.to_csv(f'{out}/report/leiden_clusters.tsv', sep='\t')
B.obs['marker_label'] = pd.Categorical(ctab['marker_label'].reindex(cl).values)
B.obs['compartment'] = pd.Categorical(ctab['compartment'].reindex(cl).values)
del M; gc.collect()

# integration metrics on a fixed subsample
rng = np.random.default_rng(seed)
sub = np.sort(rng.choice(B.n_obs, min(N_SUB, B.n_obs), replace=False))


def lisi(emb, labels, k=90):
    nn = NearestNeighbors(n_neighbors=k, n_jobs=-1).fit(emb)
    idx = nn.kneighbors(emb, return_distance=False)
    codes = pd.Categorical(labels).codes
    L = codes[idx]
    simpson = np.array([np.sum((np.bincount(r) / k) ** 2) for r in L])
    same = (L[:, 1:] == codes[:, None]).mean(1)
    return 1 / simpson, same


metrics = {}
o = B.obs.iloc[sub]
for emb in ('X_pca', 'X_pca_harmony'):
    E = B.obsm[emb][sub]
    for key in ('source', 'unit_id', 'compartment'):
        li, same = lisi(E, o[key].astype(str).values)
        d = pd.DataFrame({'lisi': li, 'same': same, 'comp': o['compartment'].astype(str).values})
        metrics[f'{emb}|{key}'] = {'lisi_median': float(np.median(li)), 'same_frac_median': float(np.median(same)),
                                   'by_compartment_lisi_median': d.groupby('comp')['lisi'].median().round(3).to_dict(),
                                   'by_compartment_same_median': d.groupby('comp')['same'].median().round(3).to_dict()}
    log(step='metrics', emb=emb)
metrics['note'] = 'iLISI (source, unit_id): higher = better mixing (max = number of batches); cLISI (compartment): close to 1 = compartments kept separate; k=90 unweighted neighbours, 200k-cell subsample'
json.dump(metrics, open(f'{out}/report/integration_metrics.json', 'w'), indent=1)

# figures
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
U = B.obsm['X_umap']; order = rng.permutation(B.n_obs)
PAL = list(plt.get_cmap('tab20').colors) + list(plt.get_cmap('tab20b').colors) + list(plt.get_cmap('Set3').colors)


def overview(key, fn, top=None):
    v = B.obs[key].astype(str).values
    cats = pd.Series(v).value_counts().index.tolist()
    if top and len(cats) > top:
        v = np.where(np.isin(v, cats[:top]), v, 'other'); cats = cats[:top] + ['other']
    cmap = {c: (PAL[i % len(PAL)] if c != 'other' else (0.8, 0.8, 0.8)) for i, c in enumerate(cats)}
    fig, ax = plt.subplots(figsize=(11, 8))
    ax.scatter(U[order, 0], U[order, 1], c=[cmap[x] for x in v[order]], s=0.05, linewidths=0, rasterized=True)
    ax.set_axis_off(); ax.set_title(f'TME v1 UMAP ({B.n_obs:,} cells): {key}')
    for c in cats:
        ax.scatter([], [], c=[cmap[c]], s=20, label=f'{c} ({(v == c).sum():,})')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=6, frameon=False, ncol=1 if len(cats) <= 30 else 2)
    fig.savefig(f'{out}/figures/{fn}.png', dpi=200, bbox_inches='tight'); plt.close(fig)


def facets(key, fn, maxn=40):
    v = B.obs[key].astype(str).values
    cats = pd.Series(v).value_counts().index.tolist()[:maxn]
    n = len(cats); nc = 6; nr = int(np.ceil(n / nc))
    fig, axs = plt.subplots(nr, nc, figsize=(nc * 3, nr * 3))
    bg = rng.choice(B.n_obs, min(300_000, B.n_obs), replace=False)
    for ax, c in zip(np.ravel(axs), cats):
        ax.scatter(U[bg, 0], U[bg, 1], c='#dddddd', s=0.02, linewidths=0, rasterized=True)
        w = np.flatnonzero(v == c)
        ax.scatter(U[w, 0], U[w, 1], c='#c0392b', s=0.05 if len(w) > 50_000 else 0.3, linewidths=0, rasterized=True)
        ax.set_title(f'{c[:30]} ({len(w):,})', fontsize=7); ax.set_axis_off()
    for ax in np.ravel(axs)[n:]:
        ax.set_axis_off()
    fig.tight_layout(); fig.savefig(f'{out}/figures/{fn}.png', dpi=150); plt.close(fig)


overview('compartment', 'umap_compartment'); overview('marker_label', 'umap_marker_label')
overview('analysis_main', 'umap_cancer'); overview('source', 'umap_source', top=30)
overview('tme_flag', 'umap_tme_flag'); overview('sampling_site', 'umap_sampling_site'); overview('treatment_status', 'umap_treatment')
facets('analysis_main', 'facet_cancer'); facets('source', 'facet_source', maxn=72); facets('compartment', 'facet_compartment')
fig, ax = plt.subplots(figsize=(9, 8))
ax.hexbin(U[:, 0], U[:, 1], gridsize=300, bins='log', cmap='viridis', mincnt=1); ax.set_axis_off(); ax.set_title('cell density (log)')
fig.savefig(f'{out}/figures/umap_density.png', dpi=200, bbox_inches='tight'); plt.close(fig)
log(step='figures_done')

# final zarr
zpath = f'{out}/tme_v1.zarr'
concat_on_disk([f'{out}/shards/{s}.h5ad' for s in sources], zpath, join='outer', max_loaded_elems=200_000_000)
g = zarr.open(zpath, 'r+')
zo = read_elem(g['obs'])
assert list(zo.index) == list(B.obs_names), 'cell order mismatch'
for k in ('leiden_tme', 'marker_label', 'compartment'):
    zo[k] = B.obs[k].values
del g['obs']; write_elem(g, 'obs', zo)
for k in ('X_pca', 'X_pca_harmony', 'X_umap'):
    write_elem(g['obsm'] if 'obsm' in g else g.create_group('obsm'), k, B.obsm[k])
var = read_elem(g['var']); var['selected_HVG'] = var.index.isin(selected)
del g['var']; write_elem(g, 'var', var)
commit = subprocess.run(['git', '-C', repo, 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
write_elem(g, 'uns', {'schema': 'pancancer_tme_v1', 'code_commit': commit, 'seed': str(seed),
                      'cells': 'integration v2 cells with malignancy_call != malignant (uncertain kept, tme_flag)',
                      'excluded_sources': json.dumps(EXCLUDE), 'excluded_units': json.dumps(EXCLUDE_UNITS),
                      'X': 'log1p(normalize_total 1e4) of decontX-corrected counts', 'layers': 'counts = raw integer; decontx = decontX-corrected',
                      'hvg': f'seurat_v3 per source on TME raw counts of {len(common)} common genes; {N_CAND} -> {N_HVG} non-technical',
                      'pca': f'{N_PC} PCs, centred, ARPACK', 'harmony': 'vars ' + ', '.join(f'{k} (theta {t:g})' for k, t in HV) + ', lambda 1, max_iter 20',
                      'graph': '15-NN on X_pca_harmony; Leiden 1.0 (leiden_tme); UMAP min_dist 0.3',
                      'compartment': 'per Leiden cluster, highest marker-set z-score >= 1 else unresolved'})
zarr.consolidate_metadata(zpath)
json.dump({'cells': int(B.n_obs), 'sources': len(sources), 'units': int(B.obs['unit_id'].nunique()), 'patients': int(B.obs['patient_id'].nunique()),
           'tme_flag': B.obs['tme_flag'].value_counts().to_dict(), 'compartment': B.obs['compartment'].value_counts().to_dict(),
           'clusters': int(B.obs['leiden_tme'].nunique()), 'common_genes': len(common), 'hvg': len(selected)},
          open(f'{out}/report/tme_summary.json', 'w'), indent=1)
log(step='done', cells=int(B.n_obs))
