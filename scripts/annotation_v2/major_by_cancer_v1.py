# Annotation v2 step A2 (2026-11-06): major cell types within ONE cancer type (analysis_main) of TME v1 run v2,
# following the design of Lodi et al. 2025 with the project's parameter choices (rationale table in
# methods/tme-annotation-v2-lodi-workflow-20261105.md):
#  - cells: all TME cells of the cancer type (non_malignant + uncertain; uncertain stays flagged)
#  - read : per dataset shard, contiguous row blocks, only needed genes (common genes + registry genes)
#  - HVG  : seurat_v3 per source on raw counts (merged ranks), 2,000, excluding technical (MT-/RPS/RPL/MALAT1),
#           haemoglobin, IG V and TCR V genes
#  - state regression before PCA (Lodi's reason: state programs otherwise form pseudo-subtypes): S and G2M
#    (Seurat cc.genes.updated.2019), IFN (MSigDB BROWNE), stress (V7 Stress_response), hypoxia (Lodi 15 genes),
#    scored like Seurat AddModuleScore (scanpy score_genes); sample ID and cancer type are NOT regressed
#  - scale (clip 10), PCA 50, number of PCs = Kneedle elbow of the variance ratio clipped to [15, 40]
#  - Harmony on source (theta 2) when the cancer type has >1 source; 15-NN; Leiden 0.3/0.6/1.0/1.5 + adjacent ARI
#  - evidence for the main resolution: Seurat-style DE (Wilcoxon on <=2,000 cells per cluster; avg_log2FC =
#    log2(mean(expm1)+1) difference; Bonferroni p_val_adj), expression detection over measured cells, and the V7
#    evidence engine on the merged registry v2 (no winner); plus QC, doublet evidence (A0), malignancy, source mix
# Outputs (OUT/<cancer_slug>/): object.h5ad, clusters.tsv, review.tsv, de.tsv, expression.tsv, evidence/, figures/
# Usage: major_by_cancer_v1.py RUN_DIR CANCER OUT_ROOT REGISTRY_DIR DOUBLET_DIR REPO SEED
import sys, os, re, json, datetime, subprocess
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc, harmonypy as hm, zarr, h5py
from anndata.io import read_elem
from scipy.stats import rankdata
from sklearn.metrics import adjusted_rand_score

run, cancer, outroot, regdir, dbldir, repo, seed = sys.argv[1:8]
seed = int(seed)
slug = re.sub(r'[^A-Za-z0-9]+', '_', cancer).strip('_')
out = f'{outroot}/{slug}'
os.makedirs(f'{out}/figures', exist_ok=True)
log = lambda **k: print(json.dumps({'t': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'), **k}, default=str), flush=True)
BLOCK, MIN_SRC, N_HVG, RES, MAIN = 50_000, 200, 2000, (0.3, 0.6, 1.0, 1.5), 1.0
excluded = lambda g: g.upper() == 'MALAT1' or bool(re.search(r'^(MT-|RP[SL][0-9]|HB[ABDEGQMZ][0-9]?$|HB[ABDEGQMZ][0-9]|IG[HKL]V|TR[ABGD]V)', g, re.I))
CONTEXT = {  # registry required_context strings that apply to each cancer type's primary organ
    'NSCLC': ['lung', 'lung (NSCLC)', 'squamous_or_airway', 'airway_or_ciliated_organ'], 'Breast cancer': ['breast'],
    'Head and neck squamous cell carcinoma': ['oesophagus / head and neck / cervix', 'squamous_or_airway'],
    'Renal cell carcinoma': ['kidney'], 'Renal medullary carcinoma': ['kidney'], 'Cutaneous melanoma': ['skin'], 'Acral melanoma': ['skin'],
    'Colorectal cancer': ['intestine', 'colon/rectum', 'GI tract or pancreas'], 'Pancreatic ductal adenocarcinoma': ['pancreas', 'GI tract or pancreas'],
    'Pancreatic neuroendocrine tumor': ['pancreas', 'GI tract or pancreas'], 'Diffuse glioma': ['brain (glioma)'],
    'Bladder cancer': ['bladder / upper urinary tract'], 'Prostate cancer': ['prostate'], 'Gastric cancer': ['stomach', 'intestine', 'GI tract or pancreas'],
    'Osteosarcoma': ['bone (osteosarcoma)'], 'Ovarian cancer': ['ovary / fallopian tube / endometrium'],
    'Endometrial cancer': ['ovary / fallopian tube / endometrium', 'endometrium'], 'Papillary thyroid carcinoma': ['thyroid'],
    'Hepatocellular carcinoma': ['liver (HCC, iCCA)', 'liver', 'hepatobiliary'], 'Intrahepatic cholangiocarcinoma': ['liver (HCC, iCCA)', 'liver', 'hepatobiliary'],
    'Esophageal cancer': ['oesophagus / head and neck / cervix', 'squamous_or_airway', 'stomach'],
    'Cervical cancer': ['oesophagus / head and neck / cervix', 'squamous_or_airway'], 'Nasopharyngeal carcinoma': ['oesophagus / head and neck / cervix', 'squamous_or_airway', 'airway_or_ciliated_organ'],
    'Basal cell carcinoma': ['skin', 'squamous_or_airway'], 'Mesothelioma': ['lung'], 'Uveal melanoma': [], 'Soft tissue sarcoma': []}
ctx = CONTEXT.get(cancer, []) + ['tumour']

# ---- registry and gene sets ----
R = json.load(open(f'{regdir}/registry.json'))
pgenes = sorted({g for p in R['panels'] for g in p['positive_genes']})
cc = pd.read_csv(f'{regdir}/cell_cycle_genes_seurat.tsv', sep='\t')
ifn = list(pd.read_csv(f'{regdir}/BROWNE_INTERFERON_RESPONSIVE_GENES.tsv', sep='\t')['gene'])
byid = {p['id']: p for p in R['panels']}
STATES = {'S': list(cc.gene[cc.phase == 'S']), 'G2M': list(cc.gene[cc.phase == 'G2M']), 'IFN': ifn,
          'stress': byid['Stress_response']['positive_genes'], 'hypoxia': byid['Lodi_hypoxia_15']['positive_genes']}

# ---- select and read ----
g = zarr.open(f'{run}/tme_v1.zarr', 'r')
obs_all = read_elem(g['obs'])
sel_all = (obs_all['analysis_main'].astype(str) == cancer).values
common = list(pd.read_csv(f'{run}/report/hvg_tme_v1.tsv', sep='\t', index_col=0).index)
genes = list(dict.fromkeys(common + pgenes + [x for v in STATES.values() for x in v]))
log(step='select', cancer=cancer, cells=int(sel_all.sum()), genes=len(genes), contexts=ctx)


def read_rows_cols(grp, rows, cols, ncol_out):
    indptr = grp['indptr'][:]; parts = []
    for a in range(0, len(indptr) - 1, BLOCK):
        b = min(a + BLOCK, len(indptr) - 1); r = rows[(rows >= a) & (rows < b)]
        if not len(r):
            continue
        lo, hi = int(indptr[a]), int(indptr[b])
        blk = sp.csr_matrix((grp['data'][lo:hi], grp['indices'][lo:hi], indptr[a:b + 1] - lo), shape=(b - a, len(cols)))[r - a].tocoo()
        k = cols[blk.col] >= 0
        parts.append(sp.csr_matrix((blk.data[k].astype(np.float32), (blk.row[k], cols[blk.col[k]])), shape=(len(r), ncol_out)))
    return sp.vstack(parts).tocsr()


src = obs_all['source'].astype(str).values
X, C, O, measured = [], [], [], {}
for s in dict.fromkeys(src):
    pos = np.flatnonzero(src == s); w = np.flatnonzero(sel_all[pos])
    if not len(w):
        continue
    with h5py.File(f'{run}/shards/{s}.h5ad', 'r') as f:
        vn = f['var'][f['var'].attrs['_index']].asstr()[:]
        colmap = pd.Index(genes).get_indexer(vn); measured[s] = np.isin(genes, vn)
        X.append(read_rows_cols(f['X'], w, colmap, len(genes))); C.append(read_rows_cols(f['layers/counts'], w, colmap, len(genes)))
    O.append(obs_all.iloc[pos[w]])
A = ad.AnnData(X=sp.vstack(X).tocsr(), obs=pd.concat(O), var=pd.DataFrame(index=genes))
A.layers['counts'] = sp.vstack(C).tocsr(); del X, C
M = pd.DataFrame(measured, index=genes).T
# doublet evidence (A0), if available
dbl = []
for s in A.obs['source'].astype(str).unique():
    fp = f'{dbldir}/{s}.tsv.gz'
    if os.path.exists(fp):
        d = pd.read_csv(fp, sep='\t', usecols=['cell_id', 'doublet_score', 'scrublet_predicted'])
        dbl.append(d[d.cell_id.isin(A.obs_names)].set_index('cell_id'))
if dbl:
    D = pd.concat(dbl).reindex(A.obs_names)
    A.obs['doublet_score'] = D['doublet_score'].values
    A.obs['scrublet_predicted'] = pd.to_numeric(D['scrublet_predicted'], errors='coerce').values
else:
    A.obs['doublet_score'] = np.nan; A.obs['scrublet_predicted'] = np.nan
log(step='loaded', shape=A.shape, doublet_scored=int(A.obs['doublet_score'].notna().sum()))

# ---- HVG (per source, merged ranks), state scores, regression, scale, PCA ----
A.obs['hvg_batch'] = A.obs['source'].astype(str)
vc = A.obs['hvg_batch'].value_counts()
A.obs.loc[A.obs['hvg_batch'].isin(vc.index[vc < MIN_SRC]), 'hvg_batch'] = 'pooled_small_sources'
cand = np.isin(genes, common) & ~np.array([excluded(x) for x in genes])
B = ad.AnnData(X=A.layers['counts'][:, cand], obs=A.obs[['hvg_batch']].copy(), var=pd.DataFrame(index=np.array(genes)[cand]))
nb = A.obs['hvg_batch'].nunique()
sc.pp.highly_variable_genes(B, n_top_genes=N_HVG, flavor='seurat_v3', batch_key='hvg_batch' if nb > 1 else None)
hvg = list(B.var_names[B.var['highly_variable']]); del B
for k, gs in STATES.items():
    gs = [x for x in gs if x in A.var_names]
    sc.tl.score_genes(A, gs, score_name=f'score_{k}', random_state=seed, use_raw=False)
Y = A[:, hvg].X.toarray().astype(np.float32)
Dm = np.column_stack([np.ones(A.n_obs)] + [A.obs[f'score_{k}'].values for k in STATES]).astype(np.float32)
beta, *_ = np.linalg.lstsq(Dm, Y, rcond=None)
Y -= Dm @ beta
Y = (Y - Y.mean(0)) / np.maximum(Y.std(0), 1e-6); np.clip(Y, -10, 10, out=Y)
Bp = ad.AnnData(X=Y, obs=A.obs[['source']].copy())
sc.pp.pca(Bp, n_comps=50, svd_solver='arpack', random_state=seed)
vr = Bp.uns['pca']['variance_ratio']
xs = np.arange(len(vr)); ys = (vr - vr.min()) / (vr.max() - vr.min())
line = 1 - xs / xs[-1]
n_pc = int(np.clip(np.argmax(line - ys) + 1, 15, 40))
A.obsm['X_pca'] = Bp.obsm['X_pca'][:, :n_pc]; del Y, Bp
if A.obs['source'].nunique() > 1:
    h = hm.run_harmony(np.asarray(A.obsm['X_pca'], np.float64), A.obs[['source']].astype(str), vars_use=['source'], theta=2.0, lamb=1.0, max_iter_harmony=20, random_state=seed)
    Z = np.asarray(h.Z_corr); A.obsm['X_emb'] = (Z.T if Z.shape[0] == n_pc else Z).astype(np.float32)
else:
    A.obsm['X_emb'] = A.obsm['X_pca']
sc.pp.neighbors(A, n_neighbors=15, use_rep='X_emb', random_state=seed)
for r in RES:
    sc.tl.leiden(A, resolution=r, key_added=f'leiden_r{r:g}', random_state=seed, flavor='igraph', n_iterations=2, directed=False)
sc.tl.umap(A, min_dist=0.3, random_state=seed)
ari = {f'r{a:g}-r{b:g}': round(adjusted_rand_score(A.obs[f'leiden_r{a:g}'], A.obs[f'leiden_r{b:g}']), 3) for a, b in zip(RES, RES[1:])}
K = f'leiden_r{MAIN:g}'
log(step='embedded', hvg=len(hvg), n_pc=n_pc, clusters={r: int(A.obs[f'leiden_r{r:g}'].nunique()) for r in RES}, adjacent_ARI=ari)

# ---- expression detection (measured cells only) and Seurat-style DE ----
cl = A.obs[K].astype(str).values; clusters = sorted(set(cl), key=int)
srcs = A.obs['source'].astype(str).values
pg = [x for x in pgenes if x in A.var_names]; gi = pd.Index(genes).get_indexer(pg)
E = A.X[:, gi].tocsc()
rows = []
for c in clusters:
    m = cl == c; n = int(m.sum())
    meas = M.loc[srcs[m], pg].values  # cells x genes measured mask
    nmeas = meas.sum(0)
    det = np.asarray((E[m] > 0).multiply(meas).sum(0)).ravel()
    lin = np.asarray(E[m].expm1().multiply(meas).sum(0)).ravel()
    rows.append(pd.DataFrame({'cluster': c, 'gene': pg, 'n_cells': n, 'n_measured': nmeas.astype(int),
                              'pct_expressing_measured': np.where(nmeas > 0, det / np.maximum(nmeas, 1), np.nan),
                              'mean_cp10k_measured': np.where(nmeas > 0, lin / np.maximum(nmeas, 1), np.nan)}))
EXP = pd.concat(rows); EXP.to_csv(f'{out}/expression.tsv', sep='\t', index=False)
rng = np.random.default_rng(seed)
take = np.sort(np.concatenate([rng.choice(np.flatnonzero(cl == c), min(2000, (cl == c).sum()), replace=False) for c in clusters]))
Dsub = A[take].copy(); cls = Dsub.obs[K].astype(str).values
common_measured = [x for x in genes if x in common]  # DE only on genes measured in every source
Dd = Dsub[:, common_measured].copy()
sc.tl.rank_genes_groups(Dd, groupby=K, method='wilcoxon', use_raw=False, tie_correct=True)  # one-vs-rest Wilcoxon p-values
Xs = Dd.X.tocsc()
de_rows = []
for c in clusters:
    m = cls == c
    xin, xout = Xs[m], Xs[~m]
    p1 = np.asarray((xin > 0).mean(0)).ravel(); p2 = np.asarray((xout > 0).mean(0)).ravel()
    fc = np.log2(np.asarray(xin.expm1().mean(0)).ravel() + 1) - np.log2(np.asarray(xout.expm1().mean(0)).ravel() + 1)  # Seurat avg_log2FC
    keep = (np.maximum(p1, p2) >= 0.10) & (np.abs(fc) >= 0.25)
    rg = sc.get.rank_genes_groups_df(Dd, group=c).set_index('names')['pvals']
    idx = np.flatnonzero(keep)
    p = rg.reindex(np.array(common_measured)[idx]).values
    de_rows.append(pd.DataFrame({'cluster': c, 'gene': np.array(common_measured)[idx], 'avg_log2FC': fc[idx], 'pct.1': p1[idx], 'pct.2': p2[idx],
                                 'p_val': p, 'p_val_adj': np.minimum(1, p * len(common_measured))}))  # Bonferroni, as Seurat
DE = pd.concat(de_rows); DE.to_csv(f'{out}/de.tsv', sep='\t', index=False)
del Dsub, Dd, Xs
contexts = {c: {'organ_contexts': ctx} for c in clusters}
json.dump(contexts, open(f'{out}/contexts.json', 'w'))
ev = f'{out}/evidence'
if os.path.exists(ev):
    import shutil; shutil.rmtree(ev)
subprocess.run([sys.executable, '-I', f'{repo}/scripts/v7_panel_repair_v1/evidence_engine.py', '--registry', f'{regdir}/registry.json',
                '--expression', f'{out}/expression.tsv', '--de', f'{out}/de.tsv', '--contexts', f'{out}/contexts.json', '--out', ev], check=True,
               stdout=subprocess.DEVNULL)
EV = pd.read_csv(f'{ev}/all_panel_evidence.tsv', sep='\t', dtype={'cluster': str})

# ---- cluster summary + compact review sheet ----
o = A.obs; grp = o.groupby(K, observed=True)
S = pd.DataFrame({'cells': grp.size(), 'non_malignant': grp['tme_flag'].agg(lambda x: int((x.astype(str) == 'non_malignant').sum())),
                  'uncertain_share': grp['tme_flag'].agg(lambda x: (x.astype(str) != 'non_malignant').mean()).round(3),
                  'global_c13_share': grp['leiden_tme'].agg(lambda x: (x.astype(str) == '13').mean()).round(3),
                  'sources': grp['source'].nunique(), 'units': grp['unit_id'].nunique(),
                  'top_source_share': grp['source'].agg(lambda x: x.value_counts(normalize=True).iloc[0]).round(3),
                  'median_genes': grp['n_genes'].median(), 'median_umi': grp['total_counts'].median(), 'median_pct_mt': grp['pct_mt'].median().round(2),
                  'median_contamination': grp['decontX_contamination'].median().round(3),
                  'doublet_score_mean': grp['doublet_score'].mean().round(3), 'scrublet_pred_frac': grp['scrublet_predicted'].mean().round(3),
                  'cnv_corr_median': grp['malig_cnv_corr'].median().round(3),
                  **{f'{k}_score_mean': grp[f'score_{k}'].mean().round(3) for k in STATES}})
S.index = S.index.astype(str)
S['genes_pct_of_object_p10'] = (S['median_genes'] < np.percentile(o['n_genes'], 10))
S['doublet_score_z'] = ((S['doublet_score_mean'] - S['doublet_score_mean'].mean()) / S['doublet_score_mean'].std()).round(2)
ok = (DE['avg_log2FC'] >= .25) & (DE['pct.1'] >= .10) & ((DE['pct.1'] - DE['pct.2']) >= .25) & (DE['p_val_adj'] < .01)
S['top_markers'] = DE[ok].sort_values('avg_log2FC', ascending=False).groupby('cluster')['gene'].agg(lambda x: ','.join(x.head(15))).reindex(S.index)
# Lodi major-type detection (max pct among the set) and V7 parent evidence (min over core components of max detection)
# Lodi sets mix alternatives (CD4 vs CD8A), so show the best-detected gene of each set, not the minimum over genes
ep = EXP.pivot(index='cluster', columns='gene', values='pct_expressing_measured')
lm = pd.DataFrame({p['id'].replace('Lodi_major_', 'Lodi:'): ep[[x for x in p['positive_genes'] if x in ep.columns]].max(1)
                   for p in R['panels'] if p['id'].startswith('Lodi_major_')})
par = EV[(EV.axis == 'parent')].pivot(index='cluster', columns='panel', values='min_component_max_detection')
par_deg = EV[(EV.axis == 'parent')].set_index(['cluster', 'panel'])['positive_DEG_genes']
S = S.join(lm.round(2))
S['V7_parents_ge0.3'] = [';'.join(f'{k}:{v:.2f}' for k, v in par.loc[c].dropna().sort_values(ascending=False).items() if v >= 0.3) for c in S.index]
S['V7_parents_with_DEG'] = [';'.join(f'{p}[{par_deg[(c, p)]}]' for p in par.columns if isinstance(par_deg.get((c, p)), str)) for c in S.index]
tis = EV[EV.axis.isin(['branch', 'tissue_specific_population', 'parent']) & EV.positive_DEG_genes.notna() & (EV.organ_context != 'organ_context_missing')]
S['context_modules_with_DEG'] = tis.groupby('cluster').apply(lambda d: ';'.join(f'{a}[{b}]' for a, b in zip(d.panel, d.positive_DEG_genes))).reindex(S.index)
S = S.sort_values('cells', ascending=False)
S.to_csv(f'{out}/clusters.tsv', sep='\t')
A.obs['cancer_major_cluster'] = A.obs[K].astype(str)

# ---- figures ----
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
U = A.obsm['X_umap']; order = rng.permutation(A.n_obs)
pal = list(plt.get_cmap('tab20').colors) + list(plt.get_cmap('tab20b').colors) + list(plt.get_cmap('Set3').colors)
def scat(key, fn, cont=False, top=40):
    fig, ax = plt.subplots(figsize=(9, 7.5))
    if cont:
        v = pd.to_numeric(A.obs[key], errors='coerce').values
        s_ = ax.scatter(U[order, 0], U[order, 1], c=v[order], s=0.5, cmap='viridis', linewidths=0, rasterized=True); fig.colorbar(s_, ax=ax, shrink=.6)
    else:
        v = A.obs[key].astype(str).values; cats = pd.Series(v).value_counts().index[:top].tolist()
        v = np.where(np.isin(v, cats), v, 'other'); cm = {c: pal[i % len(pal)] for i, c in enumerate(cats)}; cm['other'] = (.85, .85, .85)
        ax.scatter(U[order, 0], U[order, 1], c=[cm[x] for x in v[order]], s=0.5, linewidths=0, rasterized=True)
        for c in cats: ax.scatter([], [], c=[cm[c]], s=12, label=f'{c} ({(v == c).sum():,})')
        ax.legend(loc='center left', bbox_to_anchor=(1, .5), fontsize=6, frameon=False)
        if key == K:
            for c in cats: xy = np.median(U[v == c], 0); ax.text(xy[0], xy[1], c, fontsize=7, weight='bold')
    ax.set_axis_off(); ax.set_title(f'{cancer}: {key} ({A.n_obs:,} cells)', fontsize=9)
    fig.savefig(f'{out}/figures/{fn}.png', dpi=150, bbox_inches='tight'); plt.close(fig)
scat(K, 'umap_cluster'); scat('source', 'umap_source'); scat('tme_flag', 'umap_tme_flag'); scat('doublet_score', 'umap_doublet', cont=True)
scat('n_genes', 'umap_n_genes', cont=True); scat('malig_cnv_corr', 'umap_cnv_corr', cont=True); scat('compartment', 'umap_global_compartment')
# Lodi 29-gene + key parent genes dot plot by cluster
dot_genes = ['PTPRC', 'CD3D', 'CD3E', 'CD4', 'CD8A', 'NKG7', 'CD79A', 'MS4A1', 'MZB1', 'JCHAIN', 'CD68', 'LYZ', 'AIF1', 'ITGAM', 'CSF3R', 'LILRA4', 'CXCR3',
             'CLEC9A', 'XCR1', 'CD1C', 'CCR7', 'CCL17', 'CCL19', 'MS4A2', 'TPSAB1', 'CPA3', 'CLDN5', 'PECAM1', 'VWF', 'PROX1', 'COL1A1', 'BGN', 'DCN',
             'RGS5', 'MYH11', 'EPCAM', 'KRT7', 'KRT18', 'KRT5', 'MSLN', 'PMEL', 'MLANA', 'GFAP', 'PLP1', 'ALB', 'HBB', 'MKI67']
dg = [x for x in dot_genes if x in set(EXP.gene)]
pv = EXP[EXP.gene.isin(dg)].pivot(index='cluster', columns='gene', values='pct_expressing_measured').reindex(index=list(S.index), columns=dg)
mv = EXP[EXP.gene.isin(dg)].pivot(index='cluster', columns='gene', values='mean_cp10k_measured').reindex(index=list(S.index), columns=dg)
mz = np.log1p(mv) / np.log1p(mv).max().replace(0, np.nan)
fig, ax = plt.subplots(figsize=(0.22 * len(dg) + 3, 0.2 * len(S) + 2))
yy, xx = np.meshgrid(range(len(S)), range(len(dg)), indexing='ij')
ax.scatter(xx.ravel(), yy.ravel(), s=pv.values.ravel() * 40, c=mz.values.ravel(), cmap='Reds', vmin=0, vmax=1, linewidths=0)
ax.set_xticks(range(len(dg))); ax.set_xticklabels(dg, rotation=90, fontsize=6); ax.set_yticks(range(len(S)))
ax.set_yticklabels([f'{c} (n={S.loc[c, "cells"]:,}; unc {S.loc[c, "uncertain_share"]:.2f}; dbl z {S.loc[c, "doublet_score_z"]})' for c in S.index], fontsize=6)
ax.invert_yaxis(); fig.savefig(f'{out}/figures/dotplot_major.png', dpi=150, bbox_inches='tight'); plt.close(fig)

A.uns['params'] = json.dumps({'cancer': cancer, 'hvg': len(hvg), 'n_pc': n_pc, 'resolutions': RES, 'main': MAIN, 'adjacent_ARI': ari, 'contexts': ctx,
                              'regressed': list(STATES), 'harmony': 'source theta 2', 'seed': seed})
A.write_h5ad(f'{out}/object.h5ad')
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 80)
print(S[['cells', 'uncertain_share', 'global_c13_share', 'sources', 'median_genes', 'doublet_score_z', 'V7_parents_ge0.3', 'top_markers']].to_string(), flush=True)
log(step='done', clusters=len(S))
