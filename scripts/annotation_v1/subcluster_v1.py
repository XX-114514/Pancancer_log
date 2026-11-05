# De novo annotation v1 (methods/tme-de-novo-annotation-v1-20261105.md): re-cluster a subset of TME v1 run v2
# (one or more global leiden_tme clusters) and build per-subcluster evidence for annotation.
#  read   : per dataset shard, h5py CSR in contiguous row blocks, only the selected rows and the common + panel genes
#  embed  : seurat_v3 HVG per source on counts (sources < MIN_SRC cells pooled), 3,000 non-technical; PCA 50 (centred);
#           Harmony source theta 2; 15-NN; Leiden 0.5 / 1.0 / 2.0; UMAP min_dist 0.3
#  evidence per subcluster (resolution 1.0): panel scores (mean log expression, z across subclusters), Wilcoxon DE
#           (max 2,000 cells per subcluster), cross-source reproducibility of the top DE genes, composition, malignancy
#           evidence, QC medians, and an automatic call following the method's decision order
# Outputs (OUT): subset.h5ad, subclusters.tsv, de_top50.tsv, panel_scores.tsv, gene_dot.tsv, figures/*.png
# Usage: subcluster_v1.py RUN_DIR CLUSTERS(comma) OUT_DIR PANEL_TSV SEED
import sys, os, re, json, datetime
import numpy as np, pandas as pd, scipy.sparse as sp, anndata as ad, scanpy as sc, harmonypy as hm, zarr, h5py
from anndata.io import read_elem

run, clusters, out, panelf, seed = sys.argv[1], sys.argv[2].split(','), sys.argv[3], sys.argv[4], int(sys.argv[5])
os.makedirs(f'{out}/figures', exist_ok=True)
log = lambda **k: print(json.dumps({'t': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds'), **k}, default=str), flush=True)
BLOCK, MIN_SRC, N_HVG, N_PC = 50_000, 200, 3000, 50
technical = lambda g: g.upper() == 'MALAT1' or bool(re.search(r'^(MT-|RP[SL][0-9]|HB[ABDEGQMZ][0-9]|ERCC-)', g, re.I))
panel = pd.read_csv(panelf, sep='\t')
panel['gl'] = panel.genes.str.split(',')
pgenes = sorted({g for gl in panel.gl for g in gl})

g = zarr.open(f'{run}/tme_v1.zarr', 'r')
obs_all = read_elem(g['obs'])
sel_all = obs_all['leiden_tme'].astype(str).isin(clusters).values
common = list(pd.read_csv(f'{run}/report/hvg_tme_v1.tsv', sep='\t', index_col=0).index)
genes = list(dict.fromkeys(common + pgenes))
log(step='select', clusters=clusters, cells=int(sel_all.sum()), genes=len(genes))


def read_rows_cols(grp, rows, cols, ncol_out):
    indptr = grp['indptr'][:]
    out_ = []
    for a in range(0, len(indptr) - 1, BLOCK):
        b = min(a + BLOCK, len(indptr) - 1)
        r = rows[(rows >= a) & (rows < b)]
        if not len(r):
            continue
        lo, hi = int(indptr[a]), int(indptr[b])
        blk = sp.csr_matrix((grp['data'][lo:hi], grp['indices'][lo:hi], indptr[a:b + 1] - lo), shape=(b - a, len(cols)))[r - a].tocoo()
        keep = cols[blk.col] >= 0
        out_.append(sp.csr_matrix((blk.data[keep].astype(np.float32), (blk.row[keep], cols[blk.col[keep]])), shape=(len(r), ncol_out)))
    return sp.vstack(out_).tocsr()


src = obs_all['source'].astype(str).values
X, C, O, measured = [], [], [], {}
for s in dict.fromkeys(src):
    pos = np.flatnonzero(src == s)
    w = np.flatnonzero(sel_all[pos])
    if not len(w):
        continue
    with h5py.File(f'{run}/shards/{s}.h5ad', 'r') as f:
        vn = f['var'][f['var'].attrs['_index']].asstr()[:]
        assert f['X'].attrs['shape'][0] == len(pos), s
        colmap = pd.Index(genes).get_indexer(vn)
        measured[s] = np.isin(genes, vn)
        X.append(read_rows_cols(f['X'], w, colmap, len(genes)))
        C.append(read_rows_cols(f['layers/counts'], w, colmap, len(genes)))
    O.append(obs_all.iloc[pos[w]])
A = ad.AnnData(X=sp.vstack(X).tocsr(), obs=pd.concat(O), var=pd.DataFrame(index=genes))
A.layers['counts'] = sp.vstack(C).tocsr(); del X, C  # obs (incl. malig_* evidence) comes from the zarr
M =pd.DataFrame(measured, index=genes).T  # source x gene measured
log(step='loaded', shape=A.shape)

# embedding
A.obs['hvg_batch'] = A.obs['source'].astype(str)
small = A.obs['hvg_batch'].value_counts(); small = small.index[small < MIN_SRC]
A.obs.loc[A.obs['hvg_batch'].isin(small), 'hvg_batch'] = 'pooled_small_sources'
cg = np.isin(genes, common) & ~np.array([technical(x) for x in genes])
B = ad.AnnData(X=A.layers['counts'][:, cg], obs=A.obs[['hvg_batch']].copy(), var=pd.DataFrame(index=np.array(genes)[cg]))
sc.pp.highly_variable_genes(B, n_top_genes=N_HVG, flavor='seurat_v3', batch_key='hvg_batch')
hvg = list(B.var_names[B.var['highly_variable']]); del B
A.var['hvg_sub'] = A.var_names.isin(hvg)
sc.pp.pca(A, n_comps=N_PC, mask_var='hvg_sub', zero_center=True, svd_solver='arpack', random_state=seed)
h = hm.run_harmony(np.asarray(A.obsm['X_pca'], np.float64), A.obs[['source']].astype(str), vars_use=['source'], theta=2.0, lamb=1.0,
                   max_iter_harmony=20, random_state=seed)
Z = np.asarray(h.Z_corr); A.obsm['X_harmony'] = (Z.T if Z.shape[0] == N_PC else Z).astype(np.float32)
sc.pp.neighbors(A, n_neighbors=15, use_rep='X_harmony', random_state=seed)
for r in (0.5, 1.0, 2.0):
    sc.tl.leiden(A, resolution=r, key_added=f'sub_r{r:g}', random_state=seed, flavor='igraph', n_iterations=2, directed=False)
sc.tl.umap(A, min_dist=0.3, random_state=seed)
K = 'sub_r1'
log(step='embedded', hvg=len(hvg), **{f'r{r:g}': int(A.obs[f'sub_r{r:g}'].nunique()) for r in (0.5, 1.0, 2.0)})

# panel scores (per cell mean over measured panel genes), per subcluster z
cl = A.obs[K].astype(str).values
gi = {x: i for i, x in enumerate(genes)}
P = {}
for _, row in panel.iterrows():
    idx = [gi[x] for x in row.gl if x in gi]
    P[f'{row.tier}:{row.label}'] = np.asarray(A.X[:, idx].mean(1)).ravel()
P = pd.DataFrame(P, index=A.obs_names)
PS = P.groupby(cl).mean()
PZ = (PS - PS.mean()) / PS.std(ddof=0).replace(0, np.nan)
PS.add_prefix('mean|').join(PZ.add_prefix('z|')).to_csv(f'{out}/panel_scores.tsv', sep='\t')
# gene dot table (mean and fraction expressing) for panel genes
pg = [x for x in pgenes if x in gi]
E = A[:, pg].X
dot = []
for c in sorted(set(cl), key=int):
    m = cl == c; e = E[m]
    dot.append(pd.DataFrame({'subcluster': c, 'gene': pg, 'mean': np.asarray(e.mean(0)).ravel(), 'frac': np.asarray((e > 0).mean(0)).ravel()}))
pd.concat(dot).to_csv(f'{out}/gene_dot.tsv', sep='\t', index=False)

# DE (capped per subcluster)
rng = np.random.default_rng(seed)
take = np.concatenate([rng.choice(np.flatnonzero(cl == c), min(2000, (cl == c).sum()), replace=False) for c in sorted(set(cl), key=int)])
D = A[np.sort(take)].copy()
sc.tl.rank_genes_groups(D, groupby=K, method='wilcoxon', use_raw=False, pts=True)
de = sc.get.rank_genes_groups_df(D, group=None)
de = de[~de['names'].map(technical)]
de = de.groupby('group', observed=True).head(50)
de.to_csv(f'{out}/de_top50.tsv', sep='\t', index=False)
del D

# cross-source reproducibility of the top 10 DE genes per subcluster
srcs = A.obs['source'].astype(str).values
rep = {}
for c, dg in de.groupby('group', observed=True):
    top = [x for x in dg['names'].head(10) if x in gi]
    m = cl == str(c); ok, tot = 0, 0
    for s in np.unique(srcs[m]):
        ins, outs = m & (srcs == s), ~m & (srcs == s)
        if ins.sum() < 50 or outs.sum() < 50:
            continue
        ii = [gi[x] for x in top]
        lfc = np.asarray(A.X[ins][:, ii].mean(0)).ravel() - np.asarray(A.X[outs][:, ii].mean(0)).ravel()
        ok += int((lfc > 0.5).mean() >= 0.5); tot += 1
    rep[str(c)] = (ok, tot)

# per-subcluster summary + automatic call
o = A.obs
grp = o.groupby(K, observed=True)
S = pd.DataFrame({'cells': grp.size(), 'non_malignant': grp['tme_flag'].agg(lambda x: int((x.astype(str) == 'non_malignant').sum())),
                  'uncertain_share': grp['tme_flag'].agg(lambda x: (x.astype(str) != 'non_malignant').mean()).round(3),
                  'sources': grp['source'].nunique(), 'units': grp['unit_id'].nunique(),
                  'top_source': grp['source'].agg(lambda x: x.value_counts().index[0]),
                  'top_source_share': grp['source'].agg(lambda x: x.value_counts(normalize=True).iloc[0]).round(3),
                  'top_cancers': grp['analysis_main'].agg(lambda x: '; '.join(f'{k} {v:.2f}' for k, v in x.value_counts(normalize=True).head(3).items())),
                  'cnv_score_median': grp['malig_cnv_score'].median(), 'cnv_corr_median': grp['malig_cnv_corr'].median().round(3),
                  'scatomic_tumour_share': grp['malig_scATOMIC_tumour'].agg(lambda x: (x.astype(str) == 'tumour').mean()).round(3),
                  'scatomic_top': grp['malig_scATOMIC_pred'].agg(lambda x: '; '.join(f'{k} {v:.2f}' for k, v in x.astype(str).value_counts(normalize=True).head(3).items())),
                  'median_genes': grp['n_genes'].median(), 'median_umi': grp['total_counts'].median(),
                  'median_pct_mt': grp['pct_mt'].median().round(2), 'median_contamination': grp['decontX_contamination'].median().round(3)})
S.index = S.index.astype(str)
S['top_de'] = de.groupby('group', observed=True)['names'].agg(lambda x: ','.join(x.head(12))).reindex(S.index)
S['de_reproducible_sources'] = [f'{rep.get(c, (0, 0))[0]}/{rep.get(c, (0, 0))[1]}' for c in S.index]
L1 = PZ[[c for c in PZ.columns if c.startswith('L1:') and not c.endswith((':cycling', ':low_quality', ':erythroid', ':immune'))]]
L2 = PZ[[c for c in PZ.columns if c.startswith('L2:')]]
S['L1_top'] = L1.idxmax(1).str[3:]; S['L1_top_z'] = L1.max(1).round(2)
S['L1_second'] = L1.apply(lambda r: r.nlargest(2).index[-1][3:], axis=1); S['L1_second_z'] = L1.apply(lambda r: r.nlargest(2).iloc[-1], axis=1).round(2)
S['L2_top'] = L2.idxmax(1).str[3:]; S['L2_top_z'] = L2.max(1).round(2)
S['cycling_z'] = PZ['L1:cycling'].round(2); S['erythroid_z'] = PZ['L1:erythroid'].round(2)
lab_genes = {f'{r.tier}:{r.label}': set(r.gl) for r in panel.itertuples()}
de_sets = de.groupby('group', observed=True)['names'].agg(set)
de_sets.index = de_sets.index.astype(str)


def call(c, r):
    d = de_sets.get(c, set())
    if r.median_genes < 500 or (r.median_pct_mt > 20 and r.L1_top_z < 1):
        return 'low_quality', 'low'
    if r.erythroid_z >= 1.5 and r.L1_top_z < r.erythroid_z:
        return 'erythroid_contamination', 'low'
    if r.L1_top_z >= 1 and r.L1_second_z >= 1 and len(d & lab_genes['L1:' + r.L1_top]) >= 1 and len(d & lab_genes['L1:' + r.L1_second]) >= 1:
        return f'doublet_like_{r.L1_top}+{r.L1_second}', 'low'
    for key, z in (('L2:' + r.L2_top, r.L2_top_z), ('L1:' + r.L1_top, r.L1_top_z)):
        if z >= 1 and len(d & lab_genes[key]) >= 2:
            ok, tot = rep.get(c, (0, 0))
            conf = 'high' if tot and ok / tot >= 0.5 and r.sources >= 3 else 'medium'
            return key.split(':', 1)[1], conf
    for key, z in (('L2:' + r.L2_top, r.L2_top_z), ('L1:' + r.L1_top, r.L1_top_z)):
        if z >= 1 or len(d & lab_genes[key]) >= 2:
            return key.split(':', 1)[1], 'low'
    return f'unresolved_{c}', 'low'


calls = [call(c, r) for c, r in S.iterrows()]
S['auto_label'] = [x[0] for x in calls]; S['auto_confidence'] = [x[1] for x in calls]
S['source_restricted'] = (S.top_source_share >= 0.8) | (S.sources < 3)
S['residual_tumour_suspect'] = (S.uncertain_share >= 0.5) & ((S.cnv_corr_median > 0.3) | (S.scatomic_tumour_share >= 0.5))
S = S.sort_values('cells', ascending=False)
S.to_csv(f'{out}/subclusters.tsv', sep='\t')
A.obs['auto_label'] = A.obs[K].astype(str).map(S['auto_label'])

# figures
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
U = A.obsm['X_umap']; order = rng.permutation(A.n_obs)
pal = list(plt.get_cmap('tab20').colors) + list(plt.get_cmap('tab20b').colors) + list(plt.get_cmap('Set3').colors)
def scat(key, fn, top=25, cont=False):
    fig, ax = plt.subplots(figsize=(10, 8))
    if cont:
        v = A.obs[key].astype(float).values
        s_ = ax.scatter(U[order, 0], U[order, 1], c=v[order], s=0.3, cmap='viridis', linewidths=0, rasterized=True); fig.colorbar(s_, ax=ax, shrink=0.6)
    else:
        v = A.obs[key].astype(str).values; cats = pd.Series(v).value_counts().index[:top].tolist()
        v = np.where(np.isin(v, cats), v, 'other'); cats = cats + (['other'] if (v == 'other').any() else [])
        cm = {c: (pal[i % len(pal)] if c != 'other' else (.85, .85, .85)) for i, c in enumerate(cats)}
        ax.scatter(U[order, 0], U[order, 1], c=[cm[x] for x in v[order]], s=0.3, linewidths=0, rasterized=True)
        for c in cats: ax.scatter([], [], c=[cm[c]], s=15, label=f'{c} ({(v == c).sum():,})')
        ax.legend(loc='center left', bbox_to_anchor=(1, .5), fontsize=6, frameon=False)
        if key == K:
            for c in cats:
                if c != 'other':
                    xy = np.median(U[v == c], 0); ax.text(xy[0], xy[1], c, fontsize=7, weight='bold')
    ax.set_axis_off(); ax.set_title(f'{key} ({A.n_obs:,} cells; global clusters {",".join(clusters)})')
    fig.savefig(f'{out}/figures/{fn}.png', dpi=170, bbox_inches='tight'); plt.close(fig)
scat(K, 'umap_subcluster', top=60); scat('auto_label', 'umap_auto_label', top=40); scat('tme_flag', 'umap_tme_flag')
scat('source', 'umap_source'); scat('analysis_main', 'umap_cancer'); scat('malig_cnv_corr', 'umap_cnv_corr', cont=True)
scat('n_genes', 'umap_n_genes', cont=True)
# dot plot of L1 + epithelial L2 panel genes (first two genes per label)
sel = [g_ for r in panel[(panel.tier == 'L1') | (panel.compartment == 'epithelial')].itertuples() for g_ in r.gl[:2] if g_ in gi]
sel = list(dict.fromkeys(sel))
dd = pd.read_csv(f'{out}/gene_dot.tsv', sep='\t'); dd = dd[dd.gene.isin(sel)]
cls = [c for c in S.index]
fig, ax = plt.subplots(figsize=(0.22 * len(sel) + 3, 0.22 * len(cls) + 2))
mx = dd.pivot(index='subcluster', columns='gene', values='mean').reindex(index=[int(c) for c in cls], columns=sel)
fr = dd.pivot(index='subcluster', columns='gene', values='frac').reindex(index=[int(c) for c in cls], columns=sel)
mz = mx / mx.max().replace(0, np.nan)
yy, xx = np.meshgrid(range(len(cls)), range(len(sel)), indexing='ij')
ax.scatter(xx.ravel(), yy.ravel(), s=(fr.values.ravel() * 40), c=mz.values.ravel(), cmap='Reds', vmin=0, vmax=1, linewidths=0)
ax.set_xticks(range(len(sel))); ax.set_xticklabels(sel, rotation=90, fontsize=6)
ax.set_yticks(range(len(cls))); ax.set_yticklabels([f'{c} {S.loc[c, "auto_label"]} (n={S.loc[c, "cells"]:,}, unc {S.loc[c, "uncertain_share"]:.2f})' for c in cls], fontsize=6)
ax.invert_yaxis(); fig.savefig(f'{out}/figures/dotplot_panel.png', dpi=170, bbox_inches='tight'); plt.close(fig)

A.uns['subcluster_params'] = json.dumps({'clusters': clusters, 'hvg': N_HVG, 'pcs': N_PC, 'harmony': 'source theta 2', 'seed': seed})
M.to_csv(f'{out}/gene_measured_by_source.tsv.gz', sep='\t', compression='gzip')
A.write_h5ad(f'{out}/subset.h5ad')
pd.set_option('display.width', 250); pd.set_option('display.max_colwidth', 70)
print(S[['cells', 'non_malignant', 'uncertain_share', 'sources', 'top_source_share', 'auto_label', 'auto_confidence', 'L2_top', 'L2_top_z', 'cnv_corr_median',
         'median_genes', 'residual_tumour_suspect', 'top_de']].to_string(), flush=True)
log(step='done')
