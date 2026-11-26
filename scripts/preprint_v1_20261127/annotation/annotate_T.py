"""Phase A: identity clustering (stress/cycle/TCR/IG/sex genes removed), CellTypist, marker programmes.
Outputs per-cell table and per-cluster review table; final labels are assigned in phase B after review."""
import scanpy as sc, anndata as ad, numpy as np, pandas as pd, harmonypy, celltypist, sys, re, os
W = os.environ['TFINE_ROOT']
os.environ['CELLTYPIST_FOLDER'] = W + '/models'
A = ad.read_h5ad(W + '/T_cells_cap1000.h5ad')
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)
G = A.var_names
cc = ['MKI67','TOP2A','TYMS','STMN1','TUBA1B','TUBB','H2AZ1','HMGB2','PCNA','MCM2','MCM3','MCM5','MCM7','CDK1','CCNB1','CCNB2','CCNA2','UBE2C','BIRC5','CENPF','NUSAP1','TPX2','PTTG1','CDKN3','RRM2','ZWINT','CKS1B','CKS2','DUT','GINS2','CLSPN','HMGN2','TK1','KIAA0101','PCLAF','SMC4','NUCKS1','DEK','H2AZ2']
excl = re.compile(r'^(MT-|RP[SL]\d|RPLP|MRP[SL]|TR[ABDG][VJC]|IG[HKL][VCJDGAEM]|HSP|DNAJ|HIST|H1-|H2A|H2B|H3-|H4C|XIST|TSIX|RPS4Y|DDX3Y|UTY|KDM5D|EIF1AY|USP9Y|MALAT1|NEAT1|FOS|JUN|IER\d|EGR\d|DUSP1$|ZFP36|NR4A|ATF3|KLF[246]$|PPP1R15A|CDKN1A$|GADD45|HSPA|BAG3|SERPINH1|CHORDC1)')
bad = np.array([bool(excl.match(g)) or g in cc for g in G])
sc.pp.highly_variable_genes(A, n_top_genes=5000, flavor='seurat', batch_key='dataset', subset=False)
A.var.loc[bad, 'highly_variable'] = False
hv = A.var.highly_variable.values
print('HVG used', hv.sum(), flush=True)
B = A[:, hv].copy(); sc.pp.scale(B, max_value=10); sc.tl.pca(B, n_comps=40, random_state=0)
H = harmonypy.run_harmony(B.obsm['X_pca'], A.obs, ['dataset', 'sampleID'], theta=[2, 1], max_iter_harmony=20, random_state=0)
Zc = np.asarray(H.Z_corr); Zc = Zc.T if Zc.shape[0] != A.n_obs else Zc
A.obsm['X_pca_h'] = Zc.astype(np.float32); del B
sc.pp.neighbors(A, n_neighbors=30, use_rep='X_pca_h', random_state=0)
sc.tl.leiden(A, resolution=1.5, key_added='leiden_fine', flavor='igraph', n_iterations=2, random_state=0)
sc.tl.umap(A, random_state=0)
print('clusters', A.obs.leiden_fine.nunique(), flush=True)
# CellTypist on all genes, log1p CP10k; majority vote over our fine clusters
for m in ['Immune_All_Low.pkl', 'Immune_All_High.pkl']:
    p = celltypist.annotate(A, model=m, majority_voting=True, over_clustering='leiden_fine')
    k = m.split('.')[0].replace('Immune_All_', 'ct_')
    A.obs[k] = p.predicted_labels['predicted_labels'].values
    A.obs[k + '_mv'] = p.predicted_labels['majority_voting'].values
    A.obs[k + '_conf'] = p.probability_matrix.max(1).values
progs = {
 'CD4': ['CD4','CD40LG'], 'CD8': ['CD8A','CD8B'],
 'Treg': ['FOXP3','IL2RA','IKZF2','TNFRSF18','TNFRSF4','CTLA4'],
 'Tex': ['CXCL13','HAVCR2','PDCD1','LAG3','TOX','ENTPD1','LAYN','TIGIT'],
 'Tfh': ['CXCL13','CXCR5','ICOS','PDCD1','MAF','TOX2','SH2D1A','BTLA'],
 'Naive': ['CCR7','SELL','LEF1','TCF7','MAL','SOX4'],
 'Tcm_mem': ['IL7R','LTB','CD69','ANXA1','KLF2','S1PR1'],
 'GZMK_Tem': ['GZMK','CCL5','CST7','CD74','EOMES'],
 'TRM': ['ZNF683','ITGAE','CXCR6','ITGA1','CD69'],
 'Temra': ['FGFBP2','CX3CR1','GZMH','GNLY','FCGR3A','KLRG1'],
 'Cytotoxic': ['GZMB','PRF1','NKG7','GNLY','GZMA'],
 'MAIT': ['SLC4A10','TRAV1-2','KLRB1','ZBTB16','RORC'],
 'gdT': ['TRDC','TRGC1','TRGC2'],
 'NKlike': ['KLRD1','KLRF1','NCAM1','TYROBP','FCER1G'],
 'Th17': ['RORC','IL17A','CCR6','IL23R','KLRB1'],
 'Th1': ['TBX21','IFNG','CXCR3'],
 'ISG': ['ISG15','IFIT1','IFIT3','MX1','IFI6'],
 'Stress': ['HSPA1A','HSPA1B','DNAJB1','HSPH1','FOS','JUN'],
 'Cycling': ['MKI67','TOP2A','TYMS','STMN1','CENPF','UBE2C'],
}
for k, v in progs.items():
    v = [x for x in v if x in G]
    sc.tl.score_genes(A, v, score_name='sc_' + k, random_state=0)
keyg = sorted(set(sum(progs.values(), [])) | {'CD3E','CD3D','TRAC','CD8B','FOXP3','GZMB','CTLA4','SELL','CCR7','IL7R'})
keyg = [x for x in keyg if x in G]
E = A[:, keyg].X; E = E.toarray() if hasattr(E, 'toarray') else E
pos = pd.DataFrame(E > 0, columns=keyg, index=A.obs_names)
cl = A.obs.leiden_fine
summ = pos.groupby(cl).mean().round(2)
scores = A.obs[[c for c in A.obs if c.startswith('sc_')]].groupby(cl).mean().round(2)
size = cl.value_counts().rename('n')
nds = A.obs.groupby('leiden_fine').dataset.agg(lambda x: x.value_counts(normalize=True).round(2).to_dict()).rename('dataset_mix')
oldn = A.obs.groupby('leiden_fine').subcluster.agg(lambda x: x.value_counts(normalize=True).head(3).round(2).to_dict()).rename('old_nodes')
ctl = A.obs.groupby('leiden_fine').ct_Low.agg(lambda x: x.value_counts(normalize=True).head(3).round(2).to_dict()).rename('celltypist_low_top3')
pd.concat([size, ctl, oldn, nds, scores, summ], axis=1).to_csv(W + '/cluster_review_phaseA.tsv', sep='\t')
A.obs.to_csv(W + '/cells_phaseA.tsv.gz', sep='\t')
np.save(W + '/umap_phaseA.npy', A.obsm['X_umap'])
A.write_h5ad(W + '/T_cells_phaseA.h5ad', compression='gzip')
print('phase A done')
