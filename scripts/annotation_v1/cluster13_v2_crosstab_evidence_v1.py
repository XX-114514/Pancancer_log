import sys, numpy as np, pandas as pd, anndata as ad, scipy.sparse as sp
S, sub, out = sys.argv[1:4]
pd.set_option('display.width',250,'display.max_rows',200,'display.max_columns',40)
m = pd.read_pickle(S+'/m.pkl'); d = pd.read_csv(S+'/draft.tsv', sep='\t')
m['sub_r1'] = m.sub_r1.astype(int)
m = m.merge(d[['sub_r1','manual_label_draft','manual_class','residual_review']], on='sub_r1', validate='many_to_one')
def fate(r):
    a = r.major_label
    if not a.startswith('r3_'): return 'parenchymal' if r.group == 'parenchymal' else 'removed_or_set_aside'
    if 'keep_flag_low_depth' in (r.r1_action, r.r2_action): return 'retained_low_depth'
    if r.r2_action in ('keep_flag_sample_specific','keep_flag_source_restricted'): return 'retained_flagged'
    return 'retained_main_like'
m['fate'] = m.apply(fate, axis=1)
print(pd.crosstab(m.residual_review, m.fate, margins=True))
print(pd.crosstab(m.manual_class, m.fate, margins=True))
r = m[m.fate.str.startswith('retained')].copy()
A = ad.read_h5ad(sub, backed='r')
idx = pd.Index(A.obs.cell_id.astype(str)).get_indexer(r.cell_id); assert (idx >= 0).all()
genes = ['PTPRC','CD3E','CD3D','NKG7','MS4A1','CD79A','JCHAIN','MZB1','LYZ','CD14','CD68','COL1A1','COL1A2','DCN','RGS5','PECAM1','VWF','EPCAM','KRT8','KRT18','KRT19','KRT5','MALAT1']
vi = A.var_names.get_indexer([g for g in genes]); present = [g for g, i in zip(genes, vi) if i >= 0]; vi = vi[vi >= 0]
print('X', type(A.X), A.uns.get('log1p') is not None if 'uns' in dir(A) else '')
order = np.argsort(idx); X = A.X[idx[order]]; X = sp.csr_matrix(X)[:, vi]; 
p = np.empty_like(order); p[order] = np.arange(len(order)); X = X[p]
F = pd.DataFrame((X > 0).toarray(), columns=present); F['fate'] = r.fate.values; F['tier'] = r.residual_review.values; F['lin'] = r.major_label.values; F['sub_r1'] = r.sub_r1.values
o = A.obs.iloc[idx]; F['cnv_corr'] = pd.to_numeric(o.malig_cnv_corr.values, errors='coerce'); F['sca'] = o.malig_scATOMIC_tumour.astype(str).isin(['True','1','tumour']).values
epi = F[['EPCAM','KRT8','KRT18','KRT19']].sum(1) >= 2
LIN = {'r3_T_NK':['PTPRC','CD3E','CD3D','NKG7'],'r3_Myeloid':['PTPRC','LYZ','CD14','CD68'],'r3_B_Plasma':['PTPRC','MS4A1','CD79A','JCHAIN','MZB1'],'r3_Stromal':['COL1A1','COL1A2','DCN','RGS5'],'r3_Endothelial':['PECAM1','VWF'],'r3_Microglia':['PTPRC','CD68']}
F['lin_marker'] = [bool(F.loc[i, LIN[l]].any()) for i, l in zip(F.index, F.lin)]
F['epi2'] = epi
print(F.groupby('tier').agg(n=('lin','size'), lin_marker=('lin_marker','mean'), epi2=('epi2','mean'), cnv=('cnv_corr','median')).round(3))
print(F.groupby(['tier','lin']).agg(n=('lin','size'), lin_marker=('lin_marker','mean'), epi2=('epi2','mean'), cnv=('cnv_corr','median')).round(3))
suspect = F[F.epi2 & ~F.lin_marker]; print('epi2 and no lineage marker', len(suspect), suspect.tier.value_counts().to_dict(), suspect.lin.value_counts().to_dict())
g = F.groupby('sub_r1').agg(retained=('lin','size'), lin_marker=('lin_marker','mean'), epi2=('epi2','mean')).round(3).reset_index()
t = pd.crosstab(m.sub_r1, m.fate).reset_index()
t2 = pd.crosstab(m.sub_r1, m.major_label).reset_index()
o2 = d[['sub_r1','cells','non_malignant','manual_label_draft','manual_class','residual_review']].merge(t, on='sub_r1').merge(g, on='sub_r1', how='left').sort_values('non_malignant', ascending=False)
o2.to_csv(out+'/c13_subcluster_by_v2_fate.tsv', sep='\t', index=False)
t2.to_csv(out+'/c13_subcluster_by_v2_label.tsv', sep='\t', index=False)
F[F.tier!='not_flagged'].groupby(['tier','lin','fate']).agg(cells=('lin','size'), lineage_marker_share=('lin_marker','mean'), epithelial_2of4_share=('epi2','mean'), cnv_corr_median=('cnv_corr','median')).round(3).reset_index().to_csv(out+'/c13_retained_lineage_evidence.tsv', sep='\t', index=False)
print(o2[o2.retained>=50][['sub_r1','non_malignant','manual_label_draft','residual_review','retained','lin_marker','epi2']].to_string(index=False))
print(F[['PTPRC','EPCAM','KRT18']].mean().round(3).to_dict())
