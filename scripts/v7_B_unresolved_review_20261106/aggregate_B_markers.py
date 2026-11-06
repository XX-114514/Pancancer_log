import json
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad,h5py,scipy.sparse as sp
from anndata.io import read_elem

r=PROJECT_ROOT/'runs/20261106_v7_B_unresolved_review_v1';name='B_unresolved';d=r/'stage_d_B_unresolved_refit';o=r/'stage_e_B_marker_review'
assert json.loads((r/'launch/B_unresolved_refit_terminal.json').read_text())['exit_code']==0;o.mkdir(exist_ok=False)
p=json.loads((r/'review_policy_v1.json').read_text())
extra=['PCSK1N','CPE','CHGA','CHGB','NEFL','NEFM','MAPT','NCAM1','MGP','ALB','ASGR1','APOA1','HPX','TF','PRSS1','PRSS2','CPA1','CTRB2','CTRC','CEL','FOXJ1','RSPH4A','MPO','ELANE','PRTN3','CTSG','AZU1','MMP9','ACP5','CTSK','SIGLEC15','MCM4','PCNA','VPREB1','VPREB3','IGLL1','IL10','TGFB1','IL21R','CXCR5','CD40LG','MYC','BCL2']
genes=list(dict.fromkeys(p['genes']+extra));a=ad.read_h5ad(r/'stage_c_QC_views'/name/(name+'_clean_model_raw.h5ad'))
md=pd.read_csv(d/'metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id');assert a.obs_names.equals(md.index)
q=pd.read_pickle(r/'stage_c_QC_views/QC_decisions_PRIVATE.pkl.gz').loc[md.index];assert q.current_lineage_model_eligible.all()
cl=md['source_snn_res.0.4'].astype(int).to_numpy();nc=int(cl.max()+1);assert set(cl)==set(range(nc));target=q.is_annotation_target.to_numpy().astype(bool)
ii=a.var_names.get_indexer(genes);found=ii>=0;raw=np.zeros((a.n_obs,len(genes)),np.int32);raw[:,found]=a.X[:,ii[found]].toarray();total=np.asarray(a.X.sum(axis=1)).ravel();cp=raw*(10000.0/np.maximum(total,1)[:,None])
gses=sorted(q.gse_id.astype(str).unique());available=[]
for gs in gses:
 with h5py.File(PROJECT_ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad'),'r') as f:v=pd.Index(read_elem(f['var']).index.astype(str))
 available.append(pd.Index(genes).isin(v))
av=np.asarray(available);gi=pd.Index(gses).get_indexer(q.gse_id.astype(str));assert not ((raw>0)&~av[gi]).any()
for population,take in [('all_model',np.ones(a.n_obs,bool)),('annotation_targets',target)]:
 ids=np.flatnonzero(take);cc=cl[take];n=np.bincount(cc,minlength=nc);H=sp.csr_matrix((np.ones(len(ids),np.int32),(cc,np.arange(len(ids)))),shape=(nc,len(ids)));S=sp.csr_matrix((np.ones(len(ids),np.int32),(cc,gi[take])),shape=(nc,len(gses)))
 den=np.asarray(S@av.astype(np.int32));nz=np.asarray(H@(raw[take]>0).astype(np.int32));sums=np.asarray(H@cp[take])
 rows=[{'cluster':c,'gene':g,'n_cells':int(n[c]),'n_measured':int(den[c,j]),'coverage_fraction':float(den[c,j]/n[c]),'n_expressing':int(nz[c,j]),'pct_expressing_measured':float(nz[c,j]/den[c,j]) if den[c,j] else np.nan,'mean_cp10k_measured':float(sums[c,j]/den[c,j]) if den[c,j] else np.nan,'population':population} for c in range(nc) if n[c] for j,g in enumerate(genes)]
 pd.DataFrame(rows).to_csv(o/(population+'_marker_expression.tsv'),sep='\t',index=False,na_rep='NA')
keys=['T_positive_source_aware','NK_receptor_effector_positive','B_distinct_components_positive','ASC_distinct_components_positive','Myeloid_mature_components_positive','Endothelial_distinct_components_positive','Epithelial_distinct_components_positive','Fibroblast_distinct_components_positive','Ery_structural_synthesis_nonHB_positive','Platelet_secreted_structural_positive','HSPC_hematopoietic_components_positive']
summaries=[];modules=[];refs=[];origins=[]
for c in range(nc):
 t=q.iloc[np.flatnonzero(cl==c)];tt=t[t.is_annotation_target]
 z={'cluster':c,'n_model_cells':len(t),'n_target_cells':len(tt),'n_context_controls':int((~t.is_annotation_target).sum()),'n_sources':t.gse_id.nunique(),'n_target_sources':tt.gse_id.nunique(),'max_source_fraction':float(t.gse_id.astype(str).value_counts().iloc[0]/len(t)),'max_target_source_fraction':float(tt.gse_id.astype(str).value_counts().iloc[0]/len(tt)) if len(tt) else np.nan,'raw_library_median':float(t.full_raw_library_total.median()),'target_raw_library_median':float(tt.full_raw_library_total.median()) if len(tt) else np.nan,'decontx_median':float(pd.to_numeric(t.decontx_contamination,errors='coerce').median()),'target_decontx_median':float(pd.to_numeric(tt.decontx_contamination,errors='coerce').median()) if len(tt) else np.nan,'original_doublet_percentile_median':float(pd.to_numeric(t.doublet_within_sample_method_percentile,errors='coerce').median()),'dominant_old_label':str(tt.final_annotation_v7.value_counts().index[0]) if len(tt) else 'controls_only'}
 for key in keys:z['n_model_'+key]=int(t[key].sum());z['n_target_'+key]=int(tt[key].sum())
 summaries.append(z)
 for pop,x in [('all_model',t),('annotation_targets',tt)]:
  if not len(x):continue
  for panel in p['panels']:modules.append({'cluster':c,'population':pop,'module':panel['id'],'axis':panel['axis'],'n_cells':len(x),'n_component_evaluable':int(x[panel['id']+'_core_evaluable'].sum()),'n_component_positive_raw1':int(x[panel['id']+'_core_positive_raw1'].sum()),'n_component_positive_raw2':int(x[panel['id']+'_core_positive_raw2'].sum())})
  for col in ['ev_celltypist_label','scatomic_full_label']:
   if col in x:
    for label,num in x[col].astype(object).fillna('missing').astype(str).value_counts().head(5).items():refs.append({'cluster':c,'population':pop,'reference_field':col,'label':label,'n_cells':int(num),'n_cluster':len(x)})
 for label,num in t.final_annotation_v7.value_counts().items():origins.append({'cluster':c,'old_label':str(label),'n_cells':int(num)})
pd.DataFrame(summaries).to_csv(o/'cluster_marker_QC_summary.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(modules).to_csv(o/'module_component_evidence.tsv',sep='\t',index=False);pd.DataFrame(refs).to_csv(o/'cluster_reference_top5.tsv',sep='\t',index=False);pd.DataFrame(origins).to_csv(o/'old_label_by_cluster.tsv',sep='\t',index=False)
de=pd.read_csv(d/'RNA_FindAllMarkers_common_measured.tsv',sep='\t');de=de[(de.avg_log2FC>0)&(de.p_val_adj<.05)].sort_values(['cluster','avg_log2FC'],ascending=[True,False]).groupby('cluster').head(15);de.to_csv(o/'top15_RNA_DEG.tsv',sep='\t',index=False,na_rep='NA')
out={'status':'B_marker_aggregation_complete_not_identity_acceptance','UTC':datetime.now(timezone.utc).isoformat(),'n_model_cells':a.n_obs,'n_annotation_targets':int(target.sum()),'n_context_controls':int((~target).sum()),'n_diagnostic_clusters':nc,'n_clusters_with_targets':int(np.unique(cl[target]).size),'n_marker_genes':len(genes),'all_source_measured_denominators_verified':True,'controls_not_counted_as_annotation_targets':True,'canonical_V7_changed':False,'fine_identity_accepted':False};(o/'complete.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
