import json
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad,h5py,scipy.sparse as sp
from anndata.io import read_elem
from sklearn.metrics import adjusted_rand_score

r=PROJECT_ROOT/'runs/20261106_v7_B_unresolved_review_v1'
assert json.loads((r/'launch/B_unresolved_refit_terminal.json').read_text())['exit_code']==0
assert (r/'stage_e_B_marker_review/complete.json').exists()
o=r/'stage_e2_B_QC_context';o.mkdir(exist_ok=False)
md=pd.read_csv(r/'stage_d_B_unresolved_refit/metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id');q=pd.read_pickle(r/'stage_c_QC_views/QC_decisions_PRIVATE.pkl.gz').loc[md.index];cl=md['source_snn_res.0.4'].astype(int)
rows=[]
for c in sorted(cl.unique()):
 t=q.loc[cl.eq(c)]
 for population,x in [('all_model',t),('annotation_targets',t[t.is_annotation_target])]:
  if not len(x):continue
  z={'cluster':int(c),'population':population,'n_cells':len(x),'n_sources':x.gse_id.nunique(),'max_source_fraction':float(x.gse_id.astype(str).value_counts().iloc[0]/len(x))}
  for col in ['n_genes_by_counts','pct_counts_mt','full_raw_library_total','decontx_contamination','doublet_within_sample_method_percentile']:
   if col in x:
    v=pd.to_numeric(x[col],errors='coerce');z[col+'_median']=float(v.median());z[col+'_q25']=float(v.quantile(.25));z[col+'_q75']=float(v.quantile(.75));z[col+'_n_available']=int(v.notna().sum())
  rows.append(z)
pd.DataFrame(rows).to_csv(o/'cluster_original_QC_metric_context.tsv',sep='\t',index=False,na_rep='NA')
raw_cl=md['raw_snn_res.0.4'].astype(int);target=q.is_annotation_target.to_numpy()
existing=set(pd.read_csv(r/'stage_e_B_marker_review/all_model_marker_expression.tsv',sep='\t').gene)
new_context=['LAMP3','CCR7','FSCN1','CLEC10A','CADM1','CD1E','CD68','C1QB','C1QC','APOE','APOC1','SPP1','FCN1','S100A12','AIF1','GPNMB','BATF3','CD80','CENPF','FCER2','FOS','HSP90B1','JUN']
genes=[g for g in new_context if g not in existing]
a=ad.read_h5ad(r/'stage_c_QC_views/B_unresolved/B_unresolved_clean_model_raw.h5ad');assert a.obs_names.equals(md.index)
ii=a.var_names.get_indexer(genes);found=ii>=0;raw=np.zeros((a.n_obs,len(genes)),np.int32);raw[:,found]=a.X[:,ii[found]].toarray();total=np.asarray(a.X.sum(axis=1)).ravel();assert np.array_equal(total,q.full_raw_library_total.to_numpy());cp=raw*(10000.0/np.maximum(total,1)[:,None]);cc=cl.to_numpy();nc=int(cc.max()+1)
gses=sorted(q.gse_id.astype(str).unique());av=[]
for gs in gses:
 with h5py.File(PROJECT_ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad'),'r') as f:v=pd.Index(read_elem(f['var']).index.astype(str))
 av.append(pd.Index(genes).isin(v))
av=np.asarray(av);gi=pd.Index(gses).get_indexer(q.gse_id.astype(str));assert not ((raw>0)&~av[gi]).any()
for pop,take in [('all_model',np.ones(a.n_obs,bool)),('annotation_targets',target)]:
 inds=np.flatnonzero(take);c=cc[take];n=np.bincount(c,minlength=nc);H=sp.csr_matrix((np.ones(len(inds),np.int32),(c,np.arange(len(inds)))),shape=(nc,len(inds)));S=sp.csr_matrix((np.ones(len(inds),np.int32),(c,gi[take])),shape=(nc,len(gses)));den=np.asarray(S@av.astype(np.int32));nz=np.asarray(H@(raw[take]>0).astype(np.int32));sums=np.asarray(H@cp[take])
 rows=[{'cluster':c,'gene':gene,'n_cells':int(n[c]),'n_measured':int(den[c,j]),'coverage_fraction':float(den[c,j]/n[c]),'n_expressing':int(nz[c,j]),'pct_expressing_measured':float(nz[c,j]/den[c,j]) if den[c,j] else np.nan,'mean_cp10k_measured':float(sums[c,j]/den[c,j]) if den[c,j] else np.nan,'population':pop} for c in range(nc) if n[c] for j,gene in enumerate(genes)]
 pd.DataFrame(rows).to_csv(o/(pop+'_additional_marker_expression.tsv'),sep='\t',index=False,na_rep='NA')
out={'status':'B_original_QC_and_graph_context_complete','UTC':datetime.now(timezone.utc).isoformat(),'n_model_cells':len(md),'n_targets':int(target.sum()),'additional_context_markers':genes,'additional_markers_not_duplicated_in_main_table':True,'ARI_raw_source_graph_all_model':float(adjusted_rand_score(raw_cl,cl)),'ARI_raw_source_graph_targets':float(adjusted_rand_score(raw_cl[target],cl[target])),'ARI_source_clusters_vs_source_labels':float(adjusted_rand_score(cl,q.gse_id.astype(str))),'interpretation':'Exploratory raw-versus-source-corrected graph context; not a subtype stability or optimal-integration validation. Gene-count and mitochondrial metrics are inherited original QC; full raw library totals were recalculated.'}
(o/'complete.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
