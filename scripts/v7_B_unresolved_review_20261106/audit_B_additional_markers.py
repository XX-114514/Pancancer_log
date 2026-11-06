import json
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad,h5py
from anndata.io import read_elem
from sklearn.metrics import adjusted_rand_score

r=PROJECT_ROOT/'runs/20261106_v7_B_unresolved_review_v1';e=r/'stage_e2_B_QC_context'
meta=json.loads((e/'complete.json').read_text())
a=ad.read_h5ad(r/'stage_c_QC_views/B_unresolved/B_unresolved_clean_model_raw.h5ad')
md=pd.read_csv(r/'stage_d_B_unresolved_refit/metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id')
assert a.obs_names.equals(md.index)
q=pd.read_pickle(r/'stage_c_QC_views/QC_decisions_PRIVATE.pkl.gz').loc[md.index]
cl=md['source_snn_res.0.4'].astype(int);total=np.asarray(a.X.sum(axis=1)).ravel();available={}
for gs in md.gse_id.unique():
 with h5py.File(PROJECT_ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(str(gs)+'_symbol_counts.h5ad'),'r') as f:available[gs]=set(read_elem(f['var']).index.astype(str))
a.X=a.X.tocsc();errmax=0.;records=0
main=set(pd.read_csv(r/'stage_e_B_marker_review/all_model_marker_expression.tsv',sep='\t').gene)
for name,subset in [('all_model',np.ones(len(md),bool)),('annotation_targets',q.is_annotation_target.to_numpy())]:
 m=pd.read_csv(e/(name+'_additional_marker_expression.tsv'),sep='\t');assert not m.duplicated(['cluster','gene']).any();assert set(m.gene).isdisjoint(main);assert set(m.gene)==set(meta['additional_context_markers']);records+=len(m)
 for gene,t in m.groupby('gene'):
  ix=a.var_names.get_indexer([gene])[0];x=a.X[:,ix].toarray().ravel() if ix>=0 else np.zeros(len(md));measured=md.gse_id.map(lambda g:gene in available[g]).to_numpy();assert not ((x>0)&~measured).any()
  for z in t.to_dict('records'):
   take=cl.eq(z['cluster']).to_numpy()&subset;den=int((take&measured).sum());nz=int((take&(x>0)).sum());assert take.sum()==z['n_cells'] and den==z['n_measured'] and nz==z['n_expressing'];assert abs(den/take.sum()-z['coverage_fraction'])<1e-12
   if den:
    mean=float((10000.0*x[take]/np.maximum(total[take],1)).sum()/den);err=abs(mean-z['mean_cp10k_measured']);errmax=max(errmax,err);assert err<1e-7 and abs(nz/den-z['pct_expressing_measured'])<1e-12
metrics=pd.read_csv(e/'cluster_original_QC_metric_context.tsv',sep='\t')
for z in metrics.to_dict('records'):
 take=cl.eq(z['cluster']).to_numpy()&(q.is_annotation_target.to_numpy() if z['population']=='annotation_targets' else True);t=q.loc[take];assert len(t)==z['n_cells'] and t.gse_id.nunique()==z['n_sources']
 for col in ['n_genes_by_counts','pct_counts_mt','full_raw_library_total','decontx_contamination','doublet_within_sample_method_percentile']:
  key=col+'_median'
  if key in z:
   vals=pd.to_numeric(t[col],errors='coerce');assert vals.notna().sum()==z[col+'_n_available'];expected=vals.median();assert (pd.isna(expected) and pd.isna(z[key])) or abs(expected-z[key])<1e-8
raw=md['raw_snn_res.0.4'].astype(int);target=q.is_annotation_target.to_numpy()
for key,value in [('ARI_raw_source_graph_all_model',adjusted_rand_score(raw,cl)),('ARI_raw_source_graph_targets',adjusted_rand_score(raw[target],cl[target])),('ARI_source_clusters_vs_source_labels',adjusted_rand_score(cl,md.gse_id.astype(str)))]:assert abs(value-meta[key])<1e-12
out={'status':'PASS','UTC':datetime.now(timezone.utc).isoformat(),'n_model_cells':len(md),'n_marker_records':records,'n_additional_markers':len(meta['additional_context_markers']),'max_absolute_CP10K_error':errmax,'source_measured_denominators_verified':True,'original_QC_summary_medians_verified':True,'graph_context_ARI_independently_verified':True}
(e/'independent_additional_marker_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
