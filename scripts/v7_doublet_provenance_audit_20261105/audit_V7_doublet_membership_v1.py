"""Read-only all-cell V7/QC membership and legacy-call audit. No matrix loading."""
import json,hashlib,re
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,h5py,zarr
from anndata.io import read_elem
r=PROJECT_ROOT;b=r/'runs/20260818_pancancer_5M_merge_prep';u=r/'runs/20260725_120105_unified_pancancer_covarnet';g=zarr.open_group(str(r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/pancancer_core_fullgene_v1.zarr'),mode='r')
index=pd.Index(read_elem(g['obs/_index']));cols=['gse_id','source_obs_name','cohort','doublet_score','doublet_method','doublet_status'];d=pd.DataFrame({c:read_elem(g['obs/'+c]) for c in cols},index=index)
frozen=pd.read_pickle(b/'11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1/checkpoints/final_annotation_v7.pkl.gz');assert frozen.index.equals(index) and len(d)==4676787
d['major_celltype_v7']=frozen.major_celltype_v7.to_numpy();d['cluster_uid']=frozen.cluster_uid.to_numpy();d['global_position']=np.arange(len(d));d['historical_predicted_doublet']=False;d['historical_doublet_removed']=False;d['matched_QC_filtered_object']=False
old=pd.read_csv(u/'06_full_qc_integration/phase06_production_qc_validation_33datasets_20260726_051121.tsv',sep='\t');new=pd.read_csv(b/'01_new_cohort_prep/manifests/dataset_qc_summary.tsv',sep='\t');records=[]
def resolve(x):
 p=Path(str(x))
 if p.is_absolute():return p
 if str(p).startswith('Project_v3/'):return r.parent/p
 return r/p
for gs,positions in d.groupby('gse_id',observed=True).indices.items():
 positions=np.asarray(positions);t=d.iloc[positions];co=t.cohort.astype(str).unique();assert len(co)==1
 if co[0]=='new_20260818_v1':row=new[new.gse_id==gs].iloc[0];p=resolve(row.output_path)
 else:row=old[old.gse_id==gs].iloc[0];p=resolve(row.h5ad_output)
 rec={'gse_id':str(gs),'cohort':co[0],'n_V7':len(t),'QC_path_relative':str(p.relative_to(r)) if p.is_relative_to(r) else p.name,'QC_exists':p.exists()}
 if not p.exists():records.append(rec);continue
 with h5py.File(p,'r') as f:
  ob=f['obs'];qi=pd.Index(read_elem(ob[ob.attrs['_index']])).astype(str);candidates={'obs_index':qi}
  for name in ['phase06_global_cell_id','original_obs_name']:
   if name in ob:candidates[name]=pd.Index(read_elem(ob[name])).astype(str)
  core_candidates={'canonical_index':t.index.astype(str),'source_obs_name':pd.Index(t.source_obs_name.astype(str))};joins=[]
  for qname,qindex in candidates.items():
   if not qindex.is_unique:continue
   for cname,cindex in core_candidates.items():
    match=qindex.get_indexer(cindex);joins.append((int((match>=0).sum()),qname,cname,match))
  best=max(joins,key=lambda x:x[0]);matched,qname,cname,match=best;rec.update(n_QC_filtered=len(qi),n_matched=int(matched),QC_join_key=qname,V7_join_key=cname,all_V7_cells_in_QC_filtered=(matched==len(t)))
  if matched!=len(t):records.append(rec);continue
  flags={}
  for name in ['predicted_doublet','doublet_removed']:
   if name in ob:
    flags[name]=pd.Series(read_elem(ob[name])).astype(str).str.lower().isin(['true','1']).to_numpy();rec['n_QC_filtered_'+name]=int(flags[name].sum());rec['n_V7_'+name]=int(flags[name][match].sum());d.iloc[positions,d.columns.get_loc('historical_'+name)]=flags[name][match]
  if 'filter_reason' in ob:rec['V7_original_filter_reasons']=pd.Series(np.asarray(read_elem(ob['filter_reason']))[match]).astype(str).value_counts().to_dict()
  if 'doublet_score' in ob:
   original=pd.to_numeric(pd.Series(read_elem(ob['doublet_score'])),errors='coerce').to_numpy()[match];current=pd.to_numeric(t.doublet_score,errors='coerce').to_numpy();rec['n_score_exact_or_both_missing']=int(((original==current)|(np.isnan(original)&np.isnan(current))).sum());rec['n_finite_original_scores']=int(np.isfinite(original).sum())
  rec['QC_header_fields']=[c for c in ['predicted_doublet','doublet_removed','filter_reason','doublet_score'] if c in ob]
  d.iloc[positions,d.columns.get_loc('matched_QC_filtered_object')]=True
 records.append(rec)
assert int(d.matched_QC_filtered_object.sum())==len(d),records
assert not d.historical_doublet_removed.any()
# Aggregate only; do not export cell or sample IDs.
tables={}
for name,keys in [('V7_historical_doublet_by_cohort.tsv',['cohort','doublet_method','doublet_status','historical_predicted_doublet']),('V7_historical_doublet_by_major.tsv',['major_celltype_v7','doublet_status','historical_predicted_doublet']),('V7_historical_doublet_by_GSE.tsv',['gse_id','cohort','doublet_status','historical_predicted_doublet'])]:tables[name]=d.groupby(keys,observed=True,dropna=False).size().reset_index(name='n_cells').to_json(orient='records')
ref=r/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';T=pd.read_pickle(ref/'stage_n_T6_measurement_scope_v1/T6_hierarchical_measured_review_v2_PRIVATE.pkl.gz');assert T.index.isin(d.index).all();k9=T[T['T6_review_cluster_source_res0.4'].eq(9)];joint=k9[k9.joint_competitor_doublet_candidate];assert len(k9)==1027 and len(joint)==245
target=[]
for name,t in [('all_T6',T),('K09_all',k9),('K09_joint245',joint),('all_Mast',frozen[frozen.major_celltype_v7.eq('Mast')])]:
 z=d.loc[t.index];target.append({'scope':name,'n_cells':len(z),'n_original_predicted_doublet':int(z.historical_predicted_doublet.sum()),'n_original_doublet_removed':int(z.historical_doublet_removed.sum()),'original_status_counts':z.doublet_status.astype(str).value_counts().to_dict(),'original_method_counts':z.doublet_method.astype(str).value_counts().to_dict()})
out={'status':'PASS_full_V7_membership_and_legacy_doublet_call_audit','observed_UTC':datetime.now(timezone.utc).isoformat(),'n_V7':len(d),'frozen_index_exact_with_core':True,'n_sources':len(records),'n_all_cells_matched_QC_filtered':int(d.matched_QC_filtered_object.sum()),'n_historical_predicted_doublets_retained':int(d.historical_predicted_doublet.sum()),'n_historical_removed_cells_reintroduced':int(d.historical_doublet_removed.sum()),'n_not_run_status':int(d.doublet_status.astype(str).eq('not_run').sum()),'source_records':records,'targets':target,'tables':tables,'scope_limits':['Doublet_status is execution provenance, not a singlet/doublet call; recovered calls come from original filtered-object predicted_doublet.','Full membership in QC-filtered objects verifies no historically removed cells reentered this V7/core axis.','No new doublet model or cell exclusion was run.']}
print(json.dumps(out,default=str))
