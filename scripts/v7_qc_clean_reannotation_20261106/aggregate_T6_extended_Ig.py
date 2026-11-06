import os
from pathlib import Path
PROJECT_ROOT=Path(os.environ["PROJECT_ROOT"])
import os,json,zlib,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad,h5py,scipy.sparse as sp
from anndata.io import read_elem
run=PROJECT_ROOT/'runs/20261106_v7_QC_clean_reannotation_v1';name=os.environ['REVIEW_LINEAGE'];d=run/('stage_d_'+name+'_clean_refit');assert json.loads((run/'launch'/(name+'_refit_terminal.json')).read_text())['exit_code']==0 and (d/'complete.json').exists()
o=run/('stage_l_'+name+'_extended_Ig_marker_review_v1');o.mkdir(exist_ok=False)
spec=json.loads((run/'config'/(name+'_context_panels_v2.json')).read_text());genes=list(dict.fromkeys(g for p in spec['panels']+spec.get('reference_panels',spec.get('panMC_reference_panels',[])) for g in p['genes']))
a=ad.read_h5ad(run/'stage_c_clean_views'/name/(name+'_clean_raw.h5ad'));md=pd.read_csv(d/'metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id');assert a.obs_names.equals(md.index)
qp=run/'stage_c_clean_views'/(name+'_cleaning_decisions_PRIVATE.pkl.gz');alt=run/'stage_c_clean_views/T6_source_aware_parent_review_PRIVATE.pkl.gz';qp=alt if name=='T6' and alt.exists() else qp;q=pd.read_pickle(qp).loc[a.obs_names];assert q.clean_lineage_model_eligible.all()
ids=a.var_names.get_indexer(genes);raw=np.zeros((a.n_obs,len(genes)),np.int32);found=ids>=0;raw[:,found]=a.X[:,ids[found]].toarray();total=np.asarray(a.X.sum(axis=1)).ravel();cp=raw*(10000/np.maximum(total,1)[:,None]);cl=md['source_snn_res.0.4'].astype(int).to_numpy();nc=int(cl.max()+1);assert set(cl)==set(range(nc))
gses=sorted(a.obs.gse_id.astype(str).unique());av=[]
for gs in gses:
 with h5py.File(PROJECT_ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad'),'r') as f:var=pd.Index(read_elem(f['var']).index.astype(str))
 av.append(pd.Index(genes).isin(var))
av=np.asarray(av);gi=pd.Index(gses).get_indexer(a.obs.gse_id.astype(str));assert not ((raw>0)&~av[gi]).any();H=sp.csr_matrix((np.ones(a.n_obs,np.int32),(cl,np.arange(a.n_obs))),shape=(nc,a.n_obs));S=sp.csr_matrix((np.ones(a.n_obs,np.int32),(cl,gi)),shape=(nc,len(gses)));den=np.asarray(S@av.astype(np.int32));nd=np.asarray(H@(raw>0).astype(np.int32));su=np.asarray(H@cp);n=np.bincount(cl,minlength=nc)
e=pd.DataFrame([dict(cluster=c,gene=g,n_cells=int(n[c]),n_measured=int(den[c,j]),coverage_fraction=float(den[c,j]/n[c]),n_expressing=int(nd[c,j]),pct_expressing_measured=float(nd[c,j]/den[c,j]) if den[c,j] else np.nan,mean_cp10k_measured=float(su[c,j]/den[c,j]) if den[c,j] else np.nan) for c in range(nc) for j,g in enumerate(genes)]);e.to_csv(o/'measured_marker_expression.tsv',sep='\t',index=False,na_rep='NA')
rows=[];refs=[]
for c in range(nc):
 t=q.iloc[np.flatnonzero(cl==c)];z={'cluster':c,'n_cells':len(t),'n_sources':t.gse_id.nunique(),'max_source_fraction':t.gse_id.astype(str).value_counts().iloc[0]/len(t),'decontx_median':float(pd.to_numeric(t.decontx_contamination,errors='coerce').median()),'old_cluster_dominant':int(t.review_cluster_before_cleanup.value_counts().index[0]),'old_cluster_dominant_fraction':float(t.review_cluster_before_cleanup.value_counts().iloc[0]/len(t))}
 for col in ['T_parent_cross_component_raw2','Treg_core_positive_raw1','Treg_core_positive_raw2','T_parent_panel_evaluable','Treg_positive_rule_evaluable','B_review_program_review','Plasma_review_program_review','Myeloid_review_program_review','joint_review_original','distinct_competing_anchor_supported','two_CD3_genes_raw2_positive','source_aware_T_positive_evidence','additional_CD3_supported_T_evidence','source_aware_Treg_core_raw1_positive']+[key for key in t.columns if key.startswith('adjudication_')]:
  if col in t:z['n_'+col]=int(t[col].sum())
 rows.append(z)
 for col in ['ev_celltypist_label','scatomic_full_label']:
  if col in t:
   for label,num in t[col].astype(object).fillna('missing').astype(str).value_counts().head(5).items():refs.append(dict(cluster=c,reference_field=col,label=label,n_cells=int(num),n_cluster=len(t)))
pd.DataFrame(rows).to_csv(o/'cluster_marker_QC_summary.tsv',sep='\t',index=False,na_rep='NA');pd.DataFrame(refs).to_csv(o/'cluster_reference_top5.tsv',sep='\t',index=False,na_rep='NA')
de=pd.read_csv(d/'RNA_FindAllMarkers_common_measured.tsv',sep='\t');de=de[(de.avg_log2FC>0)&(de.p_val_adj<.05)].sort_values(['cluster','avg_log2FC'],ascending=[True,False]).groupby('cluster').head(12);de.to_csv(o/'top12_RNA_DEG.tsv',sep='\t',index=False,na_rep='NA')
mf={'status':'completed_clean_marker_aggregation_not_final_identity','UTC':datetime.now(timezone.utc).isoformat(),'lineage':name,'n_cells':a.n_obs,'n_clusters':nc,'n_genes':len(genes),'n_marker_records':len(e),'measured_denominators_verified':True,'reference_labels_used_in_fit':False,'final_identity_pending':True};(o/'complete.json').write_text(json.dumps(mf,indent=2))
tables={f:pd.read_csv(o/f,sep='\t').to_json(orient='records') for f in ['measured_marker_expression.tsv','cluster_marker_QC_summary.tsv','cluster_reference_top5.tsv','top12_RNA_DEG.tsv']}
for f in ['cluster_resolution_counts.tsv','cluster_source_concentration.tsv','old_to_clean_cluster_crosswalk.tsv']:tables[f]=pd.read_csv(d/f,sep='\t').to_json(orient='records')
packet={'manifest':mf,'D_complete':json.loads((d/'complete.json').read_text()),'D_terminal':json.loads((run/'launch'/(name+'_refit_terminal.json')).read_text()),'tables':tables};b=zlib.compress(json.dumps(packet).encode(),9);(o/'aggregate_transfer.zlib').write_bytes(b);print(json.dumps({'manifest':mf,'packet_bytes':len(b),'packet_sha256':hashlib.sha256(b).hexdigest()}))
