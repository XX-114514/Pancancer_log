import os,json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))));r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';o=r/'stage_m_T6_working_overlay_v1';o.mkdir(exist_ok=False)
obs=pd.read_pickle(r/'stage_l_T6_marker_review_v1/T6_marker_review_PRIVATE.pkl.gz');md=pd.read_csv(r/'stage_k_T6_source_robust_diagnostic_v1/fit/metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id');assert md.index.equals(obs.index)
labels=json.loads((r/'config/T6_working_cluster_labels_v1.json').read_text());cl=obs['T6_review_cluster_source_res0.4'].astype(int);assert set(map(str,cl.unique()))==set(labels)
obs['cluster_annotation_working']=cl.astype(str).map(labels);obs['parent_identity']=np.where(obs.T_parent_cross_component_raw2,'T_marker_supported_candidate','T_candidate_low_signal');obs['identity_subtype']=pd.NA;obs['subtype_status']='major_or_marker_defined_candidate_not_finally_assigned';obs['subtype_candidate']=np.where(obs.Treg_core_positive_raw1,'Treg_positive_RNA_candidate','not_assigned_by_this_positive_rule')
obs['qc_review']=np.select([obs.joint_competitor_doublet_candidate,obs.parent_conflicts.ne(''),~obs.T_parent_cross_component_raw2],['joint_doublet_candidate_review','competing_program_review','low_marker_support_review'],default='no_extended_competitor_detected_under_this_rule');obs['qc_pending']=obs.qc_review.ne('no_extended_competitor_detected_under_this_rule');obs['version']='V7_T6_20261105_working_review_v1'
obs.to_pickle(o/'T6_hierarchical_working_review_PRIVATE.pkl.gz');back=pd.read_pickle(o/'T6_hierarchical_working_review_PRIVATE.pkl.gz');pd.testing.assert_frame_equal(obs,back)
qc=[]
for c,t in obs.groupby('T6_review_cluster_source_res0.4',observed=True):
 z={'cluster':int(c),'working_annotation':labels[str(c)],'n_cells':len(t),'n_sources':t.gse_id.nunique()}
 for name in ['total_counts','n_genes_by_counts','pct_counts_mt','decontx_contamination']:
  if name in t:
   v=pd.to_numeric(t[name],errors='coerce');z[name+'_median']=float(v.median());z[name+'_q10']=float(v.quantile(.1));z[name+'_q90']=float(v.quantile(.9))
 qc.append(z)
pd.DataFrame(qc).to_csv(o/'T6_cluster_QC_depth.tsv',sep='\t',index=False,na_rep='NA')
rows=[]
for col in ['source_snn_res.0.2','source_snn_res.0.6','source_snn_res.0.8']:
 other=md[col].astype(int);joint=pd.crosstab(cl,other);sizes=other.value_counts()
 for c in joint.index:
  best=int(joint.loc[c].idxmax());inter=int(joint.loc[c,best]);own=int((cl==c).sum());other_n=int(sizes[best]);rows.append({'cluster_res0.4':int(c),'compared_resolution':col,'largest_overlap_cluster':best,'n_cells_res0.4':own,'n_cells_other_cluster':other_n,'intersection':inter,'fraction_current_in_other':inter/own,'fraction_other_from_current':inter/other_n,'jaccard':inter/(own+other_n-inter)})
pd.DataFrame(rows).to_csv(o/'T6_resolution_overlap_summary.tsv',sep='\t',index=False)
route=obs.groupby(['cluster_annotation_working','qc_review'],observed=True).size().reset_index(name='n_cells');assert route.n_cells.sum()==156631;route.to_csv(o/'T6_working_identity_QC_routes.tsv',sep='\t',index=False)
result={'status':'completed_T6_hierarchical_working_overlay','UTC':datetime.now(timezone.utc).isoformat(),'n_cells':len(obs),'n_clusters':len(labels),'old_fields_preserved':True,'n_final_subtypes_assigned':0,'n_Treg_positive_RNA_candidates':int(obs.Treg_core_positive_raw1.sum()),'n_Treg_positive_RNA_sources':obs.loc[obs.Treg_core_positive_raw1,'gse_id'].nunique(),'limits':'Cluster labels are marker-guided review hypotheses, not labels transferred to every cell. Tfh/Tph, resting/effector Treg and CD8 exhaustion are not forced. Multiple-resolution overlap is diagnostic, not a held-out replication.'}
(o/'complete.json').write_text(json.dumps(result,indent=2));print(json.dumps({'manifest':result,'tables':{name:pd.read_csv(o/name,sep='\t').to_json(orient='records') for name in ['T6_cluster_QC_depth.tsv','T6_resolution_overlap_summary.tsv','T6_working_identity_QC_routes.tsv']}}))
