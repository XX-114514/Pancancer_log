import os
from pathlib import Path
PROJECT_ROOT=Path(os.environ["PROJECT_ROOT"])
import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
run=PROJECT_ROOT/'runs/20261106_v7_QC_clean_reannotation_v1';o=run/'stage_h_working_review_overlays_v1';o.mkdir(exist_ok=False);summaries=[]
for name in ['Mast','T6']:
 d=run/('stage_i_Mast_NK_clean_refit_v2' if name=='Mast' else 'stage_d_T6_clean_refit');e=run/('stage_j_Mast_NK_clean_marker_review_v3' if name=='Mast' else 'stage_l_T6_extended_Ig_marker_review_v1');assert json.loads((e/'independent_postfit_audit.json').read_text())['status']=='PASS'
 qpath=run/'stage_g_Mast_NK_cleanup_v2/Mast_cleaning_decisions_PRIVATE.pkl.gz'
 if name=='T6':qpath=run/'stage_c_clean_views/T6_source_aware_parent_review_PRIVATE.pkl.gz'
 q=pd.read_pickle(qpath);md=pd.read_csv(d/'metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id');labels=pd.read_csv(run/'config'/(name+'_clean_working_labels.tsv'),sep='\t');assert labels.cluster.is_unique and set(labels.cluster)==set(md['source_snn_res.0.4'].unique()) and md.index.equals(q.index[q.clean_lineage_model_eligible])
 cl=pd.Series(pd.NA,index=q.index,dtype='Int64');cl.loc[md.index]=md['source_snn_res.0.4'].astype(int).to_numpy();q['clean_review_cluster']=cl;q['working_parent_review']='not_in_clean_lineage_fit'
 if name=='Mast':q.loc[md.index,'working_parent_review']='Mast_positive_RNA_support'
 else:
  pos=q.source_aware_T_positive_evidence&q.clean_lineage_model_eligible;q.loc[md.index,'working_parent_review']='T_candidate_not_positive_by_current_rules';q.loc[pos,'working_parent_review']='T_positive_RNA_support'
 q['working_parent_identity']=cl.map(labels.set_index('cluster').working_parent_identity).astype(object).fillna('outside_current_clean_model')
 q['broad_assignment_basis']='outside_current_clean_model';q.loc[md.index,'broad_assignment_basis']='manual_cluster_multi_marker_review_with_source_scope; direct_single_cell_positive_rules_retained_separately_not_required_for_every_cell'
 q['working_fine_subtype']='not_assigned';q['working_cluster_context']=cl.map(labels.set_index('cluster').display_annotation).astype(object).fillna('not_in_clean_lineage_fit');q['cluster_context_requires_review']=q.working_cluster_context.str.contains('[review]',regex=False);q['state_label_semantics']='cluster_marker_association_not_a_cellwise_state_call';q['review_overlay_version']='V7_literature_QC_working_20261106_v1';q['final_fine_identity_accepted']=False
 p=o/(name+'_working_annotation_QC_PRIVATE.pkl.gz');q.to_pickle(p);back=pd.read_pickle(p);pd.testing.assert_frame_equal(q,back)
 counts=q.groupby(['primary_QC_action','clean_lineage_model_eligible','working_parent_identity','working_parent_review','working_fine_subtype','working_cluster_context'],dropna=False,observed=True).size().reset_index(name='n_cells');assert counts.n_cells.sum()==len(q);counts.to_csv(o/(name+'_working_annotation_counts.tsv'),sep='\t',index=False,na_rep='NA')
 summaries.append({'lineage':name,'n_full_review_pool':len(q),'n_model':len(md),'n_primary_QC_keep':int(q.primary_QC_keep.sum()),'n_positive_parent_in_model':int(q.working_parent_review.isin(['Mast_positive_RNA_support','T_positive_RNA_support']).sum()),'n_working_broad_T':int(q.working_parent_identity.eq('T').sum()),'n_working_broad_Mast':int(q.working_parent_identity.eq('Mast').sum()),'n_T_with_competing_RNA_cluster_review':int(q.working_parent_identity.eq('T_with_competing_RNA_review').sum()),'n_T_candidate_low_information':int(q.working_parent_identity.eq('T_candidate_low_information').sum()),'n_diagnostic_clusters':len(labels),'n_cells_in_clusters_requiring_review':int(q.cluster_context_requires_review.sum()),'canonical_V7_changed':False,'fine_final_labels_assigned':0,'private_overlay_SHA256':hashlib.sha256(p.read_bytes()).hexdigest(),'readback_exact':True})
out={'status':'working_annotation_and_QC_overlays_created_readback_verified','UTC':datetime.now(timezone.utc).isoformat(),'lineages':summaries,'limits':'Working broad evidence and cluster marker contexts only. Review clusters remain explicit; no final subtype acceptance or automatic downstream CoVarNet input replacement.'};(o/'complete.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
