import json,hashlib,re
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad

r=PROJECT_ROOT/'runs/20261106_v7_B_unresolved_review_v1';o=r/'stage_f_B_working_annotations'
assert json.loads((r/'stage_e_B_marker_review/independent_postfit_audit.json').read_text())['status']=='PASS'
assert json.loads((r/'stage_e2_B_QC_context/independent_additional_marker_audit.json').read_text())['status']=='PASS'
q=pd.read_pickle(r/'stage_c_QC_views/QC_decisions_PRIVATE.pkl.gz');md=pd.read_csv(r/'stage_d_B_unresolved_refit/metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id')
lab=pd.read_csv(r/'review_decisions/B_cluster_decisions.tsv',sep='\t').set_index('cluster');assert lab.index.is_unique and set(lab.index)==set(md['source_snn_res.0.4'].astype(int))
target=q[q.is_annotation_target].copy();assert len(target)==158257
t=pd.DataFrame(index=target.index);t.index.name='cell_id';t['old_final_annotation_v7']=target.final_annotation_v7.astype(object);t['old_major_celltype_v7']=target.major_celltype_v7.astype(object);t['is_annotation_target']=True
t['primary_QC_keep']=target.current_QC_keep.to_numpy();t['in_current_B_model']=target.current_lineage_model_eligible.to_numpy();t['diagnostic_cluster']=pd.Series(pd.NA,index=t.index,dtype='Int64')
ids=t.index[t.in_current_B_model];assert set(ids)==set(md.index[q.loc[md.index,'is_annotation_target']]);t.loc[ids,'diagnostic_cluster']=md.loc[ids,'source_snn_res.0.4'].astype(int).to_numpy()
t['working_parent']=t.old_major_celltype_v7;t['working_annotation_label']=t.old_final_annotation_v7;t['working_subtype_candidate']='not_assigned';t['state_context']='not_refit_due_to_QC';t['naming_markers']='not_reassigned';t['identity_review_required']=True
t['assignment_basis']='QC_excluded_before_B_identity_refit';t['reason']='Excluded from this clean refit by the registered operational mixed-RNA QC rule; old identity retained for provenance.'
cl=t.loc[ids,'diagnostic_cluster'].astype(int)
for col in ['working_parent','working_annotation_label','working_subtype_candidate','state_context','naming_markers','reason']:
 t.loc[ids,col]=cl.map(lab[col]).to_numpy()
t.loc[ids,'identity_review_required']=~cl.map(lab.identity_review_status).eq('broad_working').to_numpy()
t.loc[ids,'assignment_basis']=np.where(t.loc[ids,'identity_review_required'],'cluster_context_requires_parent_review','manual_cluster_multi_marker_broad_working')
key={'T':'T_positive_source_aware','NK':'NK_receptor_effector_positive','B':'B_distinct_components_positive','ASC':'ASC_distinct_components_positive','Myeloid':'Myeloid_mature_components_positive','Endothelial':'Endothelial_distinct_components_positive','Epithelial':'Epithelial_distinct_components_positive','Fibroblast':'Fibroblast_distinct_components_positive','Erythroid':'Ery_structural_synthesis_nonHB_positive','Platelet':'Platelet_secreted_structural_positive','Mast':'Mast_core_positive_raw2'}
expected={'Bcell':{'B','ASC'},'B_lineage':{'B','ASC'},'Plasma':{'B','ASC'},'Tcell':{'T','NK'},'NK':{'NK','T'},'Myeloid':{'Myeloid'},'Macrophage':{'Myeloid'},'Monocyte':{'Myeloid'},'Dendritic':{'Myeloid','ASC'},'Granulocyte':{'Myeloid'},'Endothelial':{'Endothelial','Fibroblast'},'Mural':{'Endothelial','Fibroblast'},'Fibroblast':{'Fibroblast','Endothelial','Epithelial'},'Epithelial':{'Epithelial','Fibroblast'},'Melanocytic':{'Epithelial','Fibroblast'},'Mast':{'Mast'},'Hepatocyte_like':{'Epithelial'},'Platelet_Megakaryocytic':{'Platelet','Endothelial'},'HSPC_like':{'Myeloid','Erythroid','Platelet'}}
t['positive_competing_programme_review']=False;t['related_B_ASC_programme_present']=target.B_distinct_components_positive&target.ASC_distinct_components_positive
for par,allowed in expected.items():
 members=t.in_current_B_model&t.working_parent.eq(par)&~t.identity_review_required
 foreign=target[[col for k,col in key.items() if k not in allowed]].any(axis=1)&members
 t.loc[foreign,'identity_review_required']=True;t.loc[foreign,'positive_competing_programme_review']=True;t.loc[foreign,'assignment_basis']='cluster_broad_working_with_positive_competitor_review'
# Positive single-cell routes only in open parent contexts, never use a missing marker as a negative identity.
open_route=t.in_current_B_model&cl.reindex(t.index).map(lab.identity_review_status).ne('broad_working')
foreign_B=target[[v for k,v in key.items() if k not in {'B','ASC'}]].any(axis=1)
direct_B=open_route&target.B_distinct_components_positive&~foreign_B
t.loc[direct_B,'working_parent']='Bcell';t.loc[direct_B,'working_annotation_label']='Bcell_broad';t.loc[direct_B,'working_subtype_candidate']='not_assigned';t.loc[direct_B,'identity_review_required']=False;t.loc[direct_B,'assignment_basis']='positive_B_distinct_component_route'
foreign_T=target[[v for k,v in key.items() if k not in {'T','NK'}]].any(axis=1)
direct_T=open_route&target.T_positive_source_aware&~foreign_T
t.loc[direct_T,'working_parent']='Tcell';t.loc[direct_T,'working_annotation_label']='Tcell_broad';t.loc[direct_T,'working_subtype_candidate']='not_assigned';t.loc[direct_T,'identity_review_required']=False;t.loc[direct_T,'assignment_basis']='positive_T_apparatus_route'
for k,col in key.items():t['direct_'+k+'_positive_evidence']=target[col].astype(bool)
t['hematologic_disease_context']=target.cancer_type.fillna('').str.contains(r'(?<![A-Za-z])(?:AML|ALL|MDS|MPN|CML|CLL|DLBCL)(?![A-Za-z])|leuk|lymphoma|myeloma|Hodgkin',flags=re.I,regex=True)
t['fine_identity_final_accepted']=False;t['canonical_V7_changed']=False;t['state_semantics']='cluster-associated RNA context, not a per-cell functional state confirmation'
assert t.loc[~t.in_current_B_model,'identity_review_required'].all() and not t.loc[~t.in_current_B_model,'primary_QC_keep'].any()
o.mkdir(exist_ok=False);p=o/'B_targets_working_annotation_PRIVATE.pkl.gz';t.to_pickle(p);pd.testing.assert_frame_equal(t,pd.read_pickle(p))
t.groupby(['old_final_annotation_v7','working_parent','working_annotation_label','working_subtype_candidate','identity_review_required','primary_QC_keep','assignment_basis'],observed=True,dropna=False).size().rename('n_cells').reset_index().to_csv(o/'B_old_to_working_annotation_counts.tsv',sep='\t',index=False)
t.groupby(['working_parent','working_annotation_label','working_subtype_candidate','identity_review_required','primary_QC_keep'],observed=True,dropna=False).size().rename('n_cells').reset_index().to_csv(o/'B_working_identity_counts.tsv',sep='\t',index=False)
rows=[]
for c,z in lab.iterrows():
 g=t[t.in_current_B_model&t.diagnostic_cluster.eq(c).fillna(False)];rows.append({'cluster':int(c),'n_target_cells':len(g),'n_identity_review_required':int(g.identity_review_required.sum()),'n_broad_working_targets':int((~g.identity_review_required).sum()),**z.to_dict()})
pd.DataFrame(rows).to_csv(o/'B_all_cluster_working_review.tsv',sep='\t',index=False)
out={'status':'B_target_working_sidecar_written_and_readback_verified','UTC':datetime.now(timezone.utc).isoformat(),'n_target_records':len(t),'n_targets_in_clean_model':int(t.in_current_B_model.sum()),'n_QC_excluded_before_refit':int((~t.in_current_B_model).sum()),'n_broad_working_targets':int((t.in_current_B_model&~t.identity_review_required).sum()),'n_identity_review_required_in_model':int((t.in_current_B_model&t.identity_review_required).sum()),'n_direct_B_routes':int(direct_B.sum()),'n_direct_T_routes':int(direct_T.sum()),'n_related_B_ASC_programme_cells':int(t.related_B_ASC_programme_present.sum()),'n_controls_overlaid':0,'canonical_V7_changed':False,'fine_identity_final_accepted':False,'sidecar_sha256':hashlib.sha256(p.read_bytes()).hexdigest()};(o/'complete.json').write_text(json.dumps(out,indent=2));print(json.dumps(out))
