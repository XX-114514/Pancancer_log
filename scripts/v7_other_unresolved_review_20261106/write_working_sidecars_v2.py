"""Write new V7 working annotations; retain frozen labels and explicit review routes."""
import json,re,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd,anndata as ad

r=PROJECT_ROOT/'runs/20261106_v7_other_unresolved_review_v1'
o=r/'stage_f_working_annotations';o.mkdir(exist_ok=False)
q=pd.read_pickle(r/'stage_c_QC_views/QC_decisions_PRIVATE.pkl.gz')
all_markers=ad.read_h5ad(r/'stage_a_candidates/candidate_marker_raw_PRIVATE.h5ad')
assert all_markers.obs_names.equals(q.index)
ptprc=pd.Series(np.asarray(all_markers.X[:,all_markers.var_names.get_loc('PTPRC')].toarray()).ravel(),index=q.index).ge(2)
key={'T':'T_positive_source_aware','NK':'NK_receptor_effector_positive','B':'B_distinct_components_positive',
     'ASC':'ASC_distinct_components_positive','Myeloid':'Myeloid_mature_components_positive',
     'Endothelial':'Endothelial_distinct_components_positive','Epithelial':'Epithelial_distinct_components_positive',
     'Fibroblast':'Fibroblast_distinct_components_positive','Erythroid':'Ery_structural_synthesis_nonHB_positive',
     'Platelet':'Platelet_secreted_structural_positive'}
excluded_competitors={
    'T':{'T','NK'},'Myeloid':{'Myeloid'},'Mural':{'Fibroblast','Endothelial'},
    'Epithelial':{'Epithelial','Fibroblast'},'Endothelial':{'Endothelial','Fibroblast'},
    'Erythroid':{'Erythroid','Platelet'},'Platelet_Megakaryocytic':{'Platelet','Endothelial'},
    'HSPC_like':{'Myeloid','Erythroid','Platelet'},'Immature_myeloid_like':{'Myeloid'},
    'Hepatocyte_like':{'Epithelial'},'Melanocytic':{'Epithelial','Fibroblast'},
}
blood_re=r'(?<![A-Za-z])(?:AML|ALL|MDS|MPN|CML|CLL|DLBCL)(?![A-Za-z])|leuk|lymphoma|myeloma|Hodgkin'
sidecars=[];summaries=[];old_cross=[];cluster_summaries=[];seen=set()
for name in ['Erythroid','Platelet','HSPC','ILC']:
    e=r/('stage_e2_'+name+'_extended_marker_review');audit=json.loads((e/'independent_postfit_audit.json').read_text());assert audit['status']=='PASS'
    d=r/('stage_d_'+name+'_refit');md=pd.read_csv(d/'metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id')
    assert not seen.intersection(md.index);seen.update(md.index)
    x=q.loc[md.index];assert x.current_lineage_model_eligible.all()
    labels=pd.read_csv(r/'review_decisions'/(name+'_cluster_decisions.tsv'),sep='\t').set_index('cluster')
    cl=md['source_snn_res.0.4'].astype(int);assert set(cl)==set(labels.index) and labels.index.is_unique
    t=pd.DataFrame(index=x.index);t.index.name='cell_id';t['model_pool']=name;t['diagnostic_cluster']=cl
    t['old_final_annotation_v7']=x.final_annotation_v7.astype(object)
    t['old_major_cell_type']=x.major_celltype_v7.astype(object)
    for col in ['working_parent','identity_review_status','state_context','naming_markers','reason']:
        t['cluster_'+col]=cl.map(labels[col])
    t['working_parent']=t.cluster_working_parent;t['assignment_basis']=np.where(t.cluster_identity_review_status.eq('broad_working'),'RNA_supported_cluster_broad_working','cluster_context_requires_review')
    t['identity_review_required']=~t.cluster_identity_review_status.eq('broad_working')
    t['positive_competing_programme_review']=False
    for parent in t.loc[t.cluster_identity_review_status.eq('broad_working'),'cluster_working_parent'].unique():
        members=t.cluster_working_parent.eq(parent)&t.cluster_identity_review_status.eq('broad_working')
        foreign=[col for k,col in key.items() if k not in excluded_competitors[parent]]
        conflict=x[foreign].any(axis=1)&members
        t.loc[conflict,'identity_review_required']=True
        t.loc[conflict,'positive_competing_programme_review']=True
        t.loc[conflict,'assignment_basis']='cluster_broad_working_with_positive_competitor_review'
    # In a reviewed parent, positive T apparatus supports routing to T; missing T is not ILC proof.
    open_route=~t.cluster_identity_review_status.eq('broad_working')
    foreign_T=x[[v for k,v in key.items() if k not in {'T','NK'}]].any(axis=1)
    direct_T=open_route&x[key['T']]&~foreign_T
    t.loc[direct_T,'working_parent']='T';t.loc[direct_T,'assignment_basis']='positive_T_apparatus_route';t.loc[direct_T,'identity_review_required']=False
    # Direct erythroid components can recover cells inside weak/mixed erythroid partitions.
    foreign_E=x[[v for k,v in key.items() if k not in {'Erythroid','Platelet'}]].any(axis=1)
    direct_E=open_route&(name=='Erythroid')&x[key['Erythroid']]&~foreign_E
    t.loc[direct_E,'working_parent']='Erythroid';t.loc[direct_E,'assignment_basis']='positive_nonHB_structural_synthesis_route';t.loc[direct_E,'identity_review_required']=False
    # ILC-like remains a candidate, even with multi-component programme evidence.
    non_ILC=x[[v for k,v in key.items() if k not in {'Erythroid','Platelet'}]].any(axis=1)
    non_ILC=non_ILC|x['Melanocytic_core_positive_raw2']|x['Pericyte_core_positive_raw2']|x['Smooth_muscle_core_positive_raw2']
    ilc_candidate=open_route&(name=='ILC')&x.ILC2_or_ILC3_specific_program_positive&ptprc.loc[x.index]&~non_ILC
    t.loc[ilc_candidate,'working_parent']='ILC_like_candidate';t.loc[ilc_candidate,'assignment_basis']='positive_ILC_programme_candidate_no_Lin_negative_proof';t.loc[ilc_candidate,'identity_review_required']=True
    t['hematologic_disease_context']=x.cancer_type.fillna('').str.contains(blood_re,flags=re.I,regex=True)
    t['normal_HSC_identity_established']=False;t['final_fine_identity_accepted']=False;t['canonical_V7_changed']=False
    t['QC_view_version']='v3';t['primary_QC_keep']=True
    for k,col in key.items():t['direct_'+k+'_positive_evidence']=x[col].astype(bool)
    t['direct_HSPC_hematopoietic_positive_evidence']=x.HSPC_hematopoietic_components_positive.astype(bool)
    t['direct_ILC2_or_ILC3_program_positive_evidence']=x.ILC2_or_ILC3_specific_program_positive.astype(bool)
    p=o/(name+'_working_annotation_PRIVATE.pkl.gz');t.to_pickle(p);back=pd.read_pickle(p);pd.testing.assert_frame_equal(t,back)
    assert back.index.equals(md.index) and back.old_final_annotation_v7.equals(x.final_annotation_v7.astype(object).rename('old_final_annotation_v7'))
    sidecars.append(t)
    summaries.append(t.groupby(['model_pool','working_parent','assignment_basis','identity_review_required','hematologic_disease_context'],dropna=False).size().rename('n_cells').reset_index())
    old_cross.append(t.groupby(['model_pool','old_final_annotation_v7','working_parent','identity_review_required'],dropna=False).size().rename('n_cells').reset_index())
    for c,z in labels.iterrows():
        take=t.diagnostic_cluster.eq(c);row=dict(lineage=name,cluster=int(c),n_cells=int(take.sum()),**z.to_dict())
        row['n_identity_review_required']=int(t.loc[take,'identity_review_required'].sum());row['n_hematologic_disease_context']=int(t.loc[take,'hematologic_disease_context'].sum());cluster_summaries.append(row)
combined=pd.concat(sidecars);assert combined.index.is_unique and len(combined)==70778
combined.to_pickle(o/'four_models_working_annotation_PRIVATE.pkl.gz')
routes=pd.DataFrame(index=q.index);routes.index.name='cell_id';routes['old_final_annotation_v7']=q.final_annotation_v7.astype(object)
routes['review_route']='background_not_identity_refit_this_pass';routes.loc[combined.index,'review_route']='four_model_working_annotation_available'
routes.loc[q.new_QC_exclude,'review_route']='QC_v3_operational_exclusion';routes.loc[q.new_QC_holdout,'review_route']='QC_v3_hematologic_context_holdout'
assert int(routes.review_route.eq('QC_v3_operational_exclusion').sum())==913 and int(routes.review_route.eq('QC_v3_hematologic_context_holdout').sum())==2
routes.to_pickle(o/'all_candidate_review_routes_PRIVATE.pkl.gz')
pd.testing.assert_frame_equal(routes,pd.read_pickle(o/'all_candidate_review_routes_PRIVATE.pkl.gz'))
pd.concat(summaries).to_csv(o/'working_identity_counts.tsv',sep='\t',index=False)
pd.concat(old_cross).to_csv(o/'old_to_working_identity_crosswalk.tsv',sep='\t',index=False)
pd.DataFrame(cluster_summaries).to_csv(o/'all_cluster_annotation_review.tsv',sep='\t',index=False)
routes.review_route.value_counts().rename_axis('review_route').rename('n_cells').to_csv(o/'candidate_route_counts.tsv',sep='\t')
manifest={'status':'working_sidecars_written_and_readback_verified','UTC':datetime.now(timezone.utc).isoformat(),'n_review_scope':len(q),'n_model_cells':len(combined),'n_diagnostic_clusters':len(cluster_summaries),'n_identity_review_required':int(combined.identity_review_required.sum()),'n_broad_working_without_positive_competitor_review':int((~combined.identity_review_required).sum()),'canonical_V7_changed':False,'normal_HSC_identity_established':False,'final_fine_identity_accepted':False,'annotation_masks_change_global_QC':False,'files':[{'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in o.iterdir()]}
(o/'complete.json').write_text(json.dumps(manifest,indent=2));print(json.dumps({k:v for k,v in manifest.items() if k!='files'}))
