"""Build a new full-axis working annotation view; never overwrite frozen V7."""
import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd

r=PROJECT_ROOT;run=r/'runs/20261106_v7_B_unresolved_review_v1';o=run/'stage_g_consolidated_annotation'
vpath=r/'runs/20260818_pancancer_5M_merge_prep/11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1/checkpoints/final_annotation_v7.pkl.gz'
prior=r/'runs/20261106_v7_QC_clean_reannotation_v1';rare=r/'runs/20261106_v7_other_unresolved_review_v1'
bpath=run/'stage_f_B_working_annotations/B_targets_working_annotation_PRIVATE.pkl.gz'
assert bpath.exists() and json.loads((run/'stage_f_B_working_annotations/independent_B_sidecar_audit.json').read_text())['status']=='PASS'
assert json.loads((rare/'stage_f_working_annotations/independent_sidecar_audit.json').read_text())['status']=='PASS'
v=pd.read_pickle(vpath);assert len(v)==4676787 and v.index.is_unique
m=np.load(run/'stage_c_QC_views/global_QC_view_v4_masks_PRIVATE.npz');assert np.array_equal(m['primary_keep'],~m['primary_exclude']&~m['QC_holdout'])
ih=hashlib.sha256(pd.util.hash_pandas_object(v.index,index=False).to_numpy().tobytes()).hexdigest();assert ih=='3f802eb8cc52f04beb27fa2987d2ffe3ddf75f4b4f51dec8e6448b26df8d1f16'
out=pd.DataFrame(index=v.index);out.index.name='cell_id'
out['frozen_parent']=v.major_celltype_v7.astype(object);out['frozen_label']=v.final_annotation_v7.astype(object)
out['frozen_confidence']=v.annotation_confidence_v7.astype(object).fillna('missing')
out['frozen_strict_major_eligible']=v.strict_major_eligible_v7.fillna(False).astype(bool)
out['frozen_strict_subtype_eligible']=v.strict_subtype_eligible_v7.fillna(False).astype(bool)
malignancy_columns=['malignancy_call_v7','malignancy_tier_v7','is_malignant_primary_v7','is_malignant_high_confidence_v7','is_malignant_sensitive_v7','strict_malignant_downstream_eligible_v7','malignancy_evidence_version_v7']
for col in malignancy_columns:
 assert col in v;out['frozen_'+col]=v[col].copy()
out['malignancy_axis_review_status']='frozen_inherited_not_reassessed'
out['current_parent']=out.frozen_parent;out['current_annotation_label']=out.frozen_label
out['annotation_layer']='frozen_V7_inherited';out['identity_review_status']='not_reaudited_in_current_cycle'
out['candidate_subtype']='not_reassessed';out['cluster_RNA_context']='frozen_label_retained';out['naming_markers']='not_reaudited'
out['evidence_scope']='not_screened_in_current_QC_cycle';out['assignment_basis']='frozen_V7_inherited_not_new_acceptance'
out['diagnostic_model']='not_in_current_refit';out['diagnostic_cluster']=pd.Series(pd.NA,index=v.index,dtype='Int64')
out['refit_identity_reviewed_this_cycle']=False;out['working_identity_review_required']=pd.Series(pd.NA,index=v.index,dtype='boolean')
out['current_QC_action']=np.select([m['primary_exclude'],m['QC_holdout']],['exclude_primary_QC','QC_holdout'],default='primary_QC_keep')
out['current_primary_QC_keep']=m['primary_keep'];out['QC_reason_code']=m['reason_code'];out['fine_identity_final_accepted_this_cycle']=False
canon={'T':'Tcell','B':'Bcell','ASC':'Plasma'}
def parent(x):return canon.get(str(x),str(x))
seen=set();layers=[]
def claim(index,layer):
 assert pd.Index(index).is_unique and pd.Index(index).isin(v.index).all()
 collision=seen.intersection(index);assert not collision,(layer,len(collision))
 seen.update(index);out.loc[index,'annotation_layer']=layer;out.loc[index,'refit_identity_reviewed_this_cycle']=True
 assert out.loc[index,'current_primary_QC_keep'].all()

for name in ['Mast','T6']:
 q=pd.read_pickle(prior/'stage_h_working_review_overlays_v1'/(name+'_working_annotation_QC_PRIVATE.pkl.gz'))
 assert q.index.is_unique and q.index.isin(v.index).all();out.loc[q.index,'evidence_scope']='Mast_T6_QC_marker_review'
 if 'final_annotation_v7' in q:assert np.array_equal(q.final_annotation_v7.astype(str).to_numpy(),v.loc[q.index,'final_annotation_v7'].astype(str).to_numpy())
 ids=q.index[q.clean_lineage_model_eligible];x=q.loc[ids];claim(ids,'Mast_T6_working_20261106')
 outside=q.index[~q.clean_lineage_model_eligible];out.loc[outside,'identity_review_status']='reviewed_pool_not_in_clean_identity_model';out.loc[outside,'working_identity_review_required']=True
 out.loc[ids,'current_parent']=x.working_parent_identity.map(parent).to_numpy();out.loc[ids,'current_annotation_label']=x.working_parent_identity.map(parent).to_numpy()
 bad=x.working_parent_identity.isin(['T_with_competing_RNA_review','T_candidate_low_information'])
 out.loc[ids,'identity_review_status']=np.where(bad,'identity_review_required','broad_working')
 out.loc[ids,'working_identity_review_required']=bad.to_numpy();out.loc[ids,'candidate_subtype']='not_assigned'
 out.loc[ids,'cluster_RNA_context']=x.working_cluster_context.astype(str).to_numpy();out.loc[ids,'assignment_basis']=x.broad_assignment_basis.astype(str).to_numpy()
 out.loc[ids,'diagnostic_model']=name;out.loc[ids,'diagnostic_cluster']=x.clean_review_cluster.astype('Int64').to_numpy()
 out.loc[ids,'naming_markers']='TPSAB1;TPSB2;CPA3;KIT' if name=='Mast' else 'CD3D;CD3E;TRAC;TRBC2'
 layers.append({'layer':name,'n_refit_target_cells':len(ids),'n_identity_review_required':int(bad.sum()),'n_context_controls_overlaid':0})

qscope=pd.read_pickle(rare/'stage_c_QC_views/QC_decisions_PRIVATE.pkl.gz');out.loc[qscope.index,'evidence_scope']='rare_four_marker_QC_screen'
q=pd.read_pickle(rare/'stage_f_working_annotations/four_models_working_annotation_PRIVATE.pkl.gz');claim(q.index,'rare_four_working_20261106')
assert np.array_equal(q.old_final_annotation_v7.astype(str).to_numpy(),v.loc[q.index,'final_annotation_v7'].astype(str).to_numpy())
out.loc[q.index,'current_parent']=q.working_parent.map(parent).to_numpy();out.loc[q.index,'current_annotation_label']=q.working_parent.map(parent).to_numpy()
out.loc[q.index,'identity_review_status']=np.where(q.identity_review_required,'identity_review_required','broad_working');out.loc[q.index,'working_identity_review_required']=q.identity_review_required.to_numpy()
out.loc[q.index,'candidate_subtype']='not_assigned';out.loc[q.index,'cluster_RNA_context']=q.cluster_state_context.astype(str).to_numpy();out.loc[q.index,'naming_markers']=q.cluster_naming_markers.astype(str).to_numpy()
out.loc[q.index,'assignment_basis']=q.assignment_basis.astype(str).to_numpy();out.loc[q.index,'diagnostic_model']=q.model_pool.astype(str).to_numpy();out.loc[q.index,'diagnostic_cluster']=q.diagnostic_cluster.to_numpy()
layers.append({'layer':'rare_four','n_refit_target_cells':len(q),'n_identity_review_required':int(q.identity_review_required.sum()),'n_context_controls_overlaid':0})

bqscope=pd.read_pickle(run/'stage_c_QC_views/QC_decisions_PRIVATE.pkl.gz');out.loc[bqscope.index,'evidence_scope']='B_ASC_marker_QC_screen'
q=pd.read_pickle(bpath);assert q.index.is_unique and q.index.isin(v.index).all() and q.is_annotation_target.all()
assert np.array_equal(q.old_final_annotation_v7.astype(str).to_numpy(),v.loc[q.index,'final_annotation_v7'].astype(str).to_numpy())
ids=q.index[q.in_current_B_model];x=q.loc[ids];claim(ids,'B_working_20261106')
outside=q.index[~q.in_current_B_model];out.loc[outside,'identity_review_status']='QC_excluded_before_B_identity_refit';out.loc[outside,'working_identity_review_required']=True
out.loc[ids,'current_parent']=x.working_parent.map(parent).to_numpy();out.loc[ids,'current_annotation_label']=x.working_annotation_label.astype(str).to_numpy()
out.loc[ids,'identity_review_status']=np.where(x.identity_review_required,'identity_review_required','broad_working');out.loc[ids,'working_identity_review_required']=x.identity_review_required.to_numpy()
out.loc[ids,'candidate_subtype']=x.working_subtype_candidate.astype(str).to_numpy();out.loc[ids,'cluster_RNA_context']=x.state_context.astype(str).to_numpy();out.loc[ids,'naming_markers']=x.naming_markers.astype(str).to_numpy()
out.loc[ids,'assignment_basis']=x.assignment_basis.astype(str).to_numpy();out.loc[ids,'diagnostic_model']='B_unresolved_targets';out.loc[ids,'diagnostic_cluster']=x.diagnostic_cluster.astype('Int64').to_numpy()
controls=bqscope.index[bqscope.current_lineage_model_eligible&~bqscope.is_annotation_target];assert not set(controls).intersection(ids)
assert out.loc[controls,'annotation_layer'].ne('B_working_20261106').all()
layers.append({'layer':'B_unresolved_targets','n_refit_target_cells':len(ids),'n_identity_review_required':int(x.identity_review_required.sum()),'n_context_controls_overlaid':0})

changed_parent=out.refit_identity_reviewed_this_cycle&out.current_parent.ne(out.frozen_parent);out.loc[changed_parent,'malignancy_axis_review_status']='parent_changed_recheck_context_not_a_new_malignancy_call'
assert len(seen)==int(out.refit_identity_reviewed_this_cycle.sum())==sum(x['n_refit_target_cells'] for x in layers)
assert out.frozen_label.equals(v.final_annotation_v7.astype(object).rename('frozen_label'))
o.mkdir(exist_ok=False)
# One flat record per frozen cell; QC, parent, subtype candidate and cluster context remain distinct.
for col in out.select_dtypes(include='object').columns:out[col]=out[col].astype('category')
p=o/'current_annotation_all_cells_PRIVATE.pkl.gz';out.to_pickle(p);back=pd.read_pickle(p);pd.testing.assert_frame_equal(out,back)
out.to_csv(o/'current_annotation_all_cells_PRIVATE.tsv.gz',sep='\t',index=True,compression={'method':'gzip','compresslevel':1},na_rep='NA')
pd.DataFrame(layers).to_csv(o/'review_layer_counts.tsv',sep='\t',index=False)

temp=out.assign(n_cells=1)
cols=['current_parent','current_annotation_label','annotation_layer','identity_review_status','candidate_subtype','cluster_RNA_context','current_QC_action']
current=temp.groupby(cols,observed=True,dropna=False).n_cells.sum().reset_index();assert int(current.n_cells.sum())==len(v)
current.to_csv(o/'current_annotation_label_counts.tsv',sep='\t',index=False)
cross=temp.groupby(['frozen_parent','frozen_label','current_parent','current_annotation_label','identity_review_status','current_QC_action'],observed=True,dropna=False).n_cells.sum().reset_index();assert int(cross.n_cells.sum())==len(v)
cross.to_csv(o/'frozen_to_current_crosswalk.tsv',sep='\t',index=False)

rows=[]
for par,inds in out.groupby('current_parent',observed=True).indices.items():
 t=out.iloc[inds];keep=t.current_primary_QC_keep;review=t.refit_identity_reviewed_this_cycle
 rows.append({'current_parent':str(par),'n_all_records':len(t),'n_primary_QC_keep':int(keep.sum()),'n_QC_excluded':int(t.current_QC_action.eq('exclude_primary_QC').sum()),'n_QC_holdout':int(t.current_QC_action.eq('QC_holdout').sum()),'n_refit_reviewed_kept':int((keep&review).sum()),'n_broad_working_kept':int((keep&review&~t.working_identity_review_required).sum()),'n_refit_identity_review_required':int((keep&review&t.working_identity_review_required).sum()),'n_inherited_not_refit_kept':int((keep&~review).sum()),'n_current_annotation_labels':int(t.current_annotation_label.nunique())})
pd.DataFrame(rows).to_csv(o/'current_parent_summary.tsv',sep='\t',index=False)
rows=[]
for (par,lab),inds in v.groupby(['major_celltype_v7','final_annotation_v7'],observed=True).indices.items():
 t=v.iloc[inds];x=out.iloc[inds];keep=x.current_primary_QC_keep
 rows.append({'frozen_parent':str(par),'frozen_label':str(lab),'n_frozen':len(t),'n_primary_QC_keep':int(keep.sum()),'n_QC_excluded':int(x.current_QC_action.eq('exclude_primary_QC').sum()),'n_QC_holdout':int(x.current_QC_action.eq('QC_holdout').sum()),'n_refit_reviewed_this_cycle':int(x.refit_identity_reviewed_this_cycle.sum()),'n_not_refit_kept':int((keep&~x.refit_identity_reviewed_this_cycle).sum()),'n_frozen_strict_major_eligible':int(x.frozen_strict_major_eligible.sum()),'n_frozen_strict_subtype_eligible':int(x.frozen_strict_subtype_eligible.sum()),'n_sources':int(t.gse_id.nunique()),'frozen_confidence_counts':json.dumps(t.annotation_confidence_v7.astype(object).fillna('missing').astype(str).value_counts().to_dict(),ensure_ascii=False),'frozen_major_reason_counts':json.dumps(t.major_resolution_reason_v2.astype(str).value_counts().to_dict(),ensure_ascii=False),'frozen_subtype_method_counts':json.dumps(t.subtype_resolution_method_v2.astype(str).value_counts().to_dict(),ensure_ascii=False),'interpretation':'Frozen strict flags are inherited provenance, not new publication acceptance.'})
frozen=pd.DataFrame(rows);assert len(frozen)==106 and frozen.n_frozen.sum()==len(v);frozen.to_csv(o/'all_106_frozen_labels_current_status.tsv',sep='\t',index=False)
mal_cols=['frozen_parent','current_parent','frozen_malignancy_call_v7','frozen_malignancy_tier_v7','malignancy_axis_review_status','current_QC_action'];out.groupby(mal_cols,observed=True,dropna=False).size().rename('n_cells').reset_index().to_csv(o/'frozen_malignancy_axis_current_context_counts.tsv',sep='\t',index=False,na_rep='missing')
out.groupby(['frozen_label','frozen_confidence'],observed=True,dropna=False).size().rename('n_cells').reset_index().to_csv(o/'frozen_confidence_counts.tsv',sep='\t',index=False)
cluster_rows=[]
for name in ['Mast','T6']:
 z=pd.read_pickle(prior/'stage_h_working_review_overlays_v1'/(name+'_working_annotation_QC_PRIVATE.pkl.gz'));z=z[z.clean_lineage_model_eligible]
 for c,g in z.groupby('clean_review_cluster',observed=True):
  bad=g.working_parent_identity.isin(['T_with_competing_RNA_review','T_candidate_low_information'])
  cluster_rows.append({'source_model':name,'cluster_id':int(c),'cluster_working_parent':parent(g.working_parent_identity.iloc[0]),'naming_markers':'TPSAB1;TPSB2;CPA3;KIT' if name=='Mast' else 'CD3D;CD3E;TRAC;TRBC2','candidate_subtype':'not_assigned','cluster_RNA_context':str(g.working_cluster_context.iloc[0]),'n_model_cells':len(g),'n_target_cells':len(g),'n_context_controls':0,'n_identity_review_required':int(bad.sum()),'n_sources':int(v.loc[g.index,'gse_id'].nunique()),'reason':'Earlier independent marker review; see the preserved cluster RNA context and original run.','source_run':'20261106_v7_QC_clean_reannotation_v1'})
z=pd.read_pickle(rare/'stage_f_working_annotations/four_models_working_annotation_PRIVATE.pkl.gz')
for (name,c),g in z.groupby(['model_pool','diagnostic_cluster'],observed=True):
 cluster_rows.append({'source_model':name,'cluster_id':int(c),'cluster_working_parent':parent(g.cluster_working_parent.iloc[0]),'naming_markers':str(g.cluster_naming_markers.iloc[0]),'candidate_subtype':'not_assigned','cluster_RNA_context':str(g.cluster_state_context.iloc[0]),'n_model_cells':len(g),'n_target_cells':len(g),'n_context_controls':0,'n_identity_review_required':int(g.identity_review_required.sum()),'n_sources':int(v.loc[g.index,'gse_id'].nunique()),'reason':str(g.cluster_reason.iloc[0]),'source_run':'20261106_v7_other_unresolved_review_v1'})
z=pd.read_pickle(bpath);lab=pd.read_csv(run/'review_decisions/B_cluster_decisions.tsv',sep='\t').set_index('cluster');bm=pd.read_csv(run/'stage_e_B_marker_review/cluster_marker_QC_summary.tsv',sep='\t')
for row in bm.to_dict('records'):
 c=int(row['cluster']);g=z[z.in_current_B_model&z.diagnostic_cluster.eq(c).fillna(False)];d=lab.loc[c]
 cluster_rows.append({'source_model':'B_unresolved','cluster_id':c,'cluster_working_parent':parent(d.working_parent),'naming_markers':str(d.naming_markers),'candidate_subtype':str(d.working_subtype_candidate),'cluster_RNA_context':str(d.state_context),'n_model_cells':int(row['n_model_cells']),'n_target_cells':int(row['n_target_cells']),'n_context_controls':int(row['n_context_controls']),'n_identity_review_required':int(g.identity_review_required.sum()),'n_sources':int(row['n_target_sources']),'reason':str(d.reason),'source_run':'20261106_v7_B_unresolved_review_v1'})
ci=pd.DataFrame(cluster_rows);assert ci.n_target_cells.sum()==len(seen) and ci.n_context_controls.sum()==len(controls);ci.to_csv(o/'all_reviewed_cluster_inventory.tsv',sep='\t',index=False)
files=[{'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in o.iterdir() if p.is_file()]
receipt={'status':'full_axis_working_annotation_written_and_readback_verified','UTC':datetime.now(timezone.utc).isoformat(),'n_frozen_cells':len(v),'n_frozen_labels':106,'n_current_parent_categories':out.current_parent.nunique(),'n_current_annotation_labels':out.current_annotation_label.nunique(),'n_current_label_status_records':len(current),'n_primary_QC_keep':int(m['primary_keep'].sum()),'n_QC_exclude':int(m['primary_exclude'].sum()),'n_QC_holdout':int(m['QC_holdout'].sum()),'n_refit_reviewed_targets':len(seen),'n_broad_working_refit_targets':int((out.refit_identity_reviewed_this_cycle&~out.working_identity_review_required).sum()),'n_refit_identity_review_required':int((out.refit_identity_reviewed_this_cycle&out.working_identity_review_required).sum()),'n_B_context_controls_overlaid':0,'n_parent_changed_malignancy_context_recheck':int(changed_parent.sum()),'frozen_malignancy_axis_preserved_not_rerun':True,'index_hash':ih,'canonical_V7_changed':False,'final_fine_identity_accepted_this_cycle':False,'state_semantics':'Cluster-associated RNA context, not a per-cell functional state confirmation.','files':files}
(o/'complete.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({k:v for k,v in receipt.items() if k!='files'}))
