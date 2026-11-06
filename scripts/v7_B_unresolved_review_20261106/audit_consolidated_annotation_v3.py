import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd

r=PROJECT_ROOT;run=r/'runs/20261106_v7_B_unresolved_review_v1';o=run/'stage_g_consolidated_annotation'
vpath=r/'runs/20260818_pancancer_5M_merge_prep/11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1/checkpoints/final_annotation_v7.pkl.gz'
v=pd.read_pickle(vpath);t=pd.read_pickle(o/'current_annotation_all_cells_PRIVATE.pkl.gz');m=np.load(run/'stage_c_QC_views/global_QC_view_v4_masks_PRIVATE.npz')
assert t.index.equals(v.index) and len(t)==4676787 and t.index.is_unique
assert np.array_equal(t.frozen_label.astype(str).to_numpy(),v.final_annotation_v7.astype(str).to_numpy())
assert np.array_equal(t.frozen_parent.astype(str).to_numpy(),v.major_celltype_v7.astype(str).to_numpy())
assert np.array_equal(t.frozen_strict_major_eligible.to_numpy(),v.strict_major_eligible_v7.fillna(False).astype(bool).to_numpy())
assert np.array_equal(t.frozen_strict_subtype_eligible.to_numpy(),v.strict_subtype_eligible_v7.fillna(False).astype(bool).to_numpy())
assert np.array_equal(t.current_primary_QC_keep.to_numpy(),m['primary_keep'])
assert np.array_equal(t.current_QC_action.eq('exclude_primary_QC').to_numpy(),m['primary_exclude'])
assert np.array_equal(t.current_QC_action.eq('QC_holdout').to_numpy(),m['QC_holdout'])
for col in ['malignancy_call_v7','malignancy_tier_v7','is_malignant_primary_v7','is_malignant_high_confidence_v7','is_malignant_sensitive_v7','strict_malignant_downstream_eligible_v7','malignancy_evidence_version_v7']:
 pd.testing.assert_series_equal(t['frozen_'+col].astype(object).rename(col),v[col].astype(object),check_dtype=False,check_categorical=False)
mal=pd.read_csv(o/'frozen_malignancy_axis_current_context_counts.tsv',sep='\t');assert mal.n_cells.sum()==len(v)
assert not t.fine_identity_final_accepted_this_cycle.any()
expected=[];old=r/'runs/20261106_v7_QC_clean_reannotation_v1/stage_h_working_review_overlays_v1'
for name in ['Mast','T6']:
 q=pd.read_pickle(old/(name+'_working_annotation_QC_PRIVATE.pkl.gz'));expected.extend(q.index[q.clean_lineage_model_eligible])
rare=pd.read_pickle(r/'runs/20261106_v7_other_unresolved_review_v1/stage_f_working_annotations/four_models_working_annotation_PRIVATE.pkl.gz');expected.extend(rare.index)
b=pd.read_pickle(run/'stage_f_B_working_annotations/B_targets_working_annotation_PRIVATE.pkl.gz');expected.extend(b.index[b.in_current_B_model])
idx=pd.Index(expected);assert idx.is_unique and set(idx)==set(t.index[t.refit_identity_reviewed_this_cycle])
assert t.loc[idx,'current_primary_QC_keep'].all()
q=pd.read_pickle(run/'stage_c_QC_views/QC_decisions_PRIVATE.pkl.gz');controls=q.index[q.current_lineage_model_eligible&~q.is_annotation_target]
assert not t.loc[controls,'annotation_layer'].eq('B_working_20261106').any()
assert set(t.index[t.annotation_layer.eq('B_working_20261106')])==set(b.index[b.in_current_B_model])
unrefit=t.index[~t.refit_identity_reviewed_this_cycle]
assert np.array_equal(t.loc[unrefit,'current_annotation_label'].astype(str).to_numpy(),v.loc[unrefit,'final_annotation_v7'].astype(str).to_numpy())
cols=['current_parent','current_annotation_label','annotation_layer','identity_review_status','candidate_subtype','cluster_RNA_context','current_QC_action']
actual=t.groupby(cols,observed=True,dropna=False).size().rename('n_cells').reset_index()
saved=pd.read_csv(o/'current_annotation_label_counts.tsv',sep='\t')
for col in cols:actual[col]=actual[col].astype(str);saved[col]=saved[col].astype(str)
pd.testing.assert_frame_equal(actual.sort_values(cols).reset_index(drop=True),saved.sort_values(cols).reset_index(drop=True),check_dtype=False)
assert saved.n_cells.sum()==len(v)
par=pd.read_csv(o/'current_parent_summary.tsv',sep='\t');assert par.n_all_records.sum()==len(v) and par.n_primary_QC_keep.sum()==int(m['primary_keep'].sum())
assert par.n_QC_excluded.sum()==int(m['primary_exclude'].sum()) and par.n_QC_holdout.sum()==int(m['QC_holdout'].sum())
assert par.n_refit_reviewed_kept.sum()==len(idx) and np.array_equal((par.n_broad_working_kept+par.n_refit_identity_review_required).to_numpy(),par.n_refit_reviewed_kept.to_numpy())
f=pd.read_csv(o/'all_106_frozen_labels_current_status.tsv',sep='\t');assert len(f)==106 and f.n_frozen.sum()==len(v)
assert np.array_equal((f.n_primary_QC_keep+f.n_QC_excluded+f.n_QC_holdout).to_numpy(),f.n_frozen.to_numpy())
conf=t.groupby(['frozen_label','frozen_confidence'],observed=True,dropna=False).size().rename('n_cells').reset_index()
saved_conf=pd.read_csv(o/'frozen_confidence_counts.tsv',sep='\t')
for col in ['frozen_label','frozen_confidence']:conf[col]=conf[col].astype(str);saved_conf[col]=saved_conf[col].astype(str)
pd.testing.assert_frame_equal(conf.sort_values(['frozen_label','frozen_confidence']).reset_index(drop=True),saved_conf.sort_values(['frozen_label','frozen_confidence']).reset_index(drop=True),check_dtype=False)
ci=pd.read_csv(o/'all_reviewed_cluster_inventory.tsv',sep='\t');assert ci.n_target_cells.sum()==len(idx) and ci.n_context_controls.sum()==len(controls) and np.array_equal((ci.n_target_cells+ci.n_context_controls).to_numpy(),ci.n_model_cells.to_numpy())
for row in ci.to_dict('records'):
 model='B_unresolved_targets' if row['source_model']=='B_unresolved' else row['source_model'];take=t.diagnostic_model.eq(model)&t.diagnostic_cluster.eq(row['cluster_id']).fillna(False)&t.refit_identity_reviewed_this_cycle;assert int(take.sum())==row['n_target_cells']
pos=0
for chunk in pd.read_csv(o/'current_annotation_all_cells_PRIVATE.tsv.gz',sep='\t',usecols=['cell_id','frozen_label','current_annotation_label','current_QC_action'],chunksize=100000):
 n=len(chunk);want=t.iloc[pos:pos+n];assert np.array_equal(chunk.cell_id.to_numpy(),want.index.to_numpy())
 for col in ['frozen_label','current_annotation_label','current_QC_action']:assert np.array_equal(chunk[col].astype(str).to_numpy(),want[col].astype(str).to_numpy())
 pos+=n
assert pos==len(v)
manifest=json.loads((o/'complete.json').read_text())
for z in manifest['files']:
 p=o/z['name'];assert p.stat().st_size==z['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==z['sha256']
h=hashlib.sha256()
with vpath.open('rb') as f:
 for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
assert h.hexdigest()=='e4fda365b2bfa9ea6f82b9b48c4ce6009a0bf31dcfb5a57e7e94cab6d8d3981e'
receipt={'status':'PASS','UTC':datetime.now(timezone.utc).isoformat(),'n_frozen_cells':len(v),'n_frozen_labels':106,'n_refit_reviewed_targets':len(idx),'n_primary_QC_keep':int(m['primary_keep'].sum()),'n_QC_excluded':int(m['primary_exclude'].sum()),'n_QC_holdout':int(m['QC_holdout'].sum()),'all_original_labels_and_strict_flags_preserved':True,'all_frozen_malignancy_columns_preserved_not_rerun':True,'all_overlay_target_sets_disjoint_and_exact':True,'B_reference_controls_not_overlaid':len(controls),'all_public_aggregate_counts_recomputed':True,'full_private_TSV_stream_axis_and_values_verified':True,'all_private_artifact_SHA_verified':True,'frozen_V7_SHA256_unchanged':h.hexdigest(),'fine_identity_final_accepted_this_cycle':False}
(o/'independent_consolidation_audit.json').open('x').write(json.dumps(receipt,indent=2));print(json.dumps(receipt))
