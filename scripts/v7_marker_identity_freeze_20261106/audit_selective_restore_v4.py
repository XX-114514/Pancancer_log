"""Independently verify selective restoration against the last subtype-preserving source."""
import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';o=r/'stage_h_selective_subtype_restore_v4'
q=pd.read_pickle(PROJECT_ROOT/'runs/20261106_v7_B_unresolved_review_v1/stage_g_consolidated_annotation/current_annotation_all_cells_PRIVATE.pkl.gz');b=pd.read_pickle(r/'stage_c4_harmonized_identity_freeze/annotation_all_cells_with_release_status_PRIVATE.pkl.gz');t=pd.read_pickle(o/'all_axis_selective_annotation_with_history_PRIVATE.pkl.gz');f=pd.read_pickle(o/'annotation_frozen_V7_selective_subtype_20261106_v4_PRIVATE.pkl.gz');assert q.index.equals(b.index) and q.index.equals(t.index)
for col in q:
 if not col.startswith('frozen_'):continue
 x=q[col].astype(object).to_numpy();y=t[col].astype(object).to_numpy();na=pd.isna(x);assert np.array_equal(na,pd.isna(y)) and np.array_equal(x[~na],y[~na]),col
assert np.array_equal(q.current_primary_QC_keep.to_numpy(),t.QC_keep.to_numpy());assert np.array_equal(q.current_QC_action.astype(str).to_numpy(),t.QC_action.astype(str).to_numpy())
names=q.frozen_label.astype(str);unfinished=names.str.contains(r'(?i)unresolved|coarse|provisional|pending|boundary|unknown|ambiguous',regex=True).to_numpy();qc=q.current_primary_QC_keep.to_numpy(bool);completed_rare_refit=q.annotation_layer.astype(str).eq('rare_four_working_20261106').to_numpy();parent_supported=q.current_parent.astype(str).eq(q.frozen_parent.astype(str)).to_numpy() & q.identity_review_status.astype(str).eq('broad_working').to_numpy();exception=completed_rare_refit&~unfinished&~parent_supported;restore=qc&~unfinished&~exception;fallback=qc&~restore&b.annotation_frozen.to_numpy(bool);selected=restore|fallback
assert np.array_equal(selected,t.annotation_frozen.to_numpy());assert np.array_equal(restore,t.prior_named_annotation_retained.to_numpy());assert np.array_equal(exception,t.previous_refit_contradicts_named_parent.to_numpy())
assert np.array_equal(t.loc[restore,'final_annotation'].astype(str).to_numpy(),names.loc[restore].to_numpy());assert np.array_equal(t.loc[restore,'final_parent'].astype(str).to_numpy(),q.loc[restore,'frozen_parent'].astype(str).to_numpy())
assert np.array_equal(t.loc[fallback,'final_annotation'].astype(str).to_numpy(),b.loc[fallback,'final_annotation'].astype(str).to_numpy());assert np.array_equal(t.loc[fallback,'final_parent'].astype(str).to_numpy(),b.loc[fallback,'final_parent'].astype(str).to_numpy());assert not (restore&unfinished).any();assert not (fallback&~unfinished&~exception).any();assert not t.loc[selected,'final_annotation'].astype(str).str.contains('unresolved|coarse|provisional|pending|boundary|unknown|ambiguous',case=False).any()
assert f.index.equals(t.index[selected]);assert np.array_equal(f.final_annotation.astype(str).to_numpy(),t.loc[selected,'final_annotation'].astype(str).to_numpy());assert np.array_equal(f.major_cell_type.astype(str).to_numpy(),t.loc[selected,'final_parent'].astype(str).to_numpy())
assert not {'frozen_label','legacy_label','old_label','freeze_status','diagnostic_group'}.intersection(f.columns)
for col in q:
 if col.startswith('frozen_') and ('malignan' in col or 'evidence_version' in col):
  x=q.loc[selected,col].astype(object).to_numpy();y=f[col.removeprefix('frozen_')].astype(object).to_numpy();na=pd.isna(x);assert np.array_equal(na,pd.isna(y)) and np.array_equal(x[~na],y[~na]),col
for col in ['raw_IFN_programme_positive','raw_cycling_programme_positive','raw_heat_shock_programme_positive']:assert np.array_equal(b.loc[selected,col].to_numpy(),f[col].to_numpy())
pos=0
for chunk in pd.read_csv(o/'annotation_frozen_V7_selective_subtype_20261106_v4_PRIVATE.tsv.gz',sep='\t',usecols=['cell_id','final_annotation','major_cell_type'],chunksize=100000):
 n=len(chunk);want=f.iloc[pos:pos+n];assert np.array_equal(chunk.cell_id.to_numpy(),want.index.to_numpy())
 for col in ['final_annotation','major_cell_type']:assert np.array_equal(chunk[col].to_numpy(),want[col].astype(str).to_numpy())
 pos+=n
assert pos==len(f)
saved=pd.read_csv(o/'corrected_final_label_counts.tsv',sep='\t');expected=f.groupby(['major_cell_type','final_annotation','annotation_resolution','annotation_application'],observed=True,dropna=False).size().rename('n_cells').reset_index();keys=['major_cell_type','final_annotation','annotation_resolution','annotation_application'];saved=saved.sort_values(keys);expected=expected.sort_values(keys)
for col in keys:assert np.array_equal(saved[col].astype(str).to_numpy(),expected[col].astype(str).to_numpy())
assert np.array_equal(saved.n_cells.to_numpy(),expected.n_cells.to_numpy());assert pd.read_csv(o/'all_previous_labels_correction_crosswalk.tsv',sep='\t').n_cells.sum()==len(q)
receipt=json.loads((o/'complete.json').read_text())
for item in receipt['files']:
 p=o/item['name'];assert p.stat().st_size==item['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==item['sha256']
out={'status':'PASS_selective_subtype_preservation_and_complete_export','UTC':datetime.now(timezone.utc).isoformat(),'n_all_axis':len(t),'n_selected':len(f),'n_labels':int(f.final_annotation.nunique()),'n_retained_previous_names':int(restore.sum()),'n_fallback_major_marker':int(fallback.sum()),'n_QC_keep_outside':int((qc&~selected).sum()),'all_prior_names_outside_explicit_refit_exceptions_preserved_exactly':True,'fallback_restricted_to_unfinished_labels_or_explicit_prior_refit_exceptions':True,'no_broad_gate_applied_to_unchanged_named_subtypes':True,'all_QC_original_fields_CNV_state_and_axis_preserved':True,'private_TSV_all_rows_and_labels_verified':True,'all_aggregate_counts_and_artifact_SHA_verified':True,'new_subtype_validation_claims':0,'all_atlas_publication_acceptance':False};(o/'independent_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
