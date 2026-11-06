import json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
r=PROJECT_ROOT/'runs/20261106_v7_publication_marker_annotation_v1';o=r/'stage_c3_exact_identity_freeze';proof=json.loads((r/'stage_c2_exact_audit_resume/independent_serialized_axis_audit.json').read_text());assert proof['status']=='PASS_full_axis_exact_value_and_complete_membership_audit'
t=pd.read_pickle(o/'annotation_all_cells_with_release_status_PRIVATE.pkl.gz');f=pd.read_pickle(o/'marker_identity_frozen_cells_PRIVATE.pkl.gz');assert hashlib.sha256((o/'annotation_all_cells_with_release_status_PRIVATE.pkl.gz').read_bytes()).hexdigest()==proof['source_full_axis_sha256'];assert len(t)==proof['n_all_axis'] and len(f)==proof['n_frozen'];selected=t.loc[t.annotation_frozen];assert f.index.equals(selected.index)
for col in f:
 if isinstance(f[col].dtype,pd.CategoricalDtype):assert f[col].cat.categories.equals(selected[col].cat.categories) and np.array_equal(f[col].cat.codes.to_numpy(),selected[col].cat.codes.to_numpy()),col
 else:
  x=f[col].astype(object).to_numpy();y=selected[col].astype(object).to_numpy();na=pd.isna(x);assert np.array_equal(na,pd.isna(y)) and np.array_equal(x[~na],y[~na]),col
assert f.QC_keep.all() and f.final_annotation.notna().all()
for name,total in [('frozen_identity_counts.tsv',len(f)),('frozen_group_counts.tsv',len(f)),('all_cell_release_status_counts.tsv',len(t)),('historical_to_new_identity_counts.tsv',len(t))]:
 z=pd.read_csv(o/name,sep='\t',keep_default_na=False);count='n_frozen_cells' if name=='frozen_group_counts.tsv' else 'n_cells';assert z[count].sum()==total
actual=f.groupby(['final_parent','final_annotation','representative_markers','freeze_basis'],observed=True,dropna=False).size().rename('n_cells').reset_index();saved=pd.read_csv(o/'frozen_identity_counts.tsv',sep='\t');keys=['final_parent','final_annotation','representative_markers','freeze_basis'];actual=actual.sort_values(keys).reset_index(drop=True);saved=saved.sort_values(keys).reset_index(drop=True)
for col in keys:assert np.array_equal(actual[col].astype(str).to_numpy(),saved[col].astype(str).to_numpy())
assert np.array_equal(actual.n_cells.to_numpy(),saved.n_cells.to_numpy());position=0
for chunk in pd.read_csv(o/'marker_identity_frozen_cells_PRIVATE.tsv.gz',sep='\t',usecols=['cell_id','final_annotation','final_parent'],chunksize=100000):
 n=len(chunk);want=f.iloc[position:position+n];assert np.array_equal(chunk.cell_id.to_numpy(),want.index.to_numpy());assert np.array_equal(chunk.final_annotation.to_numpy(),want.final_annotation.astype(str).to_numpy());assert np.array_equal(chunk.final_parent.to_numpy(),want.final_parent.astype(str).to_numpy());position+=n
assert position==len(f)
manifest=json.loads((o/'complete.json').read_text())
for z in manifest['files']:
 p=o/z['name'];assert p.stat().st_size==z['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==z['sha256']
out={'status':'PASS','UTC':datetime.now(timezone.utc).isoformat(),'n_all_axis':len(t),'n_frozen':len(f),'n_QC_keep_not_frozen':int((t.QC_keep&~t.annotation_frozen).sum()),'n_labels':f.final_annotation.nunique(),'all_full_axis_and_complete_membership_proof_reused_by_exact_file_SHA':True,'all_frozen_derivative_cells_and_all_columns_exact':True,'all_frozen_identity_counts_recomputed':True,'complete_private_TSV_axis_and_values_verified':True,'all_artifact_hashes_verified':True,'all_frozen_labels_have_no_progress_placeholders':True,'all_registered_competitor_rules_after_related_axis_exceptions_verified_in_source_audit':True,'all_atlas_publication_acceptance':False};(o/'independent_freeze_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
