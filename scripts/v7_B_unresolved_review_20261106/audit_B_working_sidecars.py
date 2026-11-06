import json,hashlib
from datetime import datetime,timezone
import pandas as pd,numpy as np
r=PROJECT_ROOT/'runs/20261106_v7_B_unresolved_review_v1';o=r/'stage_f_B_working_annotations';t=pd.read_pickle(o/'B_targets_working_annotation_PRIVATE.pkl.gz');q=pd.read_pickle(r/'stage_c_QC_views/QC_decisions_PRIVATE.pkl.gz');md=pd.read_csv(r/'stage_d_B_unresolved_refit/metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id')
assert t.index.equals(q.index[q.is_annotation_target]) and t.index.is_unique and len(t)==158257
assert np.array_equal(t.old_final_annotation_v7.astype(str).to_numpy(),q.loc[t.index,'final_annotation_v7'].astype(str).to_numpy())
assert np.array_equal(t.primary_QC_keep.to_numpy(),q.loc[t.index,'current_QC_keep'].to_numpy())
ids=t.index[t.in_current_B_model];assert set(ids)==set(md.index[q.loc[md.index,'is_annotation_target']]) and len(ids)==157888
assert np.array_equal(t.loc[ids,'diagnostic_cluster'].astype(int).to_numpy(),md.loc[ids,'source_snn_res.0.4'].astype(int).to_numpy())
lab=pd.read_csv(r/'review_decisions/B_cluster_decisions.tsv',sep='\t').set_index('cluster')
for c in lab.index:assert int(t.diagnostic_cluster.eq(c).fillna(False).sum())==int((md['source_snn_res.0.4'].astype(int).eq(c)&q.loc[md.index,'is_annotation_target']).sum())
assert not t.fine_identity_final_accepted.any() and not t.canonical_V7_changed.any()
assert not t.loc[~t.in_current_B_model,'primary_QC_keep'].any() and t.loc[ids,'primary_QC_keep'].all()
for basis,col in [('positive_B_distinct_component_route','B_distinct_components_positive'),('positive_T_apparatus_route','T_positive_source_aware')]:
 ix=t.index[t.assignment_basis.eq(basis)];assert q.loc[ix,col].all()
saved=pd.read_csv(o/'B_working_identity_counts.tsv',sep='\t');keys=['working_parent','working_annotation_label','working_subtype_candidate','identity_review_required','primary_QC_keep'];actual=t.groupby(keys,observed=True,dropna=False).size().rename('n_cells').reset_index();pd.testing.assert_frame_equal(actual.sort_values(keys).reset_index(drop=True),saved.sort_values(keys).reset_index(drop=True),check_dtype=False)
receipt=json.loads((o/'complete.json').read_text());assert hashlib.sha256((o/'B_targets_working_annotation_PRIVATE.pkl.gz').read_bytes()).hexdigest()==receipt['sidecar_sha256']
out={'status':'PASS','UTC':datetime.now(timezone.utc).isoformat(),'n_target_records':len(t),'n_model_targets':len(ids),'n_excluded_targets':int((~t.in_current_B_model).sum()),'n_broad_working_targets':int((t.in_current_B_model&~t.identity_review_required).sum()),'n_parent_review_targets':int((t.in_current_B_model&t.identity_review_required).sum()),'n_controls_overlaid':0,'all_target_axes_and_old_labels_preserved':True,'all_public_identity_counts_recomputed':True,'fine_identity_final_accepted':False};(o/'independent_B_sidecar_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
