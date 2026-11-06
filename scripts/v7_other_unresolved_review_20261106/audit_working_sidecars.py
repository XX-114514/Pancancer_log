import json,hashlib
from datetime import datetime,timezone
import pandas as pd,numpy as np
r=PROJECT_ROOT/'runs/20261106_v7_other_unresolved_review_v1';o=r/'stage_f_working_annotations'
q=pd.read_pickle(r/'stage_c_QC_views/QC_decisions_PRIVATE.pkl.gz');t=pd.read_pickle(o/'four_models_working_annotation_PRIVATE.pkl.gz');routes=pd.read_pickle(o/'all_candidate_review_routes_PRIVATE.pkl.gz')
assert t.index.is_unique and len(t)==70778 and routes.index.equals(q.index)
assert t.old_final_annotation_v7.equals(q.loc[t.index,'final_annotation_v7'].astype(object).rename('old_final_annotation_v7'))
assert t.old_major_cell_type.equals(q.loc[t.index,'major_celltype_v7'].astype(object).rename('old_major_cell_type'))
assert q.loc[t.index,'current_lineage_model_eligible'].all()
assert not t.final_fine_identity_accepted.any() and not t.canonical_V7_changed.any() and not t.normal_HSC_identity_established.any()
assert set(t.index)==set(q.index[q.current_lineage_model_eligible])
assert routes.review_route.eq('QC_v3_operational_exclusion').equals(q.new_QC_exclude.rename('review_route'))
assert routes.review_route.eq('QC_v3_hematologic_context_holdout').equals(q.new_QC_holdout.rename('review_route'))
checks=[]
for name in ['Erythroid','Platelet','HSPC','ILC']:
 a=t[t.model_pool.eq(name)];b=pd.read_pickle(o/(name+'_working_annotation_PRIVATE.pkl.gz'));pd.testing.assert_frame_equal(a,b)
 md=pd.read_csv(r/('stage_d_'+name+'_refit')/'metadata_PRIVATE.tsv.gz',sep='\t',index_col='cell_id');assert a.index.equals(md.index)
 assert np.array_equal(a.diagnostic_cluster.to_numpy(),md['source_snn_res.0.4'].to_numpy())
 lab=pd.read_csv(r/'review_decisions'/(name+'_cluster_decisions.tsv'),sep='\t').set_index('cluster')
 for col in ['working_parent','identity_review_status','state_context','naming_markers','reason']:
  assert a['cluster_'+col].equals(a.diagnostic_cluster.map(lab[col]).rename('cluster_'+col))
 direct=a.assignment_basis.eq('positive_T_apparatus_route');assert a.loc[direct,'direct_T_positive_evidence'].all() and a.loc[direct,'working_parent'].eq('T').all()
 direct=a.assignment_basis.eq('positive_nonHB_structural_synthesis_route');assert a.loc[direct,'direct_Erythroid_positive_evidence'].all()
 checks.append({'lineage':name,'n_model_cells':len(a),'n_clusters':a.diagnostic_cluster.nunique(),'n_identity_review_required':int(a.identity_review_required.sum()),'n_broad_working_without_positive_competitor_review':int((~a.identity_review_required).sum()),'n_positive_T_routed':int(a.assignment_basis.eq('positive_T_apparatus_route').sum()),'n_direct_nonHB_Ery_routed':int(a.assignment_basis.eq('positive_nonHB_structural_synthesis_route').sum()),'n_ILC_like_candidates':int(a.working_parent.eq('ILC_like_candidate').sum()),'normal_HSC_identity_established':False,'fine_identity_accepted':False})
saved=pd.read_csv(o/'working_identity_counts.tsv',sep='\t');actual=t.groupby(['model_pool','working_parent','assignment_basis','identity_review_required','hematologic_disease_context'],dropna=False).size().rename('n_cells').reset_index();pd.testing.assert_frame_equal(saved.sort_values(list(actual.columns[:-1])).reset_index(drop=True),actual.sort_values(list(actual.columns[:-1])).reset_index(drop=True))
manifest=json.loads((o/'complete.json').read_text())
for z in manifest['files']:
 p=o/z['name'];assert hashlib.sha256(p.read_bytes()).hexdigest()==z['sha256'] and p.stat().st_size==z['bytes']
m=np.load(r/'stage_c_QC_views/global_QC_view_v3_masks_PRIVATE.npz');assert int(m['primary_exclude'].sum())==2916 and int(m['QC_holdout'].sum())==238 and int(m['primary_keep'].sum())==4673633
frozen=PROJECT_ROOT/'runs/20260818_pancancer_5M_merge_prep/11_v7_cnv_freeze/20260905_v7_cnv_rerun_freeze_v1/checkpoints/final_annotation_v7.pkl.gz'
expected='e4fda365b2bfa9ea6f82b9b48c4ce6009a0bf31dcfb5a57e7e94cab6d8d3981e';h=hashlib.sha256()
with frozen.open('rb') as f:
 for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
assert h.hexdigest()==expected
out={'status':'PASS','UTC':datetime.now(timezone.utc).isoformat(),'n_review_scope':len(q),'n_model_cells':len(t),'n_diagnostic_clusters':int(t.groupby('model_pool').diagnostic_cluster.nunique().sum()),'exact_old_labels_and_axes_preserved':True,'four_models_disjoint_and_exactly_eligible':True,'all_aggregate_counts_recomputed':True,'all_sidecar_SHA256_verified':True,'frozen_V7_SHA256_unchanged':h.hexdigest(),'QC_v3_exclude':2916,'QC_v3_holdout':238,'QC_v3_keep':4673633,'models':checks,'final_fine_annotation_accepted':False}
(o/'independent_sidecar_audit.json').open('x').write(json.dumps(out,indent=2));print(json.dumps(out))
