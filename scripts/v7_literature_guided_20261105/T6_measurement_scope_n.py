"""Separate unavailable marker panels from weak positive evidence in the review overlay."""
import os,json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd,h5py
from anndata.io import read_elem
ROOT=Path(os.environ.get('PROJECT_ROOT',str(globals().get('PROJECT_ROOT',''))));r=ROOT/'annotation_rebuild/20261105_v7_literature_guided_reannotation_v1';o=r/'stage_n_T6_measurement_scope_v1';o.mkdir(exist_ok=False)
q=pd.read_pickle(r/'stage_m_T6_working_overlay_v1/T6_hierarchical_working_review_PRIVATE.pkl.gz');gses=sorted(q.gse_id.astype(str).unique());rows=[]
for gs in gses:
 with h5py.File(ROOT/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects/harmonized_per_gse'/(gs+'_symbol_counts.h5ad'),'r') as f:v=pd.Index(read_elem(f['var']).index.astype(str))
 rows.append({'gse_id':gs,'CD3_component_measurable':any(g in v for g in ['CD3D','CD3E','CD3G']),'TCR_component_measurable':any(g in v for g in ['TRAC','TRBC1','TRBC2']),'FOXP3_measurable':'FOXP3' in v,'second_Treg_core_component_measurable':any(g in v for g in ['IL2RA','IKZF2','IKZF4'])})
s=pd.DataFrame(rows).set_index('gse_id');key=q.gse_id.astype(str)
for c in s:q[c]=key.map(s[c]).to_numpy().astype(bool)
q['T_parent_panel_evaluable']=q.CD3_component_measurable&q.TCR_component_measurable;q['Treg_positive_rule_evaluable']=q.T_parent_panel_evaluable&q.FOXP3_measurable&q.second_Treg_core_component_measurable
assert (~q.T_parent_cross_component_raw2|q.T_parent_panel_evaluable).all();assert (~q.Treg_core_positive_raw1|q.Treg_positive_rule_evaluable).all()
q['parent_support_status']=np.select([q.T_parent_cross_component_raw2,~q.T_parent_panel_evaluable],['positive_cross_component_RNA','marker_component_unavailable_at_source'],default='below_conservative_positive_rule_not_identity_negative')
q['parent_identity_evidence']='T_candidate_evaluate_with_positive_marker_and_source_scope';q['parent_identity']=np.where(q.T_parent_cross_component_raw2,'T_marker_supported_candidate','T_candidate_unsupported_by_this_positive_rule')
q['Treg_evidence_status']=np.select([q.Treg_core_positive_raw1,~q.Treg_positive_rule_evaluable],['positive_FOXP3_plus_second_core_and_T_components','positive_rule_not_evaluable_at_source'],default='not_positive_under_this_rule_not_Treg_negative')
take=q.qc_review.eq('low_marker_support_review')&~q.T_parent_panel_evaluable;q.loc[take,'qc_review']='source_marker_unavailability_review';q['version']='V7_T6_20261105_hierarchical_measured_scope_v2'
q.to_pickle(o/'T6_hierarchical_measured_review_v2_PRIVATE.pkl.gz');back=pd.read_pickle(o/'T6_hierarchical_measured_review_v2_PRIVATE.pkl.gz');pd.testing.assert_frame_equal(q,back)
counts=q.groupby(['parent_support_status','Treg_evidence_status'],observed=True).size().reset_index(name='n_cells');assert counts.n_cells.sum()==156631;counts.to_csv(o/'T6_marker_evaluability_counts.tsv',sep='\t',index=False)
manifest={'status':'completed_T6_marker_measurement_scope_and_overlay_readback','UTC':datetime.now(timezone.utc).isoformat(),'n_cells':len(q),'n_sources':len(gses),'n_T_parent_positive_rule_not_evaluable':int((~q.T_parent_panel_evaluable).sum()),'n_Treg_positive_rule_not_evaluable':int((~q.Treg_positive_rule_evaluable).sum()),'n_Treg_core_positive_raw1_unchanged':int(q.Treg_core_positive_raw1.sum()),'n_T_parent_positive_raw2_unchanged':int(q.T_parent_cross_component_raw2.sum()),'fine_final_labels_changed':0,'limits':'Unmeasured marker components are unknown, not negative. Other independent T evidence can still support major identity when this conservative cross-component rule cannot be evaluated.'}
(o/'complete.json').write_text(json.dumps(manifest,indent=2));print(json.dumps({'manifest':manifest,'counts':counts.to_dict('records')}))
