"""Review accepted full-scope cluster proposals before shared-identity modeling.

Reuse independently rebuilt native profiles; no counts copy, relabeling or
clinical fitting. Source-level summaries and anonymous unit diagnostics remain
separate from independent patients and biological accuracy.
"""
import argparse,collections,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from TME_all_lineage_cluster_policy_v1 import GATES
from TME_cluster_annotation_policy_v1 import STATE_PANELS

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')

ROUTES={'T_NK':{'T','NK','T_NK'},'B_Plasma':{'B','Plasma'},'Myeloid':{'Myeloid'},'Endothelial':{'Endothelial'},'Fibroblast':{'Fibroblast'},'Mural':{'Mural'},'Epithelial':{'Epithelial'},'Neural':{'Glial','Astrocyte','Neuron'},'Mast':{'Mast'},'Erythroid':{'Erythroid'},'Platelet':{'Platelet'},'unresolved':{'unresolved'}}
MARKER_ROLES={'erythroid_associated_leading_marker':set('HBA1 HBA2 HBB HBD HBG1 HBG2 HBM HBZ AHSP ALAS2 GYPA SLC4A1'.split()),'platelet_associated_leading_marker':set('PF4 PPBP GP9 GP1BA GP1BB ITGA2B'.split()),'heatshock_processing_or_stress_marker':set('HSPA1A HSPA1B DNAJB1 HSPH1 HSP90AA1 HSPA6'.split()),'immediate_early_program_marker':set('FOS FOSB JUN JUNB JUND DUSP1 EGR1 ATF3'.split()),'broad_RNA_content_marker':set('MALAT1 NEAT1 ACTB GAPDH TMSB4X TMSB10 EEF1A1'.split())}

def main(root):
 begin=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);A=Path(cfg['auditor']);D=Path(cfg['delivery'])
 for r in [P,A,D]:assert load(r/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json');pt=load(P/'provenance/terminal_v1.json');pc=load(P/'config/spec_v1.json');assert pc['root_protocol']=='qualified_full_TME_source_roots_v1';assert sha(P/'outputs/cluster_hierarchy_named_catalog_v1.json')==pt['catalog_SHA256'] and sha(P/'outputs/cluster_hierarchy_index_v1.json')==pt['index_SHA256']
 nodes=load(P/'outputs/cluster_hierarchy_named_catalog_v1.json')['nodes'];leaves=[r for r in nodes if r['final_leaf']];by={(r['source'],r['source_node_code']):r for r in nodes};assert len(by)==len(nodes) and sum(r['n_records'] for r in leaves)==pt['totals']['n_pool_records']
 accepted=load(Path(pc['working_delivery'])/'outputs/accepted_TME_scoped_working_annotation_v1.json');assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];w=load(accepted['producer_manifest']['path']);assert sha(w['index']['path'])==w['index']['SHA256'];wb=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])}
 root_cfg=load(Path(pc['root_model'])/'config/spec_v1.json');feature_inputs={r['source']:r for r in root_cfg['frozen_preflight']['sources']};excluded=set(g for modules in GATES.values() for module in modules for g in module)|set(g for panel in STATE_PANELS.values() for g in panel);features={};source_summary=[];crosswalk=[];group_review=[];route_counts=collections.defaultdict(dict);all_counts=collections.Counter();role_totals=collections.defaultdict(collections.Counter);states=collections.defaultdict(collections.Counter)
 for e in load(P/'outputs/cluster_hierarchy_index_v1.json'):
  source=e['source'];mp=P/'outputs'/e['metadata'];assert sha(mp)==e['metadata_SHA256'];m=load(mp);assert sha(m['arrays']['path'])==m['arrays']['SHA256'];fp=Path(feature_inputs[source]['candidate_feature_axis_path']);assert sha(fp)==feature_inputs[source]['candidate_feature_axis_SHA256'];fi=load(fp);assert fi['feature_selection_not_performed'];features[source]=set(fi['genes'])&set(m['native_measured_genes'])
  with np.load(m['arrays']['path'],allow_pickle=False) as f:pool=f['source_scope_position'];final=f['final_leaf_node']
  wr=wi[source];assert sha(wb/wr['arrays'])==wr['arrays_SHA256'];assert sha(wb/wr['metadata'])==wr['metadata_SHA256'];wm=load(wb/wr['metadata'])
  with np.load(wb/wr['arrays'],allow_pickle=False) as f:
   np.testing.assert_array_equal(pool,np.flatnonzero(f['fit_pool_candidate_not_final_QC']));old=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[f['consensus_celltype_code'][pool]];n_all=len(f['native_row'])
  proposed=np.asarray([by[source,int(j)]['cluster_RNA_proposal']['celltype'] for j in final],object);trans=collections.Counter(zip(map(str,old),map(str,proposed)))
  for (a,b),n in sorted(trans.items()):crosswalk.append(dict(source=source,baseline_v5_celltype=a,new_cluster_RNA_proposal=b,canonical_cluster_proposal={'Tcell':'T','NKcell':'NK','gamma_delta_T_RNA_candidate':'T_gamma_delta_RNA_candidate'}.get(b,b),n_fit_records=n,not_identity_accuracy=True))
  audit=load(A/'outputs'/(source+'_cluster_hierarchy_audit_v1.json'));assert audit['producer_metadata_SHA256']==sha(mp);dominance={r['source_node_code']:r for r in audit['anonymous_unit_dominance']};members=[r for r in leaves if r['source']==source];assert len(members)==m['n_final_leaves'];assert sum(r['n_records'] for r in members)==len(pool);counts=collections.Counter();typ=collections.Counter();statuses=collections.Counter();sr=collections.Counter()
  for r in members:
   call=r['cluster_RNA_proposal'];major=call['major'];counts[major]+=r['n_records'];typ[call['celltype']]+=r['n_records'];statuses[call['status']]+=r['n_records'];role=[k for k,v in MARKER_ROLES.items() if (r['marker'] or '').upper() in v];dom=dominance[r['source_node_code']];own_specific=r['marker_evidence_type']=='own_subset_expression_specificity_not_inferential_DEG';state_names=[k for k,v in r['descriptive_state_programs'].items() if v['descriptively_enriched']];unmeasured_states=[k for k,v in r['descriptive_state_programs'].items() if v['n_measured_genes']<max(2,int(np.ceil(v['n_panel_genes']*.5)))];unit_limited=dom['n_declared_units_with_at_least_20_cells']<2;concentrated=dom['maximum_declared_unit_fraction'] is not None and dom['maximum_declared_unit_fraction']>=.8
   for k in role:role_totals[k].update(n_clusters=1,n_records=r['n_records'])
   for k in state_names:states[k].update(n_clusters=1,n_records=r['n_records'])
   sr.update(n_clusters=1,n_records=r['n_records'],n_records_without_own_positive_specificity=0 if own_specific else r['n_records'],n_records_anonymous_unit_limited=r['n_records'] if unit_limited else 0,n_records_anonymous_unit_concentrated=r['n_records'] if concentrated else 0)
   group_review.append(dict(source=source,subcluster=r['display_label'],stable_membership_SHA256=r['stable_membership_SHA256'],n_records=r['n_records'],major_proposal=major,celltype_proposal=call['celltype'],annotation_status=call['status'],T_evidence_mode=r.get('assay_aware_T_evidence_mode'),own_positive_marker_specificity=own_specific,marker_role_flags=role,state_program_names=state_names,state_programs_unassessable=unmeasured_states,anonymous_expression_unit_limited=unit_limited,anonymous_expression_unit_concentrated=concentrated,does_not_diagnose_artifact_or_independent_patient_support=True))
  assert sum(counts.values())==len(pool);source_summary.append(dict(source=source,n_current_qualified_records=n_all,n_fit_records=len(pool),n_unfitted_retained=n_all-len(pool),n_leaf_groups=len(members),major_proposal_record_counts=dict(counts),celltype_proposal_record_counts=dict(typ),annotation_status_record_counts=dict(statuses),review_counts=dict(sr)))
  all_counts.update(sr)
  for route,majors in ROUTES.items():route_counts[route][source]=sum(n for k,n in counts.items() if k in majors)
 assert all_counts['n_records']==pt['totals']['n_pool_records'];assert sum(r['n_current_qualified_records'] for r in source_summary)==pc['frozen_preflight']['n_current_qualified_records'];assert sum(r['n_unfitted_retained'] for r in source_summary)==pc['frozen_preflight']['n_unfitted_retained'];assert sum(r['n_fit_records'] for r in crosswalk)==all_counts['n_records']
 feature_review=[]
 for route,counts in ROUTES.items():
  available=route_counts[route];eligible=sorted(s for s,n in available.items() if n>=cfg['minimum_source_route_records']);observed=sorted(s for s,n in available.items() if n>0);common=set.intersection(*(features[s] for s in eligible)) if eligible else set();model=sorted(common-excluded);enough=len(eligible)>=cfg['minimum_fit_sources'] and len(model)>=cfg['minimum_common_model_features'];feature_review.append(dict(route=route,n_sources_observed=len(observed),n_sources_meeting_frozen_size=len(eligible),source_record_counts=available,eligible_source_ids=eligible,n_common_unique_measured_candidate_features=len(common),n_shared_model_features_after_identity_state_holdout=len(model),model_feature_genes_private=model,identity_and_state_panel_genes_excluded=True,passes_preliminary_shared_reference_size_and_measurement_gate=enough,not_a_model_fit_or_production_choice=True,source_alias_patient_and_model_state_gates_still_separate=True))
 summary=dict(schema='full_TME_cluster_review_and_shared_reference_preflight_v1',observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope='all_current_source_qualified_fit_records_retained_unfitted_outside_model',totals=dict(all_counts),n_qualified_records=sum(r['n_current_qualified_records'] for r in source_summary),n_unfitted_retained=sum(r['n_unfitted_retained'] for r in source_summary),n_sources=len(source_summary),n_final_leaf_clusters=len(leaves),source_summary=source_summary,marker_role_counts={k:dict(v) for k,v in role_totals.items()},descriptive_state_counts={k:dict(v) for k,v in states.items()},shared_reference_preflight=[{k:v for k,v in r.items() if k not in ['model_feature_genes_private','eligible_source_ids','source_record_counts']} for r in feature_review],new_CM_started=False,current_working_v5_unchanged=True,not_final_shared_taxonomy=True,not_orthogonal_accuracy_or_artifact_diagnosis=True)
 save(root/'outputs/TME_full_cluster_review_summary_v1.json',summary);save(root/'outputs/TME_shared_native_feature_preflight_v1.json',dict(frozen_minimum_source_route_records=cfg['minimum_source_route_records'],frozen_minimum_fit_sources=cfg['minimum_fit_sources'],frozen_minimum_common_model_features=cfg['minimum_common_model_features'],routes=feature_review));save(root/'outputs/TME_full_cluster_group_review_v1.json',group_review);save(root/'outputs/TME_full_cluster_v5_crosswalk_v1.json',crosswalk)
 save(root/'provenance/terminal_v1.json',dict(status='FULL_SCOPE_CLUSTER_REVIEW_AND_MEASUREMENT_PREFLIGHT_COMPLETE_NOT_TAXONOMY',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),auditor_terminal_SHA256=sha(A/'provenance/terminal_v1.json'),delivery_terminal_SHA256=sha(D/'provenance/terminal_v1.json'),files={p.name:sha(p) for p in (root/'outputs').iterdir()},n_fit_records=all_counts['n_records'],n_leaf_groups=len(leaves),current_labels_counts_states_and_clusters_unchanged=True,no_shared_model_fit_or_new_CM=True,wall_seconds=time.monotonic()-begin))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
