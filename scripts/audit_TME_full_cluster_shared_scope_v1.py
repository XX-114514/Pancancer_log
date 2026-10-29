"""Independent denominators, version crosswalk and measurement-gate replay.

This accepts group/source arithmetic and feature eligibility, not identity
accuracy, independent patients, state structure or a selected shared model.
"""
import argparse,collections,datetime,hashlib,json,time
from pathlib import Path
import numpy as np
from TME_all_lineage_cluster_policy_v1 import GATES
from TME_cluster_annotation_policy_v1 import STATE_PANELS

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')

def main(root):
 cfg=load(root/'config/spec_v1.json');R=Path(cfg['reviewer']);assert load(R/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;rc=load(R/'config/spec_v1.json');P=Path(rc['producer']);A=Path(rc['auditor']);rt=load(R/'provenance/terminal_v1.json');assert rt['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json')
 for name,digest in rt['files'].items():assert sha(R/'outputs'/name)==digest
 summary=load(R/'outputs/TME_full_cluster_review_summary_v1.json');groups=load(R/'outputs/TME_full_cluster_group_review_v1.json');cross=load(R/'outputs/TME_full_cluster_v5_crosswalk_v1.json');pre=load(R/'outputs/TME_shared_native_feature_preflight_v1.json');nodes=load(P/'outputs/cluster_hierarchy_named_catalog_v1.json')['nodes'];leaves=[r for r in nodes if r['final_leaf']];by={r['stable_membership_SHA256']:r for r in leaves};assert len(by)==len(leaves)==len(groups)==summary['n_final_leaf_clusters'];pc=load(P/'config/spec_v1.json');root_cfg=load(Path(pc['root_model'])/'config/spec_v1.json');fi={r['source']:r for r in root_cfg['frozen_preflight']['sources']};accepted=load(Path(pc['working_delivery'])/'outputs/accepted_TME_scoped_working_annotation_v1.json');w=load(accepted['producer_manifest']['path']);wb=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};features={};record_total=0;unfitted=0;source_totals={};canonical={'Tcell':'T','NKcell':'NK','gamma_delta_T_RNA_candidate':'T_gamma_delta_RNA_candidate'}
 for e in load(P/'outputs/cluster_hierarchy_index_v1.json'):
  source=e['source'];m=load(P/'outputs'/e['metadata']);assert sha(P/'outputs'/e['metadata'])==e['metadata_SHA256'];wr=wi[source];assert sha(wb/wr['arrays'])==wr['arrays_SHA256'];wm=load(wb/wr['metadata'])
  with np.load(wb/wr['arrays'],allow_pickle=False) as f:positions=np.flatnonzero(f['fit_pool_candidate_not_final_QC']);old=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[f['consensus_celltype_code'][positions]];n_all=len(f['native_row'])
  with np.load(m['arrays']['path'],allow_pickle=False) as f:np.testing.assert_array_equal(f['source_scope_position'],positions);codes=f['final_leaf_node']
  lookup={r['source_node_code']:r for r in leaves if r['source']==source};new=[lookup[int(i)]['cluster_RNA_proposal']['celltype'] for i in codes];observed=collections.Counter(zip(map(str,old),map(str,new)));reported={(r['baseline_v5_celltype'],r['new_cluster_RNA_proposal']):r['n_fit_records'] for r in cross if r['source']==source};assert len(reported)==sum(r['source']==source for r in cross) and dict(observed)==reported
  for r in cross:
   if r['source']==source:assert r['canonical_cluster_proposal']==canonical.get(r['new_cluster_RNA_proposal'],r['new_cluster_RNA_proposal'])
  counts=collections.Counter()
  for n in lookup.values():counts[n['cluster_RNA_proposal']['major']]+=n['n_records']
  ss=next(r for r in summary['source_summary'] if r['source']==source);assert ss['major_proposal_record_counts']==dict(counts) and ss['n_fit_records']==len(positions) and ss['n_unfitted_retained']==n_all-len(positions) and ss['n_leaf_groups']==len(lookup);source_totals[source]=counts;record_total+=len(positions);unfitted+=n_all-len(positions)
  ep=Path(fi[source]['candidate_feature_axis_path']);assert sha(ep)==fi[source]['candidate_feature_axis_SHA256'];f=load(ep);features[source]=set(f['genes']).intersection(m['native_measured_genes'])
 for r in groups:
  n=by[r['stable_membership_SHA256']];assert (r['source'],r['subcluster'],r['n_records'],r['major_proposal'],r['celltype_proposal'])==(n['source'],n['display_label'],n['n_records'],n['cluster_RNA_proposal']['major'],n['cluster_RNA_proposal']['celltype']);assert r['own_positive_marker_specificity']==(n['marker_evidence_type']=='own_subset_expression_specificity_not_inferential_DEG');assert r['state_program_names']==[k for k,v in n['descriptive_state_programs'].items() if v['descriptively_enriched']]
 excluded={g for v in GATES.values() for module in v for g in module}|{g for panel in STATE_PANELS.values() for g in panel};routes={'T_NK':{'T','NK','T_NK'},'B_Plasma':{'B','Plasma'},'Myeloid':{'Myeloid'},'Endothelial':{'Endothelial'},'Fibroblast':{'Fibroblast'},'Mural':{'Mural'},'Epithelial':{'Epithelial'},'Neural':{'Glial','Astrocyte','Neuron'},'Mast':{'Mast'},'Erythroid':{'Erythroid'},'Platelet':{'Platelet'},'unresolved':{'unresolved'}}
 for r in pre['routes']:
  counts={s:sum(n for major,n in v.items() if major in routes[r['route']]) for s,v in source_totals.items()};eligible=sorted(s for s,n in counts.items() if n>=rc['minimum_source_route_records']);common=set.intersection(*(features[s] for s in eligible)) if eligible else set();model=sorted(common-excluded);assert r['source_record_counts']==counts and r['eligible_source_ids']==eligible and r['n_sources_observed']==sum(n>0 for n in counts.values());assert r['n_common_unique_measured_candidate_features']==len(common) and r['model_feature_genes_private']==model and r['n_shared_model_features_after_identity_state_holdout']==len(model);assert r['passes_preliminary_shared_reference_size_and_measurement_gate']==(len(eligible)>=rc['minimum_fit_sources'] and len(model)>=rc['minimum_common_model_features'])
 assert record_total==summary['totals']['n_records']==rt['n_fit_records']==sum(r['n_records'] for r in groups);assert unfitted==summary['n_unfitted_retained'] and record_total+unfitted==summary['n_qualified_records'];assert len(source_totals)==summary['n_sources'];save(root/'outputs/full_cluster_review_scope_acceptance_v1.json',dict(all_sources_and_records_crosswalk_and_shared_measured_feature_gates_recomputed=True,n_sources=len(source_totals),n_fit_records=record_total,n_unfitted_retained=unfitted,n_leaf_groups=len(groups),not_identity_accuracy_or_shared_model_validation=True));save(root/'provenance/terminal_v1.json',dict(status='FULL_CLUSTER_REVIEW_DENOMINATORS_CROSSWALK_AND_FEATURE_GATES_ACCEPTED_NOT_SHARED_TAXONOMY',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),reviewer_terminal_SHA256=sha(R/'provenance/terminal_v1.json'),no_labels_clusters_counts_or_states_changed=True,new_CM_started=False))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
