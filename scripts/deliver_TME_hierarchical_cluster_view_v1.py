"""After real producer/auditor success, exercise view and export group-only tables."""
import argparse,collections,csv,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def table(p,rows):
 with Path(p).open('x',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def main(root):
 begin=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);A=Path(cfg['auditor'])
 for r in [P,A]:
  while not (r/'provenance/supervisor_terminal_v1.json').exists():time.sleep(3)
  assert load(r/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 at=A/'provenance/terminal_v1.json';assert load(at)['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json');pc=load(P/'config/spec_v1.json');B=Path(pc['working_delivery'])/'outputs/accepted_TME_scoped_working_annotation_v1.json';bm=load(load(B)['producer_manifest']['path']);J=Path(cfg['joined']);sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(root/'scripts'))
 catalog=P/'outputs/cluster_hierarchy_named_catalog_v1.json';index=P/'outputs/cluster_hierarchy_index_v1.json';assert sha(catalog)==load(P/'provenance/terminal_v1.json')['catalog_SHA256'];assert sha(index)==load(P/'provenance/terminal_v1.json')['index_SHA256'];accept=dict(schema='accepted_numeric_TME_hierarchical_cluster_view_v1',accepted_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),catalog=dict(path=str(catalog),SHA256=sha(catalog)),index=dict(path=str(index),SHA256=sha(index)),audit_terminal=dict(path=str(at),SHA256=sha(at)),baseline_working_acceptance=dict(path=str(B),SHA256=sha(B)),scope='first_round_source_local_T_NK_and_frozen_unknown_controls_not_all_TME',all_baseline_fields_and_counts_preserved=True,biological_accuracy_certified=False,final_shared_taxonomy=False,new_CM_started=False)
 ap=root/'outputs/accepted_TME_hierarchical_cluster_view_v1.json';save(ap,accept)
 from TME_hierarchical_cluster_view_v1 import TMEHierarchicalClusterView
 view=TMEHierarchicalClusterView(ap);nodes=load(catalog)['nodes'];roles={}
 for source in view.sources:
  s=view.source(source);a=s['arrays'];base=view.baseline.source(source)['arrays'];pool=a['source_scope_position'];leaf=a['final_leaf_node']
  for label in ['Treg','CD8_T','CD4_T','gamma_delta_T_RNA_candidate','Tcell','unresolved']:
   key='cluster_'+label
   if key in roles:continue
   for n in np.unique(leaf):
    r=s['nodes'][int(n)]
    if r['cluster_RNA_proposal']['celltype']==label:roles[key]=(source,np.array([pool[np.flatnonzero(leaf==n)[0]]],np.int64));break
  child=np.array([s['nodes'][int(n)]['level']==2 for n in leaf])
  if 'nested_child_and_repeated_query' not in roles and child.any():
   p=int(pool[np.flatnonzero(child)[0]]);roles['nested_child_and_repeated_query']=(source,np.array([p,p],np.int64))
  outside=np.setdiff1d(np.arange(len(base['native_row'])),pool)
  if 'outside_first_round_pool' not in roles and len(outside):roles['outside_first_round_pool']=(source,outside[:1])
  nofit=np.flatnonzero(~base['fit_pool_candidate_not_final_QC'])
  if 'unfitted_retained_baseline' not in roles and len(nofit):roles['unfitted_retained_baseline']=(source,nofit[:1])
 assert {'nested_child_and_repeated_query','outside_first_round_pool','unfitted_retained_baseline','cluster_Tcell','cluster_unresolved'}.issubset(roles)
 tests=[]
 for role,(source,p) in roles.items():
  actual=view.metadata(source,p)['columns'];baseline=view.baseline.metadata(source,p)['columns']
  for k,v in baseline.items():np.testing.assert_equal(actual[k],v)
  if role in ['outside_first_round_pool','unfitted_retained_baseline']:assert (actual['cluster_leaf_number']==-1).all() and (actual['cluster_view_status']=='outside_frozen_first_round_pool').all()
  else:assert (actual['cluster_leaf_number']>=0).all() and all(x is not None for x in actual['cluster_leaf_display'])
  genes=['CD3D','FOXP3','__UNMEASURED_SENTINEL__'];x=view.read_counts(source,p,genes);y=view.baseline.read_counts(source,p,genes);assert (x['counts']!=y['counts']).nnz==0
  for k in ['native_row','native_total_counts','measured']:np.testing.assert_array_equal(x[k],y[k])
  assert not x['measured'][-1];tests.append(dict(role=role,n_query_records=len(p),baseline_fields_and_native_counts_mask_equal=True))
 nodeby={(r['source'],r['source_node_code']):r for r in nodes};rows=[];links=[]
 for r in nodes:
  call=r['cluster_RNA_proposal'];parent=nodeby[(r['source'],r['source_root_cluster'])];states=[k for k,v in r['descriptive_state_programs'].items() if v['descriptively_enriched']]
  rows.append(dict(source=r['source'],level=r['level'],root_cluster_number=r['source_root_cluster'],child_cluster_number=r['child_cluster_number'],display_number=r['global_display_number'],subcluster=r['display_label'],major=call['major'],celltype=call['celltype'],marker=r['marker'],n_records=r['n_records'],final_leaf=r['final_leaf'],annotation_status=call['status'],marker_evidence_type=r['marker_evidence_type'],marker_candidate_names=';'.join(x['gene'] for x in r['marker_candidates']),state_program_names=';'.join(states),stable_membership_SHA256=r['stable_membership_SHA256'],parent_membership_SHA256=parent['stable_membership_SHA256'],n_declared_units_not_verified_patients=r['n_declared_expression_units_not_verified_patients_or_libraries'],parent_marker_borrowed=False))
  if r['level']==2:links.append(dict(source=r['source'],parent_display=parent['display_label'],child_display=r['display_label'],child_n_records=r['n_records'],child_is_final_leaf=r['final_leaf'],parent_membership_SHA256=parent['stable_membership_SHA256'],child_membership_SHA256=r['stable_membership_SHA256']))
 table(root/'outputs/TME_hierarchical_cluster_catalog_names_v1.tsv',rows);table(root/'outputs/TME_hierarchical_cluster_leaf_names_v1.tsv',[r for r in rows if r['final_leaf']]);table(root/'outputs/TME_hierarchical_cluster_parent_child_v1.tsv',links)
 summary=load(P/'outputs/cluster_hierarchy_summary_v1.json');save(root/'outputs/TME_hierarchical_cluster_summary_v1.json',summary);coverage=[]
 for source in view.sources:
  s=view.source(source);q=view.baseline.source(source)['arrays'];coverage.append(dict(source=source,n_current_TME_records=len(q['native_row']),n_first_round_cluster_pool=len(s['arrays']['native_row']),n_outside_first_round_pool_retained=len(q['native_row'])-len(s['arrays']['native_row']),n_root_clusters=s['metadata']['n_roots'],n_final_leaf_clusters=s['metadata']['n_final_leaves']))
 table(root/'outputs/TME_hierarchical_cluster_source_coverage_v1.tsv',coverage);assert sum(r['n_first_round_cluster_pool'] for r in coverage)==summary['totals']['n_pool_records'];assert sum(r['n_current_TME_records'] for r in coverage)==load(bm['summary']['path'])['totals']['n_source_qualified_rows'];save(root/'outputs/TME_hierarchical_cluster_view_readback_v1.json',dict(tested_cases=tests,numeric_membership_and_native_marker_audit_accepted=True,baseline_v5_not_overwritten=True,all_current_TME_sources_readable_with_explicit_first_round_scope=True,individual_records_and_expression_values_not_in_exported_tables=True,new_CM_started=False))
 files={p.name:sha(p) for p in (root/'outputs').glob('*') if p.is_file()};save(root/'provenance/terminal_v1.json',dict(status='NUMERICALLY_AUDITED_CLUSTER_TREE_VIEW_AND_SAFE_GROUP_NAME_TABLES_DELIVERED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=files,totals=summary['totals'],n_real_reader_cases=len(tests),baseline_current_v5_unchanged=True,not_final_common_taxonomy=True,new_CM_started=False,wall_seconds=time.monotonic()-begin))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
