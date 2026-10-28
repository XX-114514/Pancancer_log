"""Replay alternate T support with an independent CSR gather and Boolean rule.

Reuse the frozen cluster classifier explicitly as policy replay, never as
orthogonal identity accuracy. All membership, marker and state fields are exact
invariants; source assays with measured chains must remain untouched.
"""
import argparse,collections,datetime,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from TME_cluster_annotation_policy_v1 import GATES,GATE_ORDER,annotate_cluster

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')

def independently_replay_alternate_flags(raw,genes):
 by=dict(zip(genes,range(len(genes))));n=raw.shape[0]
 assert not set(['TRAC','TRBC1','TRBC2','TRDC'])&set(by)
 def any_of(symbols):
  present=[by[g] for g in symbols if g in by]
  return np.asarray((raw[:,present]>0).getnnz(axis=1)).ravel()>0 if present else np.zeros(n,bool)
 assert {'CD3D','CD3E'}.issubset(by)
 signals=np.zeros(n,np.int16)
 for g in ['CD247','LCK','LAT','ZAP70','BCL11B']:signals+=any_of([g])
 flags={};assessed={}
 for name,modules in GATES.items():
  assessed[name]=all(any(g in by for g in m) for m in modules)
  v=np.ones(n,bool)
  for module in modules:v&=any_of(module)
  flags[name]=v if assessed[name] else np.zeros(n,bool)
 flags['T']=any_of(['CD3D'])&any_of(['CD3E'])&any_of(['CD2','CD7','TRAT1'])&(signals>=2)
 assessed['T']=True
 for name in ['Treg','CD4_T','CD8_T','gamma_delta_T_RNA_candidate']:flags[name]&=flags['T']
 flags['CD4_T']&=~flags['Treg']
 assert not flags['gamma_delta_T_RNA_candidate'].any() and not assessed['gamma_delta_T_RNA_candidate']
 return np.column_stack([flags[k] for k in GATE_ORDER]),assessed

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);H=Path(cfg['original_cluster_producer']);J=Path(cfg['joined']);EA=Path(cfg['independent_native_auditor'])
 while not (P/'provenance/supervisor_terminal_v1.json').exists():time.sleep(3)
 assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 pt=load(P/'provenance/terminal_v1.json');pc=load(P/'config/spec_v1.json');hc=load(H/'config/spec_v1.json');target=set(pc['frozen_assay_spec']['eligible_source_ids_private'])
 assert sha(P/'outputs/cluster_hierarchy_named_catalog_v1.json')==pt['catalog_SHA256'] and sha(P/'outputs/cluster_hierarchy_index_v1.json')==pt['index_SHA256']
 assert sha(P/'outputs/assay_aware_cluster_support_index_v1.json')==pt['support_index_SHA256']
 oldcat=load(H/'outputs/cluster_hierarchy_named_catalog_v1.json')['nodes'];newcat=load(P/'outputs/cluster_hierarchy_named_catalog_v1.json')['nodes'];oldby={r['stable_membership_SHA256']:r for r in oldcat};newby={r['stable_membership_SHA256']:r for r in newcat};assert list(oldby)==list(newby) and len(oldby)==len(oldcat)==len(newcat)
 oldindex={r['source']:r for r in load(H/'outputs/cluster_hierarchy_index_v1.json')};newindex={r['source']:r for r in load(P/'outputs/cluster_hierarchy_index_v1.json')};supportindex={r['source']:r for r in load(P/'outputs/assay_aware_cluster_support_index_v1.json')};assert set(oldindex)==set(newindex) and set(supportindex)==target
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(EA/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from audit_full_TME_native_evidence_v1 import gather
 jm=load(J/'outputs/joined_native_atlas_delivery_v1.json');native=NativeArchiveCollection(jm['archive_manifest']);w=load(load(Path(hc['working_delivery'])/'outputs/accepted_TME_scoped_working_annotation_v1.json')['producer_manifest']['path']);wb=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};mode='T_chain_unmeasured_alternate_CD3_accessory_two_signal_RNA_evidence';totals=collections.Counter();results=[]
 for source in oldindex:
  oi=oldindex[source];ni=newindex[source];op=H/'outputs'/oi['metadata'];npth=P/'outputs'/ni['metadata'];assert sha(op)==oi['metadata_SHA256'] and sha(npth)==ni['metadata_SHA256'];old=load(op);new=load(npth)
  nodes=[r for r in newcat if r['source']==source]
  if source not in target:
   assert op==npth and old==new
   for r in nodes:assert r==oldby[r['stable_membership_SHA256']]
   totals.update(n_untouched_sources=1);continue
  assert new['previous_cluster_source_metadata']==dict(path=str(op),SHA256=sha(op))
  for key,value in old.items():
   if key not in ['nodes','parents','gate_assessed']:assert new[key]==value,key
  assert new['previous_gate_assessed']==old['gate_assessed'] and new['original_TCR_chain_module_assessed'] is False
  assert sha(old['arrays']['path'])==old['arrays']['SHA256'] and sha(old['statistics']['path'])==old['statistics']['SHA256']
  with np.load(old['arrays']['path'],allow_pickle=False) as f:pool=f['source_scope_position'];rows=f['native_row'];roots=f['root_cluster_code'];children=f['candidate_child_node']
  with np.load(old['statistics']['path'],allow_pickle=False) as f:cols=f['native_columns'];node_n=f['n']
  wr=wi[source];assert sha(wb/wr['arrays'])==wr['arrays_SHA256']
  with np.load(wb/wr['arrays'],allow_pickle=False) as f:rna=f['native_total_RNA'][pool];np.testing.assert_array_equal(rows,f['native_row'][pool])
  csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']];support=np.zeros((len(node_n),len(GATE_ORDER)),np.int64);assessed=None
  for lo in range(0,len(rows),1024):
   take=np.arange(lo,min(lo+1024,len(rows)));raw=gather(csr,rows[take]);np.testing.assert_array_equal(np.asarray(raw.sum(axis=1)).ravel(),rna[take]);flags,a=independently_replay_alternate_flags(raw[:,cols].tocsr(),old['native_measured_genes'])
   if assessed is None:assessed=a
   else:assert assessed==a
   for mapping in [roots[take],children[take]]:
    for j in range(len(GATE_ORDER)):
     ok=(mapping>=0)&flags[:,j];support[:,j]+=np.bincount(mapping[ok],minlength=len(support))
  si=supportindex[source];assert sha(si['path'])==si['SHA256'] and si['metadata_SHA256']==sha(npth)
  with np.load(si['path'],allow_pickle=False) as f:np.testing.assert_array_equal(f['node_gate_counts'],support);np.testing.assert_array_equal(f['node_counts'],node_n)
  assert new['gate_assessed']==assessed and new['assay_aware_T_evidence_mode']==mode
  changed=0
  for r in nodes:
   before=oldby[r['stable_membership_SHA256']];node=r['source_node_code'];expected=annotate_cluster(r['n_records'],dict(zip(GATE_ORDER,support[node].tolist())),assessed,hc['minimum_named_group_cells'],hc['family_fraction']);assert r['cluster_RNA_proposal']==expected and r['previous_cluster_RNA_proposal']==before['cluster_RNA_proposal'] and r['assay_aware_T_evidence_mode']==mode and r['chain_genes_not_inferred']
   for key,value in before.items():
    if key not in ['cluster_RNA_proposal','display_label']:assert r[key]==value,key
   expected_label=expected['celltype']+'_'+(r['marker'] if r['marker'] else 'marker_unresolved')+'_'+str(r['global_display_number'])
   assert r['display_label']==expected_label
   if r['final_leaf'] and before['cluster_RNA_proposal']!=expected:changed+=r['n_records']
  assert len(old['parents'])==len(new['parents'])
  newnodes={r['source_node_code']:r for r in new['nodes']}
  for before,after in zip(old['parents'],new['parents']):
   for key,value in before.items():
    if key!='cluster_RNA_proposal':assert after[key]==value
   assert after['previous_cluster_RNA_proposal']==before['cluster_RNA_proposal'] and after['cluster_RNA_proposal']==newnodes[after['source_root_cluster']]['cluster_RNA_proposal']
  # Metadata nodes and the named catalog must refer to the same proposals.
  for r in new['nodes']:assert r==newby[r['stable_membership_SHA256']]
  totals.update(n_sources_reassessed=1,n_pool_records_directly_read=len(rows),n_final_leaf_records_with_changed_proposals=changed);results.append(dict(source=source,n_records_independently_read=len(rows),all_group_support_counts_exact=True,membership_markers_numbers_states_and_measured_masks_unchanged=True,no_gamma_chains_inferred=True));native._fields.clear()
 assert {k:totals[k] for k in pt['totals']}==pt['totals']
 out=root/'outputs/assay_aware_T_cluster_full_replay_v1.json';save(out,dict(totals=dict(totals),sources=results,original_classifier_reused_as_explicit_policy_replay_not_accuracy=True));save(root/'provenance/terminal_v1.json',dict(status='ASSAY_AWARE_T_CLUSTER_FULL_NATIVE_SUPPORT_REPLAY_ACCEPTED_PROVISIONAL_RNA_VIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),output_SHA256=sha(out),totals=dict(totals),independent_CSR_and_alternate_boolean_implementation=True,all_existing_cluster_membership_marker_state_and_number_fields_preserved=True,not_orthogonal_identity_validation=True,new_CM_started=False,wall_seconds=time.monotonic()-start))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
