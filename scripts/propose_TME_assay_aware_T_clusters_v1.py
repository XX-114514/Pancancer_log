"""Re-evaluate only source assays missing the entire TCR-chain module.

Cluster membership, expression markers, state context and display numbers are
kept. New calls are RNA proposals with an explicit alternate-anchor grade.
"""
import argparse,collections,datetime,hashlib,json,subprocess,sys,time
from pathlib import Path
import numpy as np
from TME_assay_aware_T_anchor_policy_v1 import CHAIN,ACCESSORY,SIGNAL,assay_aware_flags
from TME_cluster_annotation_policy_v1 import GATE_ORDER,annotate_cluster
from TME_native_measurement_core_v1 import canonical_native

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');H=Path(cfg['cluster_producer']);A=Path(cfg['cluster_auditor']);D=Path(cfg['cluster_delivery']);frozen=cfg['frozen_assay_spec']
 assert frozen['alternative_required_all']==['CD3D','CD3E'] and frozen['alternative_required_any']==ACCESSORY and frozen['signaling_gene_set']==SIGNAL and frozen['alternative_minimum_signaling_genes']==2
 for p in [H,A,D]:assert load(p/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(H/'provenance/terminal_v1.json');cat=H/'outputs/cluster_hierarchy_named_catalog_v1.json';assert sha(cat)==load(H/'provenance/terminal_v1.json')['catalog_SHA256'];catalog=load(cat)['nodes'];original={r['stable_membership_SHA256']:dict(r) for r in catalog};target=set(frozen['eligible_source_ids_private']);hc=load(H/'config/spec_v1.json');w=load(load(Path(hc['working_delivery'])/'outputs/accepted_TME_scoped_working_annotation_v1.json')['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};J=Path(cfg['joined']);jm=load(J/'outputs/joined_native_atlas_delivery_v1.json');sys.path.insert(0,str(J/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from native_panel_aggregate_v1 import gather_csr_rows
 from atlas_identity_release_keys_v1 import display_label
 t=subprocess.run([sys.executable,'-B',str(root/'tests/test_TME_assay_aware_T_anchor_policy_v1.py')],capture_output=True,text=True,timeout=60);(root/'provenance/test_stderr.txt').write_text(t.stderr);assert t.returncode==0 and 'Ran 6 tests' in t.stderr and 'skipped' not in t.stderr
 native=NativeArchiveCollection(jm['archive_manifest']);index=[];counts_index=[];seen=set();leaf_transitions=collections.Counter();totals=collections.Counter();all_expected=set()
 for entry in load(H/'outputs/cluster_hierarchy_index_v1.json'):
  source=entry['source'];mp=H/'outputs'/entry['metadata'];assert sha(mp)==entry['metadata_SHA256'];m=load(mp);eligible=not m['gate_assessed']['T'] and not any(g in m['native_measured_genes'] for g in CHAIN)
  if eligible:all_expected.add(source)
  if source not in target:index.append(dict(source=source,metadata=str(mp),metadata_SHA256=sha(mp)));continue
  assert eligible;seen.add(source);genes=m['native_measured_genes'];assert sha(m['arrays']['path'])==m['arrays']['SHA256']
  with np.load(m['arrays']['path'],allow_pickle=False) as f:pool=f['source_scope_position'];rows=f['native_row'];roots=f['root_cluster_code'];children=f['candidate_child_node']
  assert sha(m['statistics']['path'])==m['statistics']['SHA256']
  with np.load(m['statistics']['path'],allow_pickle=False) as f:cols=f['native_columns'];node_n=f['n']
  wr=wi[source];assert sha(base/wr['arrays'])==wr['arrays_SHA256']
  with np.load(base/wr['arrays'],allow_pickle=False) as f:rna=f['native_total_RNA'][pool];np.testing.assert_array_equal(rows,f['native_row'][pool])
  csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']];support=np.zeros((m['n_candidate_nodes'],len(GATE_ORDER)),np.int64);assessed=None;mode=None
  for lo in range(0,len(pool),2048):
   take=np.arange(lo,min(lo+2048,len(pool)));raw=canonical_native(gather_csr_rows(csr,rows[take],max_span=8192));np.testing.assert_array_equal(np.asarray(raw.sum(axis=1)).ravel(),rna[take]);flags,a,b=assay_aware_flags(raw[:,cols].tocsr(),genes)
   if assessed is None:assessed=a;mode=b
   else:assert assessed==a and mode==b
   for mapping in [roots[take],children[take]]:
    valid=mapping>=0
    for j in range(len(GATE_ORDER)):support[:,j]+=np.bincount(mapping[valid],weights=flags[valid,j],minlength=len(support)).astype(np.int64)
  assert mode=='T_chain_unmeasured_alternate_CD3_accessory_two_signal_RNA_evidence';current=[];changed=0
  for old in m['nodes']:
   r=dict(original[old['stable_membership_SHA256']]);node=r['source_node_code'];n=r['n_records'];assert n==node_n[node];call=annotate_cluster(n,dict(zip(GATE_ORDER,support[node].tolist())),assessed,hc['minimum_named_group_cells'],hc['family_fraction']);previous=r['cluster_RNA_proposal'];r.update(previous_cluster_RNA_proposal=previous,cluster_RNA_proposal=call,assay_aware_T_evidence_mode=mode,chain_genes_not_inferred=True,display_label=display_label(call['celltype'],r['marker'],r['global_display_number']));current.append(r)
   if r['final_leaf']:
    leaf_transitions[previous['celltype'],call['celltype']]+=n
    if previous!=call:changed+=n
  updated={r['source_node_code']:r for r in current};meta=dict(m,previous_cluster_source_metadata=dict(path=str(mp),SHA256=sha(mp)),nodes=current,previous_gate_assessed=m['gate_assessed'],gate_assessed=assessed,original_TCR_chain_module_assessed=False,assay_aware_T_evidence_mode=mode)
  for p in meta['parents']:p['previous_cluster_RNA_proposal']=p['cluster_RNA_proposal'];p['cluster_RNA_proposal']=updated[p['source_root_cluster']]['cluster_RNA_proposal']
  out=root/'outputs'/(source+'_assay_aware_cluster_proposals_v1.json');save(out,meta);cp=root/'outputs'/(source+'_assay_aware_cluster_support_v1.npz')
  with cp.open('xb') as f:np.savez_compressed(f,node_gate_counts=support,node_counts=node_n)
  counts_index.append(dict(source=source,path=str(cp),SHA256=sha(cp),metadata=out.name,metadata_SHA256=sha(out),n_model_pool_records_directly_read=len(pool)));index.append(dict(source=source,metadata=out.name,metadata_SHA256=sha(out)));totals.update(n_sources_reassessed=1,n_pool_records_directly_read=len(pool),n_final_leaf_records_with_changed_proposals=changed);replacement={r['stable_membership_SHA256']:r for r in current};catalog=[replacement.get(r['stable_membership_SHA256'],r) for r in catalog];native._fields.clear()
 assert seen==target==all_expected
 ip=root/'outputs/cluster_hierarchy_index_v1.json';cp=root/'outputs/cluster_hierarchy_named_catalog_v1.json';ci=root/'outputs/assay_aware_cluster_support_index_v1.json';save(ip,index);save(cp,dict(nodes=catalog,cluster_membership_markers_numbers_and_states_unchanged=True));save(ci,counts_index);summary=load(H/'outputs/cluster_hierarchy_summary_v1.json');summary.update(root_cluster_proposal_counts=dict(collections.Counter(r['cluster_RNA_proposal']['celltype'] for r in catalog if r['level']==1)),final_leaf_cluster_proposal_counts=dict(collections.Counter(r['cluster_RNA_proposal']['celltype'] for r in catalog if r['final_leaf'])),final_leaf_record_proposal_counts=dict(sum((collections.Counter({r['cluster_RNA_proposal']['celltype']:r['n_records']}) for r in catalog if r['final_leaf']),collections.Counter())),assay_aware_anchor_scope=dict(totals),assay_transition_counts=[dict(previous=a,current=b,n_records=n) for (a,b),n in sorted(leaf_transitions.items())]);save(root/'outputs/cluster_hierarchy_summary_v1.json',summary)
 save(root/'provenance/terminal_v1.json',dict(status='ASSAY_AWARE_NATIVE_T_CLUSTER_PROPOSALS_REQUIRE_INDEPENDENT_REPLAY',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),index_SHA256=sha(ip),catalog_SHA256=sha(cp),support_index_SHA256=sha(ci),original_cluster_catalog_SHA256=sha(cat),working_v5_counts_states_and_cluster_memberships_unchanged=True,not_orthogonal_T_identity_validation=True,no_gamma_chain_evidence_inferred=True,new_CM_started=False,wall_seconds=time.monotonic()-start))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
