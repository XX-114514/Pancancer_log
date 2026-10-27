"""Exact full-gene profiles for TCR candidates within unresolved parents.

Directly reread candidate native CSR rows. Derive each nonempty remainder from
its previously independently accepted exact parent, retaining missing features.
These are descriptive review profiles and do not assign new cell identities.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
from TME_native_measurement_core_v1 import canonical_native,update_group_stats,expression_candidates

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');R=Path(cfg['context_producer']);A=Path(cfg['context_auditor'])
 for p in [R,A]:assert load(p/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(R/'provenance/terminal_v1.json')
 rc=load(R/'config/spec_v1.json');K=Path(rc['TCR_producer']);pc=load(K/'config/spec_v1.json');E=Path(pc['native_evidence']);J=Path(pc['joined']);EA=Path(pc['native_audit']);assert load(EA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(E/'provenance/terminal_v1.json')
 accepted=load(rc['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];w=load(accepted['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};catalog=load(w['group_catalog']['path'])['groups']
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(root/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
 from native_panel_aggregate_v1 import gather_csr_rows
 from atlas_identity_release_keys_v1 import membership_key
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);ei={r['source']:r for r in load(E/'outputs/TME_full_native_evidence_index_v1.json')};ip=R/'outputs/unknown_parent_TCR_context_index_v1.json';assert sha(ip)==load(R/'provenance/terminal_v1.json')['index_SHA256'];totals=collections.Counter();outputs=[];resource_cfg=dict(pc,output_cap_bytes=cfg.get('output_cap_bytes',2147483648))
 for entry in load(ip):
  source=entry['source'];rp=R/'outputs'/entry['metadata'];assert sha(rp)==entry['metadata_SHA256'];review=load(rp)
  if not review['n_candidates']:continue
  storage_gate(root,resource_cfg);assert sha(review['arrays'])==review['arrays_SHA256']
  with np.load(review['arrays'],allow_pickle=False) as f:target=f['target_TCR_unknown_parent'];rows=f['native_row']
  wr=wi[source];wm=load(base/wr['metadata']);wa=base/wr['arrays'];assert sha(wa)==wr['arrays_SHA256']
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','annotation_group_code','native_total_RNA']}
  np.testing.assert_array_equal(rows,q['native_row']);er=ei[source];ep=E/'outputs'/er['metadata'];assert sha(ep)==er['metadata_SHA256'];e=load(ep);cols=np.asarray(e['native_columns']);genes=e['native_measured_gene_names'];esp=e['group_statistics'];assert sha(esp['path'])==esp['SHA256']
  with np.load(esp['path'],allow_pickle=False) as f:old={k:f[k] for k in ['n','cp_sum','log_sum','detected','native_columns']}
  np.testing.assert_array_equal(old['native_columns'],cols);oldby={g['stable_group_sha256']:g for g in e['groups']};current={g['source_annotation_group_code']:g for g in catalog if g['source']==source};parents=sorted(set(map(int,q['annotation_group_code'][target])));ids=native.native_ids(source);assert text_axis_sha256(ids)==wm['native_ids_source_SHA256']
  for parent in parents:
   g=current[parent];assert g['stable_group_sha256'] in oldby;sel=np.flatnonzero(q['annotation_group_code']==parent);assert len(sel)==old['n'][oldby[g['stable_group_sha256']]['source_annotation_group_code']]==g['n_scope_group_rows'];assert membership_key(source,ids[q['source_native_position'][sel]].astype(str).tolist())==g['stable_group_sha256']
  stat={'n':np.zeros(len(parents),np.int64),**{k:np.zeros((len(parents),len(cols)),np.float64) for k in ['cp_sum','log_sum','detected']}};take=np.flatnonzero(target);codes=np.searchsorted(parents,q['annotation_group_code'][take]);csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']]
  for lo in range(0,len(take),2048):
   use=take[lo:lo+2048];raw=canonical_native(gather_csr_rows(csr,q['native_row'][use],max_span=8192));total=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(total,q['native_total_RNA'][use]);update_group_stats(stat,raw[:,cols].tocsr(),total,codes[lo:lo+len(use)])
  groups=[];vectors={k:[] for k in stat};scope_codes=np.full(len(target),-1,np.int32);units=native.frozen_expression_units(source)[q['source_native_position']].astype(str)
  for j,parent in enumerate(parents):
   g=current[parent];oldcode=oldby[g['stable_group_sha256']]['source_annotation_group_code'];remainder={k:old[k][oldcode]-stat[k][j] for k in stat}
   for name,value in remainder.items():
    assert np.min(value)>=(0 if name in ['n','detected'] else -1e-7)
    if name in ['cp_sum','log_sum']:remainder[name]=np.maximum(value,0)
   for role,mask,moments,other in [('TCR_candidate_not_reclassified',(q['annotation_group_code']==parent)&target,{k:stat[k][j] for k in stat},remainder),('unresolved_remainder',(q['annotation_group_code']==parent)&~target,remainder,{k:stat[k][j] for k in stat})]:
    sel=np.flatnonzero(mask);assert len(sel)==moments['n']
    if not len(sel):continue
    code=len(groups);scope_codes[sel]=code;key=membership_key(source,ids[q['source_native_position'][sel]].astype(str).tolist());count=int(other['n']);mean=moments['cp_sum']/len(sel);det=moments['detected']/len(sel);candidates=expression_candidates(genes,mean,det,other['cp_sum']/count if count else None,other['detected']/count if count else None,.1,20)
    missing=np.isin(np.char.lower(np.char.strip(units[sel])),['','nan','none','na','n/a','unknown','__missing__','__source_value_missing__'])
    groups.append(dict(measurement_subset_code=code,previous_group_code=parent,previous_group_number=g['global_subcluster_number'],previous_group_SHA256=g['stable_group_sha256'],role=role,n_exact_subset=len(sel),stable_membership_SHA256=key,native_rows_SHA256=hashlib.sha256(np.asarray(q['native_row'][sel],'<i8').tobytes()).hexdigest(),native_parent_profile_row=oldcode,n_declared_expression_units_not_libraries_or_patients=len(set(units[sel][~missing])),n_records_with_missing_declared_unit=int(missing.sum()),own_full_native_gene_marker_candidates_not_DEG=candidates,comparator_scope='other_role_in_same_exact_current_unresolved_parent' if count else 'no_nonempty_parent_remainder_or_counterpart',n_comparator_cells=count,new_identity_assigned=False))
    for name,value in moments.items():vectors[name].append(value)
  arrays={name:np.asarray(value) for name,value in vectors.items()};npz=root/'outputs'/(source+'_unknown_parent_TCR_profiles_v1.npz')
  with npz.open('xb') as f:np.savez_compressed(f,**arrays,native_columns=cols,scope_position=np.flatnonzero(scope_codes>=0),native_row=q['native_row'][scope_codes>=0],measurement_subset_code=scope_codes[scope_codes>=0])
  with np.load(npz,allow_pickle=False) as f:
   for name,value in arrays.items():np.testing.assert_array_equal(f[name],value)
  for parent in parents:
   ix=[g['measurement_subset_code'] for g in groups if g['previous_group_code']==parent];oldcode=oldby[current[parent]['stable_group_sha256']]['source_annotation_group_code']
   for name in ['n','detected']:np.testing.assert_array_equal(arrays[name][ix].sum(axis=0),old[name][oldcode])
   for name in ['cp_sum','log_sum']:np.testing.assert_allclose(arrays[name][ix].sum(axis=0),old[name][oldcode],rtol=2e-10,atol=1e-7)
  mp=root/'outputs'/(source+'_unknown_parent_TCR_profiles_v1.json');out=dict(source=source,n_parent_groups=len(parents),n_exact_subsets=len(groups),n_candidate_native_rows_read=len(take),n_parent_scope_records=int((scope_codes>=0).sum()),groups=groups,native_measured_genes=genes,statistics=dict(path=str(npz),SHA256=sha(npz)),context_metadata_SHA256=sha(rp),working_arrays_SHA256=sha(wa),accepted_complete_native_metadata_SHA256=sha(ep),candidate_native_rows_directly_read=True,remainder_from_previously_accepted_exact_parent=True,all_original_working_labels_states_and_source_gate_preserved=True);save(mp,out);outputs.append(dict(source=source,metadata=mp.name,metadata_SHA256=sha(mp)));totals.update(n_sources=1,n_parent_groups=len(parents),n_exact_subsets=len(groups),n_candidate_native_rows_read=len(take),n_parent_scope_records=out['n_parent_scope_records']);native._fields.clear();del q,old,stat,arrays;gc.collect()
 ip=root/'outputs/unknown_parent_TCR_profile_index_v1.json';save(ip,outputs);assert totals['n_candidate_native_rows_read']==load(R/'provenance/terminal_v1.json')['totals']['n_candidate_records']
 save(root/'provenance/terminal_v1.json',dict(status='ALL_UNKNOWN_PARENT_TCR_FULL_NATIVE_PROFILES_REQUIRE_INDEPENDENT_AUDIT',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),index_SHA256=sha(ip),context_terminal_SHA256=sha(R/'provenance/terminal_v1.json'),context_independent_audit_SHA256=sha(A/'provenance/terminal_v1.json'),new_identity_assigned=False,new_CM_started=False,wall_seconds=time.monotonic()-start,peak_RSS_bytes_Linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
