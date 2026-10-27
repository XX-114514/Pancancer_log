"""Exact native evidence for accepted common-context identity deltas.

Directly read only proposed records; prove the unchanged-label remainder from
its accepted exact parent. Compare descriptive markers with the accepted full
source common-space pool. Current identities and all archived data stay frozen.
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
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['proposal_producer']);A=Path(cfg['proposal_auditor'])
 for s in [P,A]:assert load(s/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json');pc=load(P/'config/spec_v1.json');S=Path(pc['group_profile_producer']);sc=load(S/'config/spec_v1.json');M=Path(sc['model_producer']);mc=load(M/'config/spec_v1.json');R=Path(mc['context_producer']);rc=load(R/'config/spec_v1.json');kc=load(Path(rc['TCR_producer'])/'config/spec_v1.json');E=Path(kc['native_evidence']);EA=Path(kc['native_audit']);J=Path(kc['joined']);assert load(EA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(E/'provenance/terminal_v1.json')
 accepted=load(rc['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];w=load(accepted['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};catalog=load(w['group_catalog']['path'])['groups'];ei={r['source']:r for r in load(E/'outputs/TME_full_native_evidence_index_v1.json')};si={r['source']:r for r in load(S/'outputs/common_TNK_group_native_profile_index_v1.json')}
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(root/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
 from native_panel_aggregate_v1 import gather_csr_rows
 from atlas_identity_release_keys_v1 import membership_key
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);ip=P/'outputs/common_reentry_identity_proposal_index_v1.json';assert sha(ip)==load(P/'provenance/terminal_v1.json')['index_SHA256'];outputs=[];totals=collections.Counter();resource_cfg=dict(kc,output_cap_bytes=cfg.get('output_cap_bytes',2147483648))
 for entry in load(ip):
  source=entry['source'];pp=P/'outputs'/entry['metadata'];assert sha(pp)==entry['metadata_SHA256'];p=load(pp)
  if not p['n_proposed']:continue
  storage_gate(root,resource_cfg);assert sha(p['arrays'])==p['arrays_SHA256']
  with np.load(p['arrays'],allow_pickle=False) as f:proposed=f['proposed_T'];typecode=f['proposed_celltype_code'];rows=f['native_row']
  wr=wi[source];wm=load(base/wr['metadata']);wa=base/wr['arrays'];assert sha(wa)==wr['arrays_SHA256']
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','annotation_group_code','consensus_celltype_code','consensus_major_code','native_total_RNA']}
  np.testing.assert_array_equal(rows,q['native_row']);before=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']];major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];assert (before[proposed]=='unresolved').all() and (major[proposed]=='unresolved').all();after=before.copy();after[proposed]=np.asarray(p['celltype_dictionary'],object)[typecode[proposed]]
  ep=E/'outputs'/ei[source]['metadata'];assert sha(ep)==ei[source]['metadata_SHA256'];e=load(ep);cols=np.asarray(e['native_columns']);genes=e['native_measured_gene_names'];assert sha(e['group_statistics']['path'])==e['group_statistics']['SHA256']
  with np.load(e['group_statistics']['path'],allow_pickle=False) as f:old={k:f[k] for k in ['n','cp_sum','log_sum','detected']}
  oldby={g['stable_group_sha256']:g for g in e['groups']};current={g['source_annotation_group_code']:g for g in catalog if g['source']==source};parents=sorted(set(map(int,q['annotation_group_code'][proposed])));groups=[];code_for_scope=np.full(len(before),-1,np.int32);ids=native.native_ids(source);assert text_axis_sha256(ids)==wm['native_ids_source_SHA256']
  for parent in parents:
   oldg=current[parent];assert oldg['stable_group_sha256'] in oldby;domain=q['annotation_group_code']==parent;assert int(domain.sum())==oldg['n_scope_group_rows']==old['n'][oldby[oldg['stable_group_sha256']]['source_annotation_group_code']];assert membership_key(source,ids[q['source_native_position'][domain]].astype(str).tolist())==oldg['stable_group_sha256']
   for label in sorted(set(after[domain])):
    sel=np.flatnonzero(domain&(after==label));j=len(groups);code_for_scope[sel]=j;direct=bool(proposed[sel].all());assert direct or not proposed[sel].any();assert direct or label=='unresolved'
    groups.append(dict(measurement_subset_code=j,previous_group_code=parent,previous_group_number=oldg['global_subcluster_number'],previous_group_SHA256=oldg['stable_group_sha256'],fine_celltype=str(label),consensus_major_proposal='T' if direct else 'unresolved',n_exact_subset=len(sel),exact_membership_reused=len(sel)==int(domain.sum()),native_profile_derivation='direct_proposed_original_CSR' if direct else 'accepted_exact_parent_minus_all_new_proposals',scope_positions=sel))
  ng=len(groups);stat={'n':np.zeros(ng,np.int64),**{k:np.zeros((ng,len(cols)),np.float64) for k in ['cp_sum','log_sum','detected']}};take=np.flatnonzero(proposed);csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']]
  for lo in range(0,len(take),2048):
   sel=take[lo:lo+2048];raw=canonical_native(gather_csr_rows(csr,q['native_row'][sel],max_span=8192));rna=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(rna,q['native_total_RNA'][sel]);update_group_stats(stat,raw[:,cols].tocsr(),rna,code_for_scope[sel])
  for parent in parents:
   ix=[g['measurement_subset_code'] for g in groups if g['previous_group_code']==parent];changed=[j for j in ix if groups[j]['consensus_major_proposal']=='T'];rem=[j for j in ix if groups[j]['consensus_major_proposal']=='unresolved'];oldcode=oldby[current[parent]['stable_group_sha256']]['source_annotation_group_code'];assert len(rem)<=1
   if rem:
    for name in stat:
     value=old[name][oldcode]-stat[name][changed].sum(axis=0);assert np.min(value)>=(0 if name in ['n','detected'] else -1e-7);stat[name][rem[0]]=np.maximum(value,0) if name in ['cp_sum','log_sum'] else value
   for name in ['n','detected']:np.testing.assert_array_equal(stat[name][ix].sum(axis=0),old[name][oldcode])
   for name in ['cp_sum','log_sum']:np.testing.assert_allclose(stat[name][ix].sum(axis=0),old[name][oldcode],rtol=2e-10,atol=1e-7)
  sp=S/'outputs'/si[source]['metadata'];assert sha(sp)==si[source]['metadata_SHA256'];sm=load(sp);assert sm['native_measured_genes']==genes and sha(sm['statistics']['path'])==sm['statistics']['SHA256']
  with np.load(sm['statistics']['path'],allow_pickle=False) as f:poolN=int(f['n'].sum());poolCP=f['cp_sum'].sum(axis=0);poolDet=f['detected'].sum(axis=0);np.testing.assert_array_equal(f['native_columns'],cols)
  for g in groups:
   j=g['measurement_subset_code'];sel=g.pop('scope_positions');assert stat['n'][j]==len(sel);g['stable_membership_SHA256']=membership_key(source,ids[q['source_native_position'][sel]].astype(str).tolist());g['native_rows_SHA256']=hashlib.sha256(np.asarray(q['native_row'][sel],'<i8').tobytes()).hexdigest();count=poolN-len(sel);assert count>=0;otherCP=poolCP-stat['cp_sum'][j];otherDet=poolDet-stat['detected'][j];assert otherCP.min()>=-1e-7 and otherDet.min()>=0;otherCP=np.maximum(otherCP,0)
   g.update(own_full_native_gene_marker_candidates_not_DEG=expression_candidates(genes,stat['cp_sum'][j]/len(sel),stat['detected'][j]/len(sel),otherCP/count if count else None,otherDet/count if count else None,.1,20),comparator_scope='other_current_source_common_TNK_and_unknown_parent_controls' if count else 'no_other_cells_in_current_source_common_pool',n_comparator_cells=count,n_source_common_pool_records=poolN,new_identity_applied=False)
  npz=root/'outputs'/(source+'_reentry_annotation_delta_native_v1.npz')
  with npz.open('xb') as f:np.savez_compressed(f,**stat,native_columns=cols,scope_position=np.flatnonzero(code_for_scope>=0),native_row=q['native_row'][code_for_scope>=0],measurement_subset_code=code_for_scope[code_for_scope>=0])
  with np.load(npz,allow_pickle=False) as f:
   for name,value in stat.items():np.testing.assert_array_equal(f[name],value)
  out=dict(source=source,n_parent_groups=len(parents),n_exact_subsets=ng,n_proposed_native_records_directly_read=len(take),n_affected_parent_records=int((code_for_scope>=0).sum()),n_exact_memberships_reused=sum(g['exact_membership_reused'] for g in groups),native_measured_genes=genes,statistics=dict(path=str(npz),SHA256=sha(npz)),groups=groups,proposal_metadata_SHA256=sha(pp),working_arrays_SHA256=sha(wa),accepted_parent_native_metadata_SHA256=sha(ep),accepted_common_pool_metadata_SHA256=sha(sp),all_current_labels_states_and_source_gate_unchanged=True);mp=root/'outputs'/(source+'_reentry_annotation_delta_native_v1.json');save(mp,out);outputs.append(dict(source=source,metadata=mp.name,metadata_SHA256=sha(mp)));totals.update(n_sources=1,n_parent_groups=len(parents),n_exact_subsets=ng,n_proposed_native_records_directly_read=len(take),n_affected_parent_records=out['n_affected_parent_records'],n_exact_memberships_reused=out['n_exact_memberships_reused']);native._fields.clear();del q,old,stat;gc.collect()
 ip=root/'outputs/reentry_annotation_delta_native_index_v1.json';save(ip,outputs);assert totals['n_proposed_native_records_directly_read']==load(P/'provenance/terminal_v1.json')['totals']['n_proposed'];save(root/'provenance/terminal_v1.json',dict(status='ALL_REENTRY_ANNOTATION_DELTA_EXACT_NATIVE_PROFILES_REQUIRE_INDEPENDENT_AUDIT',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),index_SHA256=sha(ip),proposal_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),proposal_audit_terminal_SHA256=sha(A/'provenance/terminal_v1.json'),new_identity_applied=False,new_CM_started=False,wall_seconds=time.monotonic()-start,peak_RSS_bytes_Linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
