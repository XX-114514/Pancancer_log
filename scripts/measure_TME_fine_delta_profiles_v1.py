"""Own full-gene evidence for every fine-label delta, with exact-profile reuse."""
import argparse,collections,datetime,gc,hashlib,json,sys,time
from pathlib import Path
import numpy as np
from TME_accepted_native_profiles_v1 import AcceptedNativeProfiles
from TME_native_measurement_core_v1 import canonical_native,update_group_stats,expression_candidates

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def subtract_nonnegative(total,own,integer=False):
 d=total-own
 if d.min() < (0 if integer else -1e-7):raise ValueError('Comparator subtraction is not a nonnegative cell pool')
 return np.maximum(d,0)

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['proposal']);A=Path(cfg['proposal_audit']);pc=load(P/'config/spec_v1.json')
 for s in [P,A]:assert load(s/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json')
 J=Path(pc['joined']);sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(root/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
 from native_panel_aggregate_v1 import gather_csr_rows
 from atlas_identity_release_keys_v1 import membership_key
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 profiles=AcceptedNativeProfiles(pc['working_acceptance'],pc['native_evidence'],pc['native_audit'],cfg['previous_changed_profiles'],cfg['previous_changed_audit']);working=profiles.working;base=Path(working['index']['path']).parent;wi={r['source']:r for r in load(working['index']['path'])};native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest'])
 ip=P/'outputs/direct_T_fine_proposal_index_v1.json';assert sha(ip)==load(P/'provenance/terminal_v1.json')['index_SHA256'];results=[];totals=collections.Counter()
 for entry in load(ip):
  if entry['n_fine_changes']==0:continue
  storage_gate(root,pc);source=entry['source'];mp=P/'outputs'/entry['metadata'];assert sha(mp)==entry['metadata_SHA256'];m=load(mp);pa=P/'outputs'/m['arrays'];assert sha(pa)==m['arrays_SHA256']
  with np.load(pa,allow_pickle=False) as f:proposal=f['proposed_fine_code'];positions=f['source_native_position']
  wr=wi[source];wm=load(base/wr['metadata']);wa=base/wr['arrays'];assert sha(wa)==wr['arrays_SHA256']
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','annotation_group_code','consensus_celltype_code','native_total_RNA']}
  np.testing.assert_array_equal(positions,q['source_native_position']);fine=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']];after=fine.copy();changed=proposal>0;after[changed]=np.asarray(m['fine_dictionary'],object)[proposal[changed]]
  p=profiles.source(source);genes=p['genes'];cols=p['native_columns'];oldstat=p['stat'];groups=[];raw_domains=[];code_for_scope=np.full(len(fine),-1,np.int32)
  for oldcode in sorted(set(map(int,q['annotation_group_code'][changed&(q['annotation_group_code']>=0)]))):
   oldsel=q['annotation_group_code']==oldcode;parent=p['groups'][oldcode];labels=sorted(set(after[oldsel]));reuse=len(labels)==1
   for label in labels:
    sel=np.flatnonzero(oldsel&(after==label));j=len(groups);code_for_scope[sel]=j;groups.append({'measurement_subset_code':j,'previous_group_code':oldcode,'previous_group_number':parent['global_subcluster_number'],'previous_group_SHA256':parent['stable_group_sha256'],'fine_celltype':str(label),'n_exact_subset':len(sel),'exact_membership_reused':reuse,'accepted_profile_reference':p['profile_references'][oldcode] if reuse else None,'scope_positions':sel})
   if not reuse:raw_domains.append(np.flatnonzero(oldsel))
  n=len(groups);stat={'n':np.zeros(n,np.int64),**{k:np.zeros((n,len(genes)),np.float64) for k in ['cp_sum','log_sum','detected']}}
  for g in groups:
   if g['exact_membership_reused']:
    for k in stat:stat[k][g['measurement_subset_code']]=oldstat[k][g['previous_group_code']]
  take=np.sort(np.concatenate(raw_domains)) if raw_domains else np.asarray([],np.int64);ids=native.native_ids(source);assert text_axis_sha256(ids)==wm['native_ids_source_SHA256'];np.testing.assert_array_equal(native.rows(source)[positions],q['native_row'])
  key=native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key'];csr=native.store(source)[key]
  for low in range(0,len(take),2048):
   sel=take[low:low+2048];raw=canonical_native(gather_csr_rows(csr,q['native_row'][sel],max_span=8192));rna=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(rna,q['native_total_RNA'][sel]);update_group_stats(stat,raw[:,cols].tocsr(),rna,code_for_scope[sel])
  assert (stat['n']>0).all() if n else True
  for oldcode in sorted({g['previous_group_code'] for g in groups}):
   child=[g['measurement_subset_code'] for g in groups if g['previous_group_code']==oldcode]
   for k in ['n','detected']:np.testing.assert_array_equal(stat[k][child].sum(axis=0),oldstat[k][oldcode])
   for k in ['cp_sum','log_sum']:np.testing.assert_allclose(stat[k][child].sum(axis=0),oldstat[k][oldcode],rtol=2e-10,atol=1e-7)
  tn=p['T_NK_group_mask'];poolN=int(oldstat['n'][tn].sum());poolCP=oldstat['cp_sum'][tn].sum(axis=0);poolDet=oldstat['detected'][tn].sum(axis=0)
  for g in groups:
   j=g['measurement_subset_code'];sel=g.pop('scope_positions');g['stable_membership_SHA256']=membership_key(source,ids[positions[sel]].astype(str).tolist());g['native_rows_SHA256']=hashlib.sha256(np.asarray(q['native_row'][sel],'<i8').tobytes()).hexdigest()
   if g['exact_membership_reused']:assert g['stable_membership_SHA256']==g['previous_group_SHA256']
   nc=poolN-int(stat['n'][j]);cp=poolCP;det=poolDet;scope='other_current_source_T_NK_fitted_group_cells'
   if nc==0:nc=int(oldstat['n'].sum()-stat['n'][j]);cp=oldstat['cp_sum'].sum(axis=0);det=oldstat['detected'].sum(axis=0);scope='other_current_source_fitted_group_cells_no_T_NK_comparator'
   assert nc>=0;other=subtract_nonnegative(cp,stat['cp_sum'][j])/nc if nc else None;otherdet=subtract_nonnegative(det,stat['detected'][j],True)/nc if nc else None
   g['own_full_native_gene_marker_candidates_not_DEG']=expression_candidates(genes,stat['cp_sum'][j]/stat['n'][j],stat['detected'][j]/stat['n'][j],other,otherdet,.1,20);g.update(comparator_scope=scope,n_comparator_cells=nc,new_display_number_not_assigned=True)
  sp=root/'outputs'/(source+'_fine_delta_native_profiles_v1.npz')
  with sp.open('xb') as f:np.savez_compressed(f,**stat,native_columns=cols,scope_position=np.flatnonzero(code_for_scope>=0),native_row=q['native_row'][code_for_scope>=0],measurement_subset_code=code_for_scope[code_for_scope>=0])
  with np.load(sp,allow_pickle=False) as f:
   for k,v in stat.items():np.testing.assert_array_equal(f[k],v)
  out={'source':source,'n_exact_subsets':n,'n_parent_groups':len({g['previous_group_code'] for g in groups}),'n_raw_records_read':len(take),'n_exact_memberships_reused':sum(g['exact_membership_reused'] for g in groups),'n_fine_changes':int(changed.sum()),'n_unfitted_fine_changes':int((changed&(q['annotation_group_code']<0)).sum()),'native_measured_genes':genes,'statistics':{'path':str(sp),'SHA256':sha(sp)},'groups':groups,'proposal_metadata_SHA256':sha(mp),'working_arrays_SHA256':sha(wa),'source_complete_profile_metadata_SHA256':p['native_metadata_SHA256'],'current_labels_not_applied':True}
  om=root/'outputs'/(source+'_fine_delta_native_profiles_v1.json');save(om,out);results.append({'source':source,'metadata':om.name,'metadata_SHA256':sha(om)});totals.update(n_sources=1,n_parent_groups=out['n_parent_groups'],n_exact_subsets=n,n_raw_records_read=len(take),n_exact_memberships_reused=out['n_exact_memberships_reused'],n_fine_changes=int(changed.sum()),n_unfitted_fine_changes=out['n_unfitted_fine_changes']);native._fields.clear();profiles.cache=None;del p,q,oldstat,stat;gc.collect()
 op=root/'outputs/fine_delta_native_profile_index_v1.json';save(op,results);assert totals['n_fine_changes']==load(P/'provenance/terminal_v1.json')['totals']['n_fine_changes']
 save(root/'provenance/terminal_v1.json',{'status':'ALL_FINE_DELTA_OWN_NATIVE_PROFILES_READY_INDEPENDENT_AUDIT_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals),'index_SHA256':sha(op),'current_labels_states_and_CM_unchanged':True,'wall_seconds':time.monotonic()-start})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
