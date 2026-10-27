"""Full native own-group expression and existing RNA context for new groups.

The provisional common-space groups are not new cell identities. Read original
counts, preserve all measurement masks, and require an independent profile audit
plus model-native acceptance before any annotation release.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
import psutil
from TME_native_measurement_core_v1 import canonical_native,update_group_stats,expression_candidates

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def hist(v):return dict(collections.Counter(map(str,v)))

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['model_producer']);A=Path(cfg['model_scope_auditor'])
 for s in [P,A]:assert load(s/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json');pc=load(P/'config/spec_v1.json');R=Path(pc['context_producer']);rc=load(R/'config/spec_v1.json');K=Path(rc['TCR_producer']);kc=load(K/'config/spec_v1.json');J=Path(kc['joined']);E=Path(kc['native_evidence']);F=Path(kc['candidate_producer']);accepted=load(rc['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];w=load(accepted['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])}
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(root/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
 from native_panel_aggregate_v1 import gather_csr_rows
 from atlas_identity_release_keys_v1 import membership_key
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);ei={r['source']:r for r in load(E/'outputs/TME_full_native_evidence_index_v1.json')};fi={r['source']:r for r in load(F/'outputs/full_native_cell_identity_candidate_index_v1.json')};ki={r['source']:r for r in load(K/'outputs/TCR_boundary_proposal_index_v1.json')};ri={r['source']:r for r in load(R/'outputs/unknown_parent_TCR_context_index_v1.json')};ip=P/'outputs/common_TNK_reentry_index_v1.json';assert sha(ip)==load(P/'provenance/terminal_v1.json')['index_SHA256'];totals=collections.Counter();index=[];resource_cfg=dict(kc,output_cap_bytes=cfg['output_cap_bytes'])
 for item in load(ip):
  source=item['source'];storage_gate(root,resource_cfg);mp=P/'outputs'/item['metadata'];assert sha(mp)==item['metadata_SHA256'];m=load(mp);assert sha(m['arrays'])==m['arrays_SHA256']
  with np.load(m['arrays'],allow_pickle=False) as f:pool=f['source_scope_position'];codes=f['reentry_group_code'];rows=f['native_row'];known=f['baseline_known_T_NK'];control=f['unresolved_parent_control']
  n=len(pool);ng=m['n_groups'];wr=wi[source];wp=base/wr['metadata'];wa=base/wr['arrays'];assert sha(wp)==wr['metadata_SHA256'] and sha(wa)==wr['arrays_SHA256'];wm=load(wp)
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','native_total_RNA','consensus_major_code','consensus_celltype_code','native_family_bitmask','reference_family_bitmask','RNA_state_native_logCP10K']}
  np.testing.assert_array_equal(rows,q['native_row'][pool]);ep=E/'outputs'/ei[source]['metadata'];assert sha(ep)==ei[source]['metadata_SHA256'];e=load(ep);cols=np.asarray(e['native_columns']);genes=e['native_measured_gene_names'];bound=ng*len(cols)*32+cfg['RAM_reserve_bytes'];assert bound<psutil.virtual_memory().available
  stat={'n':np.zeros(ng,np.int64),**{k:np.zeros((ng,len(cols)),np.float64) for k in ['cp_sum','log_sum','detected']}};csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']];nnz=0
  for lo in range(0,n,2048):
   sel=pool[lo:lo+2048];raw=canonical_native(gather_csr_rows(csr,q['native_row'][sel],max_span=8192));rna=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(rna,q['native_total_RNA'][sel]);update_group_stats(stat,raw[:,cols].tocsr(),rna,codes[lo:lo+len(sel)]);nnz+=raw.nnz
  np.testing.assert_array_equal(stat['n'],np.bincount(codes,minlength=ng));assert int(stat['n'].sum())==n
  fr=fi[source];fp=F/'outputs'/fr['metadata'];fa=F/'outputs'/fr['arrays'];assert sha(fp)==fr['metadata_sha256'] and sha(fa)==fr['arrays_sha256'];fm=load(fp)
  with np.load(fa,allow_pickle=False) as f:
   pos=q['source_native_position'][pool];np.testing.assert_array_equal(rows,f['native_row'][pos]);own=f['own_fine_bitmask'][pos];reference=f['reference_fine_bitmask'][pos]
  direct=(own&reference).copy();direct[(direct&1)!=0]&=np.asarray(65533,dtype=direct.dtype)
  kp=K/'outputs'/ki[source]['metadata'];assert sha(kp)==ki[source]['metadata_SHA256'];km=load(kp);ka=K/'outputs'/km['arrays'];assert sha(ka)==km['arrays_SHA256']
  with np.load(ka,allow_pickle=False) as f:np.testing.assert_array_equal(q['native_row'],f['native_row']);tcore=f['TCR_core'][pool];competitor=f['competing_core_mask'][pool];unknown=f['unassessed_old_competitor_mask'][pool];lowT=f['reference_low_T'][pool]
  rp=R/'outputs'/ri[source]['metadata'];assert sha(rp)==ri[source]['metadata_SHA256'];rm=load(rp);assert sha(rm['arrays'])==rm['arrays_SHA256']
  with np.load(rm['arrays'],allow_pickle=False) as f:target=f['target_TCR_unknown_parent'][pool]
  major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code'][pool]];fine=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code'][pool]];ids=native.native_ids(source);assert text_axis_sha256(ids)==wm['native_ids_source_SHA256'];units=native.frozen_expression_units(source)[q['source_native_position'][pool]].astype(str);groups=[]
  totalCP=stat['cp_sum'].sum(axis=0);totalDet=stat['detected'].sum(axis=0)
  for code in range(ng):
   sel=np.flatnonzero(codes==code);count=n-len(sel);otherCP=totalCP-stat['cp_sum'][code];otherDet=totalDet-stat['detected'][code];assert otherCP.min()>=-1e-7 and otherDet.min()>=0;otherCP=np.maximum(otherCP,0);candidates=expression_candidates(genes,stat['cp_sum'][code]/len(sel),stat['detected'][code]/len(sel),otherCP/count if count else None,otherDet/count if count else None,.1,20)
   state={}
   for j,name in enumerate(wm['state_order']):
    values=q['RNA_state_native_logCP10K'][pool[sel],j];good=np.isfinite(values);state[name]=dict(n_measured=int(good.sum()),mean_native_logCP10K=float(values[good].mean(dtype=np.float64)) if good.any() else None,measurement_coverage=wm['state_measurement_coverage'][name],functional_identity_assigned=False)
   missing=np.isin(np.char.lower(np.char.strip(units[sel])),['','nan','none','na','n/a','unknown','__missing__','__source_value_missing__']);groups.append(dict(source_reentry_group_code=code,n_exact_group=len(sel),stable_membership_SHA256=membership_key(source,ids[q['source_native_position'][pool[sel]]].astype(str).tolist()),native_rows_SHA256=hashlib.sha256(np.asarray(rows[sel],'<i8').tobytes()).hexdigest(),n_baseline_known_T_NK=int(known[sel].sum()),n_unknown_parent_control=int(control[sel].sum()),n_strict_TCR_unknown_parent_candidate=int(target[sel].sum()),n_TCR_core=int(tcore[sel].sum()),n_TCR_core_with_coherent_competitor=int((tcore[sel]&(competitor[sel]!=0)).sum()),coherent_competitor_mask_counts=hist(competitor[sel]),unassessed_old_competitor_mask_counts=hist(unknown[sel]),n_unique_low_reference_T=int(lowT[sel].sum()),baseline_major_counts=hist(major[sel]),baseline_celltype_counts=hist(fine[sel]),baseline_labels_are_not_independent_truth=True,native_family_mask_counts=hist(q['native_family_bitmask'][pool[sel]]),reference_family_mask_counts=hist(q['reference_family_bitmask'][pool[sel]]),direct_own_and_reference_fine_cue_counts={name:int(((direct[sel]&(1<<j))!=0).sum()) for j,name in enumerate(fm['fine_candidate_order'])},n_declared_expression_units_not_libraries_or_patients=len(set(units[sel][~missing])),n_records_with_missing_declared_unit=int(missing.sum()),descriptive_existing_RNA_states=state,own_full_native_gene_marker_candidates_not_DEG=candidates,comparator_scope='other_current_source_common_TNK_and_unknown_control_groups' if count else 'no_other_group_in_current_source_common_pool',n_comparator_cells=count,new_identity_assigned=False))
  spath=root/'outputs'/(source+'_common_TNK_group_native_profiles_v1.npz')
  with spath.open('xb') as f:np.savez_compressed(f,**stat,native_columns=cols)
  with np.load(spath,allow_pickle=False) as f:
   for name,value in stat.items():np.testing.assert_array_equal(f[name],value)
  op=root/'outputs'/(source+'_common_TNK_group_native_profiles_v1.json');out=dict(source=source,n_model_pool_records=n,n_exact_groups=ng,n_native_nnz_read=nnz,native_measured_genes=genes,statistics=dict(path=str(spath),SHA256=sha(spath)),groups=groups,model_group_arrays=dict(path=m['arrays'],SHA256=m['arrays_SHA256']),model_metadata_SHA256=sha(mp),working_arrays_SHA256=sha(wa),complete_native_metadata_SHA256=sha(ep),TCR_metadata_SHA256=sha(kp),direct_fine_metadata_SHA256=sha(fp),no_current_working_identity_or_state_changed=True,independent_profile_audit_and_model_native_acceptance_required=True);save(op,out);index.append(dict(source=source,metadata=op.name,metadata_SHA256=sha(op)));totals.update(n_sources=1,n_model_pool_records=n,n_exact_groups=ng,n_native_nnz_read=nnz);native._fields.clear();del q,stat;gc.collect()
  save(root/'outputs'/(source+'_common_profile_checkpoint_v1.json'),dict(completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals)))
 assert totals['n_model_pool_records']==load(P/'provenance/terminal_v1.json')['totals']['n_fit_pool'];ip=root/'outputs/common_TNK_group_native_profile_index_v1.json';save(ip,index)
 save(root/'provenance/terminal_v1.json',dict(status='ALL_COMMON_TNK_GROUP_NATIVE_PROFILES_AND_EXISTING_CONTEXT_READY_FOR_INDEPENDENT_AUDIT',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),index_SHA256=sha(ip),model_producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),model_scope_audit_terminal_SHA256=sha(A/'provenance/terminal_v1.json'),no_identity_or_functional_state_assigned=True,new_CM_started=False,wall_seconds=time.monotonic()-start,peak_RSS_bytes_Linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
