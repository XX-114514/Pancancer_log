"""Independent full-scope join and rule replay for common-context proposals."""
import argparse,collections,datetime,gc,hashlib,json,sys,time
from pathlib import Path
import numpy as np

TYPES=['not_proposed','T','Treg','CD4_T','CD8_T','T_gamma_delta_RNA_candidate']
REASONS=['outside_original_unknown_parent_TCR_target','qualified_common_context_T_proposal','common_context_not_T_under_original_group_policy']
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def hist(a):return dict(collections.Counter(map(str,a)))

def main(root):
 start=time.monotonic();P=Path(load(root/'config/spec_v1.json')['producer']);assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 for name,digest in load(P/'config/supervisor_request_v1.json')['script_SHA256'].items():assert sha(P/name)==digest
 cfg=load(P/'config/spec_v1.json');S=Path(cfg['group_profile_producer']);SA=Path(cfg['group_profile_auditor']);NA=Path(cfg['model_native_auditor']);sc=load(S/'config/spec_v1.json');M=Path(sc['model_producer']);mc=load(M/'config/spec_v1.json');R=Path(mc['context_producer']);rc=load(R/'config/spec_v1.json');K=Path(rc['TCR_producer']);kc=load(K/'config/spec_v1.json');J=Path(kc['joined']);F=Path(kc['candidate_producer'])
 for stage in [S,SA,NA,M]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(SA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(S/'provenance/terminal_v1.json');assert load(NA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(M/'provenance/terminal_v1.json')
 w=load(load(rc['working_acceptance'])['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};oldgroups=load(w['group_catalog']['path'])['groups'];L=Path(load(J/'outputs/joined_native_atlas_delivery_v1.json')['local_review']);policy=L/'scripts/native_annotation_review_v1.py';assert sha(policy)==load(L/'config/supervisor_request_v1.json')['script_SHA256']['scripts/native_annotation_review_v1.py'];sys.path.insert(0,str(L/'scripts'))
 from native_annotation_review_v1 import classify_group
 ref=Path(load(F/'config/full_cell_candidate_spec_v1.json')['reference']);models=load(ref/'outputs/reference_models_v1.json');ri={r['source']:r for r in load(ref/'outputs/full_archive_reference_index_v1.json')};fi={r['source']:r for r in load(F/'outputs/full_native_cell_identity_candidate_index_v1.json')};ki={r['source']:r for r in load(K/'outputs/TCR_boundary_proposal_index_v1.json')};oi={r['source']:r for r in load(R/'outputs/unknown_parent_TCR_context_index_v1.json')};mi={r['source']:r for r in load(M/'outputs/common_TNK_reentry_index_v1.json')};si={r['source']:r for r in load(S/'outputs/common_TNK_group_native_profile_index_v1.json')};ip=P/'outputs/common_reentry_identity_proposal_index_v1.json';assert sha(ip)==load(P/'provenance/terminal_v1.json')['index_SHA256'];index=load(ip);assert {r['source'] for r in index}==set(wi);totals=collections.Counter();types=collections.Counter();groupcalls=collections.Counter()
 for entry in index:
  source=entry['source'];mp=P/'outputs'/entry['metadata'];assert sha(mp)==entry['metadata_SHA256'];out=load(mp);assert sha(out['arrays'])==out['arrays_SHA256']
  with np.load(out['arrays'],allow_pickle=False) as f:actual={k:f[k] for k in f.files}
  wr=wi[source];wp=base/wr['metadata'];wa=base/wr['arrays'];assert sha(wp)==wr['metadata_SHA256'] and sha(wa)==wr['arrays_SHA256'];wm=load(wp)
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','annotation_group_code','consensus_major_code','consensus_celltype_code','native_total_RNA','fit_pool_candidate_not_final_QC','reference_family_bitmask']}
  for name in ['native_row','source_native_position']:np.testing.assert_array_equal(actual[name],q[name])
  sp=S/'outputs'/si[source]['metadata'];assert sha(sp)==si[source]['metadata_SHA256'];sm=load(sp);assert sha(sm['statistics']['path'])==sm['statistics']['SHA256']
  with np.load(sm['statistics']['path'],allow_pickle=False) as f:det=f['detected'];ng=f['n']
  modelp=M/'outputs'/mi[source]['metadata'];assert sha(modelp)==mi[source]['metadata_SHA256'];model=load(modelp);assert sha(model['arrays'])==model['arrays_SHA256']
  with np.load(model['arrays'],allow_pickle=False) as f:pool=f['source_scope_position'];codes=f['reentry_group_code'];np.testing.assert_array_equal(f['native_row'],q['native_row'][pool])
  rr=ri[source];rp=ref/'outputs'/rr['metadata'];ra=ref/'outputs'/rr['arrays'];assert sha(rp)==rr['metadata_sha256'] and sha(ra)==rr['arrays_sha256'];rm=load(rp);pos=q['source_native_position'][pool]
  with np.load(ra,allow_pickle=False) as f:np.testing.assert_array_equal(q['native_row'][pool],f['native_row'][pos]);panels=f['native_panel_bitmask'][pos];low=f['immune_low__unique_code'][pos]
  n=len(q['native_row']);common=np.full(n,-1,np.int32);common[pool]=codes;localT=np.zeros(n,bool);genes=sm['native_measured_genes'];bygene={g:j for j,g in enumerate(genes)};ngroup=len(sm['groups']);panelcounts={name:np.bincount(codes,weights=(panels&(1<<j))!=0,minlength=ngroup).astype(np.int64) for j,name in enumerate(rm['native_panel_bit_order'])};group_records=[]
  for code,g in enumerate(sm['groups']):
   sel=codes==code;count=int(sel.sum());assert count==ng[code];panel={name:int(v[code])/count for name,v in panelcounts.items()};detection={gene:float(det[code,bygene[gene]]/count) for gene in rm['review_genes'] if gene in bygene};reference=hist(np.asarray(models['immune_low']['labels'],object)[low[sel][low[sel]>=0]]);call=classify_group(panel,detection,reference,count)
   if call['local_major'] in {'T','T_NK'}:localT[pool[sel]]=True
   group_records.append(dict(source_reentry_group_code=code,n_exact_group=count,stable_membership_SHA256=g['stable_membership_SHA256'],local_major_proposal=call['local_major'],local_celltype_proposal=call['local_celltype'],local_annotation_evidence_status=call['annotation_evidence_status'],local_native_supported_families=call['native_supported_families'],native_panel_fraction=panel,source_group_rule_unchanged=True,original_group_policy_SHA256=sha(policy),native_group_profile_metadata_SHA256=sha(sp),reference_unique_low_counts=reference));groupcalls[call['local_major']]+=1
  assert out['groups']==group_records;np.testing.assert_array_equal(actual['source_common_group_code'],common);np.testing.assert_array_equal(actual['new_common_local_T_context'],localT)
  oldp=R/'outputs'/oi[source]['metadata'];assert sha(oldp)==oi[source]['metadata_SHA256'];oldm=load(oldp);assert sha(oldm['arrays'])==oldm['arrays_SHA256']
  with np.load(oldm['arrays'],allow_pickle=False) as f:target=f['target_TCR_unknown_parent'];np.testing.assert_array_equal(f['native_row'],q['native_row'])
  np.testing.assert_array_equal(actual['original_unknown_parent_TCR_target'],target);kp=K/'outputs'/ki[source]['metadata'];assert sha(kp)==ki[source]['metadata_SHA256'];km=load(kp);assert sha(K/'outputs'/km['arrays'])==km['arrays_SHA256']
  with np.load(K/'outputs'/km['arrays'],allow_pickle=False) as f:assert (f['review_code'][target]==7).all() and f['TCR_core'][target].all() and f['reference_low_T'][target].all() and (f['competing_core_mask'][target]==0).all() and (f['unassessed_old_competitor_mask'][target]==0).all()
  major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];assert (major[target]=='unresolved').all() and (q['reference_family_bitmask'][target]==1).all() and q['fit_pool_candidate_not_final_QC'][target].all() and (q['native_total_RNA'][target]>0).all() and (common[target]>=0).all();proposed=np.logical_and(localT,target);np.testing.assert_array_equal(actual['proposed_T'],proposed)
  fr=fi[source];fp=F/'outputs'/fr['metadata'];fa=F/'outputs'/fr['arrays'];assert sha(fp)==fr['metadata_sha256'] and sha(fa)==fr['arrays_sha256'];assert load(fp)['fine_candidate_order'][:4]==TYPES[2:]
  with np.load(fa,allow_pickle=False) as f:pos=q['source_native_position'];np.testing.assert_array_equal(q['native_row'],f['native_row'][pos]);both=(f['own_fine_bitmask'][pos]&f['reference_fine_bitmask'][pos]&15).astype(np.uint8)
  lookup=np.asarray([v&~2 if v&1 else v for v in range(16)],np.uint8);direct=lookup[both];fine_lookup=np.asarray([({1:2,2:3,4:4,8:5}.get(int(v),1)) for v in range(16)],np.uint8);finecode=np.where(proposed,fine_lookup[direct],0).astype(np.uint8);reason=np.zeros(n,np.uint8);reason[target&~localT]=2;reason[proposed]=1
  for name,value in [('proposed_celltype_code',finecode),('direct_fine_intersection_bitmask',direct),('review_code',reason)]:np.testing.assert_array_equal(actual[name],value)
  assert out['review_dictionary']==REASONS and out['celltype_dictionary']==TYPES;before=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']];after=before.copy();after[proposed]=np.asarray(TYPES,object)[finecode[proposed]];current={g['source_annotation_group_code']:g for g in oldgroups if g['source']==source};affected=[]
  for code in sorted(set(map(int,q['annotation_group_code'][proposed]))):
   sel=q['annotation_group_code']==code;g=current[code];affected.append(dict(previous_group_code=code,previous_group_number=g['global_subcluster_number'],previous_group_SHA256=g['stable_group_sha256'],n_previous_group=int(sel.sum()),n_proposed=int((sel&proposed).sum()),after_celltype_counts=hist(after[sel])))
  counts=hist(np.asarray(TYPES,object)[finecode[proposed]]);assert affected==out['affected_previous_exact_groups'] and counts==out['proposed_celltype_counts'] and out['n_scope_records']==n and out['n_original_target']==int(target.sum()) and out['n_proposed']==int(proposed.sum())
  for k,p in [('working_metadata_SHA256',wp),('working_arrays_SHA256',wa),('common_model_metadata_SHA256',modelp),('common_native_profile_metadata_SHA256',sp),('reference_metadata_SHA256',rp),('reference_arrays_SHA256',ra),('frozen_fine_metadata_SHA256',fp),('original_unknown_parent_target_metadata_SHA256',oldp),('original_group_policy_SHA256',policy)]:assert out[k]==sha(p)
  totals.update(n_sources=1,n_scope_records=n,n_original_target=int(target.sum()),n_proposed=int(proposed.sum()),n_affected_previous_groups=len(affected),n_new_common_groups=ngroup);types.update(counts);save(root/'outputs'/(source+'_common_reentry_identity_audit_v1.json'),dict(source=source,all_memberships_contexts_frozen_rules_and_proposal_fields_replayed=True));del q,actual;gc.collect()
 sp=P/'outputs/common_reentry_identity_proposal_summary_v1.json';assert sha(sp)==load(P/'provenance/terminal_v1.json')['summary_SHA256'];summary=load(sp);assert summary['totals']==dict(totals) and summary['proposed_celltype_counts']==dict(types) and summary['new_group_original_policy_major_calls']==dict(groupcalls) and summary['source_group_policy_SHA256']==sha(policy)
 save(root/'provenance/terminal_v1.json',dict(status='ALL_COMMON_CONTEXT_IDENTITY_PROPOSAL_JOINS_AND_FROZEN_RULES_ACCEPTED_NOT_APPLIED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),proposed_celltype_counts=dict(types),producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),producer_index_SHA256=sha(ip),original_group_policy_SHA256=sha(policy),new_production_proposal_code_not_imported=True,original_accepted_group_classifier_reused_unchanged=True,not_orthogonal_biological_validation=True,exact_membership_markers_and_working_release_pending=True,new_identity_applied=False,new_CM_started=False,wall_seconds=time.monotonic()-start))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
