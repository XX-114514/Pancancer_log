"""Replace only the frozen local-context gate with the qualified common graph.

Reuse the original source-group classifier and direct fine-cue rules unchanged.
All previous TCR/reference/competitor/fit gates remain, and only the accepted
unknown-parent target is eligible. This stage writes proposals, not annotations.
"""
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
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');S=Path(cfg['group_profile_producer']);SA=Path(cfg['group_profile_auditor']);NA=Path(cfg['model_native_auditor']);sc=load(S/'config/spec_v1.json');P=Path(sc['model_producer']);pc=load(P/'config/spec_v1.json');R=Path(pc['context_producer']);rc=load(R/'config/spec_v1.json');K=Path(rc['TCR_producer']);kc=load(K/'config/spec_v1.json');J=Path(kc['joined']);F=Path(kc['candidate_producer'])
 for stage in [S,SA,NA,P,Path(sc['model_scope_auditor'])]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(SA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(S/'provenance/terminal_v1.json');assert load(NA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json')
 accepted=load(rc['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];w=load(accepted['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};catalog=load(w['group_catalog']['path'])['groups']
 L=Path(load(J/'outputs/joined_native_atlas_delivery_v1.json')['local_review']);policy=L/'scripts/native_annotation_review_v1.py';request=load(L/'config/supervisor_request_v1.json');assert sha(policy)==request['script_SHA256']['scripts/native_annotation_review_v1.py'];sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(L/'scripts'))
 from native_annotation_review_v1 import classify_group
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 ref=Path(load(F/'config/full_cell_candidate_spec_v1.json')['reference']);models=load(ref/'outputs/reference_models_v1.json');ri={r['source']:r for r in load(ref/'outputs/full_archive_reference_index_v1.json')};fi={r['source']:r for r in load(F/'outputs/full_native_cell_identity_candidate_index_v1.json')};ki={r['source']:r for r in load(K/'outputs/TCR_boundary_proposal_index_v1.json')};oi={r['source']:r for r in load(R/'outputs/unknown_parent_TCR_context_index_v1.json')};pi={r['source']:r for r in load(P/'outputs/common_TNK_reentry_index_v1.json')};ip=S/'outputs/common_TNK_group_native_profile_index_v1.json';assert sha(ip)==load(S/'provenance/terminal_v1.json')['index_SHA256'];totals=collections.Counter();types=collections.Counter();groupcalls=collections.Counter();index=[];resource_cfg=dict(kc,output_cap_bytes=cfg.get('output_cap_bytes',1073741824))
 for entry in load(ip):
  source=entry['source'];storage_gate(root,resource_cfg);sp=S/'outputs'/entry['metadata'];assert sha(sp)==entry['metadata_SHA256'];sm=load(sp);assert sha(sm['statistics']['path'])==sm['statistics']['SHA256']
  with np.load(sm['statistics']['path'],allow_pickle=False) as f:detected=f['detected'];groupn=f['n']
  wp=base/wi[source]['metadata'];wa=base/wi[source]['arrays'];assert sha(wp)==wi[source]['metadata_SHA256'] and sha(wa)==wi[source]['arrays_SHA256'];wm=load(wp)
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','annotation_group_code','consensus_major_code','consensus_celltype_code','native_total_RNA','fit_pool_candidate_not_final_QC','reference_family_bitmask']}
  modelp=P/'outputs'/pi[source]['metadata'];assert sha(modelp)==pi[source]['metadata_SHA256'];model=load(modelp);assert sha(model['arrays'])==model['arrays_SHA256']
  with np.load(model['arrays'],allow_pickle=False) as f:pool=f['source_scope_position'];codes=f['reentry_group_code'];np.testing.assert_array_equal(f['native_row'],q['native_row'][pool])
  rr=ri[source];rp=ref/'outputs'/rr['metadata'];ra=ref/'outputs'/rr['arrays'];assert sha(rp)==rr['metadata_sha256'] and sha(ra)==rr['arrays_sha256'];rm=load(rp);positions=q['source_native_position'][pool]
  with np.load(ra,allow_pickle=False) as f:
   np.testing.assert_array_equal(f['native_row'][positions],q['native_row'][pool]);panel=f['native_panel_bitmask'][positions];low=f['immune_low__unique_code'][positions];assert not np.any((low>=0)&~f['immune_reference_scope'][positions])
  genes=sm['native_measured_genes'];bygene={g:j for j,g in enumerate(genes)};groups=[];localT=np.zeros(len(q['native_row']),bool);common=np.full(len(localT),-1,np.int32)
  common[pool]=codes
  for code,g in enumerate(sm['groups']):
   sel=np.flatnonzero(codes==code);n=len(sel);assert n==groupn[code]==g['n_exact_group'];panels={name:float(((panel[sel]&(1<<j))!=0).mean()) for j,name in enumerate(rm['native_panel_bit_order'])};detection={gene:float(detected[code,bygene[gene]]/n) for gene in rm['review_genes'] if gene in bygene};lowvalid=low[sel][low[sel]>=0];reference=hist(np.asarray(models['immune_low']['labels'],object)[lowvalid]);result=classify_group(panels,detection,reference,n)
   if result['local_major'] in ['T','T_NK']:localT[pool[sel]]=True
   groups.append(dict(source_reentry_group_code=code,n_exact_group=n,stable_membership_SHA256=g['stable_membership_SHA256'],local_major_proposal=result['local_major'],local_celltype_proposal=result['local_celltype'],local_annotation_evidence_status=result['annotation_evidence_status'],local_native_supported_families=result['native_supported_families'],native_panel_fraction=panels,source_group_rule_unchanged=True,original_group_policy_SHA256=sha(policy),native_group_profile_metadata_SHA256=sha(sp),reference_unique_low_counts=reference));groupcalls[result['local_major']]+=1
  oldp=R/'outputs'/oi[source]['metadata'];assert sha(oldp)==oi[source]['metadata_SHA256'];oldmeta=load(oldp);assert sha(oldmeta['arrays'])==oldmeta['arrays_SHA256']
  with np.load(oldmeta['arrays'],allow_pickle=False) as f:target=f['target_TCR_unknown_parent'];np.testing.assert_array_equal(f['native_row'],q['native_row'])
  kp=K/'outputs'/ki[source]['metadata'];assert sha(kp)==ki[source]['metadata_SHA256'];km=load(kp);ka=K/'outputs'/km['arrays'];assert sha(ka)==km['arrays_SHA256']
  with np.load(ka,allow_pickle=False) as f:
   assert (f['review_code'][target]==7).all() and f['TCR_core'][target].all() and (f['competing_core_mask'][target]==0).all() and (f['unassessed_old_competitor_mask'][target]==0).all() and f['reference_low_T'][target].all()
  major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];assert (major[target]=='unresolved').all() and q['fit_pool_candidate_not_final_QC'][target].all() and (q['native_total_RNA'][target]>0).all() and (q['reference_family_bitmask'][target]==1).all() and (common[target]>=0).all();proposed=target&localT
  fr=fi[source];fp=F/'outputs'/fr['metadata'];fa=F/'outputs'/fr['arrays'];assert sha(fp)==fr['metadata_sha256'] and sha(fa)==fr['arrays_sha256'];assert load(fp)['fine_candidate_order'][:4]==TYPES[2:];pos=q['source_native_position']
  with np.load(fa,allow_pickle=False) as f:np.testing.assert_array_equal(f['native_row'][pos],q['native_row']);direct=(f['own_fine_bitmask'][pos]&f['reference_fine_bitmask'][pos]&15).astype(np.uint8)
  direct[(direct&1)!=0]&=13;finecode=np.zeros(len(target),np.uint8);finecode[proposed]=1
  for j in range(4):finecode[proposed&(direct==(1<<j))]=j+2
  reason=np.zeros(len(target),np.uint8);reason[target]=2;reason[proposed]=1;before=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']];after=before.copy();after[proposed]=np.asarray(TYPES,object)[finecode[proposed]];affected=[]
  current={g['source_annotation_group_code']:g for g in catalog if g['source']==source}
  for oldcode in sorted(set(map(int,q['annotation_group_code'][proposed]))):
   sel=q['annotation_group_code']==oldcode;g=current[oldcode];affected.append(dict(previous_group_code=oldcode,previous_group_number=g['global_subcluster_number'],previous_group_SHA256=g['stable_group_sha256'],n_previous_group=int(sel.sum()),n_proposed=int((sel&proposed).sum()),after_celltype_counts=hist(after[sel])))
  op=root/'outputs'/(source+'_common_reentry_identity_proposals_v1.npz')
  with op.open('xb') as f:np.savez_compressed(f,native_row=q['native_row'],source_native_position=q['source_native_position'],original_unknown_parent_TCR_target=target,source_common_group_code=common,new_common_local_T_context=localT,proposed_T=proposed,proposed_celltype_code=finecode,direct_fine_intersection_bitmask=direct,review_code=reason)
  result=dict(source=source,n_scope_records=len(target),n_original_target=int(target.sum()),n_proposed=int(proposed.sum()),proposed_celltype_counts=hist(np.asarray(TYPES,object)[finecode[proposed]]),review_dictionary=REASONS,celltype_dictionary=TYPES,arrays=str(op),arrays_SHA256=sha(op),groups=groups,affected_previous_exact_groups=affected,working_metadata_SHA256=sha(wp),working_arrays_SHA256=sha(wa),common_model_metadata_SHA256=sha(modelp),common_native_profile_metadata_SHA256=sha(sp),reference_metadata_SHA256=sha(rp),reference_arrays_SHA256=sha(ra),frozen_fine_metadata_SHA256=sha(fp),original_unknown_parent_target_metadata_SHA256=sha(oldp),original_group_policy_SHA256=sha(policy),current_working_labels_states_and_source_gate_unchanged=True,new_identity_applied=False,independent_proposal_and_exact_membership_marker_audits_required=True)
  mp=root/'outputs'/(source+'_common_reentry_identity_proposals_v1.json');save(mp,result);index.append(dict(source=source,metadata=mp.name,metadata_SHA256=sha(mp)));totals.update(n_sources=1,n_scope_records=len(target),n_original_target=int(target.sum()),n_proposed=int(proposed.sum()),n_affected_previous_groups=len(affected),n_new_common_groups=len(groups));types.update(result['proposed_celltype_counts']);del q;gc.collect()
 ip=root/'outputs/common_reentry_identity_proposal_index_v1.json';save(ip,index);summary=dict(totals=dict(totals),proposed_celltype_counts=dict(types),new_group_original_policy_major_calls=dict(groupcalls),source_group_policy_SHA256=sha(policy),old_TCR_and_fine_rules_unchanged=True,only_local_context_replaced_with_current_qualified_common_space=True,no_current_identity_or_QC_or_functional_state_applied=True,new_CM_started=False);sp=root/'outputs/common_reentry_identity_proposal_summary_v1.json';save(sp,summary)
 assert totals['n_original_target']==load(R/'provenance/terminal_v1.json')['totals']['n_candidate_records'];save(root/'provenance/terminal_v1.json',dict(status='COMMON_SPACE_CONTEXT_IDENTITY_PROPOSALS_REQUIRE_INDEPENDENT_AUDIT_AND_EXACT_MARKERS',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),proposed_celltype_counts=dict(types),index_SHA256=sha(ip),summary_SHA256=sha(sp),original_group_policy_SHA256=sha(policy),current_labels_and_new_CM_unchanged=True,wall_seconds=time.monotonic()-start))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
