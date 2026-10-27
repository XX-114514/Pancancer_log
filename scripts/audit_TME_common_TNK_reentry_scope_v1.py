"""Full-source pool/model-axis audit, separate from native feature replay.

Consumes immutable completed-source outputs while the fitter runs. A successful
terminal certifies scope, axes and preserved old annotations; it deliberately
does not certify native feature selection, PCA arithmetic or biological identity.
"""
import argparse,collections,datetime,gc,hashlib,json,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def live(root):
 p=root/'provenance/child_launch_v1.json'
 if not p.exists():return False
 c=load(p);s=Path('/proc')/str(c['PID'])/'stat'
 return s.exists() and s.read_text().rsplit(')',1)[1].split()[19]==c['startticks']

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);pc=load(P/'config/spec_v1.json');R=Path(pc['context_producer']);rc=load(R/'config/spec_v1.json');J=Path(load(Path(rc['TCR_producer'])/'config/spec_v1.json')['joined']);accepted=load(rc['working_acceptance']);w=load(accepted['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};ri={r['source']:r for r in load(R/'outputs/unknown_parent_TCR_context_index_v1.json')};frozen=pc['frozen_preflight'];expected={r['source']:r for r in frozen['sources']}
 for name,digest in load(P/'config/supervisor_request_v1.json')['script_SHA256'].items():assert sha(P/name)==digest
 assert set(wi)==set(expected);sys.path.insert(0,str(J/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection,row_axis_sha256,text_axis_sha256
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);lm=Path(load(J/'outputs/joined_native_atlas_delivery_v1.json')['local_review']);G=Path(load(lm/'config/local_identity_review_spec_v1.json')['graphs']);gcfg=load(G/'config/full_local_graph_spec_v1.json');I=Path(gcfg['local_inputs']);ii={r['source']:r for r in load(I/'outputs/full_native_local_discovery_input_index_v1.json')};params=load(gcfg['original_model_spec']);done={};totals=collections.Counter()
 while set(done)!=set(wi):
  for source in sorted(set(wi)-set(done)):
   op=P/'outputs'/(source+'_common_TNK_reentry_groups_v1.json')
   if not op.exists():continue
   # Metadata is written only after array/model readback; source files are
   # immutable afterwards. A partial JSON write is retried, not accepted.
   try:m=load(op)
   except json.JSONDecodeError:continue
   r=wi[source];mp=base/r['metadata'];arr=base/r['arrays'];assert sha(mp)==r['metadata_SHA256'] and sha(arr)==r['arrays_SHA256'];wm=load(mp)
   with np.load(arr,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','consensus_major_code','annotation_group_code','fit_pool_candidate_not_final_QC','source_qualification_code']}
   rp=R/'outputs'/ri[source]['metadata'];assert sha(rp)==ri[source]['metadata_SHA256'];review=load(rp);major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];known=np.isin(major,['T','NK','T_NK']);control=np.isin(q['annotation_group_code'],[g['previous_group_code'] for g in review['groups']]);pool=np.flatnonzero((known|control)&q['fit_pool_candidate_not_final_QC']);assert (q['source_qualification_code'][pool]==1).all();assert sha(m['arrays'])==m['arrays_SHA256']
   with np.load(m['arrays'],allow_pickle=False) as f:
    for k,v in [('source_scope_position',pool),('source_native_position',q['source_native_position'][pool]),('native_row',q['native_row'][pool]),('baseline_known_T_NK',known[pool]),('unresolved_parent_control',control[pool]),('previous_annotation_group_code',q['annotation_group_code'][pool])]:np.testing.assert_array_equal(f[k],v)
    codes=f['reentry_group_code']
   n=len(pool);ng=len(set(codes.tolist())) if n else 0;assert m['n_fit_pool']==n==expected[source]['n_new_fit_pool'] and m['n_groups']==ng and codes.shape==(n,)
   if n:assert set(codes.tolist())==set(range(ng))
   assert m['n_known_T_NK_fit']==int(known[pool].sum()) and m['n_unknown_parent_control_fit']==int(control[pool].sum()) and m['n_qualified_scope']==len(major);assert not np.any(known&control);np.testing.assert_array_equal(native.rows(source)[q['source_native_position'][pool]],q['native_row'][pool]);assert text_axis_sha256(native.native_ids(source))==wm['native_ids_source_SHA256']
   im=load(I/'outputs'/ii[source]['metadata']);fp=I/'outputs'/im['feature_candidates'];assert sha(fp)==im['feature_candidates_sha256']==expected[source]['candidate_feature_axis_SHA256'];features=load(fp);bycolumn=dict(zip(features['native_columns'],features['genes']));model=m['model']
   if model is not None:
    assert m['status']=='qualified_source_common_TNK_and_unknown_control_graph' and n>=params['minimum_fit_group'] and sha(model['path'])==model['SHA256'] and model['native_candidate_feature_metadata_SHA256']==sha(fp)
    with np.load(model['path'],allow_pickle=False) as f:
     cols=f['native_columns'];embed=f['fit_embedding'];components=f['components'];mean=f['mean'];var=f['explained_variance'];np.testing.assert_array_equal(f['source_scope_position'],pool);np.testing.assert_array_equal(f['source_native_position'],q['source_native_position'][pool]);np.testing.assert_array_equal(f['native_row'],q['native_row'][pool]);np.testing.assert_array_equal(f['reentry_group_code'],codes)
    assert len(cols)==len(set(cols))==len(model['selected_native_genes'])==model['n_selected_features']<=gcfg['n_features'];assert [bycolumn[int(c)] for c in cols]==model['selected_native_genes'];assert not set(cols)&set(features['heldout_physical_alias_columns_excluded']);assert embed.shape==(n,components.shape[0]) and components.shape[1]==len(cols) and mean.shape==(len(cols),) and var.shape==(components.shape[0],);assert all(np.isfinite(v).all() for v in [embed,components,mean,var]) and (var>=0).all();assert model['fit_native_rows_SHA256']==row_axis_sha256(q['native_row'][pool]);assert model['original_feature_axis_and_heldout_alias_exclusion_preserved'] and model['fit_is_current_scope_transductive_not_independent_validation']
   elif n:assert m['status'] in ['small_source_pool_kept_coarse_no_embedding','insufficient_variable_features_kept_coarse_no_embedding'] and ng==1
   else:assert m['status']=='empty_qualified_pool' and ng==0
   assert m['current_working_metadata_SHA256']==sha(mp) and m['current_working_arrays_SHA256']==sha(arr) and m['unknown_parent_context_metadata_SHA256']==sha(rp) and m['all_training_records_current_source_qualified'] and not m['new_identity_assigned']
   done[source]=dict(source=source,metadata=op.name,metadata_SHA256=sha(op));totals.update(n_sources=1,n_qualified_scope=len(major),n_fit_pool=n,n_known_T_NK_fit=int(known[pool].sum()),n_unknown_parent_control_fit=int(control[pool].sum()),n_groups=ng,n_native_nnz_read=m['n_native_nnz_read']);save(root/'outputs'/(source+'_common_TNK_scope_audit_v1.json'),dict(source=source,n_pool=n,n_groups=ng,old_annotation_source_scope_and_model_axes_accepted=True,native_feature_and_PCA_replay_not_yet_accepted=True));native._fields.clear();del q;gc.collect()
  terminal=P/'provenance/supervisor_terminal_v1.json'
  if terminal.exists():
   assert load(terminal)['actual_child_exit_code']==0,'Producer failed, accepted source checkpoints retained'
   if set(done)!=set(wi):raise ValueError('Producer ended without complete source outputs')
  elif not live(P):raise RuntimeError('Exact producer process missing; source checkpoints retained')
  if set(done)!=set(wi):
   if time.monotonic()-start>cfg['maximum_wait_seconds']:raise TimeoutError('Bounded producer wait exhausted')
   time.sleep(15)
 while not (P/'provenance/supervisor_terminal_v1.json').exists():
  assert live(P)
  if time.monotonic()-start>cfg['maximum_wait_seconds']:raise TimeoutError('Bounded terminal wait exhausted')
  time.sleep(2)
 assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;pt=load(P/'provenance/terminal_v1.json');ip=P/'outputs/common_TNK_reentry_index_v1.json';assert sha(ip)==pt['index_SHA256'];actual={r['source']:r for r in load(ip)};assert actual==done and pt['totals']==dict(totals)
 assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];assert totals['n_fit_pool']==frozen['n_new_fit_pool'];save(root/'provenance/terminal_v1.json',dict(status='ALL_CURRENT_COMMON_TNK_SCOPE_AND_MODEL_AXES_ACCEPTED_NATIVE_FEATURE_REPLAY_PENDING',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),producer_index_SHA256=sha(ip),all_source_gate_fit_masks_and_control_memberships_exact=True,old_annotations_and_heldout_alias_rules_preserved=True,native_feature_selection_and_PCA_arithmetic_acceptance_pending=True,new_identity_assigned=False,new_CM_started=False,wall_seconds=time.monotonic()-start))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
