"""Full native-input rebuild, independent feature ranking and PCA replay.

Re-read every model-pool record through a different CSR reader. PCA is a
deterministic replay of the same installed solver, not an independent solver or
biological validation. Graph-state preservation and identities remain separate.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
import psutil
import scipy.sparse as sp
from sklearn.decomposition import PCA

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
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);SA=Path(cfg['scope_auditor']);pc=load(P/'config/spec_v1.json');R=Path(pc['context_producer']);rc=load(R/'config/spec_v1.json');old=load(Path(rc['TCR_producer'])/'config/spec_v1.json');J=Path(old['joined']);EA=Path(old['native_audit']);w=load(load(rc['working_acceptance'])['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])}
 for name,digest in load(P/'config/supervisor_request_v1.json')['script_SHA256'].items():assert sha(P/name)==digest
 assert load(EA/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 lm=Path(load(J/'outputs/joined_native_atlas_delivery_v1.json')['local_review']);G=Path(load(lm/'config/local_identity_review_spec_v1.json')['graphs']);gcfg=load(G/'config/full_local_graph_spec_v1.json');I=Path(gcfg['local_inputs']);ii={r['source']:r for r in load(I/'outputs/full_native_local_discovery_input_index_v1.json')};params=load(gcfg['original_model_spec']);seed=gcfg['seed'];expected={r['source']:r for r in pc['frozen_preflight']['sources']}
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(EA/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from audit_full_TME_native_evidence_v1 import gather
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);done={};totals=collections.Counter();error_max=collections.defaultdict(float);resource_cfg=dict(old,output_cap_bytes=cfg.get('output_cap_bytes',1073741824))
 (root/'checkpoints').mkdir()
 while set(done)!=set(wi):
  pending=sorted(set(wi)-set(done),key=lambda s:(expected[s]['n_new_fit_pool'],s))
  for source in pending:
   op=P/'outputs'/(source+'_common_TNK_reentry_groups_v1.json');scope=SA/'outputs'/(source+'_common_TNK_scope_audit_v1.json')
   if not op.exists() or not scope.exists():continue
   try:m=load(op);scope_report=load(scope)
   except json.JSONDecodeError:continue
   assert scope_report['old_annotation_source_scope_and_model_axes_accepted'];storage_gate(root,resource_cfg);wr=wi[source];mp=base/wr['metadata'];wa=base/wr['arrays'];assert sha(mp)==wr['metadata_SHA256'] and sha(wa)==wr['arrays_SHA256'];wm=load(mp)
   with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','native_total_RNA','fit_pool_candidate_not_final_QC','source_qualification_code']}
   assert sha(m['arrays'])==m['arrays_SHA256']
   with np.load(m['arrays'],allow_pickle=False) as f:pool=f['source_scope_position'];np.testing.assert_array_equal(f['native_row'],q['native_row'][pool]);np.testing.assert_array_equal(f['source_native_position'],q['source_native_position'][pool])
   n=len(pool);assert n==m['n_fit_pool']==expected[source]['n_new_fit_pool'];assert q['fit_pool_candidate_not_final_QC'][pool].all() and (q['source_qualification_code'][pool]==1).all()
   im=load(I/'outputs'/ii[source]['metadata']);fp=I/'outputs'/im['feature_candidates'];assert sha(fp)==im['feature_candidates_sha256']==expected[source]['candidate_feature_axis_SHA256'];features=load(fp);cols=np.asarray(features['native_columns'],np.int64);genes=features['genes'];assert len(cols)==len(set(cols))==len(set(genes));assert not set(cols)&set(features['heldout_physical_alias_columns_excluded'])
   bound=n*len(cols)*16+n*min(gcfg['n_features'],len(cols))*12+cfg['RAM_reserve_bytes']
   if bound>psutil.virtual_memory().available:raise MemoryError('Fresh full native feature replay RAM gate failed')
   save(root/'checkpoints'/(str(time.time_ns())+'.json'),dict(observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage='independent_native_feature_source_started',source=source,n_pool=n,n_complete_sources=len(done)))
   chunks=[];nnz=0;csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']]
   for lo in range(0,n,1024):
    sel=pool[lo:lo+1024];raw=gather(csr,q['native_row'][sel]);total=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(total,q['native_total_RNA'][sel]);assert (total>0).all();nnz+=raw.nnz
    # Independent CSR reader coalesces repeated slots before log. Keep the
    # frozen float32 model arithmetic, with complete native RNA denominators.
    x=raw[:,cols].astype(np.float32).tocsr();scale=np.asarray(10000.0/total,np.float32);x=x.multiply(scale[:,None]).tocsr();x.data=np.log1p(x.data);chunks.append(x)
   x=sp.vstack(chunks,format='csr') if chunks else sp.csr_matrix((0,len(cols)),dtype=np.float32);del chunks
   if n>=params['minimum_fit_group']:assert nnz==m['n_native_nnz_read']
   mean=np.asarray(x.mean(axis=0)).ravel() if n else np.zeros(len(cols),np.float32);variance=np.maximum(np.asarray(x.multiply(x).mean(axis=0)).ravel()-mean*mean,0) if n else mean.copy();variable=np.flatnonzero((variance>0)&(np.asarray((x>0).sum(axis=0)).ravel()>=5));model=m['model'];errors={}
   if model is not None:
    assert n>=params['minimum_fit_group'] and len(variable)>=3 and sha(model['path'])==model['SHA256'];candidate=x[:,variable];del x
    cm=np.asarray(candidate.mean(axis=0)).ravel();cv=np.maximum(np.asarray(candidate.multiply(candidate).mean(axis=0)).ravel()-cm*cm,0);score=cv/(cm+1e-3);score[np.asarray((candidate>0).sum(axis=0)).ravel()<5]=-np.inf
    # One-source average ranks reduce to descending dispersion and original
    # variable-column order for ties; no production selection function used.
    order=np.lexsort((np.arange(len(variable)),-score))[:min(gcfg['n_features'],len(variable))];chosen=variable[order];dense=candidate[:,order].toarray().astype(np.float32,copy=False);del candidate;gc.collect()
    with np.load(model['path'],allow_pickle=False) as f:stored={k:f[k] for k in ['native_columns','mean','components','explained_variance','fit_embedding','native_row','source_scope_position']}
    np.testing.assert_array_equal(stored['native_columns'],cols[chosen]);assert model['selected_native_genes']==[genes[j] for j in chosen] and model['n_variable_candidates']==len(variable);np.testing.assert_array_equal(stored['native_row'],q['native_row'][pool]);np.testing.assert_array_equal(stored['source_scope_position'],pool)
    components=min(params['components'],n-1,len(order)-1);replay=PCA(n_components=components,svd_solver='randomized',random_state=seed,iterated_power=params['pca_iterations'],copy=False);embed=replay.fit_transform(dense).astype(np.float32,copy=False);del dense;gc.collect()
    for name,value,rtol,atol in [('mean',replay.mean_,1e-6,1e-6),('components',replay.components_,1e-5,1e-6),('explained_variance',replay.explained_variance_,2e-5,1e-6),('fit_embedding',embed,2e-5,2e-4)]:
     errors[name]=float(np.max(np.abs(stored[name].astype(np.float64)-np.asarray(value,np.float64)),initial=0));np.testing.assert_allclose(stored[name],value,rtol=rtol,atol=atol);error_max[name]=max(error_max[name],errors[name])
    gram=np.asarray(replay.components_,np.float64)@np.asarray(replay.components_,np.float64).T;np.testing.assert_allclose(gram,np.eye(components),rtol=0,atol=1e-4)
    del replay,embed,stored
   elif n and n<params['minimum_fit_group']:assert m['status']=='small_source_pool_kept_coarse_no_embedding'
   elif n:assert len(variable)<3 and m['status']=='insufficient_variable_features_kept_coarse_no_embedding'
   else:assert m['status']=='empty_qualified_pool'
   if 'x' in locals():del x
   report=dict(source=source,n_model_pool_records_directly_reread=n,n_complete_native_nnz_independently_read=nnz,n_source_measured_candidate_features=len(cols),n_native_variable_candidates=len(variable),all_native_totals_and_heldout_alias_masks_checked=True,feature_selection_rebuilt_without_production_helper=True,deterministic_same_solver_PCA_replay_checked=model is not None,maximum_absolute_PCA_replay_difference=errors,producer_source_metadata_SHA256=sha(op),not_biological_identity_or_graph_state_validation=True)
   om=root/'outputs'/(source+'_common_TNK_native_feature_audit_v1.json');save(om,report);done[source]=dict(source=source,metadata=om.name,metadata_SHA256=sha(om),producer_metadata_SHA256=sha(op));totals.update(n_sources=1,n_model_pool_records_directly_reread=n,n_complete_native_nnz_independently_read=nnz,n_sources_with_PCA_replay=int(model is not None));save(root/'checkpoints'/(str(time.time_ns())+'.json'),dict(observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage='independent_native_feature_source_accepted',source=source,totals=dict(totals)));native._fields.clear();del q;gc.collect()
  for stage in [P,SA]:
   terminal=stage/'provenance/supervisor_terminal_v1.json'
   if terminal.exists():assert load(terminal)['actual_child_exit_code']==0,'Upstream failed; completed source audits retained'
   elif not live(stage):raise RuntimeError('Exact upstream process missing; completed source audits retained')
  if set(done)!=set(wi):
   if time.monotonic()-start>cfg['maximum_wait_seconds']:raise TimeoutError('Bounded full-source audit wait exhausted')
   time.sleep(15)
 for stage in [P,SA]:
  while not (stage/'provenance/supervisor_terminal_v1.json').exists():
   assert live(stage)
   if time.monotonic()-start>cfg['maximum_wait_seconds']:raise TimeoutError('Bounded terminal wait exhausted')
   time.sleep(2)
  assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 pt=load(P/'provenance/terminal_v1.json');ip=P/'outputs/common_TNK_reentry_index_v1.json';assert sha(ip)==pt['index_SHA256'];records=load(ip);assert set(done)=={r['source'] for r in records}
 for r in records:assert done[r['source']]['producer_metadata_SHA256']==r['metadata_SHA256']
 assert totals['n_model_pool_records_directly_reread']==pt['totals']['n_fit_pool']==pc['frozen_preflight']['n_new_fit_pool'];assert load(SA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json');out=root/'outputs/common_TNK_native_feature_audit_index_v1.json';save(out,[done[s] for s in sorted(done)])
 save(root/'provenance/terminal_v1.json',dict(status='ALL_COMMON_TNK_NATIVE_INPUT_FEATURE_SELECTION_AND_DETERMINISTIC_PCA_REPLAY_ACCEPTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),maximum_absolute_PCA_replay_difference=dict(error_max),producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),scope_audit_terminal_SHA256=sha(SA/'provenance/terminal_v1.json'),index_SHA256=sha(out),production_normalization_and_feature_selection_helpers_not_imported=True,PCA_same_solver_replay_not_independent_solver_validation=True,graph_state_preservation_and_biological_identity_review_pending=True,new_identity_assigned=False,new_CM_started=False,wall_seconds=time.monotonic()-start,peak_RSS_bytes_Linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
