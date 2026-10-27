"""Reassess unknown local parents without turning a missing label into conflict.

Rebuild exact, source/branch-conditional neighborhoods from frozen embeddings.
Only qualified source-view cells enter neighbor pools. The original embedding
training exposure is reported separately, never called a new TME-only fit.
No identities, state values, source eligibility, or count matrices are changed.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
import psutil
from sklearn import config_context
from sklearn.neighbors import NearestNeighbors

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def hist(a):return dict(collections.Counter(map(str,a)))

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');K=Path(cfg['TCR_producer']);KA=Path(cfg['TCR_audit']);KD=Path(cfg['TCR_delivery']);pc=load(K/'config/spec_v1.json');J=Path(pc['joined'])
 for s in [K,KA,KD]:assert load(s/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(KA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(K/'provenance/terminal_v1.json')
 assert sha(cfg['working_acceptance'])==cfg['working_acceptance_SHA256'];accepted=load(cfg['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];w=load(accepted['producer_manifest']['path'])
 for field in ['index','group_catalog']:assert sha(w[field]['path'])==w[field]['SHA256']
 wi=load(w['index']['path']);base=Path(w['index']['path']).parent;catalog=load(w['group_catalog']['path'])['groups'];ki={r['source']:r for r in load(K/'outputs/TCR_boundary_proposal_index_v1.json')}
 assert sha(K/'outputs/TCR_boundary_proposal_index_v1.json')==load(K/'provenance/terminal_v1.json')['index_SHA256']
 lmroot=Path(load(J/'outputs/joined_native_atlas_delivery_v1.json')['local_review']);lmconfig=load(lmroot/'config/local_identity_review_spec_v1.json');graphs=Path(lmconfig['graphs']);li={r['source']:r for r in load(lmroot/'outputs/full_local_identity_review_index_v1.json')};gcfg=load(graphs/'config/full_local_graph_spec_v1.json');neighbors=int(load(gcfg['original_model_spec'])['neighbors']);assert neighbors>0
 sys.path.insert(0,str(J/'scripts'))
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 resource_cfg=dict(pc,output_cap_bytes=cfg.get('output_cap_bytes',2147483648));index=[];totals=collections.Counter();reasons=collections.Counter();routes=collections.Counter()
 for wr in wi:
  storage_gate(root,resource_cfg);source=wr['source'];mp=base/wr['metadata'];arr=base/wr['arrays'];assert sha(mp)==wr['metadata_SHA256'] and sha(arr)==wr['arrays_SHA256'];wm=load(mp)
  with np.load(arr,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','consensus_major_code','annotation_group_code','parent_source_local_group_code','fit_pool_candidate_not_final_QC']}
  kr=ki[source];kp=K/'outputs'/kr['metadata'];assert sha(kp)==kr['metadata_SHA256'];km=load(kp);ka=K/'outputs'/km['arrays'];assert sha(ka)==km['arrays_SHA256']
  with np.load(ka,allow_pickle=False) as f:t={k:f[k] for k in ['native_row','review_code','TCR_core','competing_core_mask','reference_low_T']}
  np.testing.assert_array_equal(q['native_row'],t['native_row']);lr=li[source];lp=lmroot/'outputs'/lr['output'];assert sha(lp)==lr['output_sha256'];local=load(lp);localby={g['source_local_group_code']:g for g in local['groups']}
  parent_unknown=np.asarray([localby[int(code)]['local_major_proposal']=='unresolved' if code>=0 else False for code in q['parent_source_local_group_code']],bool);major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']]
  target=(t['review_code']==7)&parent_unknown;assert (major[target]=='unresolved').all() and q['fit_pool_candidate_not_final_QC'][target].all() and (t['competing_core_mask'][target]==0).all()
  gr=local['source_graph_record'];gp=graphs/'outputs'/gr['metadata'];ga=graphs/'outputs'/gr['arrays'];assert sha(gp)==gr['metadata_sha256'] and sha(ga)==gr['arrays_sha256'];gm=load(gp)
  with np.load(ga,allow_pickle=False) as f:
   graph_rows=f['native_row'];np.testing.assert_array_equal(q['native_row'],graph_rows[q['source_native_position']]);np.testing.assert_array_equal(q['parent_source_local_group_code'],f['source_local_group_code'][q['source_native_position']])
  lookup=np.full(len(graph_rows),-1,np.int64);lookup[q['source_native_position']]=np.arange(len(major));per=np.full((len(major),4),np.nan);covered=np.zeros(len(major),bool);branch_reports=[]
  if target.any():
   for branch in gm['branches']:
    if 'model' not in branch:continue
    model=graphs/branch['model'];assert sha(model)==branch['model_sha256'];mm=graphs/branch['model_metadata'];assert sha(mm)==branch['model_metadata_sha256']
    with np.load(model,allow_pickle=False) as f:positions=f['source_row_positions'];keep=lookup[positions]>=0;pool=lookup[positions[keep]];query=np.flatnonzero(target[pool])
    if not len(query):continue
    n=len(pool);k=min(neighbors,n-1);assert k>=1
    if n*int(gcfg.get('n_features',2000))*16+2147483648>psutil.virtual_memory().available:raise MemoryError('Fresh conditional-neighbor RAM reserve failed')
    with np.load(model,allow_pickle=False) as f:X=f['fit_embedding'][keep].astype(np.float64);np.testing.assert_array_equal(f['native_rows'][keep],q['native_row'][pool])
    assert np.isfinite(X).all()
    with config_context(working_memory=256):
     distance,near=NearestNeighbors(n_neighbors=k+1,algorithm='brute',metric='euclidean',n_jobs=2).fit(X).kneighbors(X[query])
    ni=np.empty((len(query),k),np.int64);nd=np.empty((len(query),k),np.float64)
    for i,me in enumerate(query):
     sel=near[i]!=me;ni[i]=near[i][sel][:k];nd[i]=distance[i][sel][:k];assert len(set(ni[i]))==k and me not in ni[i]
    scope=pool[query];neighbors_scope=pool[ni];covered[scope]=True;per[scope]=np.column_stack([(major[neighbors_scope]=='T').mean(axis=1),(major[neighbors_scope]=='unresolved').mean(axis=1),target[neighbors_scope].mean(axis=1),t['TCR_core'][neighbors_scope].mean(axis=1)])
    op=root/'outputs'/(source+'_'+branch['branch_id']+'_conditional_neighbors_v1.npz')
    with op.open('xb') as f:np.savez_compressed(f,query_scope_position=scope,neighbor_scope_position=neighbors_scope,neighbor_euclidean_distance=nd,pool_scope_position=pool)
    branch_reports.append(dict(branch_id=branch['branch_id'],provisional_family=branch['provisional_family'],model_path=str(model),model_SHA256=sha(model),model_metadata_SHA256=sha(mm),n_original_embedding_fit_records=len(positions),n_qualified_neighbor_pool=n,n_fit_records_outside_current_source_gate=int((~keep).sum()),n_queries=len(query),requested_k=neighbors,effective_k=k,neighbors_path=str(op),neighbors_SHA256=sha(op),n_known_T_in_qualified_pool=int((major[pool]=='T').sum()),original_embedding_was_fitted_before_current_source_gate=True))
    del X,near,distance;gc.collect()
  groups=[];bycode={g['source_annotation_group_code']:g for g in catalog if g['source']==source}
  for code in np.unique(q['annotation_group_code'][target]):
   sel=q['annotation_group_code']==code;take=sel&target;g=bycode[int(code)];parents=set(map(int,q['parent_source_local_group_code'][sel]));assert len(parents)==1;parent=localby[parents.pop()];assert parent['local_major_proposal']=='unresolved';ok=take&covered
   r=dict(previous_group_code=int(code),previous_group_number=g['global_subcluster_number'],previous_group_SHA256=g['stable_group_sha256'],n_scope_group=int(sel.sum()),n_TCR_candidate=int(take.sum()),n_TCR_core_in_parent=int((sel&t['TCR_core']).sum()),n_TCR_and_competing_core_in_parent=int((sel&t['TCR_core']&(t['competing_core_mask']!=0)).sum()),n_unique_low_T_in_parent=int((sel&t['reference_low_T']).sum()),local_parent_group=parent['source_local_group_code'],local_parent_major=parent['local_major_proposal'],local_parent_reason=parent['proposal_evidence_status'],local_parent_provisional_route=parent['provisional_routing_family'],local_parent_native_supported_families=parent['native_supported_families_in_legacy_rule'],n_neighbor_evaluable=int(ok.sum()),n_no_frozen_embedding_or_comparable_neighbors=int((take&~covered).sum()))
   for j,name in enumerate(['known_T','unresolved','same_TCR_candidate','TCR_core']):r[name+'_neighbor_fraction_mean']=float(per[ok,j].mean()) if ok.any() else None
   groups.append(r);reasons[parent['proposal_evidence_status']]+=int(take.sum());routes[parent['provisional_routing_family']]+=int(take.sum())
  sp=root/'outputs'/(source+'_unknown_parent_TCR_review_v1.npz')
  with sp.open('xb') as f:np.savez_compressed(f,native_row=q['native_row'],source_native_position=q['source_native_position'],target_TCR_unknown_parent=target,neighbor_evaluable=covered,neighbor_fraction_values=per)
  om=root/'outputs'/(source+'_unknown_parent_TCR_review_v1.json');out=dict(source=source,n_scope_records=len(major),n_candidates=int(target.sum()),n_groups=len(groups),groups=groups,branches=branch_reports,arrays=str(sp),arrays_SHA256=sha(sp),working_metadata_SHA256=sha(mp),working_arrays_SHA256=sha(arr),TCR_metadata_SHA256=sha(kp),local_review_SHA256=sha(lp),source_graph_metadata_SHA256=sha(gp),n_unknown_parent_groups_with_other_coherent_cores=sum(r['n_TCR_and_competing_core_in_parent']>0 for r in groups),new_identity_assigned=False,neighbor_context_not_orthogonal_truth=True);save(om,out);index.append(dict(source=source,metadata=om.name,metadata_SHA256=sha(om)))
  totals.update(n_sources=1,n_scope_records=len(major),n_affected_sources=int(target.any()),n_candidate_records=int(target.sum()),n_exact_groups=len(groups),n_parent_scope_records=sum(g['n_scope_group'] for g in groups),n_neighbor_evaluable=int(covered.sum()),n_unavailable_neighbors=int((target&~covered).sum()));save(root/'outputs'/(source+'_context_checkpoint_v1.json'),dict(completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals)));del q,t,per;gc.collect()
 expected=load(KD/'outputs/TCR_proposal_delivery_summary_v1.json')['local_group_discordance_by_parent_local_major']['unresolved'];assert totals['n_candidate_records']==expected
 ip=root/'outputs/unknown_parent_TCR_context_index_v1.json';save(ip,index);summary=dict(totals=dict(totals),candidate_counts_by_original_parent_reason=dict(reasons),candidate_counts_by_original_provisional_route=dict(routes),all_queries_and_neighbor_pools_restricted_to_current_source_gate=True,original_embedding_training_exposure_separate=True,not_original_graph_adjacency_replay=True,new_identity_or_model_or_CM_assigned=False);op=root/'outputs/unknown_parent_TCR_context_summary_v1.json';save(op,summary)
 save(root/'provenance/terminal_v1.json',dict(status='ALL_UNKNOWN_PARENT_TCR_CONDITIONAL_CONTEXT_REVIEW_REQUIRES_INDEPENDENT_AUDIT',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),index_SHA256=sha(ip),summary_SHA256=sha(op),previous_working_contract_SHA256=sha(cfg['working_acceptance']),new_identity_assigned=False,new_CM_started=False,wall_seconds=time.monotonic()-start,peak_RSS_bytes_Linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
