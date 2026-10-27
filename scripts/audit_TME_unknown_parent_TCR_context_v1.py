"""Independent membership and exact-neighbor distance audit for unknown parents."""
import argparse,collections,datetime,gc,hashlib,json,resource,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')

def main(root):
 start=time.monotonic();P=Path(load(root/'config/spec_v1.json')['producer']);assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 cfg=load(P/'config/spec_v1.json');K=Path(cfg['TCR_producer']);pc=load(K/'config/spec_v1.json');J=Path(pc['joined']);w=load(load(cfg['working_acceptance'])['producer_manifest']['path']);base=Path(w['index']['path']).parent
 for name,digest in load(P/'config/supervisor_request_v1.json')['script_SHA256'].items():assert sha(P/name)==digest
 assert sha(cfg['working_acceptance'])==cfg['working_acceptance_SHA256'];wi={r['source']:r for r in load(w['index']['path'])};ki={r['source']:r for r in load(K/'outputs/TCR_boundary_proposal_index_v1.json')};catalog=load(w['group_catalog']['path'])['groups']
 lmroot=Path(load(J/'outputs/joined_native_atlas_delivery_v1.json')['local_review']);lmconfig=load(lmroot/'config/local_identity_review_spec_v1.json');graphs=Path(lmconfig['graphs']);li={r['source']:r for r in load(lmroot/'outputs/full_local_identity_review_index_v1.json')};gcfg=load(graphs/'config/full_local_graph_spec_v1.json');kreq=int(load(gcfg['original_model_spec'])['neighbors'])
 for b in gcfg['bindings'].values():assert sha(b['path'])==b['sha256']
 ip=P/'outputs/unknown_parent_TCR_context_index_v1.json';assert sha(ip)==load(P/'provenance/terminal_v1.json')['index_SHA256'];index=load(ip);assert {r['source'] for r in index}==set(wi)
 totals=collections.Counter();reasons=collections.Counter();routes=collections.Counter();n_distances=0;n_exposure=0
 for entry in index:
  source=entry['source'];op=P/'outputs'/entry['metadata'];assert sha(op)==entry['metadata_SHA256'];out=load(op);assert sha(out['arrays'])==out['arrays_SHA256']
  with np.load(out['arrays'],allow_pickle=False) as f:stored={k:f[k] for k in f.files}
  wr=wi[source];mp=base/wr['metadata'];wa=base/wr['arrays'];assert sha(mp)==wr['metadata_SHA256'] and sha(wa)==wr['arrays_SHA256'];wm=load(mp)
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','consensus_major_code','parent_source_local_group_code','annotation_group_code','fit_pool_candidate_not_final_QC']}
  kp=K/'outputs'/ki[source]['metadata'];assert sha(kp)==ki[source]['metadata_SHA256'];km=load(kp);ka=K/'outputs'/km['arrays'];assert sha(ka)==km['arrays_SHA256']
  with np.load(ka,allow_pickle=False) as f:t={k:f[k] for k in ['native_row','review_code','TCR_core','competing_core_mask','reference_low_T']}
  np.testing.assert_array_equal(q['native_row'],t['native_row']);lp=lmroot/'outputs'/li[source]['output'];assert sha(lp)==li[source]['output_sha256'];local=load(lp);byparent={g['source_local_group_code']:g for g in local['groups']};unknown_codes=[k for k,g in byparent.items() if g['local_major_proposal']=='unresolved'];target=np.isin(q['parent_source_local_group_code'],unknown_codes)&(t['review_code']==7)
  major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];assert (major[target]=='unresolved').all() and q['fit_pool_candidate_not_final_QC'][target].all();np.testing.assert_array_equal(stored['target_TCR_unknown_parent'],target)
  for name in ['native_row','source_native_position']:np.testing.assert_array_equal(stored[name],q[name])
  gr=local['source_graph_record'];gp=graphs/'outputs'/gr['metadata'];ga=graphs/'outputs'/gr['arrays'];assert sha(gp)==gr['metadata_sha256'] and sha(ga)==gr['arrays_sha256'];gm=load(gp);actual_branches={r['branch_id']:r for r in out['branches']};needed=set()
  with np.load(ga,allow_pickle=False) as f:n_archive=len(f['native_row']);np.testing.assert_array_equal(q['native_row'],f['native_row'][q['source_native_position']])
  lookup=np.full(n_archive,-1,np.int64);lookup[q['source_native_position']]=np.arange(len(major));covered=np.zeros(len(major),bool);fractions=np.full((len(major),4),np.nan)
  for branch in gm['branches']:
   if 'model' not in branch:continue
   model=graphs/branch['model'];assert sha(model)==branch['model_sha256']
   with np.load(model,allow_pickle=False) as f:original=f['source_row_positions'];in_scope=lookup[original]>=0;pool=lookup[original[in_scope]];query=np.flatnonzero(target[pool])
   if not len(query):continue
   needed.add(branch['branch_id']);b=actual_branches[branch['branch_id']];assert b['model_path']==str(model) and b['model_SHA256']==sha(model) and b['model_metadata_SHA256']==sha(graphs/branch['model_metadata']);assert sha(b['neighbors_path'])==b['neighbors_SHA256']
   n=len(pool);k=min(kreq,n-1);assert b['n_original_embedding_fit_records']==len(original) and b['n_qualified_neighbor_pool']==n and b['n_fit_records_outside_current_source_gate']==int((~in_scope).sum()) and b['n_queries']==len(query) and b['requested_k']==kreq and b['effective_k']==k
   assert b['n_known_T_in_qualified_pool']==int((major[pool]=='T').sum()) and b['original_embedding_was_fitted_before_current_source_gate'];n_exposure+=int((~in_scope).sum())
   with np.load(model,allow_pickle=False) as f:X=f['fit_embedding'][in_scope].astype(np.float64)
   with np.load(b['neighbors_path'],allow_pickle=False) as f:queries=f['query_scope_position'];neighbors=f['neighbor_scope_position'];dist=f['neighbor_euclidean_distance'];np.testing.assert_array_equal(pool,f['pool_scope_position'])
   np.testing.assert_array_equal(queries,pool[query]);assert neighbors.shape==dist.shape==(len(query),k);inv=np.full(len(major),-1,np.int64);inv[pool]=np.arange(n);ni=inv[neighbors];assert (ni>=0).all()
   assert not np.any(neighbors==queries[:,None]);assert all(len(set(row))==k for row in ni);assert np.isfinite(dist).all() and (dist>=0).all() and (np.diff(dist,axis=1)>=-1e-9).all()
   norms=np.einsum('ij,ij->i',X,X)
   for lo in range(0,len(query),64):
    take=query[lo:lo+64];selected=ni[lo:lo+64];delta=X[take,None,:]-X[selected];actual=np.sqrt(np.einsum('ijk,ijk->ij',delta,delta));np.testing.assert_allclose(dist[lo:lo+64],actual,rtol=1e-7,atol=1e-7)
    full=norms[take,None]+norms[None,:]-2*(X[take]@X.T);np.maximum(full,0,out=full);full[np.arange(len(take)),take]=np.inf;kth=np.partition(full,k-1,axis=1)[:,k-1];bound=np.max(actual**2,axis=1);np.testing.assert_allclose(bound,kth,rtol=1e-7,atol=1e-7);n_distances+=len(take)*n
   covered[queries]=True
   fractions[queries]=np.asarray([[np.count_nonzero(major[row]=='T')/k,np.count_nonzero(major[row]=='unresolved')/k,np.count_nonzero(target[row])/k,np.count_nonzero(t['TCR_core'][row])/k] for row in neighbors])
   del X;gc.collect()
  assert set(actual_branches)==needed;np.testing.assert_array_equal(stored['neighbor_evaluable'],covered);np.testing.assert_allclose(stored['neighbor_fraction_values'],fractions,rtol=0,atol=0,equal_nan=True)
  current={g['source_annotation_group_code']:g for g in catalog if g['source']==source};expected=[]
  for code in sorted(set(map(int,q['annotation_group_code'][target]))):
   g=current[code];sel=q['annotation_group_code']==code;take=sel&target;ok=take&covered;parents=set(map(int,q['parent_source_local_group_code'][sel]));assert len(parents)==1;parent=byparent[parents.pop()]
   r=dict(previous_group_code=code,previous_group_number=g['global_subcluster_number'],previous_group_SHA256=g['stable_group_sha256'],n_scope_group=int(sel.sum()),n_TCR_candidate=int(take.sum()),n_TCR_core_in_parent=int((sel&t['TCR_core']).sum()),n_TCR_and_competing_core_in_parent=int((sel&t['TCR_core']&(t['competing_core_mask']!=0)).sum()),n_unique_low_T_in_parent=int((sel&t['reference_low_T']).sum()),local_parent_group=parent['source_local_group_code'],local_parent_major=parent['local_major_proposal'],local_parent_reason=parent['proposal_evidence_status'],local_parent_provisional_route=parent['provisional_routing_family'],local_parent_native_supported_families=parent['native_supported_families_in_legacy_rule'],n_neighbor_evaluable=int(ok.sum()),n_no_frozen_embedding_or_comparable_neighbors=int((take&~covered).sum()))
   for j,name in enumerate(['known_T','unresolved','same_TCR_candidate','TCR_core']):r[name+'_neighbor_fraction_mean']=float(fractions[ok,j].mean()) if ok.any() else None
   expected.append(r);reasons[parent['proposal_evidence_status']]+=int(take.sum());routes[parent['provisional_routing_family']]+=int(take.sum())
  assert out['groups']==expected and out['n_scope_records']==len(major) and out['n_candidates']==int(target.sum()) and out['n_groups']==len(expected)
  for k,p in [('working_metadata_SHA256',mp),('working_arrays_SHA256',wa),('TCR_metadata_SHA256',kp),('local_review_SHA256',lp),('source_graph_metadata_SHA256',gp)]:assert out[k]==sha(p)
  totals.update(n_sources=1,n_scope_records=len(major),n_affected_sources=int(target.any()),n_candidate_records=int(target.sum()),n_exact_groups=len(expected),n_parent_scope_records=sum(r['n_scope_group'] for r in expected),n_neighbor_evaluable=int(covered.sum()),n_unavailable_neighbors=int((target&~covered).sum()));save(root/'outputs'/(source+'_context_audit_v1.json'),dict(source=source,all_memberships_full_distance_search_and_scope_rebuilt=True));del q,t,stored;gc.collect()
 summary=load(P/'outputs/unknown_parent_TCR_context_summary_v1.json');assert dict(totals)==summary['totals']==load(P/'provenance/terminal_v1.json')['totals'];assert dict(reasons)==summary['candidate_counts_by_original_parent_reason'] and dict(routes)==summary['candidate_counts_by_original_provisional_route']
 save(root/'provenance/terminal_v1.json',dict(status='ALL_UNKNOWN_PARENT_TCR_CONTEXT_MEMBERSHIP_AND_EXACT_NEIGHBORS_ACCEPTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),n_full_query_pool_distances_independently_checked=n_distances,n_branch_fit_records_outside_current_scope=n_exposure,producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),producer_index_SHA256=sha(ip),new_production_review_or_neighbor_search_library_not_imported=True,not_orthogonal_identity_accuracy=True,new_identity_assigned=False,wall_seconds=time.monotonic()-start,peak_RSS_bytes_Linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
