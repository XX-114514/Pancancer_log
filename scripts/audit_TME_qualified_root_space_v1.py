"""Full qualified-root native readback, marker/feature reconstruction and PCA replay."""
import argparse,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
import psutil
import scipy.sparse as sp
from sklearn.decomposition import PCA
from audit_TME_hierarchical_clusters_v1 import independent_markers,equal

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def live(P):
 p=P/'provenance/child_launch_v1.json'
 if not p.exists():return False
 c=load(p);f=Path('/proc')/str(c['PID'])/'stat';return f.exists() and f.read_text().rsplit(')',1)[1].split()[19]==c['startticks']

def main(root):
 begin=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);pc=load(P/'config/spec_v1.json');pre=pc['frozen_preflight'];ac=Path(pc['working_acceptance']);assert sha(ac)==pre['working_acceptance_SHA256'];a=load(ac);assert sha(a['producer_manifest']['path'])==a['producer_manifest']['SHA256'];w=load(a['producer_manifest']['path']);assert sha(w['index']['path'])==w['index']['SHA256'];base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};expected={r['source']:r for r in pre['sources']};EA=Path(pc['native_evidence_audit']);J=Path(pc['joined']);jm=load(J/'outputs/joined_native_atlas_delivery_v1.json');L=Path(jm['local_review']);G=Path(load(L/'config/local_identity_review_spec_v1.json')['graphs']);gcfg=load(G/'config/full_local_graph_spec_v1.json');params=load(gcfg['original_model_spec'])
 for name,digest in load(P/'config/supervisor_request_v1.json')['script_SHA256'].items():assert sha(P/name)==digest
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(EA/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
 from audit_full_TME_native_evidence_v1 import gather
 from atlas_identity_release_keys_v1 import membership_key
 native=NativeArchiveCollection(jm['archive_manifest']);todo=set(expected);done={};totals=collections.Counter();maximum=collections.defaultdict(float);(root/'checkpoints').mkdir()
 while todo:
  progressed=False
  for source in sorted(todo.copy(),key=lambda s:(expected[s]['n_fit_pool'],s)):
   mp=P/'outputs'/(source+'_TME_qualified_root_space_v1.json')
   if not mp.exists():continue
   try:m=load(mp)
   except json.JSONDecodeError:continue
   progressed=True;ex=expected[source];r=wi[source];wp=base/r['metadata'];ap=base/r['arrays'];assert sha(wp)==r['metadata_SHA256']==ex['working_metadata_SHA256']==m['working_metadata_SHA256'];assert sha(ap)==r['arrays_SHA256']==ex['working_arrays_SHA256']==m['working_arrays_SHA256'];wm=load(wp)
   with np.load(ap,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','native_total_RNA','source_qualification_code','fit_pool_candidate_not_final_QC']}
   assert sha(m['arrays']['path'])==m['arrays']['SHA256'];assert sha(m['statistics']['path'])==m['statistics']['SHA256']
   with np.load(m['arrays']['path'],allow_pickle=False) as f:pool=f['source_scope_position'];rows=f['native_row'];positions=f['source_native_position'];codes=f['root_cluster_code']
   np.testing.assert_array_equal(pool,np.flatnonzero(q['fit_pool_candidate_not_final_QC']));np.testing.assert_array_equal(rows,q['native_row'][pool]);np.testing.assert_array_equal(positions,q['source_native_position'][pool]);assert (q['source_qualification_code']==1).all();n=len(pool);assert n==m['n_fit_pool']==ex['n_fit_pool'];assert m['n_unfitted_retained']==len(q['native_row'])-n
   with np.load(m['statistics']['path'],allow_pickle=False) as f:stored_stats={k:f[k] for k in ['n','cp_sum','log_sum','detected']};ecols=f['native_columns']
   ep=Path(ex['native_evidence_metadata_path']);assert sha(ep)==ex['native_evidence_metadata_SHA256']==m['native_evidence_metadata_SHA256'];e=load(ep);np.testing.assert_array_equal(ecols,e['native_columns']);assert m['native_measured_genes']==e['native_measured_gene_names'];genes=m['native_measured_genes'];ng=m['n_roots'];assert np.array_equal(np.unique(codes),np.arange(ng))
   fp=Path(ex['candidate_feature_axis_path']);assert sha(fp)==ex['candidate_feature_axis_SHA256']==m['candidate_feature_metadata_SHA256'];fmeta=load(fp);assert fmeta['feature_selection_not_performed'];fcols=np.asarray(fmeta['native_columns'],np.int64);assert len(fcols)==len(set(fcols));assert not set(fcols)&set(fmeta['heldout_physical_alias_columns_excluded'])
   assert n*len(fcols)*16+ng*len(ecols)*24+pc['RAM_reserve_bytes']<psutil.virtual_memory().available
   csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']];stat=dict(n=np.zeros(ng,np.int64),**{k:np.zeros((ng,len(ecols)),np.float64) for k in ['cp_sum','log_sum','detected']});chunks=[];nnz=0
   save(root/'checkpoints'/(str(time.time_ns())+'.json'),dict(observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage='independent_complete_native_source_started',source=source,n_fit_pool=n))
   for lo in range(0,n,1024):
    p=np.arange(lo,min(n,lo+1024));raw=gather(csr,rows[p]).astype(np.float64).tocsr();raw.sum_duplicates();raw.eliminate_zeros();assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.equal(raw.data,np.floor(raw.data)).all();total=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(total,q['native_total_RNA'][pool[p]]);assert (total>0).all();nnz+=raw.nnz
    x=raw[:,ecols].tocsr();cp=x.multiply((10000/total)[:,None]).tocsr();lg=cp.copy();lg.data=np.log1p(lg.data);u,inv=np.unique(codes[p],return_inverse=True);indicator=sp.coo_matrix((np.ones(len(p)),(inv,np.arange(len(p)))),shape=(len(u),len(p))).tocsr();stat['n'][u]+=np.bincount(inv,minlength=len(u));stat['cp_sum'][u]+=(indicator@cp).toarray();stat['log_sum'][u]+=(indicator@lg).toarray();stat['detected'][u]+=(indicator@(x>0)).toarray()
    fx=raw[:,fcols].astype(np.float32).tocsr();fx=fx.multiply(np.asarray(10000/total,np.float32)[:,None]).tocsr();fx.data=np.log1p(fx.data);chunks.append(fx)
   for k in stat:
    if k in ['n','detected']:np.testing.assert_array_equal(stat[k],stored_stats[k])
    else:np.testing.assert_allclose(stat[k],stored_stats[k],rtol=2e-10,atol=1e-7)
   assert nnz==m['n_native_nnz_read_in_model_pass'];fx=sp.vstack(chunks,format='csr');del chunks;mean=np.asarray(fx.mean(axis=0)).ravel();var=np.maximum(np.asarray(fx.multiply(fx).mean(axis=0)).ravel()-mean*mean,0);valid=np.flatnonzero((var>0)&(np.asarray((fx>0).sum(axis=0)).ravel()>=5));errors={}
   if m['model'] is not None:
    candidate=fx[:,valid];del fx;mu=np.asarray(candidate.mean(axis=0)).ravel();variance=np.maximum(np.asarray(candidate.multiply(candidate).mean(axis=0)).ravel()-mu*mu,0);score=variance/(mu+1e-3);score[np.asarray((candidate>0).sum(axis=0)).ravel()<5]=-np.inf;order=np.lexsort((np.arange(len(valid)),-score))[:min(pc['n_features'],len(valid))];chosen=valid[order];model=m['model'];assert sha(model['path'])==model['SHA256']
    with np.load(model['path'],allow_pickle=False) as f:saved={k:f[k] for k in f.files}
    np.testing.assert_array_equal(saved['native_columns'],fcols[chosen]);assert model['selected_native_genes']==[fmeta['genes'][j] for j in chosen];np.testing.assert_array_equal(saved['source_scope_position'],pool);np.testing.assert_array_equal(saved['native_row'],rows);np.testing.assert_array_equal(saved['root_cluster_code'],codes)
    dense=candidate[:,order].toarray().astype(np.float32,copy=False);del candidate;gc.collect();pca=PCA(n_components=min(params['components'],n-1,len(order)-1),svd_solver='randomized',random_state=pc['seed'],iterated_power=params['pca_iterations'],copy=False);embedding=pca.fit_transform(dense).astype(np.float32,copy=False);del dense
    for key,value,rtol,atol in [('mean',pca.mean_,1e-6,1e-6),('components',pca.components_,1e-5,1e-6),('explained_variance',pca.explained_variance_,2e-5,1e-6),('fit_embedding',embedding,2e-5,2e-4)]:
     errors[key]=float(np.max(np.abs(saved[key].astype(float)-np.asarray(value,float)),initial=0));np.testing.assert_allclose(saved[key],value,rtol=rtol,atol=atol);maximum[key]=max(maximum[key],errors[key])
    del pca,embedding,saved
   else:assert n<params['minimum_fit_group'] or len(valid)<3
   if 'fx' in locals():del fx
   ids=native.native_ids(source);assert text_axis_sha256(ids)==wm['native_ids_source_SHA256'];totalcp=stat['cp_sum'].sum(axis=0);totaldet=stat['detected'].sum(axis=0)
   for code,g in enumerate(m['groups']):
    take=np.flatnonzero(codes==code);assert len(take)==g['n_exact_group']==stat['n'][code];assert membership_key(source,ids[positions[take]].astype(str).tolist())==g['stable_membership_SHA256'];assert hashlib.sha256(np.asarray(rows[take],'<i8').tobytes()).hexdigest()==g['native_rows_SHA256'];other=n-len(take);othercp=np.maximum(totalcp-stat['cp_sum'][code],0);otherdet=totaldet-stat['detected'][code];assert otherdet.min()>=0;expectedmarkers=independent_markers(genes,stat['cp_sum'][code]/len(take),stat['detected'][code]/len(take),othercp/other if other else None,otherdet/other if other else None);equal(expectedmarkers,g['own_full_native_gene_marker_candidates_not_DEG']);assert not g['new_identity_assigned']
   report=dict(source=source,n_current_qualified_records=len(q['native_row']),n_fit_pool=n,n_roots=ng,all_fit_rows_and_unfitted_retention_and_masks_accepted=True,full_native_gene_statistics_and_marker_ranking_independently_rebuilt=True,feature_universe_and_physical_heldout_alias_exclusion_accepted=True,maximum_absolute_PCA_replay_difference=errors,same_PCA_solver_reexecution_not_orthogonal_biology=True,producer_metadata_SHA256=sha(mp));op=root/'outputs'/(source+'_TME_qualified_root_audit_v1.json');save(op,report);done[source]=dict(source=source,metadata=op.name,metadata_SHA256=sha(op));todo.remove(source);totals.update(n_sources=1,n_qualified_scope_records=len(q['native_row']),n_fit_pool=n,n_unfitted_retained=len(q['native_row'])-n,n_roots=ng,n_native_nnz_independently_read=nnz);save(root/'checkpoints'/(str(time.time_ns())+'.json'),dict(observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage='source_accepted',source=source,totals=dict(totals)));native._fields.clear();del q,stat,stored_stats;gc.collect()
   sup=P/'provenance/supervisor_terminal_v1.json'
   if sup.exists():assert load(sup)['actual_child_exit_code']==0
  sup=P/'provenance/supervisor_terminal_v1.json'
  if sup.exists():assert load(sup)['actual_child_exit_code']==0
  else:assert live(P),'Producer exact handle disappeared; do not infer observation timeout'
  if todo and not progressed:time.sleep(5)
 while not (P/'provenance/supervisor_terminal_v1.json').exists():assert live(P);time.sleep(3)
 assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;pt=load(P/'provenance/terminal_v1.json');ip=P/'outputs/TME_qualified_root_space_index_v1.json';assert sha(ip)==pt['index_SHA256'];index=load(ip);assert {r['source'] for r in index}==set(done)
 for r in index:assert r['metadata_SHA256']==load(root/'outputs'/done[r['source']]['metadata'])['producer_metadata_SHA256']
 assert totals['n_fit_pool']==pre['n_fit_pool']==pt['totals']['n_fit_pool'];out=root/'outputs/TME_qualified_root_audit_index_v1.json';save(out,[done[s] for s in sorted(done)]);save(root/'provenance/terminal_v1.json',dict(status='ALL_QUALIFIED_TME_ROOT_NATIVE_FEATURE_PCA_AND_MARKER_NUMERICS_ACCEPTED_NOT_TAXONOMY',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),maximum_absolute_PCA_replay_difference=dict(maximum),producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),index_SHA256=sha(out),production_normalization_feature_selection_and_marker_ranking_helpers_not_imported=True,PCA_same_solver_reexecution=True,graph_and_biological_state_structure_validation_pending=True,no_new_working_identity_or_CM_applied=True,wall_seconds=time.monotonic()-begin,peak_RSS_bytes_Linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
