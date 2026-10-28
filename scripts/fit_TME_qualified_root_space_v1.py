"""Rebuild source-local roots on every qualified fit record, before lineage routing.

Retain original counts, author/consensus labels, QC flags and RNA states. The
feature universe is native/measurement-defined with frozen physical held-outs;
variance selection is learned again only in this exact current fit pool.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
import psutil
import scipy.sparse as sp
from TME_native_measurement_core_v1 import canonical_native,update_group_stats,expression_candidates

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def check(root,stage,**v):save(root/'checkpoints'/(str(time.time_ns())+'.json'),dict(observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage=stage,**v))

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');pre=cfg['frozen_preflight'];ac=Path(cfg['working_acceptance']);assert sha(ac)==pre['working_acceptance_SHA256'];a=load(ac);assert sha(a['producer_manifest']['path'])==a['producer_manifest']['SHA256'];w=load(a['producer_manifest']['path']);base=Path(w['index']['path']).parent;assert sha(w['index']['path'])==w['index']['SHA256'];J=Path(cfg['joined']);E=Path(cfg['native_evidence']);EA=Path(cfg['native_evidence_audit']);assert load(EA/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 jm=load(J/'outputs/joined_native_atlas_delivery_v1.json');L=Path(jm['local_review']);G=Path(load(L/'config/local_identity_review_spec_v1.json')['graphs']);gcfg=load(G/'config/full_local_graph_spec_v1.json');params=dict(load(gcfg['original_model_spec']));params['threads']=cfg['threads'];sys.path.insert(0,str(G/'scripts'));sys.path.insert(0,str(J/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
 from native_panel_aggregate_v1 import gather_csr_rows
 from annotation_representation_pilot_v1 import normalized_log_counts,select_features,pca_fit
 from compare_annotation_granularity_v1 import graph_grid
 from atlas_identity_release_keys_v1 import membership_key
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 native=NativeArchiveCollection(jm['archive_manifest']);expected={r['source']:r for r in pre['sources']};results=[];totals=collections.Counter();(root/'models').mkdir();(root/'checkpoints').mkdir();resource_cfg=dict(load(cfg['resource_policy_source']),output_cap_bytes=cfg['output_cap_bytes'])
 for item in sorted(load(w['index']['path']),key=lambda r:(expected[r['source']]['n_fit_pool'],r['source'])):
  source=item['source'];storage_gate(root,resource_cfg);ex=expected[source];mp=base/item['metadata'];ap=base/item['arrays'];assert sha(mp)==item['metadata_SHA256']==ex['working_metadata_SHA256'] and sha(ap)==item['arrays_SHA256']==ex['working_arrays_SHA256'];wm=load(mp)
  with np.load(ap,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','native_total_RNA','source_qualification_code','fit_pool_candidate_not_final_QC']}
  pool=np.flatnonzero(q['fit_pool_candidate_not_final_QC']);n=len(pool);assert n==ex['n_fit_pool'];assert (q['source_qualification_code']==1).all() and (q['native_total_RNA'][pool]>0).all();fp=Path(ex['candidate_feature_axis_path']);assert sha(fp)==ex['candidate_feature_axis_SHA256'];feat=load(fp);assert feat['feature_selection_not_performed'];cols=np.asarray(feat['native_columns'],np.int64);genes=feat['genes'];assert len(cols)==len(set(cols))==len(set(genes));assert not set(cols)&set(feat['heldout_physical_alias_columns_excluded'])
  ep=Path(ex['native_evidence_metadata_path']);assert sha(ep)==ex['native_evidence_metadata_SHA256'];e=load(ep);ecols=np.asarray(e['native_columns'],np.int64);egenes=e['native_measured_gene_names'];assert len(ecols)==len(set(ecols))==len(set(egenes));bound=n*len(cols)*16+n*min(cfg['n_features'],len(cols))*12+cfg['RAM_reserve_bytes'];assert bound<psutil.virtual_memory().available
  ids=native.native_ids(source);assert text_axis_sha256(ids)==wm['native_ids_source_SHA256'];np.testing.assert_array_equal(native.rows(source)[q['source_native_position'][pool]],q['native_row'][pool]);csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']];chunks=[];nnz=0;check(root,'source_native_model_input_started',source=source,n_fit_pool=n)
  for lo in range(0,n,cfg['read_block_size']):
   sel=pool[lo:lo+cfg['read_block_size']];raw=canonical_native(gather_csr_rows(csr,q['native_row'][sel],max_span=8192));total=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(total,q['native_total_RNA'][sel]);x,positive=normalized_log_counts(raw[:,cols],total);assert positive.all();chunks.append(x);nnz+=raw.nnz
  x=sp.vstack(chunks,format='csr');del chunks,raw;mean=np.asarray(x.mean(axis=0)).ravel();var=np.maximum(np.asarray(x.multiply(x).mean(axis=0)).ravel()-mean*mean,0);variable=np.flatnonzero((var>0)&(np.asarray((x>0).sum(axis=0)).ravel()>=5));model=None;codes=np.zeros(n,np.int32);status='small_or_low_variation_source_retained_as_one_root'
  if n>=params['minimum_fit_group'] and len(variable)>=3:
   candidate=x[:,variable];del x;rank=select_features(candidate,np.zeros(n,np.int8),cfg['n_features']);chosen=variable[rank];emb,fit=pca_fit(candidate[:,rank],cfg['seed'],params['components'],params['pca_iterations']);del candidate;gc.collect();codes=graph_grid(emb,cfg['seed'],params,[cfg['resolution']])[str(cfg['resolution'])];mf=root/'models'/(source+'_TME_qualified_root_model_v1.npz')
   with mf.open('xb') as f:np.savez_compressed(f,mean=fit.mean_,components=fit.components_,explained_variance=fit.explained_variance_,fit_embedding=emb,native_columns=cols[chosen],source_scope_position=pool,source_native_position=q['source_native_position'][pool],native_row=q['native_row'][pool],root_cluster_code=codes)
   model=dict(path=str(mf),SHA256=sha(mf),selected_native_genes=[genes[j] for j in chosen],n_variable_candidates=len(variable),n_selected_features=len(chosen),candidate_feature_axis_SHA256=sha(fp),heldout_physical_aliases_excluded=True);status='all_qualified_source_root_PCA_and_Leiden';del emb,fit
  if 'x' in locals():del x
  ng=int(codes.max())+1;stat=dict(n=np.zeros(ng,np.int64),**{k:np.zeros((ng,len(ecols)),np.float64) for k in ['cp_sum','log_sum','detected']});check(root,'source_all_native_root_profiles_started',source=source,n_fit_pool=n,n_roots=ng)
  for lo in range(0,n,cfg['read_block_size']):
   sel=pool[lo:lo+cfg['read_block_size']];raw=canonical_native(gather_csr_rows(csr,q['native_row'][sel],max_span=8192));total=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(total,q['native_total_RNA'][sel]);update_group_stats(stat,raw[:,ecols].tocsr(),total,codes[lo:lo+len(sel)])
  np.testing.assert_array_equal(stat['n'],np.bincount(codes,minlength=ng));sumcp=stat['cp_sum'].sum(axis=0);sumdet=stat['detected'].sum(axis=0);groups=[]
  for code in range(ng):
   take=np.flatnonzero(codes==code);count=n-len(take);othercp=np.maximum(sumcp-stat['cp_sum'][code],0);otherdet=sumdet-stat['detected'][code];assert otherdet.min()>=0
   groups.append(dict(source_root_cluster=code,n_exact_group=len(take),stable_membership_SHA256=membership_key(source,ids[q['source_native_position'][pool[take]]].astype(str).tolist()),native_rows_SHA256=hashlib.sha256(np.asarray(q['native_row'][pool[take]],'<i8').tobytes()).hexdigest(),own_full_native_gene_marker_candidates_not_DEG=expression_candidates(egenes,stat['cp_sum'][code]/len(take),stat['detected'][code]/len(take),othercp/count if count else None,otherdet/count if count else None,.1,20),n_comparator_records=count,comparator_scope='other_roots_in_same_current_qualified_source',new_identity_assigned=False))
  arrays=root/'outputs'/(source+'_TME_qualified_root_space_v1.npz');stats=root/'outputs'/(source+'_TME_qualified_root_profiles_v1.npz')
  with arrays.open('xb') as f:np.savez_compressed(f,source_scope_position=pool,source_native_position=q['source_native_position'][pool],native_row=q['native_row'][pool],root_cluster_code=codes)
  with stats.open('xb') as f:np.savez_compressed(f,**stat,native_columns=ecols)
  op=root/'outputs'/(source+'_TME_qualified_root_space_v1.json');metadata=dict(source=source,status=status,n_qualified_scope_records=len(q['native_row']),n_fit_pool=n,n_unfitted_retained=len(q['native_row'])-n,n_roots=ng,n_native_nnz_read_in_model_pass=nnz,arrays=dict(path=str(arrays),SHA256=sha(arrays)),statistics=dict(path=str(stats),SHA256=sha(stats)),native_measured_genes=egenes,groups=groups,model=model,working_metadata_SHA256=sha(mp),working_arrays_SHA256=sha(ap),candidate_feature_metadata_SHA256=sha(fp),native_evidence_metadata_SHA256=sha(ep),pool_selection='all_current_source_qualified_fit_records_no_lineage_or_reference_condition',heldout_alias_exclusion_is_physical_not_symbol_only=True,current_labels_QC_states_and_counts_unchanged=True,new_identity_assigned=False,new_CM_started=False);save(op,metadata);results.append(dict(source=source,metadata=op.name,metadata_SHA256=sha(op)));totals.update(n_sources=1,n_qualified_scope_records=len(q['native_row']),n_fit_pool=n,n_unfitted_retained=len(q['native_row'])-n,n_roots=ng,n_native_nnz_read_in_model_pass=nnz);check(root,'source_complete',source=source,totals=dict(totals),metadata_SHA256=sha(op));native._fields.clear();del q,stat;gc.collect()
  assert sum(p.stat().st_size for p in root.rglob('*') if p.is_file())<=cfg['output_cap_bytes']
 ip=root/'outputs/TME_qualified_root_space_index_v1.json';save(ip,results);assert totals['n_fit_pool']==pre['n_fit_pool'] and totals['n_qualified_scope_records']==pre['n_source_qualified_records'];save(root/'provenance/terminal_v1.json',dict(status='ALL_QUALIFIED_TME_SOURCE_ROOTS_AND_NATIVE_MARKERS_COMPLETE_REQUIRES_FULL_AUDIT',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),index_SHA256=sha(ip),all_unfitted_and_archive_records_retained=True,no_old_lineage_or_reference_condition_on_fit_pool=True,current_working_labels_and_RNA_states_unchanged=True,new_CM_started=False,wall_seconds=time.monotonic()-start,peak_RSS_bytes_Linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
