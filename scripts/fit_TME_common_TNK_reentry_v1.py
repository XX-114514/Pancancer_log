"""Source-local unsupervised common-space review of T/NK and unknown parents.

Use every eligible known T/NK and every cell in the audited unknown parents,
including TCR-negative and competing-core controls. No new identity is assigned.
Raw counts remain in native partitions; only model features/embeddings are saved.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
import psutil
import scipy.sparse as sp
from TME_native_measurement_core_v1 import canonical_native

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def checkpoint(root,stage,**values):save(root/'checkpoints'/(str(time.time_ns())+'.json'),dict(observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage=stage,**values))

def main(root):
 begin=time.monotonic();cfg=load(root/'config/spec_v1.json');R=Path(cfg['context_producer']);RA=Path(cfg['context_auditor']);PA=Path(cfg['native_profile_auditor']);PP=Path(cfg['native_profile_producer'])
 for d in [R,RA,PP,PA]:assert load(d/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(RA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(R/'provenance/terminal_v1.json');assert load(PA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(PP/'provenance/terminal_v1.json')
 rc=load(R/'config/spec_v1.json');pc=load(Path(rc['TCR_producer'])/'config/spec_v1.json');J=Path(pc['joined']);accepted=load(rc['working_acceptance']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];w=load(accepted['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi=load(w['index']['path']);ri={r['source']:r for r in load(R/'outputs/unknown_parent_TCR_context_index_v1.json')}
 lm=Path(load(J/'outputs/joined_native_atlas_delivery_v1.json')['local_review']);G=Path(load(lm/'config/local_identity_review_spec_v1.json')['graphs']);graph_cfg=load(G/'config/full_local_graph_spec_v1.json');I=Path(graph_cfg['local_inputs']);ii={r['source']:r for r in load(I/'outputs/full_native_local_discovery_input_index_v1.json')};params=dict(load(graph_cfg['original_model_spec']));params['threads']=cfg['threads'];seed=graph_cfg['seed'];resolution=graph_cfg['resolution'];preflight=cfg['frozen_preflight']
 assert preflight['working_contract_SHA256']==sha(rc['working_acceptance']);expected={r['source']:r for r in preflight['sources']}
 for binding in graph_cfg['bindings'].values():assert sha(binding['path'])==binding['sha256']
 sys.path.insert(0,str(G/'scripts'));sys.path.insert(0,str(J/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256,row_axis_sha256
 from native_panel_aggregate_v1 import gather_csr_rows
 from annotation_representation_pilot_v1 import normalized_log_counts,select_features,pca_fit
 from compare_annotation_granularity_v1 import graph_grid
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);(root/'models').mkdir();(root/'checkpoints').mkdir();totals=collections.Counter();results=[];resource_cfg=dict(pc,output_cap_bytes=cfg['output_cap_bytes'])
 # Execution order is source size, never endpoint, marker effect or CM ranking.
 for wr in sorted(wi,key=lambda r:(expected[r['source']]['n_new_fit_pool'],r['source'])):
  source=wr['source'];storage_gate(root,resource_cfg);mp=base/wr['metadata'];wa=base/wr['arrays'];assert sha(mp)==wr['metadata_SHA256'] and sha(wa)==wr['arrays_SHA256'];wm=load(mp)
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','consensus_major_code','annotation_group_code','fit_pool_candidate_not_final_QC','native_total_RNA']}
  rp=R/'outputs'/ri[source]['metadata'];assert sha(rp)==ri[source]['metadata_SHA256'];review=load(rp);parent_codes=[g['previous_group_code'] for g in review['groups']];major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];known=np.isin(major,['T','NK','T_NK']);control=np.isin(q['annotation_group_code'],parent_codes);pool=np.flatnonzero((known|control)&q['fit_pool_candidate_not_final_QC']);n=len(pool);assert n==expected[source]['n_new_fit_pool']
  im=load(I/'outputs'/ii[source]['metadata']);fp=I/'outputs'/im['feature_candidates'];assert sha(fp)==im['feature_candidates_sha256']==expected[source]['candidate_feature_axis_SHA256'];features=load(fp);genes=features['genes'];cols=np.asarray(features['native_columns'],np.int64);assert len(cols)==len(set(cols))==len(set(genes));assert not set(cols)&set(features['heldout_physical_alias_columns_excluded'])
  ids=native.native_ids(source);assert text_axis_sha256(ids)==wm['native_ids_source_SHA256'];np.testing.assert_array_equal(native.rows(source)[q['source_native_position'][pool]],q['native_row'][pool]);codes=np.full(n,-1,np.int32);status='empty_qualified_pool';model_record=None;actual_nnz=0;started=time.monotonic()
  if n and n<params['minimum_fit_group']:codes[:]=0;status='small_source_pool_kept_coarse_no_embedding'
  elif n:
   upper=n*len(genes)*16+n*min(graph_cfg['n_features'],len(genes))*12+cfg['RAM_reserve_bytes'];available=psutil.virtual_memory().available
   if upper>available:raise MemoryError('Fresh conservative source feature/PCA RAM bound failed')
   checkpoint(root,'native_pool_read_started',source=source,n_pool=n,estimated_peak_bound_bytes=upper,available_RAM_bytes=available)
   csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']];chunks=[]
   for lo in range(0,n,cfg['read_block_size']):
    take=pool[lo:lo+cfg['read_block_size']];raw=canonical_native(gather_csr_rows(csr,q['native_row'][take],max_span=8192));total=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(total,q['native_total_RNA'][take]);x,positive=normalized_log_counts(raw[:,cols],total);assert positive.all();chunks.append(x);actual_nnz+=raw.nnz
   x=sp.vstack(chunks,format='csr');del chunks,raw;gc.collect();mean=np.asarray(x.mean(axis=0)).ravel();variance=np.maximum(np.asarray(x.multiply(x).mean(axis=0)).ravel()-mean*mean,0);variable=np.flatnonzero((variance>0)&(np.asarray((x>0).sum(axis=0)).ravel()>=5))
   if len(variable)<3:codes[:]=0;status='insufficient_variable_features_kept_coarse_no_embedding'
   else:
    candidate=x[:,variable];del x;rank=select_features(candidate,np.zeros(n,np.int8),graph_cfg['n_features']);chosen=variable[rank];fitx=candidate[:,rank];del candidate;gc.collect();embedding,model=pca_fit(fitx,seed,params['components'],params['pca_iterations']);del fitx;gc.collect();codes=graph_grid(embedding,seed,params,[resolution])[str(resolution)];status='qualified_source_common_TNK_and_unknown_control_graph'
    mf=root/'models'/(source+'_common_TNK_reentry_model_v1.npz')
    with mf.open('xb') as f:np.savez_compressed(f,mean=model.mean_,components=model.components_,explained_variance=model.explained_variance_,fit_embedding=embedding,native_columns=cols[chosen],source_scope_position=pool,source_native_position=q['source_native_position'][pool],native_row=q['native_row'][pool],reentry_group_code=codes)
    model_record=dict(path=str(mf),SHA256=sha(mf),selected_native_genes=[genes[j] for j in chosen],native_candidate_feature_metadata_SHA256=sha(fp),n_selected_features=len(chosen),n_variable_candidates=len(variable),fit_native_rows_SHA256=row_axis_sha256(q['native_row'][pool]),feature_selection='same_source_log_dispersion_on_complete_qualified_common_pool',original_feature_axis_and_heldout_alias_exclusion_preserved=True,counts_canonicalized_before_log=True,fit_is_current_scope_transductive_not_independent_validation=True)
    del embedding,model
   if 'x' in locals():del x
  output=root/'outputs'/(source+'_common_TNK_reentry_groups_v1.npz')
  with output.open('xb') as f:np.savez_compressed(f,source_scope_position=pool,source_native_position=q['source_native_position'][pool],native_row=q['native_row'][pool],reentry_group_code=codes,baseline_known_T_NK=known[pool],unresolved_parent_control=control[pool],previous_annotation_group_code=q['annotation_group_code'][pool])
  with np.load(output,allow_pickle=False) as f:np.testing.assert_array_equal(f['native_row'],q['native_row'][pool]);np.testing.assert_array_equal(f['reentry_group_code'],codes)
  record=dict(source=source,status=status,n_qualified_scope=len(major),n_fit_pool=n,n_known_T_NK_fit=int(known[pool].sum()),n_unknown_parent_control_fit=int(control[pool].sum()),n_groups=len(set(codes.tolist())) if n else 0,n_native_nnz_read=actual_nnz,arrays=str(output),arrays_SHA256=sha(output),model=model_record,current_working_metadata_SHA256=sha(mp),current_working_arrays_SHA256=sha(wa),unknown_parent_context_metadata_SHA256=sha(rp),all_training_records_current_source_qualified=True,clinical_or_CM_or_new_identity_labels_not_used_for_fit=True,new_identity_assigned=False,wall_seconds=time.monotonic()-started)
  op=root/'outputs'/(source+'_common_TNK_reentry_groups_v1.json');save(op,record);results.append(dict(source=source,metadata=op.name,metadata_SHA256=sha(op)));totals.update(n_sources=1,n_qualified_scope=len(major),n_fit_pool=n,n_known_T_NK_fit=record['n_known_T_NK_fit'],n_unknown_parent_control_fit=record['n_unknown_parent_control_fit'],n_groups=record['n_groups'],n_native_nnz_read=actual_nnz);checkpoint(root,'source_complete',source=source,totals=dict(totals),source_metadata_SHA256=sha(op));native._fields.clear();del q,codes;gc.collect()
 assert totals['n_fit_pool']==preflight['n_new_fit_pool'] and totals['n_unknown_parent_control_fit']==preflight['n_unknown_parent_control_fit'];ip=root/'outputs/common_TNK_reentry_index_v1.json';save(ip,results)
 save(root/'provenance/terminal_v1.json',dict(status='COMPLETE_CURRENT_QUALIFIED_COMMON_TNK_AND_UNKNOWN_CONTROL_DISCOVERY_REQUIRES_AUDIT',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),index_SHA256=sha(ip),current_working_contract_SHA256=sha(rc['working_acceptance']),new_identity_assigned=False,new_CM_started=False,not_controlled_ABC_winner_comparison=True,wall_seconds=time.monotonic()-begin,peak_RSS_bytes_Linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
