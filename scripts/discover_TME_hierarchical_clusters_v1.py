"""Cluster-led annotation and one nested subclustering round on accepted roots.

Existing accepted source-local roots define the first tree level, including all
unknown/competing controls. Fit child features anew inside each eligible root.
All labels are a separate provisional cluster view; working cell labels stay.
"""
import argparse,collections,datetime,gc,hashlib,json,re,resource,subprocess,sys,time
from pathlib import Path
import numpy as np
import psutil
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score
from TME_native_measurement_core_v1 import canonical_native,update_group_stats,expression_candidates
from TME_cluster_annotation_policy_v1 import GATE_ORDER,coherent_flags,annotate_cluster,state_context

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def checkpoint(root,stage,**v):save(root/'checkpoints'/(str(time.time_ns())+'.json'),dict(observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage=stage,**v))

def child_graphs(embedding,seed,parameters,resolutions):
 import anndata as ad
 import scanpy as sc
 import numba
 sc.settings.n_jobs=parameters['threads'];numba.set_num_threads(parameters['threads'])
 a=ad.AnnData(sp.csr_matrix((len(embedding),0)));a.obsm['X_child']=embedding
 sc.pp.neighbors(a,use_rep='X_child',n_neighbors=min(parameters['neighbors'],len(embedding)-1),random_state=seed,transformer='pynndescent',metric='euclidean')
 out={};curves=[]
 for resolution in resolutions:
  pair=[]
  for offset in [0,1]:
   sc.tl.leiden(a,resolution=resolution,random_state=seed+offset,flavor='igraph',n_iterations=2,directed=False)
   _,c=np.unique(a.obs['leiden'].to_numpy().astype(str),return_inverse=True);pair.append(c.astype(np.int32))
  sizes=np.bincount(pair[0]);ari=float(adjusted_rand_score(pair[0],pair[1]));out[str(resolution)]=pair
  curves.append(dict(resolution=resolution,n_children=len(sizes),minimum_child_cells=int(sizes.min()),n_children_at_least_30=int((sizes>=30).sum()),same_graph_Leiden_seed_ARI=ari))
 return out,curves

def main(root):
 global GATE_ORDER,coherent_flags,annotate_cluster,state_context
 begin=time.monotonic();cfg=load(root/'config/spec_v1.json');M=Path(cfg['root_model']);S=Path(cfg['root_profiles']);A=Path(cfg['root_profile_auditor']);NA=Path(cfg['root_model_native_auditor']);D=Path(cfg['working_delivery']);frozen=cfg['frozen_preflight']
 full=cfg.get('root_protocol')=='qualified_full_TME_source_roots_v1'
 if full:
  from TME_all_lineage_cluster_policy_v1 import GATE_ORDER,coherent_flags,annotate_cluster,state_context
 for d in [M,S,A,NA,D]:assert load(d/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert sha(M/'provenance/terminal_v1.json')==frozen['root_model_terminal_SHA256'] and sha(A/'provenance/terminal_v1.json')==frozen['root_profile_audit_terminal_SHA256']
 ac=D/'outputs/accepted_TME_scoped_working_annotation_v1.json';assert sha(ac)==frozen['working_contract_SHA256'];accepted=load(ac);w=load(accepted['producer_manifest']['path']);assert sha(accepted['producer_manifest']['path'])==accepted['producer_manifest']['SHA256'];base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};history=load(w['display_number_registry']['path']);assert sha(w['display_number_registry']['path'])==w['display_number_registry']['SHA256']
 if full:
  assert sha(cfg['display_registry_path'])==frozen['display_registry_SHA256'];history=load(cfg['display_registry_path']);assert sha(NA/'provenance/terminal_v1.json')==frozen['root_native_audit_terminal_SHA256']
 mc=load(M/'config/spec_v1.json')
 if full:kc=load(mc['resource_policy_source']);J=Path(mc['joined'])
 else:rc=load(Path(mc['context_producer'])/'config/spec_v1.json');kc=load(Path(rc['TCR_producer'])/'config/spec_v1.json');J=Path(kc['joined'])
 jm=load(J/'outputs/joined_native_atlas_delivery_v1.json');L=Path(jm['local_review']);G=Path(load(L/'config/local_identity_review_spec_v1.json')['graphs']);gcfig=load(G/'config/full_local_graph_spec_v1.json');parameters=dict(load(gcfig['original_model_spec']));parameters['threads']=cfg['threads'];seed=cfg['seed']
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(G/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
 from native_panel_aggregate_v1 import gather_csr_rows
 from annotation_representation_pilot_v1 import select_features,pca_fit
 from atlas_identity_release_keys_v1 import membership_key,assign_display_numbers,display_label
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 profile_index='TME_qualified_root_space_index_v1.json' if full else 'common_TNK_group_native_profile_index_v1.json';model_index='TME_qualified_root_space_index_v1.json' if full else 'common_TNK_reentry_index_v1.json'
 native=NativeArchiveCollection(jm['archive_manifest']);(root/'models').mkdir();(root/'checkpoints').mkdir();spidx={r['source']:r for r in load(S/'outputs'/profile_index)};expected={r['source']:r for r in frozen['sources']};resource_cfg=dict(kc,output_cap_bytes=cfg['output_cap_bytes']);results=[];allnodes=[];totals=collections.Counter()
 test_name='test_TME_all_lineage_cluster_policy_v1.py' if full else 'test_TME_cluster_annotation_policy_v1.py'
 test=subprocess.run([sys.executable,'-B',str(root/'tests'/test_name)],capture_output=True,text=True,timeout=60);(root/'provenance/policy_test_stdout.txt').write_text(test.stdout);(root/'provenance/policy_test_stderr.txt').write_text(test.stderr);assert test.returncode==0 and 'Ran 8 tests' in test.stderr and 'skipped' not in test.stderr
 for item in sorted(load(M/'outputs'/model_index),key=lambda r:(expected[r['source']]['n_pool'],r['source'])):
  source=item['source'];storage_gate(root,resource_cfg);mp=M/'outputs'/item['metadata'];assert sha(mp)==item['metadata_SHA256'];model=load(mp);model_arrays=model['arrays']['path'] if full else model['arrays'];model_arrays_sha=model['arrays']['SHA256'] if full else model['arrays_SHA256'];assert sha(model_arrays)==model_arrays_sha
  with np.load(model_arrays,allow_pickle=False) as f:pool=f['source_scope_position'];native_rows=f['native_row'];roots=f['root_cluster_code' if full else 'reentry_group_code']
  wr=wi[source];wp=base/wr['metadata'];wa=base/wr['arrays'];assert sha(wp)==wr['metadata_SHA256'] and sha(wa)==wr['arrays_SHA256'];wm=load(wp)
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','native_total_RNA','fit_pool_candidate_not_final_QC','consensus_major_code','consensus_celltype_code']}
  np.testing.assert_array_equal(native_rows,q['native_row'][pool]);assert q['fit_pool_candidate_not_final_QC'][pool].all();assert len(pool)==expected[source]['n_pool']
  if full:np.testing.assert_array_equal(pool,np.flatnonzero(q['fit_pool_candidate_not_final_QC']))
  pr=spidx[source];pp=S/'outputs'/pr['metadata'];assert sha(pp)==pr['metadata_SHA256'];profile=load(pp);assert sha(profile['statistics']['path'])==profile['statistics']['SHA256'];genes=profile['native_measured_genes']
  with np.load(profile['statistics']['path'],allow_pickle=False) as f:rootstats={k:f[k] for k in ['n','cp_sum','log_sum','detected']};cols=f['native_columns']
  anchor_mode='T_chain_unmeasured_alternate_CD3_accessory_two_signal_RNA_evidence' if full and not any(g in genes for g in ['TRAC','TRBC1','TRBC2','TRDC']) and all(g in genes for g in ['CD3D','CD3E']) and any(g in genes for g in ['CD2','CD7','TRAT1']) and sum(g in genes for g in ['CD247','LCK','LAT','ZAP70','BCL11B'])>=2 else 'original_native_module_assay_coverage_in_proposal'
  ng=len(rootstats['n']);np.testing.assert_array_equal(rootstats['n'],np.bincount(roots,minlength=ng));assert ng==expected[source]['n_root_clusters']
  feature_indices=np.asarray([j for j,g in enumerate(genes) if not re.match(r'^(MT-|RPL|RPS|ENSG|IG[HKL][VDJC]|TR[ABGD][VJ])',g)],np.int64)
  ids=native.native_ids(source);assert text_axis_sha256(ids)==wm['native_ids_source_SHA256'];units=native.frozen_expression_units(source)[q['source_native_position'][pool]].astype(str);csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']]
  coherent=np.zeros((len(pool),len(GATE_ORDER)),bool);assessed=None;candidate_nodes=np.full(len(pool),-1,np.int32);parents=[];node_parent=list(range(ng));child_number=[None]*ng;node_cursor=ng;raw_nnz=0
  for parent in range(ng):
   take=np.flatnonzero(roots==parent);n=len(take);chunks=[];fit=n>=cfg['minimum_parent_cells'];upper=n*len(genes)*24+cfg['RAM_reserve_bytes'];assert upper<psutil.virtual_memory().available
   checkpoint(root,'parent_native_feature_read',source=source,parent=parent,n_records=n,eligible_child_fit=fit)
   for lo in range(0,n,cfg['read_block_size']):
    pos=take[lo:lo+cfg['read_block_size']];scope=pool[pos];raw=canonical_native(gather_csr_rows(csr,q['native_row'][scope],max_span=8192));total=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(total,q['native_total_RNA'][scope]);x=raw[:,cols].tocsr();flags,measured=coherent_flags(x,genes);coherent[pos]=flags
    if assessed is None:assessed=measured
    else:assert assessed==measured
    if fit:
     log=x[:,feature_indices].multiply((10000/total)[:,None]).tocsr();log.data=np.log1p(log.data);chunks.append(log.astype(np.float32))
    raw_nnz+=raw.nnz
   call=annotate_cluster(n,dict(zip(GATE_ORDER,coherent[take].sum(axis=0).tolist())),assessed,cfg['minimum_named_group_cells'],cfg['family_fraction'])
   p=dict(source_root_cluster=parent,n_parent_records=n,cluster_RNA_proposal=call,child_fit_status='below_frozen_parent_size',candidate_children=[],split_retained=False,primary_resolution=cfg['primary_resolution'],parameter_sensitivity_not_optimized_by_marker_or_clinical_result=True)
   if fit:
    x=sp.vstack(chunks,format='csr');chunks=[];mean=np.asarray(x.mean(axis=0)).ravel();var=np.maximum(np.asarray(x.multiply(x).mean(axis=0)).ravel()-mean*mean,0);valid=np.flatnonzero((var>0)&(np.asarray((x>0).sum(axis=0)).ravel()>=5))
    if len(valid)<3:p['child_fit_status']='insufficient_variable_features_keep_parent'
    else:
     rank=select_features(x[:,valid],np.zeros(n,np.uint8),cfg['n_features']);chosen=feature_indices[valid[rank]];embedding,fitmodel=pca_fit(x[:,valid[rank]],seed,parameters['components'],parameters['pca_iterations']);partitions,curves=child_graphs(embedding,seed,parameters,cfg['resolutions']);p['resolution_sensitivity']=curves;pair=partitions[str(cfg['primary_resolution'])];primary=pair[0];sizes=np.bincount(primary);primarycurve=next(c for c in curves if c['resolution']==cfg['primary_resolution']);candidate=primarycurve['n_children_at_least_30']>=2 and primarycurve['same_graph_Leiden_seed_ARI']>=cfg['minimum_seed_ARI']
     modelpath=root/'models'/(source+'_root_'+str(parent)+'_child_model_v1.npz')
     with modelpath.open('xb') as f:np.savez_compressed(f,mean=fitmodel.mean_,components=fitmodel.components_,explained_variance=fitmodel.explained_variance_,fit_embedding=embedding,native_columns=cols[chosen],pool_positions=take,native_row=native_rows[take],**{'resolution_'+str(r)+'_seed_'+str(j):v for r,pair2 in partitions.items() for j,v in enumerate(pair2)})
     p.update(child_fit_status='stable_candidate_requires_sibling_marker_contrast' if candidate else 'seed_or_size_diagnostic_keep_parent',model=dict(path=str(modelpath),SHA256=sha(modelpath)),n_variable_features=len(valid),n_selected_features=len(chosen),model_features_may_include_identity_and_RNA_state_genes=True,heldout_state_independence_not_claimed=True)
     if candidate:
      for child in range(len(sizes)):
       node=node_cursor;node_cursor+=1;candidate_nodes[take[primary==child]]=node;node_parent.append(parent);child_number.append(child);p['candidate_children'].append(node)
     del embedding,fitmodel,partitions
    del x
   parents.append(p);gc.collect()
  ngenes=len(cols);bound=node_cursor*ngenes*3*8+cfg['RAM_reserve_bytes'];assert bound<psutil.virtual_memory().available
  stat=dict(n=np.zeros(node_cursor,np.int64),**{k:np.zeros((node_cursor,ngenes),np.float64) for k in ['cp_sum','log_sum','detected']})
  for k in stat:stat[k][:ng]=rootstats[k]
  for lo in range(0,len(pool),cfg['read_block_size']):
   pos=np.arange(lo,min(len(pool),lo+cfg['read_block_size']));sel=candidate_nodes[pos]>=0
   if not sel.any():continue
   pos=pos[sel];scope=pool[pos];raw=canonical_native(gather_csr_rows(csr,q['native_row'][scope],max_span=8192));rna=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(rna,q['native_total_RNA'][scope]);update_group_stats(stat,raw[:,cols].tocsr(),rna,candidate_nodes[pos])
  node_rows=[];sourceN=int(rootstats['n'].sum());sourceCP=rootstats['cp_sum'].sum(axis=0);sourceLog=rootstats['log_sum'].sum(axis=0);sourceDet=rootstats['detected'].sum(axis=0)
  oldmajor=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code'][pool]];oldfine=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code'][pool]]
  for node in range(node_cursor):
   parent=node_parent[node];sel=np.flatnonzero((roots==node) if node<ng else (candidate_nodes==node));n=len(sel);assert n==stat['n'][node] and n>0;comparatorN=sourceN-n if node<ng else int(stat['n'][parent])-n
   cp=(sourceCP if node<ng else stat['cp_sum'][parent])-stat['cp_sum'][node];lg=(sourceLog if node<ng else stat['log_sum'][parent])-stat['log_sum'][node];de=(sourceDet if node<ng else stat['detected'][parent])-stat['detected'][node];assert cp.min()>=-1e-7 and de.min()>=0;cp=np.maximum(cp,0)
   markers=expression_candidates(genes,stat['cp_sum'][node]/n,stat['detected'][node]/n,cp/comparatorN if comparatorN else None,de/comparatorN if comparatorN else None,.1,20)
   # Clusters are named from their own statistics; never borrow a parent marker.
   best=next((m for m in markers if re.fullmatch('[A-Za-z0-9][A-Za-z0-9.-]*',m['gene'])),None)
   call=annotate_cluster(n,dict(zip(GATE_ORDER,coherent[sel].sum(axis=0).tolist())),assessed,cfg['minimum_named_group_cells'],cfg['family_fraction']);missing=np.isin(np.char.lower(np.char.strip(units[sel])),['','nan','none','na','n/a','unknown','__missing__','__source_value_missing__'])
   record=dict(source=source,source_node_code=node,source_root_cluster=parent,level=1 if node<ng else 2,child_cluster_number=child_number[node],n_records=n,stable_membership_SHA256=membership_key(source,ids[q['source_native_position'][pool[sel]]].astype(str).tolist()),native_rows_SHA256=hashlib.sha256(np.asarray(native_rows[sel],'<i8').tobytes()).hexdigest(),cluster_RNA_proposal=call,marker=best['gene'] if best else None,marker_candidates=markers,marker_evidence_type=best['marker_evidence_type'] if best else 'no_eligible_own_marker',parent_marker_borrowed=False,comparator_scope=('other_roots_in_same_qualified_source' if full else 'other_roots_in_same_source_frozen_common_pool') if node<ng else 'other_candidate_children_within_exact_same_parent',n_comparator_records=comparatorN,descriptive_state_programs=state_context(genes,stat['log_sum'][node]/n,lg/comparatorN if comparatorN else None,stat['detected'][node]/n),n_declared_expression_units_not_verified_patients_or_libraries=len(set(units[sel][~missing])),n_records_missing_declared_expression_unit=int(missing.sum()),old_v5_major_counts=dict(collections.Counter(map(str,oldmajor[sel]))),old_v5_fine_counts=dict(collections.Counter(map(str,oldfine[sel]))),old_labels_not_used_as_classifier_or_fit_input=True,final_leaf=False)
   if full:record['assay_aware_T_evidence_mode']=anchor_mode
   record['sibling_marker_separation_supported']=bool(node>=ng and n>=cfg['minimum_named_group_cells'] and any(m['descriptive_log2_mean_ratio'] is not None and m['descriptive_log2_mean_ratio']>=cfg['minimum_sibling_log2_ratio'] and m['native_detection_fraction']-m['comparator_detection_fraction']>=cfg['minimum_sibling_detection_difference'] for m in markers));node_rows.append(record)
  final=roots.copy()
  for p in parents:
   children=p['candidate_children'];supported=[c for c in children if node_rows[c]['sibling_marker_separation_supported']]
   p['n_children_with_native_sibling_marker_contrast']=len(supported);p['split_retained']=len(supported)>=2
   if p['split_retained']:
    take=np.flatnonzero(roots==p['source_root_cluster']);assert (candidate_nodes[take]>=ng).all();final[take]=candidate_nodes[take]
   elif children:p['child_fit_status']='sibling_marker_separation_insufficient_keep_parent'
  for node in set(map(int,final)):node_rows[node]['final_leaf']=True
  for parent in parents:
   child=parent['candidate_children']
   if child:
    for k in ['n','cp_sum','log_sum','detected']:
     a=stat[k][child].sum(axis=0);b=stat[k][parent['source_root_cluster']]
     if k in ['n','detected']:np.testing.assert_array_equal(a,b)
     else:np.testing.assert_allclose(a,b,rtol=2e-10,atol=1e-7)
  assert sum(r['n_records'] for r in node_rows if r['final_leaf'])==len(pool)
  arrays=root/'outputs'/(source+'_cluster_hierarchy_v1.npz')
  with arrays.open('xb') as f:np.savez_compressed(f,source_scope_position=pool,source_native_position=q['source_native_position'][pool],native_row=native_rows,root_cluster_code=roots,candidate_child_node=candidate_nodes,final_leaf_node=final)
  stats=root/'outputs'/(source+'_cluster_hierarchy_native_profiles_v1.npz')
  with stats.open('xb') as f:np.savez_compressed(f,**stat,native_columns=cols)
  meta=dict(source=source,n_pool_records=len(pool),n_roots=ng,n_candidate_nodes=node_cursor,n_final_leaves=len(set(final.tolist())),native_measured_genes=genes,gate_assessed=assessed,parents=parents,nodes=node_rows,arrays=dict(path=str(arrays),SHA256=sha(arrays)),statistics=dict(path=str(stats),SHA256=sha(stats)),root_model_metadata_SHA256=sha(mp),root_profile_metadata_SHA256=sha(pp),current_working_metadata_SHA256=sha(wp),current_working_arrays_SHA256=sha(wa),existing_current_cell_annotations_unchanged=True,new_CM_started=False)
  op=root/'outputs'/(source+'_cluster_hierarchy_v1.json');save(op,meta);results.append(dict(source=source,metadata=op.name,metadata_SHA256=sha(op)));allnodes+=node_rows;totals.update(n_sources=1,n_pool_records=len(pool),n_root_clusters=ng,n_candidate_nodes=node_cursor,n_final_leaf_clusters=meta['n_final_leaves'],n_roots_with_retained_child_split=sum(p['split_retained'] for p in parents),n_raw_native_nnz_in_feature_pass=raw_nnz);checkpoint(root,'source_complete',source=source,totals=dict(totals),metadata_SHA256=sha(op));native._fields.clear();del q,stat,rootstats,coherent,node_rows;gc.collect()
 assigned=assign_display_numbers([r['stable_membership_SHA256'] for r in allnodes],history)
 for r in allnodes:r.update(global_display_number=assigned['current'][r['stable_membership_SHA256']],display_label=display_label(r['cluster_RNA_proposal']['celltype'],r['marker'],assigned['current'][r['stable_membership_SHA256']]))
 catalog=root/'outputs/cluster_hierarchy_named_catalog_v1.json';save(catalog,dict(nodes=allnodes,previous_historical_numbers_preserved=True,numbers_do_not_define_biological_types=True));ip=root/'outputs/cluster_hierarchy_index_v1.json';save(ip,results);save(root/'outputs/cluster_hierarchy_display_registry_v1.json',assigned['registry'])
 summary=dict(totals=dict(totals),root_cluster_proposal_counts=dict(collections.Counter(r['cluster_RNA_proposal']['celltype'] for r in allnodes if r['level']==1)),final_leaf_cluster_proposal_counts=dict(collections.Counter(r['cluster_RNA_proposal']['celltype'] for r in allnodes if r['final_leaf'])),final_leaf_record_proposal_counts=dict(sum((collections.Counter({r['cluster_RNA_proposal']['celltype']:r['n_records']}) for r in allnodes if r['final_leaf']),collections.Counter())),root_protocol=cfg.get('root_protocol','frozen_TNK_and_unknown_pool_v1'),no_existing_working_labels_or_RNA_states_or_counts_changed=True,independent_marker_and_membership_audit_required=True,new_CM_started=False)
 save(root/'outputs/cluster_hierarchy_summary_v1.json',summary);assert totals['n_pool_records']==frozen['n_pool'] and totals['n_root_clusters']==frozen['n_root_clusters']
 save(root/'provenance/terminal_v1.json',dict(status='CLUSTER_HIERARCHY_AND_NATIVE_MARKER_PROPOSALS_COMPLETE_REQUIRES_AUDIT',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),index_SHA256=sha(ip),catalog_SHA256=sha(catalog),all_candidate_children_partition_original_parent_and_native_statistics_recompose=True,current_working_annotation_unchanged=True,new_CM_started=False,wall_seconds=time.monotonic()-begin,peak_RSS_bytes_Linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
