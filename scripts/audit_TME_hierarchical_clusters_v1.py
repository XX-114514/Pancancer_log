"""Full independent CSR statistics and tree/member audit of cluster proposals.

Replays the frozen proposal policy, not biological ground truth. Checks all
native moments and coherent modules, marker contrasts, final leaf coverage,
saved model axes and seed diagnostics; does not claim a second PCA solver.
"""
import argparse,collections,datetime,gc,hashlib,json,re,sys,time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from sklearn.metrics import adjusted_rand_score
from TME_cluster_annotation_policy_v1 import GATES,GATE_ORDER,annotate_cluster,state_context

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):
 with Path(p).open('x') as f:json.dump(v,f,ensure_ascii=False,indent=2,allow_nan=False);f.write('\n')
def equal(a,b):
 if isinstance(a,dict):
  assert set(a)==set(b)
  for k in a:equal(a[k],b[k])
 elif isinstance(a,list):
  assert len(a)==len(b)
  for x,y in zip(a,b):equal(x,y)
 elif isinstance(a,float):np.testing.assert_allclose(a,b,rtol=2e-10,atol=1e-7)
 else:assert a==b,(a,b)

def independent_markers(genes,cp,det,otherCP,otherDet):
 ok=np.asarray([not re.match(r'^(MT-|RPL|RPS|ENSG|IG[HKL][VDJC])',g) for g in genes])&(det>=.1)&(cp>0)
 effect=np.log2((cp+.1)/(otherCP+.1)) if otherCP is not None else None;order=[];kind='own_subset_expression_without_comparator_not_inferential_DEG'
 if effect is not None:
  valid=np.flatnonzero(ok&(det>otherDet)&(effect>0));score=effect*np.maximum(det-otherDet,0);order=sorted(valid,key=lambda j:(-score[j],genes[j]))[:20];kind='own_subset_expression_specificity_not_inferential_DEG'
 if not order:
  order=sorted(np.flatnonzero(ok),key=lambda j:(-cp[j]*det[j],genes[j]))[:20]
  if effect is not None:kind='own_subset_expression_without_positive_specificity_not_inferential_DEG'
 return [dict(gene=genes[j],native_mean_CP10K=float(cp[j]),native_detection_fraction=float(det[j]),comparator_detection_fraction=float(otherDet[j]) if otherDet is not None else None,descriptive_log2_mean_ratio=float(effect[j]) if effect is not None else None,marker_evidence_type=kind) for j in order]

def main(root):
 global GATES,GATE_ORDER,annotate_cluster,state_context
 begin=time.monotonic();cfg=load(root/'config/spec_v1.json');P=Path(cfg['producer']);pc=load(P/'config/spec_v1.json');S=Path(pc['root_profiles']);mc=load(Path(pc['root_model'])/'config/spec_v1.json');full=pc.get('root_protocol')=='qualified_full_TME_source_roots_v1'
 if full:
  from TME_all_lineage_cluster_policy_v1 import GATES,GATE_ORDER,annotate_cluster,state_context
  kc=load(mc['resource_policy_source']);J=Path(mc['joined'])
 else:rc=load(Path(mc['context_producer'])/'config/spec_v1.json');kc=load(Path(rc['TCR_producer'])/'config/spec_v1.json');J=Path(kc['joined'])
 jm=load(J/'outputs/joined_native_atlas_delivery_v1.json');D=Path(pc['working_delivery']);accepted=load(D/'outputs/accepted_TME_scoped_working_annotation_v1.json');w=load(accepted['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};EA=Path(cfg['accepted_native_audit']);assert load(EA/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 sys.path.insert(0,str(EA/'scripts'));sys.path.insert(0,str(J/'scripts'))
 from audit_full_TME_native_evidence_v1 import gather
 from native_archive_collection_v1 import NativeArchiveCollection,text_axis_sha256
 from atlas_identity_release_keys_v1 import membership_key,assign_display_numbers,display_label
 native=NativeArchiveCollection(jm['archive_manifest']);expected={r['source']:r for r in pc['frozen_preflight']['sources']};todo=set(expected);checked={};totals=collections.Counter();allnodes=[]
 while todo:
  progress=False
  for source in sorted(todo.copy()):
   mp=P/'outputs'/(source+'_cluster_hierarchy_v1.json')
   if not mp.exists():continue
   try:m=load(mp)
   except json.JSONDecodeError:continue
   progress=True;assert m['current_working_metadata_SHA256']==sha(base/wi[source]['metadata']) and m['current_working_arrays_SHA256']==sha(base/wi[source]['arrays']);wm=load(base/wi[source]['metadata']);assert sha(m['arrays']['path'])==m['arrays']['SHA256'] and sha(m['statistics']['path'])==m['statistics']['SHA256']
   with np.load(m['arrays']['path'],allow_pickle=False) as f:pool=f['source_scope_position'];positions=f['source_native_position'];rows=f['native_row'];roots=f['root_cluster_code'];children=f['candidate_child_node'];final=f['final_leaf_node']
   with np.load(base/wi[source]['arrays'],allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','native_total_RNA','fit_pool_candidate_not_final_QC','consensus_major_code','consensus_celltype_code']}
   np.testing.assert_array_equal(rows,q['native_row'][pool]);np.testing.assert_array_equal(positions,q['source_native_position'][pool]);assert q['fit_pool_candidate_not_final_QC'][pool].all() and len(set(rows.tolist()))==len(rows)==expected[source]['n_pool']
   if full:np.testing.assert_array_equal(pool,np.flatnonzero(q['fit_pool_candidate_not_final_QC']))
   with np.load(m['statistics']['path'],allow_pickle=False) as f:produced={k:f[k] for k in ['n','cp_sum','log_sum','detected']};cols=f['native_columns']
   if full:
    root_entry=next(r for r in load(Path(pc['root_model'])/'outputs/TME_qualified_root_space_index_v1.json') if r['source']==source);root_meta_path=Path(pc['root_model'])/'outputs'/root_entry['metadata'];assert sha(root_meta_path)==root_entry['metadata_SHA256']==m['root_model_metadata_SHA256']==m['root_profile_metadata_SHA256']==expected[source]['metadata_SHA256'];up=load(root_meta_path);assert sha(up['arrays']['path'])==up['arrays']['SHA256'] and sha(up['statistics']['path'])==up['statistics']['SHA256']
    with np.load(up['arrays']['path'],allow_pickle=False) as f:np.testing.assert_array_equal(f['source_scope_position'],pool);np.testing.assert_array_equal(f['root_cluster_code'],roots);np.testing.assert_array_equal(f['native_row'],rows)
    with np.load(up['statistics']['path'],allow_pickle=False) as f:np.testing.assert_array_equal(f['native_columns'],cols)
    assert up['native_measured_genes']==m['native_measured_genes']
   nnode=m['n_candidate_nodes'];genes=m['native_measured_genes'];by={g:j for j,g in enumerate(genes)};assert len(cols)==len(set(cols.tolist()))==len(set(genes));stat=dict(n=np.zeros(nnode,np.int64),**{k:np.zeros((nnode,len(cols)),np.float64) for k in ['cp_sum','log_sum','detected']});flags=np.zeros((len(pool),len(GATE_ORDER)),bool);assessed={name:all(any(g in by for g in module) for module in modules) for name,modules in GATES.items()}
   alternate=full and not any(g in by for g in ['TRAC','TRBC1','TRBC2','TRDC']) and all(g in by for g in ['CD3D','CD3E']) and any(g in by for g in ['CD2','CD7','TRAT1']) and sum(g in by for g in ['CD247','LCK','LAT','ZAP70','BCL11B'])>=2
   if alternate:assessed['T']=True
   for name in ['Treg','CD4_T','CD8_T','gamma_delta_T_RNA_candidate']:assessed[name]&=assessed['T']
   if full:
    from TME_all_lineage_cluster_policy_v1 import FINE_BY_PARENT
    for parent,names in FINE_BY_PARENT.items():
     for name in names:assessed[name]&=assessed[parent]
   equal(assessed,m['gate_assessed']);csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']]
   for lo in range(0,len(pool),1024):
    pos=np.arange(lo,min(lo+1024,len(pool)));raw=gather(csr,rows[pos]).astype(np.float64).tocsr();raw.sum_duplicates();raw.eliminate_zeros();assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.equal(raw.data,np.floor(raw.data)).all();total=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(total,q['native_total_RNA'][pool[pos]]);x=raw[:,cols].tocsr();cp=x.multiply((10000/total)[:,None]).tocsr();lg=cp.copy();lg.data=np.log1p(lg.data)
    for mapping in [roots[pos],children[pos]]:
     keep=mapping>=0
     if not keep.any():continue
     c=mapping[keep];u,inv=np.unique(c,return_inverse=True);ind=sp.coo_matrix((np.ones(len(c)),(inv,np.arange(len(c)))),shape=(len(u),len(c))).tocsr();stat['n'][u]+=np.bincount(inv,minlength=len(u))
     stat['cp_sum'][u]+=(ind@cp[keep]).toarray();stat['log_sum'][u]+=(ind@lg[keep]).toarray();stat['detected'][u]+=(ind@(x[keep]>0)).toarray()
    local={}
    for name,modules in GATES.items():
     v=np.ones(len(pos),bool)
     for module in modules:
      physical=[by[g] for g in module if g in by]
      if physical:v&=np.asarray((x[:,physical]>0).sum(axis=1)).ravel()>0
      else:v[:]=False
     local[name]=v
    if alternate:
     def observed(g):return np.asarray((x[:,[by[g]]]>0).sum(axis=1)).ravel()>0 if g in by else np.zeros(len(pos),bool)
     auxiliary=np.zeros(len(pos),bool);signal=np.zeros(len(pos),np.int16)
     for g in ['CD2','CD7','TRAT1']:auxiliary|=observed(g)
     for g in ['CD247','LCK','LAT','ZAP70','BCL11B']:signal+=observed(g)
     local['T']=observed('CD3D')&observed('CD3E')&auxiliary&(signal>=2)
    for name in ['Treg','CD4_T','CD8_T','gamma_delta_T_RNA_candidate']:local[name]&=local['T']
    if full:
     from TME_all_lineage_cluster_policy_v1 import FINE_BY_PARENT
     for parent,names in FINE_BY_PARENT.items():
      for name in names:local[name]&=local[parent]
    local['CD4_T']&=~local['Treg']
    for j,k in enumerate(GATE_ORDER):flags[pos,j]=local[k]
   for k in stat:
    if k in ['n','detected']:np.testing.assert_array_equal(stat[k],produced[k])
    else:np.testing.assert_allclose(stat[k],produced[k],rtol=2e-10,atol=1e-7)
   ids=native.native_ids(source);assert text_axis_sha256(ids)==wm['native_ids_source_SHA256'];units=native.frozen_expression_units(source)[positions].astype(str);ng=m['n_roots'];N=int(stat['n'][:ng].sum());sumcp=stat['cp_sum'][:ng].sum(axis=0);sumlg=stat['log_sum'][:ng].sum(axis=0);sumdet=stat['detected'][:ng].sum(axis=0);node_by={r['source_node_code']:r for r in m['nodes']};assert set(node_by)==set(range(nnode));sibling_support={}
   oldmajor=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code'][pool]];oldfine=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code'][pool]]
   for code,r in node_by.items():
    if full:assert r['assay_aware_T_evidence_mode']==('T_chain_unmeasured_alternate_CD3_accessory_two_signal_RNA_evidence' if alternate else 'original_native_module_assay_coverage_in_proposal')
    sel=np.flatnonzero(roots==code) if r['level']==1 else np.flatnonzero(children==code);n=len(sel);assert n==r['n_records']==stat['n'][code] and n>0;parent=r['source_root_cluster'];assert (roots[sel]==parent).all();assert r['stable_membership_SHA256']==membership_key(source,ids[positions[sel]].astype(str).tolist());assert r['native_rows_SHA256']==hashlib.sha256(np.asarray(rows[sel],'<i8').tobytes()).hexdigest()
    call=annotate_cluster(n,dict(zip(GATE_ORDER,flags[sel].sum(axis=0).tolist())),assessed,pc['minimum_named_group_cells'],pc['family_fraction']);equal(call,r['cluster_RNA_proposal']);assert not r['parent_marker_borrowed']
    otherN=N-n if code<ng else int(stat['n'][parent])-n;otherCP=((sumcp if code<ng else stat['cp_sum'][parent])-stat['cp_sum'][code])/otherN if otherN else None;otherLog=((sumlg if code<ng else stat['log_sum'][parent])-stat['log_sum'][code])/otherN if otherN else None;otherDet=((sumdet if code<ng else stat['detected'][parent])-stat['detected'][code])/otherN if otherN else None
    if otherCP is not None:otherCP=np.maximum(otherCP,0)
    markers=independent_markers(genes,stat['cp_sum'][code]/n,stat['detected'][code]/n,otherCP,otherDet);equal(markers,r['marker_candidates']);best=next((x for x in markers if re.fullmatch('[A-Za-z0-9][A-Za-z0-9.-]*',x['gene'])),None);assert r['marker']==(best['gene'] if best else None);equal(state_context(genes,stat['log_sum'][code]/n,otherLog,stat['detected'][code]/n),r['descriptive_state_programs'])
    support=code>=ng and n>=pc['minimum_named_group_cells'] and any(x['descriptive_log2_mean_ratio'] is not None and x['descriptive_log2_mean_ratio']>=pc['minimum_sibling_log2_ratio'] and x['native_detection_fraction']-x['comparator_detection_fraction']>=pc['minimum_sibling_detection_difference'] for x in markers);assert support==r['sibling_marker_separation_supported'];sibling_support[code]=support
    equal(dict(collections.Counter(map(str,oldmajor[sel]))),r['old_v5_major_counts']);equal(dict(collections.Counter(map(str,oldfine[sel]))),r['old_v5_fine_counts']);missing=np.isin(np.char.lower(np.char.strip(units[sel])),['','nan','none','na','n/a','unknown','__missing__','__source_value_missing__']);assert len(set(units[sel][~missing]))==r['n_declared_expression_units_not_verified_patients_or_libraries'];assert int(missing.sum())==r['n_records_missing_declared_expression_unit']
   expectedfinal=roots.copy()
   for p in m['parents']:
    parent=p['source_root_cluster'];child=p['candidate_children'];retained=sum(sibling_support[c] for c in child)>=2;assert retained==p['split_retained'];take=np.flatnonzero(roots==parent)
    if 'model' in p:
     assert sha(p['model']['path'])==p['model']['SHA256']
     with np.load(p['model']['path'],allow_pickle=False) as f:
      np.testing.assert_array_equal(f['pool_positions'],take);np.testing.assert_array_equal(f['native_row'],rows[take]);com=f['components'];emb=f['fit_embedding'];assert com.shape[1]==len(f['native_columns']) and emb.shape==(len(take),len(com));assert set(f['native_columns'].tolist()).issubset(set(cols.tolist()))
      np.testing.assert_allclose(com@com.T,np.eye(len(com)),rtol=0,atol=2e-5);np.testing.assert_allclose(np.square(emb.astype(float)).sum(axis=0)/(len(take)-1),f['explained_variance'],rtol=2e-5,atol=1e-6)
      for curve in p['resolution_sensitivity']:
       key=str(curve['resolution']);a=f['resolution_'+key+'_seed_0'];b=f['resolution_'+key+'_seed_1'];assert float(adjusted_rand_score(a,b))==curve['same_graph_Leiden_seed_ARI'];sizes=np.bincount(a);assert len(sizes)==curve['n_children'] and int(sizes.min())==curve['minimum_child_cells']
      if child:
       label=f['resolution_'+str(pc['primary_resolution'])+'_seed_0'];assert len(child)==len(set(label.tolist()));primary=next(c for c in p['resolution_sensitivity'] if c['resolution']==pc['primary_resolution']);assert primary['n_children_at_least_30']>=2 and primary['same_graph_Leiden_seed_ARI']>=pc['minimum_seed_ARI'];assert sorted(set(map(int,children[take])))==child
       for j,c in enumerate(child):np.testing.assert_array_equal(np.flatnonzero(children==c),take[label==j])
    if retained:expectedfinal[take]=children[take]
   np.testing.assert_array_equal(final,expectedfinal);assert sum(r['n_records'] for r in m['nodes'] if r['final_leaf'])==len(pool)
   for r in m['nodes']:assert r['final_leaf']==(r['source_node_code'] in set(final.tolist()))
   # Raw unit identifiers stay server-side; only anonymous dominance diagnostics.
   dominance=[]
   for r in m['nodes']:
    if not r['final_leaf']:continue
    sel=final==r['source_node_code'];u=units[sel];known=~np.isin(np.char.lower(np.char.strip(u)),['','nan','none','na','n/a','unknown','__missing__','__source_value_missing__']);c=collections.Counter(u[known]);dominance.append(dict(source_node_code=r['source_node_code'],maximum_declared_unit_fraction=max(c.values())/int(sel.sum()) if c else None,n_declared_units_with_at_least_20_cells=sum(v>=20 for v in c.values()),declared_units_not_verified_independent_patients=True))
   op=root/'outputs'/(source+'_cluster_hierarchy_audit_v1.json');save(op,dict(source=source,n_records=len(pool),n_nodes=nnode,n_final_leaves=m['n_final_leaves'],all_full_native_stats_markers_coherent_modules_member_axes_checked=True,model_axes_seed_curves_checked_not_full_PCA_or_graph_independent_reexecution=True,anonymous_unit_dominance=dominance,producer_metadata_SHA256=sha(mp)));checked[source]=dict(source=source,metadata=op.name,metadata_SHA256=sha(op));allnodes+=m['nodes'];todo.remove(source);totals.update(n_sources=1,n_pool_records=len(pool),n_nodes=nnode,n_final_leaves=m['n_final_leaves']);native._fields.clear();del q,stat,produced,flags;gc.collect()
   save(root/'checkpoints'/(str(time.time_ns())+'.json'),dict(observed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),remaining_sources=len(todo)))
  sup=P/'provenance/supervisor_terminal_v1.json'
  if sup.exists() and load(sup)['actual_child_exit_code']!=0:raise RuntimeError('Producer failed; preserve partial evidence')
  if todo and not progress:time.sleep(3)
 while not (P/'provenance/supervisor_terminal_v1.json').exists():time.sleep(3)
 assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;terminal=load(P/'provenance/terminal_v1.json');assert terminal['totals']['n_pool_records']==totals['n_pool_records'];index=load(P/'outputs/cluster_hierarchy_index_v1.json')
 for r in index:assert r['metadata_SHA256']==load(root/'outputs'/checked[r['source']]['metadata'])['producer_metadata_SHA256']
 catalog=load(P/'outputs/cluster_hierarchy_named_catalog_v1.json')['nodes'];history_path=pc['display_registry_path'] if full else w['display_number_registry']['path'];historical=load(history_path)
 if full:assert sha(history_path)==pc['frozen_preflight']['display_registry_SHA256']
 registry=assign_display_numbers([r['stable_membership_SHA256'] for r in allnodes],historical);assert load(P/'outputs/cluster_hierarchy_display_registry_v1.json')==registry['registry'];bykey={r['stable_membership_SHA256']:r for r in allnodes};assert len(catalog)==len(bykey)==len(allnodes)
 for r in catalog:
  expectedrow=dict(bykey[r['stable_membership_SHA256']],global_display_number=registry['current'][r['stable_membership_SHA256']],display_label=display_label(r['cluster_RNA_proposal']['celltype'],r['marker'],registry['current'][r['stable_membership_SHA256']]));equal(expectedrow,r)
 ip=root/'outputs/cluster_hierarchy_audit_index_v1.json';save(ip,[checked[s] for s in sorted(checked)]);save(root/'provenance/terminal_v1.json',dict(status='FULL_CLUSTER_HIERARCHY_NATIVE_MEMBER_MARKER_AND_POLICY_REPLAY_ACCEPTED_PROVISIONAL_VIEW',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),producer_index_SHA256=sha(P/'outputs/cluster_hierarchy_index_v1.json'),index_SHA256=sha(ip),native_CSR_read_and_marker_ranking_implementation_independent=True,cluster_rule_replay_not_orthogonal_accuracy=True,model_audit_scope='all_axes_covariance_and_saved_seed_diagnostics_not_full_PCA_or_knn_reexecution',current_working_cell_labels_unchanged=True,new_CM_started=False,wall_seconds=time.monotonic()-begin))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output);(root/'checkpoints').mkdir()
 try:main(root)
 except BaseException:
  import traceback
  save(root/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
