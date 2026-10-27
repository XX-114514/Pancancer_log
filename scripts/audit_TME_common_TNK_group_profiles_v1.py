"""Complete group-expression proof by accepted exact unions and direct splits.

Every model-pool record is accounted for. Whole, unchanged accepted memberships
can contribute their already audited additive native moments; all other selected
records are independently reread. This is not sampling or biological validation.
"""
import argparse,collections,datetime,gc,hashlib,json,re,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def hist(a):return dict(collections.Counter(map(str,a)))
def membership(source,ids):
 assert len(ids)==len(set(ids));h=hashlib.sha256()
 for v in ['atlas_exact_membership_v1',source,*sorted(ids)]:
  b=v.encode();h.update(len(b).to_bytes(8,'big'));h.update(b)
 return h.hexdigest()
def live(root):
 p=root/'provenance/child_launch_v1.json'
 if not p.exists():return False
 c=load(p);s=Path('/proc')/str(c['PID'])/'stat'
 return s.exists() and s.read_text().rsplit(')',1)[1].split()[19]==c['startticks']

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');S=Path(cfg['producer']);NA=Path(cfg['model_native_auditor']);sc=load(S/'config/spec_v1.json');P=Path(sc['model_producer']);pc=load(P/'config/spec_v1.json');R=Path(pc['context_producer']);rc=load(R/'config/spec_v1.json');K=Path(rc['TCR_producer']);kc=load(K/'config/spec_v1.json');J=Path(kc['joined']);E=Path(kc['native_evidence']);EA=Path(kc['native_audit']);F=Path(kc['candidate_producer']);w=load(load(rc['working_acceptance'])['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};catalog=load(w['group_catalog']['path'])['groups']
 for name,digest in load(S/'config/supervisor_request_v1.json')['script_SHA256'].items():assert sha(S/name)==digest
 assert load(EA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(E/'provenance/terminal_v1.json')
 roots=[]
 for idxname,auditname in [('changed_subset_index','changed_subset_audit'),('fine_delta_native_index','fine_delta_native_audit')]:
  b=w['bindings'][idxname];ip=Path(b['path']);assert sha(ip)==b['SHA256'];b=w['bindings'][auditname];ap=Path(b['path']);assert sha(ap)==b['SHA256'];rp=ip.parent.parent;assert load(ap)['producer_terminal_SHA256']==sha(rp/'provenance/terminal_v1.json');roots.append((rp,{r['source']:r for r in load(ip)}))
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(EA/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from audit_full_TME_native_evidence_v1 import gather
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);ei={r['source']:r for r in load(E/'outputs/TME_full_native_evidence_index_v1.json')};fi={r['source']:r for r in load(F/'outputs/full_native_cell_identity_candidate_index_v1.json')};ki={r['source']:r for r in load(K/'outputs/TCR_boundary_proposal_index_v1.json')};ri={r['source']:r for r in load(R/'outputs/unknown_parent_TCR_context_index_v1.json')};pi={r['source']:r for r in load(P/'outputs/common_TNK_reentry_index_v1.json')};done={};totals=collections.Counter()
 while set(done)!=set(wi):
  for source in sorted(set(wi)-set(done)):
   mp=S/'outputs'/(source+'_common_TNK_group_native_profiles_v1.json')
   if not mp.exists():continue
   try:m=load(mp)
   except json.JSONDecodeError:continue
   assert sha(m['statistics']['path'])==m['statistics']['SHA256']
   with np.load(m['statistics']['path'],allow_pickle=False) as f:actual={k:f[k] for k in f.files}
   wr=wi[source];wa=base/wr['arrays'];assert sha(wa)==wr['arrays_SHA256'];wm=load(base/wr['metadata'])
   with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','annotation_group_code','native_total_RNA','consensus_major_code','consensus_celltype_code','native_family_bitmask','reference_family_bitmask','RNA_state_native_logCP10K']}
   modelp=P/'outputs'/pi[source]['metadata'];assert sha(modelp)==pi[source]['metadata_SHA256'];model=load(modelp);assert sha(model['arrays'])==model['arrays_SHA256']
   with np.load(model['arrays'],allow_pickle=False) as f:pool=f['source_scope_position'];codes=f['reentry_group_code'];rows=f['native_row'];known=f['baseline_known_T_NK'];control=f['unresolved_parent_control']
   np.testing.assert_array_equal(rows,q['native_row'][pool]);n=len(pool);ng=model['n_groups'];ep=E/'outputs'/ei[source]['metadata'];assert sha(ep)==ei[source]['metadata_SHA256'];e=load(ep);genes=e['native_measured_gene_names'];cols=np.asarray(e['native_columns']);assert genes==m['native_measured_genes'];np.testing.assert_array_equal(actual['native_columns'],cols)
   p=e['group_statistics'];assert sha(p['path'])==p['SHA256']
   with np.load(p['path'],allow_pickle=False) as f:old={k:f[k] for k in ['n','cp_sum','log_sum','detected']}
   accepted={g['stable_group_sha256']:(old,g['source_annotation_group_code']) for g in e['groups']}
   for rp,index in roots:
    if source not in index:continue
    r=index[source];sp=rp/'outputs'/r['metadata'];assert sha(sp)==r['metadata_SHA256'];sm=load(sp);assert sm['native_measured_genes']==genes and sha(sm['statistics']['path'])==sm['statistics']['SHA256']
    with np.load(sm['statistics']['path'],allow_pickle=False) as f:a={k:f[k] for k in ['n','cp_sum','log_sum','detected','native_columns']}
    np.testing.assert_array_equal(a['native_columns'],cols)
    for g in sm['groups']:accepted[g['stable_membership_SHA256']]=(a,g['measurement_subset_code'])
   expected={'n':np.zeros(ng,np.int64),**{k:np.zeros((ng,len(cols)),np.float64) for k in ['cp_sum','log_sum','detected']}};current={g['source_annotation_group_code']:g for g in catalog if g['source']==source};oldcodes=q['annotation_group_code'][pool];ids=native.native_ids(source);direct=np.ones(n,bool);union_records=union_groups=0
   for code in np.unique(oldcodes):
    assert code>=0;g=current[int(code)];sel=np.flatnonzero(oldcodes==code);dest=np.unique(codes[sel])
    if len(sel)!=g['n_scope_group_rows'] or len(dest)!=1:continue
    assert membership(source,ids[q['source_native_position'][pool[sel]]].astype(str).tolist())==g['stable_group_sha256'];a,j=accepted[g['stable_group_sha256']];assert a['n'][j]==len(sel)
    for name in expected:expected[name][int(dest[0])]+=a[name][j]
    direct[sel]=False;union_records+=len(sel);union_groups+=1
   take=np.flatnonzero(direct);csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']]
   for lo in range(0,len(take),1024):
    use=take[lo:lo+1024];raw=gather(csr,rows[use]);rna=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(rna,q['native_total_RNA'][pool[use]]);x=raw[:,cols].tocsr()
    for code in np.unique(codes[use]):
     sub=x[codes[use]==code];den=rna[codes[use]==code];cp=sub.multiply((10000/den)[:,None]).tocsr();logged=cp.copy();logged.data=np.log1p(logged.data);expected['n'][code]+=len(den)
     for name,value in [('cp_sum',cp),('log_sum',logged),('detected',sub>0)]:expected[name][code]+=np.asarray(value.sum(axis=0)).ravel()
   assert union_records+len(take)==n
   for name in ['n','detected']:np.testing.assert_array_equal(actual[name],expected[name])
   for name in ['cp_sum','log_sum']:np.testing.assert_allclose(actual[name],expected[name],rtol=2e-10,atol=1e-7)
   fr=fi[source];fp=F/'outputs'/fr['metadata'];fa=F/'outputs'/fr['arrays'];assert sha(fp)==fr['metadata_sha256'] and sha(fa)==fr['arrays_sha256'];fm=load(fp)
   with np.load(fa,allow_pickle=False) as f:pos=q['source_native_position'][pool];np.testing.assert_array_equal(f['native_row'][pos],rows);directfine=f['own_fine_bitmask'][pos]&f['reference_fine_bitmask'][pos]
   directfine[(directfine&1)!=0]&=np.asarray(65533,dtype=directfine.dtype);kp=K/'outputs'/ki[source]['metadata'];assert sha(kp)==ki[source]['metadata_SHA256'];km=load(kp);ka=K/'outputs'/km['arrays'];assert sha(ka)==km['arrays_SHA256']
   with np.load(ka,allow_pickle=False) as f:np.testing.assert_array_equal(f['native_row'],q['native_row']);tcore=f['TCR_core'][pool];comp=f['competing_core_mask'][pool];unknown=f['unassessed_old_competitor_mask'][pool];low=f['reference_low_T'][pool]
   rp=R/'outputs'/ri[source]['metadata'];assert sha(rp)==ri[source]['metadata_SHA256'];rm=load(rp);assert sha(rm['arrays'])==rm['arrays_SHA256']
   with np.load(rm['arrays'],allow_pickle=False) as f:target=f['target_TCR_unknown_parent'][pool]
   major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code'][pool]];fine=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code'][pool]];units=native.frozen_expression_units(source)[q['source_native_position'][pool]].astype(str);assert len(m['groups'])==ng
   for code,g in enumerate(m['groups']):
    sel=np.flatnonzero(codes==code);assert g['source_reentry_group_code']==code and g['n_exact_group']==len(sel);assert g['stable_membership_SHA256']==membership(source,ids[q['source_native_position'][pool[sel]]].astype(str).tolist()) and g['native_rows_SHA256']==hashlib.sha256(np.asarray(rows[sel],'<i8').tobytes()).hexdigest()
    for name,value in [('n_baseline_known_T_NK',known[sel]),('n_unknown_parent_control',control[sel]),('n_strict_TCR_unknown_parent_candidate',target[sel]),('n_TCR_core',tcore[sel]),('n_TCR_core_with_coherent_competitor',tcore[sel]&(comp[sel]!=0)),('n_unique_low_reference_T',low[sel])]:assert g[name]==int(value.sum())
    for name,value in [('coherent_competitor_mask_counts',comp[sel]),('unassessed_old_competitor_mask_counts',unknown[sel]),('baseline_major_counts',major[sel]),('baseline_celltype_counts',fine[sel]),('native_family_mask_counts',q['native_family_bitmask'][pool[sel]]),('reference_family_mask_counts',q['reference_family_bitmask'][pool[sel]])]:assert g[name]==hist(value)
    assert g['direct_own_and_reference_fine_cue_counts']=={name:int(((directfine[sel]&(1<<j))!=0).sum()) for j,name in enumerate(fm['fine_candidate_order'])};missing=np.isin(np.char.lower(np.char.strip(units[sel])),['','nan','none','na','n/a','unknown','__missing__','__source_value_missing__']);assert g['n_declared_expression_units_not_libraries_or_patients']==len(set(units[sel][~missing])) and g['n_records_with_missing_declared_unit']==int(missing.sum())
    for j,name in enumerate(wm['state_order']):
     values=q['RNA_state_native_logCP10K'][pool[sel],j];ok=np.isfinite(values);assert g['descriptive_existing_RNA_states'][name]==dict(n_measured=int(ok.sum()),mean_native_logCP10K=float(values[ok].mean(dtype=np.float64)) if ok.any() else None,measurement_coverage=wm['state_measurement_coverage'][name],functional_identity_assigned=False)
    others=np.arange(ng)!=code;count=int(actual['n'][others].sum());own=actual['cp_sum'][code]/len(sel);det=actual['detected'][code]/len(sel);othercp=actual['cp_sum'][others].sum(axis=0)/count if count else None;otherdet=actual['detected'][others].sum(axis=0)/count if count else None;effect=np.log2((own+.1)/(othercp+.1)) if count else None;eligible=np.asarray([not re.match(r'^(MT-|RPL|RPS|ENSG|IG[HKL][VDJC])',name) for name in genes])&(own>0)&(det>=.1)
    if count:
     score=effect*np.maximum(det-otherdet,0);order=sorted(np.flatnonzero(eligible&(effect>0)&(det>otherdet)),key=lambda j:(-score[j],genes[j]))[:20];kind='own_subset_expression_specificity_not_inferential_DEG'
    else:order=[];kind='own_subset_expression_without_comparator_not_inferential_DEG'
    if not order:
     order=sorted(np.flatnonzero(eligible),key=lambda j:(-own[j]*det[j],genes[j]))[:20]
     if count:kind='own_subset_expression_without_positive_specificity_not_inferential_DEG'
    cand=g['own_full_native_gene_marker_candidates_not_DEG'];assert [v['gene'] for v in cand]==[genes[j] for j in order]
    for j,c in zip(order,cand):
     assert c['native_mean_CP10K']==float(own[j]) and c['native_detection_fraction']==float(det[j]) and c['marker_evidence_type']==kind
     if count:np.testing.assert_allclose(c['descriptive_log2_mean_ratio'],float(effect[j]),rtol=1e-10,atol=1e-10);np.testing.assert_allclose(c['comparator_detection_fraction'],float(otherdet[j]),rtol=1e-12,atol=1e-12)
     else:assert c['descriptive_log2_mean_ratio'] is None and c['comparator_detection_fraction'] is None
    assert g['n_comparator_cells']==count and g['comparator_scope']==('other_current_source_common_TNK_and_unknown_control_groups' if count else 'no_other_group_in_current_source_common_pool') and not g['new_identity_assigned']
   assert m['model_metadata_SHA256']==sha(modelp) and m['working_arrays_SHA256']==sha(wa) and m['complete_native_metadata_SHA256']==sha(ep) and m['TCR_metadata_SHA256']==sha(kp) and m['direct_fine_metadata_SHA256']==sha(fp)
   out=root/'outputs'/(source+'_common_TNK_group_profile_audit_v1.json');save(out,dict(source=source,n_model_pool_records=n,n_groups=ng,n_direct_native_reread=len(take),n_proven_by_accepted_exact_unions=union_records,n_accepted_old_groups_reused=union_groups,all_records_and_native_statistics_and_context_and_markers_checked=True));done[source]=dict(source=source,metadata=out.name,metadata_SHA256=sha(out),producer_metadata_SHA256=sha(mp));totals.update(n_sources=1,n_model_pool_records=n,n_exact_groups=ng,n_direct_native_reread=len(take),n_proven_by_accepted_exact_unions=union_records,n_accepted_old_groups_reused=union_groups);native._fields.clear();del q,actual,expected,accepted;gc.collect()
  for stage in [S,NA]:
   p=stage/'provenance/supervisor_terminal_v1.json'
   if p.exists():assert load(p)['actual_child_exit_code']==0,'Upstream failed; completed source proofs retained'
   elif not live(stage):raise RuntimeError('Exact upstream process missing')
  if set(done)!=set(wi):
   if time.monotonic()-start>cfg['maximum_wait_seconds']:raise TimeoutError('Bounded source wait exhausted')
   time.sleep(15)
 for stage in [S,NA]:
  while not (stage/'provenance/supervisor_terminal_v1.json').exists():
   assert live(stage)
   if time.monotonic()-start>cfg['maximum_wait_seconds']:raise TimeoutError('Bounded terminal wait exhausted')
   time.sleep(2)
  assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 st=load(S/'provenance/terminal_v1.json');ip=S/'outputs/common_TNK_group_native_profile_index_v1.json';assert sha(ip)==st['index_SHA256'];assert set(done)=={r['source'] for r in load(ip)}
 for r in load(ip):assert done[r['source']]['producer_metadata_SHA256']==r['metadata_SHA256']
 assert totals['n_model_pool_records']==st['totals']['n_model_pool_records']==totals['n_direct_native_reread']+totals['n_proven_by_accepted_exact_unions'];assert totals['n_exact_groups']==st['totals']['n_exact_groups'];op=root/'outputs/common_TNK_group_profile_audit_index_v1.json';save(op,[done[s] for s in sorted(done)])
 save(root/'provenance/terminal_v1.json',dict(status='ALL_COMMON_TNK_GROUP_NATIVE_EXACT_UNION_AND_DIRECT_SPLIT_PROOFS_ACCEPTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),producer_terminal_SHA256=sha(S/'provenance/terminal_v1.json'),model_native_audit_terminal_SHA256=sha(NA/'provenance/terminal_v1.json'),index_SHA256=sha(op),all_model_pool_records_accounted_for_not_sampling=True,production_measurement_and_marker_ranking_helpers_not_imported=True,biological_identity_and_graph_state_preservation_not_certified=True,new_identity_assigned=False,wall_seconds=time.monotonic()-start))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
