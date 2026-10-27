"""Independent candidate CSR replay and accepted-parent remainder proof."""
import argparse,collections,datetime,gc,hashlib,json,re,resource,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def key(source,ids):
 assert len(ids)==len(set(ids));h=hashlib.sha256()
 for v in ['atlas_exact_membership_v1',source,*sorted(ids)]:
  b=v.encode();h.update(len(b).to_bytes(8,'big'));h.update(b)
 return h.hexdigest()

def main(root):
 start=time.monotonic();P=Path(load(root/'config/spec_v1.json')['producer']);assert load(P/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 cfg=load(P/'config/spec_v1.json');R=Path(cfg['context_producer']);rc=load(R/'config/spec_v1.json');pc=load(Path(rc['TCR_producer'])/'config/spec_v1.json');E=Path(pc['native_evidence']);EA=Path(pc['native_audit']);J=Path(pc['joined']);assert load(EA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(E/'provenance/terminal_v1.json')
 for name,digest in load(P/'config/supervisor_request_v1.json')['script_SHA256'].items():assert sha(P/name)==digest
 assert load(Path(cfg['context_auditor'])/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(R/'provenance/terminal_v1.json')
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(EA/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from audit_full_TME_native_evidence_v1 import gather
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);w=load(load(rc['working_acceptance'])['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};catalog=load(w['group_catalog']['path'])['groups'];ri={r['source']:r for r in load(R/'outputs/unknown_parent_TCR_context_index_v1.json')};ei={r['source']:r for r in load(E/'outputs/TME_full_native_evidence_index_v1.json')}
 ip=P/'outputs/unknown_parent_TCR_profile_index_v1.json';assert sha(ip)==load(P/'provenance/terminal_v1.json')['index_SHA256'];index=load(ip);totals=collections.Counter()
 for entry in index:
  source=entry['source'];mp=P/'outputs'/entry['metadata'];assert sha(mp)==entry['metadata_SHA256'];m=load(mp);sp=m['statistics'];assert sha(sp['path'])==sp['SHA256']
  with np.load(sp['path'],allow_pickle=False) as f:actual={k:f[k] for k in f.files}
  wr=wi[source];wa=base/wr['arrays'];assert sha(wa)==wr['arrays_SHA256']
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','annotation_group_code','native_total_RNA']}
  rp=R/'outputs'/ri[source]['metadata'];assert sha(rp)==ri[source]['metadata_SHA256'];review=load(rp);assert sha(review['arrays'])==review['arrays_SHA256']
  with np.load(review['arrays'],allow_pickle=False) as f:target=f['target_TCR_unknown_parent'];np.testing.assert_array_equal(f['native_row'],q['native_row'])
  ep=E/'outputs'/ei[source]['metadata'];assert sha(ep)==ei[source]['metadata_SHA256'];e=load(ep);genes=e['native_measured_gene_names'];cols=np.asarray(e['native_columns']);assert genes==m['native_measured_genes'];np.testing.assert_array_equal(actual['native_columns'],cols);assert sha(e['group_statistics']['path'])==e['group_statistics']['SHA256']
  with np.load(e['group_statistics']['path'],allow_pickle=False) as f:old={k:f[k] for k in ['n','cp_sum','log_sum','detected']}
  oldby={g['stable_group_sha256']:g for g in e['groups']};current={g['source_annotation_group_code']:g for g in catalog if g['source']==source};parents=sorted(set(map(int,q['annotation_group_code'][target])));expected_groups=[];codes=np.full(len(target),-1,np.int32);candidate_code={}
  for parent in parents:
   for role,sel in [('TCR_candidate_not_reclassified',(q['annotation_group_code']==parent)&target),('unresolved_remainder',(q['annotation_group_code']==parent)&~target)]:
    if not sel.any():continue
    code=len(expected_groups);codes[sel]=code;expected_groups.append((parent,role,np.flatnonzero(sel)))
    if role=='TCR_candidate_not_reclassified':candidate_code[parent]=code
  take=np.flatnonzero(codes>=0)
  for name,v in [('scope_position',take),('native_row',q['native_row'][take]),('measurement_subset_code',codes[take])]:np.testing.assert_array_equal(actual[name],v)
  n=len(expected_groups);expected={'n':np.zeros(n,np.int64),**{k:np.zeros((n,len(cols))) for k in ['cp_sum','log_sum','detected']}};tpos=np.flatnonzero(target);csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']]
  for lo in range(0,len(tpos),1024):
   use=tpos[lo:lo+1024];raw=gather(csr,q['native_row'][use]);rna=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(rna,q['native_total_RNA'][use]);assert (rna>0).all();x=raw[:,cols].tocsr()
   for code in np.unique(codes[use]):
    sub=x[codes[use]==code];den=rna[codes[use]==code];cp=sub.multiply((10000/den)[:,None]).tocsr();logged=cp.copy();logged.data=np.log1p(logged.data);expected['n'][code]+=len(den)
    for name,values in [('cp_sum',cp),('log_sum',logged),('detected',sub>0)]:expected[name][code]+=np.asarray(values.sum(axis=0)).ravel()
  for j,(parent,role,sel) in enumerate(expected_groups):
   if role!='unresolved_remainder':continue
   parent_code=oldby[current[parent]['stable_group_sha256']]['source_annotation_group_code'];own=candidate_code[parent]
   for name in expected:
    values=old[name][parent_code]-expected[name][own];assert np.min(values)>=(0 if name in ['n','detected'] else -1e-7);expected[name][j]=np.maximum(values,0) if name in ['cp_sum','log_sum'] else values
  for name in ['n','detected']:np.testing.assert_array_equal(actual[name],expected[name])
  for name in ['cp_sum','log_sum']:np.testing.assert_allclose(actual[name],expected[name],rtol=2e-10,atol=1e-7)
  assert len(m['groups'])==n;ids=native.native_ids(source);units=native.frozen_expression_units(source)[q['source_native_position']].astype(str)
  for j,(parent,role,sel) in enumerate(expected_groups):
   g=m['groups'][j];pg=current[parent];oldcode=oldby[pg['stable_group_sha256']]['source_annotation_group_code'];assert g['measurement_subset_code']==j and g['previous_group_code']==parent and g['previous_group_number']==pg['global_subcluster_number'] and g['previous_group_SHA256']==pg['stable_group_sha256'] and g['native_parent_profile_row']==oldcode and g['role']==role and g['n_exact_subset']==len(sel)
   assert g['stable_membership_SHA256']==key(source,ids[q['source_native_position'][sel]].astype(str).tolist()) and g['native_rows_SHA256']==hashlib.sha256(np.asarray(q['native_row'][sel],'<i8').tobytes()).hexdigest();missing=np.isin(np.char.lower(np.char.strip(units[sel])),['','nan','none','na','n/a','unknown','__missing__','__source_value_missing__']);assert g['n_declared_expression_units_not_libraries_or_patients']==len(set(units[sel][~missing])) and g['n_records_with_missing_declared_unit']==int(missing.sum())
   siblings=[k for k,c in enumerate(expected_groups) if c[0]==parent];other=[k for k in siblings if k!=j]
   for name in ['n','detected']:np.testing.assert_array_equal(actual[name][siblings].sum(axis=0),old[name][oldcode])
   for name in ['cp_sum','log_sum']:np.testing.assert_allclose(actual[name][siblings].sum(axis=0),old[name][oldcode],rtol=2e-10,atol=1e-7)
   count=int(actual['n'][other].sum());own=actual['cp_sum'][j]/len(sel);det=actual['detected'][j]/len(sel);othercp=actual['cp_sum'][other].sum(axis=0)/count if count else None;otherdet=actual['detected'][other].sum(axis=0)/count if count else None;eligible=np.asarray([not re.match(r'^(MT-|RPL|RPS|ENSG|IG[HKL][VDJC])',g) for g in genes])&(own>0)&(det>=.1);effect=np.log2((own+.1)/(othercp+.1)) if count else None
   if count:
    score=effect*np.maximum(det-otherdet,0);order=sorted(np.flatnonzero(eligible&(effect>0)&(det>otherdet)),key=lambda k:(-score[k],genes[k]))[:20];kind='own_subset_expression_specificity_not_inferential_DEG'
   else:order=[];kind='own_subset_expression_without_comparator_not_inferential_DEG'
   if not order:
    order=sorted(np.flatnonzero(eligible),key=lambda k:(-own[k]*det[k],genes[k]))[:20]
    if count:kind='own_subset_expression_without_positive_specificity_not_inferential_DEG'
   candidates=g['own_full_native_gene_marker_candidates_not_DEG'];assert [v['gene'] for v in candidates]==[genes[k] for k in order]
   for k,c in zip(order,candidates):
    assert c['native_mean_CP10K']==float(own[k]) and c['native_detection_fraction']==float(det[k]) and c['marker_evidence_type']==kind and c['comparator_detection_fraction']==(float(otherdet[k]) if count else None)
    assert c['descriptive_log2_mean_ratio']==(float(effect[k]) if count else None)
   assert g['n_comparator_cells']==count and g['comparator_scope']==('other_role_in_same_exact_current_unresolved_parent' if count else 'no_nonempty_parent_remainder_or_counterpart') and not g['new_identity_assigned']
  assert m['context_metadata_SHA256']==sha(rp) and m['working_arrays_SHA256']==sha(wa) and m['accepted_complete_native_metadata_SHA256']==sha(ep)
  totals.update(n_sources=1,n_parent_groups=len(parents),n_exact_subsets=n,n_candidate_native_rows_read=len(tpos),n_parent_scope_records=len(take));save(root/'outputs'/(source+'_unknown_parent_TCR_profile_audit_v1.json'),dict(source=source,all_candidate_native_rows_directly_reread=True,all_parent_remainder_keys_markers_and_units_checked=True));native._fields.clear();del q,actual,expected,old;gc.collect()
 assert dict(totals)==load(P/'provenance/terminal_v1.json')['totals']
 save(root/'provenance/terminal_v1.json',dict(status='ALL_UNKNOWN_PARENT_TCR_NATIVE_PROFILES_AND_REMAINDER_PROOFS_ACCEPTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),producer_terminal_SHA256=sha(P/'provenance/terminal_v1.json'),producer_index_SHA256=sha(ip),new_production_measurement_and_ranking_not_imported=True,new_identity_assigned=False,not_orthogonal_biological_accuracy=True,wall_seconds=time.monotonic()-start,peak_RSS_bytes_Linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
