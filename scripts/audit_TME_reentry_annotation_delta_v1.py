"""Independent direct proposed-subtype counts and accepted-parent remainder proof."""
import argparse,collections,datetime,gc,hashlib,json,re,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def membership(source,ids):
 assert len(ids)==len(set(ids));h=hashlib.sha256()
 for v in ['atlas_exact_membership_v1',source,*sorted(ids)]:
  b=v.encode();h.update(len(b).to_bytes(8,'big'));h.update(b)
 return h.hexdigest()

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');S=Path(cfg['producer']);assert load(S/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;sc=load(S/'config/spec_v1.json');P=Path(sc['proposal_producer']);pc=load(P/'config/spec_v1.json');C=Path(pc['group_profile_producer']);cc=load(C/'config/spec_v1.json');M=Path(cc['model_producer']);mc=load(M/'config/spec_v1.json');R=Path(mc['context_producer']);rc=load(R/'config/spec_v1.json');kc=load(Path(rc['TCR_producer'])/'config/spec_v1.json');E=Path(kc['native_evidence']);EA=Path(kc['native_audit']);J=Path(kc['joined'])
 assert load(Path(sc['proposal_auditor'])/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(P/'provenance/terminal_v1.json');assert load(EA/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(E/'provenance/terminal_v1.json')
 for name,digest in load(S/'config/supervisor_request_v1.json')['script_SHA256'].items():assert sha(S/name)==digest
 w=load(load(rc['working_acceptance'])['producer_manifest']['path']);base=Path(w['index']['path']).parent;wi={r['source']:r for r in load(w['index']['path'])};catalog=load(w['group_catalog']['path'])['groups'];pi={r['source']:r for r in load(P/'outputs/common_reentry_identity_proposal_index_v1.json')};ei={r['source']:r for r in load(E/'outputs/TME_full_native_evidence_index_v1.json')};ci={r['source']:r for r in load(C/'outputs/common_TNK_group_native_profile_index_v1.json')}
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(EA/'scripts'))
 from native_archive_collection_v1 import NativeArchiveCollection
 from audit_full_TME_native_evidence_v1 import gather
 native=NativeArchiveCollection(load(J/'outputs/joined_native_atlas_delivery_v1.json')['archive_manifest']);ip=S/'outputs/reentry_annotation_delta_native_index_v1.json';assert sha(ip)==load(S/'provenance/terminal_v1.json')['index_SHA256'];totals=collections.Counter()
 for entry in load(ip):
  source=entry['source'];mp=S/'outputs'/entry['metadata'];assert sha(mp)==entry['metadata_SHA256'];m=load(mp);assert sha(m['statistics']['path'])==m['statistics']['SHA256']
  with np.load(m['statistics']['path'],allow_pickle=False) as f:actual={k:f[k] for k in f.files}
  pp=P/'outputs'/pi[source]['metadata'];assert sha(pp)==pi[source]['metadata_SHA256'];p=load(pp);assert sha(p['arrays'])==p['arrays_SHA256']
  with np.load(p['arrays'],allow_pickle=False) as f:proposed=f['proposed_T'];typecode=f['proposed_celltype_code'];rows=f['native_row']
  wr=wi[source];wa=base/wr['arrays'];assert sha(wa)==wr['arrays_SHA256'];wm=load(base/wr['metadata'])
  with np.load(wa,allow_pickle=False) as f:q={k:f[k] for k in ['native_row','source_native_position','native_total_RNA','annotation_group_code','consensus_celltype_code','consensus_major_code']}
  np.testing.assert_array_equal(rows,q['native_row']);before=np.asarray(wm['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']];major=np.asarray(wm['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];assert (before[proposed]=='unresolved').all() and (major[proposed]=='unresolved').all();after=before.copy();after[proposed]=np.asarray(p['celltype_dictionary'],object)[typecode[proposed]]
  ep=E/'outputs'/ei[source]['metadata'];assert sha(ep)==ei[source]['metadata_SHA256'];e=load(ep);genes=e['native_measured_gene_names'];cols=np.asarray(e['native_columns']);assert m['native_measured_genes']==genes;np.testing.assert_array_equal(actual['native_columns'],cols);assert sha(e['group_statistics']['path'])==e['group_statistics']['SHA256']
  with np.load(e['group_statistics']['path'],allow_pickle=False) as f:old={k:f[k] for k in ['n','cp_sum','log_sum','detected']}
  oldby={g['stable_group_sha256']:g for g in e['groups']};current={g['source_annotation_group_code']:g for g in catalog if g['source']==source};parents=sorted(set(map(int,q['annotation_group_code'][proposed])));expected_groups=[];codes=np.full(len(proposed),-1,np.int32);reused=0
  for parent in parents:
   domain=q['annotation_group_code']==parent
   for label in sorted(set(after[domain])):
    sel=np.flatnonzero(domain&(after==label));direct=bool(proposed[sel].all());assert direct or not proposed[sel].any();codes[sel]=len(expected_groups);expected_groups.append((parent,str(label),sel,direct,len(sel)==int(domain.sum())));reused+=int(len(sel)==int(domain.sum()))
  take=np.flatnonzero(codes>=0)
  for name,value in [('scope_position',take),('native_row',q['native_row'][take]),('measurement_subset_code',codes[take])]:np.testing.assert_array_equal(actual[name],value)
  n=len(expected_groups);expected={'n':np.zeros(n,np.int64),**{k:np.zeros((n,len(cols)),np.float64) for k in ['cp_sum','log_sum','detected']}};directpos=np.flatnonzero(proposed);csr=native.store(source)[native.legacy_config['counts_key'] if native.is_legacy(source) else native.primary.partitions[source]['counts_key']]
  for lo in range(0,len(directpos),1024):
   sel=directpos[lo:lo+1024];raw=gather(csr,q['native_row'][sel]);rna=np.asarray(raw.sum(axis=1)).ravel();np.testing.assert_array_equal(rna,q['native_total_RNA'][sel]);x=raw[:,cols].tocsr()
   for code in np.unique(codes[sel]):
    sub=x[codes[sel]==code];den=rna[codes[sel]==code];cp=sub.multiply((10000/den)[:,None]).tocsr();log=cp.copy();log.data=np.log1p(log.data);expected['n'][code]+=len(den)
    for name,value in [('cp_sum',cp),('log_sum',log),('detected',sub>0)]:expected[name][code]+=np.asarray(value.sum(axis=0)).ravel()
  for j,(parent,label,sel,direct,reuse) in enumerate(expected_groups):
   if direct:continue
   oldcode=oldby[current[parent]['stable_group_sha256']]['source_annotation_group_code'];children=[k for k,v in enumerate(expected_groups) if v[0]==parent and v[3]]
   for name in expected:
    value=old[name][oldcode]-expected[name][children].sum(axis=0);assert np.min(value)>=(0 if name in ['n','detected'] else -1e-7);expected[name][j]=np.maximum(value,0) if name in ['cp_sum','log_sum'] else value
  for name in ['n','detected']:np.testing.assert_array_equal(actual[name],expected[name])
  for name in ['cp_sum','log_sum']:np.testing.assert_allclose(actual[name],expected[name],rtol=2e-10,atol=1e-7)
  cp=C/'outputs'/ci[source]['metadata'];assert sha(cp)==ci[source]['metadata_SHA256'];common=load(cp);assert common['native_measured_genes']==genes and sha(common['statistics']['path'])==common['statistics']['SHA256']
  with np.load(common['statistics']['path'],allow_pickle=False) as f:poolN=int(f['n'].sum());poolCP=f['cp_sum'].sum(axis=0);poolDet=f['detected'].sum(axis=0);np.testing.assert_array_equal(f['native_columns'],cols)
  ids=native.native_ids(source);assert len(m['groups'])==n
  for j,(parent,label,sel,direct,reuse) in enumerate(expected_groups):
   g=m['groups'][j];pg=current[parent];assert g['measurement_subset_code']==j and g['previous_group_code']==parent and g['previous_group_number']==pg['global_subcluster_number'] and g['previous_group_SHA256']==pg['stable_group_sha256'] and g['fine_celltype']==label and g['consensus_major_proposal']==('T' if direct else 'unresolved') and g['n_exact_subset']==len(sel) and g['exact_membership_reused']==reuse
   assert g['native_profile_derivation']==('direct_proposed_original_CSR' if direct else 'accepted_exact_parent_minus_all_new_proposals');assert g['stable_membership_SHA256']==membership(source,ids[q['source_native_position'][sel]].astype(str).tolist()) and g['native_rows_SHA256']==hashlib.sha256(np.asarray(q['native_row'][sel],'<i8').tobytes()).hexdigest()
   siblings=[k for k,v in enumerate(expected_groups) if v[0]==parent];oldcode=oldby[pg['stable_group_sha256']]['source_annotation_group_code']
   for name in ['n','detected']:np.testing.assert_array_equal(actual[name][siblings].sum(axis=0),old[name][oldcode])
   for name in ['cp_sum','log_sum']:np.testing.assert_allclose(actual[name][siblings].sum(axis=0),old[name][oldcode],rtol=2e-10,atol=1e-7)
   count=poolN-len(sel);own=actual['cp_sum'][j]/len(sel);det=actual['detected'][j]/len(sel);dcp=poolCP-actual['cp_sum'][j];ddet=poolDet-actual['detected'][j];assert count>=0 and dcp.min()>=-1e-7 and ddet.min()>=0;otherCP=np.maximum(dcp,0)/count if count else None;otherDet=ddet/count if count else None;eligible=np.asarray([not re.match(r'^(MT-|RPL|RPS|ENSG|IG[HKL][VDJC])',gene) for gene in genes])&(own>0)&(det>=.1);effect=np.log2((own+.1)/(otherCP+.1)) if count else None
   if count:
    score=effect*np.maximum(det-otherDet,0);order=sorted(np.flatnonzero(eligible&(effect>0)&(det>otherDet)),key=lambda k:(-score[k],genes[k]))[:20];kind='own_subset_expression_specificity_not_inferential_DEG'
   else:order=[];kind='own_subset_expression_without_comparator_not_inferential_DEG'
   if not order:
    order=sorted(np.flatnonzero(eligible),key=lambda k:(-own[k]*det[k],genes[k]))[:20]
    if count:kind='own_subset_expression_without_positive_specificity_not_inferential_DEG'
   candidates=g['own_full_native_gene_marker_candidates_not_DEG'];assert [v['gene'] for v in candidates]==[genes[k] for k in order]
   for k,c in zip(order,candidates):assert c['native_mean_CP10K']==float(own[k]) and c['native_detection_fraction']==float(det[k]) and c['marker_evidence_type']==kind and c['descriptive_log2_mean_ratio']==(float(effect[k]) if count else None) and c['comparator_detection_fraction']==(float(otherDet[k]) if count else None)
   assert g['n_comparator_cells']==count and g['n_source_common_pool_records']==poolN and g['comparator_scope']==('other_current_source_common_TNK_and_unknown_parent_controls' if count else 'no_other_cells_in_current_source_common_pool') and not g['new_identity_applied']
  assert m['proposal_metadata_SHA256']==sha(pp) and m['working_arrays_SHA256']==sha(wa) and m['accepted_parent_native_metadata_SHA256']==sha(ep) and m['accepted_common_pool_metadata_SHA256']==sha(cp)
  totals.update(n_sources=1,n_parent_groups=len(parents),n_exact_subsets=n,n_proposed_native_records_directly_read=len(directpos),n_affected_parent_records=len(take),n_exact_memberships_reused=reused);save(root/'outputs'/(source+'_reentry_annotation_delta_native_audit_v1.json'),dict(source=source,n_direct_proposed_native_records=len(directpos),n_exact_subsets=n,all_subtype_counts_parent_remainders_memberships_and_markers_checked=True));native._fields.clear();del q,old,expected,actual;gc.collect()
 assert dict(totals)==load(S/'provenance/terminal_v1.json')['totals'];save(root/'provenance/terminal_v1.json',dict(status='ALL_REENTRY_ANNOTATION_DELTA_NATIVE_MEMBERSHIPS_AND_MARKERS_ACCEPTED_NOT_APPLIED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),producer_terminal_SHA256=sha(S/'provenance/terminal_v1.json'),producer_index_SHA256=sha(ip),new_production_measurement_and_ranking_helpers_not_imported=True,new_identity_applied=False,working_v5_delta_audit_and_reader_required=True,wall_seconds=time.monotonic()-start))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
