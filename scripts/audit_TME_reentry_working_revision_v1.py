"""Complete working-v5 identity/marker/history delta audit over all sources."""
import argparse,collections,datetime,hashlib,json,re,time
from pathlib import Path
import numpy as np

REVIEW6='additional_TCR_core_reference_local_RNA_coarse_proposal_applied'
REVIEW7='TCR_core_reference_and_common_reentry_context_RNA_identity_applied'
FINE_REVIEW=['not_revised_by_this_direct_fine_node','direct_native_and_unique_reference_fine_RNA','gamma_RNA_candidate_requires_orthogonal_review']
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def hist(a):return dict(collections.Counter(map(str,a)))

def main(root):
 start=time.monotonic();U=Path(load(root/'config/spec_v1.json')['producer']);assert load(U/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;uc=load(U/'config/spec_v1.json');S=Path(uc['native_subsets']);SA=Path(uc['native_subset_audit']);sc=load(S/'config/spec_v1.json');P=Path(sc['proposal_producer']);PA=Path(sc['proposal_auditor'])
 for s in [S,SA,P,PA]:assert load(s/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 for p,a in [(S,SA),(P,PA)]:assert load(a/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(p/'provenance/terminal_v1.json')
 mp=U/'outputs/TME_reentry_working_annotation_manifest_v5.json';assert sha(mp)==load(U/'provenance/terminal_v1.json')['manifest_SHA256'];m=load(mp);assert m['label_version']=='TME_source_qualified_working_annotation_20261027_v5_common_reentry'
 for b in m['bindings'].values():assert sha(b['path'])==b['SHA256']
 old=load(m['bindings']['previous_working_annotation']['path'])
 for v in [old,m]:
  for k in ['index','group_catalog','summary','display_number_registry']:assert sha(v[k]['path'])==v[k]['SHA256']
 oi={r['source']:r for r in load(old['index']['path'])};ni={r['source']:r for r in load(m['index']['path'])};assert set(oi)==set(ni);ob=Path(old['index']['path']).parent;nb=Path(m['index']['path']).parent;oldgroups=load(old['group_catalog']['path'])['groups'];oldby={(g['source'],g['source_annotation_group_code']):g for g in oldgroups};groups=load(m['group_catalog']['path'])['groups'];history=load(old['display_number_registry']['path']);assert all(history[g['stable_group_sha256']]==g['global_subcluster_number'] for g in oldgroups)
 expected=dict(history);nxt=max(history.values())+1
 for k in sorted({g['stable_group_sha256'] for g in groups}-set(history)):expected[k]=nxt;nxt+=1
 assert load(m['display_number_registry']['path'])==expected;sip=S/'outputs/reentry_annotation_delta_native_index_v1.json';assert sha(sip)==load(S/'provenance/terminal_v1.json')['index_SHA256'];subsets={}
 for r in load(sip):
  p=S/'outputs'/r['metadata'];assert sha(p)==r['metadata_SHA256']
  for g in load(p)['groups']:subsets[(r['source'],g['stable_membership_SHA256'])]=(g,sha(p))
 pi={r['source']:r for r in load(P/'outputs/common_reentry_identity_proposal_index_v1.json')};totals=collections.Counter();types=collections.Counter();cross=[]
 for source,oe in oi.items():
  ne=ni[source]
  for entry,base in [(oe,ob),(ne,nb)]:assert sha(base/entry['metadata'])==entry['metadata_SHA256'] and sha(base/entry['arrays'])==entry['arrays_SHA256']
  om=load(ob/oe['metadata']);nm=load(nb/ne['metadata'])
  with np.load(ob/oe['arrays'],allow_pickle=False) as f:o={k:f[k] for k in f.files}
  with np.load(nb/ne['arrays'],allow_pickle=False) as f:n={k:f[k] for k in f.files}
  pp=P/'outputs'/pi[source]['metadata'];assert sha(pp)==pi[source]['metadata_SHA256'];pm=load(pp);assert sha(pm['arrays'])==pm['arrays_SHA256']
  with np.load(pm['arrays'],allow_pickle=False) as f:proposed=f['proposed_T'];typecode=f['proposed_celltype_code'];np.testing.assert_array_equal(f['native_row'],o['native_row'])
  assert set(n)==set(o)|({'fine_revision_review_code'} if proposed.any() else set());allowedarrays={'consensus_major_code','consensus_celltype_code','annotation_review_code','annotation_group_code','fine_revision_review_code'}
  for k in o:
   if k not in allowedarrays:np.testing.assert_array_equal(o[k],n[k])
  dictionaries={k:list(v) for k,v in om['candidate_dictionaries'].items()};names=np.asarray(pm['celltype_dictionary'],object)[typecode];major=o['consensus_major_code'].copy();fine=o['consensus_celltype_code'].copy();review=o['annotation_review_code'].copy();fine_review=o.get('fine_revision_review_code',np.zeros(len(proposed),np.uint8)).copy()
  if proposed.any():
   if 'T' not in dictionaries['consensus_major']:dictionaries['consensus_major'].append('T')
   major[proposed]=dictionaries['consensus_major'].index('T')
   for label in sorted(set(names[proposed])):
    if label not in dictionaries['consensus_celltype']:dictionaries['consensus_celltype'].append(label)
    fine[proposed&(names==label)]=dictionaries['consensus_celltype'].index(label)
   if len(dictionaries['annotation_review'])==6:dictionaries['annotation_review'].append(REVIEW6)
   assert len(dictionaries['annotation_review'])==7 and dictionaries['annotation_review'][6]==REVIEW6;dictionaries['annotation_review'].append(REVIEW7);review[proposed]=7
   if 'fine_revision_review' in dictionaries:assert dictionaries['fine_revision_review']==FINE_REVIEW
   else:dictionaries['fine_revision_review']=FINE_REVIEW.copy()
   fine_review[proposed]=np.where(typecode[proposed]==5,2,np.where(typecode[proposed]>=2,1,0)).astype(np.uint8)
  assert nm['candidate_dictionaries']==dictionaries
  for k,value in [('consensus_major_code',major),('consensus_celltype_code',fine),('annotation_review_code',review)]:np.testing.assert_array_equal(n[k],value)
  if 'fine_revision_review_code' in n:np.testing.assert_array_equal(n['fine_revision_review_code'],fine_review)
  before=np.asarray(om['candidate_dictionaries']['consensus_celltype'],object)[o['consensus_celltype_code']];beforemajor=np.asarray(om['candidate_dictionaries']['consensus_major'],object)[o['consensus_major_code']];assert (before[proposed]=='unresolved').all() and (beforemajor[proposed]=='unresolved').all();after=np.asarray(dictionaries['consensus_celltype'],object)[fine];aftermajor=np.asarray(dictionaries['consensus_major'],object)[major];types.update(after[proposed]);updates=int(proposed.sum());fineupdates=int((proposed&(typecode>=2)).sum())
  assert nm['major_counts']==hist(aftermajor) and nm['celltype_counts']==hist(after) and nm['annotation_review_counts']==hist(review) and nm['fine_revision_review_counts']==hist(fine_review);assert nm['n_reentry_TCR_updates']==updates and nm['n_reentry_fine_updates']==fineupdates;assert nm['identity_state_and_row_arrays_unchanged']==(updates==0) and nm['RNA_state_native_rows_source_gate_and_major_unchanged']==(updates==0) and nm['RNA_states_native_rows_source_gate_and_author_fields_unchanged']
  assert nm['previous_label_version']==om['label_version'] and nm['label_version']==m['label_version'] and nm['group_catalog']==m['group_catalog'];assert nm['arrays']==ne['arrays'] and nm['arrays_SHA256']==ne['arrays_SHA256'];assert nm['old_to_new_celltype_counts']==[dict(previous_working_celltype=a,current_celltype=b,n_records=count) for (a,b),count in sorted(collections.Counter(zip(before,after)).items())]
  allowed={'previous_label_version','label_version','arrays','arrays_SHA256','candidate_dictionaries','n_source_annotation_groups','n_reentry_TCR_updates','n_reentry_fine_updates','major_counts','celltype_counts','annotation_review_counts','fine_revision_review_counts','old_to_new_celltype_counts','group_catalog','identity_state_and_row_arrays_unchanged','RNA_state_native_rows_source_gate_and_major_unchanged','RNA_states_native_rows_source_gate_and_author_fields_unchanged'}
  for k,v in om.items():
   if k not in allowed:assert nm[k]==v,(source,k)
  sg=sorted([g for g in groups if g['source']==source],key=lambda g:g['source_annotation_group_code']);assert [g['source_annotation_group_code'] for g in sg]==list(range(len(sg)))
  for g in sg:
   sel=n['annotation_group_code']==g['source_annotation_group_code'];parents=set(map(int,o['annotation_group_code'][sel]));assert len(parents)==1;parent=parents.pop();og=oldby[(source,parent)];assert int(sel.sum())==g['n_scope_group_rows'];assert g['native_rows_SHA256']==hashlib.sha256(np.asarray(n['native_row'][sel],'<i8').tobytes()).hexdigest();assert g['global_subcluster_number']==expected[g['stable_group_sha256']] and g['global_subcluster_id']=='TME_annotation_group_'+str(g['global_subcluster_number'])
   if (source,g['stable_group_sha256']) in subsets:
    s,h=subsets[(source,g['stable_group_sha256'])];np.testing.assert_array_equal(sel,(o['annotation_group_code']==parent)&(after==s['fine_celltype']));assert s['previous_group_SHA256']==og['stable_group_sha256'] and s['n_exact_subset']==int(sel.sum()) and s['native_rows_SHA256']==g['native_rows_SHA256'] and g['consensus_celltype']==s['fine_celltype'] and g['major_distribution']==hist(aftermajor[sel]);assert g['previous_reentry_group_SHA256']==og['stable_group_sha256'] and g['previous_reentry_group_number']==og['global_subcluster_number']
    choices=[v for v in s['own_full_native_gene_marker_candidates_not_DEG'] if re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9.-]*',v['gene']) and v['gene'].lower() not in ['nan','none','unknown']];best=choices[0] if choices else None;assert g['marker']==(best['gene'] if best else None) and g['marker_origin']==(g['stable_group_sha256'] if best else 'unavailable') and g['marker_evidence_type']==(best['marker_evidence_type'] if best else 'no_accepted_current_exact_subset_marker') and g['marker_expression_evidence_metadata_SHA256']==h
    for field,key in [('marker_native_mean_CP10K','native_mean_CP10K'),('marker_detection_fraction','native_detection_fraction'),('marker_descriptive_log2_mean_ratio','descriptive_log2_mean_ratio')]:assert g[field]==(best[key] if best else None)
    assert g['marker_comparator_scope']==s['comparator_scope'] and g['marker_comparator_cell_count']==s['n_comparator_cells'] and g['marker_evidence_membership_matches_exact_group'] and g['marker_DEG_status']=='no_qualified_replicate_aware_DEG_supplied';relation='same_membership_identity_revision' if s['exact_membership_reused'] else 'split'
    if s['exact_membership_reused']:assert g['global_subcluster_number']==og['global_subcluster_number']
    else:assert g['global_subcluster_number']>max(history.values())
    for k,v in og.items():
     if k=='marker' or k.startswith('marker_') or k in ['parent_marker_borrowed','borrow_reason','subcluster']:assert g['previous_reentry__'+k]==v
     elif k not in {'source_annotation_group_code','stable_group_sha256','n_scope_group_rows','consensus_celltype','major_distribution','native_rows_SHA256','membership_identical_to_archive_parent','global_subcluster_number','global_subcluster_id'}:assert g[k]==v,(source,k)
   else:assert g==og;np.testing.assert_array_equal(sel,o['annotation_group_code']==parent);relation='same_membership'
   assert not g['parent_marker_borrowed'] and g['subcluster']==g['consensus_celltype']+'_'+(g['marker'] or 'marker_unresolved')+'_'+str(g['global_subcluster_number']);cross.append(dict(source=source,old_group_SHA256=og['stable_group_sha256'],old_group_number=og['global_subcluster_number'],old_display=og['subcluster'],new_group_SHA256=g['stable_group_sha256'],n_records=int(sel.sum()),relation=relation,new_celltype=g['consensus_celltype'],new_group_number=g['global_subcluster_number'],new_display=g['subcluster']))
  assert nm['n_source_annotation_groups']==len(sg);totals.update(n_sources=1,n_source_qualified_rows=len(proposed),n_fit_candidates_not_final_QC=int(n['fit_pool_candidate_not_final_QC'].sum()),n_without_fitted_group=int((~n['fit_pool_candidate_not_final_QC']).sum()),n_TCR_coarse_updates=nm.get('n_TCR_coarse_updates',0),n_direct_fine_updates=nm.get('n_direct_fine_updates',0),n_reentry_TCR_updates=updates,n_reentry_fine_updates=fineupdates,n_exact_scoped_groups=len(sg))
 xp=m['annotation_membership_crosswalk'];assert sha(xp['path'])==xp['SHA256'];actual=load(xp['path']);order=lambda r:(r['source'],r['new_group_SHA256']);assert sorted(actual,key=order)==sorted(cross,key=order);assert sum(r['n_records'] for r in actual)==totals['n_fit_candidates_not_final_QC'];assert load(m['summary']['path'])['totals']==dict(totals)==load(U/'provenance/terminal_v1.json')['totals'];assert dict(types)==load(P/'provenance/terminal_v1.json')['proposed_celltype_counts']
 save(root/'provenance/terminal_v1.json',dict(status='ALL_TME_COMMON_REENTRY_WORKING_IDENTITIES_MARKERS_HISTORY_AND_NATIVE_INVARIANTS_ACCEPTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),new_celltype_counts=dict(types),producer_manifest_SHA256=sha(mp),producer_terminal_SHA256=sha(U/'provenance/terminal_v1.json'),all_source_identity_deltas_and_crosswalk_checked=True,all_inactive_historical_numbers_preserved=True,not_biological_accuracy_or_final_taxonomy_acceptance=True,new_CM_started=False,final_module_input_ready=False,wall_seconds=time.monotonic()-start))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',dict(traceback=traceback.format_exc()));raise
