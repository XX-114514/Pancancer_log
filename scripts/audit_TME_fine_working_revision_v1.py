"""All-source fine-only delta audit over accepted proposals and native profiles."""
import argparse,collections,datetime,hashlib,json,re,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')

def main(root):
 start=time.monotonic();U=Path(load(root/'config/spec_v1.json')['producer']);assert load(U/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 uc=load(U/'config/spec_v1.json');S=Path(uc['native_subsets']);SA=Path(uc['native_subset_audit']);sc=load(S/'config/spec_v1.json');P=Path(sc['proposal']);PA=Path(sc['proposal_audit'])
 for stage in [S,SA,P,PA]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 for producer,auditor in [(S,SA),(P,PA)]:assert load(auditor/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(producer/'provenance/terminal_v1.json')
 mp=U/'outputs/TME_direct_fine_working_annotation_manifest_v4.json';assert sha(mp)==load(U/'provenance/terminal_v1.json')['manifest_SHA256'];m=load(mp)
 for b in m['bindings'].values():assert sha(b['path'])==b['SHA256']
 old=load(m['bindings']['previous_working_annotation']['path'])
 for c in [old,m]:
  for k in ['index','group_catalog','summary','display_number_registry']:assert sha(c[k]['path'])==c[k]['SHA256']
 oi={r['source']:r for r in load(old['index']['path'])};ni={r['source']:r for r in load(m['index']['path'])};assert set(oi)==set(ni);ob=Path(old['index']['path']).parent;nb=Path(m['index']['path']).parent
 oldgroups=load(old['group_catalog']['path'])['groups'];oldby={(g['source'],g['source_annotation_group_code']):g for g in oldgroups};groups=load(m['group_catalog']['path'])['groups']
 history=load(old['display_number_registry']['path']);assert all(history[g['stable_group_sha256']]==g['global_subcluster_number'] for g in oldgroups);expected=dict(history);nextnumber=max(history.values())+1
 for key in sorted({g['stable_group_sha256'] for g in groups}-set(history)):expected[key]=nextnumber;nextnumber+=1
 assert load(m['display_number_registry']['path'])==expected
 sip=S/'outputs/fine_delta_native_profile_index_v1.json';assert sha(sip)==load(S/'provenance/terminal_v1.json')['index_SHA256'];subsets={}
 for r in load(sip):
  p=S/'outputs'/r['metadata'];assert sha(p)==r['metadata_SHA256'];sm=load(p)
  for g in sm['groups']:subsets[(r['source'],g['stable_membership_SHA256'])]=(g,sha(p))
 pi={r['source']:r for r in load(P/'outputs/direct_T_fine_proposal_index_v1.json')};totals=collections.Counter();newfine=collections.Counter();cross=[]
 for source,oldentry in oi.items():
  ne=ni[source]
  for entry,base in [(oldentry,ob),(ne,nb)]:assert sha(base/entry['metadata'])==entry['metadata_SHA256'] and sha(base/entry['arrays'])==entry['arrays_SHA256']
  om=load(ob/oldentry['metadata']);nm=load(nb/ne['metadata'])
  with np.load(ob/oldentry['arrays'],allow_pickle=False) as f:o={k:f[k] for k in f.files}
  with np.load(nb/ne['arrays'],allow_pickle=False) as f:n={k:f[k] for k in f.files}
  pr=pi[source];pp=P/'outputs'/pr['metadata'];assert sha(pp)==pr['metadata_SHA256'];pm=load(pp);pa=P/'outputs'/pm['arrays'];assert sha(pa)==pm['arrays_SHA256']
  with np.load(pa,allow_pickle=False) as f:proposal=f['proposed_fine_code'];reason=f['fine_review_code'];np.testing.assert_array_equal(f['native_row'],o['native_row'])
  changed=proposal>0;assert set(n)==set(o)|({'fine_revision_review_code'} if changed.any() else set())
  for k in o:
   if k not in ['consensus_celltype_code','annotation_group_code']:np.testing.assert_array_equal(o[k],n[k])
  before=np.asarray(om['candidate_dictionaries']['consensus_celltype'],object)[o['consensus_celltype_code']];after=np.asarray(nm['candidate_dictionaries']['consensus_celltype'],object)[n['consensus_celltype_code']];want=before.copy();want[changed]=np.asarray(pm['fine_dictionary'],object)[proposal[changed]];np.testing.assert_array_equal(after,want)
  assert (before[changed]=='T').all();major=np.asarray(om['candidate_dictionaries']['consensus_major'],object)[o['consensus_major_code']];assert (major[changed]=='T').all();newfine.update(after[changed])
  expected_dict={k:list(v) for k,v in om['candidate_dictionaries'].items()}
  if changed.any():
   for label in sorted(set(want[changed])):
    if label not in expected_dict['consensus_celltype']:expected_dict['consensus_celltype'].append(label)
   expected_dict['fine_revision_review']=['not_revised_by_this_direct_fine_node',pm['review_dictionary'][1],pm['review_dictionary'][2]];np.testing.assert_array_equal(n['fine_revision_review_code'],np.where(changed,reason,0).astype(np.uint8))
  assert nm['candidate_dictionaries']==expected_dict and nm['celltype_counts']==dict(collections.Counter(want)) and nm['major_counts']==om['major_counts'] and nm['n_direct_fine_updates']==int(changed.sum())
  allowed={'previous_label_version','label_version','arrays','arrays_SHA256','candidate_dictionaries','n_source_annotation_groups','n_direct_fine_updates','major_counts','celltype_counts','old_to_new_celltype_counts','group_catalog','identity_state_and_row_arrays_unchanged','RNA_state_native_rows_source_gate_and_major_unchanged'}
  for k,v in om.items():
   if k not in allowed:assert nm[k]==v,(source,k)
  assert nm['RNA_state_native_rows_source_gate_and_major_unchanged'] and nm['identity_state_and_row_arrays_unchanged']==(not changed.any())
  sg=sorted([g for g in groups if g['source']==source],key=lambda g:g['source_annotation_group_code']);assert [g['source_annotation_group_code'] for g in sg]==list(range(len(sg)))
  for g in sg:
   code=g['source_annotation_group_code'];sel=n['annotation_group_code']==code;assert int(sel.sum())==g['n_scope_group_rows'];parents=set(map(int,o['annotation_group_code'][sel]));assert len(parents)==1;parent=parents.pop();og=oldby[(source,parent)]
   assert g['native_rows_SHA256']==hashlib.sha256(np.asarray(n['native_row'][sel],'<i8').tobytes()).hexdigest();assert g['global_subcluster_number']==expected[g['stable_group_sha256']] and g['global_subcluster_id']=='TME_annotation_group_'+str(g['global_subcluster_number'])
   if (source,g['stable_group_sha256']) in subsets:
    s,h=subsets[(source,g['stable_group_sha256'])];assert s['previous_group_SHA256']==og['stable_group_sha256'];np.testing.assert_array_equal(sel,(o['annotation_group_code']==parent)&(want==s['fine_celltype']))
    assert g['consensus_celltype']==s['fine_celltype'] and g['major_distribution']=={'T':int(sel.sum())} and g['native_rows_SHA256']==s['native_rows_SHA256'] and g['previous_fine_revision_group_SHA256']==og['stable_group_sha256'] and g['previous_fine_revision_group_number']==og['global_subcluster_number']
    choices=[v for v in s['own_full_native_gene_marker_candidates_not_DEG'] if re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9.-]*',v['gene']) and v['gene'].lower() not in ['nan','none','unknown']];best=choices[0] if choices else None
    assert g['marker']==(best['gene'] if best else None) and g['marker_origin']==(g['stable_group_sha256'] if best else 'unavailable') and g['marker_expression_evidence_metadata_SHA256']==h
    for field,k in [('marker_native_mean_CP10K','native_mean_CP10K'),('marker_detection_fraction','native_detection_fraction'),('marker_descriptive_log2_mean_ratio','descriptive_log2_mean_ratio')]:assert g[field]==(best[k] if best else None)
    assert g['marker_evidence_type']==(best['marker_evidence_type'] if best else 'no_accepted_current_exact_subset_marker') and g['marker_comparator_scope']==s['comparator_scope'] and g['marker_comparator_cell_count']==s['n_comparator_cells']
    assert g['marker_evidence_membership_matches_exact_group'] and g['marker_DEG_status']=='no_qualified_replicate_aware_DEG_supplied'
    relation='same_membership_fine_label_revision' if s['exact_membership_reused'] else 'split'
    if s['exact_membership_reused']:assert g['global_subcluster_number']==og['global_subcluster_number']
    else:assert g['global_subcluster_number']>max(history.values())
    for k,v in og.items():
     if k=='marker' or k.startswith('marker_') or k in ['parent_marker_borrowed','borrow_reason','subcluster']:assert g['previous_fine_revision__'+k]==v
     elif k not in {'source_annotation_group_code','stable_group_sha256','n_scope_group_rows','consensus_celltype','major_distribution','native_rows_SHA256','membership_identical_to_archive_parent','global_subcluster_number','global_subcluster_id'}:assert g[k]==v,(source,k)
   else:
    assert g==og;np.testing.assert_array_equal(sel,o['annotation_group_code']==parent);relation='same_membership'
   assert not g['parent_marker_borrowed'] and g['subcluster']==g['consensus_celltype']+'_'+(g['marker'] or 'marker_unresolved')+'_'+str(g['global_subcluster_number'])
   cross.append(dict(source=source,old_group_SHA256=og['stable_group_sha256'],old_group_number=og['global_subcluster_number'],old_display=og['subcluster'],new_group_SHA256=g['stable_group_sha256'],n_records=int(sel.sum()),relation=relation,new_celltype=g['consensus_celltype'],new_group_number=g['global_subcluster_number'],new_display=g['subcluster']))
  assert nm['n_source_annotation_groups']==len(sg);totals.update(n_sources=1,n_source_qualified_rows=len(n['native_row']),n_fit_candidates_not_final_QC=int(n['fit_pool_candidate_not_final_QC'].sum()),n_without_fitted_group=int((~n['fit_pool_candidate_not_final_QC']).sum()),n_TCR_coarse_updates=nm.get('n_TCR_coarse_updates',0),n_direct_fine_updates=int(changed.sum()),n_exact_scoped_groups=len(sg))
 xp=m['annotation_membership_crosswalk'];assert sha(xp['path'])==xp['SHA256'];actual=load(xp['path']);order=lambda r:(r['source'],r['new_group_SHA256']);assert sorted(actual,key=order)==sorted(cross,key=order)
 assert sum(r['n_records'] for r in actual)==totals['n_fit_candidates_not_final_QC'];assert load(m['summary']['path'])['totals']==dict(totals)==load(U/'provenance/terminal_v1.json')['totals'];assert dict(newfine)==load(P/'provenance/terminal_v1.json')['new_fine_counts']
 save(root/'provenance/terminal_v1.json',dict(status='ALL_TME_DIRECT_FINE_WORKING_DELTAS_HISTORY_STATES_AND_OWN_MARKERS_ACCEPTED',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),new_fine_counts=dict(newfine),producer_manifest_SHA256=sha(mp),producer_terminal_SHA256=sha(U/'provenance/terminal_v1.json'),all_scope_deltas_and_crosswalks_checked=True,all_inactive_historical_numbers_retained=True,not_orthogonal_identity_accuracy=True,final_module_input_ready=False,wall_seconds=time.monotonic()-start))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
