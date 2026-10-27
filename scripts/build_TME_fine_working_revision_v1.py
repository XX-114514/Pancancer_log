"""Apply accepted direct T fine cues with exact own profiles and full history."""
import argparse,collections,datetime,hashlib,json,re,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):
 with Path(p).open('x') as f:json.dump(d,f,indent=2,allow_nan=False);f.write('\n')
def hist(a):return dict(collections.Counter(map(str,a)))

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');S=Path(cfg['native_subsets']);SA=Path(cfg['native_subset_audit']);sc=load(S/'config/spec_v1.json');P=Path(sc['proposal']);PA=Path(sc['proposal_audit']);pc=load(P/'config/spec_v1.json');J=Path(pc['joined'])
 for stage in [S,SA,P,PA]:assert load(stage/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 for producer,auditor in [(S,SA),(P,PA)]:assert load(auditor/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(producer/'provenance/terminal_v1.json')
 accepted=load(pc['working_acceptance']);oldmp=Path(accepted['producer_manifest']['path']);assert sha(oldmp)==accepted['producer_manifest']['SHA256'];old=load(oldmp)
 for name in ['index','group_catalog','summary','display_number_registry']:assert sha(old[name]['path'])==old[name]['SHA256']
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(root/'scripts'))
 from atlas_identity_release_keys_v1 import assign_display_numbers,display_label
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 oldgroups=load(old['group_catalog']['path'])['groups'];history=load(old['display_number_registry']['path'])
 assert all(history[g['stable_group_sha256']]==g['global_subcluster_number'] for g in oldgroups)
 sip=S/'outputs/fine_delta_native_profile_index_v1.json';assert sha(sip)==load(S/'provenance/terminal_v1.json')['index_SHA256'];subsets={};subset_hash={}
 for r in load(sip):
  mp=S/'outputs'/r['metadata'];assert sha(mp)==r['metadata_SHA256'];subsets[r['source']]=load(mp);subset_hash[r['source']]=sha(mp)
 pip=P/'outputs/direct_T_fine_proposal_index_v1.json';assert sha(pip)==load(P/'provenance/terminal_v1.json')['index_SHA256'];pi={r['source']:r for r in load(pip)}
 base=Path(old['index']['path']).parent;catalog=[];metas=[];cross=[];totals=collections.Counter();sources=[];version='TME_source_qualified_working_annotation_20261027_v4_direct_T_fine'
 for item in load(old['index']['path']):
  storage_gate(root,pc);source=item['source'];mp=base/item['metadata'];arr=base/item['arrays'];assert sha(mp)==item['metadata_SHA256'] and sha(arr)==item['arrays_SHA256'];m=load(mp)
  with np.load(arr,allow_pickle=False) as f:q={k:f[k] for k in f.files}
  before=np.asarray(m['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']];oldcodes=q['annotation_group_code'].copy();oldgs=sorted([dict(g) for g in oldgroups if g['source']==source],key=lambda g:g['source_annotation_group_code']);groups=[dict(g) for g in oldgs];changed=np.zeros(len(before),bool);n_changed=0;relations={g['stable_group_sha256']:'same_membership' for g in groups}
  if source in subsets:
   pr=pi[source];pm=P/'outputs'/pr['metadata'];assert sha(pm)==pr['metadata_SHA256'];p=load(pm);pa=P/'outputs'/p['arrays'];assert sha(pa)==p['arrays_SHA256']
   with np.load(pa,allow_pickle=False) as f:proposed=f['proposed_fine_code'];reason=f['fine_review_code'];np.testing.assert_array_equal(f['native_row'],q['native_row'])
   changed=proposed>0;n_changed=int(changed.sum());assert (before[changed]=='T').all();fine=np.asarray(p['fine_dictionary'],object)[proposed]
   for label in sorted(set(fine[changed])):
    if label not in m['candidate_dictionaries']['consensus_celltype']:m['candidate_dictionaries']['consensus_celltype'].append(label)
    q['consensus_celltype_code'][changed&(fine==label)]=m['candidate_dictionaries']['consensus_celltype'].index(label)
   assert 'fine_revision_review_code' not in q;q['fine_revision_review_code']=np.where(changed,reason,0).astype(np.uint8);m['candidate_dictionaries']['fine_revision_review']=['not_revised_by_this_direct_fine_node',p['review_dictionary'][1],p['review_dictionary'][2]]
   after=np.asarray(m['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']];byparent=collections.defaultdict(list)
   for s in subsets[source]['groups']:byparent[s['previous_group_code']].append(s)
   for parent,children in sorted(byparent.items()):
    oldg=oldgs[parent];assert all(s['previous_group_SHA256']==oldg['stable_group_sha256'] for s in children)
    for j,s in enumerate(sorted(children,key=lambda s:s['fine_celltype'])):
     mask=(oldcodes==parent)&(after==s['fine_celltype']);assert int(mask.sum())==s['n_exact_subset'];g=dict(oldg)
     for name in list(g):
      if name=='marker' or name.startswith('marker_') or name in ['parent_marker_borrowed','borrow_reason','subcluster']:g['previous_fine_revision__'+name]=g.pop(name)
     choices=[v for v in s['own_full_native_gene_marker_candidates_not_DEG'] if re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9.-]*',v['gene']) and v['gene'].lower() not in ['nan','none','unknown']];best=choices[0] if choices else None
     code=parent if j==0 else len(groups)
     g.update(source_annotation_group_code=code,stable_group_sha256=s['stable_membership_SHA256'],n_scope_group_rows=s['n_exact_subset'],consensus_celltype=s['fine_celltype'],major_distribution={'T':s['n_exact_subset']},native_rows_SHA256=s['native_rows_SHA256'],membership_identical_to_archive_parent=s['stable_membership_SHA256']==oldg['parent_archive_group_sha256'],previous_fine_revision_group_SHA256=oldg['stable_group_sha256'],previous_fine_revision_group_number=oldg['global_subcluster_number'],marker=best['gene'] if best else None,marker_origin=s['stable_membership_SHA256'] if best else 'unavailable',marker_evidence_type=best['marker_evidence_type'] if best else 'no_accepted_current_exact_subset_marker',marker_native_mean_CP10K=best['native_mean_CP10K'] if best else None,marker_detection_fraction=best['native_detection_fraction'] if best else None,marker_descriptive_log2_mean_ratio=best['descriptive_log2_mean_ratio'] if best else None,marker_comparator_scope=s['comparator_scope'],marker_comparator_cell_count=s['n_comparator_cells'],marker_expression_evidence_metadata_SHA256=subset_hash[source],marker_evidence_membership_matches_exact_group=True,marker_selection_rule='accepted_exact_own_native_gene_descriptive_ranking_source_T_NK_comparator_not_DEG',marker_DEG_status='no_qualified_replicate_aware_DEG_supplied',marker_expression_unit_support_status='not_newly_estimated_no_patient_independence_claim',marker_small_group_status='small_group_descriptive_only' if s['n_exact_subset']<30 else 'descriptive_expression_not_inferential_DEG',parent_marker_borrowed=False,borrow_reason=None)
     q['annotation_group_code'][mask]=code
     if j==0:groups[parent]=g
     else:groups.append(g)
     relations[g['stable_group_sha256']]='same_membership_fine_label_revision' if s['exact_membership_reused'] else 'split'
   out=root/'outputs'/(source+'_TME_working_annotation_v4.npz')
   with out.open('xb') as f:np.savez_compressed(f,**q)
   with np.load(out,allow_pickle=False) as f:
    for k,v in q.items():np.testing.assert_array_equal(f[k],v)
   arrays=out
  else:arrays=arr.resolve()
  major=np.asarray(m['candidate_dictionaries']['consensus_major'],object)[q['consensus_major_code']];fine=np.asarray(m['candidate_dictionaries']['consensus_celltype'],object)[q['consensus_celltype_code']]
  for g in groups:
   sel=q['annotation_group_code']==g['source_annotation_group_code'];parents=set(map(int,oldcodes[sel]));assert len(parents)==1;og=oldgs[parents.pop()]
   cross.append(dict(source=source,old_group_SHA256=og['stable_group_sha256'],old_group_number=og['global_subcluster_number'],old_display=og['subcluster'],new_group_SHA256=g['stable_group_sha256'],n_records=int(sel.sum()),relation=relations[g['stable_group_sha256']],new_celltype=g['consensus_celltype']))
  m.update(previous_label_version=m['label_version'],label_version=version,arrays=str(arrays),arrays_SHA256=sha(arrays),n_source_annotation_groups=len(groups),n_direct_fine_updates=n_changed,major_counts=hist(major),celltype_counts=hist(fine),identity_state_and_row_arrays_unchanged=n_changed==0,RNA_state_native_rows_source_gate_and_major_unchanged=True,old_to_new_celltype_counts=[dict(previous_working_celltype=a,current_celltype=b,n_records=n) for (a,b),n in sorted(collections.Counter(zip(before,fine)).items())])
  metas.append(m);catalog+=groups;totals.update(n_sources=1,n_source_qualified_rows=len(major),n_fit_candidates_not_final_QC=int(q['fit_pool_candidate_not_final_QC'].sum()),n_without_fitted_group=int((~q['fit_pool_candidate_not_final_QC']).sum()),n_TCR_coarse_updates=m.get('n_TCR_coarse_updates',0),n_direct_fine_updates=n_changed,n_exact_scoped_groups=len(groups));sources.append(dict(source=source,major_counts=hist(major),celltype_counts=hist(fine)))
 assigned=assign_display_numbers([g['stable_group_sha256'] for g in catalog],history)
 for g in catalog:g.update(global_subcluster_number=assigned['current'][g['stable_group_sha256']],global_subcluster_id='TME_annotation_group_'+str(assigned['current'][g['stable_group_sha256']]),subcluster=display_label(g['consensus_celltype'],g['marker'],assigned['current'][g['stable_group_sha256']]))
 by={g['stable_group_sha256']:g for g in catalog}
 for r in cross:r.update(new_group_number=by[r['new_group_SHA256']]['global_subcluster_number'],new_display=by[r['new_group_SHA256']]['subcluster'])
 gp=root/'outputs/TME_working_group_catalog_v4.json';xp=root/'outputs/TME_annotation_membership_crosswalk_v4.json';rp=root/'outputs/TME_display_number_registry_v4.json';save(gp,{'groups':catalog});save(xp,cross);save(rp,assigned['registry']);index=[]
 for m in metas:
  m['group_catalog']={'path':str(gp),'SHA256':sha(gp)};op=root/'outputs'/(m['source']+'_TME_working_annotation_v4.json');save(op,m);index.append(dict(source=m['source'],metadata=op.name,metadata_SHA256=sha(op),arrays=m['arrays'],arrays_SHA256=m['arrays_SHA256'],n_source_qualified_rows=m['n_source_qualified_rows']))
 ip=root/'outputs/TME_working_annotation_index_v4.json';sp=root/'outputs/TME_working_annotation_summary_v4.json';save(ip,index);save(sp,dict(totals=dict(totals),sources=sources,final_module_input_ready=False))
 bindings=dict(old['bindings']);bindings.update({k:dict(path=str(p),SHA256=sha(p)) for k,p in dict(previous_working_annotation=oldmp,previous_display_number_registry=Path(old['display_number_registry']['path']),direct_fine_proposal_index=pip,direct_fine_proposal_audit=PA/'provenance/terminal_v1.json',fine_delta_native_index=sip,fine_delta_native_audit=SA/'provenance/terminal_v1.json').items()})
 manifest=dict(old,label_version=version,bindings=bindings,index=dict(path=str(ip),SHA256=sha(ip)),group_catalog=dict(path=str(gp),SHA256=sha(gp)),summary=dict(path=str(sp),SHA256=sha(sp)),annotation_membership_crosswalk=dict(path=str(xp),SHA256=sha(xp)),display_number_registry=dict(path=str(rp),SHA256=sha(rp)),release_status='WORKING_DIRECT_FINE_REVISION_INDEPENDENT_DELTA_AUDIT_REQUIRED_NOT_FINAL_TAXONOMY');manifest.pop('marker_display_crosswalk',None)
 mp=root/'outputs/TME_direct_fine_working_annotation_manifest_v4.json';save(mp,manifest);assert totals['n_direct_fine_updates']==load(P/'provenance/terminal_v1.json')['totals']['n_fine_changes']
 save(root/'provenance/terminal_v1.json',dict(status='FULL_TME_DIRECT_FINE_WORKING_REVISION_REQUIRES_INDEPENDENT_AUDIT',completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),totals=dict(totals),manifest_SHA256=sha(mp),new_CM_started=False,final_module_input_ready=False,wall_seconds=time.monotonic()-start))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);r=Path(p.parse_args().output)
 try:main(r)
 except BaseException:
  import traceback
  save(r/'provenance/failure_v1.json',{'traceback':traceback.format_exc()});raise
