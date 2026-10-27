"""Audit where fine identities are lost between native evidence and references.

Reads accepted detection sidecars for the complete archive. This stage neither
releases labels nor selects a rule by clinical associations or reduced conflict.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,subprocess,sys,time
from pathlib import Path
import numpy as np

def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def hist(x):
 names,n=np.unique(np.asarray(x).astype(str),return_counts=True)
 return dict(zip(names.tolist(),map(int,n)))

def main(root):
 began=time.monotonic();cfg=load(root/'config/spec_v1.json')
 for b in cfg['bindings'].values():assert sha(b['path'])==b['sha256']
 J=Path(cfg['joined']);frozen=Path(cfg['frozen_candidates'])
 sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(frozen/'scripts'))
 from joined_native_atlas_delivery_v1 import JoinedNativeAtlas
 from scoped_low_reference_identity_v1 import scoped_low_identity
 from cell_identity_confidence_v2 import FINE_TYPES
 from native_lineage_routing_v1 import FAMILIES
 from native_fine_resolution_v1 import evaluate_patterns,FINE_ANCHOR_PATTERNS,ANCHOR_FAMILIES
 from review_source_qualified_identity_context_v1 import source_status
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 test=subprocess.run([sys.executable,'-B',str(root/'tests/test_native_fine_resolution_v1.py')],capture_output=True,text=True,timeout=60)
 (root/'provenance/test_stderr_v1.txt').write_text(test.stderr)
 assert test.returncode==0 and 'Ran 7 tests' in test.stderr and 'skipped' not in test.stderr.lower()
 atlas=JoinedNativeAtlas(J/'outputs/joined_native_atlas_delivery_v1.json')
 view=load(cfg['source_view']);core=load(view['work_units']['path'])['units'];assert sha(view['work_units']['path'])==view['work_units']['sha256']
 models=load(atlas.evidence/'outputs/reference_models_v1.json');labels=models['immune_low']['labels'];fine_by_label=[scoped_low_identity(s)[1] for s in labels]
 vocabulary={t:[s for s,v in zip(labels,fine_by_label) if t in v] for t in FINE_ANCHOR_PATTERNS}
 summary=[];index=[];total=collections.Counter()
 for ds in atlas.sources:
  storage_gate(root,cfg);s=atlas.source(ds);ca,ea,em=s['ca'],s['ea'],s['em'];n=len(s['rows'])
  result=evaluate_patterns(ea['native_gene_detection_packbits'],em['review_genes'],np.asarray(em['gene_measured'],bool))
  status=source_status(s,core,view);fit=ca['fit_pool_candidate'].astype(bool);active=status==1
  scopes={'all_archive':np.ones(n,bool),'fit_pool_not_final_QC':fit,'source_qualified_archive':active,'source_qualified_fit_not_final_QC':active&fit}
  major=np.asarray(s['cm']['candidate_dictionaries']['consensus_major'],object)[ca['consensus_major_code']]
  fine=np.asarray(s['cm']['candidate_dictionaries']['consensus_celltype'],object)[ca['consensus_celltype_code']]
  reason=np.asarray(s['cm']['candidate_dictionaries']['coarse_evidence'],object)[ca['coarse_evidence_code']]
  ref_lookup=np.asarray([sum(1<<FINE_TYPES.index(t) for t in ts) for ts in fine_by_label],np.uint16)
  ref_expected=np.zeros(n,np.uint16);have_unique=ea['immune_low__unique_code']>=0
  ref_expected[have_unique]=ref_lookup[ea['immune_low__unique_code'][have_unique]]
  np.testing.assert_array_equal(ref_expected,ca['reference_fine_bitmask'])
  type_reports=[];group_reports=[]
  for j,t in enumerate(result['pattern_order']):
   positive=(result['positive_anchor_bitmask']&(1<<j))!=0;unknown=(result['unassessable_anchor_bitmask']&(1<<j))!=0
   family=ANCHOR_FAMILIES[t];family_bit=1<<FAMILIES.index(family);single_native=ca['native_family_bitmask']==family_bit
   compatible=np.isin(major,{'T_NK':['T','T_NK'],'Myeloid':['Myeloid'],'Endothelial':['Endothelial']}[family])
   k=FINE_TYPES.index(t) if t in FINE_TYPES else None
   old_own=(ca['own_fine_bitmask']&(1<<k))!=0 if k is not None else np.zeros(n,bool)
   old_ref=(ca['reference_fine_bitmask']&(1<<k))!=0 if k is not None else np.zeros(n,bool)
   target_codes=[i for i,v in enumerate(fine_by_label) if t in v]
   best=np.isin(ea['immune_low__best_code'],target_codes);unique=np.isin(ea['immune_low__unique_code'],target_codes)
   organ_best=np.zeros(n,bool);organ_unique=np.zeros(n,bool);organ_vocabulary={}
   for model in em['source_reference_models']:
    if model not in ['breast','lung']:continue
    def compatible_organ_label(label):
     z=label.lower()
     if t=='Neutrophil':return model=='breast' and z=='neutrophil'
     if t=='Lymphatic_endothelial_related':return (model=='breast' and z.startswith('lymph-')) or (model=='lung' and z.startswith('lymphatic ec'))
     return False
    match=[i for i,v in enumerate(models[model]['labels']) if compatible_organ_label(v)]
    organ_vocabulary[model]=[models[model]['labels'][i] for i in match]
    organ_best|=np.isin(ea[model+'__best_code'],match);organ_unique|=np.isin(ea[model+'__unique_code'],match)
   scope_reports={}
   for scope,mask in scopes.items():
    p=positive&mask;c=fine==t
    scope_reports[scope]={'n_scope_rows':int(mask.sum()),'n_positive_anchor_pattern':int(p.sum()),'n_unassessable_required_pattern':int((unknown&mask).sum()),
     'n_old_native_fine_panel_support':int((old_own&mask).sum()),'n_old_unique_reference_fine_support':int((old_ref&mask).sum()),'n_current_fine_label':int((c&mask).sum()),
     'n_anchor_positive_current_fine_matches':int((p&c).sum()),'n_anchor_positive_current_coarse_compatible':int((p&compatible).sum()),'n_anchor_positive_single_native_family_matches':int((p&single_native).sum()),
     'n_anchor_positive_current_major_unresolved':int((p&(major=='unresolved')).sum()),'n_anchor_positive_best_reference_fine_matches':int((p&best).sum()),
     'n_anchor_positive_unique_reference_fine_matches':int((p&unique).sum()),'n_anchor_positive_best_fine_match_unique_rejected':int((p&best&~unique).sum()),
     'n_anchor_positive_scoped_organ_best_fine_matches':int((p&organ_best).sum()),'n_anchor_positive_scoped_organ_unique_fine_matches':int((p&organ_unique).sum()),
     'n_anchor_positive_no_eligible_unique_reference_fine':int((p&(ca['reference_fine_bitmask']==0)).sum()),
     'anchor_positive_coarse_candidate_distribution':hist(major[p]),'anchor_positive_coarse_evidence_distribution':hist(reason[p]),
     'anchor_positive_low_reference_n_positive_labels':hist(ea['immune_low__n_positive_labels'][p]),
     'anchor_positive_low_reference_best_label_distribution':hist(np.asarray(labels,object)[ea['immune_low__best_code'][p& (ea['immune_low__best_code']>=0)]]),
     'anchor_positive_native_total_RNA_quantiles':np.quantile(ea['native_total_RNA'][p],[0,.25,.5,.75,1]).tolist() if p.any() else None}
   type_reports.append({'type_pattern':t,'family_context':family,'reference_low_vocab_labels':vocabulary[t],'scoped_organ_reference_vocab_labels':organ_vocabulary,'all_requested_genes_measured':result['coverage'][t]['all_requested_genes_measured'],'scopes':scope_reports})
   if positive[fit].any():
    codes=ca['guarded_source_subset_code'];ng=len([v for v in atlas.groups.values() if v['source']==ds]);z=np.bincount(codes[fit],minlength=ng);p=np.bincount(codes[fit&positive],minlength=ng);q=np.bincount(codes[fit&positive&active],minlength=ng)
    for group in np.flatnonzero(p):group_reports.append({'guarded_subset_code':int(group),'pattern':t,'n_original_fit_members':int(z[group]),'n_anchor_positive_fit_members':int(p[group]),'n_source_qualified_anchor_positive_fit_members':int(q[group])})
  npz=root/'outputs'/(ds+'_fine_resolution_anchor_evidence_v1.npz')
  arrays={'native_row':s['rows'],'positive_anchor_bitmask':result['positive_anchor_bitmask'],'unassessable_anchor_bitmask':result['unassessable_anchor_bitmask'],'source_qualification_code':status}
  with npz.open('xb') as f:np.savez_compressed(f,**arrays)
  with np.load(npz,allow_pickle=False) as f:
   assert set(f.files)==set(arrays)
   for k,v in arrays.items():np.testing.assert_array_equal(f[k],v)
  report={'source':ds,'n_archive_rows':n,'native_ids_SHA256':em['native_ids_sha256'],'n_source_qualified_rows':int(active.sum()),'n_source_qualified_fit_not_final_QC':int((active&fit).sum()),'n_original_fit_rows':int(fit.sum()),'gene_axis':em['review_genes'],'pattern_order':result['pattern_order'],'pattern_measurement_coverage':result['coverage'],'type_reports':type_reports,'exact_group_pattern_support':group_reports,'arrays':npz.name,'arrays_SHA256':sha(npz),'candidate_metadata_SHA256':atlas.ci[ds]['metadata_sha256'],'evidence_metadata_SHA256':atlas.ei[ds]['metadata_sha256'],'new_identity_QC_or_malignancy_assigned':False,'binary_detection_patterns_do_not_establish_RNA_amount_or_validate_identity':True,'native_patterns_and_reference_features_not_independent_validation':True}
  mp=npz.with_suffix('.json');save(mp,report);index.append({'source':ds,'n_archive_rows':n,'path':str(mp),'SHA256':sha(mp),'arrays_SHA256':sha(npz)})
  summary.append({k:report[k] for k in ['source','n_archive_rows','n_source_qualified_rows','n_source_qualified_fit_not_final_QC','n_original_fit_rows','type_reports']})
  total.update(n_sources=1,n_archive_rows=n,n_source_qualified_rows=int(active.sum()),n_source_qualified_fit_not_final_QC=int((active&fit).sum()))
  save(root/'checkpoints'/('source_'+ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_complete_sources':len(index),'source':ds,'n_complete_archive_rows':total['n_archive_rows']})
  del s,ca,ea,major,fine,reason,result,arrays,report,type_reports,group_reports;atlas._cache.clear();atlas.native._fields.clear();gc.collect()
 assert total['n_archive_rows']==sum(v['n_archive_rows'] for v in atlas.native.archive_partitions.values())
 assert total['n_source_qualified_rows']==sum(u['n_archive_rows'] for u in core if u['source_qualification']=='来源资格通过')
 ip=root/'outputs/fine_resolution_gap_index_v1.json';sp=root/'outputs/fine_resolution_gap_summary_v1.json';save(ip,index);save(sp,{'totals':dict(total),'reference_fine_vocabulary':vocabulary,'source_summaries':summary,'no_labels_or_QC_changed':True,'predefined_patterns_not_selected_by_yield_or_clinical_association':True})
 save(root/'provenance/terminal_v1.json',{'status':'ALL_ARCHIVE_NATIVE_FINE_PATTERNS_REFERENCE_GAPS_AND_CURRENT_RETENTION_MEASURED_INDEPENDENT_AUDIT_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(total),'summary_SHA256':sha(sp),'index_SHA256':sha(ip),'source_view_SHA256':sha(cfg['source_view']),'output_bytes':sum(p.stat().st_size for p in root.rglob('*') if p.is_file()),'native_count_matrix_read':False,'no_identity_or_QC_release':True,'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException as e:
  import traceback
  save(root/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e),'traceback':traceback.format_exc()});raise
