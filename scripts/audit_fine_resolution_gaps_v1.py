"""Independent complete-array replay of the native fine-resolution diagnostic.

The diagnostic producer and its pattern-evaluation helper are not imported.
Accepted native/reference sidecars remain evidence, not biological truth.
"""
import argparse,ast,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def hist(x):
 a,b=np.unique(np.asarray(x).astype(str),return_counts=True);return dict(zip(a.tolist(),map(int,b)))

def main(root):
 started=time.monotonic();R=Path(load(root/'config/spec_v1.json')['producer']);cfg=load(R/'config/spec_v1.json')
 terminal=load(R/'provenance/terminal_v1.json');supervisor=load(R/'provenance/supervisor_terminal_v1.json');assert supervisor['actual_child_exit_code']==0 and not supervisor['wall_timeout']
 for b in cfg['bindings'].values():assert sha(b['path'])==b['sha256']
 ip=R/'outputs/fine_resolution_gap_index_v1.json';sp=R/'outputs/fine_resolution_gap_summary_v1.json';assert sha(ip)==terminal['index_SHA256'] and sha(sp)==terminal['summary_SHA256']
 summary=load(sp);index=load(ip);J=Path(cfg['joined']);sys.path.insert(0,str(J/'scripts'));sys.path.insert(0,str(Path(cfg['frozen_candidates'])/'scripts'))
 from joined_native_atlas_delivery_v1 import JoinedNativeAtlas
 from cell_identity_confidence_v2 import FINE_TYPES
 from native_lineage_routing_v1 import FAMILIES
 atlas=JoinedNativeAtlas(J/'outputs/joined_native_atlas_delivery_v1.json');assert [i['source'] for i in index]==list(atlas.sources)
 tree=ast.parse((R/'scripts/native_fine_resolution_v1.py').read_text())
 patterns=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='FINE_ANCHOR_PATTERNS' for t in n.targets))
 order=list(patterns);assert len(order)==11
 view=load(cfg['source_view']);core=load(view['work_units']['path'])['units'];bykey={u['qualification_unit_key_sha256']:u for u in core}
 inherited=load(view['base_v4_manifest']['path']);origin=inherited['bindings']['author_bridge_registry'];assert sha(origin['path'])==origin['sha256'];origin_ds=load(origin['path'])['source'];overrides={origin_ds:{'path':inherited['override_arrays_path'],'sha256':inherited['override_arrays_sha256']},**view['new_source_overrides']}
 models=load(atlas.evidence/'outputs/reference_models_v1.json');labels=models['immune_low']['labels'];total=collections.Counter();scope_limits=[]
 for entry,ss in zip(index,summary['source_summaries']):
  ds=entry['source'];assert sha(entry['path'])==entry['SHA256'];report=load(entry['path']);ap=Path(entry['path']).parent/report['arrays'];assert sha(ap)==entry['arrays_SHA256']==report['arrays_SHA256']
  with np.load(ap,allow_pickle=False) as f:d={k:f[k] for k in f.files}
  a=atlas.source(ds);ca,ea,em=a['ca'],a['ea'],a['em'];n=len(a['rows']);np.testing.assert_array_equal(d['native_row'],a['rows']);assert report['native_ids_SHA256']==em['native_ids_sha256']
  assert report['pattern_order']==order and report['gene_axis']==em['review_genes'];genes=em['review_genes'];measurement=np.asarray(em['gene_measured'],bool);lookup={g:i for i,g in enumerate(genes)}
  detection=np.unpackbits(ea['native_gene_detection_packbits'],axis=1,count=len(genes),bitorder='little').astype(bool)
  if ds in overrides:
   o=overrides[ds];assert sha(o['path'])==o['sha256']
   with np.load(o['path'],allow_pickle=False) as f:np.testing.assert_array_equal(f['native_row'],a['rows']);status=f['source_qualification_code'].copy()
  else:
   decode={'排除':0,'来源资格通过':1,'待核实':2};u=[bykey[r['native_unit_key_sha256']] for r in a['um']['unit_records']];codes=np.asarray([decode[r['source_qualification']] for r in u],np.int8);status=codes[a['ua']['native_expression_unit_code']]
  np.testing.assert_array_equal(d['source_qualification_code'],status)
  major=np.asarray(a['cm']['candidate_dictionaries']['consensus_major'],object)[ca['consensus_major_code']];fine=np.asarray(a['cm']['candidate_dictionaries']['consensus_celltype'],object)[ca['consensus_celltype_code']];reason=np.asarray(a['cm']['candidate_dictionaries']['coarse_evidence'],object)[ca['coarse_evidence_code']]
  fit=ca['fit_pool_candidate'];active=status==1;scopes={'all_archive':np.ones(n,bool),'fit_pool_not_final_QC':fit,'source_qualified_archive':active,'source_qualified_fit_not_final_QC':active&fit}
  positive_bits=np.zeros(n,np.uint16);unknown_bits=np.zeros(n,np.uint16);all_group=[]
  for j,(t,clauses) in enumerate(patterns.items()):
   state=np.ones(n,np.uint8);descriptions=[]
   for clause in clauses:
    available=[g for g in clause if g in lookup and measurement[lookup[g]]];clause_state=np.full(n,0 if len(available)==len(clause) else 2,np.uint8)
    if available:clause_state[detection[:,[lookup[g] for g in available]].any(axis=1)]=1
    state=np.where((state==0)|(clause_state==0),0,np.where((state==2)|(clause_state==2),2,1)).astype(np.uint8)
    descriptions.append({'requested':list(clause),'measured':available,'fully_measured':len(available)==len(clause)})
   positive=state==1;unknown=state==2;positive_bits|=positive.astype(np.uint16)<<j;unknown_bits|=unknown.astype(np.uint16)<<j
   assert report['pattern_measurement_coverage'][t]=={'AND_of_OR_clauses':descriptions,'all_requested_genes_measured':all(x['fully_measured'] for x in descriptions)}
   record=report['type_reports'][j];assert record==ss['type_reports'][j] and record['type_pattern']==t
   family='T_NK' if j<4 else 'Endothelial' if t=='Lymphatic_endothelial_related' else 'Myeloid';assert record['family_context']==family
   native=ca['native_family_bitmask']==(1<<FAMILIES.index(family));compatible=np.isin(major,{'T_NK':['T','T_NK'],'Myeloid':['Myeloid'],'Endothelial':['Endothelial']}[family])
   k=FINE_TYPES.index(t) if t in FINE_TYPES else None;own=(ca['own_fine_bitmask']&(1<<k))!=0 if k is not None else np.zeros(n,bool);ref=(ca['reference_fine_bitmask']&(1<<k))!=0 if k is not None else np.zeros(n,bool)
   vocab=summary['reference_fine_vocabulary'][t];assert record['reference_low_vocab_labels']==vocab
   ref_codes=[labels.index(label) for label in vocab];best=np.isin(ea['immune_low__best_code'],ref_codes);unique=np.isin(ea['immune_low__unique_code'],ref_codes)
   # The reported fine vocabulary must reproduce the already accepted cell-level
   # reference bit for every observed unique prediction, not just a few rows.
   if k is not None:np.testing.assert_array_equal(ref,unique)
   org_best=np.zeros(n,bool);org_unique=np.zeros(n,bool)
   for m,vocab in record['scoped_organ_reference_vocab_labels'].items():
    assert m in em['source_reference_models'];codes=[models[m]['labels'].index(v) for v in vocab]
    org_best|=np.isin(ea[m+'__best_code'],codes);org_unique|=np.isin(ea[m+'__unique_code'],codes)
    assert not (org_best&~ea['organ_reference_scope']).any() and not (org_unique&~ea['organ_reference_scope']).any()
   for name,mask in scopes.items():
    p=positive&mask;c=fine==t;valid_ref=p&(ea['immune_low__best_code']>=0)
    expected={'n_scope_rows':int(mask.sum()),'n_positive_anchor_pattern':int(p.sum()),'n_unassessable_required_pattern':int((unknown&mask).sum()),'n_old_native_fine_panel_support':int((own&mask).sum()),'n_old_unique_reference_fine_support':int((ref&mask).sum()),'n_current_fine_label':int((c&mask).sum()),'n_anchor_positive_current_fine_matches':int((p&c).sum()),'n_anchor_positive_current_coarse_compatible':int((p&compatible).sum()),'n_anchor_positive_single_native_family_matches':int((p&native).sum()),'n_anchor_positive_current_major_unresolved':int((p&(major=='unresolved')).sum()),'n_anchor_positive_best_reference_fine_matches':int((p&best).sum()),'n_anchor_positive_unique_reference_fine_matches':int((p&unique).sum()),'n_anchor_positive_best_fine_match_unique_rejected':int((p&best&~unique).sum()),'n_anchor_positive_scoped_organ_best_fine_matches':int((p&org_best).sum()),'n_anchor_positive_scoped_organ_unique_fine_matches':int((p&org_unique).sum()),'n_anchor_positive_no_eligible_unique_reference_fine':int((p&(ca['reference_fine_bitmask']==0)).sum()),'anchor_positive_coarse_candidate_distribution':hist(major[p]),'anchor_positive_coarse_evidence_distribution':hist(reason[p]),'anchor_positive_low_reference_n_positive_labels':hist(ea['immune_low__n_positive_labels'][p]),'anchor_positive_low_reference_best_label_distribution':hist(np.asarray(labels,object)[ea['immune_low__best_code'][valid_ref]]),'anchor_positive_native_total_RNA_quantiles':np.quantile(ea['native_total_RNA'][p],[0,.25,.5,.75,1]).tolist() if p.any() else None}
    assert expected==record['scopes'][name],(ds,t,name,'diagnostic scalar/distribution differs')
    if t in ['Neutrophil','Lymphatic_endothelial_related']:
     scope_limits.append({'source':ds,'pattern':t,'scope':name,'n_anchor_positive':int(p.sum()),'n_anchor_positive_inside_existing_organ_reference_scope':int((p&ea['organ_reference_scope']).sum()),'n_anchor_positive_inside_existing_immune_reference_scope':int((p&ea['immune_reference_scope']).sum()),'zero_organ_reference_matches_outside_scope_are_unavailable_not_negative':True})
   groupcode=ca['guarded_source_subset_code'];ng=len([g for g in atlas.groups.values() if g['source']==ds]);group_indices,group_n=np.unique(groupcode[fit],return_counts=True);lookup_n=dict(zip(group_indices.tolist(),map(int,group_n)))
   pos_indices,pos_n=np.unique(groupcode[fit&positive],return_counts=True)
   qi,qn=np.unique(groupcode[fit&positive&active],return_counts=True);qualified_counts=dict(zip(qi.tolist(),map(int,qn)))
   for g,count in zip(pos_indices,pos_n):all_group.append({'guarded_subset_code':int(g),'pattern':t,'n_original_fit_members':lookup_n[int(g)],'n_anchor_positive_fit_members':int(count),'n_source_qualified_anchor_positive_fit_members':qualified_counts.get(int(g),0)})
  np.testing.assert_array_equal(d['positive_anchor_bitmask'],positive_bits);np.testing.assert_array_equal(d['unassessable_anchor_bitmask'],unknown_bits);assert report['exact_group_pattern_support']==all_group
  assert report['n_archive_rows']==ss['n_archive_rows']==n;assert report['n_source_qualified_rows']==int(active.sum()) and report['n_source_qualified_fit_not_final_QC']==int((active&fit).sum())
  total.update(n_sources=1,n_archive_rows=n,n_source_qualified_rows=int(active.sum()),n_source_qualified_fit_not_final_QC=int((active&fit).sum()))
  save(root/'checkpoints'/('source_'+ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source':ds,'n_sources_checked':total['n_sources'],'n_archive_rows_checked':total['n_archive_rows']})
  del a,ca,ea,em,detection,d,major,fine,reason,positive_bits,unknown_bits;atlas._cache.clear();atlas.native._fields.clear();gc.collect()
 assert dict(total)==summary['totals']==terminal['totals']
 out={'status':'ALL_ARCHIVE_FINE_ANCHOR_MASKS_REFERENCE_RETENTION_SCALARS_AND_GROUP_COUNTS_INDEPENDENTLY_ACCEPTED','totals':dict(total),'producer_rule_or_summary_helpers_imported':False,'native_count_matrix_read':False,'not_biological_type_validation_or_label_release':True,'existing_organ_reference_scope_diagnostics':scope_limits,'producer_terminal_SHA256':sha(R/'provenance/terminal_v1.json'),'summary_SHA256':terminal['summary_SHA256']};op=root/'outputs/independent_fine_resolution_gap_audit_v1.json';save(op,out)
 save(root/'provenance/terminal_v1.json',{'status':out['status'],'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'output_SHA256':sha(op),'producer_terminal_SHA256':out['producer_terminal_SHA256'],'totals':dict(total),'wall_seconds':time.monotonic()-started,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException as e:
  import traceback
  save(root/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e),'traceback':traceback.format_exc()});raise
