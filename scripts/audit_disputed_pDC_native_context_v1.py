"""Independently replay every all-archive disputed-pDC raw count row."""
import argparse,ast,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def hist(x):
 a,b=np.unique(np.asarray(x).astype(str),return_counts=True);return dict(zip(a.tolist(),map(int,b)))
def main(root):
 began=time.monotonic();R=Path(load(root/'config/spec_v1.json')['producer']);cfg=load(R/'config/spec_v1.json');assert load(R/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0;terminal=load(R/'provenance/terminal_v1.json')
 for b in cfg['bindings'].values():assert sha(b['path'])==b['sha256']
 ip=R/'outputs/disputed_pDC_context_index_v1.json';sp=R/'outputs/disputed_pDC_context_summary_v1.json';assert sha(ip)==terminal['index_SHA256'] and sha(sp)==terminal['summary_SHA256'];index=load(ip);by_source={r['source']:r for r in index};summary=load(sp)
 tree=ast.parse((R/'scripts/measure_disputed_pDC_native_context_v1.py').read_text());panels=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='CONTEXT_PANELS' for t in n.targets));genes=sorted({g for v in panels.values() for g in v})
 J=Path(cfg['joined']);D=Path(cfg['fine_diagnostic']);sys.path.insert(0,str(J/'scripts'));from joined_native_atlas_delivery_v1 import JoinedNativeAtlas
 atlas=JoinedNativeAtlas(J/'outputs/joined_native_atlas_delivery_v1.json');di={r['source']:r for r in load(D/'outputs/fine_resolution_gap_index_v1.json')};pdc=load(atlas.evidence/'outputs/reference_models_v1.json')['immune_low']['labels'].index('pDC');totals=collections.Counter()
 for ds in atlas.sources:
  a=atlas.source(ds);ca,ea,em=a['ca'],a['ea'],a['em'];de=di[ds];assert sha(de['path'])==de['SHA256'];dm=load(de['path']);ap=Path(de['path']).parent/dm['arrays'];assert sha(ap)==dm['arrays_SHA256']
  with np.load(ap,allow_pickle=False) as f:np.testing.assert_array_equal(f['native_row'],a['rows']);pattern=(f['positive_anchor_bitmask']&(1<<dm['pattern_order'].index('pDC')))!=0;status=f['source_qualification_code'].copy()
  major=np.asarray(a['cm']['candidate_dictionaries']['consensus_major'],object)[ca['consensus_major_code']];fine=np.asarray(a['cm']['candidate_dictionaries']['consensus_celltype'],object)[ca['consensus_celltype_code']];cases=pattern&(ea['immune_low__unique_code']==pdc)&(major=='unresolved');count=int(cases.sum());totals.update(n_sources_examined=1,n_archive_rows_examined=len(cases),n_disputed_case_rows=count)
  if not count:
   assert ds not in by_source;save(root/'checkpoints'/('source_'+ds+'.json'),{'source':ds,'n_sources_checked':totals['n_sources_examined'],'no_disputed_cases':True});atlas._cache.clear();atlas.native._fields.clear();gc.collect();continue
  e=by_source[ds];assert sha(e['path'])==e['SHA256'];m=load(e['path']);npz=Path(e['path']).parent/m['arrays'];assert sha(npz)==m['arrays_SHA256'];assert m['genes']==genes and m['context_panels']==panels
  with np.load(npz,allow_pickle=False) as f:d={k:f[k] for k in f.files}
  masks=[cases,(fine=='pDC')&ca['fit_pool_candidate'],(major=='B')&(ca['native_family_bitmask']==2)&ca['fit_pool_candidate']&~pattern];positions=[];strata=[]
  for j,mask in enumerate(masks):
   members=np.flatnonzero(mask)
   if j==0:take=members
   else:
    k=min(cfg['max_controls_per_source_stratum'],len(members));take=members[np.asarray([i*(len(members)-1)//(k-1) for i in range(k)],np.int64)] if k>1 else members[:k]
   assert m['frames'][j]['n_in_frame']==len(members) and m['frames'][j]['n_measured']==len(take);positions.extend(take.tolist());strata.extend([j]*len(take))
  positions=np.asarray(positions,np.int64);assert len(set(positions))==len(positions);np.testing.assert_array_equal(d['source_native_position'],positions);np.testing.assert_array_equal(d['native_row'],a['rows'][positions]);np.testing.assert_array_equal(d['stratum_code'],strata)
  np.testing.assert_array_equal(d['source_qualification_code'],status[positions])
  for field in ['fit_pool_candidate','parent_source_local_group_code','guarded_source_subset_code','native_family_bitmask','reference_family_bitmask']:np.testing.assert_array_equal(d[field],ca[field][positions])
  lookup=atlas.native.feature_lookup(ds);physical=[next(iter(lookup[g])) if len(lookup.get(g,set()))==1 else None for g in genes];measured=np.asarray([v is not None for v in physical]);np.testing.assert_array_equal(d['measured'],measured);x=np.full_like(d['native_gene_values'],np.nan);rna=np.zeros(len(positions));nnz=0
  key=atlas.native.legacy_config['counts_key'] if atlas.native.is_legacy(ds) else atlas.native.primary.partitions[ds]['counts_key'];csr=atlas.native.store(ds)[key];ncolumns=int(csr.attrs['shape'][1])
  for j,row in enumerate(d['native_row']):
   lo,hi=int(csr['indptr'][int(row)]),int(csr['indptr'][int(row)+1]);v=np.asarray(csr['data'][lo:hi],np.float64);columns=np.asarray(csr['indices'][lo:hi],np.int64);assert np.isfinite(v).all() and (v>=0).all() and np.equal(v,np.floor(v)).all();nnz+=len(v);rna[j]=v.sum();counts=np.bincount(columns,weights=v,minlength=ncolumns);x[j,measured]=counts[np.asarray([v for v in physical if v is not None],np.int64)]
  np.testing.assert_array_equal(x,d['native_gene_values']);np.testing.assert_array_equal(rna,d['native_total_RNA']);np.testing.assert_array_equal(rna,ea['native_total_RNA'][positions])
  shared=[g for g in genes if measured[genes.index(g)] and g in em['review_genes'] and em['gene_measured'][em['review_genes'].index(g)]]
  for g in shared:
   k=em['review_genes'].index(g);np.testing.assert_array_equal(x[:,genes.index(g)]>0,((ea['native_gene_detection_packbits'][positions,k//8]>>(k%8))&1).astype(bool))
  assert m['n_all_disputed_cases_measured']==count and m['n_total_measured']==len(positions) and m['n_controls']==len(positions)-count;assert m['case_source_gate_distribution']==hist(status[cases]);assert m['n_selected_native_nnz']==nnz and m['n_shared_detection_checks']==len(shared)*len(positions)
  totals.update(n_sources_with_disputed_cases=1,n_total_measured=len(positions),n_controls=len(positions)-count,n_selected_native_nnz=nnz,n_shared_detection_checks=len(shared)*len(positions));save(root/'checkpoints'/('source_'+ds+'.json'),{'source':ds,'n_sources_checked':totals['n_sources_examined'],'n_all_disputed_case_rows_checked':totals['n_disputed_case_rows'],'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat()});del a,ca,ea,x,rna,d,major,fine;atlas._cache.clear();atlas.native._fields.clear();gc.collect()
 assert dict(totals)==terminal['totals']==summary['totals'];out={'status':'COMPLETE_ALL_ARCHIVE_DISPUTED_pDC_CASE_SELECTION_NATIVE_COUNTS_AND_BINDINGS_INDEPENDENTLY_ACCEPTED','totals':dict(totals),'producer_selection_and_CSR_helpers_imported':False,'no_identity_QC_or_eligibility_release':True,'producer_terminal_SHA256':sha(R/'provenance/terminal_v1.json'),'summary_SHA256':terminal['summary_SHA256']};op=root/'outputs/independent_disputed_pDC_native_context_audit_v1.json';save(op,out);save(root/'provenance/terminal_v1.json',{'status':out['status'],'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'output_SHA256':sha(op),'producer_terminal_SHA256':out['producer_terminal_SHA256'],'totals':dict(totals),'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException as e:
  import traceback
  save(root/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e),'traceback':traceback.format_exc()});raise
