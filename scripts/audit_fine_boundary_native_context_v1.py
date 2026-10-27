"""Replay every sampled original CSR row without producer aggregation helpers."""
import argparse,ast,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def hist(x):
 a,b=np.unique(np.asarray(x).astype(str),return_counts=True);return dict(zip(a.tolist(),map(int,b)))
def main(root):
 start=time.monotonic();R=Path(load(root/'config/spec_v1.json')['producer']);cfg=load(R/'config/spec_v1.json');assert load(R/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 for b in cfg['bindings'].values():assert sha(b['path'])==b['sha256']
 term=load(R/'provenance/terminal_v1.json');ip=R/'outputs/fine_boundary_context_index_v1.json';sp=R/'outputs/fine_boundary_context_summary_v1.json';assert sha(ip)==term['index_SHA256'] and sha(sp)==term['summary_SHA256']
 index=load(ip);summary=load(sp);assert [r['source'] for r in index]==cfg['sources'];tree=ast.parse((R/'scripts/measure_fine_boundary_native_context_v1.py').read_text());panels=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='PANELS' for t in n.targets));genes=sorted({g for v in panels.values() for g in v})
 J=Path(cfg['joined']);D=Path(cfg['fine_diagnostic']);sys.path.insert(0,str(J/'scripts'));from joined_native_atlas_delivery_v1 import JoinedNativeAtlas
 atlas=JoinedNativeAtlas(J/'outputs/joined_native_atlas_delivery_v1.json');di={r['source']:r for r in load(D/'outputs/fine_resolution_gap_index_v1.json')};labels=load(atlas.evidence/'outputs/reference_models_v1.json')['immune_low']['labels'];pdc_code=labels.index('pDC');totals=collections.Counter()
 for e,sm in zip(index,summary['source_summaries']):
  ds=e['source'];assert sha(e['path'])==e['SHA256'];m=load(e['path']);p=Path(e['path']).parent/m['arrays'];assert sha(p)==m['arrays_SHA256'];assert m['genes']==genes
  with np.load(p,allow_pickle=False) as z:d={k:z[k] for k in z.files}
  a=atlas.source(ds);ca,ea,em=a['ca'],a['ea'],a['em'];de=di[ds];assert sha(de['path'])==de['SHA256'];dm=load(de['path']);dp=Path(de['path']).parent/dm['arrays'];assert sha(dp)==dm['arrays_SHA256']
  with np.load(dp,allow_pickle=False) as z:anchor=z['positive_anchor_bitmask'];status=z['source_qualification_code'];np.testing.assert_array_equal(z['native_row'],a['rows'])
  has=lambda name:(anchor&(1<<dm['pattern_order'].index(name)))!=0
  major=np.asarray(a['cm']['candidate_dictionaries']['consensus_major'],object)[ca['consensus_major_code']];fine=np.asarray(a['cm']['candidate_dictionaries']['consensus_celltype'],object)[ca['consensus_celltype_code']];base=(status==1)&ca['fit_pool_candidate'];endo=(ea['native_panel_bitmask']&(1<<em['native_panel_bit_order'].index('Endothelial')))!=0
  frames=[base&has('pDC')&(ea['immune_low__unique_code']==pdc_code)&(major=='unresolved'),base&(fine=='pDC'),base&(major=='B')&(ca['native_family_bitmask']==2)&~has('pDC'),base&has('Neutrophil')&(major=='Myeloid')&(ca['native_family_bitmask']==4),base&(fine=='Monocyte')&~has('Neutrophil'),base&has('Lymphatic_endothelial_related')&(major=='Endothelial'),base&has('Lymphatic_endothelial_related')&(major=='unresolved'),base&(major=='Endothelial')&endo&~has('Lymphatic_endothelial_related')]
  expected=[];strata=[]
  for j,mask in enumerate(frames):
   members=np.flatnonzero(mask);k=min(len(members),cfg['max_per_stratum']);sel=members[np.asarray([i*(len(members)-1)//(k-1) for i in range(k)],np.int64)] if k>1 else members[:k]
   assert m['strata'][j]['n_source_qualified_fit_records_in_frame']==len(members) and m['strata'][j]['n_sampled']==k;expected.extend(sel.tolist());strata.extend([j]*k)
  np.testing.assert_array_equal(d['source_native_position'],expected);np.testing.assert_array_equal(d['stratum_code'],strata);np.testing.assert_array_equal(d['native_row'],a['rows'][expected]);assert m['n_unique_native_rows']==len(set(expected))
  lookup=atlas.native.feature_lookup(ds);physical=[next(iter(lookup[g])) if len(lookup.get(g,set()))==1 else None for g in genes];measured=np.asarray([v is not None for v in physical]);np.testing.assert_array_equal(d['measured'],measured)
  key=atlas.native.legacy_config['counts_key'] if atlas.native.is_legacy(ds) else atlas.native.primary.partitions[ds]['counts_key'];csr=atlas.native.store(ds)[key];ncolumns=int(csr.attrs['shape'][1]);values=np.full_like(d['native_gene_values'],np.nan);rna=np.zeros(len(expected));nnz=0
  for j,row in enumerate(d['native_row']):
   lo,hi=int(csr['indptr'][int(row)]),int(csr['indptr'][int(row)+1]);v=np.asarray(csr['data'][lo:hi],np.float64);idx=np.asarray(csr['indices'][lo:hi],np.int64);assert np.isfinite(v).all() and (v>=0).all() and np.equal(v,np.floor(v)).all();rna[j]=v.sum();nnz+=len(v);gene_sum=np.bincount(idx,weights=v,minlength=ncolumns);values[j,measured]=gene_sum[np.asarray([i for i in physical if i is not None],np.int64)]
  np.testing.assert_array_equal(d['native_gene_values'],values);np.testing.assert_array_equal(d['native_total_RNA'],rna);np.testing.assert_array_equal(rna,ea['native_total_RNA'][expected]);np.testing.assert_array_equal(d['native_family_bitmask'],ca['native_family_bitmask'][expected]);np.testing.assert_array_equal(d['reference_family_bitmask'],ca['reference_family_bitmask'][expected])
  shared=[g for g in genes if measured[genes.index(g)] and g in em['review_genes'] and em['gene_measured'][em['review_genes'].index(g)]]
  for g in shared:
   i=em['review_genes'].index(g);np.testing.assert_array_equal(values[:,genes.index(g)]>0,((ea['native_gene_detection_packbits'][expected,i//8]>>(i%8))&1).astype(bool))
  author=atlas.native.author_labels(ds);av=np.asarray(author['values'],object)[expected];assert m['author_label_field']==author['field'] and m['author_label_status']==author['status']
  for j,r in enumerate(m['strata']):
   mask=d['stratum_code']==j;v=values[mask];den=rna[mask];assert r==sm['strata'][j]
   assert r['current_coarse_counts']==hist(major[np.asarray(expected)[mask]]) and r['current_fine_counts']==hist(fine[np.asarray(expected)[mask]]) and r['carried_author_label_counts']==hist(av[mask])
   assert r['native_total_RNA_quantiles']==(np.quantile(den,[0,.25,.5,.75,1]).tolist() if mask.any() else None)
   for name,gs in panels.items():
    observed=[g for g in gs if measured[genes.index(g)]];x=v[:,[genes.index(g) for g in observed]]
    calc={'measured_genes':observed,'n_requested':len(gs),'n_sampled':int(mask.sum()),'fraction_at_least2_genes_detected':float(((x>0).sum(axis=1)>=2).mean()) if mask.any() and len(observed)>=2 else None,'fraction_at_least2_genes_each_native_value_ge2':float(((x>=2).sum(axis=1)>=2).mean()) if mask.any() and len(observed)>=2 else None,'per_gene_detection_fraction':{g:float((v[:,genes.index(g)]>0).mean()) for g in observed} if mask.any() else {},'per_gene_median_native_value':{g:float(np.median(v[:,genes.index(g)])) for g in observed} if mask.any() else {},'median_panel_native_RNA_fraction':float(np.median(x.sum(axis=1)/den)) if mask.any() and observed else None}
    assert calc==r['panels'][name],(ds,j,name)
  assert m['n_sampled_native_nnz_observations']==nnz and m['n_original_detection_gene_checks']==len(shared)*len(expected)
  totals.update(n_sources=1,n_sampled_stratum_observations=len(expected),n_unique_native_rows=len(set(expected)),n_sampled_native_nnz_observations=nnz,n_shared_gene_detection_checks=len(shared)*len(expected));save(root/'checkpoints'/('source_'+ds+'.json'),{'source':ds,'n_sources_checked':totals['n_sources'],'n_observations_checked':totals['n_sampled_stratum_observations']});del a,ca,ea,values,rna,d,major,fine;atlas._cache.clear();atlas.native._fields.clear();gc.collect()
 assert dict(totals)==term['totals']==summary['totals'];out={'status':'ALL_FINE_BOUNDARY_SAMPLE_ROWS_ORIGINAL_NATIVE_COUNTS_MASKS_AXES_AND_SUMMARIES_INDEPENDENTLY_ACCEPTED','totals':dict(totals),'producer_sampling_or_CSR_helpers_imported':False,'no_population_accuracy_or_biological_identity_claim':True,'producer_terminal_SHA256':sha(R/'provenance/terminal_v1.json'),'summary_SHA256':term['summary_SHA256']};op=root/'outputs/independent_fine_boundary_native_context_audit_v1.json';save(op,out);save(root/'provenance/terminal_v1.json',{'status':out['status'],'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'output_SHA256':sha(op),'producer_terminal_SHA256':out['producer_terminal_SHA256'],'totals':dict(totals),'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException as e:
  import traceback
  save(root/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e),'traceback':traceback.format_exc()});raise
