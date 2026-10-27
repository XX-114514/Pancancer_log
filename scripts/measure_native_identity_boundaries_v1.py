"""Bounded quantitative native-RNA boundary checks, with fixed sampling rules.

Panel coexpression is descriptive and cannot diagnose doublets, contamination,
malignancy or a cell identity. No labels, QC calls or clinical models change.
"""
import argparse, collections, datetime, hashlib, json, resource, sys, time, gc
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from review_source_qualified_identity_context_v1 import source_status, load, save, sha

PANELS={
 'T':['CD3D','CD3E','TRAC','CD247','TRBC1','TRBC2','LCK','LAT','ZAP70'],
 'NK':['NKG7','GNLY','KLRD1','FCGR3A','NCR1','KLRF1'],
 'Myeloid':['LST1','LYZ','AIF1','FCER1G','CTSS','TYROBP','SPI1','CD68','MS4A7'],
 'Epithelial':['EPCAM','KRT8','KRT18','KRT19','KRT7','KRT14','KRT5','SFN'],
 'Fibroblast':['COL1A1','COL1A2','DCN','LUM','PDGFRA','COL3A1','COL6A1','COL6A2'],
 'Mural':['RGS5','CSPG4','PDGFRB','MCAM','MYH11','CNN1','TAGLN','ACTA2'],
}

def main(root):
 start=time.monotonic();cfg=load(root/'config/spec_v1.json');J=Path(cfg['joined'])
 for b in cfg['bindings'].values():assert sha(b['path'])==b['sha256']
 sys.path.insert(0,str(J/'scripts'))
 from joined_native_atlas_delivery_v1 import JoinedNativeAtlas
 from native_panel_aggregate_v1 import gather_csr_rows
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 atlas=JoinedNativeAtlas(J/'outputs/joined_native_atlas_delivery_v1.json');view=load(cfg['source_view']);core=load(view['work_units']['path'])['units']
 source_n=collections.Counter()
 for u in core:
  if u['source_qualification']=='来源资格通过':source_n[u['source']]+=u['n_archive_rows']
 ordered=sorted(source_n,key=lambda ds:(source_n[ds],ds));selected=[ordered[i] for i in np.linspace(0,len(ordered)-1,min(8,len(ordered)),dtype=int)]
 assert selected==cfg['selected_sources_by_size_spaced_not_phenotype']
 genes=sorted(set(g for panel in PANELS.values() for g in panel));index=[];summary=[];counts=collections.Counter()
 for ds in selected:
  storage_gate(root,cfg);a=atlas.source(ds);ca=a['ca'];eligible=(source_status(a,core,view)==1)&ca['fit_pool_candidate'];native=ca['native_family_bitmask'];ref=ca['reference_family_bitmask'];dic=a['cm']['candidate_dictionaries'];major=np.asarray(dic['consensus_major'],object)[ca['consensus_major_code']]
  strata={'T_NK_Myeloid_native_conflict':eligible&(native==5),
    'Epithelial_Myeloid_native_conflict':eligible&(native==132),
    'Epithelial_Fibroblast_native_conflict':eligible&(native==160),
    'T_NK_within_family_ambiguity':eligible&(major=='T_NK'),
    'T_concordant_control':eligible&(native==1)&(ref==1)&(major=='T'),
    'Myeloid_concordant_control':eligible&(native==4)&(ref==4)&(major=='Myeloid'),
    'Epithelial_concordant_control':eligible&(native==128)&(ref==128)&(major=='Epithelial')}
  positions=[];codes=[];descriptors=[]
  for j,(name,mask) in enumerate(strata.items()):
   choices=np.flatnonzero(mask);pick=choices[np.linspace(0,len(choices)-1,min(32,len(choices)),dtype=int)] if len(choices) else choices
   positions.extend(pick.tolist());codes.extend([j]*len(pick));descriptors.append({'stratum':name,'n_source_qualified_fit_records_in_frame':len(choices),'n_sampled':len(pick)})
  p=np.asarray(positions,np.int64);code=np.asarray(codes,np.int16)
  assert len(p)==len(set(p)) and len(p)<=224
  key=atlas.native.legacy_config['counts_key'] if atlas.native.is_legacy(ds) else atlas.native.primary.partitions[ds]['counts_key'];store=atlas.native.store(ds)
  raw=gather_csr_rows(store[key],a['rows'][p],max_span=128);tot=np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel()
  np.testing.assert_allclose(tot,a['ea']['native_total_RNA'][p],rtol=1e-10,atol=1e-6)
  assert (tot>0).all() and np.isfinite(raw.data).all() and (raw.data>=0).all()
  lookup=atlas.native.feature_lookup(ds);cols=[next(iter(lookup[g])) if len(lookup.get(g,set()))==1 else None for g in genes]
  measured=np.asarray([c is not None for c in cols]);physical=[c for c in cols if c is not None];assert len(set(physical))==len(physical)
  panel=sp.hstack([raw[:,c] if c is not None else sp.csr_matrix((len(p),1)) for c in cols],format='csr').toarray()
  arrays={'source_native_position':p,'native_row':a['rows'][p],'stratum_code':code,'native_total_RNA':tot,
          'measured':measured,'native_gene_values':panel,'native_family_bitmask':native[p],'reference_family_bitmask':ref[p]}
  for name,gs in PANELS.items():
   q=[genes.index(g) for g in gs if measured[genes.index(g)]];counts_here=panel[:,q]
   arrays[name+'__n_genes_detected']=(counts_here>0).sum(axis=1)
   arrays[name+'__n_genes_value_ge2']=(counts_here>=2).sum(axis=1)
   arrays[name+'__count_slot_sum']=counts_here.sum(axis=1)
   arrays[name+'__fraction_of_native_RNA']=counts_here.sum(axis=1)/tot
  out=root/'outputs'/(ds+'_native_identity_boundary_probe_v1.npz')
  with out.open('xb') as h:np.savez_compressed(h,**arrays)
  with np.load(out,allow_pickle=False) as f:
   for k,v in arrays.items():np.testing.assert_array_equal(f[k],v)
  result=[]
  for j,d in enumerate(descriptors):
   m=code==j;record={**d,'native_total_RNA_quantiles':np.quantile(tot[m],[0,.25,.5,.75,1]).tolist() if m.any() else None,'panels':{}}
   for name,gs in PANELS.items():
    available=[g for g in gs if measured[genes.index(g)]]
    record['panels'][name]={'measured_genes':available,'n_measured':len(available),'n_requested':len(gs),
       'fraction_with_at_least2_detected_genes':float((arrays[name+'__n_genes_detected'][m]>=2).mean()) if m.any() and len(available)>=2 else None,
       'fraction_with_at_least2_genes_each_value_ge2':float((arrays[name+'__n_genes_value_ge2'][m]>=2).mean()) if m.any() and len(available)>=2 else None,
       'median_native_RNA_fraction':float(np.median(arrays[name+'__fraction_of_native_RNA'][m])) if m.any() and available else None}
   result.append(record)
  metadata={'source':ds,'sampling_rule':'eight size-spaced source-qualified sources; up to32 axis-spaced rows per fixed stratum; no clinical labels',
      'genes':genes,'measured':measured.tolist(),'strata':descriptors,'panels':PANELS,'arrays':out.name,'arrays_SHA256':sha(out),
      'count_slot_values_integer_in_this_probe':bool(np.equal(raw.data,np.floor(raw.data)).all()),
      'not_independent_truth_and_no_doublet_or_ambient_diagnosis':True,'identity_and_QC_unchanged':True,'summaries':result}
  mp=out.with_suffix('.json');save(mp,metadata);index.append({'source':ds,'metadata':str(mp),'SHA256':sha(mp)});summary.append({'source':ds,'results':result});counts.update(n_sources=1,n_sampled_rows=len(p))
  save(root/'checkpoints'/(ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_completed_sources':len(index)})
  del a,ca,raw,panel,arrays;atlas._cache.clear();atlas.native._fields.clear();gc.collect()
 save(root/'outputs/native_identity_boundary_probe_index_v1.json',index);save(root/'outputs/native_identity_boundary_probe_summary_v1.json',{'counts':dict(counts),'source_results':summary,'fixed_diagnostic_thresholds_not_new_labels':True})
 save(root/'provenance/terminal_v1.json',{'status':'BOUNDED_NATIVE_IDENTITY_BOUNDARY_QUANTITATIVE_PROBE_COMPLETE_NO_RELABELING','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'counts':dict(counts),'index_SHA256':sha(root/'outputs/native_identity_boundary_probe_index_v1.json'),'summary_SHA256':sha(root/'outputs/native_identity_boundary_probe_summary_v1.json'),'wall_seconds':time.monotonic()-start,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException as e:
  import traceback
  save(root/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e),'traceback':traceback.format_exc()});raise
