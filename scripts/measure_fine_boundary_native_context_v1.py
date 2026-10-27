"""Frozen source/stratum native-count context for fine-resolution boundaries.

This targeted diagnostic preserves full-archive rows and labels. The design is
not a prevalence or accuracy sample. Added genes may also be model features.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
PANELS={
 'pDC_context':['CLEC4C','IL3RA','TCF4','GZMB','LILRA4','PLD4','IRF7','CCDC50','TSPAN13'],
 'B_context':['CD79A','CD79B','CD19','MS4A1','PAX5','BANK1','CD22','CD37','CD74'],
 'Plasma_context':['JCHAIN','MZB1','SDC1','XBP1','DERL3'],
 'transitional_DC_context':['AXL','SIGLEC6','CD1C','FCER1A','CLEC10A'],
 'Neutrophil_context':['CSF3R','FCGR3B','S100A8','S100A9','CEACAM8','CXCR2','MMP8','CAMP','LCN2','S100A12'],
 'Monocyte_context':['FCN1','LYZ','CD14','FCGR1A','LST1','CTSS','HLA-DRA','TYROBP'],
 'Lymphatic_context':['PROX1','FLT4','PDPN','CCL21','LYVE1','MMRN1'],
 'pan_endothelial_context':['PECAM1','VWF','KDR','EMCN','CDH5','CLDN5','ESAM','ERG'],
 'glial_context':['OLIG1','OLIG2','SOX10','PLP1','GFAP','AQP4','SLC1A3'],
 'fibroblast_context':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA'],
 'T_context':['CD3D','CD3E','TRAC','CD247'],
 'Macrophage_context':['C1QA','C1QB','APOE','CD68','MS4A7'],
}
def load(p):return json.loads(Path(p).read_text())
def save(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def hist(x):
 a,b=np.unique(np.asarray(x).astype(str),return_counts=True);return dict(zip(a.tolist(),map(int,b)))

def main(root):
 started=time.monotonic();cfg=load(root/'config/spec_v1.json');P=root.parent;J=Path(cfg['joined']);D=Path(cfg['fine_diagnostic']);A=Path(cfg['fine_audit'])
 for b in cfg['bindings'].values():assert sha(b['path'])==b['sha256']
 assert load(A/'provenance/supervisor_terminal_v1.json')['actual_child_exit_code']==0
 assert load(A/'provenance/terminal_v1.json')['producer_terminal_SHA256']==sha(D/'provenance/terminal_v1.json')
 sys.path.insert(0,str(J/'scripts'));from joined_native_atlas_delivery_v1 import JoinedNativeAtlas
 from native_panel_aggregate_v1 import gather_csr_rows
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 atlas=JoinedNativeAtlas(J/'outputs/joined_native_atlas_delivery_v1.json');index={r['source']:r for r in load(D/'outputs/fine_resolution_gap_index_v1.json')}
 genes=sorted(set(g for v in PANELS.values() for g in v));summaries=[];outputs=[];total=collections.Counter()
 for ds in cfg['sources']:
  storage_gate(root,cfg);s=atlas.source(ds);ca,ea,em=s['ca'],s['ea'],s['em'];n=len(s['rows']);entry=index[ds];assert sha(entry['path'])==entry['SHA256'];meta=load(entry['path']);p=Path(entry['path']).parent/meta['arrays'];assert sha(p)==meta['arrays_SHA256']
  with np.load(p,allow_pickle=False) as f:anchor=f['positive_anchor_bitmask'];q=f['source_qualification_code'];np.testing.assert_array_equal(f['native_row'],s['rows'])
  bit=lambda name:(anchor&(1<<meta['pattern_order'].index(name)))!=0
  major=np.asarray(s['cm']['candidate_dictionaries']['consensus_major'],object)[ca['consensus_major_code']];fine=np.asarray(s['cm']['candidate_dictionaries']['consensus_celltype'],object)[ca['consensus_celltype_code']]
  ref_models=load(atlas.evidence/'outputs/reference_models_v1.json');pdc_code=ref_models['immune_low']['labels'].index('pDC');base=(q==1)&ca['fit_pool_candidate']
  endo_panel=(ea['native_panel_bitmask']&(1<<em['native_panel_bit_order'].index('Endothelial')))!=0
  frames={
   'pDC_pattern_unique_reference_current_unresolved':base&bit('pDC')&(ea['immune_low__unique_code']==pdc_code)&(major=='unresolved'),
   'pDC_current_fine_control':base&(fine=='pDC'),
   'B_concordant_control':base&(major=='B')&(ca['native_family_bitmask']==2)&~bit('pDC'),
   'Neutrophil_pattern_current_Myeloid':base&bit('Neutrophil')&(major=='Myeloid')&(ca['native_family_bitmask']==4),
   'Monocyte_current_fine_control':base&(fine=='Monocyte')&~bit('Neutrophil'),
   'Lymphatic_pattern_current_Endothelial':base&bit('Lymphatic_endothelial_related')&(major=='Endothelial'),
   'Lymphatic_pattern_current_unresolved':base&bit('Lymphatic_endothelial_related')&(major=='unresolved'),
   'pan_Endothelial_control_without_Lymphatic_pattern':base&(major=='Endothelial')&endo_panel&~bit('Lymphatic_endothelial_related')}
  positions=[];strata=[];records=[]
  for j,(name,mask) in enumerate(frames.items()):
   members=np.flatnonzero(mask);take=members[np.linspace(0,len(members)-1,min(cfg['max_per_stratum'],len(members)),dtype=np.int64)] if len(members) else members
   positions.extend(take.tolist());strata.extend([j]*len(take));records.append({'stratum':name,'n_source_qualified_fit_records_in_frame':len(members),'n_sampled':len(take)})
  positions=np.asarray(positions,np.int64);strata=np.asarray(strata,np.uint8)
  # Overlapping biological patterns may select the same cell in two strata.
  # Preserve that relation; observation counts are not independent-cell counts.
  n_unique_positions=len(np.unique(positions))
  native_rows=s['rows'][positions];lookup=atlas.native.feature_lookup(ds);physical=[next(iter(lookup[g])) if len(lookup.get(g,set()))==1 else None for g in genes];measured=np.asarray([v is not None for v in physical]);values=np.full((len(positions),len(genes)),np.nan,np.float64);RNA=np.empty(len(positions),np.float64)
  key=atlas.native.legacy_config['counts_key'] if atlas.native.is_legacy(ds) else atlas.native.primary.partitions[ds]['counts_key'];csr=atlas.native.store(ds)[key];nnz=0
  for low in range(0,len(positions),cfg['read_block_size']):
   high=min(len(positions),low+cfg['read_block_size']);raw=gather_csr_rows(csr,native_rows[low:high],max_span=cfg['max_native_row_span']);assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.equal(raw.data,np.floor(raw.data)).all();nnz+=raw.nnz;RNA[low:high]=np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel();values[low:high,measured]=raw[:,np.asarray([p for p in physical if p is not None],np.int64)].toarray()
  np.testing.assert_array_equal(RNA,ea['native_total_RNA'][positions]);assert np.isnan(values[:,~measured]).all()
  old_genes=em['review_genes'];old_measure=np.asarray(em['gene_measured'],bool);shared=[g for g,p in zip(genes,physical) if p is not None and g in old_genes and old_measure[old_genes.index(g)]]
  for g in shared:
   k=old_genes.index(g);expected=((ea['native_gene_detection_packbits'][positions,k//8]>>(k%8))&1).astype(bool);np.testing.assert_array_equal(values[:,genes.index(g)]>0,expected)
  author=atlas.native.author_labels(ds);author_values=np.asarray(author['values'],object)[positions]
  for j,record in enumerate(records):
   z=strata==j;v=values[z];rna=RNA[z];panels={}
   for name,panel in PANELS.items():
    observed=[g for g in panel if measured[genes.index(g)]];cols=[genes.index(g) for g in observed];x=v[:,cols]
    panels[name]={'measured_genes':observed,'n_requested':len(panel),'n_sampled':int(z.sum()),'fraction_at_least2_genes_detected':float(((x>0).sum(axis=1)>=2).mean()) if z.any() and len(cols)>=2 else None,'fraction_at_least2_genes_each_native_value_ge2':float(((x>=2).sum(axis=1)>=2).mean()) if z.any() and len(cols)>=2 else None,'per_gene_detection_fraction':{g:float((v[:,genes.index(g)]>0).mean()) for g in observed} if z.any() else {},'per_gene_median_native_value':{g:float(np.median(v[:,genes.index(g)])) for g in observed} if z.any() else {},'median_panel_native_RNA_fraction':float(np.median(x.sum(axis=1)/rna)) if z.any() and len(cols) else None}
   record.update(native_total_RNA_quantiles=np.quantile(rna,[0,.25,.5,.75,1]).tolist() if z.any() else None,current_coarse_counts=hist(major[positions[z]]),current_fine_counts=hist(fine[positions[z]]),carried_author_label_counts=hist(author_values[z]),panels=panels)
  arrays={'source_native_position':positions,'native_row':native_rows,'stratum_code':strata,'native_gene_values':values,'measured':measured,'native_total_RNA':RNA,'native_family_bitmask':ca['native_family_bitmask'][positions],'reference_family_bitmask':ca['reference_family_bitmask'][positions]};npz=root/'outputs'/(ds+'_fine_boundary_context_v1.npz')
  with npz.open('xb') as f:np.savez_compressed(f,**arrays)
  with np.load(npz,allow_pickle=False) as f:
   for k,v in arrays.items():np.testing.assert_array_equal(f[k],v)
  mp=npz.with_suffix('.json');report={'source':ds,'genes':genes,'n_source_archive_rows':n,'n_sampled_stratum_observations':len(positions),'n_unique_native_rows':n_unique_positions,'overlapping_strata_not_counted_as_independent_cells':True,'stratum_order':list(frames),'strata':records,'arrays':npz.name,'arrays_SHA256':sha(npz),'author_label_field':author['field'],'author_label_status':author['status'],'n_original_detection_gene_checks':len(shared)*len(positions),'n_sampled_native_nnz_observations':nnz,'max_native_read_span':cfg['max_native_row_span'],'added_genes_may_be_reference_model_features_not_independent_validation':True,'sampling_is_targeted_not_prevalence_or_accuracy_estimation':True,'no_identity_QC_or_source_qualification_change':True};save(mp,report);outputs.append({'source':ds,'path':str(mp),'SHA256':sha(mp)});summaries.append({'source':ds,'strata':records});total.update(n_sources=1,n_sampled_stratum_observations=len(positions),n_unique_native_rows=n_unique_positions,n_sampled_native_nnz_observations=nnz,n_shared_gene_detection_checks=len(shared)*len(positions));save(root/'checkpoints'/('source_'+ds+'.json'),{'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_completed_sources':len(outputs),'source':ds})
  del s,ca,ea,major,fine,values,RNA,arrays,author;atlas._cache.clear();atlas.native._fields.clear();gc.collect()
 ip=root/'outputs/fine_boundary_context_index_v1.json';sp=root/'outputs/fine_boundary_context_summary_v1.json';save(ip,outputs);save(sp,{'totals':dict(total),'source_summaries':summaries,'no_labels_changed':True})
 save(root/'provenance/terminal_v1.json',{'status':'TARGETED_FINE_BOUNDARY_NATIVE_CONTEXT_MEASURED_INDEPENDENT_RAW_READ_AUDIT_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(total),'index_SHA256':sha(ip),'summary_SHA256':sha(sp),'wall_seconds':time.monotonic()-started,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'no_identity_QC_or_source_qualification_change':True})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException as e:
  import traceback
  save(root/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e),'traceback':traceback.format_exc()});raise
