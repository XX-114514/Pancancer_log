"""Measure the entire frozen pDC disagreement subset, with bounded controls.

Selection is across the full archive, independently of TME source eligibility.
No label, QC, patient or composition eligibility is assigned by this worker.
"""
import argparse,collections,datetime,gc,hashlib,json,resource,sys,time
from pathlib import Path
import numpy as np
CONTEXT_PANELS={
 'pDC':['CLEC4C','IL3RA','TCF4','GZMB','LILRA4','PLD4','TSPAN13','CCDC50','IRF7'],
 'B_specific':['MS4A1','CD19','CD79B','PAX5','BANK1','CD22'],
 'B_shared_context':['CD79A','CD37','CD74'],
 'Plasma_secretory':['MZB1','JCHAIN','SDC1','DERL3','XBP1','TNFRSF17'],
 'immunoglobulin_context':['IGKC','IGLC2','IGHG1','IGHM'],
 'T':['CD3D','CD3E','CD247','TRAC','TRBC1','TRBC2'],
 'NK':['NKG7','GNLY','KLRD1','NCR1','KLRF1','FCGR3A'],
 'Myeloid':['LYZ','LST1','AIF1','FCER1G','CTSS','FCN1','C1QA','C1QB'],
 'Endothelial':['PECAM1','CDH5','CLDN5','VWF','EMCN'],
 'Fibroblast':['COL1A1','COL1A2','DCN','LUM','PDGFRA','COL3A1'],
 'Glial':['GFAP','AQP4','OLIG2','SOX10','PLP1'],
}
def load(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):Path(p).write_text(json.dumps(d,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def hist(x):
 a,b=np.unique(np.asarray(x).astype(str),return_counts=True);return dict(zip(a.tolist(),map(int,b)))
def main(root):
 began=time.monotonic();cfg=load(root/'config/spec_v1.json');J=Path(cfg['joined']);D=Path(cfg['fine_diagnostic'])
 for b in cfg['bindings'].values():assert sha(b['path'])==b['sha256']
 sys.path.insert(0,str(J/'scripts'));from joined_native_atlas_delivery_v1 import JoinedNativeAtlas
 from native_panel_aggregate_v1 import gather_csr_rows
 from extract_full_source_local_marker_evidence_v1 import storage_gate
 atlas=JoinedNativeAtlas(J/'outputs/joined_native_atlas_delivery_v1.json');di={r['source']:r for r in load(D/'outputs/fine_resolution_gap_index_v1.json')};labels=load(atlas.evidence/'outputs/reference_models_v1.json')['immune_low']['labels'];pdc=labels.index('pDC');genes=sorted({g for v in CONTEXT_PANELS.values() for g in v});index=[];source_rows=[];totals=collections.Counter()
 for ds in atlas.sources:
  storage_gate(root,cfg);s=atlas.source(ds);ca,ea,em=s['ca'],s['ea'],s['em'];record=di[ds];assert sha(record['path'])==record['SHA256'];meta=load(record['path']);ap=Path(record['path']).parent/meta['arrays'];assert sha(ap)==meta['arrays_SHA256']
  with np.load(ap,allow_pickle=False) as a:np.testing.assert_array_equal(a['native_row'],s['rows']);pattern=(a['positive_anchor_bitmask']&(1<<meta['pattern_order'].index('pDC')))!=0;status=a['source_qualification_code'].copy()
  major=np.asarray(s['cm']['candidate_dictionaries']['consensus_major'],object)[ca['consensus_major_code']];fine=np.asarray(s['cm']['candidate_dictionaries']['consensus_celltype'],object)[ca['consensus_celltype_code']]
  cases=pattern&(ea['immune_low__unique_code']==pdc)&(major=='unresolved');count=int(cases.sum());totals.update(n_sources_examined=1,n_archive_rows_examined=len(cases),n_disputed_case_rows=count)
  if count==0:
   save(root/'checkpoints'/('source_'+ds+'.json'),{'source':ds,'n_sources_examined':totals['n_sources_examined'],'n_disputed_case_rows':0});atlas._cache.clear();atlas.native._fields.clear();gc.collect();continue
  frames=[cases,(fine=='pDC')&ca['fit_pool_candidate'],(major=='B')&(ca['native_family_bitmask']==2)&ca['fit_pool_candidate']&~pattern];names=['all_disputed_pDC_pattern_unique_reference_unresolved','current_pDC_fine_controls','B_concordant_controls'];positions=[];strata=[];frame_records=[]
  for j,mask in enumerate(frames):
   members=np.flatnonzero(mask);take=members if j==0 else (members[np.linspace(0,len(members)-1,min(cfg['max_controls_per_source_stratum'],len(members)),dtype=np.int64)] if len(members) else members);positions.extend(take.tolist());strata.extend([j]*len(take));frame_records.append({'name':names[j],'n_in_frame':len(members),'n_measured':len(take),'all_frame_rows_measured':j==0})
  positions=np.asarray(positions,np.int64);strata=np.asarray(strata,np.uint8);rows=s['rows'][positions];assert len(set(positions))==len(positions)
  lookup=atlas.native.feature_lookup(ds);physical=[next(iter(lookup[g])) if len(lookup.get(g,set()))==1 else None for g in genes];measured=np.asarray([v is not None for v in physical]);x=np.full((len(positions),len(genes)),np.nan,np.float64);rna=np.zeros(len(positions));nnz=0
  key=atlas.native.legacy_config['counts_key'] if atlas.native.is_legacy(ds) else atlas.native.primary.partitions[ds]['counts_key'];csr=atlas.native.store(ds)[key]
  for low in range(0,len(positions),cfg['read_block_size']):
   hi=min(len(positions),low+cfg['read_block_size']);raw=gather_csr_rows(csr,rows[low:hi],max_span=cfg['max_native_row_span']);assert np.isfinite(raw.data).all() and (raw.data>=0).all() and np.equal(raw.data,np.floor(raw.data)).all();rna[low:hi]=np.asarray(raw.sum(axis=1,dtype=np.float64)).ravel();x[low:hi,measured]=raw[:,np.asarray([v for v in physical if v is not None],np.int64)].toarray();nnz+=raw.nnz
  np.testing.assert_array_equal(rna,ea['native_total_RNA'][positions]);assert (rna>0).all()
  shared=[g for g in genes if measured[genes.index(g)] and g in em['review_genes'] and em['gene_measured'][em['review_genes'].index(g)]]
  for g in shared:
   k=em['review_genes'].index(g);np.testing.assert_array_equal(x[:,genes.index(g)]>0,((ea['native_gene_detection_packbits'][positions,k//8]>>(k%8))&1).astype(bool))
  arrays={'source_native_position':positions,'native_row':rows,'stratum_code':strata,'native_gene_values':x,'measured':measured,'native_total_RNA':rna,'source_qualification_code':status[positions],'fit_pool_candidate':ca['fit_pool_candidate'][positions],'parent_source_local_group_code':ca['parent_source_local_group_code'][positions],'guarded_source_subset_code':ca['guarded_source_subset_code'][positions],'native_family_bitmask':ca['native_family_bitmask'][positions],'reference_family_bitmask':ca['reference_family_bitmask'][positions]};npz=root/'outputs'/(ds+'_disputed_pDC_context_v1.npz')
  with npz.open('xb') as f:np.savez_compressed(f,**arrays)
  with np.load(npz,allow_pickle=False) as f:
   for k,v in arrays.items():np.testing.assert_array_equal(f[k],v)
  desc={'source':ds,'n_archive_rows':len(cases),'n_all_disputed_cases_measured':count,'n_controls':len(positions)-count,'n_total_measured':len(positions),'frames':frame_records,'genes':genes,'context_panels':CONTEXT_PANELS,'arrays':npz.name,'arrays_SHA256':sha(npz),'n_shared_detection_checks':len(shared)*len(positions),'n_selected_native_nnz':nnz,'case_source_gate_distribution':hist(status[cases]),'all_case_rows_including_archive_model_exclusions_measured':True,'source_gate_QC_old_labels_and_old_CM_unchanged':True,'extra_genes_not_independent_of_reference_models':True};mp=npz.with_suffix('.json');save(mp,desc);index.append({'source':ds,'path':str(mp),'SHA256':sha(mp)});source_rows.append({k:desc[k] for k in ['source','n_archive_rows','n_all_disputed_cases_measured','n_controls','n_total_measured','case_source_gate_distribution']});totals.update(n_sources_with_disputed_cases=1,n_total_measured=len(positions),n_controls=len(positions)-count,n_selected_native_nnz=nnz,n_shared_detection_checks=len(shared)*len(positions));save(root/'checkpoints'/('source_'+ds+'.json'),{'source':ds,'completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'n_sources_examined':totals['n_sources_examined'],'n_disputed_case_rows':count,'total_measured':totals['n_total_measured']});del s,ca,ea,x,rna,arrays,major,fine;atlas._cache.clear();atlas.native._fields.clear();gc.collect()
  storage_gate(root,cfg)
 assert totals['n_archive_rows_examined']==sum(r['n_archive_rows'] for r in atlas.native.archive_partitions.values())
 ip=root/'outputs/disputed_pDC_context_index_v1.json';sp=root/'outputs/disputed_pDC_context_summary_v1.json';save(ip,index);save(sp,{'totals':dict(totals),'sources':source_rows,'case_selection_entire_archive_not_only_TME_qualified':True,'no_new_labels_or_QC':True});save(root/'provenance/terminal_v1.json',{'status':'COMPLETE_ALL_ARCHIVE_DISPUTED_pDC_NATIVE_CONTEXT_MEASURED_INDEPENDENT_AUDIT_REQUIRED','completed_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'totals':dict(totals),'index_SHA256':sha(ip),'summary_SHA256':sha(sp),'wall_seconds':time.monotonic()-began,'peak_RSS_bytes_Linux':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,'no_identity_or_QC_or_source_gate_release':True})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',required=True);root=Path(p.parse_args().output)
 try:main(root)
 except BaseException as e:
  import traceback
  save(root/'provenance/failure_v1.json',{'error_class':type(e).__name__,'error':str(e),'traceback':traceback.format_exc()});raise
