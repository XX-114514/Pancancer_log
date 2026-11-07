"""Freeze same-V4 reference and all existing qualified human Visium sections."""
import json,hashlib,gc
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd,scipy.sparse as sp,zarr,h5py,anndata as ad
from anndata.io import read_elem

r=PROJECT_ROOT;run=r/'runs/20261107_v4_Visium_author_validation_v2';run.mkdir(exist_ok=False)
for name in ['reference','inputs','scripts','config','launch','sections','scores','audit','author']:(run/name).mkdir()
root=r/'annotation_rebuild/20260923_core_anndata_zarr_v1/objects';g=zarr.open_group(str(root/'pancancer_core_fullgene_v1.zarr'),mode='r');axis=pd.Index(read_elem(g['obs/_index']));var=pd.Index(read_elem(g['var/_index']));src=np.asarray(read_elem(g['obs/gse_id'])).astype(str)
annpath=r/'runs/20261106_v7_publication_marker_annotation_v1/stage_h_selective_subtype_restore_v4/annotation_frozen_V7_selective_subtype_20261106_v4_PRIVATE.pkl.gz';ann=pd.read_pickle(annpath);assert hashlib.sha256(annpath.read_bytes()).hexdigest()=='2a4bc82cb495d0fc2ffb40ecda6b1cbb53802eb2bf87efd2392746371653244c';pos=axis.get_indexer(ann.index);assert (pos>=0).all();available={}
core={'IGHD','TCL1A','MS4A1','CD79A','PECAM1','VWF','XCL1','KLRD1','GNLY','SPP1','EPCAM','KRT18'}
for source in sorted(set(src)):
 with h5py.File(root/'harmonized_per_gse'/(source+'_symbol_counts.h5ad'),'r') as f:available[source]=set(map(str,read_elem(f['var']).index))
allowed=sorted(s for s,a in available.items() if core<=a);genes=sorted(set.intersection(*(available[s] for s in allowed))-set(x for x in var if str(x).startswith('MT-')));assert core<=set(genes);take=np.isin(src[pos],allowed);f=ann.loc[take].copy();f['source']=src[pos[take]];assert f.final_annotation.nunique()==ann.final_annotation.nunique()==50
rng=np.random.default_rng(2026110701);selected=[];meta=[];catalog=[]
for label,t in f.groupby('final_annotation',observed=True,sort=True):
 n=len(t);nref=min(1000,n);nhold=min(200,n-nref);chosen=rng.choice(n,nref+nhold,replace=False);p=axis.get_indexer(t.index[chosen]);selected.extend(p);split=['reference']*nref+['heldout']*nhold
 for loc,s in zip(p,split):meta.append({'cell_id':axis[loc],'row':int(loc),'node':str(label),'source':src[loc],'split':s,'major':str(t.major_cell_type.iloc[0])})
 catalog.append({'node':str(label),'major':str(t.major_cell_type.iloc[0]),'n_V4_cells_in_measured_reference_sources':n,'n_reference_cells':nref,'n_heldout_cells':nhold,'n_reference_sources':len(set(src[p[:nref]]))})
meta=pd.DataFrame(meta).sort_values('row').set_index('cell_id');selected=meta.row.to_numpy();columns=var.get_indexer(genes);assert (columns>=0).all();ptr=g['X/indptr'][:];blocks=[];N=len(axis)
for left in range(0,N,20000):
 right=min(N,left+20000);lo=np.searchsorted(selected,left);hi=np.searchsorted(selected,right)
 if hi==lo:continue
 a,b=int(ptr[left]),int(ptr[right]);data=g['X/data'][a:b];ind=g['X/indices'][a:b];assert np.isfinite(data).all() and (data>=0).all() and np.allclose(data,np.rint(data));x=sp.csr_matrix((data,ind,ptr[left:right+1]-a),shape=(right-left,len(var)));blocks.append(x[selected[lo:hi]-left][:,columns].astype(np.float32));(run/'reference/progress.json').write_text(json.dumps({'UTC':datetime.now(timezone.utc).isoformat(),'raw_axis_rows_scanned':right,'n_selected_rows_completed':int(hi),'n_selected_total':len(meta)}));gc.collect()
X=sp.vstack(blocks,format='csr');assert X.shape==(len(meta),len(genes));assert (np.asarray(X.sum(axis=1)).ravel()>0).all();a=ad.AnnData(X,obs=meta[['node','source','split','major']].copy(),var=pd.DataFrame(index=pd.Index(genes,name='gene')));a.write_h5ad(run/'reference/reference_and_heldout_PRIVATE.h5ad',compression='gzip');meta.to_csv(run/'reference/selected_cell_rows_PRIVATE.tsv.gz',sep='\t');pd.DataFrame(catalog).to_csv(run/'reference/node_reference_coverage.tsv',sep='\t',index=False);pd.DataFrame({'source':sorted(available),'core_genes_measured':[s in allowed for s in sorted(available)],'native_gene_count':[len(available[s]) for s in sorted(available)]}).to_csv(run/'reference/source_measurement_qualification.tsv',sep='\t',index=False)
# Freeze observed canonical inputs. Missing expected catalog files are not fabricated.
pan=r.parent;qcroot=pan/'Project_v2/Downstream/0716_Data/raw_geo_supplements/QC_adata';jobs=[];scope=[];private=[]
def add_file(path,dataset):
 with h5py.File(path,'r') as h:
  obs=read_elem(h['obs']);v=read_elem(h['var']);layer='counts' if 'layers/counts' in h else 'X';keys=list(h['obsm'].keys());coord='obsm/spatial' if 'spatial' in keys else 'obs/pxl_col_in_fullres,pxl_row_in_fullres' if {'pxl_col_in_fullres','pxl_row_in_fullres'}<=set(obs) else 'obs/pxl_col,pxl_row' if {'pxl_col','pxl_row'}<=set(obs) else None
  assert v.index.is_unique
  group='gsm_id' if 'gsm_id' in obs else 'sample_id' if 'sample_id' in obs else None;groups=sorted(obs[group].astype(str).unique()) if group else ['single_section']
  for unit in groups:
   mask=obs[group].astype(str).eq(unit).to_numpy() if group else np.ones(len(obs),bool);t=obs.loc[mask];text=' '.join([path.name]+[str(x) for c in ['sample_id','sample_label','cancer_type'] if c in t for x in t[c].astype(str).unique()]);reason='eligible' if coord and len(t)>60 else 'coordinates_missing_or_fewer_than_61_spots'
   if dataset=='GSE203612' and 'PDX' in text.upper():reason='PDX_not_human_native_TME'
   key='S%03d'%(len(private)+1);entry={'section_code':key,'dataset':dataset,'input_path':str(path),'group_column':group,'group_value':unit,'counts_layer':layer,'coordinates':coord,'n_spots_before_new_filter':len(t),'input_genes':len(v),'n_common_reference_genes':len(set(v.index)&set(genes)),'status':reason,'input_file_bytes':path.stat().st_size};private.append(entry)
   scope.append({k:entry[k] for k in ['section_code','dataset','n_spots_before_new_filter','input_genes','n_common_reference_genes','status']})
   if reason=='eligible':jobs.append(entry)
for dataset in ['GSE203612','GSE175540']:
 folder=pan/'Project_v5/spatial_expansion_20260903/processed'/dataset
 for path in sorted(folder.glob('*.h5ad')):add_file(path,dataset)
for dataset in ['GSE273952','GSE238264','GSE301720']:
 path=qcroot/(dataset+'_QC_adata.h5ad')
 if path.exists():add_file(path,dataset)
 else:scope.append({'section_code':'not_available','dataset':dataset,'status':'canonical_QC_h5ad_not_present'})
assert jobs and len({x['section_code'] for x in jobs})==len(jobs);(run/'inputs/section_jobs_PRIVATE.json').write_text(json.dumps(jobs,indent=2));(run/'inputs/all_input_scope_PRIVATE.json').write_text(json.dumps(private,indent=2));pd.DataFrame(scope).to_csv(run/'inputs/all_section_technical_scope.tsv',sep='\t',index=False);pd.DataFrame(scope).groupby(['dataset','status'],dropna=False).agg(n_sections=('section_code','size'),n_spots=('n_spots_before_new_filter','sum')).reset_index().to_csv(run/'inputs/source_scope_summary.tsv',sep='\t',index=False)
receipt={'status':'SAME_V4_REFERENCE_AND_OBSERVED_VISIUM_SCOPE_PREPARED_PENDING_AUDIT','UTC':datetime.now(timezone.utc).isoformat(),'annotation_SHA256':hashlib.sha256(annpath.read_bytes()).hexdigest(),'reference_labels':len(catalog),'reference_cells':sum(x['n_reference_cells'] for x in catalog),'heldout_cells':sum(x['n_heldout_cells'] for x in catalog),'reference_sources':len(allowed),'reference_common_measured_genes':len(genes),'reference_source_rule':'core identity genes all measured; full common measured gene set; no spatial effects used','core_identity_genes':sorted(core),'CM_nodes_or_original_annotation_changed':False,'eligible_sections':len(jobs),'eligible_sources':len(set(x['dataset'] for x in jobs)),'patient_independence_certified':False,'all_label_identifiability_claim':False,'files':{str(p.relative_to(run)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in [run/'reference/reference_and_heldout_PRIVATE.h5ad',run/'reference/selected_cell_rows_PRIVATE.tsv.gz',run/'inputs/section_jobs_PRIVATE.json']}};(run/'inputs/complete.json').write_text(json.dumps(receipt,indent=2));print(json.dumps({k:v for k,v in receipt.items() if k!='files'}))
