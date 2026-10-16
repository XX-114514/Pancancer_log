"""Existing whole-transcriptome spatial asset, raw counts/axes only, no clinical reads."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd,h5py
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);project=AR.parent.parent;old=project/'Project_v2/Downstream/0716_Data/immunotherapy_4cohort_analysis/20260803_v1';assert old.is_dir(),old
# annotation_rebuild.parent is Project_v3; its parent is PanCancer.
if not old.exists():raise RuntimeError('Historical root not found')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
contractp=old/'outputs/spatial/GSE273952_spatial_contract_audit.json';contract=json.loads(contractp.read_text());path=Path(contract['input_h5ad']);assert path.is_file();O=V/'existing_GSE273952_Visium_raw_preflight_20261016_v2';assert not O.exists();O.mkdir(mode=0o700);dictionary=V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json';assert sha(dictionary)=='6a3c4311d89a9a02d55a19f2d249c960c508b5973368c06c2e5c43bb605edec9';recs=json.loads(dictionary.read_text())['programmes'];T=V/'joint_NAMPT_DDIT3_spatial_measurement_20261016_v1/private_frozen_target_and_background_definition.json';definition=json.loads(T.read_text());targets=set(definition['outside_target_weights']);background=set(definition['background_genes'])
def decode(v):
 if isinstance(v,h5py.Dataset):
  x=v[()];return np.asarray([z.decode() if isinstance(z,bytes) else str(z) for z in np.ravel(x)],dtype=str)
 if isinstance(v,h5py.Group) and 'codes' in v and 'categories' in v:
  codes=v['codes'][()].astype(int);cats=decode(v['categories']);return np.array([cats[c] if c>=0 else '' for c in codes])
 raise TypeError('Unsupported feature metadata encoding')
with h5py.File(path,'r') as f:
 ob=f['obs'];va=f['var'];index=va.attrs.get('_index','_index');index=index.decode() if isinstance(index,bytes) else index;names=decode(va[index]);nvar=len(names);oi=ob.attrs.get('_index','_index');oi=oi.decode() if isinstance(oi,bytes) else oi;nobs=len(ob[oi]);lookup={}
 for j,name in enumerate(names):lookup.setdefault(name,set()).add(j)
 for col in ['gene_symbol','gene_symbols','gene_name','gene_names','feature_name','ensembl_id','gene_id','gene_name_original','gene_id_original']:
  if col in va:
   values=decode(va[col]);assert len(values)==nvar
   for j,name in enumerate(values):
    if name not in ['','nan','None']:lookup.setdefault(name,set()).add(j)
 measured={g for g,js in lookup.items() if len(js)==1};counts=f['layers/counts'] if 'layers' in f and 'counts' in f['layers'] else None;assert counts is not None,'No explicit counts layer: review source contract before using X';encoding=counts.attrs.get('encoding-type','');encoding=encoding.decode() if isinstance(encoding,bytes) else encoding;assert encoding=='csr_matrix',encoding;shape=tuple(int(x) for x in counts.attrs['shape']);assert shape==(nobs,nvar);indptr=counts['indptr'][()].astype(np.int64);assert len(indptr)==nobs+1 and indptr[0]==0 and np.all(np.diff(indptr)>=0) and indptr[-1]==len(counts['data'])==len(counts['indices']);values=0;max_count=0;total=0
 for start in range(0,len(counts['data']),1048576):
  x=counts['data'][start:start+1048576];ix=counts['indices'][start:start+1048576];assert np.isfinite(x).all() and (x>=0).all() and np.array_equal(x,np.floor(x)) and (ix>=0).all() and (ix<nvar).all();values+=len(x);total+=float(x.sum(dtype=np.float64));max_count=max(max_count,float(x.max()) if len(x) else 0)
 obs_fields=list(ob.keys());allowed=['library_id','sample_id','sample_label','sample','slide_id','batch','orig.ident'];candidates=[];private_units={}
 for key in allowed:
  if key in ob:
   a=decode(ob[key]);assert len(a)==nobs;codes,levels=pd.factorize(a,sort=True);candidates.append(dict(field=key,unique_nonmissing=len(set(a)-{''}),missing=int((a=='').sum())))
   if len(levels)==19 and (a!='').all():private_units[key]=dict(codes=codes,levels=levels)
 assert private_units,'No19-library field found among allowed nonclinical identifiers';field=next(k for k in allowed if k in private_units);units=private_units[field]['codes'];assert 'pxl_col_in_fullres' in ob and 'pxl_row_in_fullres' in ob;coords=np.column_stack([decode(ob['pxl_col_in_fullres']).astype(float),decode(ob['pxl_row_in_fullres']).astype(float)]);assert coords.shape==(nobs,2) and np.isfinite(coords).all();np.savez_compressed(O/'private_library_rows_and_coordinates.npz',library_codes=units,spatial_coordinates=coords);spatial_uns_libraries=len(f['uns/spatial']) if 'uns' in f and 'spatial' in f['uns'] else None;matrix_metadata=dict(counts_layer='layers/counts',encoding=encoding,shape=list(shape),stored_values=values,all_stored_values_nonnegative_integers=True,maximum_single_gene_count=max_count,total_counts=total)
 genes=[]
 for r in recs:
  hit=set(r['primary_genes'])&measured;genes.append(dict(programme_id=r['programme_id'],signature_genes=len(r['primary_genes']),unambiguous_measured_genes=len(hit),measured_fraction=len(hit)/len(r['primary_genes']),original80pct10gene_measurement_gate=len(hit)>=10 and len(hit)/len(r['primary_genes'])>=.8,formal_state_or_spatial_validation=False))
 measured_lr={k:k in measured for k in ['NAMPT','ITGA5','ITGB1','DDIT3']};target_hit=targets&measured;bg_hit=background&measured;alias=sorted(measured&{'NAMPT','PBEF','PBEF1','VF','VISFATIN','ENSG00000105835'});schema=dict(feature_lookup={g:next(iter(lookup[g])) for g in sorted(measured)},library_field=field,nonclinical_library_field_choices=candidates,uns_spatial_libraries=spatial_uns_libraries,counts=matrix_metadata,coordinate_source='raw obs pixel col/row; x=column,y=row');(O/'private_native_gene_and_library_schema.json').write_text(json.dumps(schema,indent=2)+'\n')
pd.DataFrame(genes).to_csv(O/'public_all814_Visium_signature_measurement.tsv',sep='\t',index=False);pd.DataFrame([dict(dataset='GSE273952',spots=nobs,genes=nvar,libraries=19,**measured_lr,DDIT3_frozen_targets=74,DDIT3_targets_measured=len(target_hit),DDIT3_target_fraction=len(target_hit)/74,frozen_background=12657,frozen_background_measured=len(bg_hit),all_joint_gene_measurement_gates=all(measured_lr[k] for k in ['NAMPT','ITGA5','ITGB1']) and len(target_hit)/74>=.8 and len(target_hit)>=10 and len(bg_hit)>=5000,exact_original172_endpoint_identity_validated=False,spot_is_single_cell=False)]).to_csv(O/'public_Visium_joint_axis_measurement.tsv',sep='\t',index=False);out=dict(status='EXISTING_VISIUM_RAW_COUNTS_NATIVE_AXES_COORDINATES_AND_PROGRAMME_COVERAGE_CHECKED_NOT_BIOLOGICAL_VALIDATION',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),historical_technical_contract_sha256=sha(contractp),input_file_sha256=sha(path),input_file_bytes=path.stat().st_size,dataset='GSE273952',source_design='19 TLS-positive ccRCC specimens per GEO; patient/clinical bridge not reviewed here',official_design_source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273952',spots=nobs,genes=nvar,library_count=19,full_signature_measurement_gate_programmes=sum(x['original80pct10gene_measurement_gate'] for x in genes),NAMPT_identifier_matches=alias,DDIT3_targets_measured=len(target_hit),frozen_background_measured=len(bg_hit),case_level_clinical_outcome_tables_or_old_cell2location_results_used=False,old_pipeline_final_acceptance_not_inherited=True,raw_counts_read_only=True,spatial_biological_validation=False,files={p.name:sha(p) for p in O.iterdir() if p.is_file()});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='files'}))
