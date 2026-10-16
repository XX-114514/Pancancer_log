"""Independent raw spot/count/coordinate and gene coverage audit, no clinical values."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd,h5py
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);I=V/'existing_GSE273952_Visium_raw_preflight_20261016_v2';iv=json.loads((I/'inventory.json').read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();old=AR.parent.parent/'Project_v2/Downstream/0716_Data/immunotherapy_4cohort_analysis/20260803_v1';cp=old/'outputs/spatial/GSE273952_spatial_contract_audit.json';p=Path(json.loads(cp.read_text())['input_h5ad']);assert sha(p)==iv['input_file_sha256'];sp=I/'private_native_gene_and_library_schema.json';assert sha(sp)==iv['files'][sp.name];schema=json.loads(sp.read_text());gp=I/'private_library_rows_and_coordinates.npz';assert sha(gp)==iv['files'][gp.name];z=np.load(gp,allow_pickle=False);counts_rows=[];independent_lookup={};raw_nnz=0

def de(g):
 if isinstance(g,h5py.Group):
  cats=de(g['categories']);codes=g['codes'][()];return np.asarray([cats[k] if k>=0 else '' for k in codes])
 a=g[()];return np.asarray([x.decode() if isinstance(x,bytes) else str(x) for x in a])
with h5py.File(p,'r') as f:
 X=f['layers/counts'];n,m=map(int,X.attrs['shape']);assert n==71398 and m==24146;obs=f['obs'];labels=de(obs['sample_label']);categories=sorted(set(labels));assert len(categories)==19;codes=np.array([categories.index(k) for k in labels]);assert np.array_equal(codes,z['library_codes']);x=de(obs['pxl_col_in_fullres']).astype(float);y=de(obs['pxl_row_in_fullres']).astype(float);coords=np.stack([x,y],axis=1);assert np.array_equal(coords,z['spatial_coordinates']) and np.isfinite(coords).all();ar=de(obs['array_row']).astype(float);ac=de(obs['array_col']).astype(float);assert np.array_equal(ar,np.floor(ar)) and np.array_equal(ac,np.floor(ac));duplicates=0
 for c in range(19):
  ii=np.flatnonzero(codes==c);duplicates+=len(ii)-len(np.unique(coords[ii],axis=0));assert len(np.unique(np.stack([ar[ii],ac[ii]],axis=1),axis=0))==len(ii)
 assert duplicates==0
 for begin in range(0,n,257):
  end=min(n,begin+257);ip=X['indptr'][begin:end+1].astype(np.int64);data=X['data'][int(ip[0]):int(ip[-1])];ix=X['indices'][int(ip[0]):int(ip[-1])];assert np.isfinite(data).all() and (data>=0).all() and np.array_equal(data,np.floor(data)) and (ix>=0).all() and (ix<m).all();raw_nnz+=len(data);rid=np.repeat(np.arange(end-begin),np.diff(ip));sums=np.bincount(rid,weights=data,minlength=end-begin);det=np.bincount(rid,weights=data>0,minlength=end-begin);counts_rows.append(np.column_stack([sums,det]))
 q=np.vstack(counts_rows);assert len(q)==n and (q[:,0]>0).all();var=f['var'];key=var.attrs['_index'];key=key.decode() if isinstance(key,bytes) else key
 for col in [key,'gene_symbol','gene_symbols','gene_name','gene_names','feature_name','ensembl_id','gene_id','gene_name_original','gene_id_original']:
  if col in var:
   for j,name in enumerate(de(var[col])):
    if name not in ['','nan','None']:independent_lookup.setdefault(name,set()).add(j)
 lookup={g:next(iter(js)) for g,js in independent_lookup.items() if len(js)==1};assert lookup==schema['feature_lookup'];dictionary=json.loads((V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json').read_text())['programmes'];T=pd.read_csv(I/'public_all814_Visium_signature_measurement.tsv',sep='\t').set_index('programme_id');checked=0
 for r in dictionary:
  present=set(r['primary_genes'])&set(lookup);assert len(present)==T.loc[r['programme_id'],'unambiguous_measured_genes'];checked+=1
 assert checked==814
O=V/'existing_GSE273952_Visium_input_audit_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);np.savez_compressed(O/'private_all_spot_derived_QC_and_grid.npz',RNA_total=q[:,0],detected_genes=q[:,1],array_row=ar,array_col=ac,library_codes=codes);out=dict(status='PASS_ALL71398_RAW_COUNT_ROWS_NATIVE_GENE_LOOKUP19_LIBRARIES_UNIQUE_COORDINATES814_GENE_COVERAGES',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(I/'inventory.json'),input_file_sha256=iv['input_file_sha256'],spots=71398,genes=24146,specimen_library_units=19,reviewed_patient_identity_bridge=False,raw_sparse_values_checked=raw_nnz,duplicate_within_library_coordinates=0,empty_count_spots=0,signature_measurement_rows_checked=814,all_count_rows_streamed_independently=True,clinical_outcome_values_used=False,old_deconvolution_values_used=False,spatial_biological_validation=False,private_QC_sha256=sha(O/'private_all_spot_derived_QC_and_grid.npz'));(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
