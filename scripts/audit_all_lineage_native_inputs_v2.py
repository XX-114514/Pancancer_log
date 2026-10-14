"""Fresh canonical COO remapping readback, independent sparse matrix assembly."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd,anndata as ad,zarr
from scipy.sparse import csr_matrix,coo_matrix
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);I=R/'all_lineage_state_singlecell_inputs_20261014_v2';iv=json.loads((I/'inventory.json').read_text());assert iv['status']=='EXTRACTED_ALL100766_COUNTS_PENDING_NATIVE_READBACK_AUDIT'
D=R/'all_lineage_state_discovery_row_selection_20261014_v1';metadata=pd.read_csv(D/'private_discovery_metadata.tsv.gz',sep='\t',keep_default_na=False);gg=json.loads((D/'fixed_common_training_genes_by_major.json').read_text());U=AR/'20261010_uniform_expression_annotation_v1';contract=json.loads((U/'input_contract.json').read_text());g=zarr.open_group(contract['core_zarr'],mode='r');x=g['X'];checked=0;native_values=0;common_values=0;checks=[]
for item in iv['files']:
 p=Path(item['path']);assert sha(p)==item['sha256'];a=ad.read_h5ad(p);major=item['major'];expected=metadata[metadata.major.eq(major)].reset_index(drop=True)
 assert a.n_obs==len(expected) and list(a.var_names)==gg[major]
 for col in expected.columns:
  left=a.obs[col].astype(str).to_numpy();right=expected[col].astype(str).to_numpy();assert np.array_equal(left,right),col
 assert np.array_equal(a.obs.analysis_block.astype(str),expected.effective_resampling_block.astype(str))
 assert np.array_equal(a.obs.patient_identity_reviewed,expected.unit_kind.eq('reviewed_patient'))
 rows=expected.core_position.to_numpy();names=g['obs'][g['obs'].attrs['_index']].get_orthogonal_selection((rows,)).astype(str);assert np.array_equal(a.obs_names,names)
 maps=[]
 for sid in sorted(expected.gse_id.unique()):
  tab=pd.read_csv(U/'by_gse'/sid/'genes.tsv',sep='\t').set_index('gene');maps.append(tab.loc[a.var_names,'core_column'].to_numpy())
 assert all(np.array_equal(maps[0],v) for v in maps)
 colmap=np.full(contract['raw_shape'][1],-1,np.int32);colmap[maps[0]]=np.arange(a.n_vars)
 first=0
 while first<len(rows):
  last=min(first+128,int(np.searchsorted(rows,rows[first]+4096)),len(rows));assert last>first;rr=rows[first:last];lo=int(rr[0]);hi=int(rr[-1])+1;ptr=x['indptr'][lo:hi+1].astype(np.int64);start,end=int(ptr[0]),int(ptr[-1]);full=csr_matrix((x['data'][start:end],x['indices'][start:end],ptr-start),shape=(hi-lo,contract['raw_shape'][1]))[rr-lo]
  co=full.tocoo();native_values+=len(co.data);total=np.bincount(co.row,weights=co.data,minlength=len(rr));detected=np.bincount(co.row,weights=co.data>0,minlength=len(rr))
  assert np.array_equal(total,a.obs.full_RNA_total.to_numpy()[first:last]) and np.array_equal(detected,a.obs.full_detected_genes.to_numpy()[first:last])
  cc=colmap[co.col];keep=cc>=0;fresh=coo_matrix((co.data[keep],(co.row[keep],cc[keep])),shape=(len(rr),a.n_vars)).tocsr();difference=fresh-a.X[first:last];difference.eliminate_zeros();assert difference.nnz==0
  common_values+=int(keep.sum());checked+=len(rr);first=last
 checks.append(dict(major=major,cells=a.n_obs,genes=a.n_vars,sources=expected.gse_id.nunique(),exact_raw_and_native_denominator=True));del a;gc.collect()
assert checked==100766 and len(checks)==18
out=dict(status='PASS_ALL100766_NATIVE_RAW_COUNTS_FULLRNA_DENOMINATORS_IDENTITIES_AND_GENE_AXES',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_inventory_sha256=sha(I/'inventory.json'),audit_code_sha256=sha(__file__),cells=checked,majors=18,native_sparse_values_read=native_values,common_sparse_values_compared=common_values,all_cells_independently_reassembled_via_COO=True,formal_states=False,models_not_trained_by_this_audit=True,checks=checks)
O=R/'all_lineage_state_native_input_audit_20261014_v2';assert not O.exists();O.mkdir();(O/'independent_audit.json').write_text(json.dumps(out,indent=2)+'\n')
accepted=dict(status='ACCEPTED_ALL18_NATIVE_DISCOVERY_INPUTS_NOT_FORMAL_STATES',observed_UTC=out['observed_UTC'],producer_inventory_sha256=sha(I/'inventory.json'),audit_sha256=sha(O/'independent_audit.json'),files=iv['files'],cells=100766,majors=18);(O/'accepted_input_inventory.json').write_text(json.dumps(accepted,indent=2)+'\n');print(json.dumps(out),flush=True)
