"""All input IDs/quantities plus fresh canonical counts and UCell for three cells/file."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd,anndata as ad,zarr
from scipy.sparse import csr_matrix,coo_matrix,vstack
import pyucell as uc

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);I=V/'state_spatial_proxy_native_cells_20261016_v1';iv=json.loads((I/'inventory.json').read_text());J=V/'state_spatial_proxy_calibration_inputs_20261016_v1';sel=np.load(J/'private_fixed_calibration_rows.npz',allow_pickle=False);expected=dict(zip(sel['global_rows'],sel['unit_index']));U=AR/'20261010_uniform_expression_annotation_v1';contract=json.loads((U/'input_contract.json').read_text());g=zarr.open_group(contract['core_zarr'],mode='r');x=g['X'];F=R/'all_CM_state_native_input_preflight_20261014_v1';units=pd.read_csv(F/'private_all_frozen_sample_node_units.tsv.gz',sep='\t',keep_default_na=False);dictionary=json.loads((V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json').read_text());d={r['programme_id']:r for r in dictionary['programmes']};O=V/'state_spatial_native_calibration_input_audit_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);seen=set();nc=0;nr=0;nscore=0;maxerr=0;checks=[]
for info in iv['private_native_files']:
 p=I/info['relative_file'];assert sha(p)==info['sha256'];a=ad.read_h5ad(p);rows=a.obs.core_position.to_numpy(int);uu=a.obs.unit_index.to_numpy(int);assert len(set(rows))==len(rows) and not seen.intersection(rows);seen.update(rows.tolist());assert all(expected[int(r)]==u for r,u in zip(rows,uu));names=np.asarray(g['obs'][g['obs'].attrs['_index']].get_orthogonal_selection((rows,)).tolist(),dtype=str);assert np.array_equal(names,a.obs_names.to_numpy(str));meta=units.iloc[uu]
 for key in ['gse_id','major','node_id','stratum','effective_resampling_block','unit_kind']:assert np.array_equal(meta[key].to_numpy(str),a.obs[key].to_numpy(str))
 assert np.array_equal(np.asarray(a.X.sum(1)).ravel(),a.obs.full_RNA_total) and np.array_equal(np.asarray((a.X>0).sum(1)).ravel(),a.obs.full_detected_genes);gt=pd.read_csv(U/'by_gse'/info['gse_id']/'genes.tsv',sep='\t');assert list(a.var_names)==list(gt.gene);mapping=dict(zip(gt.core_column.to_numpy(int),np.arange(len(gt))));pick=np.unique(np.linspace(0,a.n_obs-1,min(3,a.n_obs),dtype=int));fresh=[]
 for i in pick:
  row=int(rows[i]);ip=x['indptr'][row:row+2].astype(np.int64);data=x['data'][int(ip[0]):int(ip[1])];cols=x['indices'][int(ip[0]):int(ip[1])];converted=np.array([mapping[int(j)] for j in cols],int);v=coo_matrix((data,(np.zeros(len(data),int),converted)),shape=(1,a.n_vars)).tocsr();diff=v-a.X[i];diff.eliminate_zeros();assert diff.nnz==0;fresh.append(v);nr+=1
 aa=ad.AnnData(vstack(fresh,format='csr'),var=pd.DataFrame(index=a.var_names.copy()));pids=list(a.uns['programme_ids']);sig={}
 for pid in pids:
  genes=d[pid]['primary_genes'];present=[z for z in genes if z in a.var_names]
  if len(present)>=10 and len(present)/len(genes)>=.8:sig[pid]=present
 if sig:uc.compute_ucell_scores(aa,signatures=sig,max_rank=1500,ties_method='average',missing_genes='skip',chunk_size=3,n_jobs=1,device='cpu',suffix='_audit')
 err=0
 for j,pid in enumerate(pids):
  old=a.obsm['accepted_full_native_UCell1500'][pick,j]
  if pid not in sig:assert np.isnan(old).all();continue
  val=aa.obs[pid+'_audit'].to_numpy();e=float(np.max(abs(val-old)));assert e<1e-6;err=max(err,e);nscore+=len(pick)
 checks.append(dict(gse_id=info['gse_id'],major=info['major'],all_cells_identity_quantity_checked=a.n_obs,fresh_canonical_count_cells=len(pick),full_native_UCell_score_values=len(sig)*len(pick),max_full_native_reference_error=err));maxerr=max(maxerr,err);nc+=a.n_obs;(O/'progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),all_ID_quantity_cells=nc,fresh_canonical_count_cells=nr,full_reference_scores=nscore))+'\n');del a,aa,fresh;gc.collect()
assert seen==set(expected) and nc==iv['cells']==91447;pd.DataFrame(checks).to_csv(O/'public_all709_native_input_readback_checks.tsv',sep='\t',index=False);out=dict(status='PASS_ALL91447_ID_NATIVE_AXES_QUANTITIES_AND_THREE_REAL_CELL_CANONICAL_COUNTS_UCELL_PER709_FILE',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),native_input_inventory_sha256=sha(I/'inventory.json'),all_ID_quantity_cells=nc,fresh_canonical_count_cells=nr,full_native_reference_values_checked=nscore,max_full_reference_score_error=maxerr,all_selected_gene_values_freshly_reassembled=False,public_output_sha256=sha(O/'public_all709_native_input_readback_checks.tsv'),spatial_calibration_or_biological_validation=False);(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
