"""Independent raw-count sparse ranks for every projected mixed-spot score and quantity."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,time,gc
import numpy as np,pandas as pd,h5py
from scipy.sparse import csr_matrix
from scipy.stats import rankdata
V=Path(os.environ['PANCANCER_STATE_MODEL_OUTPUT_ROOT']);AR=Path(os.environ['PANCANCER_ANNOTATION_ROOT']);D=V/'GSE273952_frozen_programme_projection_20261016_v1';sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for tick in range(480):
 if (D/'inventory.json').exists():break
 e=V/'GSE273952_frozen_programme_projection_20261016_v1.execution.json'
 if e.exists() and json.loads(e.read_text()).get('exit_code') is not None:raise RuntimeError('Projection terminal without inventory')
 time.sleep(30)
else:raise TimeoutError('Four-hour score dependency wait')
di=json.loads((D/'inventory.json').read_text());I=V/'existing_GSE273952_Visium_raw_preflight_20261016_v2';lookup=json.loads((I/'private_native_gene_and_library_schema.json').read_text())['feature_lookup'];old=AR.parent.parent/'Project_v2/Downstream/0716_Data/immunotherapy_4cohort_analysis/20260803_v1';path=Path(json.loads((old/'outputs/spatial/GSE273952_spatial_contract_audit.json').read_text())['input_h5ad']);assert sha(path)==di['raw_input_sha256'];recs=json.loads((V/'all_lineage_state_candidate_dictionary_20261014_v4/frozen_programme_dictionary.json').read_text())['programmes'];frozen=json.loads((V/'joint_NAMPT_DDIT3_spatial_measurement_20261016_v1/private_frozen_target_and_background_definition.json').read_text())['outside_target_weights'];pos=[g for g,w in frozen.items() if w==1];neg=[g for g,w in frozen.items() if w==-1];ids=[r['programme_id'] for r in recs]+['RNA_DDIT3_TARGETS_POS','RNA_DDIT3_TARGETS_NEG','RNA_DDIT3_TARGETS_SIGNED'];indices=[list(dict.fromkeys(lookup[g] for g in r['primary_genes'] if g in lookup)) for r in recs]+[[lookup[g] for g in pos],[lookup[g] for g in neg],[lookup[g] for g in sorted(frozen)]];counts=np.array([len(x) for x in indices]);gate=np.array([len(j)>=10 and len(j)/len(r['primary_genes'])>=.8 for r,j in zip(recs,indices[:814])]+[True,True,True]);assert len(ids)==817;rows=[];cols=[]
for j,ii in enumerate(indices):
 for k in ii:rows.append(k);cols.append(j)
O=V/'GSE273952_frozen_programme_independent_audit_20261016_v1';assert not O.exists();O.mkdir(mode=0o700);seen=np.zeros(71398,bool);manual_values=0;maxerr=0;quantity=0;checks=[];allscore=np.full((71398,817),np.nan,np.float32)
with h5py.File(path,'r') as f:
 raw=f['layers/counts'];m=int(raw.attrs['shape'][1]);inc=csr_matrix((np.ones(len(rows)),(rows,cols)),shape=(m,817));base=counts*(counts+1)/2;den=counts*1500-base
 for info in di['private_outputs']:
  p=D/info['relative_file'];assert sha(p)==info['sha256'];z=np.load(p,allow_pickle=False);rr=z['spot_rows'];assert np.array_equal(z['programme_ids'],np.array(ids)) and np.all(np.diff(rr)==1) and not seen[rr].any();lo,hi=int(rr[0]),int(rr[-1])+1;ip=raw['indptr'][lo:hi+1].astype(np.int64);x=csr_matrix((raw['data'][int(ip[0]):int(ip[-1])],raw['indices'][int(ip[0]):int(ip[-1])],ip-ip[0]),shape=(len(rr),m));x.sum_duplicates();x.sort_indices();d=np.zeros(x.nnz,float)
  for j in range(len(rr)):
   a,b=x.indptr[j:j+2];v=x.data[a:b];nz=np.flatnonzero(v>0);d[a+nz]=np.minimum(1500,rankdata(-v[nz].astype(float),method='average'))-1500
  rank_delta=csr_matrix((d,x.indices.copy(),x.indptr.copy()),shape=x.shape);delta=(rank_delta@inc).toarray();num=counts*1500+delta-base;expected=1-np.divide(num,den,out=np.full(num.shape,np.nan),where=den>0);expected[:,~gate]=np.nan
  for j in [814,815]:
   if counts[j]==0:expected[:,j]=0
  expected[:,816]=np.maximum(expected[:,814]-expected[:,815],0);assert np.allclose(expected,z['scores'],atol=1e-6,equal_nan=True);error=float(np.nanmax(abs(expected-z['scores'])));maxerr=max(maxerr,error);manual_values+=int(np.isfinite(expected).sum());total=np.asarray(x.sum(1)).ravel();det=((x>0).astype(int)@inc).toarray();fraction=(x@inc).toarray()/total[:,None];assert np.array_equal(det,z['detected_signature_genes']) and np.allclose(fraction,z['signature_RNA_fraction'],atol=1e-7);assert np.array_equal(total,z['RNA_total']) and np.array_equal(np.asarray((x>0).sum(1)).ravel(),z['detected_genes']);lr=np.column_stack([np.asarray(x[:,lookup[g]].toarray()).ravel() for g in ['NAMPT','ITGA5','ITGB1','DDIT3']]);assert np.array_equal(lr,z['raw_NAMPT_ITGA5_ITGB1_DDIT3_counts']);quantity+=det.size*2+lr.size+len(rr)*2;allscore[rr]=z['scores'];seen[rr]=True;checks.append(dict(chunk=info['relative_file'],spots=len(rr),manual_values=int(np.isfinite(expected).sum()),max_error=error));(O/'progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),spots_checked=int(seen.sum()),manual_score_values=manual_values,max_error=maxerr))+'\n');del x,rank_delta,expected;gc.collect()
assert seen.all();s=pd.read_csv(D/'public_all817_cohort_RNA_score_summary.tsv',sep='\t').set_index('programme_id');assert set(s.index)==set(ids)
for j,key in enumerate(ids):
 v=allscore[:,j];v=v[np.isfinite(v)];r=s.loc[key];assert len(v)==r.valid_spots
 for field,val in [('mean_score',v.astype(float).mean() if len(v) else np.nan),('q10',np.quantile(v,.1) if len(v) else np.nan),('median_score',np.median(v) if len(v) else np.nan),('q90',np.quantile(v,.9) if len(v) else np.nan)]:assert np.isclose(val,r[field],atol=1e-7,equal_nan=True)
result=dict(status='PASS_ALL71398_MIXED_SPOT_FROZEN_AND_SIGNED_UCELL_RAW_COUNTS_SUPPORT_FRACTIONS_AND_COHORT_SUMMARIES',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),projection_inventory_sha256=sha(D/'inventory.json'),spots_checked=71398,columns=817,manual_score_values=manual_values,max_manual_error=maxerr,quantity_values_checked=quantity,summary_rows_checked=817,clinical_outcomes_used=False,cell_identity_or_protein_activity_or_biological_validation=False);(O/'inventory.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
