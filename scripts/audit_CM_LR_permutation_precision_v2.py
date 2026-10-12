"""Separate full-X least-squares, explicit shuffled response audit."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';O=P/'CM_LR_permutation_precision_v2';A=P/'CM_LR_full_control_association_v1';L=P/'CM_edge_LR_extraction_v2';planpath=Path(os.environ['PANCANCER_PRECISION_PLAN']);plan=json.loads(planpath.read_text());assert sha(__file__)==plan['auditor_sha256'];inv=json.loads((O/'inventory.json').read_text());assert inv['plan_sha256']==sha(planpath);assert not (O/'independent_audit.json').exists()
 for n,h in inv['outputs'].items():assert sha(O/n)==h
 d=pd.read_csv(O/'all26994_precision20000_associations.tsv.gz',sep='\t');old=pd.read_csv(A/'all26994_full_control_CM_LR_associations.tsv.gz',sep='\t');assert len(d)==len(old)==26994
 for c in old:
  if pd.api.types.is_numeric_dtype(old[c]):assert np.allclose(d[c],old[c],equal_nan=True,atol=1e-12,rtol=1e-12),c
  else:assert d[c].fillna('NA').equals(old[c].fillna('NA')),c
 estimated=d[d.status=='ESTIMABLE_FULL_CONTROL_SPAN_ASSOCIATION'];assert estimated.precision_permutations.eq(20000).all();assert np.allclose((estimated.permutation_exceedances_20000+1)/20001,estimated.precision_block_permutation_p)
 for variant,g in d.groupby('variant'):
  p=g.precision_block_permutation_p;assert np.allclose(g.precision_full_family_BH_q[p.notna()],multipletests(p.fillna(1),method='fdr_bh')[1][p.notna()])
 # Preserve producer's first six real fixtures; add first positive/negative for both variants.
 keys=[]
 for file in sorted((A/'private_resampling_examples').glob('*.npz')):
  b=dict(np.load(file,allow_pickle=False));keys.append((str(b['case_id']),str(b['cm']),str(b['variant'])))
 for variant,g in estimated.groupby('variant'):
  for sign in [-1,1]:
   row=g[g.beta*sign>0].iloc[0];keys.append((row.case_id,row.cm,row.variant))
 examples=[]
 for key in dict.fromkeys(keys):
  row=d[(d.case_id==key[0])&(d.cm==key[1])&(d.variant==key[2])].iloc[0];b=dict(np.load(L/'private_group_arrays'/(row.group_id+'.npz'),allow_pickle=False));j=list(zip(b['ligands'],b['receptors'])).index((row.ligand,row.receptor));ci=list(b['CMs']).index(row.cm);mask=np.isfinite(b['scores'][:,j,0]);blocks=b['patient_blocks'][mask];raw=np.column_stack([b['endpoint_excluded_activity' if row.variant=='endpoint_excluded_library' else 'original_H'][mask,ci],b['source_within_major'][mask],b['target_within_major'][mask],np.log1p(b['source_n'][mask]),np.log1p(b['target_n'][mask])]);sd=raw.std(axis=0,ddof=1);X=np.column_stack([np.ones(mask.sum()),np.divide(raw-raw.mean(axis=0),sd,out=np.zeros_like(raw),where=sd>1e-12)]);Y=np.log1p(b['scores'][mask,j,0]);Y=(Y-Y.mean())/Y.std(ddof=1);Z=X[:,[0,2,3,4,5]];fit=Z@np.linalg.lstsq(Z,Y,rcond=None)[0];err=Y-fit;observed=np.linalg.lstsq(X,Y,rcond=None)[0][1];assert np.isclose(observed,row.beta);order=list(dict.fromkeys(blocks));ids=[np.flatnonzero(blocks==k) for k in order];bins={}
  for i,v in enumerate(ids):bins.setdefault(len(v),[]).append(i)
  rng=np.random.default_rng(int(row.resampling_seed));hits=0;first=0
  for draw in range(20000):
   vector=err.copy()
   for members in bins.values():
    shuffle=rng.permutation(members)
    for target,source in zip(members,shuffle):vector[ids[target]]=err[ids[source]]
   coefficient=np.linalg.lstsq(X,fit+vector,rcond=None)[0][1];hits+=abs(coefficient)>=abs(observed)-1e-12
   if draw==1999:first=hits
  assert hits==int(row.permutation_exceedances_20000) and np.isclose((first+1)/2001,row.residual_block_permutation_p);examples.append(dict(case_id=row.case_id,cm=row.cm,variant=row.variant,beta_sign='positive' if observed>0 else 'negative',hits=int(hits),status='PASS_EXPLICIT20000_FULLX'))
 report=dict(status='PASS_ALL26994_UNCHANGED_MODELS_AND_BH_WITH_EXPLICIT20000_BOTH_SIGN_EXAMPLES',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),n_all_rows=26994,n_extended_models=19127,all_original_columns_unchanged=True,explicit_examples=examples,scientific_validation=False);(O/'independent_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
