"""Uniform 20,000 whole-block residual permutations; frozen models/family/CI.
Extends the exact first 2,000 permutations for EVERY estimated model.
Exploratory precision analysis, no thresholds or biological labels changed.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
from statsmodels.stats.multitest import multipletests
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';A=P/'CM_LR_full_control_association_v1';L=P/'CM_edge_LR_extraction_v2';O=P/'CM_LR_permutation_precision_v2';planpath=Path(os.environ['PANCANCER_PRECISION_PLAN']);plan=json.loads(planpath.read_text());assert sha(__file__)==plan['producer_sha256'] and not O.exists();assert sha(A/'inventory.json')==plan['association_inventory_sha256'] and sha(A/'independent_audit.json')==plan['association_audit_sha256'];inv=json.loads((A/'inventory.json').read_text());f=A/'all26994_full_control_CM_LR_associations.tsv.gz';assert sha(f)==inv['outputs'][f.name];d=pd.read_csv(f,sep='\t');O.mkdir();start=datetime.now(timezone.utc).isoformat();results={};checked=0;estimated=d[d.status=='ESTIMABLE_FULL_CONTROL_SPAN_ASSOCIATION'];requested=plan['permutations'];examples=[]
 for gid,g in estimated.groupby('group_id',sort=True):
  b=dict(np.load(L/'private_group_arrays'/(gid+'.npz'),allow_pickle=False));lr={(l,r):i for i,(l,r) in enumerate(zip(b['ligands'],b['receptors']))};cms={c:i for i,c in enumerate(b['CMs'])};groups={}
  for r in g.itertuples():
   j=lr[r.ligand,r.receptor];mask=np.isfinite(b['scores'][:,j,0]);key=(r.variant,r.cm,np.packbits(mask).tobytes());groups.setdefault(key,[]).append((r,j,mask))
  for (variant,cm,_),entries in groups.items():
   mask=entries[0][2];rows=[v[0] for v in entries];columns=[v[1] for v in entries];blocks=b['patient_blocks'][mask];x=b['endpoint_excluded_activity' if variant=='endpoint_excluded_library' else 'original_H'][mask,cms[cm]];raw=np.column_stack([x,b['source_within_major'][mask],b['target_within_major'][mask],np.log1p(b['source_n'][mask]),np.log1p(b['target_n'][mask])]);sd=raw.std(axis=0,ddof=1);z=np.divide(raw-raw.mean(axis=0),sd,out=np.zeros_like(raw),where=sd>1e-12);Z=np.column_stack([np.ones(mask.sum()),z[:,1:]]);u,s,_=np.linalg.svd(Z,full_matrices=False);rank=(s>np.finfo(float).eps*max(Z.shape)*s[0]).sum();Q=u[:,:rank];v=z[:,0]-Q@(Q.T@z[:,0]);a=v/(v@v);Y=np.log1p(b['scores'][mask][:,columns,0]);Y=(Y-Y.mean(axis=0))/Y.std(axis=0,ddof=1);e=Y-Q@(Q.T@Y);beta=a@Y;assert np.allclose(beta,[r.beta for r in rows],atol=1e-8,rtol=1e-6);seeds={int(r.resampling_seed) for r in rows};assert len(seeds)==1;rng=np.random.default_rng(seeds.pop());order=list(dict.fromkeys(blocks));ids=[np.flatnonzero(blocks==k) for k in order];bins={}
   for i,index in enumerate(ids):bins.setdefault(len(index),[]).append(i)
   hits=np.zeros(len(rows),dtype=np.int64);hits2k=None
   for begin in range(0,requested,1000):
    size=min(1000,requested-begin);pl=np.empty((size,len(blocks)))
    for draw in range(size):
     p=np.arange(len(blocks))
     for members in bins.values():
      shuffle=rng.permutation(members)
      for target,source in zip(members,shuffle):p[ids[target]]=ids[source]
     pl[draw,p]=a
    hits+=(abs(pl@e)>=abs(beta)[None,:]-1e-12).sum(axis=0)
    if begin+size==2000:hits2k=hits.copy()
   assert np.allclose((hits2k+1)/2001,[r.residual_block_permutation_p for r in rows],atol=1e-12);checked+=len(rows)
   for i,r in enumerate(rows):results[r.Index]=dict(permutation_exceedances_20000=int(hits[i]),precision_block_permutation_p=float((hits[i]+1)/(requested+1)),precision_permutations=requested)
  print('completed_group',gid,'models',checked,flush=True)
 assert checked==len(estimated)==19127
 for i,r in results.items():
  for k,v in r.items():d.loc[i,k]=v
 for variant,g in d.groupby('variant'):
  p=g.precision_block_permutation_p.fillna(1).to_numpy();q=multipletests(p,method='fdr_bh')[1];d.loc[g.index,'precision_full_family_BH_q']=np.where(g.precision_block_permutation_p.notna(),q,np.nan)
 file=O/'all26994_precision20000_associations.tsv.gz';d.to_csv(file,sep='\t',index=False,na_rep='NA');summary={v:dict(models=len(g),estimated=int(g.precision_block_permutation_p.notna().sum()),q05=int((g.precision_full_family_BH_q<.05).sum()),q05_stable_positive_CI=int(((g.precision_full_family_BH_q<.05)&(g.bootstrap_ci_low>0)).sum()),q05_stable_negative_CI=int(((g.precision_full_family_BH_q<.05)&(g.bootstrap_ci_high<0)).sum())) for v,g in d.groupby('variant')};out=dict(status='COMPLETE_PRECISION_EXTENSION_PENDING_INDEPENDENT_AUDIT',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planpath),all_first2000_permutation_p_reproduced=checked,frozen_model_coefficients_and_bootstrap=True,summary=summary,scientific_validation=False,outputs={file.name:sha(file)});(O/'inventory.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)
if __name__=='__main__':main()
