"""All selection/seed/count/BH checks and independent explicit bootstrap refits.
Algorithm fixtures verify numerical bootstrap transform, not all 20000 draws replayed.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib
import numpy as np,pandas as pd
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bh(p):
 i=np.argsort(p,kind='stable');q=np.array(p)[i]*len(p)/np.arange(1,len(p)+1);out=np.empty(len(p));out[i]=np.minimum(1,np.minimum.accumulate(q[::-1])[::-1]);return out

def main():
 R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';S=P/'major_LR_receiver_TF_screen_v1';O=P/'major_LR_receiver_TF_wild_v1';plan=json.loads(Path(os.environ['PANCANCER_LR_TF_WILD_PLAN']).read_text());assert sha(__file__)==plan['auditor_sha256'] and not (O/'independent_audit.json').exists();inv=json.loads((O/'inventory.json').read_text())
 for n,h in inv['outputs'].items():assert sha(O/n)==h
 original=pd.read_csv(S/'all177393_LR_union_models.tsv.gz',sep='\t');table=pd.read_csv(O/'all177393_LR_wild_models.tsv.gz',sep='\t');records=pd.read_csv(O/'all6069_selected_recheck_records.tsv.gz',sep='\t',dtype={'seed':str});selected=original[(original.variant=='abundance_complexity')&(original.joint_p_exploratory<.05)];keys=set(zip(selected.context_id,selected.LR_index));assert len(keys)==2023 and len(records)==6069 and len(table)==177393;assert not records.duplicated(['context_id','LR_index','variant']).any();assert set(zip(records.context_id,records.LR_index))==keys;assert list(zip(table.context_id,table.LR_index,table.variant))==list(zip(original.context_id,original.LR_index,original.variant));lookup=table.set_index(['context_id','LR_index','variant']);rechecked=0
 for x in records.itertuples(index=False):
  r=lookup.loc[(x.context_id,x.LR_index,x.variant)];assert r.wild_status==x.state
  if x.state=='wild_rechecked_exploratory':
   seed=int(hashlib.sha256((x.context_id+'|'+str(x.LR_index)+'|'+x.variant+'|wildmax20261012').encode()).hexdigest()[:16],16);assert int(x.seed)==seed;B=int(x.draws);ext=int(x.extreme);assert 0<=ext<=B;enumerated=str(x.enumerated_all_sign_pairs).lower()=='true';p=ext/B if enumerated else (ext+1)/(B+1);assert np.isclose(r.wild_joint_p,p,atol=1e-14);assert B==r.wild_draws and ext==r.wild_extreme
   if not enumerated:assert B==plan['bootstrap_draws']
   rechecked+=1
  else:assert pd.isna(r.wild_joint_p)
 not_selected=~table.apply(lambda r:(r.context_id,r.LR_index) in keys,axis=1);assert table.loc[not_selected].wild_status.eq('not_selected_exploratory_screen').all() and table.loc[not_selected].wild_joint_p.isna().all()
 for v,g in table.groupby('variant'):
  assert len(g)==59131;p=g.wild_joint_p.fillna(1).to_numpy();q=bh(p);valid=g.wild_status.eq('wild_rechecked_exploratory').to_numpy();assert np.allclose(g.wild_BH_exploratory.to_numpy()[valid],q[valid],atol=1e-12);by=np.minimum(1,q*sum(1/i for i in range(1,len(g)+1)));assert np.allclose(g.wild_BY_exploratory.to_numpy()[valid],by[valid],atol=1e-12)
 fixtures=list((O/'private_weight_audit_cases').glob('*.npz'));checked=0;maximum=0
 for f in fixtures:
  a=np.load(f);Z=a['Z'];x=a['x'];Y=a['Y'];blocks=a['blocks'];U,s,_=np.linalg.svd(Z,full_matrices=False);Q=U[:,s>s[0]*1e-10];Hz=Q@Q.T;e0=Y-Q@(Q.T@Y);groups=[np.flatnonzero(blocks==g) for g in np.unique(blocks)];adj=np.zeros_like(Y)
  for ix in groups:
   vals,V=np.linalg.eigh(np.eye(len(ix))-Hz[np.ix_(ix,ix)]);A=(V/np.sqrt(vals))@V.T;adj[ix]=A@e0[ix]
  X=np.column_stack([Z,x]);pinv=np.linalg.pinv(X,rcond=1e-10);H=X@pinv;values=[]
  for w in a['weights']:
   ys=adj.copy()
   for j,ix in enumerate(groups):ys[ix]*=w[j]
   beta=(pinv@ys)[-1];res=ys-X@(pinv@ys);scores=[]
   for ix in groups:
    vals,V=np.linalg.eigh(np.eye(len(ix))-H[np.ix_(ix,ix)]);A=(V/np.sqrt(vals))@V.T;scores.append(pinv[-1,ix]@A@res[ix])
   se=np.sqrt(np.sum(np.array(scores)**2,axis=0));ts=np.divide(abs(beta),se,out=np.zeros_like(beta),where=se>1e-14);ts[(se<=1e-14)&(abs(beta)>1e-12)]=np.inf;values.append(ts.max())
  assert np.allclose(values,a['fast_max_statistics'],rtol=2e-6,atol=2e-6,equal_nan=True);maximum=max(maximum,float(np.max(abs(np.array(values)-a['fast_max_statistics']))));checked+=len(values)
 result=dict(status='PASS_ALL_SELECTION_SEEDS_PROBABILITY_COUNTS_FAMILY_BH_AND_EXPLICIT_BOOTSTRAP_FIXTURES',observed_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),selected_LR_cases=len(keys),all_recheck_records_checked=len(records),wild_rechecked_records=rechecked,all_original_family_rows_checked=len(table),independent_fullX_bootstrap_refits=checked,fixtures=len(fixtures),max_absolute_refit_difference=maximum,all20000_MC_draws_not_independently_replayed=True,same_data_selection_not_independent_validation=True,no_causal_or_binding_claim=True);(O/'independent_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if __name__=='__main__':main()
