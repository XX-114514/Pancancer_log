"""All fixed context x MCP slots; internal holdout only, no selection by discovery p."""
from pathlib import Path
import sys,json,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.stats import rankdata
P=Path(sys.argv[1]);M=Path(sys.argv[2]);O=Path(sys.argv[3]);assert not O.exists();O.mkdir(mode=0o700)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
scope=pd.read_csv(P/'public_slot_scope.tsv',sep='\t');status=pd.read_csv(M/'public_fit_status.tsv',sep='\t');path=P/'private_block_mean_scores.tsv';s=pd.read_csv(path,sep='\t') if path.exists() else pd.DataFrame();rows=[];rng=np.random.default_rng(2026101603)
def rho(x,y):
 a=rankdata(x);b=rankdata(y);a-=a.mean();b-=b.mean();d=np.linalg.norm(a)*np.linalg.norm(b);return float(a@b/d) if d>0 else np.nan
for ctx in sorted(status.context):
 for mcp in ['MCP1','MCP2']:
  r=dict(context=int(ctx),program=mcp,status='NOT_ESTIMABLE',reason='',blocks=0,discovery_rho=np.nan,holdout_rho=np.nan,p=np.nan,ci_low=np.nan,ci_high=np.nan,valid_bootstrap=0,permutations=20000,bootstrap=2000,independent_cohort_validation=False)
  q=scope[(scope.context==ctx)&(scope.program==mcp)]
  if len(q)!=2 or not q.status.eq('FROZEN_RECONSTRUCTION_PASS').all():r['reason']='BOTH_GENE_SUPPORTED_RECONSTRUCTION_PASSES_REQUIRED';rows.append(r);continue
  for part in ['discovery','holdout']:
   v=s[(s.context==ctx)&s.program.eq(mcp)&s.partition.eq(part)];a=v[v.type.eq('typeA')][['sample','score']];b=v[v.type.eq('typeB')][['sample','score']];j=a.merge(b,on='sample',suffixes=('_A','_B'),validate='one_to_one');assert len(j)==len(a)==len(b)
   if part=='discovery':r['discovery_rho']=rho(j.score_A,j.score_B);continue
   r['blocks']=len(j)
   if len(j)<10:r['reason']='FEWER_THAN_10_INTERNAL_HELDOUT_BLOCKS';continue
   x=j.score_A.to_numpy();y=j.score_B.to_numpy();obs=rho(x,y)
   if not np.isfinite(obs):r['reason']='CONSTANT_HELDOUT_SCORE';continue
   a=rankdata(x);b=rankdata(y);a-=a.mean();b-=b.mean();den=np.linalg.norm(a)*np.linalg.norm(b);exceed=0
   for start in range(0,20000,1000):
    perm=np.array([rng.permutation(len(x)) for _ in range(min(1000,20000-start))]);null=b[perm]@a/den;exceed+=int((np.abs(null)>=abs(obs)-1e-12).sum())
   boots=[]
   for _ in range(2000):
    ix=rng.integers(0,len(x),len(x));z=rho(x[ix],y[ix])
    if np.isfinite(z):boots.append(z)
   r.update(status='ESTIMABLE_INTERNAL_HOLDOUT',reason='',holdout_rho=obs,p=(1+exceed)/20001,valid_bootstrap=len(boots))
   if len(boots)>=1900:r['ci_low'],r['ci_high']=map(float,np.quantile(boots,[.025,.975]))
  rows.append(r)
t=pd.DataFrame(rows);assert len(t)==28 and not t.duplicated(['context','program']).any();t['BH_q']=np.nan;ix=np.flatnonzero(np.isfinite(t.p));order=ix[np.argsort(t.p.iloc[ix].to_numpy(),kind='stable')];n=len(order)
if n:t.loc[order,'BH_q']=np.minimum(1,np.minimum.accumulate((t.p.iloc[order].to_numpy()*n/np.arange(1,n+1))[::-1])[::-1])
t.to_csv(O/'public_internal_holdout_results.tsv',sep='\t',index=False)
summary=dict(status='ALL28_FIXED_SLOTS_ACCOUNTED_INTERNAL_HOLDOUT_NOT_INDEPENDENT_VALIDATION',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(Path(__file__)),slot_scope_sha256=sha(P/'public_slot_scope.tsv'),fixed_slots=28,estimable_tests=n,BH_below_05=int((t.BH_q<.05).sum()),permutations=20000,bootstrap=2000,seed=2026101603,missing_not_zero=True,holdout_transform='Training gene normalization, signed NNLS weights, and all four training quality/abundance coefficients frozen; block means equal weight.',results_sha256=sha(O/'public_internal_holdout_results.tsv'))
(O/'inventory.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary))
