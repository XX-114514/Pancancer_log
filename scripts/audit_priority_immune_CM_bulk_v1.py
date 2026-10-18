#!/usr/bin/env python3
"""Refit all accepted designs; independently evaluate logistic/Efron likelihoods and full FDR."""
import argparse,hashlib,json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.special import expit
from scipy.stats import chi2
from statsmodels.duration.hazard_regression import PHReg
import statsmodels.api as sm

def direct_efron(t,e,X,b):
 event_times=np.unique(t[e==1]);risk=(t[None,:]>=event_times[:,None]).astype(float)
 deaths=((t[None,:]==event_times[:,None])&(e[None,:]==1)).astype(float)
 counts=deaths.sum(1).astype(int);eta=X@b;offset=eta.max();w=np.exp(eta-offset)
 R0=risk@w;D0=deaths@w;R1=risk@(w[:,None]*X);D1=deaths@(w[:,None]*X)
 xx=(X[:,:,None]*X[:,None,:]).reshape(len(X),-1)
 R2=(risk@(w[:,None]*xx)).reshape(-1,X.shape[1],X.shape[1]);D2=(deaths@(w[:,None]*xx)).reshape(-1,X.shape[1],X.shape[1])
 group=np.repeat(np.arange(len(counts)),counts);f=np.concatenate([np.arange(n)/n for n in counts])
 den=R0[group]-f*D0[group];mean=(R1[group]-f[:,None]*D1[group])/den[:,None]
 ll=float(eta[e==1].sum()-(np.log(den)+offset).sum());gradient=X[e==1].sum(0)-mean.sum(0)
 info=((R2[group]-f[:,None,None]*D2[group])/den[:,None,None]-mean[:,:,None]*mean[:,None,:]).sum(0)
 return ll,gradient,info

def qvalues(values,by=False):
 values=np.asarray(values);order=np.argsort(values);m=len(values);adjust=np.sum(1/np.arange(1,m+1)) if by else 1
 out=np.ones(m);cap=1.
 for pos in range(m-1,-1,-1):
  cap=min(cap,values[order[pos]]*m/(pos+1)*adjust);out[order[pos]]=cap
 return out

def main(a):
 assert not a.output.exists();a.output.mkdir(mode=0o700)
 terminal=json.loads((a.run/'terminal_public.json').read_text());assert terminal['n_slots']==304
 rows=pd.read_csv(a.run/'all_frozen_model_slots_private.tsv',sep='\t');designs=json.loads((a.run/'fitted_designs_private.json').read_text())
 assert len(rows)==304 and len(designs)==sum(rows.status.eq('exploratory_model_numeric_audit_pass'))
 errors=[]
 for i,entry in enumerate(designs):
  row=rows.iloc[entry['record_index']];d=pd.DataFrame(entry['data']);X=d[entry['predictors']].to_numpy(float);X=(X-X.mean(0))/X.std(0)
  if entry['kind']=='Cox':
   t=d.t.to_numpy();e=d.e.to_numpy();model=PHReg(t,X,status=e,ties='efron');fit=model.fit(method='newton',disp=False)
   null=PHReg(t,X[:,1:],status=e,ties='efron').fit(method='newton',disp=False)
   ll,g,info=direct_efron(t,e,X,fit.params);cov=np.linalg.inv(info);coef=fit.params[0];se=np.sqrt(cov[0,0])
  else:
   X=np.column_stack([np.ones(len(X)),X]);y=d.y.to_numpy();b=np.zeros(X.shape[1])
   for _ in range(100):
    pr=expit(X@b);g=X.T@(y-pr);info=X.T@((pr*(1-pr))[:,None]*X);step=np.linalg.solve(info,g);b+=step
    if np.max(np.abs(step))<1e-10:break
   eta=X@b;pr=expit(eta);g=X.T@(y-pr);info=X.T@((pr*(1-pr))[:,None]*X);ll=float((y*eta-np.logaddexp(0,eta)).sum());coef=b[1];se=np.sqrt(np.linalg.inv(info)[1,1])
   base=X[:,[0]+list(range(2,X.shape[1]))];bb=np.zeros(base.shape[1])
   for _ in range(100):
    pp=expit(base@bb);step=np.linalg.solve(base.T@((pp*(1-pp))[:,None]*base),base.T@(y-pp));bb+=step
    if np.max(np.abs(step))<1e-10:break
   etan=base@bb;null_ll=float((y*etan-np.logaddexp(0,etan)).sum())
  if entry['kind']=='Cox':null_ll=float(null.llf)
  p=float(chi2.sf(max(0,2*(ll-null_ll)),1))
  assert np.max(np.abs(g))<1e-4 and np.linalg.eigvalsh(info).min()>0
  assert abs(coef-row.log_effect)<1e-5 and abs(se-row.se)<1e-5 and abs(ll-row.loglik)<1e-5 and abs(null_ll-row.null_loglik)<1e-5 and abs(p-row.p)<1e-5
  assert np.allclose([np.exp(coef),np.exp(coef-1.96*se),np.exp(coef+1.96*se)],[row.effect,row.CI_low,row.CI_high],atol=1e-5,rtol=1e-5)
  errors.append(max(abs(coef-row.log_effect),abs(se-row.se),abs(ll-row.loglik),abs(p-row.p)))
  if i%20==0:print(json.dumps(dict(phase='independent_design_checked',done=i+1,total=len(designs))),flush=True)
 for family,n in [('TCGA_survival',264),('ICI_response',10),('ICI_survival',20),('ICI_increment',10)]:
  ss=rows[rows.family.eq(family)];assert len(ss)==n;p=ss.p.fillna(1).to_numpy()
  assert np.allclose(qvalues(p),ss.BH_q,atol=1e-12,rtol=1e-12) and np.allclose(qvalues(p,True),ss.BY_q,atol=1e-12,rtol=1e-12)
 result=dict(status='ALL304_SLOTS_AND_ALL_FITTED_NUMERIC_DESIGNS_INDEPENDENTLY_REPLAYED',finished_utc=datetime.now(timezone.utc).isoformat(),n_slots=304,n_fitted_replayed=len(designs),maximum_numeric_difference=max(errors) if errors else 0.,full_family_BH_BY_pass=True,
  producer_terminal_sha256=hashlib.sha256((a.run/'terminal_public.json').read_bytes()).hexdigest(),auditor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
  scope='Independent numeric refit/likelihood/score/covariance and FDR. Original clinical source fields reviewed separately; not independent data collection or clinical replication.',no_scientific_promotion=True)
 (a.output/'independent_audit_public.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);p.add_argument('--output',type=Path,required=True);main(p.parse_args())
