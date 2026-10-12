"""Synthetic formula fixtures and predeclared null calibration, no human expression."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,importlib.util
import numpy as np,pandas as pd
from scipy.stats import beta as beta_dist

def main():
 C=Path(os.environ['PANCANCER_CR2_CORE']);spec=importlib.util.spec_from_file_location('core',C);core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core);O=Path(os.environ['PANCANCER_CR2_CALIBRATION']);assert not O.exists();O.mkdir();(O/'synthetic_fixtures').mkdir();planp=Path(os.environ['PANCANCER_CR2_CALIBRATION_PLAN']);plan=json.loads(planp.read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();assert sha(__file__)==plan['calibration_sha256'] and sha(C)==plan['core_sha256'];start=datetime.now(timezone.utc).isoformat();rng=np.random.default_rng(plan['seed']);results=[];fixtures=[]
 specs=[('singleton30',np.ones(30,int),4),('paired12',np.full(12,2),6),('unequal10',np.array([1,1,2,2,2,3,3,3,4,5]),6),('paired20',np.full(20,2),7),('unequal20',np.tile([1,2,3,4],5),7),('many_controls12',np.full(12,2),12)]
 for name,sizes,k in specs:
  n=int(sizes.sum());blocks=np.repeat(np.arange(len(sizes)),sizes);Z=np.column_stack([np.ones(n),rng.normal(size=(n,k))]);x=core.standardize(.3*Z[:,1]+rng.normal(size=n));Y=.4*x+Z@rng.normal(size=k+1)+rng.normal(size=n);prep=core.prepare_cr2(Z,x,blocks,min_df=0);assert 'coef_vector' in prep;out=core.estimate_cr2(prep,Y);fixture=pd.DataFrame({'block':blocks,'x':x,'y':Y});
  for j in range(Z.shape[1]):fixture['Z'+str(j)]=Z[:,j]
  fixture.to_csv(O/'synthetic_fixtures'/(name+'.tsv'),sep='\t',index=False);fixtures.append(dict(case=name,n=n,blocks=len(sizes),rank=prep['rank'],beta=float(out['beta'][0]),se=float(out['se'][0]),df=float(prep['df']),p=float(out['p'][0]),lower=float(out['lower'][0]),upper=float(out['upper'][0])))
  B=plan['null_draws'];K=plan['joint_TF_dimension'];G=len(sizes);rho=.7;block_scale=np.exp(.5*core.standardize(np.bincount(blocks,weights=x)/sizes));eta=rng.normal(size=(n,B));shared=rng.normal(size=(G,B));errors={
   'iid_normal':eta,
   'cluster_intercept_rho07':np.sqrt(rho)*shared[blocks]+np.sqrt(1-rho)*eta,
   'cluster_heteroskedastic_x':block_scale[blocks,None]*(np.sqrt(rho)*shared[blocks]+np.sqrt(1-rho)*eta),
   'cluster_heavytail_t5':np.sqrt(rho)*rng.standard_t(5,size=(G,B))[blocks]/np.sqrt(5/3)+np.sqrt(1-rho)*eta}
  for error,eps in errors.items():
   p=core.estimate_cr2(prep,eps)['p'];rejections=int((p<.05).sum());lo=float(beta_dist.ppf(.025,rejections,B-rejections+1)) if rejections else 0;hi=float(beta_dist.ppf(.975,rejections+1,B-rejections)) if rejections<B else 1;results.append(dict(case=name,error=error,n=n,blocks=G,rank=prep['rank'],df=prep['df'],inference_gate_df4=prep['df']>=4,draws=B,rejections=rejections,rejection_rate=rejections/B,binomial_lower=lo,binomial_upper=hi,joint=False))
  # Ten correlated receiver TF outcomes; whole joint vector retained per null draw.
  common=errors['cluster_intercept_rho07'];ps=[]
  for j in range(K):
   extra=np.sqrt(rho)*rng.normal(size=(G,B))[blocks]+np.sqrt(1-rho)*rng.normal(size=(n,B));eps=np.sqrt(.6)*common+np.sqrt(.4)*extra;ps.append(core.estimate_cr2(prep,eps)['p'])
  omni=np.minimum(1,K*np.min(ps,axis=0));r=int((omni<.05).sum());results.append(dict(case=name,error='correlated_TF_Bonferroni_union',n=n,blocks=G,rank=prep['rank'],df=prep['df'],inference_gate_df4=prep['df']>=4,draws=B,rejections=r,rejection_rate=r/B,binomial_lower=float(beta_dist.ppf(.025,r,B-r+1)) if r else 0,binomial_upper=float(beta_dist.ppf(.975,r+1,B-r)) if r<B else 1,joint=True))
 pd.DataFrame(fixtures).to_csv(O/'synthetic_python_reference_results.tsv',sep='\t',index=False);pd.DataFrame(results).to_csv(O/'null_calibration_summary.tsv',sep='\t',index=False);j=dict(status='SYNTHETIC_PYTHON_CALIBRATION_COMPLETE_R_REFERENCE_PENDING',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),core_sha256=sha(C),calibration_sha256=sha(__file__),plan_sha256=sha(planp),seed=plan['seed'],null_draws=plan['null_draws'],joint_dimension=K,synthetic_cases=len(specs),summary_rows=len(results),real_expression_p_computed=False,formula_R_reference_pending=True,outputs={str(p.relative_to(O)):sha(p) for p in O.rglob('*') if p.is_file()});(O/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps({k:v for k,v in j.items() if k!='outputs'}));print(pd.DataFrame(results)[['case','error','df','inference_gate_df4','rejection_rate']].to_string(index=False))
if __name__=='__main__':main()
