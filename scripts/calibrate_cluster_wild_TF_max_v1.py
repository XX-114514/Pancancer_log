"""Joint TF covariance-preserving bootstrap calibration, synthetic data only."""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,importlib.util
import numpy as np,pandas as pd
from scipy.stats import beta as bd
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def module(name,p):
 s=importlib.util.spec_from_file_location(name,p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

def main():
 core=module('cr2',Path(os.environ['PANCANCER_CR2_CORE']));wild=module('wild',Path(os.environ['PANCANCER_WILD_CORE']));S=Path(os.environ['PANCANCER_CR2_CALIBRATION']);O=Path(os.environ['PANCANCER_WILD_CALIBRATION']);planp=Path(os.environ['PANCANCER_WILD_CALIBRATION_PLAN']);plan=json.loads(planp.read_text());assert sha(__file__)==plan['producer_sha256'] and sha(Path(os.environ['PANCANCER_WILD_CORE']))==plan['wild_core_sha256'] and not O.exists();O.mkdir();started=datetime.now(timezone.utc).isoformat();rng=np.random.default_rng(plan['seed']);rows=[];checks=[]
 for f in sorted((S/'synthetic_fixtures').glob('*.tsv')):
  d=pd.read_csv(f,sep='\t');name=f.stem;Z=d[[c for c in d if c.startswith('Z')]].to_numpy();x=d.x.to_numpy();blocks=d.block.to_numpy();n=len(x);G=len(set(blocks));cr=core.prepare_cr2(Z,x,blocks,min_df=0);pr=wild.prepare(Z,x,blocks,cr);assert pr['status']=='prepared';Y=np.column_stack([d.y.to_numpy(),.6*d.y.to_numpy()+rng.normal(size=n),rng.normal(size=n)]);weights=rng.choice(np.array([-1.,1.]),size=(64,G));fast=wild.bootstrap_max_statistics(pr,Y,weights);Q=pr['Q'];e0=Y-Q@(Q.T@Y);adjusted=np.zeros_like(e0)
  for ix,A in zip(pr['groups'],pr['adjustments']):adjusted[ix]=A@e0[ix]
  X=np.column_stack((Z,x));inverse=np.linalg.pinv(X,rcond=1e-10);H=X@inverse;slow=[]
  for w in weights:
   ys=adjusted*w[blocks,None];coef=(inverse@ys)[-1];res=ys-X@(inverse@ys);score=[]
   for ix in pr['groups']:
    vals,V=np.linalg.eigh(np.eye(len(ix))-H[np.ix_(ix,ix)]);A=(V/np.sqrt(vals))@V.T;score.append(inverse[-1,ix]@A@res[ix])
   se=np.sqrt(np.sum(np.array(score)**2,axis=0));slow.append(np.max(np.divide(abs(coef),se,out=np.zeros_like(coef),where=se>1e-14)))
  assert np.allclose(fast,slow,rtol=2e-7,atol=2e-7);checks.append(dict(case=name,bootstrap_refits=64,max_difference=float(np.max(abs(fast-slow)))));B=plan['null_draws'];K=plan['joint_TFs'];scale=np.exp(.5*core.standardize(np.bincount(blocks,weights=x)/np.bincount(blocks)));rho=.7
  for kind in ['iid','cluster_rho07','cluster_heteroskedastic_x','cluster_t5','joint_correlated_heteroskedastic_TFs']:
   rejects=0;dimension=K if kind.startswith('joint_') else 1
   for rep in range(B):
    if kind=='iid':eps=rng.normal(size=(n,dimension))
    elif kind=='cluster_t5':eps=np.sqrt(rho)*rng.standard_t(5,size=(G,dimension))[blocks]/np.sqrt(5/3)+np.sqrt(1-rho)*rng.normal(size=(n,dimension))
    else:
     common=np.sqrt(rho)*rng.normal(size=(G,1))[blocks]+np.sqrt(1-rho)*rng.normal(size=(n,1));extra=np.sqrt(rho)*rng.normal(size=(G,dimension))[blocks]+np.sqrt(1-rho)*rng.normal(size=(n,dimension));eps=np.sqrt(.6)*common+np.sqrt(.4)*extra
     if 'heteroskedastic' in kind:eps*=scale[blocks,None]
    estimate=core.estimate_cr2(cr,eps);result=wild.test(pr,eps,estimate['stat'],rng,draws=plan['bootstrap_draws'],exact_G_threshold=0,batch=128);rejects+=result['p']<.05
   rows.append(dict(case=name,error=kind,n=n,blocks=G,df=cr['df'],df4_eligible=cr['df']>=4,TFs=dimension,null_draws=B,wild_draws=plan['bootstrap_draws'],rejections=rejects,rate=rejects/B,binomial_lower=float(bd.ppf(.025,rejects,B-rejects+1)) if rejects else 0,binomial_upper=float(bd.ppf(.975,rejects+1,B-rejects)) if rejects<B else 1));print(json.dumps(rows[-1]),flush=True)
 pd.DataFrame(rows).to_csv(O/'null_calibration_summary.tsv',sep='\t',index=False);pd.DataFrame(checks).to_csv(O/'independent_naive_refit_checks.tsv',sep='\t',index=False);j=dict(status='SYNTHETIC_WILD_MAX_CALIBRATION_COMPLETE',started_UTC=started,finished_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),plan_sha256=sha(planp),all_naive_fullX_refits_pass=True,naive_refits=6*64,summary_rows=len(rows),synthetic_only_no_RNA_p=True,actual_source_or_receptor_validation=False,outputs={p.name:sha(p) for p in O.iterdir() if p.is_file()});(O/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps(j))
if __name__=='__main__':main()
