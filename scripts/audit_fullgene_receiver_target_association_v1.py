"""Independent QR association audit, prior-score replay, and two resampling examples.

No producer import. All1860 model statuses/coefficient/CR1/BH checks;
all325 score recipes versus accepted fullgene units and official prior;
positive/negative examples independently replay 2000case-bootstrap/permutations.
"""
from pathlib import Path
from datetime import datetime,timezone
import os,json,hashlib,gc
import numpy as np,pandas as pd,pyreadr
from scipy.stats import t as tdist
from scipy.linalg import qr


def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()


def zscore(a):
 sd=np.std(a,axis=0,ddof=1);v=np.zeros_like(a,float);good=sd>1e-12;v[:,good]=(a[:,good]-a[:,good].mean(axis=0))/sd[good]
 return v,sd


def qfit(a,y):
 n,k=len(y),a.shape[1]+1
 if not np.isfinite(a).all():return 'MISSING_REQUIRED_COVARIATE',None
 if n<=k:return 'INSUFFICIENT_RESIDUAL_DF',None
 if np.std(y,ddof=1)<=1e-12:return 'TARGET_EXPRESSION_SCORE_CONSTANT',None
 v,sd=zscore(a);X=np.column_stack([np.ones(n),v]);Y=(y-y.mean())/y.std(ddof=1);q,r,piv=qr(X,mode='economic',pivoting=True);tol=np.linalg.norm(X,2)*1e-10;rank=int((abs(np.diag(r))>tol).sum())
 if rank<k:return 'REQUIRED_DESIGN_NOT_FULL_RANK',None
 if np.linalg.cond(X)>=1e8:return 'REQUIRED_DESIGN_ILL_CONDITIONED',None
 leverage=(q*q).sum(axis=1)
 if leverage.max()>.95:return 'HIGH_INFLUENCE_DESIGN',None
 coef=np.zeros(k);coef[piv]=np.linalg.solve(r,q.T@Y);res=Y-X@coef;bread=np.linalg.solve(X.T@X,np.eye(k));meat=np.zeros((k,k))
 for i in range(n):meat+=np.outer(X[i]*res[i],X[i]*res[i])
 cov=bread@meat@bread*(n/(n-k));se=np.sqrt(max(cov[1,1],0));return 'ESTIMABLE_EXPLORATORY_ASSOCIATION',(X,Y,coef,se)


def main():
 start=datetime.now(timezone.utc).isoformat();R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'immune_same_name_family_source_review_v1';O=P/'fullgene_receiver_target_association_v1';F=P/'fullgene_receiver_target_readiness_v1';assert not (O/'independent_audit.json').exists();inv=json.loads((O/'inventory.json').read_text());plan=json.loads((R/'fullgene_receiver_target_code_v1/fullgene_receiver_target_association_plan_v1.json').read_text())
 assert sha(__file__)==plan['code_sha256'][Path(__file__).name]
 for n,h in inv['outputs'].items():assert sha(O/n)==h
 d=pd.read_csv(O/'all1860_target_association_models.tsv',sep='\t',keep_default_na=False);assert len(d)==1860 and d.test_id.nunique()==620 and not d[['test_id','variant']].duplicated().any()
 details={r['test_id']:r for r in json.loads((O/'score_provenance.json').read_text())};cases=pd.read_csv(P/'fullgene_receiver_response_preflight_v1/all_interface_original_node_preflight.tsv',sep='\t');expected=['NODE%04d'%(i+1) for i in range(620)];assert set(d.test_id)==set(expected)
 prior=next(iter(pyreadr.read_r(str(F/'verified_prior_inputs_v2/ligand_target_matrix_nsga2r_final.rds')).values()));prior.index=prior.index.astype(str);prior.columns=prior.columns.astype(str)
 score_checks=0;numeric=0;unestimated=0
 for sid in cases.gse_id.unique():
  b=dict(np.load(F/'private_source_unit_arrays'/(sid+'.npz'),allow_pickle=False));nodes=list(b['nodes']);genes=pd.Index(b['genes']);samples=pd.Index(b['samples']);lg=b['mean_log1p_CP10K'].reshape(len(nodes),len(samples),len(genes))
  for test in d[d.gse_id.eq(sid)].test_id.unique():
   sub=d[d.test_id.eq(test)];first=sub.iloc[0];idx=int(test[4:])-1;original=cases.iloc[idx]
   for col in ['interface_id','gse_id','stratum','sender_family','receiver_family','ligand','receptor','receiver_original_node','LR_status']:assert first[col]==original[col]
   path=O/'private_designs'/(test+'.npz')
   if not path.exists():
    assert all(sub.status.eq('PREFLIGHT_'+original.status));continue
   z=dict(np.load(path,allow_pickle=False));si=samples.get_indexer(z['sample_keys']);assert (si>=0).all();gi=genes.get_indexer(z['target_genes']);assert (gi>=0).all();raw=lg[nodes.index(first.receiver_original_node),si][:,gi]
   assert np.array_equal(raw,z['raw_target_gene_expression'])
   w=prior.loc[z['target_genes'],first.ligand].to_numpy(float);w=w/w.sum();assert np.allclose(w,z['target_prior_weights'],rtol=1e-12,atol=1e-14)
   standardized,sd=zscore(raw);y=sum(standardized[:,j]*w[j] for j in range(len(w)));assert np.allclose(y,z['raw_target_expression_score'],rtol=1e-10,atol=1e-10)
   assert len(set(z['patient_blocks']))==len(z['patient_blocks']);score_checks+=1
   for row in sub.itertuples():
    a=z['raw_a_'+row.variant];status,values=qfit(a,y);assert row.status==status,(test,row.variant,status,row.status)
    if values is None:unestimated+=1;continue
    X,Y,coef,se=values;assert np.allclose(X,z['X_'+row.variant],rtol=1e-8,atol=1e-8);assert np.allclose(Y,z['y_'+row.variant],rtol=1e-10,atol=1e-10)
    assert np.isclose(coef[1],float(row.sender_beta),rtol=1e-6,atol=1e-7);assert np.isclose(se,float(row.CR1_se),rtol=1e-5,atol=1e-7)
    crit=tdist.ppf(.975,len(y)-1);assert np.isclose(coef[1]-crit*se,float(row.CR1_ci_low),rtol=1e-5,atol=1e-6);numeric+=1
  del b,lg;gc.collect()
 # Diagnostic BH maintains all620 rows per variant, unavailable tests still NA.
 for variant,g in d.groupby('variant'):
  p=pd.to_numeric(g.residual_permutation_p,errors='coerce').fillna(1).to_numpy();order=np.argsort(p,kind='stable');qs=np.minimum.accumulate((p[order]*620/np.arange(1,621))[::-1])[::-1];adjusted=np.empty(620);adjusted[order]=np.minimum(qs,1)
  for j,(i,r) in enumerate(g.iterrows()):
   if r.residual_permutation_p=='NA':assert r.diagnostic_BH_q=='NA'
   else:assert np.isclose(float(r.diagnostic_BH_q),adjusted[j],rtol=1e-12,atol=1e-12)
 # Independently replay one positive/one negative estimated case; use M1 when possible.
 good=d[d.status.eq('ESTIMABLE_EXPLORATORY_ASSOCIATION')].copy();good['beta']=pd.to_numeric(good.sender_beta);examples=[]
 for sign in [1,-1]:
  pool=good[good.beta*sign>0];preferred=pool[pool.variant.eq('M1_abundance_identity')];pool=preferred if len(preferred) else pool
  if not len(pool):examples.append(dict(sign=sign,status='NO_ESTIMABLE_CASE_WITH_THIS_DIRECTION'));continue
  row=pool.iloc[(pool.beta.abs()).argmax()];z=dict(np.load(O/'private_designs'/(row.test_id+'.npz'),allow_pickle=False));X=z['X_'+row.variant];Y=z['y_'+row.variant];n,k=X.shape;B=plan['resampling']['repeats'];rng=np.random.default_rng(int(row.resampling_seed));restricted=np.delete(X,1,axis=1);nullfit=restricted@np.linalg.lstsq(restricted,Y,rcond=None)[0];res=Y-nullfit;beta=float(row.sender_beta);extreme=0
  v=np.linalg.lstsq(X,np.eye(n),rcond=None)[0][1]
  for j in range(B):
   perm=rng.permutation(n);stat=np.dot(v,nullfit+res[perm]);extreme+=abs(stat)>=abs(beta)-1e-12
  pp=(1+extreme)/(B+1);assert np.isclose(pp,float(row.residual_permutation_p),atol=1e-12)
  pick=rng.integers(0,n,(B,n));boots=[]
  for idx in pick:
   xx=X[idx];yy=Y[idx];cond=np.linalg.cond(xx)
   if cond>=1e8:continue
   q,r,piv=qr(xx,mode='economic',pivoting=True)
   if ((q*q).sum(axis=1)).max()>.95:continue
   bb=np.zeros(k);bb[piv]=np.linalg.solve(r,q.T@yy);boots.append(bb[1])
  assert len(boots)==int(float(row.bootstrap_valid)),(row.test_id,len(boots),row.bootstrap_valid)
  if len(boots)>=plan['resampling']['minimum_valid_bootstrap']:
   ci=np.quantile(boots,[.025,.975]);assert np.allclose(ci,[float(row.bootstrap_ci_low),float(row.bootstrap_ci_high)],rtol=1e-4,atol=1e-5)
  examples.append(dict(test_id=row.test_id,variant=row.variant,direction=sign,status='PASS_INDEPENDENT_QR_PERMUTATION_AND_CASE_BOOTSTRAP_REPLAY',permutation_p=pp,valid_bootstraps=len(boots)))
 failures=d[~d.status.eq('ESTIMABLE_EXPLORATORY_ASSOCIATION')];example_failure=failures[~failures.status.str.startswith('PREFLIGHT_')].head(1)[['test_id','variant','status']].to_dict('records')
 result=dict(status='PASS_ALL1860_EXPLORATORY_MODELS_TARGET_SCORE_QR_CR1_BH_AND_SIGNED_RESAMPLING_EXAMPLES',started_UTC=start,finished_UTC=datetime.now(timezone.utc).isoformat(),auditor_sha256=sha(__file__),producer_inventory_sha256=sha(O/'inventory.json'),score_replays=score_checks,estimated_models_replayed=numeric,unestimated_designs_replayed=unestimated,model_rows=1860,resampling_examples=examples,unestimable_examples=example_failure,input_controls_construction_not_yet_independently_replayed=True,independent_biological_validation=False,scientific_release=False,axis_validation_complete=False)
 (O/'independent_numerical_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
