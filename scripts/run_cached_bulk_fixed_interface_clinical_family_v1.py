"""Fixed candidate bulk proxy, all declared clinical slots; real arrays RAM only."""
import os,sys,json,hashlib,importlib.util,warnings
from pathlib import Path
from collections import Counter
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.special import expit
from scipy.stats import chi2,rankdata,spearmanr
from scipy.linalg import qr
from statsmodels.stats.multitest import multipletests
from statsmodels.duration.hazard_regression import PHReg
import statsmodels.api as sm

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):h.update(b)
 return h.hexdigest()
def close(a,b,tol=1e-6):assert np.allclose(a,b,atol=tol,rtol=tol,equal_nan=True)
def clean(v):
 if isinstance(v,np.generic):return v.item()
 raise TypeError(type(v).__name__)
def load_module(p):
 spec=importlib.util.spec_from_file_location('accepted_bulk_readers',str(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def newton(fun,p):
 b=np.zeros(p);ll,g,I=fun(b)
 for it in range(80):
  step=np.linalg.solve(I,g);scale=1.
  for _ in range(35):
   bb=b+scale*step;lll,gg,II=fun(bb)
   if np.isfinite(lll) and lll>=ll-1e-10:break
   scale*=.5
  else:raise ArithmeticError('NEWTON_LINESEARCH_FAILED')
  b,ll,g,I=bb,lll,gg,II
  if np.max(np.abs(g))<1e-7 and np.max(np.abs(scale*step))<1e-7:break
 else:raise ArithmeticError('NEWTON_NO_CONVERGENCE')
 assert np.max(np.abs(g))<1e-5 and np.linalg.eigvalsh(I).min()>0
 return b,ll,I

def logistic_function(y,X):
 def fun(b):
  eta=X@b;mu=expit(eta);ll=float((y*eta-np.logaddexp(0,eta)).sum())
  return ll,X.T@(y-mu),X.T@((mu*(1-mu))[:,None]*X)
 return fun

def efron_function(t,e,X):
 # Explicit risk/death sets differ from statsmodels/cumulative producer risk sets.
 event_times=np.unique(t[e==1]);risk=(t[None,:]>=event_times[:,None]).astype(float)
 death=((t[None,:]==event_times[:,None])&(e[None,:]==1)).astype(float)
 d=death.sum(1).astype(int);group=np.repeat(np.arange(len(d)),d)
 frac=np.concatenate([np.arange(n)/n for n in d]);xx=(X[:,:,None]*X[:,None,:]).reshape(len(X),-1)
 def fun(b):
  eta=X@b;offset=eta.max();w=np.exp(eta-offset)
  R0=risk@w;D0=death@w;R1=risk@(w[:,None]*X);D1=death@(w[:,None]*X)
  R2=(risk@(w[:,None]*xx)).reshape(-1,X.shape[1],X.shape[1]);D2=(death@(w[:,None]*xx)).reshape(-1,X.shape[1],X.shape[1])
  den=R0[group]-frac*D0[group];mu=(R1[group]-frac[:,None]*D1[group])/den[:,None]
  ll=float(eta[e==1].sum()-(np.log(den)+offset).sum());g=X[e==1].sum(0)-mu.sum(0)
  I=((R2[group]-frac[:,None,None]*D2[group])/den[:,None,None]-mu[:,:,None]*mu[:,None,:]).sum(0)
  return ll,g,I
 return fun

def fit_clinical(data,ep,tcga,plan):
 kind='logistic' if ep=='response' else 'Cox';predictors=['interface','pan_T','pan_myeloid','interferon_disjoint','proliferation','pan_leukocyte']+(['age','sex','stage'] if tcga else [])
 cols=(['y'] if kind=='logistic' else ['t','e'])+predictors
 d=data[cols].replace([np.inf,-np.inf],np.nan).dropna()
 d=d[d.y.isin([0,1])] if kind=='logistic' else d[d.t.gt(0)&d.e.isin([0,1])]
 n=len(d);events=int(d.y.eq(1).sum()) if kind=='logistic' else int(d.e.sum());minor=min(events,n-events) if kind=='logistic' else events
 r=dict(status='UNESTIMABLE',n_units=n,events=events,reason='',effect=None,CI_low=None,CI_high=None,P=None)
 if tcga and float(data.stage.notna().mean())<.70:r['reason']='REQUIRED_STAGE_COVERAGE_LT70_PERCENT';return r
 if n<60 or minor<max(20,10*len(predictors)):r['reason']='FIXED_N_OR_EVENTS_PER_PREDICTOR_GATE';return r
 X=d[predictors].to_numpy(float);sd=X.std(0)
 if (sd<1e-10).any():r['reason']='CONSTANT_REQUIRED_PREDICTOR_NO_TERMS_DROPPED';return r
 X=(X-X.mean(0))/sd;Z=np.c_[np.ones(n),X] if kind=='logistic' else X
 if np.linalg.matrix_rank(Z)!=Z.shape[1] or np.linalg.cond(Z)>1e6:r['reason']='UNIDENTIFIABLE_COMPLETE_DESIGN';return r
 try:
  with warnings.catch_warnings(record=True) as w:
   warnings.simplefilter('always')
   if kind=='logistic':
    y=d.y.to_numpy(float);model=sm.GLM(y,Z,family=sm.families.Binomial());f=model.fit(maxiter=100)
    Zn=Z[:,[0]+list(range(2,Z.shape[1]))];f0=sm.GLM(y,Zn,family=sm.families.Binomial()).fit(maxiter=100)
    assert f.converged and f0.converged
    b,ll,I=newton(logistic_function(y,Z),Z.shape[1]);b0,ll0,I0=newton(logistic_function(y,Zn),Zn.shape[1]);index=1
   else:
    t=d.t.to_numpy(float);e=d.e.to_numpy(float);model=PHReg(t,Z,status=e,ties='efron');f=model.fit(disp=False)
    Zn=Z[:,1:];f0=PHReg(t,Zn,status=e,ties='efron').fit(disp=False)
    b,ll,I=newton(efron_function(t,e,Z),Z.shape[1]);b0,ll0,I0=newton(efron_function(t,e,Zn),Zn.shape[1]);index=0
   if not np.isfinite(f.params).all() or np.max(np.abs(f.params))>15:r['reason']='NONFINITE_OR_EXTREME_COEFFICIENT';return r
   close(b,f.params);close(b0,f0.params);close(ll,f.llf);close(ll0,f0.llf);close(np.linalg.inv(I),f.cov_params());close(np.linalg.inv(I0),f0.cov_params())
   se=float(np.sqrt(np.linalg.inv(I)[index,index]));coef=float(b[index]);lr=max(0.,2*(ll-ll0));P=float(chi2.sf(lr,1))
   assert np.isfinite([coef,se,P]).all() and se>0
   close(P,chi2.sf(max(0.,2*(f.llf-f0.llf)),1))
   r.update(status='ESTIMABLE_ASYMPTOTIC_SOURCE_EXPLORATION_INDEPENDENT_NUMERIC_REPLAY',reason='',log_effect=coef,se=se,effect=float(np.exp(coef)),CI_low=float(np.exp(coef-1.96*se)),CI_high=float(np.exp(coef+1.96*se)),P=P,LR=float(lr),loglik=float(ll),null_loglik=float(ll0),condition=float(np.linalg.cond(Z)),maximum_refit_coefficient_error=float(np.max(np.abs(b-f.params))),warning_classes=sorted({type(v.message).__name__ for v in w}),n_predictors=len(predictors))
   if kind=='Cox':
    sc=np.asarray(f.schoenfeld_residuals)@f.cov_params();take=(e==1)&np.isfinite(sc[:,0]);rho=spearmanr(sc[take,0],np.log(t[take])).statistic
    r['descriptive_schoenfeld_logtime_rho']=float(rho) if np.isfinite(rho) else None
    r['PH_assumption_not_proved']=True
   return r
 except AssertionError:raise
 except Exception as ex:r['reason']='MODEL_FAILURE_'+type(ex).__name__;return r

def coordinate(scores):
 names=['pan_T','pan_myeloid','interferon_disjoint','proliferation','pan_leukocyte'];d=scores[['ligand','receiver']+names].dropna();n=len(d)
 r=dict(n_units=n,status='UNESTIMABLE',reason='',partial_rank_effect=None,no_new_P_or_CI=True)
 if n<20 or n-len(names)-2<10:r['reason']='FIXED_N_OR_RESIDUAL_DF_GATE';return r
 values=d.to_numpy(float)
 if (np.ptp(values,axis=0)<1e-12).any():r['reason']='CONSTANT_REQUIRED_READOUT_OR_BACKGROUND';return r
 rank=np.column_stack([rankdata(values[:,j]) for j in range(values.shape[1])]);close(rank,d.rank(method='average').to_numpy())
 C=np.c_[np.ones(n),rank[:,2:]];C[:,1:]=(C[:,1:]-C[:,1:].mean(0))/C[:,1:].std(0)
 if np.linalg.matrix_rank(C)!=C.shape[1] or np.linalg.cond(C)>1e6:r['reason']='UNIDENTIFIABLE_BACKGROUND';return r
 Y=rank[:,:2];B=np.linalg.lstsq(C,Y,rcond=None)[0];E=Y-C@B
 Q,R=qr(C,mode='economic');EE=Y-Q@(Q.T@Y);close(E,EE)
 rho=float(E[:,0]@E[:,1]/np.linalg.norm(E[:,0])/np.linalg.norm(E[:,1]));rr=float(np.corrcoef(EE.T)[0,1]);close(rho,rr,1e-10)
 r.update(status='DESCRIPTIVE_PARTIAL_RANK_NO_INFERENTIAL_FAMILY',partial_rank_effect=rho,residual_df=n-C.shape[1]-1,condition=float(np.linalg.cond(C)))
 return r

def make_scores(expr,plan):
 panels=plan['control_panels'];genes=plan['interface_symbols']
 required=sorted(set(genes+sum(panels.values(),[])))
 if not set(required)<=set(expr.index):return None,'REQUIRED_FEATURE_NOT_AVAILABLE'
 e=expr.loc[required]
 if not np.isfinite(e.to_numpy(float)).all() or (e.nunique(axis=1)<=1).any():return None,'NONFINITE_OR_NONVARIABLE_REQUIRED_GENE'
 ranks=e.rank(axis=1,pct=True,method='average');s=pd.DataFrame(index=e.columns)
 s['ligand']=np.minimum(ranks.loc[genes[0]],ranks.loc[genes[1]]);s['receiver']=ranks.loc[genes[2]];s['interface']=s.ligand*s.receiver
 for name,g in panels.items():s[name]=ranks.loc[g].mean()
 # Full alternative percentile implementation for every source score, ties included.
 alt={g:rankdata(e.loc[g].to_numpy())/len(e.columns) for g in required}
 close(s.interface,np.minimum(alt[genes[0]],alt[genes[1]])*alt[genes[2]],1e-12)
 for name,g in panels.items():close(s[name],np.mean([alt[a] for a in g],axis=0),1e-12)
 return s,''

def main():
 plan=json.load(sys.stdin);base=Path(os.environ['PANCANCER_BULK_ROOT']);old=Path(os.environ['PANCANCER_OLD_BULK_RUN']);qual=Path(os.environ['PANCANCER_BULK_QUALIFICATION']);out=Path(os.environ['PANCANCER_RAM_OUTPUT'])
 assert out.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
 def save(n,v):(out/n).write_text(json.dumps(v,indent=2,default=clean,allow_nan=False)+'\n')
 assert sha(qual/'private_terminal.json')==plan['accepted_gene_axis_terminal_SHA256'] and sha(qual/'private_full_metadata_independent_audit.json')==plan['accepted_gene_axis_audit_SHA256']
 qt=json.loads((qual/'private_terminal.json').read_text())
 for name,v in qt['files'].items():assert (qual/name).stat().st_size==v['bytes'] and sha(qual/name)==v['SHA256']
 module_path=Path(__file__).with_name('explore_priority_immune_CM_bulk_v1.py');M=load_module(module_path)
 save('private_plan_before_expression_and_effects.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),reader_module_SHA256=sha(module_path),fixed_plan=plan))
 axes={r['cohort']:r for r in json.loads((qual/'private_all7_physical_bulk_feature_axis_qualification.json').read_text())}
 oldcounts=json.loads((qual/'private_all6_ICI_and33_TCGA_old_clinical_unit_counts.json').read_text())
 rows=[];desc=[];proof=[];allgenes=set(plan['interface_symbols']+sum(plan['control_panels'].values(),[]));helper=M.helpers(base)
 def process(co,expr,clin,tcga=False,unavailable=None):
  slots=[v for v in plan['clinical_slots'] if v['cohort']==co]
  if unavailable:
   for sl in slots:rows.append(dict(**sl,status='UNESTIMABLE',reason=unavailable,n_units=None,events=None,P=None,effect=None,CI_low=None,CI_high=None))
   desc.append(dict(cohort=co,status='UNESTIMABLE',reason=unavailable,no_new_P_or_CI=True));return
  cp=old/('SU2C_analysis_v1' if co=='SU2C_MARK' else 'analysis_v2')/(co+'_patient_scores_private.tsv')
  if not tcga:assert sha(cp)==next(v for v in oldcounts if v['cohort']==co)['prior_clinical_input_SHA256']
  prior=pd.read_csv(cp,sep='\t',index_col=0)
  assert set(prior.index)==set(clin.index) and expr.columns.is_unique and clin.index.is_unique
  prior=prior.reindex(clin.index)
  fields=['OS','OS.time','PFI','PFI.time','age','sex','stage'] if tcga else ['response','OS','OS_event','PFS','PFS_event']
  for field in fields:close(pd.to_numeric(clin[field],errors='coerce'),pd.to_numeric(prior[field],errors='coerce'),1e-12)
  sc,reason=make_scores(expr,plan)
  if sc is None:return process(co,None,None,tcga,reason)
  for name,previous in [('pan_T','pan_T'),('interferon_disjoint','interferon'),('proliferation','proliferation')]:close(sc[name],prior[previous],1e-10)
  proof.append(dict(cohort=co,old_patient_score_SHA256=sha(cp),exact_units_and_clinical_fields_restored=True,all_old_shared_background_scores_equal=True,all_new_scores_Pandas_vs_Scipy_rank_equal=True,n_units=len(sc)))
  desc.append(dict(cohort=co,**coordinate(sc)))
  for sl in slots:
   d=sc.copy();ep=sl['endpoint']
   if tcga:
    for nm in ['age','sex','stage']:d[nm]=clin[nm]
   if ep=='response':d['y']=clin.response
   else:
    d['t']=clin[ep+'.time' if tcga else ep];d['e']=clin[ep if tcga else ep+'_event']
   fit=fit_clinical(d,ep,tcga,plan);rows.append(dict(**sl,**fit))
  save('private_progress.json',dict(stage='COMPLETE_SOURCE_CLINICAL_FAMILY',source=co,complete_slots=len(rows),descriptive_sources=len(desc)))
 for co in plan['cohorts']:
  a=axes[co]
  reason='INTERFACE_NOT_IN_CURRENT_SERVER_PANEL' if not a['all_interface_roles_unique'] else 'REQUIRED_BACKGROUND_GENE_MISSING' if any(v['unique_physical_columns']!=v['required'] for v in a['controls'].values()) else None
  if reason:process(co,None,None,False,reason);continue
  expr,clin,meta=M.load_ici(base,co,helper)
  assert len(clin)==next(v for v in oldcounts if v['cohort']==co)['old_mapped_units']
  process(co,expr,clin)
  del expr,clin
 expr,clin,meta=M.load_tcga(base,allgenes)
 for cancer in plan['tcga_cancers']:
  cc=clin[clin.cancer.eq(cancer)];co='TCGA_'+cancer
  if len(cc)==0:process(co,None,None,True,'NO_DECLARED_PRIMARY_TUMOUR_UNITS');continue
  process(co,expr.loc[:,cc.index],cc,True)
 assert len(rows)==84 and len(desc)==39 and {(r['cohort'],r['endpoint']) for r in rows}=={(r['cohort'],r['endpoint']) for r in plan['clinical_slots']}
 p=np.array([r['P'] if r['P'] is not None else 1. for r in rows]);bh=M.bh(p);by=M.bh(p,True)
 close(bh,multipletests(p,method='fdr_bh')[1],1e-12);close(by,multipletests(p,method='fdr_by')[1],1e-12)
 for k,r in enumerate(rows):r.update(BH84=float(bh[k]) if r['P'] is not None else None,BY84=float(by[k]) if r['P'] is not None else None)
 save('private_complete84_source_clinical_results.json',rows);save('private_complete39_source_RNA_coordination_descriptions.json',desc);save('private_all_source_native_and_legacy_mapping_proof.json',proof)
 files={p.name:dict(bytes=p.stat().st_size,SHA256=sha(p)) for p in out.iterdir() if p.is_file()};assert sum(v['bytes'] for v in files.values())<4*1024**2
 save('private_terminal.json',dict(status='COMPLETE_ALL84_CLINICAL_AND39_DESCRIPTIVE_SOURCE_POSITIONS_FULL_NUMERIC_REPLAY_PENDING_SCIENTIFIC_REVIEW',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),files=files,registered_clinical_slots=84,registered_descriptive_positions=39,
  all_fitted_full_null_models_replayed=True,no_real_expression_patient_scores_or_model_arrays_persisted=True,new_raw_RNA_or_denied_clinical_CSV_download=False,clinical_model_assumptions='Asymptotic likelihood and Wald CI; not exact finite, prospective prediction, causal or cell-specific validation',goal_status='active'))
 print(json.dumps(dict(status='COMPLETE_ALL84_AND39',terminal_SHA256=sha(out/'private_terminal.json'))))
if __name__=='__main__':main()
