"""Whole fixed CM proxy evidence, exact balance meaning and internal CV diagnostics."""
import os,sys,json,hashlib,warnings,importlib.util
from pathlib import Path
from datetime import datetime,timezone
from collections import Counter
import numpy as np,pandas as pd
from scipy.special import expit
from scipy.stats import percentileofscore,chi2
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score,log_loss,brier_score_loss
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def close(a,b,tol=1e-8):assert np.allclose(a,b,atol=tol,rtol=tol,equal_nan=True)
def seed(s):return int(hashlib.sha256(s.encode()).hexdigest()[:8],16)
def load_module(p):
 spec=importlib.util.spec_from_file_location('accepted_CM_proxy_bulk_readers',str(p));m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def fit(X,y):
 X=np.asarray(X,float);y=np.asarray(y,float)
 if len(y)<60 or min((y==1).sum(),(y==0).sum())<10:return None,'TRAIN_SIZE_OR_CLASS_GATE'
 mu=X.mean(0);sd=X.std(0)
 if not np.isfinite(X).all() or (sd<1e-10).any():return None,'NONFINITE_OR_CONSTANT_REQUIRED_PREDICTOR'
 Z=np.c_[np.ones(len(X)),(X-mu)/sd]
 if np.linalg.matrix_rank(Z)!=Z.shape[1] or np.linalg.cond(Z)>1e6:return None,'UNIDENTIFIABLE_FULL_DESIGN'
 try:
  with warnings.catch_warnings(record=True) as w:
   warnings.simplefilter('always');f=sm.GLM(y,Z,family=sm.families.Binomial()).fit(maxiter=200,tol=1e-12)
  if not f.converged:return None,'GLM_NOT_CONVERGED'
  b=np.zeros(Z.shape[1])
  for it in range(100):
   eta=Z@b;pr=expit(eta);g=Z.T@(y-pr);I=Z.T@((pr*(1-pr))[:,None]*Z);step=np.linalg.solve(I,g)
   ll=float((y*eta-np.logaddexp(0,eta)).sum());scale=1.
   for _ in range(35):
    bb=b+scale*step;ee=Z@bb;newll=float((y*ee-np.logaddexp(0,ee)).sum())
    if np.isfinite(newll) and newll>=ll-1e-10:break
    scale*=.5
   else:return None,'NEWTON_LINE_SEARCH_FAILED'
   b=bb
   if np.max(np.abs(scale*step))<1e-8:break
  else:return None,'NEWTON_NOT_CONVERGED'
  eta=Z@b;pr=expit(eta);g=Z.T@(y-pr);I=Z.T@((pr*(1-pr))[:,None]*Z)
  if not np.isfinite(b).all() or np.linalg.eigvalsh(I).min()<=0:return None,'NONFINITE_OR_SINGULAR_INFORMATION'
  close(b,f.params,1e-6);close(np.linalg.inv(I),f.cov_params(),1e-6);close(pr,f.fittedvalues,1e-7)
  ll=float((y*eta-np.logaddexp(0,eta)).sum());close(ll,f.llf,1e-7);assert np.max(np.abs(g))<1e-5
  return dict(beta=b,beta_GLM=np.asarray(f.params),mu=mu,sd=sd,cov=np.linalg.inv(I),ll=ll,condition=float(np.linalg.cond(Z)),maximum_coefficient_error=float(np.max(np.abs(b-f.params))),warnings=sorted({type(x.message).__name__ for x in w})),''
 except AssertionError:raise
 except Exception as e:return None,'MODEL_FAILURE_'+type(e).__name__

def predict(model,X):
 Z=np.c_[np.ones(len(X)),(np.asarray(X)-model['mu'])/model['sd']]
 p=expit(Z@model['beta']);p2=expit(Z@model['beta_GLM']);close(p,p2,1e-7);return p

def cdf_values(train,query):
 train=np.sort(np.asarray(train,float));query=np.asarray(query,float)
 assert np.isfinite(train).all() and np.isfinite(query).all() and len(train)>0
 left=np.searchsorted(train,query,'left');right=np.searchsorted(train,query,'right')
 out=(left+right+(left<right))/(2*len(train));alt=np.asarray(percentileofscore(train,query,kind='rank'))/100
 close(out,alt,1e-12);return out

def raw_to_scores(expr,pool,panels,arms):
 required=sorted(set(panels[arms[0]]+panels[arms[1]]+panels['pan_T']))
 if not set(required)<=set(expr.index):return None,'REQUIRED_GENE_NOT_MEASURED'
 raw=expr.loc[required]
 if not np.isfinite(raw.to_numpy(float)).all() or (raw.loc[:,pool].nunique(axis=1)<=1).any():return None,'NONFINITE_OR_CONSTANT_REQUIRED_TRAIN_GENE'
 rank={g:cdf_values(raw.loc[g,pool],raw.loc[g]) for g in required}
 score=pd.DataFrame(index=expr.columns)
 score['armA']=np.mean([rank[g] for g in panels[arms[0]]],axis=0);score['armB']=np.mean([rank[g] for g in panels[arms[1]]],axis=0)
 score['pan_T']=np.mean([rank[g] for g in panels['pan_T']],axis=0);score['joint']=np.minimum(score.armA,score.armB)
 # Every training-RNA score also independently matches pandas average rank.
 direct=raw.loc[:,pool].rank(axis=1,pct=True,method='average')
 for g in required:close(rank[g][expr.columns.get_indexer(pool)],direct.loc[g],1e-12)
 return score,''

def metrics(y,p):
 eps=np.finfo(float).eps;pc=np.clip(p,eps,1-eps);loss=-y*np.log(pc)-(1-y)*np.log1p(-pc)
 ll=float(loss.mean());auc=float(roc_auc_score(y,p));brier=float(((y-p)**2).mean())
 close(ll,log_loss(y,p,labels=[0,1]),1e-10);close(brier,brier_score_loss(y,p),1e-12)
 return dict(logloss=ll,AUC=auc,Brier=brier),loss

def read_text_index(p):
 d=pd.read_csv(p,sep='\t',dtype=str);return d.set_index(d.columns[0],verify_integrity=True)

def SU2C(root,M):
 src=root/'SU2C_MARK_panel_ready_v1';expr=pd.read_csv(src/'fixed_panel_TPM_private.tsv.gz',sep='\t',index_col=0)
 proof=json.loads((src/'panel_derivation_receipt_public.json').read_text());assert sha(src/'fixed_panel_TPM_private.tsv.gz')==proof['selected_panel_sha256']
 c=pd.read_excel(src/'Table_S1_Clinical_Annotations.xlsx',header=2);pid='Harmonized_SU2C_Participant_ID_v2';sid='Harmonized_SU2C_RNA_Tumor_Sample_ID_v2'
 assert c[pid].is_unique;cc=c[c.RNA_All.eq(1)&c['Pre-treatment_RNA_Sample_QC'].eq('Keep')].copy().set_index(sid)
 common=expr.columns.intersection(cc.index);assert len(common)==len(cc);expr=expr[common];cc=cc.loc[common]
 cl=pd.DataFrame(index=cc.index);cl['response']=cc.Harmonized_Confirmed_BOR.map({'CR':1,'PR':1,'SD':0,'PD':0})
 for ep in ['OS','PFS']:
  cl[ep]=pd.to_numeric(cc['Harmonized_'+ep+'_Days'],errors='coerce');e=pd.to_numeric(cc['Harmonized_'+ep+'_Event'],errors='coerce');cl[ep+'_event']=e.where(e.isin([0,1]))
 return M.collapse(expr,cl,cc[pid])[:2]

def exact_balance(sc,y,record):
 names=['joint','armA','armB','pan_T'];X=sc[names].to_numpy(float);model,reason=fit(X,y);assert model is not None
 b=model['beta'];se=float(np.sqrt(model['cov'][1,1]));close([np.exp(b[1]),np.exp(b[1]-1.96*se),np.exp(b[1]+1.96*se)],[record.effect,record.CI_low,record.CI_high],1e-7)
 null,why=fit(sc[['armA','armB','pan_T']].to_numpy(float),y);assert null is not None
 P=float(chi2.sf(max(0.,2*(model['ll']-null['ll'])),1));close(P,float(record.p),1e-8)
 A=sc.armA.to_numpy();B=sc.armB.to_numpy();J=sc.joint.to_numpy();mean=(A+B)/2;direction=(A-B)/2;gap=abs(A-B);close(J,mean-gap/2,1e-12)
 other=np.c_[mean,direction,gap,sc.pan_T];mean2=other.mean(0);sd2=other.std(0);Z2=np.c_[np.ones(len(y)),(other-mean2)/sd2];Z=np.c_[np.ones(len(y)),(X-model['mu'])/model['sd']]
 transfer=np.linalg.lstsq(Z2,Z,rcond=None)[0];close(Z2@transfer,Z,1e-10);newbeta=transfer@b;newcov=transfer@model['cov']@transfer.T;close(Z2@newbeta,Z@b,1e-10)
 close(newbeta[3],-b[1]*sd2[2]/(2*model['sd'][0]),1e-10)
 otherfit,why=fit(other,y);assert otherfit is not None;close(otherfit['beta'],newbeta,1e-6);close(otherfit['cov'],newcov,1e-6)
 ss=float(np.sqrt(newcov[3,3]));v=float(newbeta[3])
 return dict(meaning='Conditional absolute arm mismatch, same model not a molecular interaction coefficient',n_units=len(y),mismatch_OR_per_SD=float(np.exp(v)),CI_low=float(np.exp(v-1.96*ss)),CI_high=float(np.exp(v+1.96*ss)),maximum_logit_reparameterization_error=float(np.max(np.abs(Z2@newbeta-Z@b))),original_P_unchanged=P,new_P_generated=False)

def cv_case(expr,clin,panels,arms,cohort,tag,plan):
 target=clin.index[clin.response.isin([0,1])];y=clin.loc[target,'response'].to_numpy(int);n=len(y);minor=min((y==1).sum(),(y==0).sum());base=dict(cohort=cohort,module_proxy=tag,n_units=n,minority=int(minor))
 if n<80 or minor<20:return dict(**base,status='UNESTIMABLE',reason='PREDECLARED_N80_MINORITY20'),[],[],[]
 required=set(panels[arms[0]]+panels[arms[1]]+panels['pan_T'])
 if not required<=set(expr.index):return dict(**base,status='UNESTIMABLE',reason='ORIGINAL_REQUIRED_PANEL_GENE_NOT_MEASURED'),[],[],[]
 models=plan['models'];pred={k:np.full((10,n),np.nan) for k in models};repeat=[];fail=[];proof=[]
 for rep in range(10):
  folds=StratifiedKFold(n_splits=5,shuffle=True,random_state=seed('CM_proxy_CV_20261020|'+cohort+'|'+str(rep)))
  for fold,(tr,te) in enumerate(folds.split(np.zeros(n),y)):
   train_ids=target[tr];test_ids=target[te];pool=expr.columns[~expr.columns.isin(test_ids)]
   assert set(pool).isdisjoint(test_ids) and set(train_ids)<=set(pool)
   score,reason=raw_to_scores(expr,pool,panels,arms)
   if score is None:
    fail.append(dict(repeat=rep,fold=fold,model='ALL',reason=reason));continue
   for name,columns in models.items():
    m,reason=fit(score.loc[train_ids,columns],y[tr])
    if m is None:fail.append(dict(repeat=rep,fold=fold,model=name,reason=reason));continue
    pred[name][rep,te]=predict(m,score.loc[test_ids,columns])
    proof.append(dict(repeat=rep,fold=fold,model=name,train_n=len(tr),test_n=len(te),unlabeled_RNA_in_CDF_pool=len(pool)-len(tr),test_units_excluded_from_CDF=True,max_independent_coefficient_error=m['maximum_coefficient_error']))
  for name,p in pred.items():
   if np.isfinite(p[rep]).all():repeat.append(dict(repeat=rep,model=name,**metrics(y,p[rep])[0]))
 if fail or any(not np.isfinite(p).all() for p in pred.values()):return dict(**base,status='UNESTIMABLE_INCOMPLETE_CV',reason='NO_FAVOURABLE_FOLD_SUBSET',fold_failures=fail),repeat,[],proof
 pooled={k:metrics(y,p.mean(0)) for k,p in pred.items()};comparisons=[]
 for candidate,reference in plan['contrasts']:
  delta=pooled[reference][1]-pooled[candidate][1]
  rng=np.random.default_rng(seed('CM_proxy_loss_boot|'+cohort+'|'+tag+'|'+candidate));take=rng.integers(0,n,size=(2000,n));boot=delta[take].mean(1)
  independent=np.array([np.bincount(v,minlength=n)@delta/n for v in take]);close(boot,independent,1e-12)
  per_repeat=[next(r['logloss'] for r in repeat if r['repeat']==rep and r['model']==reference)-next(r['logloss'] for r in repeat if r['repeat']==rep and r['model']==candidate) for rep in range(10)]
  ci=np.quantile(boot,[.025,.975])
  comparisons.append(dict(cohort=cohort,module_proxy=tag,candidate=candidate,reference=reference,n_units=n,pooled_paired_logloss_improvement=float(delta.mean()),conditional_CI_low=float(ci[0]),conditional_CI_high=float(ci[1]),repeat_improvement_min=float(min(per_repeat)),repeat_improvement_max=float(max(per_repeat)),repeat_improvement_mean=float(np.mean(per_repeat)),all10_repeat_improvements_positive=bool(min(per_repeat)>0),no_new_P_or_FDR=True,interval_scope='Pointwise fixed-OOF-prediction patient bootstrap, not full-training or independent validation uncertainty'))
 return dict(**base,status='COMPLETE_ALL_DECLARED_FOLDS_MODELS',pooled_OOF_metrics={k:v[0] for k,v in pooled.items()},repeats_not_independent_patients=True),repeat,comparisons,proof

def fixture():
 train=np.array([1.,1.,3.,7.]);query=np.array([0.,1.,2.,3.,9.]);close(cdf_values(train,query),[0.,.375,.5,.75,1.],1e-12)
 frozen=cdf_values(train,train);cdf_values(train,query+1e6);assert np.array_equal(cdf_values(train,train),frozen)
 rng=np.random.default_rng(2026102027);x=rng.normal(size=(400,2));y=rng.binomial(1,expit(.9*x[:,0]-.6*x[:,1]));m,r=fit(x,y);assert m is not None and m['beta'][1]>0 and m['beta'][2]<0
 bad,r=fit(np.c_[x[:,0],np.ones(400)],y);assert bad is None and r=='NONFINITE_OR_CONSTANT_REQUIRED_PREDICTOR'
 return dict(status='ACCEPTED_TIED_CDF_HELDOUT_INVARIANCE_AND_INDEPENDENT_FIT_FIXTURE',real_inputs_used=False)

def main():
 cfg=json.load(sys.stdin);plan=cfg['diagnostic_plan'];oldplan=cfg['old_plan'];base=Path(os.environ['PANCANCER_BULK_ROOT']);root=Path(os.environ['PANCANCER_OLD_BULK_RUN']);out=Path(os.environ['PANCANCER_RAM_OUTPUT']);assert out.parent==Path('/dev/shm') and not out.exists();out.mkdir(mode=0o700)
 def save(n,v):(out/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
 save('private_frozen_plan_before_current_CV.json',dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),diagnostic_plan=plan,original_fixed_panel_plan=oldplan,fixture=fixture()))
 p=root/'research_results_v2/all_model_slots_private.tsv';assert sha(p)==plan['original_family_file_SHA256'];prior=pd.read_csv(p,sep='\t',float_precision='round_trip');assert len(prior)==312
 M=load_module(Path(__file__).with_name('explore_priority_immune_CM_bulk_v1.py'));helper=M.helpers(base)
 families=[]
 for family,g in prior.groupby('family'):
  P=g.p.fillna(1).to_numpy();close(M.bh(P),g.BH_q,1e-10);close(M.bh(P,True),g.BY_q,1e-10)
  close(M.bh(P),multipletests(P,method='fdr_bh')[1],1e-12)
  families.append(dict(family=family,registered=len(g),all_BH_BY_equal=True))
 evidence=[]
 fields=['family','cohort','endpoint','model','status','reason','n','events','effect','CI_low','CI_high','p','BH_q','BY_q']
 for name in ['research_results_v2/all_model_slots_private.tsv','diagnostics_v1/ICI_common_program_sensitivity_descriptive_private.tsv','SU2C_clinical_controls_v1/all8_clinical_controls_descriptive_private.tsv']:
  p=root/name;d=pd.read_csv(p,sep='\t');d['module_proxy']=d.axis.map(lambda s:'CM02_proxy' if str(s).startswith('pCM02') else 'CM04_proxy' if str(s).startswith('pCM04') else 'UNKNOWN');assert not d.module_proxy.eq('UNKNOWN').any()
  evidence.append(dict(kind=name.split('/')[0],source_SHA256=sha(p),records=json.loads(d[['module_proxy']+[f for f in fields if f in d]].to_json(orient='records',double_precision=15))))
 save('private_all312_and18_existing_source_evidence.json',evidence)
 results=[];rep_all=[];compare=[];proof_all=[];balance=[];source_proof=[]
 for co in plan['all_ICI_sources']:
  if co=='SU2C_MARK':expr,clin=SU2C(root,M)
  else:expr,clin,meta=M.load_ici(base,co,helper)
  oldfile=root/('SU2C_analysis_v1' if co=='SU2C_MARK' else 'analysis_v2')/(co+'_patient_scores_private.tsv');sc0=read_text_index(oldfile);assert set(sc0.index)==set(clin.index)
  sc0=sc0.reindex(clin.index);close(pd.to_numeric(sc0.response,errors='coerce'),clin.response,1e-12)
  sc,coverage=M.scores(expr,oldplan)
  for col in sc.columns:close(sc[col],pd.to_numeric(sc0[col],errors='coerce'),1e-10)
  source_proof.append(dict(cohort=co,exact_RNA_clinical_units_and_all_old_scores_equal=True,old_scores_SHA256=sha(oldfile),n_RNA_units=len(clin)))
  for tag in plan['module_proxies']:
   ax=next(a for a in oldplan['axes'] if a.startswith('pCM02' if tag=='CM02_proxy' else 'pCM04'));arms=oldplan['axes'][ax]
   result,reps,comparisons,proofs=cv_case(expr,clin,oldplan['panels'],arms,co,tag,plan);results.append(result)
   rep_all.extend([dict(cohort=co,module_proxy=tag,**v) for v in reps]);compare.extend(comparisons);proof_all.extend([dict(cohort=co,module_proxy=tag,**v) for v in proofs])
   expected=prior[prior.cohort.eq(co)&prior.axis.eq(ax)&prior.endpoint.eq('response_increment')];assert len(expected)==1
   if expected.p.notna().iloc[0]:
    ss=pd.DataFrame(dict(joint=sc[ax],armA=sc[arms[0]],armB=sc[arms[1]],pan_T=sc.pan_T,y=clin.response)).dropna();ss=ss[ss.y.isin([0,1])]
    balance.append(dict(cohort=co,module_proxy=tag,**exact_balance(ss,ss.y.to_numpy(int),expected.iloc[0])))
   save('private_progress.json',dict(stage='ALL_FIXED_CM_PROXY_CV',completed_source_module_positions=len(results),registered=12))
  del expr,clin,sc,sc0
 assert len(results)==12 and len(compare)<=24
 save('private_all12_internal_CV_summary.json',results);save('private_all_declared_CV_repeat_metrics.json',rep_all);save('private_all24_performance_contrast_positions.json',[next((r for r in compare if r['cohort']==co and r['module_proxy']==tag and r['candidate']==pair[0]),dict(cohort=co,module_proxy=tag,candidate=pair[0],reference=pair[1],status='UNESTIMABLE_MATCH_PARENT_CONTEXT',no_new_P_or_FDR=True)) for co in plan['all_ICI_sources'] for tag in plan['module_proxies'] for pair in plan['contrasts']])
 save('private_all_existing_balance_reparameterizations.json',balance)
 save('private_all_fold_and_source_numeric_proofs.json',dict(source_proofs=source_proof,fold_proofs=proof_all,old_full_family_qvalues=families,no_patient_ID_or_individual_predictions_saved=True))
 files={p.name:dict(bytes=p.stat().st_size,SHA256=sha(p)) for p in out.iterdir() if p.is_file()};assert sum(x['bytes'] for x in files.values())<4*1024**2
 save('private_terminal.json',dict(status='COMPLETE_EXISTING312_EVIDENCE_ALL12_CV_ALL24_CONTRASTS_EXACT_BALANCE_NUMERIC_REPLAY_PENDING_SCIENTIFIC_REVIEW',completed_UTC=datetime.now(timezone.utc).isoformat(),source_SHA256=sha(__file__),files=files,new_P_or_FDR_family=False,all_real_expression_patient_scores_fold_predictions_designs_and_bootstrap_arrays_RAM_only=True,no_new_RNA_or_denied_clinical_CSV_download=True,goal_status='active'))
 print(json.dumps(dict(status='COMPLETE',terminal_SHA256=sha(out/'private_terminal.json'))))
if __name__=='__main__':main()
