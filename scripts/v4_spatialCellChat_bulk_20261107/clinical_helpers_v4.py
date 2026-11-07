"""Reused audited clinical readers and model arithmetic; no old CM signatures."""
import ast,gzip,json,hashlib,re,warnings
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.special import expit
from scipy.stats import chi2,norm
from statsmodels.duration.hazard_regression import PHReg
import statsmodels.api as sm

def now():return datetime.now(timezone.utc).isoformat()

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1024**2),b''):h.update(b)
 return h.hexdigest()

def event(phase,**kw):print(json.dumps(dict(utc=now(),phase=phase,**kw)),flush=True)

def writej(p,x):
 with Path(p).open('x') as f:json.dump(x,f,indent=2,allow_nan=False)

def bh(p,by=False):
 p=np.asarray(p,float);n=len(p);o=np.argsort(p);factor=np.sum(1/np.arange(1,n+1)) if by else 1
 q=np.minimum.accumulate((p[o]*n*factor/np.arange(1,n+1))[::-1])[::-1];out=np.empty(n);out[o]=np.clip(q,0,1);return out

def numeric(s):return pd.to_numeric(s,errors='coerce')

def binary(s,mapping):return s.astype(str).str.strip().str.lower().map(mapping)

def stage(s):
 def one(x):
  x=str(x).strip().upper().replace('STAGE','').replace(' ','')
  z=re.fullmatch('(IV|III|II|I)[ABC]?',x)
  return {'I':1,'II':2,'III':3,'IV':4}.get(z.group(1),np.nan) if z else np.nan
 return s.map(one)

def helpers(base):
 source=base/'ici_bulk_analysis.py';tree=ast.parse(source.read_text())
 names={'parse_series_matrix','load_riaz','load_jung','load_imvigor210'}
 body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names]
 scope=dict(BASE=base,ICI=base/'ici_bulk',Path=Path,pd=pd,np=np,gzip=gzip)
 exec(compile(ast.Module(body=body,type_ignores=[]),str(source),'exec'),scope);return scope

def collapse(expr,clin,patient):
 assert expr.columns.is_unique and clin.index.is_unique and set(expr.columns)<=set(clin.index)
 clin=clin.reindex(expr.columns).copy();clin['patient']=patient.reindex(clin.index)
 assert clin.patient.notna().all()
 conflicts=clin.groupby('patient')[['response','OS','OS_event','PFS','PFS_event']].nunique(dropna=True).gt(1).any(axis=1)
 bad=set(conflicts.index[conflicts]);good=~clin.patient.isin(bad);clin=clin.loc[good];expr=expr.loc[:,clin.index]
 out=expr.T.groupby(clin.patient).mean().T
 c=clin.groupby('patient').first().reindex(out.columns);assert c.index.is_unique
 return out,c,dict(conflicting_patients_excluded=len(bad),eligible_samples=len(clin),unique_patients=len(c))

def load_ici(base,name,h):
 root=base/'ici_bulk';meta={};time_status='pretreatment_source_design'
 if name=='Riaz2017':
  expr,c,pre,_=h['load_riaz']();c=c.loc[pre];expr=expr.loc[:,c.index]
  cc=pd.DataFrame(index=c.index);cc['response']=c.response.map({'PRCR':1,'PD':0,'SD':0})
  pt=pd.Series(c.index,index=c.index).str.extract(r'^(Pt\d+)_',expand=False)
 elif name=='Hugo2016':
  raw=pd.read_excel(root/'hugo2016/GSE78220_FPKM.xlsx');raw=raw.set_index(raw.columns[0]);raw=raw.apply(pd.to_numeric,errors='coerce')
  c=h['parse_series_matrix'](root/'hugo2016/GSE78220_series_matrix.txt.gz');c=c[c['biopsy time'].str.lower().eq('pre-treatment')].copy()
  rows=[]
  for col in raw.columns:
   key=str(col).split('.')[0];hit=c[c.title.eq(key)]
   if len(hit)==0:hit=c[c['patient id'].eq(key)]
   if len(hit)==1:rows.append((col,hit.iloc[0]))
   elif len(hit)>1:raise ValueError('Ambiguous Hugo expression/metadata mapping')
  c=pd.DataFrame([x[1] for x in rows],index=[x[0] for x in rows]);expr=np.log2(raw.loc[:,c.index].clip(lower=0)+1)
  cc=pd.DataFrame(index=c.index);cc['response']=c['anti-pd-1 response'].map({'Complete Response':1,'Partial Response':1,'Progressive Disease':0,'Stable Disease':0})
  cc['OS']=numeric(c['overall survival (days)']);cc['OS_event']=binary(c['vital status'],{'dead':1,'alive':0});pt=c['patient id']
 elif name=='Jung2019':
  expr,c,t,e=h['load_jung']();cc=pd.DataFrame(index=c.index);cc['PFS']=t;cc['PFS_event']=e.where(e.isin([0,1]));pt=pd.Series(c.index,index=c.index)
 elif name=='IMvigor210':
  expr,c,t,e,r=h['load_imvigor210']();cc=pd.DataFrame(index=c.index);cc['response']=r.map({'R':1,'NR':0})
  cc['OS']=t;cc['OS_event']=binary(c.OS_STATUS,{'1:deceased':1,'0:living':0});pt=pd.Series(c.index,index=c.index)
  assert c.ICI_RX.eq('Atezolizumab').all()
 else:
  expr=pd.read_csv(root/'rose2021_urothelial/GSE176307_salmon_tpm.tsv.gz',sep='\t',index_col=0);expr=expr[~expr.index.duplicated()];expr=np.log2(expr.clip(lower=0)+1)
  key=pd.read_csv(root/'rose2021_urothelial/sample_key.csv.gz');assert key.shape[1]==2;key.columns=['sample_id','rs_id']
  c=h['parse_series_matrix'](root/'rose2021_urothelial/series_matrix.txt.gz');c['sample_id']=c.title.str.replace('Patient sample ','',regex=False)
  c=c.merge(key,on='sample_id',validate='one_to_one').set_index('rs_id')
  common=expr.columns.intersection(c.index);unmapped=len(expr.columns)-len(common);expr=expr[common];c=c.loc[common]
  assert c.sample_id.notna().all();cc=pd.DataFrame(index=c.index);pt=c.sample_id
  cc['response']=c['io.response'].str.upper().map({'CR':1,'PR':1,'SD':0,'PD':0})
  cc['OS']=numeric(c['overall survival']);cc['OS_event']=binary(c.alive,{'no':1,'yes':0})
  cc['PFS']=numeric(c.pfs);cc['PFS_event']=binary(c.progressed,{'yes':1,'no':0,'1':1,'0':0})
  time_status='archival_tumor_individual_sampling_time_not_verified_no_baseline_prediction_claim'
 for col in ['response','OS','OS_event','PFS','PFS_event']:
  if col not in cc:cc[col]=np.nan
 expr,c,meta=collapse(expr,cc,pt)
 meta.update(cohort=name,time_status=time_status,age_sex_adjustment_available=False)
 if name=='Rose2021':meta['expression_samples_without_exact_clinical_key_excluded']=unmapped
 return expr,c,meta

def load_tcga(base,genes):
 r=base/'tcga';pm=pd.read_csv(r/'gencode.v23.annotation.gene.probemap',sep='\t',dtype=str)
 pm['id']=pm.id.str.split('.').str[0];idmap=pm.drop_duplicates('id').set_index('id').gene.to_dict()
 duplicate=pm.groupby('gene').id.nunique();ambiguous=set(duplicate.index[duplicate>1]);rows=[];names=[]
 with gzip.open(r/'tcga_RSEM_gene_tpm.gz','rt') as f:
  header=f.readline().rstrip('\n').split('\t')[1:]
  for line in f:
   ident,_,values=line.rstrip('\n').partition('\t');g=idmap.get(ident.split('.')[0])
   if g in genes and g not in ambiguous:
    rows.append(np.fromstring(values,sep='\t'));names.append(g)
 assert len(names)==len(set(names))
 expr=pd.DataFrame(rows,index=names,columns=header);cols=[c for c in expr if re.match(r'^TCGA-[^-]+-[^-]+-01',c)]
 expr=expr[cols];expr=expr.T.groupby([c[:12] for c in cols]).mean().T
 c=pd.read_csv(r/'TCGA_survival_S1.tsv',sep='\t',dtype=str).rename(columns={'_PATIENT':'patient','cancer type abbreviation':'cancer','age_at_initial_pathologic_diagnosis':'age'})
 c.patient=c.patient.str[:12]
 for k in ['age','OS','OS.time','PFI','PFI.time']:c[k]=numeric(c[k])
 c['sex']=binary(c.gender,{'male':1,'female':0});c['stage']=stage(c.ajcc_pathologic_tumor_stage).fillna(stage(c.clinical_stage))
 required=['cancer','age','sex','stage','OS','OS.time','PFI','PFI.time']
 bad=c.groupby('patient')[required].nunique(dropna=True).gt(1).any(axis=1);bad=set(bad.index[bad])
 c=c[~c.patient.isin(bad)].groupby('patient').first();common=c.index.intersection(expr.columns);c=c.loc[common];expr=expr[common]
 return expr,c[required],dict(cohort='TCGA',unique_patients=len(c),clinical_conflicting_patients_excluded=len(bad),sample_type='01 primary solid tumour',ambiguous_requested_symbols=sorted(genes&ambiguous))

def independent_cox(t,e,X,beta):
 order=np.argsort(t,kind='stable');t=t[order];e=e[order];X=X[order];eta=X@beta;offset=eta.max();w=np.exp(eta-offset)
 r0=np.cumsum(w[::-1])[::-1];r1=np.cumsum((w[:,None]*X)[::-1],axis=0)[::-1]
 r2=np.cumsum((w[:,None,None]*X[:,:,None]*X[:,None,:])[::-1],axis=0)[::-1]
 ll=float(eta[e==1].sum());g=X[e==1].sum(0);info=np.zeros((X.shape[1],X.shape[1]))
 for time in np.unique(t[e==1]):
  deaths=np.flatnonzero((t==time)&(e==1));j=np.searchsorted(t,time);d=len(deaths);f=np.arange(d)/d
  d0=w[deaths].sum();d1=(w[deaths,None]*X[deaths]).sum(0);d2=(w[deaths,None,None]*X[deaths,:,None]*X[deaths,None,:]).sum(0)
  den=r0[j]-f*d0;mu=(r1[j]-f[:,None]*d1)/den[:,None]
  ll-=float((np.log(den)+offset).sum());g-=mu.sum(0)
  info+=((r2[j]-f[:,None,None]*d2)/den[:,None,None]-mu[:,:,None]*mu[:,None,:]).sum(0)
 return ll,g,info

def fit_model(data,kind,predictors,tcga=False):
 cols=['t','e'] if kind=='Cox' else ['y'];d=data[cols+predictors].replace([np.inf,-np.inf],np.nan).dropna()
 if kind=='Cox':d=d[(d.t>0)&d.e.isin([0,1])]
 else:d=d[d.y.isin([0,1])]
 result=dict(status='unavailable',n=len(d),events=int(d.e.sum()) if kind=='Cox' else int(d.y.sum()),reason='')
 n=len(d);p=len(predictors)
 gate=n>=60 and result['events']>=max(20,10*p) if tcga else n>=30 and (result['events']>=max(10,5*p) if kind=='Cox' else min(result['events'],n-result['events'])>=10)
 if not gate:result['reason']='insufficient_patients_or_events';return result,None
 X=d[predictors].to_numpy(float);sd=X.std(axis=0)
 if (sd<1e-10).any():result['reason']='constant_required_predictor';return result,None
 X=(X-X.mean(0))/sd;design=X if kind=='Cox' else sm.add_constant(X)
 if np.linalg.matrix_rank(design)!=design.shape[1] or np.linalg.cond(design)>1e6:
  result['reason']='unidentifiable_design';return result,None
 with warnings.catch_warnings(record=True) as caught:
  warnings.simplefilter('always')
  try:
   if kind=='Cox':
    model=PHReg(d.t.to_numpy(),X,status=d.e.to_numpy(),ties='efron');fit=model.fit(disp=False)
    null=PHReg(d.t.to_numpy(),X[:,1:],status=d.e.to_numpy(),ties='efron').fit(disp=False)
    b=fit.params[0];se=fit.bse[0];ll,g,info=independent_cox(d.t.to_numpy(),d.e.to_numpy(),X,fit.params)
    assert abs(ll-fit.llf)<1e-5 and np.max(np.abs(g))<1e-4 and np.max(np.abs(info+model.hessian(fit.params)))<1e-4
    assert np.allclose(np.linalg.inv(info),fit.cov_params(),atol=1e-5,rtol=1e-5)
    score_index=0
   else:
    model=sm.GLM(d.y.to_numpy(),design,family=sm.families.Binomial());fit=model.fit(maxiter=100)
    null=sm.GLM(d.y.to_numpy(),design[:,[0]+list(range(2,design.shape[1]))],family=sm.families.Binomial()).fit(maxiter=100)
    assert fit.converged and null.converged
    b=fit.params[1];se=fit.bse[1];eta=design@fit.params;mu=expit(eta)
    ll=float((d.y.to_numpy()*eta-np.logaddexp(0,eta)).sum());g=design.T@(d.y.to_numpy()-mu);info=design.T@((mu*(1-mu))[:,None]*design)
    assert abs(ll-fit.llf)<1e-6 and np.max(np.abs(g))<1e-5 and np.allclose(np.linalg.inv(info),fit.cov_params(),atol=1e-5,rtol=1e-5)
    if np.max(np.abs(fit.params))>15:raise ValueError('separation_or_extreme_coefficients')
    score_index=1
   assert np.isfinite(fit.params).all() and np.isfinite(fit.bse).all() and se>0
   lr=max(0,float(2*(fit.llf-null.llf)));pv=float(chi2.sf(lr,1))
   result.update(status='exploratory_model_numeric_audit_pass',reason='',log_effect=float(b),se=float(se),effect=float(np.exp(b)),CI_low=float(np.exp(b-1.96*se)),CI_high=float(np.exp(b+1.96*se)),p=pv,loglik=float(fit.llf),null_loglik=float(null.llf),LR=lr,maximum_independent_score_residual=float(np.max(np.abs(g))),predictors=';'.join(predictors),warnings=';'.join(sorted({str(x.message) for x in caught})))
   return result,d
  except Exception as ex:
   if isinstance(ex,AssertionError):raise
   result['reason']=type(ex).__name__+':'+str(ex)[:150];return result,None
