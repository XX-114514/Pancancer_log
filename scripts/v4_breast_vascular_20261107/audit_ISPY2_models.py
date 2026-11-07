"""Independent native-key/rank audit and NumPy likelihood refit of R glm models."""
from pathlib import Path
import os,json,gzip,re,collections,math,hashlib
from datetime import datetime,timezone
import numpy as np,pandas as pd
w=Path(os.environ['BREAST_WORK']);r=w/'clinical_v2';src=Path(os.environ['ISPY2_SOURCE']);signature=Path(os.environ['V4_SIGNATURES'])
meta=[]
for block in re.split(r'\^SAMPLE = ',(src/'GSE194040_sample_metadata_PRIVATE.txt').read_text())[1:]:
 d={}
 for line in block.splitlines():
  if line.startswith('!Sample_characteristics_ch1 = '):
   k,_,v=line.split(' = ',1)[1].partition(': ');d[k]=v
 meta.append(d)
ct=collections.Counter(x['patient id'] for x in meta);chosen=[d for d in meta if ct[d['patient id']]==1 and d['her2']=='0' and d['hr'] in ['0','1'] and d['pcr'] in ['0','1'] and d['arm'] in ['Paclitaxel','Paclitaxel + Pembrolizumab']];keys=[d['patient id'] for d in chosen]
sig=pd.read_csv(signature,sep='\t');panels={cm:sig.loc[sig.cm.eq(cm),'gene'].tolist() for cm in ['CM03','CM04','CM07','CM09','CM11']};panels.update(B=['IGHD','MS4A1','CD79A'],EC=['PECAM1','VWF','CDH5'],NK=['XCL1','NCR1','KLRD1'],panT=['CD3D','CD3E'])
need=set(sum(panels.values(),[]));values={};occ=collections.Counter()
with gzip.open(src/'GSE194040_geneLevel_PRIVATE.txt.gz','rt') as f:
 header=f.readline().rstrip('\r\n').split('\t');idx=[header.index(k)+1 for k in keys];assert len(header)==len(set(header))
 for line in f:
  fields=line.rstrip('\r\n').split('\t');gene=fields[0]
  if gene in need:
   occ[gene]+=1;values[gene]=np.array([float(fields[i]) if fields[i] not in ['NA','NaN','nan',''] else np.nan for i in idx])
def ranks(v):
 result=np.full(len(v),np.nan);ii=np.where(np.isfinite(v))[0];ss=ii[np.argsort(v[ii],kind='stable')];j=0
 while j<len(ss):
  end=j+1
  while end<len(ss) and v[ss[end]]==v[ss[j]]:end+=1
  result[ss[j:end]]=(j+1+end)/2/len(ss);j=end
 return result
scores={}
for name,genes in panels.items():
 if not all(occ[g]==1 for g in genes):scores[name]=np.full(len(keys),np.nan);continue
 v=np.column_stack([values[g] for g in genes]);scores[name]=np.column_stack([ranks(values[g]) for g in genes]).mean(axis=1);scores[name][~np.isfinite(v).all(axis=1)]=np.nan
s=np.column_stack([scores[k] for k in ['B','EC','NK']]);scores['Vascular_B_NK_equal_mean']=s.mean(axis=1);scores['Vascular_B_NK_min_component']=s.min(axis=1);scores['mean_B_NK']=np.column_stack([scores['B'],scores['NK']]).mean(axis=1)
raw=pd.read_csv(r/'unscaled_scores_PRIVATE.csv',index_col=0,float_precision='round_trip');raw.index=raw.index.astype(str);assert list(raw.index)==keys
max_measurement=0
for k,v in scores.items():
 a=raw[k].to_numpy();assert np.array_equal(np.isnan(a),np.isnan(v));max_measurement=max(max_measurement,float(np.max(np.abs(a[np.isfinite(v)]-v[np.isfinite(v)]))) if np.isfinite(v).any() else 0)
scaled={}
for k,v in scores.items():
 scaled[k]=(v-np.nanmean(v))/np.nanstd(v,ddof=1) if np.isfinite(v).any() else v
scaled.update(HR=np.array([int(d['hr']) for d in chosen]),pCR=np.array([int(d['pcr']) for d in chosen]),treatment=np.array([int(d['arm']=='Paclitaxel + Pembrolizumab') for d in chosen]))
lookup={k:i for i,k in enumerate(keys)};t=pd.read_csv(r/'all24_prespecified_target_terms.tsv',sep='\t',float_precision='round_trip');max_numeric=0;max_gradient=0;n_design=0
for f in sorted((r/'models_PRIVATE').glob('*.csv')):
 a=pd.read_csv(f,dtype={'patient_key':str},float_precision='round_trip');z=json.loads(f.with_suffix('.json').read_text());family,score,scope=f.stem.split('__');row=t[(t.family==family)&(t.score==score)&(t.target==scope)];assert len(row)==1;row=row.iloc[0]
 expected=np.isfinite(scaled[score])&np.isfinite(scaled['panT'])
 if family=='secondary_incremental3':expected&=np.isfinite(scaled['EC'])&np.isfinite(scaled['mean_B_NK'])
 if scope!='interaction':expected&=scaled['treatment']==(scope=='pembro')
 assert a.patient_key.tolist()==[k for i,k in enumerate(keys) if expected[i]]
 ids=np.array([lookup[k] for k in a.patient_key]);y=scaled['pCR'][ids];assert np.array_equal(y,a.outcome)
 cols=[]
 for name in a.columns[2:]:
  if name=='Intercept':v=np.ones(len(ids))
  elif name.endswith('_treatment'):
   base=name[:-len('_treatment')];v=scaled[score if base=='score' else base][ids]*scaled['treatment'][ids]
  else:v=scaled[score if name=='score' else name][ids]
  assert np.allclose(a[name],v,atol=2e-14,rtol=2e-14),name;cols.append(v)
 x=np.column_stack(cols);b=np.zeros(x.shape[1])
 def loss(bb):eta=x@bb;return float(np.sum(np.logaddexp(0,eta)-y*eta))
 for iteration in range(150):
  eta=x@b;p=1/(1+np.exp(-eta));grad=x.T@(p-y);h=x.T@(x*(p*(1-p))[:,None]);step=np.linalg.solve(h,grad)
  if np.max(np.abs(grad))<1e-10:break
  alpha=1.;old=loss(b)
  while loss(b-alpha*step)>old+1e-12 and alpha>1e-10:alpha/=2
  b=b-alpha*step
 assert iteration<149
 p=1/(1+np.exp(-(x@b)));cov=np.linalg.inv(x.T@(x*(p*(1-p))[:,None]));j=z['target_index']-1;se=np.sqrt(cov[j,j]);pval=math.erfc(abs(b[j]/se)/math.sqrt(2))
 errors=[np.max(np.abs(b-np.array(z['beta']))),np.max(np.abs(cov-np.array(z['covariance']))),abs(b[j]-row.beta),abs(se-row.SE),abs(pval-row.P),abs(math.exp(b[j])-row.OR),abs(math.exp(b[j]-1.959963984540054*se)-row.CI_low),abs(math.exp(b[j]+1.959963984540054*se)-row.CI_high)]
 max_numeric=max(max_numeric,*[float(v) for v in errors]);max_gradient=max(max_gradient,float(np.max(np.abs(x.T@(p-y)))));assert max(errors)<1e-5,(f.name,errors);n_design+=1
# Independent full-family FDR, including every NA as P1.
def adjust(p,by=False):
 order=np.argsort(p);rank=np.arange(1,len(p)+1);q=p[order]*len(p)/rank
 if by:q*=np.sum(1/rank)
 q=np.minimum(1,np.minimum.accumulate(q[::-1])[::-1]);out=np.empty(len(p));out[order]=q;return out
for family,n in [('primary21',21),('secondary_incremental3',3)]:
 z=t[t.family==family];assert len(z)==n;p=z.P.fillna(1).to_numpy();valid=z.P.notna().to_numpy()
 for label,by in [('BH',False),('BY',True)]:assert np.allclose(z[label].to_numpy()[valid],adjust(p,by)[valid],atol=1e-13)
assert max_measurement<1e-12
receipt={'UTC':datetime.now(timezone.utc).isoformat(),'status':'PASS_INDEPENDENT_RAW_KEYS_RANKS_DESIGNS_NUMPY_LIKELIHOOD_AND_FULL_FDR','n_eligible_patients':len(keys),'n_designs_refitted':n_design,'n_primary_slots':21,'n_secondary_slots':3,'max_native_rank_score_error':max_measurement,'max_coefficient_covariance_OR_CI_P_error':max_numeric,'max_abs_gradient':max_gradient,'plan_sha256':hashlib.sha256((w/'frozen_breast_plan_v1.json').read_bytes()).hexdigest(),'measurement_addendum_sha256':hashlib.sha256((w/'measurement_addendum_before_estimable_effects_v2.json').read_bytes()).hexdigest(),'source_patient_design_and_expression_not_public':True}
(r/'independent_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt))
