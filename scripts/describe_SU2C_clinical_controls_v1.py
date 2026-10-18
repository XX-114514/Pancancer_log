#!/usr/bin/env python3
"""Complete descriptive age/sex/line/combination controls; no new significance family."""
import argparse,json
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from explore_priority_immune_CM_bulk_v1 import fit_model

def main(a):
 assert json.loads((a.audit/'independent_audit_public.json').read_text())['n_slots']==8
 assert not a.output.exists();a.output.mkdir(mode=0o700)
 c=pd.read_csv(a.run/'SU2C_MARK_patient_scores_private.tsv',sep='\t',index_col=0)
 v=pd.read_csv(a.run/'clinical_descriptive_covariates_private.tsv',sep='\t',index_col=0).reindex(c.index)
 age=pd.to_numeric(v.Patient_Age_at_Diagnosis,errors='coerce');sex=v.Patient_Sex.map({'M':1,'F':0})
 line=pd.to_numeric(v.Line_of_Therapy,errors='coerce');combo=v.Agent_PD1_Category.map({'PD(L)1':0,'PD(L)1 + CTLA4':1,'PD(L)1 + Other':1})
 rows=[]
 for ax in ['pCM02_regulatory_cytotoxic','pCM04_memory_GZMK_effector']:
  armA,armB=('regulatory_associated','cytotoxic') if ax.startswith('pCM02') else ('memory_associated','GZMK_effector')
  d=pd.DataFrame(dict(joint=c[ax],armA=c[armA],armB=c[armB],pan_T=c.pan_T,age=age,sex=sex,line=line,combination=combo,y=c.response))
  for ep in ['response','OS','PFS','response_increment']:
   d['t']=c.get(ep,pd.Series(np.nan,index=c.index));d['e']=c.get(ep+'_event',pd.Series(np.nan,index=c.index));kind='logistic' if ep.startswith('response') else 'Cox'
   predictors=['joint','pan_T','age','sex','line','combination']+(['armA','armB'] if ep=='response_increment' else [])
   used=d.dropna(subset=predictors+(['y'] if kind=='logistic' else ['t','e']))
   minority=min(used.y.eq(1).sum(),used.y.eq(0).sum()) if kind=='logistic' else used.e.eq(1).sum()
   if len(used)<(100 if ep=='response_increment' else 60) or minority<(40 if ep=='response_increment' else 30):fit=dict(status='unavailable',n=len(used),reason='prespecified_descriptive_clinical_design_size_gate')
   else:fit,_=fit_model(d,kind,predictors)
   for k in ['p','LR','null_loglik']:fit.pop(k,None)
   rows.append(dict(cohort='SU2C_MARK',axis=ax,endpoint=ep,interpretation='Descriptive clinical adjustment; not a replacement primary significance test',**fit))
 pd.DataFrame(rows).to_csv(a.output/'all8_clinical_controls_descriptive_private.tsv',sep='\t',index=False)
 terminal=dict(status='ALL8_DESCRIPTIVE_CLINICAL_CONTROLS_COMPLETE',finished_utc=datetime.now(timezone.utc).isoformat(),n_descriptive_slots=8,no_new_significance_family=True,all_primary_8_slots_unchanged=True)
 (a.output/'terminal_public.json').write_text(json.dumps(terminal,indent=2));print(json.dumps(terminal))

if __name__=='__main__':
 p=argparse.ArgumentParser()
 for name in ['run','audit','output']:p.add_argument('--'+name,type=Path,required=True)
 main(p.parse_args())
