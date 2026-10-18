#!/usr/bin/env python3
"""All-fitted-model PH warnings and prespecified descriptive common-program controls."""
import argparse,json,warnings
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from statsmodels.duration.hazard_regression import PHReg
from explore_priority_immune_CM_bulk_v1 import fit_model

def main(a):
 proof=json.loads((a.audit/'independent_audit_public.json').read_text());assert proof['full_family_BH_BY_pass']
 assert not a.output.exists();a.output.mkdir(mode=0o700)
 results=pd.read_csv(a.run/'all_frozen_model_slots_private.tsv',sep='\t');designs=json.loads((a.run/'fitted_designs_private.json').read_text());ph=[]
 for entry in designs:
  if entry['kind']!='Cox':continue
  rec=results.iloc[entry['record_index']];d=pd.DataFrame(entry['data']);X=d[entry['predictors']].to_numpy();X=(X-X.mean(0))/X.std(0)
  with warnings.catch_warnings():
   warnings.simplefilter('ignore');fit=PHReg(d.t.to_numpy(),X,status=d.e.to_numpy(),ties='efron').fit(disp=False)
  scaled=np.asarray(fit.schoenfeld_residuals)@fit.cov_params();mask=d.e.eq(1).to_numpy()&np.isfinite(scaled[:,0]);rho,p=spearmanr(scaled[mask,0],np.log(d.t.to_numpy()[mask]))
  ph.append(dict(record_index=int(entry['record_index']),cohort=rec.cohort,axis=rec.axis,endpoint=rec.endpoint,model=rec.model,n_events=int(mask.sum()),diagnostic_spearman_rho=float(rho),diagnostic_p=float(p),interpretation='Descriptive scaled-Schoenfeld/time correlation warning, not proof of PH or a new clinical endpoint.'))
 pd.DataFrame(ph).to_csv(a.output/'all_Cox_PH_diagnostic_private.tsv',sep='\t',index=False)
 rows=[]
 if a.mode=='existing':
  for co in ['IMvigor210','Riaz2017','Hugo2016','Jung2019','Rose2021']:
   c=pd.read_csv(a.run/(co+'_patient_scores_private.tsv'),sep='\t',index_col=0)
   for axis in ['pCM02_regulatory_cytotoxic','pCM04_memory_GZMK_effector']:
    d=pd.DataFrame(dict(y=c.response,joint=c[axis],pan_T=c.pan_T,interferon=c.interferon,proliferation=c.proliferation))
    nminor=min(d.y.eq(1).sum(),d.y.eq(0).sum())
    if nminor<20:rec=dict(status='unavailable',reason='min20_minority_for4predictor_descriptive_model')
    else:rec,_=fit_model(d,'logistic',['joint','pan_T','interferon','proliferation'])
    for k in ['p','LR','null_loglik']:rec.pop(k,None)
    rows.append(dict(cohort=co,axis=axis,**rec))
  pd.DataFrame(rows).to_csv(a.output/'ICI_common_program_sensitivity_descriptive_private.tsv',sep='\t',index=False)
 result=dict(status='ALL_FITTED_COX_PH_DIAGNOSTICS_AND_DECLARED_DESCRIPTIVE_CONTROLS_COMPLETE',finished_utc=datetime.now(timezone.utc).isoformat(),n_Cox_diagnostics=len(ph),n_descriptive_controls=len(rows),no_new_significance_family=True,primary_results_unchanged=True)
 (a.output/'diagnostics_public.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))

if __name__=='__main__':
 p=argparse.ArgumentParser()
 for name in ['run','audit','output']:p.add_argument('--'+name,type=Path,required=True)
 p.add_argument('--mode',choices=['existing','external'],default='existing');main(p.parse_args())
