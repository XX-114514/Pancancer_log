#!/usr/bin/env python3
"""Fixed external NSCLC proxy analysis; public author RNA, exact clinical/QC matching."""
import argparse,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
from explore_priority_immune_CM_bulk_v1 import scores,collapse,fit_model,bh

def main(a):
 assert not a.output.exists();a.output.mkdir(mode=0o700)
 plan=json.loads(a.plan.read_text());external=json.loads(a.external_plan.read_text())
 assert external['no_new_outcomes_read_yet'] and external['axes_identical_to_original_frozen_plan']
 expr=pd.read_csv(a.source/'fixed_panel_TPM_private.tsv.gz',sep='\t',index_col=0)
 receipt=json.loads((a.source/'panel_derivation_receipt_public.json').read_text())
 assert hashlib.sha256((a.source/'fixed_panel_TPM_private.tsv.gz').read_bytes()).hexdigest()==receipt['selected_panel_sha256']
 c=pd.read_excel(a.source/'Table_S1_Clinical_Annotations.xlsx',header=2)
 pid='Harmonized_SU2C_Participant_ID_v2';sid='Harmonized_SU2C_RNA_Tumor_Sample_ID_v2'
 assert c[pid].is_unique and c.loc[c.RNA_All.eq(1),sid].is_unique
 cc=c[c.RNA_All.eq(1)&c['Pre-treatment_RNA_Sample_QC'].eq('Keep')].copy().set_index(sid)
 common=expr.columns.intersection(cc.index);assert len(common)==len(cc),'All author Keep samples require exact RNA key'
 cc=cc.loc[common];expr=expr[common]
 clin=pd.DataFrame(index=cc.index);clin['response']=cc.Harmonized_Confirmed_BOR.map({'CR':1,'PR':1,'SD':0,'PD':0})
 for ep in ['OS','PFS']:
  clin[ep]=pd.to_numeric(cc['Harmonized_'+ep+'_Days'],errors='coerce');event=pd.to_numeric(cc['Harmonized_'+ep+'_Event'],errors='coerce');clin[ep+'_event']=event.where(event.isin([0,1]))
 expr,clin,merge_audit=collapse(expr,clin,cc[pid]);sc,coverage=scores(expr,plan)
 sc.join(clin).to_csv(a.output/'SU2C_MARK_patient_scores_private.tsv',sep='\t')
 rows=[];designs=[]
 for axis,arms in plan['axes'].items():
  d=pd.DataFrame(dict(joint=sc[axis],pan_T=sc.pan_T,armA=sc[arms[0]],armB=sc[arms[1]],y=clin.response))
  for ep in ['response','OS','PFS','response_increment']:
   family='SU2C_increment' if ep=='response_increment' else 'SU2C_response' if ep=='response' else 'SU2C_survival'
   kind='logistic' if ep.startswith('response') else 'Cox';d['t']=clin.get(ep,pd.Series(np.nan,index=clin.index));d['e']=clin.get(ep+'_event',pd.Series(np.nan,index=clin.index))
   predictors=['joint','pan_T'] if ep!='response_increment' else ['joint','armA','armB','pan_T']
   if ep=='response_increment' and (len(d.dropna(subset=['y','joint','pan_T','armA','armB']))<80 or min(clin.response.eq(1).sum(),clin.response.eq(0).sum())<20):
    fit=dict(status='unavailable',n=len(clin),events=int(clin.response.eq(1).sum()),reason='increment_min80_minority20_gate');used=None
   else:fit,used=fit_model(d,kind,predictors)
   rows.append(dict(family=family,cohort='SU2C_MARK',axis=axis,endpoint=ep,model=kind,**fit))
   if used is not None:designs.append(dict(record_index=len(rows)-1,kind=kind,predictors=predictors,data=used.to_dict(orient='list')))
 result=pd.DataFrame(rows);assert len(result)==8
 for family,n in [('SU2C_response',2),('SU2C_survival',4),('SU2C_increment',2)]:
  ix=result.index[result.family.eq(family)];assert len(ix)==n;p=result.loc[ix,'p'].fillna(1).to_numpy();result.loc[ix,'BH_q']=bh(p);result.loc[ix,'BY_q']=bh(p,True)
 result.to_csv(a.output/'all_frozen_model_slots_private.tsv',sep='\t',index=False)
 (a.output/'fitted_designs_private.json').write_text(json.dumps(designs,indent=2))
 qc=dict(n_full_author_patients=len(c),n_RNA_source_columns=receipt['n_source_RNA_samples'],n_author_RNA_Keep=int(c['Pre-treatment_RNA_Sample_QC'].eq('Keep').sum()),n_author_RNA_Flag_excluded=int(c['Pre-treatment_RNA_Sample_QC'].eq('Flag').sum()),exact_RNA_clinical_mapping=True,coverage=coverage,merge_audit=merge_audit,expression_full_source_retained_local=True,server_analysis_uses_exact_frozen_gene_subset=True,no_assumption_393_RNA_patients=True)
 (a.output/'cohort_acceptance_public.json').write_text(json.dumps(qc,indent=2))
 cols=[pid,'Patient_Age_at_Diagnosis','Patient_Sex','Histology_Harmonized','Line_of_Therapy','Agent_PD1_Category']
 cc[cols].set_index(pid).to_csv(a.output/'clinical_descriptive_covariates_private.tsv',sep='\t')
 terminal=dict(status='ALL8_FROZEN_EXTERNAL_SLOTS_COMPLETE_INDEPENDENT_AUDIT_PENDING',finished_utc=datetime.now(timezone.utc).isoformat(),n_slots=8,n_fitted=int(result.status.eq('exploratory_model_numeric_audit_pass').sum()),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),external_plan_sha256=hashlib.sha256(a.external_plan.read_bytes()).hexdigest(),original_panels_and_thresholds_unchanged=True)
 (a.output/'terminal_public.json').write_text(json.dumps(terminal,indent=2));print(json.dumps(terminal),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser()
 for k in ['source','plan','external-plan','output']:p.add_argument('--'+k,type=Path,required=True)
 main(p.parse_args())
