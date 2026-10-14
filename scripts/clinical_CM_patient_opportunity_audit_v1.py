"""Independent patient and capture opportunity audit; no new outcome tests.
Only aggregate outputs may be public. Patient and library rows stay private.
"""
import argparse,json,hashlib
from pathlib import Path
from datetime import datetime,timezone
import pandas as pd
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--sidecar',type=Path,required=True);p.add_argument('--counts',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();assert not a.output.exists();s=pd.read_csv(a.sidecar,sep='\t',dtype=str,keep_default_na=False);c=pd.read_csv(a.counts,sep='\t');assert not s.duplicated(['dataset','library_id']).any();assert (c.n_cells>=0).all();totals=c.groupby(['dataset','library_id']).n_cells.sum();rows=[];details=[];opportunities=[]
 for ds,g in s.groupby('dataset'):
  assert all((ds,l) in totals.index for l in g.library_id);baseline=g[g.timepoint.isin(['Pre','PRE'])];solid=baseline[(baseline.tissue!='blood')&baseline.response_evaluable.eq('True')];patient=g.groupby('patient_id').agg(response_count=('response','nunique'),arm_count=('treatment_arm_normalized','nunique'));assert (patient.response_count<=1).all() and (patient.arm_count<=1).all();rows.append(dict(dataset=ds,clinical_rows=len(g),exact_count_library_keys=len(g),source_patient_tokens=g.patient_id.nunique(),baseline_rows=len(baseline),baseline_patient_tokens=baseline.patient_id.nunique(),baseline_solid_evaluable_rows=len(solid),baseline_solid_evaluable_patient_tokens=solid.patient_id.nunique(),patient_response_conflicts=0,patient_treatment_arm_conflicts=0,independent_validation_accepted=False))
  for (arm,tissue,phase),v in g.groupby(['treatment_arm_normalized','tissue','timepoint']):details.append(dict(dataset=ds,arm=arm,tissue=tissue,phase=phase,libraries=len(v),patient_tokens=v.patient_id.nunique(),evaluable_rows=int(v.response_evaluable.eq('True').sum())))
  if ds!='GSE169246':continue
  assert not solid.patient_id.duplicated().any();tc=c[c.dataset.eq(ds)];wide=tc.groupby(['library_id','major_class']).n_cells.sum().unstack(fill_value=0).reindex(solid.library_id,fill_value=0);wide['Macrophage']=tc[tc.annotation_final.str.contains('Macro',na=False)].groupby('library_id').n_cells.sum().reindex(solid.library_id,fill_value=0);wide.index=solid.index
  for groupname,base in [('baseline_all',solid),('baseline_solid',solid),('baseline_blood',baseline[(baseline.tissue=='blood')&baseline.response_evaluable.eq('True')])]:
   if groupname=='baseline_all':base=baseline[baseline.response_evaluable.eq('True')]
   if groupname=='baseline_all':continue
   for arm,v in base.groupby('treatment_arm_normalized'):
    unique=v.drop_duplicates('patient_id');rows.append(dict(dataset=ds,subset=groupname,arm=arm,libraries=len(v),patient_tokens=unique.patient_id.nunique(),PR=int(unique.response.eq('PR').sum()),SD=int(unique.response.eq('SD').sum()),PD=int(unique.response.eq('PD').sum())))
  for sender,receiver in [('T','Myeloid'),('T','Macrophage'),('T','B'),('Myeloid','B')]:
   for n in [20,50]:
    valid=(wide[sender]>=n)&(wide[receiver]>=n)
    for arm in ['ICI_plus_chemo','chemo']:
     v=solid[valid&solid.treatment_arm_normalized.eq(arm)];opportunities.append(dict(dataset=ds,subset='baseline_solid_evaluable',sender=sender,receiver=receiver,min_endpoint_cells=n,arm=arm,n_patient_tokens=v.patient_id.nunique(),n_libraries=len(v),PR=int(v.response.eq('PR').sum()),SD_PD=int(v.response.isin(['SD','PD']).sum()),group_definition='legacy_QC_major; Macrophage=legacy_Macro_name; not_new_172_mapping'))
 a.output.mkdir(parents=True);pd.DataFrame(rows).to_csv(a.output/'patient_opportunity_summary.tsv',sep='\t',index=False);pd.DataFrame(details).to_csv(a.output/'material_timepoint_arm_coverage.tsv',sep='\t',index=False);pd.DataFrame(opportunities).to_csv(a.output/'all16_baseline_pair_opportunities.tsv',sep='\t',index=False)
 j=dict(status='PATIENT_KEY_AND_CAPTURE_OPPORTUNITY_CHECKED_NOT_CLINICAL_MODEL_ACCEPTANCE',observed_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(__file__),clinical_sidecar_sha256=sha(a.sidecar),counts_sha256=sha(a.counts),exact_library_joins_verified=True,source_patient_tokens_not_cross_study_independence=True,no_new_p_values=True,no_expression_export=True,no_new_CM_projection=True,legacy_cell_annotation_pending_new_identity_review=True,outputs={f.name:sha(f) for f in a.output.iterdir()});(a.output/'inventory.json').write_text(json.dumps(j,indent=2)+'\n');print(json.dumps({k:v for k,v in j.items() if k!='outputs'}))
if __name__=='__main__':main()
