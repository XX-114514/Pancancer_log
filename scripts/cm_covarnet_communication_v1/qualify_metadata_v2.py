from pathlib import Path
import pandas as pd,numpy as np,json
r=Path(__file__).resolve().parents[1];(r/'audit/metadata_v2').mkdir(exist_ok=True);v=r.parent/'20260818_pancancer_5M_merge_prep/10_final_v5_malignancy_covarnet_liana/20260825_v5_final_v1'
s=pd.read_csv(r/'private/sample_eligibility.tsv',sep='\t').set_index('sampleID');m=pd.read_csv(v/'metadata_audit_20260829/manifests/final_sample_metadata.tsv',sep='\t').set_index('sample_uid');assert m.index.is_unique
print('JOIN',len(s.index.intersection(m.index)),len(s),flush=True)
q=s.join(m.drop(columns=['gse_id','cohort'],errors='ignore'),how='left');q.to_csv(r/'private/sample_metadata_joined_v2.tsv',sep='\t')
for col in ['patient_id_source','metadata_evidence_level','geo_match_status','immunotherapy_timing','sample_class','manual_review_required']:
 print(col,q[col].fillna('NA').value_counts().to_dict(),flush=True)
unknown={'unknown','na','nan','none','not_available','not_reported','unclear',''}
q['patient_known']=(q.patient_id_source.isin(['GEO_GSM','local_h5ad_obs'])) & q.patient_id.fillna('').astype(str).str.lower().map(lambda x:x not in unknown|{'unavailable','not available','not reported'})
q['timing_known']=~q.immunotherapy_timing.fillna('').astype(str).str.lower().isin(unknown)
q['material_known']=~q.sample_class.fillna('').astype(str).str.lower().isin(unknown)
q['patient_key']=q.gse_id.astype(str)+'|'+q.patient_id.astype(str)
rows=[];strata=[]
for cm in sorted(q.CMT.dropna().unique()):
 for group,mask in [('High',q.CMT==cm),('Low',q.CMT.notna()&(q.CMT!=cm))]:
  x=q[mask];rows.append(dict(cm=cm,group=group,n_samples=len(x),n_sources=x.gse_id.nunique(),n_cancers=x.cancer_type.nunique(),n_patient_known=int(x.patient_known.sum()),n_known_patient_keys=x.loc[x.patient_known,'patient_key'].nunique(),n_timing_known=int(x.timing_known.sum()),n_material_known=int(x.material_known.sum()),median_cells=float(x.n_cells.median())))
  for col in ['sample_class','immunotherapy_timing','gse_id','cancer_type']:
   x[col].fillna('NA').value_counts().rename_axis(col).reset_index(name='n_samples').assign(cm=cm,group=group).to_csv(r/'audit/metadata_v2'/f'{cm}_{group}_{col}.tsv',sep='\t',index=False)
 x=q[q.patient_known&q.timing_known&q.material_known].copy();x['High']=x.CMT==cm
 for key,g in x.groupby(['gse_id','cancer_type','sample_class','immunotherapy_timing'],dropna=False):
  nh=g.loc[g.High,'patient_key'].nunique();nl=g.loc[~g.High,'patient_key'].nunique();strata.append(dict(cm=cm,source=key[0],cancer=key[1],material=key[2],timing=key[3],n_high_patients=nh,n_low_patients=nl,status='candidate_needs_patient_bridge_validation' if min(nh,nl)>=10 else 'insufficient_or_confounded'))
pd.DataFrame(rows).to_csv(r/'audit/metadata_v2/HighLow_eligibility.tsv',sep='\t',index=False);pd.DataFrame(strata).to_csv(r/'audit/metadata_v2/patient_strata_eligibility.tsv',sep='\t',index=False)
print('candidate strata',sum(x['status'].startswith('candidate') for x in strata),flush=True)
