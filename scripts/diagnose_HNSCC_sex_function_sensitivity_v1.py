"""Sex/age follow-up prompted by complete Hallmark diagnostics, no new P."""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from diagnose_HNSCC_clinical_composition_v1 import beta
from diagnose_HNSCC_full50_orthogonal_functions_v1 import rho


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inputs=root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1'
    sub=root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1'
    functions=root/'HNSCC_full50_orthogonal_function_diagnostics_20261014_v1'
    out=root/'HNSCC_sex_function_descriptive_20261014_v1';assert not out.exists();out.mkdir()
    (out/'plan_before_sensitivity.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),
      reason='all50 orthogonal diagnostics included androgen/spermatogenesis among strongest descriptors',
      scope='all3 primary endpoints and all3 subtype clinical scores, all50 functions within each subtype',
      clinical_variants=['treatment_HPV_sex','treatment_HPV_sex_age'],function_variant='RNA_genes_cellcount_sex_age',
      no_new_P_BH_or_primary_replacement=True,min_residual_df=10,full_rank_required=True),indent=2)+'\n')
    md=pd.read_csv(inputs/'private_eligible_libraries_clinical.tsv',sep='\t')
    baseline=md[md.timepoint.eq('pre')].set_index('patient_id')
    author=pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t')
    pre=author[author.time.eq('Pre')].set_index('patient_registry')
    gender=pre.Gender.astype(str).str.strip().str.lower()
    sex=gender.map({'male':1,'female':0,'m':1,'f':0});assert sex.notna().all()
    baseline['male']=sex.loc[baseline.index];baseline['age']=pd.to_numeric(pre.loc[baseline.index,'Age'])
    scores=pd.read_csv(inputs/'private_patient_library_frozen_scores.tsv',sep='\t')
    s=scores[scores.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='score')
    z=scores[scores.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='CRC_training_z')
    s['equal_mean_CRC_training_z']=z.mean(axis=1);baseline=baseline.loc[s.index]
    groups=[(name,s[name].to_numpy(),baseline) for name in s.columns]
    subtype=pd.read_csv(sub/'private_patient_subtype_scores.tsv',sep='\t')
    function_means=pd.read_csv(functions/'private_patient_subtype_function_means.tsv',sep='\t').set_index(['library_id','subtype_candidate'])
    functional=[]
    for name,g in subtype.groupby('subtype_candidate',sort=True):
        b=baseline.loc[g.patient_id].copy();b.index=g.index
        groups.append((name,g.score.to_numpy(),b))
        f=function_means.xs(name,level='subtype_candidate').loc[g.library_id]
        cov=np.column_stack((g.mean_logRNA,g.mean_log_genes,np.log1p(g.native_subtype_cells),b.male,b.age))
        x=g.score.to_numpy()
        for reference in f.columns:
            y=f[reference].to_numpy();r=rho(x,y,cov);loo=[]
            for i in range(len(x)):
                keep=np.arange(len(x))!=i;v=rho(x[keep],y[keep],cov[keep])
                if np.isfinite(v):loo.append(v)
            functional.append(dict(subtype_candidate=name,reference=reference,n_patients=len(g),rho=r,
              status='descriptive_estimated' if np.isfinite(r) else 'rank_or_df_insufficient',
              LOO_valid=len(loo),LOO_min=min(loo) if loo else np.nan,LOO_max=max(loo) if loo else np.nan))
    clinical=[]
    for name,x,b in groups:
        y=b.regression.to_numpy(float)
        C=np.column_stack((np.ones(len(b)),b.treatment_x.eq('D+T').to_numpy(float),b.HPV.eq('positive').to_numpy(float),b.male.to_numpy(float)))
        for variant,cov in [('treatment_HPV_sex',np.empty((len(b),0))),('treatment_HPV_sex_age',b.age.to_numpy(float)[:,None])]:
            point=beta(x,y,C,cov);loo=[]
            for i in range(len(x)):
                keep=np.arange(len(x))!=i;v=beta(x[keep],y[keep],C[keep],cov[keep])
                if np.isfinite(v):loo.append(v)
            clinical.append(dict(endpoint=name,variant=variant,n_patients=len(b),sex_counts=b.male.value_counts().to_dict(),
              beta=point,status='descriptive_estimated' if np.isfinite(point) else 'rank_or_df_insufficient',
              LOO_valid=len(loo),LOO_min=min(loo) if loo else np.nan,LOO_max=max(loo) if loo else np.nan))
    c=pd.DataFrame(clinical);f=pd.DataFrame(functional);assert len(c)==12 and len(f)==150
    c.to_csv(out/'all12_private_sex_clinical_diagnostics.tsv',sep='\t',index=False,na_rep='NA')
    f.to_csv(out/'all150_private_sex_function_diagnostics.tsv',sep='\t',index=False,na_rep='NA')
    report=dict(status='SEX_AGE_CLINICAL_AND_ALL50_FUNCTION_SENSITIVITY_COMPLETE_NO_NEW_P',
      finished_UTC=datetime.now(timezone.utc).isoformat(),baseline_sex_counts=baseline.male.value_counts().to_dict(),
      primary_clinical_tests_replaced=False,formal_state_labels_changed=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report));print(c.to_json(orient='records'))
    print(f[f.subtype_candidate.eq('Macro_C1QC')&f.reference.str.contains('INTERFERON|ANDROGEN|SPERMATOGENESIS',regex=True)].to_json(orient='records'))


if __name__=='__main__':main()
