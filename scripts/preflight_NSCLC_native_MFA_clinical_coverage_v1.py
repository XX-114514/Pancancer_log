"""Verify original clinical fields and design coverage; no score-outcome effects."""
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import load, quota_gate, sha, VIEWS

COMMON=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']
FIELDS=['sampleID','age','gender','smoking_history','cancer_type','pre_treatment_staging',
        'anti-PD1_therapy','chemotherapy','targeted_therapy','cycles','pathological_response',
        'pathological_response_rate','radiological_response']


def main():
    data=load();root=data['root'];proof=root/'NSCLC_native_MFA_interpretation_20261016_v1/private_independent_interpretation_audit.json'
    assert json.loads(proof.read_text())['status'].startswith('PASS_')
    out=root/'NSCLC_native_MFA_clinical_coverage_20261016_v1';assert not out.exists();quota=quota_gate(root);out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),basis_interpretation_audit_sha256=sha(proof),
        question='Before new native-factor clinical inference, verify original author outcomes/controls and all discovery, holdout and complete-cohort design coverage. No factor-specific clinical association, score-effect coefficient, P or outcome-based selection.',
        scope='All208 frozen four-view cases. Existing discovery145 and holdout63 remain unchanged. All4 view designs retain all8 clinical controls and all original own-view quality variables.',
        outcome='Same prior author pCR/MPR=1 versus non-MPR=0; no new pathology-rate threshold. Post-neoadjuvant tissue, association only, not baseline prediction.',
        coverage_gates='At least50 complete author cases, >=10 per response group, all fixed controls varying, full rank, condition<1e8, residualdf>=10, leverage<.999999. Gates describe estimability, not scientific null.',
        new_clinical_tests=False,new_model_fit=False,quota=quota)
    (out/'plan_before_native_factor_clinical_coverage.json').write_text(json.dumps(plan,indent=2)+'\n')
    source=root/'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz';records=[]
    for chunk in pd.read_csv(source,usecols=FIELDS,dtype=str,keep_default_na=False,chunksize=200000):records.append(chunk.drop_duplicates())
    raw=pd.concat(records).drop_duplicates();assert raw.sampleID.nunique()==243 and len(raw)==243
    raw=raw.set_index('sampleID',drop=False)
    old=root/'NSCLC_frozen_signature_clinical_preflight_20261016_v1/private_all243_clinical_eligibility.tsv'
    c=pd.read_csv(old,sep='\t',dtype=str,keep_default_na=False).set_index('sampleID',drop=False)
    assert len(c)==243 and c.index.is_unique and set(c.index)==set(raw.index)
    for name in FIELDS:assert c.loc[raw.index,name].tolist()==raw[name].tolist(), 'Original clinical values changed'
    expected=pd.DataFrame(index=raw.index)
    expected['response']=raw.pathological_response.map({'pCR':1.,'MPR':1.,'non-MPR':0.})
    expected['age_numeric']=pd.to_numeric(raw.age,errors='coerce');expected['male']=raw.gender.map({'M':1.,'F':0.})
    expected['smoker']=raw.smoking_history.map({'Y':1.,'N':0.});expected['LUSC']=raw.cancer_type.map({'LUSC':1.,'LUAD':0.})
    def stage(s):
        return next((v for prefix,v in [('IV',4.),('III',3.),('II',2.),('I',1.)] if s.startswith(prefix)),np.nan)
    expected['stage_ordinal']=raw.pre_treatment_staging.map(stage)
    expected['cycles_numeric']=pd.to_numeric(raw.cycles,errors='coerce')
    def chemo(s):
        if s=='No':return 0.
        if s=='Yes' or any(x in s for x in ['Carboplatin','Cisplatin','Nedaplatin','Lobaplatin','Abraxane','Pemetrexed','Paclitaxel','Gemcitabine','Docetaxel']):return 1.
        return np.nan
    expected['chemo_binary']=raw.chemotherapy.map(chemo)
    expected['targeted_binary']=raw.targeted_therapy.map(lambda s:0. if s=='No' else 1. if any(x in s for x in ['Apatinib','Bevacizumab','Avastin','Endostar','Anlotinib']) else np.nan)
    agents=['Tislelizumab','Pembrolizumab','Camrelizumab','Sintilimab','Nivolumab','Toripalimab']
    expected['PD1_recorded']=raw['anti-PD1_therapy'].map(lambda s:(s=='Yes' or any(a in s for a in agents)) and 'Placebo' not in s)
    expected['clinical_eligible']=expected.PD1_recorded&expected[COMMON+['response']].notna().all(axis=1)
    for name in COMMON+['response']:
        values=pd.to_numeric(c.loc[expected.index,name],errors='coerce')
        assert np.allclose(values,expected[name],equal_nan=True)
    for name in ['PD1_recorded','clinical_eligible']:
        assert np.array_equal(c.loc[expected.index,name].map({'True':True,'False':False}),expected[name])
    axis=data['cases'];selected=expected.loc[axis.sampleID].copy();selected['partition']=axis.partition.to_numpy()
    selected['missing_reason']=[';'.join(['MISSING_'+name for name in COMMON+['response'] if pd.isna(getattr(r,name))]) for r in selected.itertuples()]
    assert selected.PD1_recorded.all()
    rows=[]
    for part in ['discovery','internal_holdout','all_frozen_four_view_cases']:
        take=np.ones(len(axis),bool) if part=='all_frozen_four_view_cases' else axis.partition.eq(part).to_numpy()
        eligible=take&selected.clinical_eligible.to_numpy(bool);clinical=selected.loc[eligible]
        for vi,view in enumerate(VIEWS):
            p=data['prepared'][vi];controls=np.column_stack([clinical[COMMON].to_numpy(float),p['quality'][eligible]])
            n=len(clinical);sd=controls.std(axis=0) if n else np.zeros(controls.shape[1]);names=COMMON+p['covariate_names'].astype(str).tolist()
            x=np.column_stack([np.ones(n),clinical.response.to_numpy(),(controls-controls.mean(axis=0))/np.where(sd>1e-12,sd,1)]) if n else np.empty((0,controls.shape[1]+2))
            rank=int(np.linalg.matrix_rank(x)) if n else 0;condition=float(np.linalg.cond(x)) if n else np.inf
            group=clinical.response.value_counts();n0=int(group.get(0.,0));n1=int(group.get(1.,0))
            status='FULL_FIXED_DESIGN_READY_NO_CLINICAL_EFFECTS'
            if n<50 or min(n0,n1)<10:status='INSUFFICIENT_CASE_OR_RESPONSE_GROUP_COVERAGE'
            elif (sd<=1e-12).any() or rank<x.shape[1]:status='CONSTANT_OR_COLLINEAR_FIXED_CONTROLS'
            elif n-rank<10 or condition>=1e8:status='LOW_DF_OR_ILL_CONDITIONED'
            elif np.max(np.sum(x*np.linalg.pinv(x).T,axis=1))>=.999999:status='EXTREME_LEVERAGE'
            rows.append(dict(partition=part,view=view,original_cases=int(take.sum()),complete_cases=n,response0=n0,response1=n1,
                design_columns=x.shape[1],design_rank=rank,residual_df=n-rank,condition=condition,
                constant_controls=';'.join(np.array(names)[sd<=1e-12]),status=status))
    table=pd.DataFrame(rows);assert len(table)==12
    table.to_csv(out/'private_all12_partition_view_clinical_design_coverage.tsv',sep='\t',index=False)
    selected.to_csv(out/'private_all208_clinical_eligibility.tsv',sep='\t',index=True)
    (out/'private_all_raw_pathology_rate_values.json').write_text(json.dumps(dict(values=raw.pathological_response_rate.value_counts().to_dict(),not_reinterpreted_as_new_outcome=True),indent=2)+'\n')
    terminal=dict(status='COMPLETE_NATIVE_MFA_CLINICAL_INPUT_RECHECK_AND_ALL12_DESIGN_COVERAGE_NO_CLINICAL_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),original_raw_author_case_fields_exact=True,prior_all243_codings_reproduced=True,
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
        clinical_tests_executed=False,new_model_fit=False,independent_patient_identity_reviewed=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
