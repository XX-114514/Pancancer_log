"""Freeze complete native-lineage signature clinical family before effects."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    R=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);P=R/'NSCLC_context_native_programme_projection_20261016_v1'
    assert json.loads((P/'private_independent_audit.json').read_text())['status'].startswith('PASS_ALL56_SCOPE_')
    report=json.loads((P/'private_terminal.json').read_text());source=P/'private_native_author_subtype_programme_scores.tsv.gz'
    assert sha(source)==report['files'][source.name]
    C=R/'NSCLC_all_case_clinical_scope_20261016_v1';ct=json.loads((C/'private_terminal.json').read_text())
    assert ct['within_case_conflict_records']==0
    clinical_path=C/'private_all243_author_case_fields.tsv';assert sha(clinical_path)==ct['files'][clinical_path.name]
    O=R/'NSCLC_frozen_signature_clinical_preflight_20261016_v1';assert not O.exists();O.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),source_score_sha256=sha(source),clinical_table_sha256=sha(clinical_path),
        scope='Every actually measured native-lineage programme slot and every original compatible author subtype, not just the prioritised LR or one programme. All parent56 scope exclusions retained separately.',
        question='Post-treatment pathological-response association with frozen weighted gene expression within author subtypes, beyond subtype abundance, quality and available treatment/clinical factors. Not baseline prediction or proof of an LR mechanism.',
        outcome='Author pCR or MPR ->1; non-MPR ->0; unknown unavailable. Do not infer numerical pathology-rate semantics.',
        exposure='Recorded named PD1 agent or literal Yes, exclude No, unknown and placebo-ambiguous. Chemotherapy/targeted binary only where literal exposure is identifiable.',
        analysis_units='One author case per subtype per model; case blocks, no assertion of reviewed independent patient linkage.',
        response_model='Y=raw fixed-weight signature, standardized within model eligible cases; X=author response0/1. Controls age, sex, smoking, histology, ordinal stage I-IV, cycles, chemotherapy, targeted treatment, subtype within-major fraction, log cell count, full native mean logRNA/loggenes.',
        controls='All12 fixed controls retained, no automatic removal. Continuous controls standardized. Complete clinical cases, >=20 subtype cells, >=50 cases, >=10 each response group, residualdf>=10, full rank and condition<1e8, maximum leverage<.999999.',
        primary_inference='20000 Rademacher null wild bootstrap studentized HC3 t, nuisance-only HC2 residuals; case-block bootstrap2000 coefficients and CI only >=1900 valid; all planned subtype/signature slots BH with unestimable p=1 retained.',
        sensitivity='Unadjusted two-sided Mann-Whitney/rank-biserial; original frozen-residual score with identical new-cohort controls. No new clinical definition or rank/threshold variants.',
        algebra='Old frozen residual correction is linear in the four included quality/composition controls. Check correction lies in nuisance span before effects; adjusted raw/residual coefficient and t must agree on same raw-score SD.',
        interpretation='A frozen gene signature can have within-cohort value despite old correction extrapolation. This does not establish original fine-identity transfer, same state, independent project validation, clinical prediction or causality.',
        seed=2026101604,permutations=0,wild_bootstrap=20000,case_bootstrap=2000,minimum_valid_bootstrap=1900,
        finite_stop='One frozen full-family run plus its audit. Do not trim gene weights, covariates, response definitions or family based on these effects.')
    (O/'plan_before_any_clinical_effect.json').write_text(json.dumps(plan,indent=2)+'\n')
    c=pd.read_csv(clinical_path,sep='\t',dtype=str,keep_default_na=False).set_index('sampleID',drop=False);assert c.index.is_unique and len(c)==243
    c['response']=c.pathological_response.map({'pCR':1.,'MPR':1.,'non-MPR':0.})
    agents=['Tislelizumab','Pembrolizumab','Camrelizumab','Sintilimab','Nivolumab','Toripalimab']
    c['PD1_recorded']=c['anti-PD1_therapy'].apply(lambda x:(x=='Yes' or any(a in x for a in agents)) and 'Placebo' not in x)
    c['age_numeric']=pd.to_numeric(c.age,errors='coerce');c['male']=c.gender.map({'M':1.,'F':0.});c['smoker']=c.smoking_history.map({'Y':1.,'N':0.})
    c['LUSC']=c.cancer_type.map({'LUSC':1.,'LUAD':0.})
    def stage(x):
        for prefix,value in [('IV',4.),('III',3.),('II',2.),('I',1.)]:
            if x.startswith(prefix):return value
        return np.nan
    c['stage_ordinal']=c.pre_treatment_staging.map(stage);c['cycles_numeric']=pd.to_numeric(c.cycles,errors='coerce')
    def chemo(x):
        if x=='No':return 0.
        if x=='Yes' or any(s in x for s in ['Carboplatin','Cisplatin','Nedaplatin','Lobaplatin','Abraxane','Pemetrexed','Paclitaxel','Gemcitabine','Docetaxel']):return 1.
        return np.nan
    c['chemo_binary']=c.chemotherapy.map(chemo)
    c['targeted_binary']=c.targeted_therapy.map(lambda x:0. if x=='No' else 1. if any(s in x for s in ['Apatinib','Bevacizumab','Avastin','Endostar','Anlotinib']) else np.nan)
    clinical_controls=['age_numeric','male','smoker','LUSC','stage_ordinal','cycles_numeric','chemo_binary','targeted_binary']
    complete=c[clinical_controls+['response']].notna().all(axis=1);c['clinical_eligible']=c.PD1_recorded&complete
    c['exclusion_reasons']=[';'.join((['PD1_EXPOSURE_NOT_CONFIRMED'] if not r.PD1_recorded else [])+['MISSING_'+n for n in clinical_controls+['response'] if pd.isna(getattr(r,n))]) for r in c.itertuples()]
    c.to_csv(O/'private_all243_clinical_eligibility.tsv',sep='\t',index=False,na_rep='NA')
    scores=pd.read_csv(source,sep='\t');keys=['context','type','program','sub_cell_type'];scope=[];design=O/'private_frozen_designs';design.mkdir(mode=0o700)
    for number,(key,f) in enumerate(scores.groupby(keys,sort=True),1):
        case=c.reindex(f.sampleID);use=f.n_cells.ge(20).to_numpy()&case.clinical_eligible.to_numpy(bool)
        f=f.loc[use].copy();case=case.loc[use].copy();assert f.sampleID.is_unique and f.sampleID.tolist()==case.sampleID.tolist()
        controls=np.column_stack([case[clinical_controls].to_numpy(float),f[['fraction','logN','mean_logRNA','mean_loggenes']].to_numpy(float)])
        x=case.response.to_numpy(float);y=f.raw_weighted_score.to_numpy(float);old=f.score.to_numpy(float)
        n=len(x);sd=controls.std(0) if n else np.zeros(12);Z=np.column_stack([np.ones(n),(controls-controls.mean(0))/np.where(sd>1e-12,sd,1)]) if n else np.empty((0,13))
        D=np.column_stack([Z[:,0],x,Z[:,1:]])
        rank=int(np.linalg.matrix_rank(D)) if n else 0;df=n-rank
        condition=float(np.linalg.cond(D)) if n else np.inf
        groups=np.bincount(x.astype(int),minlength=2);status='DESIGN_READY_NO_EFFECT_COMPUTED';ysd=float(y.std(ddof=1)) if n>1 else np.nan
        if n<50 or groups.min()<10:status='INSUFFICIENT_CASE_OR_RESPONSE_GROUP_COVERAGE'
        elif (sd<=1e-12).any() or rank!=14:status='CONSTANT_OR_COLLINEAR_FIXED_CONTROLS'
        elif df<10 or condition>=1e8:status='LOW_RESIDUAL_DF_OR_ILL_CONDITIONED'
        elif not np.isfinite(ysd) or ysd<=1e-12:status='NO_SIGNATURE_VARIATION'
        else:
            h=np.sum(D*np.linalg.pinv(D).T,axis=1)
            if h.max()>=.999999:status='FULL_DESIGN_EXTREME_LEVERAGE'
        identity='CLIN'+str(number).zfill(4)
        record=dict(zip(keys,key));record.update(record=identity,source_family=f.source_family.iloc[0] if n else scores[(scores.context==key[0])&scores.type.eq(key[1])&scores.program.eq(key[2])&scores.sub_cell_type.eq(key[3])].source_family.iloc[0],
            eligible_cases=n,response0_cases=int(groups[0]),response1_cases=int(groups[1]),design_rank=rank,design_columns=14,residual_df=df,condition=condition,status=status)
        if status=='DESIGN_READY_NO_EFFECT_COMPUTED':
            delta=y-old;remainder=delta-Z@np.linalg.lstsq(Z,delta,rcond=None)[0];err=float(np.max(np.abs(remainder)))
            assert err<1e-7,'Frozen correction not in planned nuisance span'
            record['maximum_raw_residual_correction_span_error']=err
            np.savez_compressed(design/(identity+'.npz'),D=D,Z=Z,X=x,Y=(y-y.mean())/ysd,Y_original_frozen_residual=(old-old.mean())/ysd,
                case_blocks=f.sampleID.to_numpy(str),raw_signature_mean=y.mean(),raw_signature_sd=ysd)
        scope.append(record)
    table=pd.DataFrame(scope);assert len(table)>0 and not table.duplicated(keys).any()
    table.to_csv(O/'private_all_planned_signature_subtype_designs.tsv',sep='\t',index=False,na_rep='NA')
    terminal=dict(status='COMPLETE_FIXED_SIGNATURE_AUTHOR_SUBTYPE_CLINICAL_FAMILY_DESIGNS_NO_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        planned_models=len(table),status_counts=table.status.value_counts().to_dict(),clinical_eligible_cases=int(c.clinical_eligible.sum()),
        design_scope_sha256=sha(O/'private_all_planned_signature_subtype_designs.tsv'),source_projection_audit_sha256=sha(P/'private_independent_audit.json'),
        effects_computed=False,original_clinical_fields_changed=False,patient_linkage_reviewed=False)
    (O/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
