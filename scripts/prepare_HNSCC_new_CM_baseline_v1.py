"""Fixed original model projection of baseline-only composition; no effects."""
import hashlib, json, os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import nnls

SCOPES=['whole_TME','T_lineage_biopsy','whole_TME_primary_tumor','T_whole_TME_primary_tumor']


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def checked_nnls(A,y):
    h,residual=nnls(A,y,maxiter=10000)
    g=A.T@(A@h-y);tol=1e-7*(1+np.linalg.norm(A,2)*np.linalg.norm(y))
    assert np.isfinite(h).all() and (h>=0).all() and g.min()>=-tol
    assert np.max(np.abs(h*g))<=tol*(1+np.linalg.norm(h))
    return h,float(residual)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    original=root/'HNSCC_new_CM_frozen_model_export_20261016_v2'
    terminal=json.loads((original/'private_terminal.json').read_text())
    assert terminal['status']=='ALL4_ORIGINAL_MINMAX_AND_NSNMF_W_S_H_PARITY_EXPORTED_NO_CLINICAL_EFFECTS'
    index=json.loads((atlas/'manifests/CM_recovery_v3_public.json').read_text())
    out=root/'HNSCC_new_CM_baseline_preparation_20261016_v1';assert not out.exists();out.mkdir()
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=sha(Path(__file__)),scopes=SCOPES,primary_slots=30,
        original_export_terminal_sha256=sha(original/'private_terminal.json'),
        baseline_only='Exact GSE286827::original library key; pre only, average normalized baseline frequencies within reviewed author patient, then NNLS.',
        model='All original W and smoothing S fixed. Primary A=W S reproduces original fitted model; raw W-only is descriptive sensitivity, no second significance family.',
        scaling='Use full-discovery fixed node minimum/range; no baseline scaling fit.',
        identifiability='Full effective-basis column rank and condition<=1e8; unavailable scores retain all planned slots.',
        clinical_outcome='Author continuous Tumor.regression; no binary label repair or derived threshold primary.',
        primary_model='rank-standardized continuous regression ~ rank-standardized module activity + D_plus_T + HPV_positive',
        model_gate='At least20 reviewed patients,10 residual df,full design rank; unknown mandatory covariates excluded explicitly, no automatic covariate deletion.',
        inference='20000 patient residual permutations,2000 patient bootstrap reranking scores/outcome; at least1800 valid for interval.',
        primary_BH='All30 frozen scope/component slots; missing internally P1, shown unavailable.',
        descriptive_sensitivities=['unadjusted','raw_W_only','age_and_sex','leave_one_patient_all_variants'],
        interpretation='Internal clinical association; W/scaling discovery included post-treatment data. Not independent prediction or original v20261010 9CM validation.',
        no_new_expression_download=True,no_clinical_effects_computed=True)
    (out/'plan_before_projection_or_effects.json').write_text(json.dumps(plan,indent=2)+'\n')
    author_path=root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    author=pd.read_csv(author_path,sep='\t',keep_default_na=False)
    author['sampleID']='GSE286827::'+author.library_id.astype(str)
    assert author.sampleID.is_unique
    author=author.set_index('sampleID')
    scores=[];patients=[];coverage=[];projection=[];manifest=[]
    for scope in SCOPES:
        od=original/scope;A=pd.read_csv(od/'A_private.tsv',sep='\t',index_col=0);W=pd.read_csv(od/'W_private.tsv',sep='\t',index_col=0)
        H=pd.read_csv(od/'H_private.tsv',sep='\t',index_col=0);F=pd.read_csv(od/'F_private.tsv',sep='\t',index_col=0)
        normalized=pd.read_csv(od/'norm_private.tsv',sep='\t',index_col=0)
        meta=pd.read_csv(atlas/'CM'/scope/'inputs/sample_metadata_private.tsv',sep='\t',keep_default_na=False)
        fq=pd.read_csv(atlas/'CM'/scope/'inputs/frequency_private.tsv',sep='\t',index_col=0)
        assert A.index.equals(W.index) and A.index.equals(normalized.index) and A.index.equals(fq.index)
        assert list(A.columns)==list(H.index)
        npA=A.to_numpy();npW=W.to_numpy();rank=int(np.linalg.matrix_rank(npA));condition=float(np.linalg.cond(npA))
        identifiable=rank==A.shape[1] and condition<=1e8
        recovered=[]
        if identifiable:
            for j in range(H.shape[1]):
                h,r=checked_nnls(npA,F.iloc[:,j].to_numpy())
                assert np.allclose(h,H.iloc[:,j],atol=1e-6,rtol=1e-7) and r<1e-6
                recovered.append(float(np.max(np.abs(h-H.iloc[:,j].to_numpy()))))
        take=meta[meta.dataset.eq('GSE286827') & meta.timepoint.eq('pre')]
        assert take.sampleID.is_unique and set(take.sampleID)<=set(author.index)
        joined=take.merge(author[['patient_registry','time','Tumor.regression','neo_type','HPV/p16','Age','Gender']],left_on='sampleID',right_index=True,validate='one_to_one')
        assert joined.time.eq('Pre').all()
        assert joined.groupby('independent_unit').patient_registry.nunique().max()==1
        for patient,group in joined.groupby('patient_registry',sort=True):
            for field in ['Tumor.regression','neo_type','HPV/p16','Age','Gender']:
                assert group[field].nunique()==1
            first=group.iloc[0];regression=float(first['Tumor.regression']);assert 0<=regression<=100
            hpv=str(first['HPV/p16']).strip().lower();treatment=str(first.neo_type).strip()
            hpv_value={'positive':1,'negative':0,'pos':1,'neg':0,'+':1,'-':0}.get(hpv,np.nan)
            treatment_value={'D':0,'D+T':1}.get(treatment,np.nan)
            gender=str(first.Gender).strip().lower();male={'m':1,'male':1,'f':0,'female':0}.get(gender,np.nan)
            clinical_ok=np.isfinite(hpv_value) and np.isfinite(treatment_value)
            patient_row=dict(scope=scope,patient_id=patient,regression=regression,HPV_positive=hpv_value,D_plus_T=treatment_value,
                age=pd.to_numeric(first.Age,errors='coerce'),male=male,baseline_libraries=len(group),mandatory_covariates_complete=clinical_ok)
            patients.append(patient_row)
            y=normalized.loc[:,group.sampleID].mean(axis=1).to_numpy()
            if identifiable:
                h,residual=checked_nnls(npA,y);raw,rawres=checked_nnls(npW,y)
            else:
                h=np.full(A.shape[1],np.nan);raw=h.copy();residual=rawres=np.nan
            total=float(np.nansum(h))
            projection.append(dict(scope=scope,patient_id=patient,relative_residual=residual/max(np.linalg.norm(y),1e-12),
                raw_W_relative_residual=rawres/max(np.linalg.norm(y),1e-12),effective_basis_identifiable=identifiable))
            for j,component in enumerate(A.columns):
                scores.append(dict(scope=scope,patient_id=patient,component=component,activity=h[j],raw_W_activity=raw[j],
                    activity_fraction=h[j]/total if total>0 else np.nan,
                    status='original_frozen_model_baseline_projection' if identifiable else 'effective_basis_not_identifiable'))
        accepted=next(d for d in index['results'] if d['scope']==scope)
        members=pd.read_csv(atlas/'CM'/scope/accepted['discovery_attempt']/'each_node.tsv',sep='\t')
        for component in A.columns:
            weight=W[component];by_major=weight.groupby(fq.iloc[:,0]).sum();by_major=by_major/by_major.sum()
            subset=members[members.cm.astype(str).eq(component)] if 'cm' in members else pd.DataFrame()
            coverage.append(dict(scope=scope,component=component,nodes=W.shape[0],rank=A.shape[1],basis_identifiable=identifiable,
                effective_condition=condition,reviewed_baseline_patients=joined.patient_registry.nunique(),
                supported_network_members=len(subset),n_major_with_positive_weight=int((by_major>0).sum()),
                weighted_major_composition=json.dumps({str(k):float(v) for k,v in by_major.items()},sort_keys=True)))
        manifest.append(dict(scope=scope,maximum_original_fitted_H_replay_error=max(recovered) if recovered else None,
            effective_basis_sha256=sha(od/'A_private.tsv'),normalized_matrix_sha256=sha(od/'norm_private.tsv')))
    assert len(coverage)==30 and len({(r['scope'],r['component']) for r in coverage})==30
    pd.DataFrame(scores).to_csv(out/'private_patient_baseline_component_scores.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(patients).to_csv(out/'private_patient_clinical_covariates.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(projection).to_csv(out/'private_patient_projection_diagnostics.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(coverage).to_csv(out/'private_all30_component_scope.tsv',sep='\t',index=False,na_rep='NA')
    report=dict(status='ALL30_FROZEN_MODEL_BASELINE_COMPONENTS_PREPARED_NO_CLINICAL_EFFECTS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        scopes=manifest,components=30,unique_author_patients=len({r['patient_id'] for r in patients}),
        missing_mandatory_covariate_scope_patient_rows=sum(not r['mandatory_covariates_complete'] for r in patients),
        original_H_not_used_as_baseline=True,no_W_S_or_scaling_refit=True,
        scores_sha256=sha(out/'private_patient_baseline_component_scores.tsv'),clinical_sha256=sha(out/'private_patient_clinical_covariates.tsv'))
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
