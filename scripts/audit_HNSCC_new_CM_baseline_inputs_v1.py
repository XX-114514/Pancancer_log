"""Independent raw-frequency/BVLS projection and exact clinical-key audit."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.optimize import lsq_linear


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);project=Path(os.environ['PANCANCER_PROJECT_ROOT'])
    atlas=project/'Project_v5/runs/20261013_immunotherapy_atlas_v1'
    original=root/'HNSCC_new_CM_frozen_model_export_20261016_v2'
    inputs=root/'HNSCC_new_CM_baseline_preparation_20261016_v1'
    terminal=json.loads((inputs/'private_terminal.json').read_text())
    scores_path=inputs/'private_patient_baseline_component_scores.tsv';clinical_path=inputs/'private_patient_clinical_covariates.tsv'
    assert hashlib.sha256(scores_path.read_bytes()).hexdigest()==terminal['scores_sha256']
    assert hashlib.sha256(clinical_path.read_bytes()).hexdigest()==terminal['clinical_sha256']
    scores=pd.read_csv(scores_path,sep='\t');clinical=pd.read_csv(clinical_path,sep='\t')
    assert not scores.duplicated(['scope','patient_id','component']).any()
    assert not clinical.duplicated(['scope','patient_id']).any()
    author=pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t')
    author=author[author.time.eq('Pre')].copy();author['sampleID']='GSE286827::'+author.library_id.astype(str)
    assert author.sampleID.is_unique and author.patient_registry.is_unique
    checked=0;errors=[];modules=0
    for scope in sorted(scores.scope.unique()):
        raw=pd.read_csv(atlas/'CM'/scope/'inputs/frequency_private.tsv',sep='\t',index_col=0)
        X=raw.iloc[:,1:].to_numpy(float)
        lower=np.array([min(row) for row in X]);upper=np.array([max(row) for row in X])
        normalized=(X-lower[:,None])/(upper-lower)[:,None]
        A=pd.read_csv(original/scope/'A_private.tsv',sep='\t',index_col=0)
        W=pd.read_csv(original/scope/'W_private.tsv',sep='\t',index_col=0)
        original_norm=pd.read_csv(original/scope/'norm_private.tsv',sep='\t',index_col=0)
        assert raw.index.equals(A.index) and np.allclose(normalized,original_norm,atol=1e-12,rtol=1e-12)
        meta=pd.read_csv(atlas/'CM'/scope/'inputs/sample_metadata_private.tsv',sep='\t')
        chosen=meta[meta.dataset.eq('GSE286827')&meta.timepoint.eq('pre')]
        keys=chosen[['sampleID']].merge(author,on='sampleID',validate='one_to_one')
        assert len(keys)==len(chosen)
        scope_clinical=clinical[clinical.scope.eq(scope)].set_index('patient_id')
        assert set(scope_clinical.index)==set(keys.patient_registry)
        modules+=A.shape[1]
        for patient,group in keys.groupby('patient_registry'):
            record=scope_clinical.loc[patient];first=group.iloc[0]
            assert float(first['Tumor.regression'])==record.regression
            assert record.D_plus_T==int(first.neo_type=='D+T')
            assert str(first['HPV/p16']).strip().lower() in ['positive','negative']
            assert record.HPV_positive==int(str(first['HPV/p16']).strip().lower()=='positive')
            positions=[list(raw.columns[1:]).index(k) for k in group.sampleID]
            y=np.mean(normalized[:,positions],axis=1)
            expected=scores[scores.scope.eq(scope)&scores.patient_id.eq(patient)].set_index('component').loc[A.columns]
            for matrix,column in [(A,'activity'),(W,'raw_W_activity')]:
                solution=lsq_linear(matrix.to_numpy(),y,bounds=(0,np.inf),method='bvls',tol=1e-12,max_iter=10000)
                assert solution.success
                error=float(np.max(np.abs(solution.x-expected[column].to_numpy())))
                assert np.allclose(solution.x,expected[column],atol=1e-7,rtol=1e-7)
                errors.append(error)
            checked+=1
    assert modules==30 and checked==len(clinical)
    path=inputs/'private_independent_input_audit.json';assert not path.exists()
    report=dict(status='PASS_ALL30_BASELINE_RAW_FREQUENCY_BVLS_AND_EXACT_CLINICAL_KEYS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),scope_patient_rows=checked,
        maximum_independent_projection_error=max(errors),scores_sha256=terminal['scores_sha256'],clinical_sha256=terminal['clinical_sha256'],
        original_mixed_time_H_used=False,clinical_effects_computed=False)
    path.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
