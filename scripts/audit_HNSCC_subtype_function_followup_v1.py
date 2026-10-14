"""Independent numerical audit of all subtype/functional/sex diagnostics."""
import json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from audit_HNSCC_fixed_MCP_clinical_v1 import ranks


def effect(x,y,clinical,cov):
    x=ranks(x);y=ranks(y)
    if min(x.std(),y.std())<1e-12:return np.nan
    x=(x-x.mean())/x.std();y=(y-y.mean())/y.std()
    cv=np.column_stack([ranks(cov[:,i]) for i in range(cov.shape[1])]) if cov.shape[1] else cov
    N=np.column_stack((clinical,cv));X=np.column_stack((N,x))
    if len(x)-X.shape[1]<10 or np.linalg.matrix_rank(X)<X.shape[1]:return np.nan
    Q,_=np.linalg.qr(N,mode='reduced');xr=x-Q@(Q.T@x)
    return float(xr@y/(xr@xr))


def partial(x,y,cov):
    N=np.column_stack((np.ones(len(x)),*[ranks(cov[:,i]) for i in range(cov.shape[1])]))
    if len(x)-N.shape[1]<10 or np.linalg.matrix_rank(N)<N.shape[1]:return np.nan
    Q,_=np.linalg.qr(N,mode='reduced');x=ranks(x);y=ranks(y)
    ex=x-Q@(Q.T@x);ey=y-Q@(Q.T@y)
    if min(ex.std(),ey.std())<1e-12:return np.nan
    return float(ex@ey/(np.linalg.norm(ex)*np.linalg.norm(ey)))


def same(a,b):assert np.isclose(a,b,atol=1e-9,rtol=1e-9,equal_nan=True),(a,b)


def check_loo(fn,x,y,args,row):
    same(fn(x,y,*args),getattr(row,'beta',getattr(row,'rho',np.nan)))
    values=[]
    for i in range(len(x)):
        keep=np.arange(len(x))!=i;v=fn(x[keep],y[keep],*[a[keep] for a in args])
        if np.isfinite(v):values.append(v)
    assert len(values)==row.LOO_valid
    same(min(values) if values else np.nan,row.LOO_min);same(max(values) if values else np.nan,row.LOO_max)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inputs=root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1'
    sub=root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1'
    funcs=root/'HNSCC_full50_orthogonal_function_diagnostics_20261014_v1'
    sexroot=root/'HNSCC_sex_function_descriptive_20261014_v1'
    out=root/'HNSCC_subtype_function_independent_audit_20261014_v1';assert not out.exists()
    baseline=pd.read_csv(inputs/'private_eligible_libraries_clinical.tsv',sep='\t');baseline=baseline[baseline.timepoint.eq('pre')].set_index('patient_id')
    author=pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv',sep='\t')
    pre=author[author.time.eq('Pre')].set_index('patient_registry')
    male=pre.Gender.str.strip().str.lower().map({'m':1,'male':1,'f':0,'female':0})
    assert male.notna().all();baseline['male']=male.loc[baseline.index];baseline['age']=pre.Age.loc[baseline.index]
    programme=pd.read_csv(sub/'private_patient_subtype_scores.tsv',sep='\t')
    assert not programme.duplicated(['patient_id','subtype_candidate']).any()
    assert programme.native_subtype_cells.ge(20).all()
    clinical=pd.read_csv(sub/'all_private_subtype_descriptive_diagnostics.tsv',sep='\t')
    functions=pd.read_csv(funcs/'private_patient_subtype_function_means.tsv',sep='\t').set_index(['library_id','subtype_candidate'])
    full=pd.read_csv(funcs/'all300_private_descriptive_functions.tsv',sep='\t')
    sexclinical=pd.read_csv(sexroot/'all12_private_sex_clinical_diagnostics.tsv',sep='\t')
    sexfunction=pd.read_csv(sexroot/'all150_private_sex_function_diagnostics.tsv',sep='\t')
    for subtype,g in programme.groupby('subtype_candidate',sort=True):
        b=baseline.loc[g.patient_id];x=g.score.to_numpy();y=b.regression.to_numpy(float)
        C=np.column_stack((np.ones(len(g)),b.treatment_x.eq('D+T').to_numpy(float),b.HPV.eq('positive').to_numpy(float)))
        cov=dict(clinical_arm_HPV=np.empty((len(g),0)),clinical_arm_HPV_age=b.age.to_numpy(float)[:,None],
                 clinical_arm_HPV_native_quality_count=np.column_stack((g.mean_logRNA,g.mean_log_genes,np.log1p(g.native_subtype_cells))))
        for row in clinical[clinical.subtype_candidate.eq(subtype)].itertuples():check_loo(effect,x,y,[C,cov[row.variant]],row)
        for row in sexclinical[sexclinical.endpoint.eq(subtype)].itertuples():
            cv=b.age.to_numpy(float)[:,None] if row.variant.endswith('_age') else np.empty((len(g),0))
            check_loo(effect,x,y,[np.column_stack((C,b.male.to_numpy(float))),cv],row)
        f=functions.xs(subtype,level='subtype_candidate').loc[g.library_id]
        q=np.column_stack((g.mean_logRNA,g.mean_log_genes,np.log1p(g.native_subtype_cells)))
        for row in full[full.subtype_candidate.eq(subtype)].itertuples():
            cv=q if row.variant=='native_RNA_detected_cellcount_rank_residual' else np.empty((len(g),0))
            check_loo(partial,x,f[row.reference].to_numpy(),[cv],row)
        for row in sexfunction[sexfunction.subtype_candidate.eq(subtype)].itertuples():
            check_loo(partial,x,f[row.reference].to_numpy(),[np.column_stack((q,b.male,b.age))],row)
    score=pd.read_csv(inputs/'private_patient_library_frozen_scores.tsv',sep='\t')
    s=score[score.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='score')
    z=score[score.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='CRC_training_z');s['equal_mean_CRC_training_z']=z.mean(axis=1)
    b=baseline.loc[s.index];C=np.column_stack((np.ones(len(b)),b.treatment_x.eq('D+T').to_numpy(float),b.HPV.eq('positive').to_numpy(float),b.male.to_numpy(float)))
    for row in sexclinical[sexclinical.endpoint.isin(s.columns)].itertuples():
        cv=b.age.to_numpy(float)[:,None] if row.variant.endswith('_age') else np.empty((len(b),0))
        check_loo(effect,s[row.endpoint].to_numpy(),b.regression.to_numpy(float),[C,cv],row)
    assert len(clinical)==9 and len(full)==300 and len(sexclinical)==12 and len(sexfunction)==150
    report=dict(status='PASS_ALL471_DESCRIPTIVE_EFFECTS_AND_ALL_LOO_INCLUDING_UNESTIMABLE',
      finished_UTC=datetime.now(timezone.utc).isoformat(),no_new_P_or_BH=True,formal_states_changed=False)
    out.mkdir();(out/'private_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
