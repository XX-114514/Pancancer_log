"""All30 fixed new-atlas baseline components; original clinical family retained."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata


def seed(key):return 20261016+int(hashlib.sha256(key.encode()).hexdigest()[:8],16)


def rankz(x):
    x=np.asarray(x,float)
    if not np.isfinite(x).all() or np.ptp(x)<=1e-10*max(1,float(np.max(np.abs(x)))):
        raise ValueError('constant_or_numerically_negligible_activity')
    r=rankdata(x);return (r-r.mean())/r.std()


def fitted(x,y,N):
    X=rankz(x);Y=rankz(y);D=np.column_stack((N,X));df=len(y)-D.shape[1]
    if len(y)<20 or df<10 or np.linalg.matrix_rank(D)<D.shape[1]:raise ValueError('coverage_or_design_rank')
    inverse=np.linalg.inv(D.T@D);solver=inverse@D.T;b=solver@Y;e=Y-D@b
    se=np.sqrt(max(float(e@e)/df*inverse[-1,-1],0))
    if se<1e-12:raise ValueError('zero_residual_variance')
    return float(b[-1]),float(b[-1]/se),D,Y,solver,float(inverse[-1,-1])


def test(x,y,N,key):
    b,t,D,Y,solver,v=fitted(x,y,N);df=len(y)-D.shape[1]
    reduced=N@np.linalg.lstsq(N,Y,rcond=None)[0];residual=Y-reduced;rng=np.random.default_rng(seed(key));extreme=0
    for start in range(0,20000,250):
        P=np.column_stack([reduced+residual[rng.permutation(len(y))] for _ in range(min(250,20000-start))])
        pb=solver@P;pe=P-D@pb;se=np.sqrt(np.maximum(np.sum(pe*pe,axis=0)/df*v,0))
        pt=np.divide(pb[-1],se,out=np.full(len(se),np.inf),where=se>1e-12)
        extreme+=np.count_nonzero(np.abs(pt)>=abs(t)-1e-12)
    boots=[]
    for _ in range(2000):
        ix=rng.integers(0,len(y),len(y))
        try:boots.append(fitted(x[ix],y[ix],N[ix])[0])
        except (ValueError,np.linalg.LinAlgError):pass
    lo,hi=np.quantile(boots,[.025,.975]) if len(boots)>=1800 else (np.nan,np.nan)
    return dict(status='estimated_internal_baseline_CM_association',beta=b,t=t,p=(extreme+1)/20001,
        lower=float(lo),upper=float(hi),valid_bootstrap=len(boots),residual_df=df)


def bh(p):
    p=np.asarray(p,float);order=np.argsort(p);q=np.empty(len(p))
    q[order]=np.minimum(1,np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q


def descriptor(x,y,N):
    point=np.nan;loo=[];status='rank_or_coverage_insufficient'
    try:point=fitted(x,y,N)[0];status='descriptive_estimated'
    except (ValueError,np.linalg.LinAlgError):pass
    for i in range(len(y)):
        keep=np.arange(len(y))!=i
        try:loo.append(fitted(x[keep],y[keep],N[keep])[0])
        except (ValueError,np.linalg.LinAlgError):pass
    return dict(beta=point,LOO_valid=len(loo),LOO_min=min(loo) if loo else np.nan,LOO_max=max(loo) if loo else np.nan,status=status)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);inputs=root/'HNSCC_new_CM_baseline_preparation_20261016_v1'
    audit=json.loads((inputs/'private_independent_input_audit.json').read_text());assert audit['status']=='PASS_ALL30_BASELINE_RAW_FREQUENCY_BVLS_AND_EXACT_CLINICAL_KEYS'
    for filename,key in [('private_patient_baseline_component_scores.tsv','scores_sha256'),('private_patient_clinical_covariates.tsv','clinical_sha256')]:
        assert hashlib.sha256((inputs/filename).read_bytes()).hexdigest()==audit[key]
    scores=pd.read_csv(inputs/'private_patient_baseline_component_scores.tsv',sep='\t')
    clinical=pd.read_csv(inputs/'private_patient_clinical_covariates.tsv',sep='\t')
    scope=pd.read_csv(inputs/'private_all30_component_scope.tsv',sep='\t');assert len(scope)==30
    out=root/'HNSCC_new_CM_clinical_results_20261016_v1';assert not out.exists();out.mkdir()
    contract=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_audit_sha256=hashlib.sha256((inputs/'private_independent_input_audit.json').read_bytes()).hexdigest(),
        primary_family=30,model='rank regression ~ ranked fixed-model baseline activity + D_plus_T + HPV',
        fixed_smoothing='Primary uses original W S; raw-W projection only descriptive, no added P.',
        numerical_variation_gate='range>1e-10*max(1,maxabs(raw activity)), fixed before first effects',
        inference='20000 residual permutations,2000 patient bootstrap; invalid counts kept,1800 minimum for CI',
        descriptive_family='All30 * primary/unadjusted/rawW/age_sex, all LOO, no extra significance tests',
        internal_exploration_not_independent_prediction=True)
    (out/'execution_contract_before_effects.json').write_text(json.dumps(contract,indent=2)+'\n')
    rows=[];descriptors=[]
    for record in scope.itertuples():
        mask=scores.scope.eq(record.scope)&scores.component.eq(record.component)
        frame=scores[mask].merge(clinical[clinical.scope.eq(record.scope)],on=['scope','patient_id'],validate='one_to_one').sort_values('patient_id')
        row=dict(scope=record.scope,component=record.component,n_input_patients=len(frame),n_patients=0,
            status='missing_or_unidentifiable_projection',beta=np.nan,t=np.nan,p=np.nan,lower=np.nan,upper=np.nan,valid_bootstrap=0)
        complete=frame[['D_plus_T','HPV_positive','regression']].notna().all(axis=1)
        f=frame[complete];row['n_patients']=len(f);row['missing_mandatory_covariates']=int((~complete).sum())
        x=f.activity.to_numpy(float);y=f.regression.to_numpy(float)
        N=np.column_stack((np.ones(len(f)),f.D_plus_T,f.HPV_positive))
        if record.basis_identifiable and np.isfinite(x).all():
            try:row.update(test(x,y,N,record.scope+'|'+record.component))
            except (ValueError,np.linalg.LinAlgError) as exc:row['status']='unestimable_'+str(exc)
        rows.append(row)
        for variant,z,design,subset in [('primary',x,N,f),('unadjusted',x,np.ones((len(f),1)),f),('raw_W_only',f.raw_W_activity.to_numpy(float),N,f)]:
            descriptors.append(dict(scope=record.scope,component=record.component,variant=variant,n_patients=len(subset),**descriptor(z,y,design)))
        valid=f[['age','male']].notna().all(axis=1);demographic=f[valid];z=demographic.activity.to_numpy(float);yy=demographic.regression.to_numpy(float)
        try:
            age=rankz(demographic.age.to_numpy(float));D=np.column_stack((np.ones(len(demographic)),demographic.D_plus_T,demographic.HPV_positive,demographic.male,age))
            description=descriptor(z,yy,D)
        except ValueError:description=dict(beta=np.nan,LOO_valid=0,LOO_min=np.nan,LOO_max=np.nan,status='age_missing_or_constant')
        descriptors.append(dict(scope=record.scope,component=record.component,variant='age_sex',n_patients=len(demographic),**description))
    result=pd.DataFrame(rows);assert len(result)==30
    result['BH_all30']=np.where(result.p.notna(),bh(result.p.fillna(1)),np.nan)
    result.to_csv(out/'all30_private_clinical_statistics.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(descriptors).to_csv(out/'all120_private_descriptive_and_LOO.tsv',sep='\t',index=False,na_rep='NA')
    report=dict(status='ALL30_NEW_CM_BASELINE_CLINICAL_TESTS_AND120_DESCRIPTORS_COMPLETE_PRIVATE',finished_UTC=datetime.now(timezone.utc).isoformat(),
        original_20261010_9CM_changed=False,internal_clinical_association_not_prediction=True,
        result_sha256=hashlib.sha256((out/'all30_private_clinical_statistics.tsv').read_bytes()).hexdigest(),
        direct_interaction_verified=False,new_expression_download=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    print(result.to_json(orient='records'))


if __name__=='__main__':main()
