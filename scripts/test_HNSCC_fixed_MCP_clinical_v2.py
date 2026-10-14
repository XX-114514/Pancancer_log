"""Patient baseline clinical validation under the pre-score frozen contract.
All patient-level tables and biological statistics remain on the server.
"""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr


def seed(key):return 20261014+int(hashlib.sha256(key.encode()).hexdigest()[:8],16)


def bh(p):
    p=np.asarray(p,float);ix=np.argsort(p);q=np.empty(len(p))
    q[ix]=np.minimum(1,np.minimum.accumulate((p[ix]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1])
    return q


def rank_standardize(x):
    r=rankdata(x);sd=r.std()
    if sd<1e-12:raise ValueError('Constant score')
    return (r-r.mean())/sd


def fitted_effect(x,y,nuisance):
    sx=rank_standardize(x);sy=rank_standardize(y)
    design=np.column_stack((nuisance,sx));df=len(y)-design.shape[1]
    if np.linalg.matrix_rank(design)<design.shape[1] or df<10:raise ValueError('Design rank or residual df')
    coef=np.linalg.lstsq(design,sy,rcond=None)[0]
    residual=sy-design@coef
    variance=(residual@residual)/df
    se=np.sqrt(variance*np.linalg.inv(design.T@design)[-1,-1])
    if se<1e-12:raise ValueError('Zero residual variance')
    return float(coef[-1]),float(coef[-1]/se),design,sy


def adjusted_test(x,y,nuisance,key):
    beta,t,design,sy=fitted_effect(x,y,nuisance)
    n=len(x);df=n-design.shape[1]
    reduced=nuisance@np.linalg.lstsq(nuisance,sy,rcond=None)[0];residual=sy-reduced
    inv=np.linalg.inv(design.T@design);solver=inv@design.T
    rng=np.random.default_rng(seed(key));extreme=0
    for _ in range(20000):
        permuted=reduced+residual[rng.permutation(n)]
        b=solver@permuted;r=permuted-design@b
        v=(r@r)/df;se=np.sqrt(max(v,0)*inv[-1,-1])
        perm_t=abs(b[-1]/se) if se>1e-12 else np.inf
        extreme+=perm_t>=abs(t)-1e-12
    boots=[]
    for _ in range(2000):
        ix=rng.integers(0,n,n)
        try:boots.append(fitted_effect(x[ix],y[ix],nuisance[ix])[0])
        except (ValueError,np.linalg.LinAlgError):pass
    low,high=map(float,np.quantile(boots,[.025,.975])) if len(boots)>=1800 else (np.nan,np.nan)
    return dict(beta=beta,t=t,p=(extreme+1)/20001,lower=low,upper=high,
                valid_bootstrap=len(boots),invalid_bootstrap=2000-len(boots),residual_df=df)


def raw_test(x,y,key,binary=False):
    n=len(x);rng=np.random.default_rng(seed(key))
    def effect(a,b):
        if binary:
            u=a[b==1];v=a[b==0]
            if min(len(u),len(v))<1:return np.nan
            return float(np.mean(np.sign(u[:,None]-v[None,:])))
        if np.std(a)<1e-12 or np.std(b)<1e-12:return np.nan
        return float(spearmanr(a,b).statistic)
    observed=effect(x,y)
    if not np.isfinite(observed):raise ValueError('Constant or missing outcome')
    extreme=0
    if binary:
        for _ in range(20000):extreme+=abs(effect(x,y[rng.permutation(n)]))>=abs(observed)-1e-12
    else:
        rx=rankdata(x);ry=rankdata(y)
        rx=(rx-rx.mean())/np.linalg.norm(rx-rx.mean())
        ry=(ry-ry.mean())/np.linalg.norm(ry-ry.mean())
        for _ in range(20000):extreme+=abs(rx@ry[rng.permutation(n)])>=abs(observed)-1e-12
    boots=[]
    for _ in range(2000):
        if binary:
            ix=np.concatenate([rng.choice(np.flatnonzero(y==group),int((y==group).sum()),replace=True) for group in [0,1]])
        else:ix=rng.integers(0,n,n)
        value=effect(x[ix],y[ix])
        if np.isfinite(value):boots.append(value)
    low,high=map(float,np.quantile(boots,[.025,.975])) if len(boots)>=1800 else (np.nan,np.nan)
    return dict(effect=observed,p=(extreme+1)/20001,lower=low,upper=high,
                valid_bootstrap=len(boots),invalid_bootstrap=2000-len(boots))


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inputs=root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1'
    out=root/'HNSCC_fixed_MCP_clinical_results_20261014_v2';assert not out.exists()
    terminal=json.loads((inputs/'private_terminal.json').read_text())
    assert terminal['status']=='CLINICAL_DEFINITIONS_AND_FIXED_PARTIAL_SCORES_PREPARED_NO_EFFECT_TESTS'
    frozen=json.loads((inputs/'plan_before_scores_or_effects.json').read_text())
    assert frozen['primary_tests']==6 and frozen['original_author_label_conflicts_resolved'] is False
    out.mkdir()
    (out/'execution_contract_before_effects.json').write_text(json.dumps(dict(
      created_UTC=datetime.now(timezone.utc).isoformat(),pre_score_plan_sha256=hashlib.sha256((inputs/'plan_before_scores_or_effects.json').read_bytes()).hexdigest(),
      producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      extra_author_sensitivity='same six-slot rank-biserial family on predefined concordant patients only',
      two_sided=True,not_response_prediction_validation=True),indent=2)+'\n')
    scores=pd.read_csv(inputs/'private_patient_library_frozen_scores.tsv',sep='\t')
    md=pd.read_csv(inputs/'private_eligible_libraries_clinical.tsv',sep='\t',keep_default_na=False)
    assert md.treatment_x.eq(md.treatment_y).all(), 'GEO/author treatment conflict'
    md['treatment']=md.treatment_x
    baseline=md[md.timepoint.eq('pre')].copy();assert baseline.patient_id.is_unique
    score=scores[scores.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='score')
    z=scores[scores.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='CRC_training_z')
    assert set(score.index)==set(baseline.patient_id)
    baseline=baseline.set_index('patient_id').loc[score.index]
    score['equal_mean_CRC_training_z']=z[['Tnode9','Macnode18']].mean(axis=1)
    # Author encoding is audited against the actual categorical values.
    hpv=baseline.HPV.str.strip().str.lower()
    known={'positive':1,'negative':0,'pos':1,'neg':0,'+':1,'-':0}
    hpv_value=hpv.map(known)
    complete=hpv_value.notna()&baseline.treatment.isin(['D','D+T'])
    nuisance=np.column_stack((np.ones(int(complete.sum())),baseline.loc[complete,'treatment'].eq('D+T').to_numpy(dtype=float),hpv_value[complete].to_numpy(dtype=float)))
    rows=[]
    for program in ['MCP1','MCP2']:
        for endpoint in frozen['endpoints']:
            for method in ['primary_adjusted_continuous','secondary_unadjusted_continuous','secondary_derived_binary','sensitivity_author_concordant_binary']:
                row=dict(program=program,endpoint=endpoint,method=method,status='unsupported_original_program_slot',n_patients=0,
                         beta=np.nan,effect=np.nan,p=np.nan,lower=np.nan,upper=np.nan)
                if program=='MCP1':
                    x=score[endpoint].to_numpy();y=baseline.regression.to_numpy(float)
                    if method=='primary_adjusted_continuous':
                        xx=x[complete];yy=y[complete];row['n_patients']=len(xx)
                        row['missing_covariate_patients']=int((~complete).sum())
                        if len(xx)<20:row['status']='coverage_insufficient'
                        else:
                            try:row.update(adjusted_test(xx,yy,nuisance,method+endpoint),status='estimated_partial_measurement_transfer')
                            except (ValueError,np.linalg.LinAlgError) as exc:row['status']='unestimable_'+str(exc)
                    else:
                        mask=np.ones(len(x),dtype=bool)
                        binary=method!='secondary_unadjusted_continuous'
                        yy=y if not binary else baseline.derived_response.to_numpy(int)
                        if method=='sensitivity_author_concordant_binary':
                            mask=baseline.author_binary_concordant.astype(str).str.lower().eq('true').to_numpy()
                        xx=x[mask];yy=yy[mask];row['n_patients']=len(xx)
                        if binary:row.update(responders=int((yy==1).sum()),nonresponders=int((yy==0).sum()))
                        if len(xx)<10 or (binary and min(int((yy==1).sum()),int((yy==0).sum()))<5):row['status']='coverage_insufficient'
                        else:
                            try:row.update(raw_test(xx,yy,method+endpoint,binary),status='estimated_partial_measurement_sensitivity')
                            except ValueError as exc:row['status']='unestimable_'+str(exc)
                rows.append(row)
    result=pd.DataFrame(rows);assert len(result)==24
    result['BH_six_slot_family']=np.nan
    for method,index in result.groupby('method').groups.items():
        assert len(index)==6
        q=bh(result.loc[index,'p'].fillna(1))
        result.loc[index,'BH_six_slot_family']=np.where(result.loc[index,'p'].notna(),q,np.nan)
    result.to_csv(out/'all24_private_clinical_statistics.tsv',sep='\t',index=False,na_rep='NA')
    (out/'private_terminal.json').write_text(json.dumps(dict(status='FIXED_MCP_BASELINE_CLINICAL_FAMILIES_COMPLETE_PRIVATE',
      finished_UTC=datetime.now(timezone.utc).isoformat(),clinical_outcome='pathological regression, not RECIST',
      baseline_patients=len(baseline),complete_primary_covariate_patients=int(complete.sum()),
      original_author_binary_overwritten=False,not_verified_direct_interaction=True),indent=2)+'\n')
    print('FIXED_MCP_BASELINE_CLINICAL_FAMILIES_COMPLETE_RESULTS_REMAIN_SERVER')


if __name__=='__main__':main()
