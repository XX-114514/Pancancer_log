"""Independent rank, residual-projection and BH replay of real clinical tests."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd


def ranks(x):
    x=np.asarray(x);order=np.argsort(x,kind='stable');out=np.empty(len(x),float);i=0
    while i<len(x):
        j=i+1
        while j<len(x) and x[order[j]]==x[order[i]]:j+=1
        out[order[i:j]]=(i+1+j)/2;i=j
    return out


def unit(x):
    x=ranks(x)
    if x.std()<1e-12:raise ValueError('constant')
    return (x-x.mean())/x.std()


def coefficient(x,y,N):
    X=np.column_stack((N,unit(x)));Y=unit(y)
    if len(y)-X.shape[1]<10 or np.linalg.matrix_rank(X)<X.shape[1]:raise ValueError('rank')
    Q,_=np.linalg.qr(N,mode='reduced');xr=X[:,-1]-Q@(Q.T@X[:,-1])
    beta=float(xr@Y/(xr@xr));res=Y-X@np.linalg.lstsq(X,Y,rcond=None)[0]
    se=np.sqrt((res@res)/(len(y)-X.shape[1])/(xr@xr))
    if se<1e-12:raise ValueError('variance')
    return beta,float(beta/se),X,Y


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inputs=root/'HNSCC_frozen_MCP_clinical_preparation_20261014_v1'
    out=root/'HNSCC_fixed_MCP_clinical_results_20261014_v2'
    result=pd.read_csv(out/'all24_private_clinical_statistics.tsv',sep='\t')
    assert len(result)==24 and not result.duplicated(['program','endpoint','method']).any()
    scores=pd.read_csv(inputs/'private_patient_library_frozen_scores.tsv',sep='\t')
    md=pd.read_csv(inputs/'private_eligible_libraries_clinical.tsv',sep='\t')
    b=md[md.timepoint.eq('pre')].set_index('patient_id')
    s=scores[scores.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='score')
    z=scores[scores.timepoint.eq('pre')].pivot(index='patient_id',columns='endpoint',values='CRC_training_z')
    s['equal_mean_CRC_training_z']=z.mean(axis=1);b=b.loc[s.index]
    assert b.index.is_unique and b.treatment_x.eq(b.treatment_y).all()
    N=np.column_stack((np.ones(len(b)),b.treatment_x.eq('D+T').astype(float),b.HPV.eq('positive').astype(float)))
    assert set(b.HPV)=={'positive','negative'}
    checks=[]
    for row in result.itertuples():
        if row.program=='MCP2':assert pd.isna(row.p) and row.status=='unsupported_original_program_slot';continue
        x=s[row.endpoint].to_numpy();y=b.regression.to_numpy(float)
        if row.method=='primary_adjusted_continuous':
            beta,t,X,Y=coefficient(x,y,N)
            assert np.isclose(beta,row.beta,rtol=1e-10,atol=1e-10)
            assert np.isclose(t,row.t,rtol=1e-10,atol=1e-10)
            Q,_=np.linalg.qr(N,mode='reduced');reduced=Q@(Q.T@Y);residual=Y-reduced
            QX,_=np.linalg.qr(X,mode='reduced');xr=X[:,-1]-Q@(Q.T@X[:,-1]);den=xr@xr
            rng=np.random.default_rng(20261014+int(hashlib.sha256((row.method+row.endpoint).encode()).hexdigest()[:8],16))
            extreme=0
            for _ in range(20000):
                yp=reduced+residual[rng.permutation(len(Y))]
                p_beta=xr@yp/den;r=yp-QX@(QX.T@yp)
                se=np.sqrt(max(float(r@r),0)/(len(Y)-X.shape[1])/den)
                p_t=abs(p_beta/se) if se>1e-12 else np.inf
                extreme+=p_t>=abs(t)-1e-12
            assert np.isclose((1+extreme)/20001,row.p,rtol=0,atol=1e-14)
            boots=[]
            for _ in range(2000):
                ix=rng.integers(0,len(x),len(x))
                try:boots.append(coefficient(x[ix],y[ix],N[ix])[0])
                except (ValueError,np.linalg.LinAlgError):pass
            assert len(boots)==row.valid_bootstrap
            if len(boots)>=1800:
                ci=np.quantile(boots,[.025,.975]);assert np.allclose(ci,[row.lower,row.upper],atol=1e-9,rtol=1e-9)
            checks.append(dict(endpoint=row.endpoint,primary_beta_t_all_permutations_and_bootstrap=True))
        elif pd.notna(row.p):
            if row.method=='secondary_unadjusted_continuous':effect=float(np.corrcoef(ranks(x),ranks(y))[0,1])
            else:
                if row.method=='sensitivity_author_concordant_binary':
                    mask=b.author_binary_concordant.to_numpy(bool);x=x[mask];y=b.derived_response.to_numpy(int)[mask]
                else:y=b.derived_response.to_numpy(int)
                effect=sum(int(a>c)-int(a<c) for a in x[y==1] for c in x[y==0])/int((y==1).sum()*(y==0).sum())
            assert np.isclose(effect,row.effect,atol=1e-12)
    for _,family in result.groupby('method'):
        assert len(family)==6
        values=family.p.fillna(1).to_numpy()
        expected=[]
        for value in values:
            expected.append(min(1,min(p*6/(values<=p).sum() for p in values if p>=value)))
        mask=family.p.notna().to_numpy()
        assert np.allclose(np.asarray(expected)[mask],family.BH_six_slot_family.to_numpy()[mask])
    report=dict(status='PASS_REAL_PRIMARY_ALL_PERMUTATIONS_BOOTSTRAP_SECONDARY_POINTS_AND_ALL_BH',
      finished_UTC=datetime.now(timezone.utc).isoformat(),checks=checks,rows=24,
      unsupported_original_slots_kept=True,result_sha256=hashlib.sha256((out/'all24_private_clinical_statistics.tsv').read_bytes()).hexdigest())
    assert len(checks)==3
    (out/'private_independent_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))
    print(result[['program','endpoint','method','status','n_patients','beta','effect','lower','upper','p','BH_six_slot_family']].to_json(orient='records'))


if __name__=='__main__':main()
