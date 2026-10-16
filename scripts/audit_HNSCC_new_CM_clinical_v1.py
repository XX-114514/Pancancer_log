"""Independent manual-rank/full-QR replay of all30 tests and120 descriptors."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd


def rankz(x):
    x=np.asarray(x,float)
    if not np.isfinite(x).all() or max(x)-min(x)<=1e-10*max(1,max(abs(x))):raise ValueError('constant')
    r=np.array([sum(x<v)+(sum(x==v)+1)/2 for v in x],float)
    return (r-r.mean())/r.std()


def fit(x,y,N):
    D=np.column_stack((N,rankz(x)));Y=rankz(y);df=len(y)-D.shape[1]
    if len(y)<20 or df<10 or np.linalg.matrix_rank(D)!=D.shape[1]:raise ValueError('rank_or_coverage')
    Q,R=np.linalg.qr(D,mode='reduced');b=np.linalg.solve(R,Q.T@Y);e=Y-D@b
    v=np.linalg.inv(D.T@D)[-1,-1];se=np.sqrt(float(e@e)/df*v)
    if se<1e-12:raise ValueError('variance')
    return float(b[-1]),float(b[-1]/se),D,Y,Q,R,float(v)


def close(a,b):assert np.isclose(a,b,atol=1e-9,rtol=2e-7),(a,b)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);inputs=root/'HNSCC_new_CM_baseline_preparation_20261016_v1'
    out=root/'HNSCC_new_CM_clinical_results_20261016_v1';report_path=out/'private_independent_audit.json';assert not report_path.exists()
    source=json.loads((inputs/'private_independent_input_audit.json').read_text());assert source['status'].startswith('PASS_ALL30')
    scores=pd.read_csv(inputs/'private_patient_baseline_component_scores.tsv',sep='\t')
    clinical=pd.read_csv(inputs/'private_patient_clinical_covariates.tsv',sep='\t')
    result=pd.read_csv(out/'all30_private_clinical_statistics.tsv',sep='\t')
    descriptors=pd.read_csv(out/'all120_private_descriptive_and_LOO.tsv',sep='\t')
    assert len(result)==30 and len(descriptors)==120
    for row in result.itertuples():
        f=scores[scores.scope.eq(row.scope)&scores.component.eq(row.component)].merge(clinical[clinical.scope.eq(row.scope)],on=['scope','patient_id'],validate='one_to_one').sort_values('patient_id')
        f=f[f[['D_plus_T','HPV_positive','regression']].notna().all(axis=1)]
        assert len(f)==row.n_patients and f.patient_id.is_unique
        x=f.activity.to_numpy();y=f.regression.to_numpy();N=np.column_stack((np.ones(len(f)),f.D_plus_T,f.HPV_positive))
        try:b,t,D,Y,Q,R,v=fit(x,y,N)
        except (ValueError,np.linalg.LinAlgError):
            assert pd.isna(row.p)
            continue
        close(b,row.beta);close(t,row.t);df=len(y)-D.shape[1];assert df==row.residual_df
        reduced=N@np.linalg.lstsq(N,Y,rcond=None)[0];residual=Y-reduced
        key=row.scope+'|'+row.component;rng=np.random.default_rng(20261016+int(hashlib.sha256(key.encode()).hexdigest()[:8],16))
        P=np.column_stack([reduced+residual[rng.permutation(len(y))] for _ in range(20000)])
        pb=np.linalg.solve(R,Q.T@P);error=P-D@pb;se=np.sqrt(np.sum(error*error,axis=0)/df*v)
        pt=np.divide(pb[-1],se,out=np.full(20000,np.inf),where=se>1e-12)
        close((1+np.count_nonzero(abs(pt)>=abs(t)-1e-12))/20001,row.p)
        boots=[]
        for _ in range(2000):
            ix=rng.integers(0,len(y),len(y))
            try:boots.append(fit(x[ix],y[ix],N[ix])[0])
            except (ValueError,np.linalg.LinAlgError):pass
        assert len(boots)==row.valid_bootstrap
        if len(boots)>=1800:
            lo,hi=np.quantile(boots,[.025,.975]);close(lo,row.lower);close(hi,row.upper)
        else:assert pd.isna(row.lower) and pd.isna(row.upper)
        sub=descriptors[descriptors.scope.eq(row.scope)&descriptors.component.eq(row.component)]
        assert set(sub.variant)=={'primary','unadjusted','raw_W_only','age_sex'}
        for diagnostic in sub.itertuples():
            g=f;z=x;C=N
            if diagnostic.variant=='unadjusted':C=np.ones((len(f),1))
            elif diagnostic.variant=='raw_W_only':z=f.raw_W_activity.to_numpy()
            elif diagnostic.variant=='age_sex':
                g=f[f[['age','male']].notna().all(axis=1)];z=g.activity.to_numpy()
                # The producer uses fixed full-cohort ranked age in LOO; this is
                # descriptive sensitivity, not heldout prediction or a new test.
                C=np.column_stack((np.ones(len(g)),g.D_plus_T,g.HPV_positive,g.male,rankz(g.age.to_numpy())))
            yy=g.regression.to_numpy();assert len(g)==diagnostic.n_patients
            try:close(fit(z,yy,C)[0],diagnostic.beta)
            except (ValueError,np.linalg.LinAlgError):assert pd.isna(diagnostic.beta)
            loo=[]
            for i in range(len(g)):
                keep=np.arange(len(g))!=i
                try:loo.append(fit(z[keep],yy[keep],C[keep])[0])
                except (ValueError,np.linalg.LinAlgError):pass
            assert len(loo)==diagnostic.LOO_valid
            if loo:close(min(loo),diagnostic.LOO_min);close(max(loo),diagnostic.LOO_max)
    p=result.p.fillna(1).to_numpy()
    q=[min(1,min(v*len(p)/(p<=v).sum() for v in p if v>=own)) for own in p]
    assert np.allclose(np.asarray(q)[result.p.notna()],result.loc[result.p.notna(),'BH_all30'],atol=1e-12,rtol=1e-12)
    report=dict(status='PASS_ALL30_CLINICAL_TESTS_FULL_RESAMPLING_BH_AND_ALL120_CONTROL_LOO',
        finished_UTC=datetime.now(timezone.utc).isoformat(),point_engine='independent manual ranks and full QR',
        result_sha256=hashlib.sha256((out/'all30_private_clinical_statistics.tsv').read_bytes()).hexdigest(),
        original_20261010_9CM_changed=False,direct_interaction_verified=False)
    report_path.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
