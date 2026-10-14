"""Independent OLS, shared permutation, bootstrap and BH replay."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from audit_HNSCC_fixed_MCP_clinical_v1 import ranks


def ranked(x):
    if np.asarray(x).ndim==1:
        r=ranks(x);sd=r.std();return (r-r.mean())/sd if sd>1e-12 else r*0
    return np.column_stack([ranked(x[:,i]) for i in range(x.shape[1])])


def nuisance(cov):return np.column_stack((np.ones(len(cov)),*[ranks(cov[:,i]) for i in range(cov.shape[1])]))


def same(x,y):assert np.isclose(x,y,atol=1e-8,rtol=1e-8,equal_nan=True),(x,y)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);inputs=root/'HNSCC_frozen_T_sender_states_20261014_v1'
    out=root/'HNSCC_frozen_T_sender_associations_20261014_v1'
    flags=json.loads((root/'HNSCC_T_sender_signature_flag_audit_20261014_v1/private_flag_audit.json').read_text())
    assert all(not x['admissible_flag_conflicts'] for x in flags['groupings'])
    result=pd.read_csv(out/'all_private_sender_state_associations.tsv',sep='\t')
    candidates=pd.read_csv(out/'private_signature_dedup_contract.tsv',sep='\t')
    md=pd.read_csv(inputs/'private_patient_receiver_covariates.tsv',sep='\t',index_col=0)
    raw=pd.read_csv(inputs/'private_patient_T_state_means.tsv',sep='\t',index_col=0).loc[md.index]
    orth=pd.read_csv(inputs/'private_patient_T_orthogonal_state_means.tsv',sep='\t',index_col=0).loc[md.index]
    maincols=['T_logRNA','T_loggenes','Mac_logRNA','Mac_loggenes','T_fraction','Mac_fraction','T_logcells','Mac_logcells']
    models=[('primary_RNA_capture',raw,md[maincols].to_numpy(float)),('sensitivity_unadjusted',raw,np.empty((len(md),0))),
      ('sensitivity_demographic',raw,md[maincols+['D_plus_T','HPV_positive','male','age']].to_numpy(float)),
      ('sensitivity_sender_gene_overlap_removed',orth,md[maincols].to_numpy(float))]
    n=len(md);Y=ranked(md.receiver_programme.to_numpy());checks=[]
    for index,(name,frame,cov) in enumerate(models):
        f=result[result.model.eq(name)].copy();assert len(f)==48
        N=nuisance(cov);df=n-N.shape[1]-1
        if np.linalg.matrix_rank(N)<N.shape[1] or df<10:
            assert f.p.isna().all();continue
        fitted=N@np.linalg.lstsq(N,Y,rcond=None)[0];residual=Y-fitted
        rng=np.random.default_rng(20261014+index)
        permuted=np.column_stack([fitted+residual[rng.permutation(n)] for _ in range(20000)])
        columns=[ids.split('|')[0] for ids in f.programme_ids]
        V=frame[columns].to_numpy();X=ranked(V)
        bootstrap_indices=[rng.integers(0,n,n) for _ in range(2000)]
        boots=[[] for _ in range(len(f))]
        for ix in bootstrap_indices:
            BN=nuisance(cov[ix])
            if np.linalg.matrix_rank(BN)<BN.shape[1]:continue
            BX=ranked(V[ix]);BY=ranked(md.receiver_programme.to_numpy()[ix])
            if np.std(BY)<1e-12:continue
            BR=BX-BN@np.linalg.lstsq(BN,BX,rcond=None)[0];BE=BY-BN@np.linalg.lstsq(BN,BY,rcond=None)[0]
            norm=np.sum(BR*BR,axis=0)
            for j in range(len(f)):
                if np.std(BX[:,j])>1e-12 and norm[j]>1e-10:boots[j].append(float(BR[:,j]@BE/norm[j]))
        for j,row in enumerate(f.itertuples()):
            assert df==row.residual_df
            if pd.isna(row.p):continue
            full=np.column_stack((N,X[:,j]));assert np.linalg.matrix_rank(full)==full.shape[1]
            Q,R=np.linalg.qr(full,mode='reduced')
            b=np.linalg.solve(R,Q.T@Y);e=Y-full@b
            covdiag=np.linalg.inv(full.T@full)[-1,-1]
            se=np.sqrt(float(e@e)/df*covdiag);t=b[-1]/se
            same(b[-1],row.beta);same(t,row.t)
            pb=np.linalg.solve(R,Q.T@permuted);pe=permuted-full@pb
            pse=np.sqrt(np.maximum(np.sum(pe*pe,axis=0)/df*covdiag,0))
            pt=np.divide(pb[-1],pse,out=np.full(20000,np.inf),where=pse>1e-12)
            p=(1+np.count_nonzero(np.abs(pt)>=abs(t)-1e-12))/20001
            assert np.isclose(p,row.p,atol=1e-14,rtol=0)
            assert len(boots[j])==row.valid_bootstrap
            if len(boots[j])>=1800:
                ci=np.quantile(boots[j],[.025,.975]);same(ci[0],row.lower);same(ci[1],row.upper)
            else:assert pd.isna(row.lower) and pd.isna(row.upper)
            # All original IDs within one group really have identical scores.
            for pid in row.programme_ids.split('|'):assert np.allclose(frame[pid],frame[columns[j]],atol=1e-12)
        for field,mask in [('BH_all_candidate_unique_family',np.ones(len(f),bool)),('BH_admissible_unique_family',f.primary_admissible.to_numpy(bool))]:
            values=f.loc[mask,'p'].fillna(1).to_numpy();expected=[]
            for own in values:expected.append(min(1,min(p*len(values)/(values<=p).sum() for p in values if p>=own)))
            valid=f.loc[mask,'p'].notna().to_numpy()
            assert np.allclose(np.asarray(expected)[valid],f.loc[mask,field].to_numpy()[valid],atol=1e-12)
        checks.append(dict(model=name,groups=len(f),estimated_groups=int(f.p.notna().sum()),primary_admissible_groups=int(f.primary_admissible.sum()),
          all_permutations_bootstraps_coefficients_BH_replayed=True,
          admissible_BH05=int((f.primary_admissible&f.BH_admissible_unique_family.lt(.05)).sum())))
    assert len(result)==192 and len(checks)==4
    report=dict(status='PASS_ALL192_SENDER_MODEL_GROUPS_REAL_STATISTICS_AND_BH',finished_UTC=datetime.now(timezone.utc).isoformat(),checks=checks,
      result_sha256=hashlib.sha256((out/'all_private_sender_state_associations.tsv').read_bytes()).hexdigest(),new_clinical_tests=False)
    (out/'private_independent_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    primary=result[result.model.eq('primary_RNA_capture')&result.primary_admissible]
    print(primary[['programme_ids','status','beta','lower','upper','p','BH_admissible_unique_family','valid_bootstrap']].to_json(orient='records'))


if __name__=='__main__':main()
