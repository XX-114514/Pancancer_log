"""Independent full-group numerical replay and complete-family correction."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from audit_HNSCC_frozen_T_sender_associations_v1 import ranked,nuisance,same


def replay(f,frame,md,cov,seed):
    n=len(md);N=nuisance(cov);df=n-N.shape[1]-1;Y=ranked(md.receiver_programme.to_numpy())
    assert n>=20 and df>=10 and np.linalg.matrix_rank(N)==N.shape[1]
    cols=[ids.split('|')[0] for ids in f.programme_ids];V=frame[cols].to_numpy();X=ranked(V)
    fit=N@np.linalg.lstsq(N,Y,rcond=None)[0];e=Y-fit;rng=np.random.default_rng(seed)
    perm=np.column_stack([fit+e[rng.permutation(n)] for _ in range(20000)])
    boots=[[] for _ in range(len(f))]
    for _ in range(2000):
        ix=rng.integers(0,n,n);BN=nuisance(cov[ix])
        if np.linalg.matrix_rank(BN)<BN.shape[1]:continue
        BX=ranked(V[ix]);BY=ranked(md.receiver_programme.to_numpy()[ix])
        if BY.std()<1e-12:continue
        RX=BX-BN@np.linalg.lstsq(BN,BX,rcond=None)[0];RY=BY-BN@np.linalg.lstsq(BN,BY,rcond=None)[0]
        norm=np.sum(RX*RX,axis=0)
        for j in range(len(f)):
            if BX[:,j].std()>1e-12 and norm[j]>1e-10:boots[j].append(float(RX[:,j]@RY/norm[j]))
    for j,row in enumerate(f.itertuples()):
        if pd.isna(row.p):continue
        full=np.column_stack((N,X[:,j]));Q,R=np.linalg.qr(full,mode='reduced')
        b=np.linalg.solve(R,Q.T@Y);res=Y-full@b;variance=np.linalg.inv(full.T@full)[-1,-1]
        t=b[-1]/np.sqrt((res@res)/df*variance);same(b[-1],row.beta);same(t,row.t)
        pb=np.linalg.solve(R,Q.T@perm);pe=perm-full@pb;pse=np.sqrt(np.maximum(np.sum(pe*pe,axis=0)/df*variance,0))
        pt=np.divide(pb[-1],pse,out=np.full(20000,np.inf),where=pse>1e-12)
        p=(1+np.count_nonzero(np.abs(pt)>=abs(t)-1e-12))/20001
        assert np.isclose(p,row.p,rtol=0,atol=1e-14)
        assert len(boots[j])==row.valid_bootstrap
        if len(boots[j])>=1800:
            lo,hi=np.quantile(boots[j],[.025,.975]);same(lo,row.lower);same(hi,row.upper)
        else:assert pd.isna(row.lower) and pd.isna(row.upper)
        for pid in row.programme_ids.split('|'):assert np.allclose(frame[pid],frame[cols[j]],atol=1e-12)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);inputs=root/'HNSCC_all_T_group_sender_states_20261014_v1'
    out=root/'HNSCC_all_T_group_sender_associations_20261014_v1'
    result=pd.read_csv(out/'all576_private_T_group_associations.tsv',sep='\t')
    oldroot=root/'HNSCC_frozen_T_sender_associations_20261014_v1'
    assert json.loads((oldroot/'private_independent_audit.json').read_text())['status'].startswith('PASS_ALL192')
    old=pd.read_csv(oldroot/'all_private_sender_state_associations.tsv',sep='\t')
    cd8=result[result.T_candidate_group.eq('CD8_Effector')].set_index(['model','signature_group']).sort_index()
    original=old.set_index(['model','signature_group']).sort_index()
    for c in ['beta','t','p','lower','upper','valid_bootstrap','BH_admissible_unique_family','BH_all_candidate_unique_family']:
        assert np.allclose(cd8[c],original[c],equal_nan=True,rtol=1e-12,atol=1e-12)
    columns=['T_logRNA','T_loggenes','Mac_logRNA','Mac_loggenes','T_fraction','Mac_fraction','T_logcells','Mac_logcells']
    for group in sorted(set(result.T_candidate_group)-{'CD8_Effector'}):
        directory=inputs/group;md=pd.read_csv(directory/'private_patient_receiver_covariates.tsv',sep='\t',index_col=0)
        raw=pd.read_csv(directory/'private_patient_T_state_means.tsv',sep='\t',index_col=0).loc[md.index]
        orth=pd.read_csv(directory/'private_patient_T_orthogonal_state_means.tsv',sep='\t',index_col=0).loc[md.index]
        models=[('primary_RNA_capture',raw,md[columns].to_numpy(float)),('sensitivity_unadjusted',raw,np.empty((len(md),0))),
          ('sensitivity_demographic',raw,md[columns+['D_plus_T','HPV_positive','male','age']].to_numpy(float)),
          ('sensitivity_sender_gene_overlap_removed',orth,md[columns].to_numpy(float))]
        for model,frame,cov in models:
            f=result[result.T_candidate_group.eq(group)&result.model.eq(model)]
            assert len(f)==48
            seed=20261014+int(hashlib.sha256((group+'|'+model).encode()).hexdigest()[:8],16)
            replay(f,frame,md,cov,seed)
    summaries=[]
    for model,f in result.groupby('model'):
        assert len(f)==144
        for field,mask in [('BH_all_T_groups_all_candidate_family',np.ones(len(f),bool)),('BH_all_T_groups_admissible_family',f.primary_admissible.to_numpy(bool))]:
            p=f.loc[mask,'p'].fillna(1).to_numpy()
            q=[min(1,min(value*len(p)/(p<=value).sum() for value in p if value>=own)) for own in p]
            valid=f.loc[mask,'p'].notna().to_numpy()
            assert np.allclose(np.asarray(q)[valid],f.loc[mask,field].to_numpy()[valid],rtol=1e-12,atol=1e-12)
        summaries.append(dict(model=model,unique_tests=144,admissible_unique_tests=int(f.primary_admissible.sum()),
          admissible_BH05=int((f.primary_admissible&f.BH_all_T_groups_admissible_family.lt(.05)).sum()),
          admissible_min_BH=float(f.loc[f.primary_admissible,'BH_all_T_groups_admissible_family'].min())))
    assert len(result)==576
    report=dict(status='PASS_ALL576_T_GROUP_ASSOCIATIONS_FULL_NUMERICAL_REPLAY_AND_GLOBAL_BH',finished_UTC=datetime.now(timezone.utc).isoformat(),
      families=summaries,original_CD8_unchanged=True,new_clinical_tests=False,
      result_sha256=hashlib.sha256((out/'all576_private_T_group_associations.tsv').read_bytes()).hexdigest())
    (out/'private_independent_audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
    primary=result[result.model.eq('primary_RNA_capture')&result.primary_admissible]
    print(primary.sort_values('p').head(8)[['T_candidate_group','programme_ids','beta','lower','upper','p','BH_all_T_groups_admissible_family','valid_bootstrap']].to_json(orient='records'))


if __name__=='__main__':main()
