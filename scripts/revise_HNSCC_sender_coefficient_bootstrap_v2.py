"""Repair bootstrap estimand eligibility; preserve every original primary P/BH."""
import hashlib,json,os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import rankdata
from threadpoolctl import threadpool_limits


def coefficient(data):
    n=len(data)
    if n<20 or n-14<10 or not np.isfinite(data).all():raise ValueError('n_df_missing')
    r=rankdata(np.column_stack((data[:,:2],data[:,4:])),axis=0,method='average');sd=r.std(axis=0)
    if (sd<1e-12).any():raise ValueError('constant_required_field')
    r=(r-r.mean(axis=0))/sd;x,y=r[:,0],r[:,1]
    N=np.column_stack((np.ones(n),data[:,2:4],r[:,2:]));D=np.column_stack((N,x))
    if np.linalg.matrix_rank(D)!=14 or np.linalg.cond(D)>1e8:raise ValueError('rank_condition')
    Q,_=np.linalg.qr(N,mode='reduced');rx=x-Q@(Q.T@x);ry=y-Q@(Q.T@y)
    xx=float(rx@rx)
    if xx<1e-12:raise ValueError('zero_residual_exposure')
    return float(rx@ry/xx)


def describe(data):
    try:b=coefficient(data)
    except (ValueError,np.linalg.LinAlgError):return dict(beta=np.nan,LOO_valid=0,LOO_min=np.nan,LOO_max=np.nan,status='unestimable_full_control_descriptor')
    values=[]
    for i in range(len(data)):
        try:values.append(coefficient(np.delete(data,i,axis=0)))
        except (ValueError,np.linalg.LinAlgError):pass
    return dict(beta=b,LOO_valid=len(values),LOO_min=min(values) if values else np.nan,LOO_max=max(values) if values else np.nan,status='estimated_descriptive_no_new_P')


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);old=root/'HNSCC_frozen_sender_association_results_20261016_v1'
    prep=root/'HNSCC_frozen_sender_association_preflight_20261016_v1';out=root/'HNSCC_frozen_sender_association_results_20261016_v2'
    assert not out.exists();out.mkdir()
    original=json.loads((old/'private_terminal.json').read_text());old_path=old/'all400_private_primary_association_statistics.tsv'
    assert hashlib.sha256(old_path.read_bytes()).hexdigest()==original['result_sha256']
    result=pd.read_csv(old_path,sep='\t');fixed=result[['ligand','sender_group','beta','t','p','BH']].copy()
    config=json.loads((prep/'contract_before_any_association_effects.json').read_text());data_path=prep/'private_all_planned_patient_designs.tsv'
    accepted=json.loads((prep/'private_independent_audit.json').read_text());assert hashlib.sha256(data_path.read_bytes()).hexdigest()==accepted['design_sha256']
    frames=pd.read_csv(data_path,sep='\t');fields=['exposure','response']+config['categorical_nuisance']+config['numeric_nuisance']
    inputs=['private_original_receiver_gene_contributions.tsv','private_all_interface_contribution_removals.tsv']
    for name in inputs:os.link(old/name,out/name)
    contract=dict(created_UTC=datetime.now(timezone.utc).isoformat(),producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        correction='Coefficient bootstrap/LOO eligibility does not require a nonzero regression-error variance needed only for studentized t. Keep all original coverage/rank/condition/zero-residual-exposure guards.',
        original_primary_result_sha256=original['result_sha256'],original_descriptor_sha256=original['descriptor_sha256'],
        old_version_and_failed_audit_preserved=True,primary_points_t_P_BH_unchanged=True,all400_family_and_covariates_unchanged=True,
        randomization='Same keyed RNG stream, advance every original20000 stratum permutation draw, then the same2000 bootstrap index draws; no seed changes.',
        interval='Coefficient-only full-rank bootstrap, >=1800 valid coefficients, same percentile interval; no new P family.',
        native_interface_inputs={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in inputs})
    (out/'contract_before_bootstrap_revision.json').write_text(json.dumps(contract,indent=2)+'\n')
    removal=pd.read_csv(out/'private_all_interface_contribution_removals.tsv',sep='\t').set_index(['ligand','patient_id'])
    prior_desc=pd.read_csv(old/'all800_private_control_and_interface_descriptors.tsv',sep='\t').set_index(['ligand','sender_group','variant'])
    descriptors=[]
    with threadpool_limits(limits=1):
        for i,row in enumerate(result.itertuples()):
            frame=frames[frames.ligand.eq(row.ligand)&frames.sender_group.eq(row.sender_group)&frames.included_in_full_design].sort_values('patient_id');data=frame[fields].to_numpy(float)
            if pd.notna(row.p):
                assert np.isclose(coefficient(data),row.beta,atol=1e-10,rtol=1e-9)
                key=row.ligand+'|'+row.sender_group;rng=np.random.default_rng(20261016+int(hashlib.sha256(key.encode()).hexdigest()[:8],16))
                strata=data[:,2]+2*data[:,3];groups=[np.flatnonzero(strata==k) for k in sorted(set(strata))]
                for _ in range(80):
                    for group in groups:
                        for _ in range(250):rng.permutation(group)
                draws=[]
                for _ in range(2000):
                    ix=np.concatenate([rng.choice(group,len(group),replace=True) for group in groups])
                    try:draws.append(coefficient(data[ix]))
                    except (ValueError,np.linalg.LinAlgError):pass
                lo,hi=np.quantile(draws,[.025,.975]) if len(draws)>=1800 else (np.nan,np.nan)
                result.loc[i,['lower','upper','valid_bootstrap']]=[lo,hi,len(draws)]
                for variant in ['primary_full_control','interface_disjoint_full_control']:
                    altered=data.copy()
                    if variant.startswith('interface_'):
                        altered[:,1]-=removal.loc[[(row.ligand,p) for p in frame.patient_id],'removed_original_contribution'].to_numpy()
                    desc=prior_desc.loc[(row.ligand,row.sender_group,variant)].to_dict();desc.update(describe(altered))
                    descriptors.append(dict(ligand=row.ligand,sender_group=row.sender_group,variant=variant,**desc))
            else:
                for variant in ['primary_full_control','interface_disjoint_full_control']:
                    descriptors.append(dict(ligand=row.ligand,sender_group=row.sender_group,variant=variant,**prior_desc.loc[(row.ligand,row.sender_group,variant)].to_dict()))
            (out/'private_progress.json').write_text(json.dumps(dict(completed_slots=i+1,observed_UTC=datetime.now(timezone.utc).isoformat())))
    assert result[['ligand','sender_group','beta','t','p','BH']].equals(fixed)
    assert len(result)==400 and len(descriptors)==800
    result.to_csv(out/'all400_private_primary_association_statistics.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(descriptors).to_csv(out/'all800_private_control_and_interface_descriptors.tsv',sep='\t',index=False,na_rep='NA')
    report=dict(status='ALL400_PRIMARY_UNCHANGED_COEFFICIENT_BOOTSTRAP_AND800_DESCRIPTORS_REVISED_V2_PRIVATE',
        finished_UTC=datetime.now(timezone.utc).isoformat(),planned_family=400,estimated=int(result.p.notna().sum()),
        primary_BH_below05=int(result.BH.lt(.05).sum()),primary_beta_t_P_BH_unchanged=True,
        result_sha256=hashlib.sha256((out/'all400_private_primary_association_statistics.tsv').read_bytes()).hexdigest(),
        descriptor_sha256=hashlib.sha256((out/'all800_private_control_and_interface_descriptors.tsv').read_bytes()).hexdigest(),
        independent_audit_pending=True,direct_interaction_verified=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))


if __name__=='__main__':main()
