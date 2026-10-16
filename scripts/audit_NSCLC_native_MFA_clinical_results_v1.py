"""Replay every actual endpoint random draw, CI, LOO and full pair family."""
import json
from datetime import datetime,timezone

import numpy as np
import pandas as pd
from scipy.linalg import qr
from scipy.special import betaincinv

from NSCLC_native_MFA_clinical_inputs_v1 import read,matched_quota,sha
from NSCLC_native_MFA_clinical_audit_math_v1 import fit,wild,bootstrap,loo,adjustment
from HNSCC_author_state_clinical_kernel_v2 import seed


def agent(y,a):
    z=np.column_stack([a['Zraw'],*(a['strata']==s for s in sorted(set(a['strata'])))])
    q,r,p=qr(z,mode='economic',pivoting=True);s=np.linalg.svd(z,compute_uv=False);rank=int((s>s[0]*1e-12).sum())
    if s[rank-1]/s[0]<1e-8:return np.nan
    q=q[:,:rank];rx=a['X']-q@(q.T@a['X']);den=rx@rx
    return float(rx@y/den) if den>1e-16*len(rx) and len(rx)-rank-1>=10 else np.nan


def main():
    root,pre,plan,scope,pairs,arrays=read();src=root/'NSCLC_native_MFA_clinical_results_20261016_v1'
    calproof=root/'NSCLC_native_MFA_clinical_calibration_audit_20261016_v1/private_terminal.json'
    cp=json.loads(calproof.read_text());assert cp['original_null_gate_passed'] is True
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    proof=src/'private_independent_clinical_result_audit.json';assert not proof.exists();matched_quota(root)
    singles=pd.read_csv(src/'private_all40_native_factor_clinical_endpoints.tsv',sep='\t').set_index('record',drop=False)
    joint=pd.read_csv(src/'private_all60_native_factor_clinical_pairs.tsv',sep='\t').set_index('slot',drop=False)
    assert len(singles)==40 and len(joint)==60
    with np.load(src/'private_all80000_bootstrap_and2480_LOO_values.npz',allow_pickle=False) as z:random={k:z[k] for k in z.files}
    assert np.array_equal(random['records'],scope.record.to_numpy(str))
    errors={};wildchecked=0;bootchecked=0;loochecked=0
    def close(name,a,b):
        a=np.asarray(a,float);b=np.asarray(b,float);assert a.shape==b.shape
        assert np.allclose(a,b,atol=1e-8,rtol=1e-8,equal_nan=True),name
        good=np.isfinite(a)&np.isfinite(b);v=float(np.max(abs(a[good]-b[good]))) if good.any() else 0.
        errors[name]=max(errors.get(name,0.),v)
    for i,row in enumerate(scope.itertuples()):
        observed=singles.loc[row.record];a=arrays[('internal_holdout',row.view)];k=row.factor_index;y=a['Y'][:,k];d=a['D'];g=fit(y,d)
        close('all40_HC3',[observed.beta,observed.HC3_t,observed.HC3_SE],[g['beta'],g['t'],g['se']])
        count=wild(g,plan['wild_draws'],seed(row.record,'wild'));wildchecked+=plan['wild_draws'];assert count==observed.wild_exceedances
        close('all40_P',observed.P,(count+1)/(plan['wild_draws']+1))
        upper=1. if count==plan['wild_draws'] else betaincinv(count+1,plan['wild_draws']-count,1-.01/40)
        upper=min(1.,(1+plan['wild_draws']*upper)/(1+plan['wild_draws']))
        close('all40_MC_upper',observed.MC_upper_P,upper)
        b=bootstrap(y,d,a['strata'],2000,seed(row.record,'bootstrap'));l=loo(y,d)
        close('all80000_bootstrap',random['bootstrap'][i],b);close('all2480_LOO',random['LOO'][i],l);bootchecked+=len(b);loochecked+=len(l)
        valid=int(np.isfinite(b).sum());assert valid==observed.valid_bootstraps==row.bootstrap_design_valid
        ci=np.quantile(b[np.isfinite(b)],[.025,.975]) if valid>=1900 else [np.nan,np.nan]
        close('all40_CI',[observed.CI_low,observed.CI_high],ci)
        assert int(np.isfinite(l).sum())==observed.LOO_valid
        assert bool(observed.LOO_same_sign)==bool(np.isfinite(l).all() and (l*g['beta']>0).all())
        close('all40_LOO_limits',[observed.LOO_min,observed.LOO_max],[np.nanmin(l),np.nanmax(l)])
        x=a['X'];n1=int((x==1).sum());n0=int((x==0).sum());ranks=pd.Series(y).rank(method='average').to_numpy();u=ranks[x==1].sum()-n1*(n1+1)/2
        close('all40_rank_biserial',observed.unadjusted_rank_biserial,2*u/(n1*n0)-1)
        close('all40_mixture_beta',observed.mixture_residual_beta,fit(a['Y_mixture'][:,k],d)['beta'])
        close('all40_agent_beta',observed.agent_controlled_beta,agent(y,a))
        for part in ['discovery','all_frozen_four_view_cases']:
            a2=arrays[(part,row.view)];close('all80_nonprimary_descriptive',observed[part+'_descriptive_beta'],fit(a2['Y'][:,k],a2['D'])['beta'])
    for table,name in [(singles,'endpoint40'),(joint,'pair60')]:
        close(name+'_BH',table.BH_q,adjustment(table.P.to_numpy()))
        close(name+'_BY',table.BY_q,adjustment(table.P.to_numpy(),by=True))
        close(name+'_MCupper_BH',table.MC_upper_BH_q,adjustment(table.MC_upper_P.to_numpy()))
    for row in pairs.itertuples():
        a=singles.loc[row.record_A];b=singles.loc[row.record_B];r=joint.loc[row.slot]
        close('all60_conjunction_P',r.P,max(a.P,b.P));close('all60_conjunction_MC',r.MC_upper_P,max(a.MC_upper_P,b.MC_upper_P))
        for label,expected in [('same_clinical_direction',a.beta*b.beta>0),('both_CI_exclude_zero',a.CI_low*a.CI_high>0 and b.CI_low*b.CI_high>0),('both_LOO_same_sign',a.LOO_same_sign and b.LOO_same_sign)]:assert bool(r[label])==bool(expected)
        robust=r.BH_q<.05 and r.MC_upper_BH_q<.05 and r.same_clinical_direction and r.both_CI_exclude_zero and r.both_LOO_same_sign
        assert bool(r.robust_internal_coassociation)==bool(robust)
        close('all60_pair_endpoints',[r.beta_A,r.beta_B,r.CI_A_low,r.CI_A_high,r.CI_B_low,r.CI_B_high,r.mixture_beta_A,r.mixture_beta_B,r.agent_beta_A,r.agent_beta_B],
            [a.beta,b.beta,a.CI_low,a.CI_high,b.CI_low,b.CI_high,a.mixture_residual_beta,b.mixture_residual_beta,a.agent_controlled_beta,b.agent_controlled_beta])
    result=dict(status='PASS_ALL40_NATIVE_CLINICAL_ENDPOINTS_799960_WILD_80000_BOOTSTRAP_2480_LOO_AND60_FULL_FAMILY_CORRECTIONS',finished_UTC=datetime.now(timezone.utc).isoformat(),
        audit_source_sha256=sha(__file__),producer_terminal_sha256=sha(src/'private_terminal.json'),calibration_audit_sha256=sha(calproof),
        wild_draws_checked=wildchecked,bootstrap_values_checked=bootchecked,LOO_values_checked=loochecked,maximum_absolute_errors=errors,
        all_planned_endpoints_and_pairs_retained=True,independent_patient_or_direct_interaction_claim=False)
    proof.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],finished_UTC=result['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
