"""Frozen all40 endpoint and60 conjunction tests, gated by full input/calibration."""
import json
from datetime import datetime,timezone

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_distribution,rankdata
from statsmodels.stats.multitest import multipletests

from NSCLC_native_MFA_clinical_inputs_v1 import read,matched_quota,sha,VIEWS
from HNSCC_author_state_clinical_kernel_v2 import fit,wild_count,bootstrap,loo,seed


def agent_contrast(y,a):
    z=np.column_stack([a['Zraw'],*(a['strata']==s for s in sorted(set(a['strata'])))])
    u,s,_=np.linalg.svd(z,full_matrices=False);rank=int((s>s[0]*1e-12).sum())
    if s[rank-1]/s[0]<1e-8:return np.nan
    q=u[:,:rank];r=a['X']-q@(q.T@a['X']);ss=r@r
    if ss<=1e-16*len(r) or len(r)-rank-1<10:return np.nan
    return float(r@y/ss)


def corrections(table,p='P',prefix=''):
    values=table[p].fillna(1).to_numpy()
    table[prefix+'BH_q']=multipletests(values,method='fdr_bh')[1]
    table[prefix+'BY_q']=multipletests(values,method='fdr_by')[1]


def main():
    root,pre,plan,scope,pairs,arrays=read();inputproof=pre/'private_independent_clinical_input_audit.json'
    assert json.loads(inputproof.read_text())['status'].startswith('PASS_')
    cal=root/'NSCLC_native_MFA_clinical_calibration_20261016_v1/private_terminal.json'
    assert json.loads(cal.read_text())['status']=='PASS_ALL12_NATIVE_CLINICAL_WORKING_MODEL_NULL_SCENARIOS'
    out=root/'NSCLC_native_MFA_clinical_results_20261016_v1';assert not out.exists();matched_quota(root);out.mkdir(mode=0o700)
    contract=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_audit_sha256=sha(inputproof),
        calibration_terminal_sha256=sha(cal),original_plan_sha256=sha(pre/'plan_before_any_native_factor_clinical_effect.json'),
        actual_family='All40 endpoints and60 original same-factor pairs; primary internal holdout only. Preserve all rows.',
        extra_pre_effect_diagnostic='For every endpoint, coefficient with complete original nuisance span PLUS all original histology x agent block indicators, on identical holdout cases. Descriptive only, no P or alternative primary selection.',
        clinical_prediction_or_LR_claim=False)
    (out/'plan_before_actual_native_factor_clinical_effects.json').write_text(json.dumps(contract,indent=2)+'\n')
    rows=[];boot=np.full((40,2000),np.nan);leave=np.full((40,62),np.nan)
    for i,row in enumerate(scope.itertuples()):
        a=arrays[('internal_holdout',row.view)];k=row.factor_index;y=a['Y'][:,k];d=a['D'];g=fit(y,d)
        count=wild_count(g,plan['wild_draws'],seed(row.record,'wild'));p=(1+count)/(1+plan['wild_draws'])
        upper=1. if count==plan['wild_draws'] else float(beta_distribution.ppf(1-.01/40,count+1,plan['wild_draws']-count))
        upper=min(1.,(1+plan['wild_draws']*upper)/(1+plan['wild_draws']))
        boot[i]=bootstrap(y,d,a['strata'],2000,seed(row.record,'bootstrap'));leave[i]=loo(y,d)
        valid=int(np.isfinite(boot[i]).sum());assert valid==row.bootstrap_design_valid
        ci=np.quantile(boot[i][np.isfinite(boot[i])],[.025,.975]) if valid>=1900 else [np.nan,np.nan]
        x=a['X'];n1=int((x==1).sum());n0=int((x==0).sum());u=rankdata(y)[x==1].sum()-n1*(n1+1)/2
        r=dict(record=row.record,factor=row.factor,view=row.view,n=len(y),status='ESTIMABLE',beta=g['beta'],HC3_t=g['t'],HC3_SE=float(np.sqrt(g['variance'])),
            P=p,MC_upper_P=upper,wild_exceedances=count,valid_bootstraps=valid,CI_low=float(ci[0]),CI_high=float(ci[1]),
            LOO_valid=int(np.isfinite(leave[i]).sum()),LOO_same_sign=bool(np.isfinite(leave[i]).all() and (leave[i]*g['beta']>0).all()),
            LOO_min=float(np.nanmin(leave[i])),LOO_max=float(np.nanmax(leave[i])),unadjusted_rank_biserial=float(2*u/(n1*n0)-1),
            mixture_residual_beta=fit(a['Y_mixture'][:,k],d,hc3=False)['beta'],agent_controlled_beta=agent_contrast(y,a))
        for part in ['discovery','all_frozen_four_view_cases']:
            b=arrays[(part,row.view)];r[part+'_descriptive_beta']=fit(b['Y'][:,k],b['D'],hc3=False)['beta']
        rows.append(r)
        (out/'private_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),completed_models=i+1,total_models=40))+'\n')
    singles=pd.DataFrame(rows);corrections(singles);singles['MC_upper_BH_q']=multipletests(singles.MC_upper_P,method='fdr_bh')[1]
    lookup=singles.set_index('record');joined=[]
    for row in pairs.itertuples():
        a=lookup.loc[row.record_A];b=lookup.loc[row.record_B]
        ci=(a.CI_low*a.CI_high>0)&(b.CI_low*b.CI_high>0)
        joined.append(dict(slot=row.slot,factor=row.factor,view_A=row.view_A,view_B=row.view_B,
            beta_A=a.beta,beta_B=b.beta,P=max(a.P,b.P),MC_upper_P=max(a.MC_upper_P,b.MC_upper_P),
            same_clinical_direction=bool(a.beta*b.beta>0),both_CI_exclude_zero=bool(ci),both_LOO_same_sign=bool(a.LOO_same_sign and b.LOO_same_sign),
            CI_A_low=a.CI_low,CI_A_high=a.CI_high,CI_B_low=b.CI_low,CI_B_high=b.CI_high,
            mixture_beta_A=a.mixture_residual_beta,mixture_beta_B=b.mixture_residual_beta,
            agent_beta_A=a.agent_controlled_beta,agent_beta_B=b.agent_controlled_beta))
    joint=pd.DataFrame(joined);corrections(joint);joint['MC_upper_BH_q']=multipletests(joint.MC_upper_P,method='fdr_bh')[1]
    joint['robust_internal_coassociation']=joint.BH_q.lt(.05)&joint.MC_upper_BH_q.lt(.05)&joint.same_clinical_direction&joint.both_CI_exclude_zero&joint.both_LOO_same_sign
    singles.to_csv(out/'private_all40_native_factor_clinical_endpoints.tsv',sep='\t',index=False)
    joint.to_csv(out/'private_all60_native_factor_clinical_pairs.tsv',sep='\t',index=False)
    np.savez_compressed(out/'private_all80000_bootstrap_and2480_LOO_values.npz',bootstrap=boot,LOO=leave,records=scope.record.to_numpy(str))
    terminal=dict(status='COMPLETE_ALL40_NATIVE_FACTOR_CLINICAL_AND60_PAIR_FAMILIES_FULL_RANDOM_VALUES_AUDIT_PENDING',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},
        unchanged_prior_models=True,independent_patient_or_direct_interaction_claim=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
