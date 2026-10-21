"""Paired exploratory background sensitivity, no new P/BH or clinical tests.

Same eligible cases and complete original controls are used before/after.
The background is an a priori equal raw-score mean, not a fitted predictor.
Conditioning on shared biology is not a causal or technical-artifact test.
"""
from collections import Counter
import numpy as np
import RNA_programme_partial_rank_v1 as stats
import audit_RNA_programme_partial_rank_v1 as audit


def paired_fit(x,y,c,cat,b,m,levels,independent=False):
    joint=(np.asarray(b)+np.asarray(m))/2
    if independent:
        base=audit.reference_fit(x,y,c,cat,levels)[0]
        adjusted=audit.reference_fit(x,y,np.column_stack([c,joint]),cat,levels)[0]
    else:
        base=stats.fit(x,y,c,cat,levels)['rho']
        adjusted=stats.fit(x,y,np.column_stack([c,joint]),cat,levels)['rho']
    return np.array([base,adjusted,adjusted-base])


def analyze(x,y,c,cat,b,m,slot,draws=2000,required=1800,independent=False):
    x=np.asarray(x,float);y=np.asarray(y,float);c=np.asarray(c,float);cat=np.asarray(cat,str);b=np.asarray(b,float);m=np.asarray(m,float)
    assert len(x)==len(y)==len(c)==len(cat)==len(b)==len(m) and np.isfinite(b).all() and np.isfinite(m).all()
    output=dict(status='NOT_ESTIMABLE',n_cases=len(x),bootstrap_draws=draws,bootstrap_valid=0,no_new_P_or_BH=True)
    if len(x)<20:output['reason']='INSUFFICIENT_ELIGIBLE_CASES';return output
    levels=stats.categorical_spec(cat)
    try:point=paired_fit(x,y,c,cat,b,m,levels,independent)
    except ValueError as error:output['reason']=str(error);return output
    groups={}
    for i,key in enumerate(map(tuple,cat)):groups.setdefault(key,[]).append(i)
    strata=[np.array(groups[k],int) for k in sorted(groups)]
    rng=np.random.default_rng(np.random.SeedSequence([2026102206,int(slot)]));values=[];failure=Counter()
    for unused in range(draws):
        chosen=np.concatenate([rng.choice(g,size=len(g),replace=True) for g in strata])
        try:value=paired_fit(x[chosen],y[chosen],c[chosen],cat[chosen],b[chosen],m[chosen],levels,independent)
        except ValueError as error:failure[str(error)]+=1;continue
        values.append(value)
    intervals=np.quantile(values,[.025,.975],axis=0).T.tolist() if len(values)>=required else [None]*3
    output.update(status='DESCRIPTIVE_PAIRED_BACKGROUND_SENSITIVITY',reason=None,
        same_case_original_rho=float(point[0]),with_joint_background_rho=float(point[1]),delta_rho=float(point[2]),
        same_case_original_bootstrap95=intervals[0],joint_background_bootstrap95=intervals[1],delta_bootstrap95=intervals[2],
        bootstrap_valid=len(values),bootstrap_failure_reasons=dict(failure),stratum_sizes=sorted(map(len,strata)),
        no_independent_validation_or_causal_claim=True)
    if not independent:
        sensitivity={}
        for name,background in [('author_concordant_B_only',b),('author_concordant_myeloid_only',m)]:
            try:
                fit=stats.fit(x,y,np.column_stack([c,background]),cat,levels)
                reference=audit.reference_fit(x,y,np.column_stack([c,background]),cat,levels)[0]
                np.testing.assert_allclose(fit['rho'],reference,atol=1e-10,rtol=1e-10)
                sensitivity[name]=dict(status='DESCRIPTIVE_POINT_ONLY',rho=fit['rho'],condition=fit['condition'],residual_df=fit['residual_df'])
            except ValueError as error:sensitivity[name]=dict(status='NOT_ESTIMABLE',reason=str(error))
        fit=stats.fit(x,y,np.column_stack([c,(b+m)/2]),cat,levels)
        output.update(point_sensitivities=sensitivity,full_design_condition=fit['condition'],full_residual_df=fit['residual_df'],full_design_columns=fit['columns'])
    return output


def assert_replay(produced,independent):
    assert produced['status']==independent['status'] and produced.get('reason')==independent.get('reason')
    if produced['status']=='NOT_ESTIMABLE':return
    assert produced['bootstrap_valid']==independent['bootstrap_valid'] and produced['bootstrap_failure_reasons']==independent['bootstrap_failure_reasons']
    for name in ['same_case_original_rho','with_joint_background_rho','delta_rho','same_case_original_bootstrap95','joint_background_bootstrap95','delta_bootstrap95']:
        if produced[name] is None:assert independent[name] is None
        else:np.testing.assert_allclose(produced[name],independent[name],atol=1e-9,rtol=1e-9)
