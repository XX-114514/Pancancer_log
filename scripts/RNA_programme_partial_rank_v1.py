"""Fixed partial-rank association, case bootstrap and residual permutation.

All arrays are caller-owned RAM values. No patient rows or random vectors are
saved. Categorical levels/strata are fixed before resampling; continuous
covariates are never dropped to make an otherwise singular design estimable.
Freedman-Lane is a conditional residual-exchange approximation, not a causal
test or proof that all unmeasured source/capture effects have been removed.
"""
from itertools import product
import numpy as np


def fixed_slots():
    keys=['left_parent','left_programme','right_programme','phase','context']
    dimensions=[['author_regulatory_reference','author_conventional_CD4_reference'],
                ['s.genes','g2m.genes'],['s.genes','g2m.genes'],['pre','post'],
                ['current_T_author_reference','CM02_network_member_union_author_reference']]
    return [dict(zip(keys,values),right_parent='author_CD8_reference') for values in product(*dimensions)]


def average_rank(value):
    x=np.asarray(value,float)
    assert x.ndim==1 and np.isfinite(x).all()
    order=np.argsort(x,kind='stable');v=x[order]
    starts=np.r_[0,np.flatnonzero(v[1:]!=v[:-1])+1];ends=np.r_[starts[1:],len(x)]
    result=np.empty(len(x),float)
    for start,end in zip(starts,ends):result[order[start:end]]=(start+1+end)/2
    return result


def correlation(x,y):
    x=np.asarray(x,float);y=np.asarray(y,float)
    x=x-x.mean();y=y-y.mean();scale=np.linalg.norm(x)*np.linalg.norm(y)
    if scale<=1e-12:return None
    return float(np.clip(x@y/scale,-1,1))


def categorical_spec(categories):
    c=np.asarray(categories,str);assert c.ndim==2
    return [sorted(set(c[:,j])) for j in range(c.shape[1])]


def design_matrix(continuous,categories,levels):
    c=np.asarray(continuous,float);cat=np.asarray(categories,str)
    assert c.ndim==2 and cat.ndim==2 and len(c)==len(cat) and np.isfinite(c).all()
    assert len(levels)==cat.shape[1]
    columns=[np.ones(len(c))];constant=[]
    for j in range(c.shape[1]):
        r=average_rank(c[:,j]);sd=float(r.std())
        if sd<=1e-12:raise ValueError('CONSTANT_CONTINUOUS_CONTROL')
        columns.append((r-r.mean())/sd)
    for j,allowed in enumerate(levels):
        assert set(cat[:,j])<=set(allowed) and len(allowed)>0
        if len(allowed)==1:constant.append(j)
        for name in allowed[1:]:columns.append((cat[:,j]==name).astype(float))
    return np.column_stack(columns),constant


def fit(x,y,continuous,categories,levels,min_df=10,max_condition=1e8):
    x=np.asarray(x,float);y=np.asarray(y,float)
    if not np.isfinite(x).all() or not np.isfinite(y).all():raise ValueError('NONFINITE_SCORE')
    if len(np.unique(x))<2 or len(np.unique(y))<2:raise ValueError('CONSTANT_PROGRAMME_SCORE')
    z,constant=design_matrix(continuous,categories,levels)
    rank=int(np.linalg.matrix_rank(z));condition=float(np.linalg.cond(z));df=len(x)-rank
    if rank<z.shape[1]:raise ValueError('SINGULAR_FULL_DESIGN')
    if condition>max_condition:raise ValueError('ILL_CONDITIONED_FULL_DESIGN')
    if df<min_df:raise ValueError('INSUFFICIENT_RESIDUAL_DF')
    q=np.linalg.qr(z,mode='reduced')[0];rx=average_rank(x);ry=average_rank(y)
    ex=rx-q@(q.T@rx);ey=ry-q@(q.T@ry)
    if min(np.linalg.norm(ex),np.linalg.norm(ey))<=1e-10:raise ValueError('NO_RESIDUAL_SCORE_VARIATION')
    rho=correlation(ex,ey)
    assert rho is not None
    return dict(rho=rho,unadjusted=correlation(rx,ry),rank=rank,columns=z.shape[1],
                condition=condition,residual_df=df,constant_factors=constant,
                q=q,ex=ex,ey=ey,yhat=ry-ey,z=z)


def permuted_statistic(fitted,permutation):
    """Permute ranked-Y residuals, refit controls; never rerank pseudo-Y."""
    q=fitted['q'];pseudo=fitted['yhat']+fitted['ey'][permutation]
    ey=pseudo-q@(q.T@pseudo)
    return correlation(fitted['ex'],ey)


def analyze(x,y,continuous,categories,slot,draws=2000,minimum=20,min_df=10,required_bootstrap=1800):
    x=np.asarray(x,float);y=np.asarray(y,float);c=np.asarray(continuous,float);cat=np.asarray(categories,str)
    assert x.ndim==y.ndim==1 and len(x)==len(y)==len(c)==len(cat)
    result=dict(n_cases=len(x),status='NOT_ESTIMABLE',P=None,CI95=None,bootstrap_valid=0,
                bootstrap_draws=draws,permutation_draws=draws)
    if len(x)<minimum:result['reason']='INSUFFICIENT_REVIEWED_CASES';return result
    levels=categorical_spec(cat)
    try:base=fit(x,y,c,cat,levels,min_df=min_df)
    except ValueError as error:result['reason']=str(error);return result
    result.update(status='ESTIMABLE',reason=None,rho=base['rho'],unadjusted_rho=base['unadjusted'],
                  design_rank=base['rank'],design_columns=base['columns'],design_condition=base['condition'],
                  residual_df=base['residual_df'],fixed_factor_indices=base['constant_factors'])
    groups={}
    for i,key in enumerate(map(tuple,cat)):groups.setdefault(key,[]).append(i)
    strata=[np.array(groups[key],int) for key in sorted(groups)]
    result['stratum_sizes']=sorted(map(len,strata));result['categorical_level_counts']=list(map(len,levels))
    seed=np.random.SeedSequence([2026102204,int(slot)])
    bootseed,permseed=seed.spawn(2);boot=np.random.default_rng(bootseed);rng=np.random.default_rng(permseed)
    effects=[];failed={}
    for unused in range(draws):
        chosen=np.concatenate([boot.choice(g,size=len(g),replace=True) for g in strata])
        try:current=fit(x[chosen],y[chosen],c[chosen],cat[chosen],levels,min_df=min_df)
        except ValueError as error:
            reason=str(error);failed[reason]=failed.get(reason,0)+1;continue
        effects.append(current['rho'])
    result.update(bootstrap_valid=len(effects),bootstrap_invalid_reasons=failed)
    if len(effects)>=required_bootstrap:result['CI95']=np.quantile(effects,[.025,.975]).tolist()
    else:result['interval_status']='INSUFFICIENT_VALID_BOOTSTRAP_REPLICATES'
    exceed=0
    for unused in range(draws):
        order=np.arange(len(x))
        for g in strata:order[g]=rng.permutation(g)
        value=permuted_statistic(base,order)
        if value is None:raise ValueError('DEGENERATE_PERMUTED_RESIDUAL_SCORE')
        exceed+=abs(value)>=abs(base['rho'])-1e-12
    result['P']=(1+exceed)/(draws+1);result['permutation_exceedances']=int(exceed)
    return result


def bh_fixed(results):
    """Unestimable opportunities stay in family at P=1, not negative results."""
    p=np.array([1. if r['P'] is None else r['P'] for r in results]);order=np.argsort(p,kind='stable')
    q=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    answer=np.empty(len(p));answer[order]=np.minimum(q,1)
    for r,value in zip(results,answer):r['BH_fixed_family']=float(value);r['family_size']=len(results)
