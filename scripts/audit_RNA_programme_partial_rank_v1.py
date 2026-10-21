"""Independent pandas-rank/least-squares replay of fixed-slot random inference.

Uses the declared random streams and recomputes every refit. Only scalar audit
discrepancies/replicate counts return; no case, permutation or score arrays
are persisted. This is a numerical replay, not independent biological data.
"""
import numpy as np
import pandas as pd


def reference_fit(x,y,continuous,categories,levels):
    if len(np.unique(x))<2 or len(np.unique(y))<2:raise ValueError('CONSTANT_PROGRAMME_SCORE')
    columns=[np.ones(len(x))]
    for j in range(continuous.shape[1]):
        rank=pd.Series(continuous[:,j]).rank(method='average').to_numpy();scale=rank.std()
        if scale<=1e-12:raise ValueError('CONSTANT_CONTINUOUS_CONTROL')
        columns.append((rank-rank.mean())/scale)
    for j,allowed in enumerate(levels):
        for name in allowed[1:]:columns.append((categories[:,j]==name).astype(float))
    design=np.column_stack(columns);rank=np.linalg.matrix_rank(design)
    if rank<design.shape[1]:raise ValueError('SINGULAR_FULL_DESIGN')
    if np.linalg.cond(design)>1e8:raise ValueError('ILL_CONDITIONED_FULL_DESIGN')
    if len(x)-rank<10:raise ValueError('INSUFFICIENT_RESIDUAL_DF')
    rx=pd.Series(x).rank(method='average').to_numpy();ry=pd.Series(y).rank(method='average').to_numpy()
    ex=rx-design@np.linalg.lstsq(design,rx,rcond=None)[0];fitted=design@np.linalg.lstsq(design,ry,rcond=None)[0];ey=ry-fitted
    if min(np.linalg.norm(ex),np.linalg.norm(ey))<=1e-10:raise ValueError('NO_RESIDUAL_SCORE_VARIATION')
    return float(np.corrcoef(ex,ey)[0,1]),design,ex,ey,fitted


def replay(x,y,continuous,categories,slot,produced,draws=2000,required_bootstrap=1800):
    levels=[sorted(set(categories[:,j])) for j in range(categories.shape[1])]
    if len(x)<20:
        assert produced['status']=='NOT_ESTIMABLE' and produced['reason']=='INSUFFICIENT_REVIEWED_CASES';return dict(status='NONESTIMABLE_GATE_REPLAYED')
    try:rho,design,ex,ey,fitted=reference_fit(x,y,continuous,categories,levels)
    except ValueError as error:
        assert produced['status']=='NOT_ESTIMABLE' and produced['reason']==str(error);return dict(status='NONESTIMABLE_GATE_REPLAYED')
    assert produced['status']=='ESTIMABLE';np.testing.assert_allclose(rho,produced['rho'],atol=1e-10,rtol=1e-10)
    groups={}
    for i,key in enumerate(map(tuple,categories)):groups.setdefault(key,[]).append(i)
    strata=[np.array(groups[k],int) for k in sorted(groups)]
    bootstrap_seed,permutation_seed=np.random.SeedSequence([2026102204,int(slot)]).spawn(2)
    boot=np.random.default_rng(bootstrap_seed);rng=np.random.default_rng(permutation_seed);effects=[];failures={}
    for unused in range(draws):
        chosen=np.concatenate([boot.choice(g,size=len(g),replace=True) for g in strata])
        try:value=reference_fit(x[chosen],y[chosen],continuous[chosen],categories[chosen],levels)[0]
        except ValueError as error:
            reason=str(error);failures[reason]=failures.get(reason,0)+1;continue
        effects.append(value)
    assert len(effects)==produced['bootstrap_valid'] and failures==produced['bootstrap_invalid_reasons']
    if len(effects)>=required_bootstrap:np.testing.assert_allclose(np.quantile(effects,[.025,.975]),produced['CI95'],atol=1e-9,rtol=1e-9)
    else:assert produced['CI95'] is None
    exceed=0
    for unused in range(draws):
        order=np.arange(len(x))
        for group in strata:order[group]=rng.permutation(group)
        pseudo=fitted+ey[order];residual=pseudo-design@np.linalg.lstsq(design,pseudo,rcond=None)[0]
        value=float(np.corrcoef(ex,residual)[0,1]);assert np.isfinite(value)
        exceed+=abs(value)>=abs(rho)-1e-12
    assert exceed==produced['permutation_exceedances'];np.testing.assert_allclose((1+exceed)/(draws+1),produced['P'],atol=1e-15)
    return dict(status='ALL_BOOTSTRAP_AND_PERMUTATION_REPLAYED',bootstrap_draws=draws,permutation_draws=draws,
                bootstrap_valid=len(effects),observed_sign=int(np.sign(rho)))
