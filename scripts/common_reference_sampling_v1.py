"""Balanced reference exposure with explicit local-stratum minimums."""
import numpy as np
from annotation_pilot_sampling_v1 import balanced_allocations


def reserve_then_balance(sizes,reserves,total,seed):
    sizes,reserves=np.asarray(sizes,np.int64),np.asarray(reserves,np.int64)
    if sizes.ndim!=1 or reserves.shape!=sizes.shape or np.any(sizes<0) or np.any(reserves<0) or np.any(reserves>sizes) or total<0:raise ValueError('Invalid reference capacities or reserves')
    budget=min(int(total),int(sizes.sum()))
    if reserves.sum()>budget:raise ValueError('Reference budget cannot cover fixed local-stratum reserves')
    return reserves+balanced_allocations(sizes-reserves,budget-int(reserves.sum()),np.random.default_rng(seed))


def stratum_reserve(strata,eligible,minimum=2):
    strata,eligible=np.asarray(strata),np.asarray(eligible)
    if strata.ndim!=1 or eligible.shape!=strata.shape or eligible.dtype!=bool or minimum<0:raise ValueError('Invalid stratum reserve input')
    _,counts=np.unique(strata[eligible].astype(str),return_counts=True)
    return int(np.minimum(counts,minimum).sum())


def select_reference(unit,strata,eligible,cap,seed,minimum=2):
    unit,strata,eligible=np.asarray(unit),np.asarray(strata),np.asarray(eligible)
    if unit.ndim!=1 or strata.shape!=unit.shape or eligible.shape!=unit.shape or eligible.dtype!=bool or cap<0 or minimum<0:raise ValueError('Invalid reference row axes or quota')
    positions=np.flatnonzero(eligible);budget=min(int(cap),len(positions));rng=np.random.default_rng(seed)
    if not len(positions):return np.array([],np.int64)
    if any(str(x).strip().lower() in ['', 'nan','none','unknown','__missing__','__source_value_missing__'] for x in unit[eligible]):raise ValueError('Eligible reference row lacks expression unit')
    _,code,count=np.unique(strata[positions].astype(str),return_inverse=True,return_counts=True)
    reserve=np.minimum(count,minimum)
    if reserve.sum()>budget:raise ValueError('Reference quota cannot preserve fixed local strata')
    selected=[]
    for s,amount in enumerate(reserve):
        if not amount:continue
        rows=positions[code==s];names,inverse,counts=np.unique(unit[rows].astype(str),return_inverse=True,return_counts=True)
        quotas=balanced_allocations(counts,int(amount),rng)
        for u,k in enumerate(quotas):
            if k:selected.extend(rng.choice(rows[inverse==u],int(k),replace=False).tolist())
    chosen=np.zeros(len(unit),bool);chosen[np.asarray(selected,np.int64)]=True;remaining=positions[~chosen[positions]];left=budget-len(selected)
    if left:
        names,inverse,counts=np.unique(unit[remaining].astype(str),return_inverse=True,return_counts=True);quotas=balanced_allocations(counts,left,rng)
        for u,k in enumerate(quotas):
            if k:selected.extend(rng.choice(remaining[inverse==u],int(k),replace=False).tolist())
    result=np.sort(np.asarray(selected,np.int64))
    if len(result)!=budget or len(np.unique(result))!=budget or not eligible[result].all():raise ValueError('Reference sampler changed budget, eligibility or uniqueness')
    labels,selected_count=np.unique(strata[result].astype(str),return_counts=True)
    expected_labels,expected_count=np.unique(strata[positions].astype(str),return_counts=True)
    observed=dict(zip(labels,selected_count))
    if any(observed.get(label,0)<min(n,minimum) for label,n in zip(expected_labels,expected_count)):raise ValueError('Reference rare-stratum reserve lost')
    return result
