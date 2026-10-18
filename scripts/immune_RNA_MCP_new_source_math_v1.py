"""Prespecified small-cohort wrappers; old finite families remain unchanged."""
import hashlib
import numpy as np
from scipy.stats import rankdata
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design as _design, estimate, wild_values, adjust
from NSCLC_immune_RNA_MCP_cohort_audit_math_v1 import prepare as _prepare, evaluate, wild as independent_wild


def seed(slot,stream):
    return int.from_bytes(hashlib.sha256((str(slot)+'|'+stream+'|2026101904').encode()).digest()[:8],'little')


def design(x,y,C,names,kind,expected_rank=None):
    if len(x)<20:return None,'FEWER_THAN20_CASE_UNITS'
    m,status=_design(x,y,C,names,kind,expected_rank=expected_rank,minimum=20)
    if m is not None and np.linalg.norm(m['y']-m['Q']@(m['Q'].T@m['y']))<=1e-8*max(1.,np.linalg.norm(m['y'])):return None,'READOUT_EXPLAINED_BY_CONTROL_SPAN'
    return m,status


def independent_design(x,y,C,names,kind,expected_rank=None):
    if len(x)<20:return None,'FEWER_THAN20_CASE_UNITS'
    m,status=_prepare(x,y,C,names,kind,expected_rank=expected_rank,minimum=20)
    if m is not None and np.linalg.norm(m['y']-m['Q']@(m['Q'].T@m['y']))<=1e-8*max(1.,np.linalg.norm(m['y'])):return None,'READOUT_EXPLAINED_BY_CONTROL_SPAN'
    return m,status


def bootstrap_indices(strata,slot,n=2000):
    rng=np.random.default_rng(seed(slot,'case_bootstrap'));pieces=[]
    for label in sorted(set(strata)):
        ids=np.flatnonzero(strata==label);pieces.append(ids[rng.integers(0,len(ids),size=(n,len(ids)))])
    return np.column_stack(pieces).astype(np.int16)


def bootstrap(x,y,C,names,kind,indices,rank,independent=False):
    values=np.full(len(indices),np.nan);reasons=[]
    for k,ix in enumerate(indices):
        m,status=(independent_design if independent else design)(x[ix],y[ix],C[ix],names,kind,expected_rank=rank)
        if m is not None:
            try:values[k]=evaluate(m)['effect'] if independent else estimate(m)[1]
            except ValueError as e:status=str(e)
        reasons.append(status)
    return values,np.asarray(reasons)


def leave_one(x,y,C,names,kind,rank,independent=False):
    result=np.full(len(x),np.nan);reasons=[]
    for k in range(len(x)):
        take=np.arange(len(x))!=k;m,status=(independent_design if independent else design)(x[take],y[take],C[take],names,kind,expected_rank=rank)
        if m is not None:
            try:result[k]=evaluate(m)['effect'] if independent else estimate(m)[1]
            except ValueError as e:status=str(e)
        reasons.append(status)
    return result,np.asarray(reasons)


def rank_effect(x,y,binary=False):
    if min(np.ptp(x),np.ptp(y))==0:return np.nan
    if binary:
        n1=int((x==1).sum());n0=int((x==0).sum())
        return float(2*(rankdata(y)[x==1].sum()-n1*(n1+1)/2)/(n1*n0)-1) if n1 and n0 else np.nan
    a=rankdata(x);b=rankdata(y);a-=a.mean();b-=b.mean()
    return float(a@b/np.linalg.norm(a)/np.linalg.norm(b))


def rank_permutations(x,y,slot,binary=False,draws=19999):
    b=rankdata(y);rng=np.random.default_rng(seed(slot,'rank_permutation'));values=np.empty(draws)
    if binary:
        n1=int((x==1).sum());n0=int((x==0).sum())
        for lo in range(0,draws,1000):
            idx=np.argsort(rng.random((min(1000,draws-lo),len(x))),axis=1);s=b[idx][:,x==1].sum(axis=1);values[lo:lo+len(idx)]=2*(s-n1*(n1+1)/2)/(n1*n0)-1
    else:
        a=rankdata(x);a-=a.mean();b-=b.mean();den=np.linalg.norm(a)*np.linalg.norm(b)
        for lo in range(0,draws,1000):
            idx=np.argsort(rng.random((min(1000,draws-lo),len(x))),axis=1);values[lo:lo+len(idx)]=(b[idx]@a)/den
    return values
