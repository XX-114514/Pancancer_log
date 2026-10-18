"""Prespecified paired median/sign and full-span partial-rank cochange math."""
import hashlib
from math import comb

import numpy as np
from scipy.stats import rankdata

from HNSCC_author_state_clinical_kernel_v2 import fit
from NSCLC_immune_RNA_MCP_cohort_math_v1 import adjust, wild_values


def seed(slot,stream,base=2026101903):
    return int.from_bytes(hashlib.sha256(('%s|%s|%s'%(base,slot,stream)).encode()).digest()[:8],'little')


def median_shift(delta):
    d=np.asarray(delta,float)
    if not np.isfinite(d).all():return None,'NONFINITE_PAIRED_DELTA'
    if len(d)<20:return None,'FEWER_THAN20_COMPLETE_PAIRED_PATIENTS'
    used=d[np.abs(d)>1e-10];n=len(used);pos=int((used>0).sum());neg=n-pos
    if n<10:return None,'FEWER_THAN10_NONZERO_PAIRED_DELTAS'
    p=min(1.,2*sum(comb(n,j) for j in range(min(pos,neg)+1))/2**n)
    ordered=np.sort(d);N=len(d);ks=[k for k in range(1,N//2+1) if 2*sum(comb(N,j) for j in range(k))/2**N<=.05]
    k=max(ks);ranks=rankdata(np.abs(used));rb=float(np.sum(ranks*np.sign(used))/ranks.sum())
    return dict(effect=float(np.median(d)),CI_low=float(ordered[k-1]),CI_high=float(ordered[N-k]),
        p=float(p),patients=N,positive=pos,negative=neg,numerical_zero=N-n,nonzero=n,
        exact_interval_minimum_coverage=float(1-2*sum(comb(N,j) for j in range(k))/2**N),
        mean_delta=float(d.mean()),paired_rank_biserial=rb),'ESTIMABLE'


def design(x,y,C,names,expected_rank=None):
    x=np.asarray(x,float);y=np.asarray(y,float);C=np.asarray(C,float);n=len(x)
    if n<20:return None,'FEWER_THAN20_COMPLETE_PAIRED_PATIENTS'
    if not np.isfinite(np.c_[x,y,C]).all():return None,'NONFINITE_INPUT'
    if min(np.ptp(x),np.ptp(y))==0:return None,'CONSTANT_PAIRED_TARGET_OR_READOUT'
    controls=C.copy();continuous=np.array([s not in ['treatment_DplusT','HPV_positive'] for s in names])
    controls[:,continuous]=rankdata(controls[:,continuous],axis=0)
    x=rankdata(x);x=(x-x.mean())/x.std(ddof=0)
    y=rankdata(y);y=(y-y.mean())/y.std(ddof=0)
    mean=controls.mean(axis=0);sd=controls.std(axis=0,ddof=0)
    Z=np.c_[np.ones(n),(controls-mean)/np.where(sd>0,sd,1)]
    Zraw=np.c_[np.ones(n),controls];u,s,_=np.linalg.svd(Z,full_matrices=False)
    rank=int((s>s[0]*1e-12).sum());machine=int((s>s[0]*np.finfo(float).eps*max(Z.shape)).sum())
    if rank!=machine or s[rank-1]/s[0]<1e-8:return None,'UNSTABLE_COMPLETE_NUISANCE_SPAN'
    if expected_rank is not None and rank!=expected_rank:return None,'BOOTSTRAP_NUISANCE_RANK_CHANGED'
    Q=u[:,:rank];error=float(np.max(np.abs(Zraw-Q@(Q.T@Zraw))))
    if error>1e-8*max(1,float(np.max(np.abs(Zraw)))):return None,'NUISANCE_SPAN_RECONSTRUCTION_FAILED'
    rx=x-Q@(Q.T@x);ry=y-Q@(Q.T@y)
    if np.linalg.norm(rx)<=1e-8*max(1,np.linalg.norm(x)):return None,'TARGET_NOT_IDENTIFIABLE'
    if np.linalg.norm(ry)<=1e-8*max(1,np.linalg.norm(y)):return None,'READOUT_NOT_IDENTIFIABLE'
    D=np.c_[Q,x];sing=np.linalg.svd(D,compute_uv=False)
    if n-D.shape[1]<10:return None,'FEWER_THAN10_RESIDUAL_DF'
    if sing[0]/sing[-1]>=1e8:return None,'ILL_CONDITIONED_TARGET_DESIGN'
    inv=np.linalg.pinv(D,rcond=1e-12);h=np.einsum('ij,ji->i',D,inv)
    if h.max()>=.999999:return None,'EXTREME_LEVERAGE'
    return dict(x=x,y=y,C=C,names=np.asarray(names),Q=Q,D=D,nuisance_rank=rank,
        condition=float(sing[0]/sing[-1]),
        max_leverage=float(h.max()),residual_df=n-D.shape[1],span_error=error),'ESTIMABLE'


def estimate(model):
    return fit(model['y'],model['D'])


def effect(model):
    Q=model['Q'];x=model['x'];y=model['y'];rx=x-Q@(Q.T@x);ry=y-Q@(Q.T@y)
    return float(rx@ry/np.linalg.norm(rx)/np.linalg.norm(ry))


def patient_takes(n,slot):
    return np.random.default_rng(seed(slot,'patient_bootstrap')).integers(0,n,size=(2000,n),dtype=np.int16)


def bootstrap(x,y,C,names,indices,rank):
    values=np.full(len(indices),np.nan);reasons=[]
    for i,ix in enumerate(indices):
        m,status=design(x[ix],y[ix],C[ix],names,expected_rank=rank)
        if m is not None:values[i]=effect(m)
        reasons.append(status)
    return values,np.asarray(reasons)


def leave_one(x,y,C,names,rank):
    values=np.full(len(x),np.nan);reasons=[]
    for i in range(len(x)):
        keep=np.arange(len(x))!=i;m,status=design(x[keep],y[keep],C[keep],names,expected_rank=rank)
        if m is not None:values[i]=effect(m)
        reasons.append(status)
    return values,np.asarray(reasons)
