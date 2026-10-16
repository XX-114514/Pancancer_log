"""Frozen block-rank statistics and author-case resampling, no fitted expression."""
import hashlib
import numpy as np
from scipy.stats import rankdata


def seed(identifier,stream):return int.from_bytes(hashlib.sha256((identifier+'|'+stream+'|NSCLC_MFA_v1').encode()).digest()[:8],'little')


def centre_rank(x,blocks):
    ranks=rankdata(x,axis=-1,method='average')
    for level in sorted(set(blocks)):
        ids=np.flatnonzero(blocks==level);ranks[...,ids]-=ranks[...,ids].mean(axis=-1,keepdims=True)
    return ranks


def point(a,b,blocks,min_active=30):
    active=sum(np.sum(blocks==s) for s in set(blocks) if np.sum(blocks==s)>=2)
    if active<min_active:raise ValueError('INSUFFICIENT_CASES_IN_NON_SINGLETON_BLOCKS')
    u=centre_rank(a,blocks);v=centre_rank(b,blocks);denom=np.sqrt(np.sum(u*u)*np.sum(v*v))
    if not np.isfinite(denom) or denom<1e-12:raise ValueError('CONSTANT_OR_NONFINITE_BLOCK_RANK_SCORE')
    return float(np.sum(u*v)/denom),u,v,float(denom),active


def permutation_count(u,v,denom,observed,blocks,draws,random_seed):
    rng=np.random.default_rng(random_seed);count=0
    for start in range(0,draws,256):
        m=min(256,draws-start);take=np.broadcast_to(np.arange(len(u)),(m,len(u))).copy()
        for level in sorted(set(blocks)):
            ids=np.flatnonzero(blocks==level);take[:,ids]=ids[np.argsort(rng.random((m,len(ids))),axis=1)]
        values=v[take]@u/denom;count+=int((abs(values)>=abs(observed)-1e-12).sum())
    return count


def bootstrap(a,b,blocks,draws,random_seed):
    rng=np.random.default_rng(random_seed);take=np.empty((draws,len(a)),int);labels=[];offset=0
    for level in sorted(set(blocks)):
        ids=np.flatnonzero(blocks==level);take[:,offset:offset+len(ids)]=ids[rng.integers(0,len(ids),(draws,len(ids)))];offset+=len(ids);labels.extend([level]*len(ids))
    levels=np.array(labels);u=centre_rank(a[take],levels);v=centre_rank(b[take],levels);denom=np.sqrt(np.sum(u*u,axis=1)*np.sum(v*v,axis=1));out=np.full(draws,np.nan);good=np.isfinite(denom)&(denom>1e-12);out[good]=np.sum(u[good]*v[good],axis=1)/denom[good];return out
