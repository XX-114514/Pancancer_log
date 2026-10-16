"""Fixed-design single-target HC3 inference and stratified patient resampling."""
import hashlib
import numpy as np


def seed(model_id,stream):return int.from_bytes(hashlib.sha256((str(model_id)+'|'+stream+'|20261016').encode()).digest()[:8],'little')


def fit(y,d,min_df=10,hc3=True):
    y=np.asarray(y,float);d=np.asarray(d,float)
    if not np.isfinite(y).all() or not np.isfinite(d).all():raise ValueError('NONFINITE')
    s=np.linalg.svd(d,compute_uv=False);rank=int(np.sum(s>s[0]*1e-12));machine=int(np.sum(s>s[0]*np.finfo(float).eps*max(d.shape)))
    if rank!=machine or rank!=d.shape[1] or len(y)-rank<min_df or s[0]/s[-1]>=1e8:raise ValueError('RANK_DF_OR_CONDITION')
    inverse=np.linalg.pinv(d,rcond=1e-12);b=inverse@y;error=y-d@b;h=np.einsum('ij,ji->i',d,inverse);variance=np.nan;t=np.nan
    if hc3:
        if h.max()>=.999999:raise ValueError('EXTREME_LEVERAGE')
        variance=float(np.square(inverse[-1]*error/(1-h)).sum())
        if not np.isfinite(variance) or variance<=0:raise ValueError('DEGENERATE_HC3_VARIANCE')
        t=float(b[-1]/np.sqrt(variance))
    return dict(y=y,d=d,inverse=inverse,beta=float(b[-1]),variance=variance,t=t,h=h,condition=float(s[0]/s[-1]),df=len(y)-rank)


def wild_count(g,draws,random_seed):
    z=g['d'][:,:-1];pinv=np.linalg.pinv(z,rcond=1e-12);mu=z@(pinv@g['y']);h0=np.einsum('ij,ji->i',z,pinv);e0=(g['y']-mu)/np.sqrt(1-h0);row=g['inverse'][-1];residual_operator=np.eye(len(mu))-g['d']@g['inverse'];v=row/(1-g['h']);rng=np.random.default_rng(random_seed);count=0
    for start in range(0,draws,1000):
        signs=np.where(rng.random((min(1000,draws-start),len(mu)))<.5,-1.,1.).T;noise=e0[:,None]*signs;errors=residual_operator@noise;statistic=(row@noise)/np.sqrt(np.sum(np.square(v[:,None]*errors),axis=0));assert np.isfinite(statistic).all();count+=int((np.abs(statistic)>=abs(g['t'])-1e-12).sum())
    return count


def bootstrap(y,d,strata,draws,random_seed):
    rng=np.random.default_rng(random_seed);take=np.empty((draws,len(y)),int);offset=0
    for level in sorted(set(strata)):
        ids=np.flatnonzero(strata==level);take[:,offset:offset+len(ids)]=ids[rng.integers(0,len(ids),size=(draws,len(ids)))];offset+=len(ids)
    assert offset==len(y);out=np.full(draws,np.nan)
    for start in range(0,draws,128):
        end=min(draws,start+128);indices=take[start:end];u,s,vt=np.linalg.svd(d[indices],full_matrices=False);rank=(s>s[:,:1]*1e-12).sum(axis=1);machine=(s>s[:,:1]*np.finfo(float).eps*max(d.shape)).sum(axis=1);good=(rank==d.shape[1])&(rank==machine)&(s[:,0]/np.maximum(s[:,-1],np.finfo(float).tiny)<1e8)
        projection=np.einsum('bni,bn->bi',u[good],y[indices[good]]);coefficient=np.sum(vt[good,:,-1]*projection/s[good],axis=1);out[np.arange(start,end)[good]]=coefficient
    return out


def loo(y,d):
    out=np.full(len(y),np.nan)
    for j in range(len(y)):
        keep=np.arange(len(y))!=j
        try:out[j]=fit(y[keep],d[keep],min_df=1,hc3=False)['beta']
        except ValueError:pass
    return out
