"""Fixed-scale full-rank linear conditional association, binary clinical outcome."""
import numpy as np


def fit(y, design, min_df=1, check_hc3=True):
    y=np.asarray(y,float);d=np.asarray(design,float)
    if d.ndim!=2 or len(d)!=len(y) or not np.isfinite(d).all() or not np.isfinite(y).all():
        raise ValueError('NONFINITE_INPUT')
    s=np.linalg.svd(d,compute_uv=False)
    rank=int(np.sum(s>s[0]*1e-12))
    machine=int(np.sum(s>s[0]*np.finfo(float).eps*max(d.shape)))
    if rank!=machine or rank!=d.shape[1]:raise ValueError('NONFULL_OR_AMBIGUOUS_RANK')
    if s[0]/s[-1]>=1e8 or len(y)-rank<min_df:raise ValueError('LOW_DF_OR_CONDITION')
    inverse=np.linalg.pinv(d,rcond=1e-12);beta=inverse@y;e=y-d@beta
    h=np.sum(d*inverse.T,axis=1)
    if check_hc3 and h.max()>=.999999:raise ValueError('EXTREME_LEVERAGE')
    if check_hc3:
        variance=np.sum((inverse*e[None,:]/(1-h[None,:]))**2,axis=1)
        if not np.isfinite(variance[-2:]).all() or np.any(variance[-2:]<=0):raise ValueError('DEGENERATE_VARIANCE')
    else:variance=np.full(len(beta),np.nan)
    return dict(y=y,d=d,inverse=inverse,beta=beta,residual=e,h=h,variance=variance,t=beta/np.sqrt(variance),condition=s[0]/s[-1],rank=rank)


def wild_statistics(g,target,signs):
    d=g['d'];z=np.delete(d,target,axis=1);p=np.linalg.pinv(z,rcond=1e-12)
    mu=z@(p@g['y']);h0=np.sum(z*p.T,axis=1)
    yy=mu[:,None]+((g['y']-mu)/np.sqrt(1-h0))[:,None]*signs
    b=g['inverse']@yy;e=yy-d@b
    variance=np.sum((g['inverse'][target,:,None]*e/(1-g['h'][:,None]))**2,axis=0)
    t=b[target]/np.sqrt(variance)
    assert np.isfinite(t).all()
    return t


def wild_count(g,target,draws,seed):
    rng=np.random.default_rng(seed);exceed=0
    # Cache the exact restricted fit; do not pseudoinvert for every batch.
    z=np.delete(g['d'],target,axis=1);p=np.linalg.pinv(z,rcond=1e-12)
    mu=z@(p@g['y']);e0=(g['y']-mu)/np.sqrt(1-np.sum(z*p.T,axis=1))
    v=g['inverse'][target]/(1-g['h'])
    for start in range(0,draws,1000):
        signs=np.where(rng.random((min(1000,draws-start),len(mu)))<.5,-1.,1.).T
        yy=mu[:,None]+e0[:,None]*signs;b=g['inverse']@yy;error=yy-g['d']@b
        t=b[target]/np.sqrt(np.sum((v[:,None]*error)**2,axis=0))
        assert np.isfinite(t).all()
        exceed+=int((np.abs(t)>=abs(g['t'][target])-1e-12).sum())
    return exceed


def resample(y,d,draws,seed):
    rng=np.random.default_rng(seed);out=np.full((draws,2),np.nan)
    for i in range(draws):
        take=rng.integers(0,len(y),len(y))
        try:out[i]=fit(y[take],d[take],check_hc3=False)['beta'][-2:]
        except ValueError:pass
    return out


def loo(y,d):
    out=np.full((len(y),2),np.nan)
    for i in range(len(y)):
        keep=np.arange(len(y))!=i
        try:out[i]=fit(y[keep],d[keep],check_hc3=False)['beta'][-2:]
        except ValueError:pass
    return out
