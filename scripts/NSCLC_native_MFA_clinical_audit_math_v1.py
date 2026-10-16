"""Independent QR clinical inference and resampling, separate from producer SVD."""
import numpy as np


def fit(y,d):
    q,r=np.linalg.qr(d,mode='reduced');coef=np.linalg.solve(r,q.T@y)
    unit=np.zeros(d.shape[1]);unit[-1]=1
    row=np.linalg.solve(r.T,unit)@q.T
    residual=y-q@(q.T@y);h=(q*q).sum(axis=1)
    variance=np.sum((row*residual/(1-h))**2)
    return dict(beta=float(coef[-1]),t=float(coef[-1]/np.sqrt(variance)),se=float(np.sqrt(variance)),q=q,row=row,h=h,y=y,d=d)


def wild(g,draws,random_seed):
    z=g['d'][:,:-1];qn,_=np.linalg.qr(z,mode='reduced');e=(g['y']-qn@(qn.T@g['y']))/(1-(qn*qn).sum(axis=1))
    rng=np.random.default_rng(random_seed);count=0
    for start in range(0,draws,1000):
        noise=e[:,None]*np.where(rng.random((min(1000,draws-start),len(e)))<.5,-1.,1.).T
        residual=noise-g['q']@(g['q'].T@noise)
        value=(g['row']@noise)/np.sqrt(np.sum((g['row'][:,None]*residual/(1-g['h'][:,None]))**2,axis=0))
        count+=int((abs(value)>=abs(g['t'])-1e-12).sum())
    return count


def bootstrap(y,d,strata,draws,random_seed):
    rng=np.random.default_rng(random_seed);take=np.empty((draws,len(y)),int);offset=0
    for level in sorted(set(strata)):
        ids=np.flatnonzero(strata==level);take[:,offset:offset+len(ids)]=ids[rng.integers(0,len(ids),(draws,len(ids)))];offset+=len(ids)
    result=np.full(draws,np.nan)
    for start in range(0,draws,128):
        ix=take[start:start+128];ds=d[ix];s=np.linalg.svd(ds,compute_uv=False)
        rank=(s>s[:,[0]]*1e-12).sum(axis=1);machine=(s>s[:,[0]]*np.finfo(float).eps*max(d.shape)).sum(axis=1)
        good=(rank==d.shape[1])&(rank==machine)&(s[:,0]/np.maximum(s[:,-1],np.finfo(float).tiny)<1e8)
        q,r=np.linalg.qr(ds[good],mode='reduced');rhs=np.einsum('bni,bn->bi',q,y[ix[good]])
        b=np.linalg.solve(r,rhs[:,:,None])[:,-1,0]
        result[np.arange(start,start+len(ix))[good]]=b
    return result


def loo(y,d):
    result=np.full(len(y),np.nan)
    for j in range(len(y)):
        take=np.arange(len(y))!=j;ds=d[take];s=np.linalg.svd(ds,compute_uv=False)
        rank=int((s>s[0]*1e-12).sum());machine=int((s>s[0]*np.finfo(float).eps*max(ds.shape)).sum())
        if rank!=d.shape[1] or rank!=machine or s[0]/s[-1]>=1e8:continue
        q,r=np.linalg.qr(ds,mode='reduced');result[j]=np.linalg.solve(r,q.T@y[take])[-1]
    return result


def adjustment(p,by=False):
    p=np.asarray(p);order=np.argsort(p,kind='stable');m=len(p);values=p[order]*m/np.arange(1,m+1)
    if by:values*=np.sum(1/np.arange(1,m+1))
    values=np.minimum.accumulate(values[::-1])[::-1];result=np.empty(m);result[order]=np.minimum(values,1.)
    return result
