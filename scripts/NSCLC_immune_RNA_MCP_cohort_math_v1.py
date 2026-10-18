"""Full nuisance-span HC3 models for the fixed six programme-pair bundles."""
import hashlib
import numpy as np
from scipy.stats import rankdata

from HNSCC_author_state_clinical_kernel_v2 import fit


def seed(slot,component,stream):
    return int.from_bytes(hashlib.sha256(('%s|%s|%s|2026101901'%(slot,component,stream)).encode()).digest()[:8],'little')


def design(x,y,C,names,kind,expected_rank=None,minimum=50):
    x=np.asarray(x,float);y=np.asarray(y,float);C=np.asarray(C,float)
    n=len(x)
    if n<minimum:return None,'FEWER_THAN50_AUTHOR_CASES'
    if not np.isfinite(np.c_[x,y,C]).all():return None,'NONFINITE_INPUT'
    if np.ptp(x)==0 or np.ptp(y)==0:return None,'CONSTANT_EXPOSURE_OR_READOUT'
    controls=C.copy()
    if kind=='coordination':
        take=np.array([not s.startswith('agent::') for s in names])
        controls[:,take]=rankdata(controls[:,take],axis=0)
        x=rankdata(x);x=(x-x.mean())/x.std(ddof=0)
        y=rankdata(y);y=(y-y.mean())/y.std(ddof=0)
    elif kind!='clinical':raise ValueError('UNDECLARED_MODEL_KIND')
    Zraw=np.c_[np.ones(n),controls]
    mean=controls.mean(axis=0);sd=controls.std(axis=0,ddof=0)
    Z=np.c_[np.ones(n),(controls-mean)/np.where(sd>0,sd,1)]
    u,s,_=np.linalg.svd(Z,full_matrices=False)
    rank=int((s>s[0]*1e-12).sum());machine=int((s>s[0]*np.finfo(float).eps*max(Z.shape)).sum())
    if rank!=machine or s[rank-1]/s[0]<1e-8:return None,'UNSTABLE_COMPLETE_NUISANCE_SPAN'
    if expected_rank is not None and rank!=expected_rank:return None,'BOOTSTRAP_NUISANCE_RANK_CHANGED'
    Q=u[:,:rank]
    error=float(np.max(np.abs(Zraw-Q@(Q.T@Zraw))))
    if error>1e-8*max(1.,float(np.max(np.abs(Zraw)))):return None,'NUISANCE_SPAN_RECONSTRUCTION_FAILED'
    rx=x-Q@(Q.T@x)
    if np.linalg.norm(rx)<=1e-8*np.sqrt(n):return None,'TARGET_NOT_IDENTIFIABLE'
    D=np.c_[Q,x];singular=np.linalg.svd(D,compute_uv=False)
    if n-D.shape[1]<10:return None,'FEWER_THAN10_RESIDUAL_DF'
    if singular[0]/singular[-1]>=1e8:return None,'ILL_CONDITIONED_TARGET_DESIGN'
    inverse=np.linalg.pinv(D,rcond=1e-12);h=np.einsum('ij,ji->i',D,inverse)
    if h.max()>=.999999:return None,'EXTREME_LEVERAGE'
    return dict(x=x,y=y,controls=C,names=np.asarray(names),kind=kind,D=D,Q=Q,Zraw=Zraw,nuisance_rank=rank,
        condition=float(singular[0]/singular[-1]),df=n-D.shape[1],max_leverage=float(h.max()),span_error=error),'ESTIMABLE'


def estimate(model):
    g=fit(model['y'],model['D'])
    if model['kind']=='coordination':
        rx=model['x']-model['Q']@(model['Q'].T@model['x']);ry=model['y']-model['Q']@(model['Q'].T@model['y'])
        effect=float(rx@ry/(np.linalg.norm(rx)*np.linalg.norm(ry)))
    else:effect=g['beta']
    return g,effect


def wild_values(g,draws,random_seed):
    Z=g['d'][:,:-1];pinv=np.linalg.pinv(Z,rcond=1e-12);mu=Z@(pinv@g['y']);h0=np.einsum('ij,ji->i',Z,pinv)
    e=(g['y']-mu)/(1-h0);row=g['inverse'][-1];v=row/(1-g['h']);rng=np.random.default_rng(random_seed);values=np.empty(draws)
    for start in range(0,draws,1000):
        signs=np.where(rng.random((min(1000,draws-start),len(mu)))<.5,-1.,1.).T
        noise=e[:,None]*signs;residual=noise-g['d']@(g['inverse']@noise)
        values[start:start+noise.shape[1]]=(row@noise)/np.sqrt(np.sum((v[:,None]*residual)**2,axis=0))
    if not np.isfinite(values).all():raise ValueError('NONFINITE_WILD_STATISTIC')
    return values


def takes(strata,draws,random_seed):
    rng=np.random.default_rng(random_seed);out=np.empty((draws,len(strata)),np.int16);offset=0
    for label in sorted(set(strata)):
        ix=np.flatnonzero(strata==label)
        out[:,offset:offset+len(ix)]=ix[rng.integers(0,len(ix),size=(draws,len(ix)))];offset+=len(ix)
    assert offset==len(strata)
    return out


def bootstrap(x,y,C,names,kind,indices,expected_rank):
    values=np.full(len(indices),np.nan);reasons=[]
    for i,ix in enumerate(indices):
        m,status=design(x[ix],y[ix],C[ix],names,kind,expected_rank=expected_rank)
        if m is not None:
            try:_,values[i]=estimate(m)
            except ValueError as e:status=str(e)
        reasons.append(status)
    return values,np.asarray(reasons)


def leave_one(x,y,C,names,kind,expected_rank):
    values=np.full(len(x),np.nan)
    for i in range(len(x)):
        keep=np.arange(len(x))!=i;m,_=design(x[keep],y[keep],C[keep],names,kind,expected_rank=expected_rank)
        if m is not None:
            try:_,values[i]=estimate(m)
            except ValueError:pass
    return values


def adjust(p,by=False):
    p=np.asarray(p,float);m=len(p);order=np.argsort(p,kind='stable')
    values=p[order]*m/np.arange(1,m+1)
    if by:values*=np.sum(1/np.arange(1,m+1))
    out=np.empty(m);out[order]=np.minimum(1,np.minimum.accumulate(values[::-1])[::-1]);return out
