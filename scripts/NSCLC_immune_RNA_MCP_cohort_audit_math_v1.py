"""Independent pandas ranks/pivoted QR for all fixed cohort model quantities."""
import numpy as np
import pandas as pd
from scipy.linalg import qr


def prepare(x,y,C,names,kind,expected_rank=None,minimum=50):
    x=np.asarray(x,float);y=np.asarray(y,float);C=np.asarray(C,float);n=len(x)
    if n<minimum:return None,'FEWER_THAN50_AUTHOR_CASES'
    if not np.isfinite(np.c_[x,y,C]).all():return None,'NONFINITE_INPUT'
    if np.ptp(x)==0 or np.ptp(y)==0:return None,'CONSTANT_EXPOSURE_OR_READOUT'
    c=C.copy()
    if kind=='coordination':
        for j,name in enumerate(names):
            if not name.startswith('agent::'):c[:,j]=pd.Series(c[:,j]).rank().to_numpy()
        x=pd.Series(x).rank().to_numpy();x=(x-x.mean())/x.std(ddof=0)
        y=pd.Series(y).rank().to_numpy();y=(y-y.mean())/y.std(ddof=0)
    elif kind!='clinical':raise ValueError('UNDECLARED_MODEL_KIND')
    raw=np.c_[np.ones(n),c];sd=c.std(axis=0,ddof=0);scaled=np.c_[np.ones(n),(c-c.mean(axis=0))/np.where(sd>0,sd,1)]
    singular=np.linalg.svd(scaled,compute_uv=False);rank=int((singular>singular[0]*1e-12).sum())
    machine=int((singular>singular[0]*np.finfo(float).eps*max(scaled.shape)).sum())
    if rank!=machine or singular[rank-1]/singular[0]<1e-8:return None,'UNSTABLE_COMPLETE_NUISANCE_SPAN'
    if expected_rank is not None and rank!=expected_rank:return None,'BOOTSTRAP_NUISANCE_RANK_CHANGED'
    q,r,_=qr(scaled,mode='economic',pivoting=True);Q=q[:,:rank]
    assert int((np.abs(np.diag(r))>singular[0]*1e-12).sum())==rank
    error=float(np.max(np.abs(raw-Q@(Q.T@raw))))
    if error>1e-8*max(1,float(np.max(np.abs(raw)))):return None,'NUISANCE_SPAN_RECONSTRUCTION_FAILED'
    if np.linalg.norm(x-Q@(Q.T@x))<=1e-8*np.sqrt(n):return None,'TARGET_NOT_IDENTIFIABLE'
    D=np.c_[Q,x];s=np.linalg.svd(D,compute_uv=False)
    if n-D.shape[1]<10:return None,'FEWER_THAN10_RESIDUAL_DF'
    if s[0]/s[-1]>=1e8:return None,'ILL_CONDITIONED_TARGET_DESIGN'
    fullq,fullr=np.linalg.qr(D,mode='reduced');h=np.sum(fullq**2,axis=1)
    if h.max()>=.999999:return None,'EXTREME_LEVERAGE'
    return dict(x=x,y=y,C=C,names=names,kind=kind,D=D,Q=Q,Zraw=raw,nuisance_rank=rank,q=fullq,r=fullr,h=h,df=n-D.shape[1],condition=float(s[0]/s[-1]),max_leverage=float(h.max()),span_error=error),'ESTIMABLE'


def evaluate(model,y=None):
    yy=model['y'] if y is None else y;coef=np.linalg.solve(model['r'],model['q'].T@yy)
    target=np.zeros(model['D'].shape[1]);target[-1]=1
    row=np.linalg.solve(model['r'].T,target)@model['q'].T
    residual=yy-model['q']@(model['q'].T@yy);se=float(np.sqrt(np.sum((row*residual/(1-model['h']))**2)))
    if not np.isfinite(se) or se<=0:raise ValueError('DEGENERATE_HC3_VARIANCE')
    beta=float(coef[-1]);stat=float(beta/se)
    if model['kind']=='coordination' and y is None:
        rx=model['x']-model['Q']@(model['Q'].T@model['x']);ry=yy-model['Q']@(model['Q'].T@yy)
        effect=float(np.sum(rx*ry)/np.sqrt(np.sum(rx**2)*np.sum(ry**2)))
    else:effect=beta
    return dict(beta=beta,t=stat,se=se,effect=effect,row=row,y=yy)


def wild(model,draws,random_seed,y=None):
    observed=evaluate(model,y=y);yy=observed['y'];qn,_=np.linalg.qr(model['D'][:,:-1],mode='reduced')
    fitted=qn@(qn.T@yy);e=(yy-fitted)/(1-np.sum(qn**2,axis=1));row=observed['row'];rng=np.random.default_rng(random_seed);values=np.empty(draws)
    for start in range(0,draws,1000):
        signs=np.where(rng.random((min(1000,draws-start),len(yy)))<.5,-1.,1.).T
        Y=fitted[:,None]+e[:,None]*signs;residual=Y-model['q']@(model['q'].T@Y)
        b=row@Y;se=np.sqrt(np.sum((row[:,None]*residual/(1-model['h'][:,None]))**2,axis=0));values[start:start+Y.shape[1]]=b/se
    assert np.isfinite(values).all()
    return observed,values


def boot(x,y,C,names,kind,takes,expected_rank):
    values=np.full(len(takes),np.nan);reasons=[]
    for i,ix in enumerate(takes):
        m,status=prepare(x[ix],y[ix],C[ix],names,kind,expected_rank=expected_rank)
        if m is not None:
            try:values[i]=evaluate(m)['effect']
            except ValueError as e:status=str(e)
        reasons.append(status)
    return values,np.asarray(reasons)
