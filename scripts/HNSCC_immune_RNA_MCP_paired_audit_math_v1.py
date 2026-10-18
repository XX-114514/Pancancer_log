"""Separate pandas ranks, pivoted QR and statsmodels HC3 numerical reference."""
import numpy as np
import pandas as pd
from scipy.linalg import qr
import statsmodels.api as sm


def reference_model(x,y,C,names,expected_rank=None):
    x=np.asarray(x,float);y=np.asarray(y,float);C=np.asarray(C,float);n=len(x)
    if n<20:return None,'FEWER_THAN20_COMPLETE_PAIRED_PATIENTS'
    if not np.isfinite(np.column_stack([x,y,C])).all():return None,'NONFINITE_INPUT'
    if min(np.ptp(x),np.ptp(y))==0:return None,'CONSTANT_PAIRED_TARGET_OR_READOUT'
    cc=pd.DataFrame(C).copy()
    for j,name in enumerate(names):
        if name not in ['treatment_DplusT','HPV_positive']:cc[j]=cc[j].rank(method='average')
    cc=cc.to_numpy(float);x=pd.Series(x).rank().to_numpy();y=pd.Series(y).rank().to_numpy();x=(x-x.mean())/x.std(ddof=0);y=(y-y.mean())/y.std(ddof=0)
    sd=np.std(cc,axis=0);Z=np.column_stack([np.ones(n),(cc-cc.mean(axis=0))/np.where(sd>0,sd,1)])
    singular=np.linalg.svd(Z,compute_uv=False);r=int(np.linalg.matrix_rank(Z));stable=int((singular>singular[0]*1e-12).sum())
    if r!=stable or singular[r-1]/singular[0]<1e-8:return None,'UNSTABLE_COMPLETE_NUISANCE_SPAN'
    if expected_rank is not None and r!=expected_rank:return None,'BOOTSTRAP_NUISANCE_RANK_CHANGED'
    Q,R,p=qr(Z,mode='economic',pivoting=True);Q=Q[:,:r]
    raw=np.column_stack([np.ones(n),cc]);err=float(np.abs(raw-Q@Q.T@raw).max())
    if err>1e-8*max(1,np.abs(raw).max()):return None,'NUISANCE_SPAN_RECONSTRUCTION_FAILED'
    rx=x-Q@Q.T@x;ry=y-Q@Q.T@y
    if np.linalg.norm(rx)<=1e-8*max(1,np.linalg.norm(x)):return None,'TARGET_NOT_IDENTIFIABLE'
    if np.linalg.norm(ry)<=1e-8*max(1,np.linalg.norm(y)):return None,'READOUT_NOT_IDENTIFIABLE'
    D=np.column_stack([Q,x]);s=np.linalg.svd(D,compute_uv=False)
    if n-D.shape[1]<10:return None,'FEWER_THAN10_RESIDUAL_DF'
    if s[0]/s[-1]>=1e8:return None,'ILL_CONDITIONED_TARGET_DESIGN'
    qd,rd=np.linalg.qr(D,mode='reduced');h=np.sum(qd**2,axis=1)
    if h.max()>=.999999:return None,'EXTREME_LEVERAGE'
    return dict(D=D,Q=Q,x=x,y=y,nuisance_rank=r,condition=float(s[0]/s[-1]),max_leverage=float(h.max()),residual_df=n-D.shape[1]),'ESTIMABLE'


def reference_estimate(m):
    g=sm.OLS(m['y'],m['D']).fit(cov_type='HC3')
    rx=m['x']-m['Q']@m['Q'].T@m['x'];ry=m['y']-m['Q']@m['Q'].T@m['y']
    return float(g.params[-1]),float(g.tvalues[-1]),float(np.corrcoef(rx,ry)[0,1])


def reference_wild(y,D,draws,random_seed):
    Q,R=np.linalg.qr(D,mode='reduced');inverse=np.linalg.solve(R,Q.T);h=np.sum(Q*Q,axis=1)
    q0,r0=np.linalg.qr(D[:,:-1],mode='reduced');mu=q0@(q0.T@y);e=(y-mu)/(1-np.sum(q0*q0,axis=1))
    signs=np.where(np.random.default_rng(random_seed).random((draws,len(y)))<.5,-1.,1.).T
    noise=e[:,None]*signs;res=noise-Q@(Q.T@noise)
    return (inverse[-1]@noise)/np.sqrt(np.sum(((inverse[-1]/(1-h))[:,None]*res)**2,axis=0))


def reference_bootstrap(x,y,C,names,indices,rank):
    vals=np.full(len(indices),np.nan);statuses=[]
    for k,ix in enumerate(indices):
        m,status=reference_model(x[ix],y[ix],C[ix],names,expected_rank=rank)
        if m is not None:
            rx=m['x']-m['Q']@m['Q'].T@m['x'];ry=m['y']-m['Q']@m['Q'].T@m['y'];vals[k]=np.corrcoef(rx,ry)[0,1]
        statuses.append(status)
    return vals,np.asarray(statuses)
