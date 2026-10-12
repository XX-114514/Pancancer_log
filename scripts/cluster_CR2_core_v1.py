"""Unweighted OLS CR2, identity working covariance, coefficient-specific Satt df.
Original nuisance span retained. Small-cluster approximation, not exact inference.
"""
import numpy as np
from scipy.stats import t as student_t

def standardize(a):
 a=np.asarray(a,float);s=a.std(axis=0);return (a-a.mean(axis=0))/np.where(s>1e-12,s,1)

def prepare_cr2(Z,x,blocks,min_df=4.0):
 Z=np.asarray(Z,float);x=np.asarray(x,float);blocks=np.asarray(blocks);U,sv,_=np.linalg.svd(Z,full_matrices=False);rank=int((sv>sv[0]*1e-10).sum());Q=U[:,:rank];xr=x-Q@(Q.T@x);ss=float(xr@xr);n=len(x)
 if ss<=1e-10*n:return dict(status='exposure_unidentified')
 q=xr/np.sqrt(ss);fullQ=np.column_stack([Q,q]);H=fullQ@fullQ.T;M=np.eye(n)-H;c=xr/ss;levels=np.unique(blocks);B=np.zeros((n,len(levels)));max_cluster_hat=0.
 for j,b in enumerate(levels):
  ix=np.flatnonzero(blocks==b);Hgg=H[np.ix_(ix,ix)];e,V=np.linalg.eigh(np.eye(len(ix))-Hgg);max_cluster_hat=max(max_cluster_hat,float(1-e.min()))
  if e.min()<=1e-8:return dict(status='cluster_adjustment_singular',max_cluster_hat=max_cluster_hat)
  Ag=(V/np.sqrt(e))@V.T;B[ix,j]=Ag@c[ix]
 P=M@B;C=P.T@P;trace=float(np.trace(C));df=trace*trace/float(np.sum(C*C));return dict(status='inference_estimable' if df>=min_df else 'effective_df_insufficient',df=df,max_cluster_hat=max_cluster_hat,coef_vector=c,residual_operator=M,adjusted_cluster_vectors=P,rank=rank+1,n=n,n_blocks=len(levels))

def estimate_cr2(prep,Y):
 Y=np.asarray(Y,float)
 if Y.ndim==1:Y=Y[:,None]
 beta=prep['coef_vector']@Y;scores=prep['adjusted_cluster_vectors'].T@Y;variance=np.sum(scores*scores,axis=0);se=np.sqrt(variance);df=prep['df'];stats=np.divide(beta,se,out=np.full_like(beta,np.nan),where=se>1e-14);p=2*student_t.sf(abs(stats),df);critical=student_t.ppf(.975,df);return dict(beta=beta,se=se,df=df,stat=stats,p=p,lower=beta-critical*se,upper=beta+critical*se)

def bh(p):
 p=np.asarray(p,float);m=len(p);order=np.argsort(p,kind='stable');s=p[order]*m/np.arange(1,m+1);q=np.empty(m);q[order]=np.minimum(1,np.minimum.accumulate(s[::-1])[::-1]);return q
