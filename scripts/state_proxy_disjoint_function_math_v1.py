"""Two implementations of block-weighted, identity-stratified rank descriptions."""
import numpy as np
import pandas as pd
from scipy.stats import rankdata


def correlate(a,b,w):
    cov=np.sum(w[:,None]*a*b,axis=0)
    den=np.sqrt(np.sum(w[:,None]*a*a,axis=0)*np.sum(w[:,None]*b*b,axis=0))
    return np.divide(cov,den,out=np.full(b.shape[1],np.nan),where=den>1e-14)


def fit(proxy,functions,quality,blocks):
    _,inv,n=np.unique(blocks,return_inverse=True,return_counts=True);w=1/(len(n)*n[inv]);sw=np.sqrt(w)
    q=rankdata(quality,axis=0);q-=np.sum(w[:,None]*q,axis=0);sd=np.sqrt(np.sum(w[:,None]*q*q,axis=0))
    if (sd<=1e-12).any():return None,'CONSTANT_REQUIRED_QUALITY'
    q/=sd;X=np.column_stack([np.ones(len(w)),q]);Q,R=np.linalg.qr(sw[:,None]*X,mode='reduced')
    if np.linalg.matrix_rank(R)<X.shape[1] or np.linalg.cond(R)>=1e8:return None,'REQUIRED_FOUR_QUALITY_DESIGN_UNIDENTIFIABLE'
    a=rankdata(proxy)[:,None];b=rankdata(functions,axis=0)
    ac=a-np.sum(w[:,None]*a,axis=0);bc=b-np.sum(w[:,None]*b,axis=0);raw=correlate(ac,bc,w)
    ar=(sw[:,None]*a-Q@(Q.T@(sw[:,None]*a)))/sw[:,None];br=(sw[:,None]*b-Q@(Q.T@(sw[:,None]*b)))/sw[:,None]
    av=np.sum(w[:,None]*ar*ar,axis=0);bv=np.sum(w[:,None]*br*br,axis=0)
    keep=(np.sum(w[:,None]*ac*ac,axis=0)>1e-12)&(np.sum(w[:,None]*bc*bc,axis=0)>1e-12)&(av>np.sum(w[:,None]*ac*ac,axis=0)*1e-10)&(bv>np.sum(w[:,None]*bc*bc,axis=0)*1e-10)
    cell=correlate(ar,br,w);cell[~keep]=np.nan
    am=np.array([ar[inv==i].mean(0) for i in range(len(n))]);bm=np.array([br[inv==i].mean(0) for i in range(len(n))])
    varying=(am.var(0)>av*1e-10)&(bm.var(0)>bv*1e-10)
    am=rankdata(am,axis=0);bm=rankdata(bm,axis=0);am-=am.mean(0);bm-=bm.mean(0)
    block=correlate(am,bm,np.full(len(n),1/len(n)));block[(~keep)|(~varying)]=np.nan
    return dict(raw=raw,cell=cell,block=block),'ESTIMABLE_DESCRIPTION'


def independent_fit(proxy,functions,quality,blocks):
    labels=pd.Series(blocks);sizes=labels.value_counts();w=labels.map(1/(len(sizes)*sizes)).to_numpy();sw=np.sqrt(w)
    q=pd.DataFrame(quality).rank().to_numpy();q-=np.average(q,axis=0,weights=w);sd=np.sqrt(np.average(q*q,axis=0,weights=w))
    if np.any(sd<=1e-12):return None,'CONSTANT_REQUIRED_QUALITY'
    X=np.column_stack([np.ones(len(w)),q/sd]);Z=sw[:,None]*X
    if np.linalg.matrix_rank(Z)<X.shape[1] or np.linalg.cond(Z)>=1e8:return None,'REQUIRED_FOUR_QUALITY_DESIGN_UNIDENTIFIABLE'
    a=pd.Series(proxy).rank().to_numpy();b=pd.DataFrame(functions).rank().to_numpy();joined=np.column_stack([a,b])
    residual=joined-X@np.linalg.lstsq(Z,sw[:,None]*joined,rcond=None)[0]
    ar=residual[:,0];br=residual[:,1:];av=np.average(ar*ar,weights=w);bv=np.average(br*br,axis=0,weights=w)
    ac=a-np.average(a,weights=w);bc=b-np.average(b,axis=0,weights=w);oa=np.average(ac*ac,weights=w);ob=np.average(bc*bc,axis=0,weights=w)
    keep=(oa>1e-12)&(ob>1e-12)&(av>oa*1e-10)&(bv>ob*1e-10)
    raw=np.array([np.average(ac*bc[:,j],weights=w)/np.sqrt(oa*ob[j]) if oa*ob[j]>1e-28 else np.nan for j in range(b.shape[1])])
    cell=np.array([np.average(ar*br[:,j],weights=w)/np.sqrt(av*bv[j]) if keep[j] else np.nan for j in range(b.shape[1])])
    group=pd.DataFrame(residual).groupby(np.asarray(blocks),sort=True).mean().to_numpy();vary=(group[:,0].var()>av*1e-10)&(group[:,1:].var(0)>bv*1e-10)
    ranked=pd.DataFrame(group).rank().to_numpy();ranked-=ranked.mean(0);block=np.full(b.shape[1],np.nan)
    for j in np.flatnonzero(keep&vary):block[j]=np.dot(ranked[:,0],ranked[:,j+1])/np.sqrt(np.dot(ranked[:,0],ranked[:,0])*np.dot(ranked[:,j+1],ranked[:,j+1]))
    return dict(raw=raw,cell=cell,block=block),'ESTIMABLE_DESCRIPTION'


def fixture():
    rng=np.random.default_rng(2026102006);blocks=np.repeat(np.arange(12),np.arange(5,17));n=len(blocks);q=rng.normal(size=(n,4));p=rng.normal(size=n)+.8*q[:,0]
    y=np.column_stack([p+rng.normal(size=n)*.2,-p+rng.normal(size=n)*.2,np.ones(n)])
    a,sa=fit(p,y,q,blocks);b,sb=independent_fit(p,y,q,blocks);assert sa==sb
    for key in a:np.testing.assert_allclose(a[key],b[key],atol=2e-9,rtol=1e-8,equal_nan=True)
    assert a['cell'][0]>.5 and a['cell'][1]<-.5 and np.isnan(a['cell'][2])
    assert fit(p,y,np.ones((n,4)),blocks)[1]=='CONSTANT_REQUIRED_QUALITY'
    assert fit(p,y,np.column_stack([q[:,:3],q[:,0]]),blocks)[1]=='REQUIRED_FOUR_QUALITY_DESIGN_UNIDENTIFIABLE'
    return dict(status='PASS_POSITIVE_NEGATIVE_CONSTANT_DUPLICATE_COVARIATE_AND_INDEPENDENT_WEIGHTED_LS',synthetic_only=True)
