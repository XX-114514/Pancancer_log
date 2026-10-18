"""Prospective paired-tissue rank statistics, exact signs and shared case bootstrap."""
import hashlib
import numpy as np
from scipy.stats import rankdata,binomtest
from scipy.stats import beta as beta_distribution


def seed(label):return int.from_bytes(hashlib.sha256(('2026101906|'+label).encode()).digest()[:8],'little')


def wilcoxon_signed(delta,exact_limit=20,draws=19999):
    d=np.asarray(delta,float);nonzero=d[d!=0];r=rankdata(np.abs(nonzero));observed=float(np.sign(nonzero)@r)
    if not len(r):return {'effect':float(np.median(d)),'statistic':0.,'P':1.,'exact':True,'statistics':np.zeros(1),'nonzero':0,'sign_P':1.}
    if len(r)<=exact_limit:
        code=np.arange(2**len(r),dtype=np.uint64)[:,None];signs=2*((code>>np.arange(len(r),dtype=np.uint64))&1).astype(np.int8)-1;statistics=signs@r;p=float((np.abs(statistics)>=abs(observed)-1e-12).mean());exact=True
    else:
        signs=2*np.random.default_rng(seed('long_signedrank')).integers(0,2,size=(draws,len(r)))-1;statistics=signs@r;p=float((1+(np.abs(statistics)>=abs(observed)-1e-12).sum())/(draws+1));exact=False
    return {'effect':float(np.median(d)),'statistic':observed,'P':p,'exact':exact,'statistics':statistics,'nonzero':len(r),'sign_P':float(binomtest(int((nonzero>0).sum()),len(r),.5).pvalue)}


def rho(x,y):
    if min(np.ptp(x),np.ptp(y))==0:return np.nan
    a=rankdata(x);b=rankdata(y);return float(np.corrcoef(a,b)[0,1])


def permutations(x,y,label,draws=19999):
    a=rankdata(x);b=rankdata(y);a-=a.mean();b-=b.mean();den=np.linalg.norm(a)*np.linalg.norm(b);assert den>0;rng=np.random.default_rng(seed(label+'|cochange'));out=[]
    for start in range(0,draws,1000):
        ix=np.argsort(rng.random((min(1000,draws-start),len(x))),axis=1);out.append(b[ix]@a/den)
    return np.concatenate(out)


def mc_probability(statistics,observed,streams=2):
    n=len(statistics);count=int((np.abs(statistics)>=abs(observed)-1e-12).sum());p=(1+count)/(1+n);u=1. if count==n else beta_distribution.ppf(1-.01/streams,count+1,n-count);return float(p),min(1.,float((1+n*u)/(1+n))),count


def bootstrap_indices(n):return np.random.default_rng(seed('all_primary_paired_cases')).integers(0,n,size=(2000,n)).astype(np.int16)


def bootstrap_effects(x,y,ix,cochange=False):return np.array([rho(x[j],y[j]) if cochange else np.median(y[j]) for j in ix])


def adjust(p,BY=False):
    p=np.asarray(p,float);order=np.argsort(p);n=len(p);v=p[order]*n/np.arange(1,n+1)
    if BY:v*=np.sum(1/np.arange(1,n+1))
    result=np.empty(n);result[order]=np.minimum(1,np.minimum.accumulate(v[::-1])[::-1]);return result
