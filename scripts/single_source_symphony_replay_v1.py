"""Independent sufficient-statistic replay for one complete query source.

The query design has intercept and one source indicator (both columns one).
A closed-form ridge coefficient is fitted once from all included rows. Output
chunks apply the same coefficient and never perform independent batch refits.
"""
import numpy as np


def responsibilities(x,centers,sigma):
    x=np.asarray(x);centers=np.asarray(centers)
    if x.ndim!=2 or centers.ndim!=2 or x.shape[1]!=centers.shape[1] or not np.isfinite(x).all():raise ValueError('Invalid query/center axes')
    maximum=x.max(axis=1,keepdims=True);norms=np.linalg.norm(centers,axis=1,keepdims=True)
    if np.any(maximum<=0) or np.any(norms<=0):raise ValueError('Unsupported cosine projection/center')
    z=x/maximum;z/=np.linalg.norm(z,axis=1,keepdims=True)
    scales=np.asarray(sigma,dtype=np.float32).ravel()
    if len(scales)==1:scales=np.repeat(scales,len(centers))
    if scales.shape!=(len(centers),) or np.any(scales<=0):raise ValueError('Invalid soft-cluster bandwidth')
    logits=(-2*(1-(centers/norms)@z.T)/scales[:,None]);logits-=logits.max(axis=0,keepdims=True);weights=np.exp(logits);weights/=weights.sum(axis=0,keepdims=True)
    return weights


def fit_single_source_offsets(raw,centers,reference_mass,sigma,lamb=1.,block_rows=8192):
    raw=np.asarray(raw);centers=np.asarray(centers);nr=np.asarray(reference_mass,float)
    if block_rows<1 or len(raw)<1 or lamb<=0 or nr.shape!=(len(centers),):raise ValueError('Invalid complete-source ridge domain')
    mass=np.zeros(len(centers),np.float64);weighted=np.zeros(centers.shape,np.float64)
    for low in range(0,len(raw),block_rows):
        x=raw[low:low+block_rows];r=responsibilities(x,centers,sigma).astype(np.float64);mass+=r.sum(axis=1);weighted+=r@x
    denominator=nr*mass+lamb*(nr+mass)
    if np.any(denominator<=0):raise ValueError('Unidentified source ridge coefficient')
    offsets=(nr[:,None]*weighted-mass[:,None]*centers)/denominator[:,None]
    return offsets


def apply_single_source_offsets(raw,centers,sigma,offsets,block_rows=8192):
    out=np.empty_like(raw,dtype=np.float64)
    for low in range(0,len(raw),block_rows):
        x=raw[low:low+block_rows];r=responsibilities(x,centers,sigma);out[low:low+len(x)]=x-r.T@offsets
    return out
