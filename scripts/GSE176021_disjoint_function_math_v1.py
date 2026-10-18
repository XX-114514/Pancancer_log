"""Prospective compact permutation mathematics for the complete 400-slot family."""
import hashlib
import numpy as np
from scipy.stats import rankdata, beta
from GSE176021_paired_tissue_math_v1 import adjust


def seed(label):
    return int.from_bytes(hashlib.sha256(('2026101908|'+label).encode()).digest()[:8],'little')


def doubled_centered_ranks(x,axis=None):
    x = np.asarray(x,float)
    n = x.size if axis is None else x.shape[axis]
    return np.rint(2*rankdata(x,method='average',axis=axis)-(n+1)).astype(np.int16)


def observed(x,y):
    a,b = doubled_centered_ranks(x),doubled_centered_ranks(y)
    numerator = int(np.dot(a.astype(np.int64),b.astype(np.int64)))
    denominator = float(np.linalg.norm(a)*np.linalg.norm(b))
    return (numerator/denominator if denominator else np.nan),numerator,denominator


def permutation_numerators(x,y,label,draws=99999):
    assert len(x) == len(y) == 10
    a,b = doubled_centered_ranks(x).astype(np.int64),doubled_centered_ranks(y).astype(np.int64)
    rng = np.random.default_rng(seed(label+'|full_case_permutation'))
    out = np.empty(draws,np.int16)
    for lo in range(0,draws,1000):
        hi = min(draws,lo+1000)
        order = np.argsort(rng.random((hi-lo,10)),axis=1)
        values = b[order]@a
        assert np.max(np.abs(values)) <= 330
        out[lo:hi] = values
    return out


def probability(values,numerator):
    n = len(values)
    count = int(np.count_nonzero(np.abs(values.astype(np.int64)) >= abs(numerator)))
    p = (count+1)/(n+1)
    upper = 1. if count == n else float(beta.ppf(1-.01/400,count+1,n-count))
    return float(p),min(1.,float((1+n*upper)/(n+1))),count


def bootstrap_indices():
    return np.random.default_rng(seed('all10_shared_case_pairs')).integers(0,10,size=(2000,10)).astype(np.int16)


def bootstrap(x,y,indices):
    a = doubled_centered_ranks(np.asarray(x)[indices],axis=1).astype(np.int64)
    b = doubled_centered_ranks(np.asarray(y)[indices],axis=1).astype(np.int64)
    num = np.sum(a*b,axis=1)
    den = np.sqrt(np.sum(a*a,axis=1)*np.sum(b*b,axis=1))
    return np.divide(num,den,out=np.full(len(indices),np.nan),where=den>0)


def response_quality_point(x,y,quality):
    """Descriptive partial Spearman with both response quality controls mandatory."""
    quality = np.asarray(quality,float)
    assert quality.shape == (10,2)
    D = np.c_[np.ones(10),rankdata(quality,axis=0)]
    if np.linalg.matrix_rank(D) != 3:
        return np.nan,'COLLINEAR_RESPONSE_QUALITY_NO_COVARIATE_DROPPED'
    a,b = rankdata(x),rankdata(y)
    ra = a-D@np.linalg.lstsq(D,a,rcond=None)[0]
    rb = b-D@np.linalg.lstsq(D,b,rcond=None)[0]
    if min(np.linalg.norm(ra),np.linalg.norm(rb)) < 1e-10:
        return np.nan,'ZERO_RESIDUAL_VARIANCE_NO_POINT'
    return float(np.dot(ra,rb)/(np.linalg.norm(ra)*np.linalg.norm(rb))),'POINT_ONLY_TWO_RESPONSE_QUALITY_CONTROLS'
