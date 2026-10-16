"""Sparse binary-membership ULM score, positive RNA support and mean diagnostics."""
import numpy as np
from scipy.sparse import csr_matrix


def score(expression,targets):
    x=np.asarray(expression,dtype=np.float64)
    assert x.ndim==2 and np.isfinite(x).all();n,g=x.shape
    rr=[];cc=[]
    for i,indices in enumerate(targets):
        indices=np.asarray(indices,dtype=int);assert len(indices)>0 and len(np.unique(indices))==len(indices)
        assert indices.min()>=0 and indices.max()<g;rr.extend([i]*len(indices));cc.extend(indices)
    w=csr_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(targets),g));k=np.asarray(w.sum(axis=1)).ravel()
    assert np.all((k>0)&(k<g));mu=x.mean(axis=1);centered=x-mu[:,None];ss=np.einsum('ij,ij->i',centered,centered)
    numerator=np.asarray(w@centered.T);denominator=np.sqrt((k-k*k/g)[:,None]*ss[None,:])
    with np.errstate(invalid='ignore',divide='ignore'):
        corr=numerator/denominator
        values=corr*np.sqrt((g-2)/((1-corr+1e-16)*(1+corr+1e-16)))
    sums=np.asarray(w@x.T)
    target_mean=(sums/k[:,None]).T
    background_mean=((g*mu[None,:]-sums)/(g-k[:,None])).T
    detected=np.asarray(w@(x>0).T).T.astype(np.int16)
    return dict(ULM_score=values.T,target_mean=target_mean,nontarget_mean=background_mean,expressed_target_genes=detected)
