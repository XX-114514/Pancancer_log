"""Conservative parent compatibility bounds from existing unit marginals.

No pseudo-cell reconstruction or exact confusion matrix from marginal counts.
The reference is an author's named group, not established biological truth.
"""
import numpy as np


def regulatory_bounds(n, foxp3, il2ra, cd3d, cd3e, trac, cd8a):
    arrays=[np.asarray(v,np.int64) for v in [n,foxp3,il2ra,cd3d,cd3e,trac,cd8a]]
    n,foxp3,il2ra,cd3d,cd3e,trac,cd8a=arrays
    assert all(a.shape==n.shape for a in arrays)
    assert (n>=0).all() and all(((a>=0)&(a<=n)).all() for a in arrays[1:])
    t_lower=np.maximum.reduce([cd3d,cd3e,trac])
    t_upper=np.minimum(n,cd3d+cd3e+trac)
    lower=np.maximum(0,foxp3+il2ra+t_lower-cd8a-2*n)
    upper=np.minimum.reduce([foxp3,il2ra,t_upper,n-cd8a])
    assert (lower<=upper).all()
    return lower,upper


def reference_compatibility(lower,upper,reference,observed_total):
    lower=np.asarray(lower,np.int64);upper=np.asarray(upper,np.int64);reference=np.asarray(reference,bool)
    assert lower.shape==upper.shape==reference.shape and (0<=lower).all() and (lower<=upper).all()
    assert lower.sum()<=observed_total<=upper.sum()
    lo=max(0,int(lower[reference].sum()),int(observed_total-upper[~reference].sum()))
    hi=min(int(observed_total),int(upper[reference].sum()),int(observed_total-lower[~reference].sum()))
    assert 0<=lo<=hi<=observed_total
    return dict(observed_parent_total=int(observed_total),compatible_reference_count_lower=lo,
                compatible_reference_count_upper=hi,
                compatible_reference_fraction_lower=lo/observed_total if observed_total else None,
                compatible_reference_fraction_upper=hi/observed_total if observed_total else None,
                exact_same_cell_confusion_matrix_reconstructed=False)
