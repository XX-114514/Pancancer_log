"""Fixed-loading own-view ridge projection; no other cell view is an input."""
import hashlib
from pathlib import Path
import numpy as np


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for b in iter(lambda:stream.read(16777216),b''):h.update(b)
    return h.hexdigest()


def ridge_score(expression,weights,precision):
    x=np.asarray(expression,float);w=np.asarray(weights,float);t=np.asarray(precision,float)
    assert x.ndim==w.ndim==2 and x.shape[1]==w.shape[0] and t.shape==(w.shape[0],)
    assert np.isfinite(x).all() and np.isfinite(w).all() and np.isfinite(t).all() and (t>0).all()
    gram=np.eye(w.shape[1])+w.T@(t[:,None]*w)
    return np.linalg.solve(gram,w.T@(t[:,None]*x.T)).T
