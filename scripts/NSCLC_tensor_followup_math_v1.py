"""Supplementary fixed-design full-model wild-t intervals; no primary P change."""
import hashlib
import numpy as np


def seed(slot,stream):
    return int.from_bytes(hashlib.sha256(('2026101926|'+slot+'|'+stream).encode()).digest()[:8],'little')


def full_model_wild_t(g,draws,random_seed):
    D=g['d'];inverse=g['inverse'];residual=g['y']-D@(inverse@g['y']);adjusted=residual/(1-g['h']);target=inverse[-1];student=target/(1-g['h']);rng=np.random.default_rng(random_seed);out=np.empty(draws)
    for start in range(0,draws,1000):
        signs=np.where(rng.random((min(1000,draws-start),len(residual)))<.5,-1.,1.).T
        noise=adjusted[:,None]*signs;errors=noise-D@(inverse@noise);se=np.sqrt(np.sum((student[:,None]*errors)**2,axis=0));out[start:start+noise.shape[1]]=(target@noise)/se
    assert np.isfinite(out).all()
    return out


def studentized_interval(g,values):
    lo,hi=np.quantile(values,[.025,.975]);se=np.sqrt(g['variance'])
    return np.array([g['beta']-hi*se,g['beta']-lo*se])
