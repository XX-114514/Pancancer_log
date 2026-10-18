"""Prespecified native interface statistic, full nuisance span and streams."""
import hashlib

import numpy as np

from NSCLC_external_programme_coordination_math_v2 import prepare as rank_prepare,effect,permutation_count
from preflight_NSCLC_continuous_prior_native_association_v1 import design_gate


def seed(record,stream):
    return int.from_bytes(hashlib.sha256((record+'|'+stream+'|native_continuous_prior_association_20261016_v1').encode()).digest()[:8],'little')


def prepare(a,b,quality,blocks,minimum=50):
    p=rank_prepare(a,b,quality,blocks,min_n=minimum,min_df=10)
    design_gate(p)
    return p


def bootstrap(p,draws,random_seed):
    rng=np.random.default_rng(random_seed);groups=[np.flatnonzero(p['blocks']==b) for b in sorted(set(p['blocks']))]
    out=np.full(draws,np.nan)
    for i in range(draws):
        ix=np.concatenate([g[rng.integers(len(g),size=len(g))] for g in groups])
        try:out[i]=effect(prepare(p['a'][ix],p['b'][ix],p['quality'][ix],p['blocks'][ix]))
        except ValueError:pass
    return out


def leave_one(p):
    out=np.full(len(p['a']),np.nan)
    for i in range(len(out)):
        keep=np.arange(len(out))!=i
        try:out[i]=effect(prepare(p['a'][keep],p['b'][keep],p['quality'][keep],p['blocks'][keep],minimum=3))
        except ValueError:pass
    return out


def secondary(a,b,quality,blocks):
    try:return effect(prepare(a,b,quality,blocks)),'ESTIMABLE'
    except ValueError as error:return np.nan,str(error)
