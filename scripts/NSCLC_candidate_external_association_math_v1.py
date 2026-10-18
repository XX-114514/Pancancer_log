"""Fixed external rank statistic and streams, complete nuisance-span guards."""
import hashlib
import numpy as np

import NSCLC_external_programme_coordination_math_v2 as base
from preflight_NSCLC_continuous_prior_native_association_v1 import design_gate


def prepare(x,y,q,blocks,minimum=20):
    p=base.prepare(x,y,q,blocks,min_n=minimum,min_df=8)
    design_gate(p)
    return p


def seed(model,stream):
    return int.from_bytes(hashlib.sha256((model+'|'+stream+'|external_fixed_native_ligand_programme_20261018_v1').encode()).digest()[:8],'little')


def bootstrap(p,draws,random_seed):
    rng=np.random.default_rng(random_seed);groups=[np.flatnonzero(p['blocks']==g) for g in sorted(set(p['blocks']))]
    values=np.full(draws,np.nan)
    for j in range(draws):
        ix=np.concatenate([g[rng.integers(len(g),size=len(g))] for g in groups])
        try:values[j]=base.effect(prepare(p['a'][ix],p['b'][ix],p['quality'][ix],p['blocks'][ix]))
        except ValueError:pass
    return values


def leave_one(p):
    values=np.full(len(p['a']),np.nan)
    for j in range(len(values)):
        ix=np.arange(len(values))!=j
        try:values[j]=base.effect(prepare(p['a'][ix],p['b'][ix],p['quality'][ix],p['blocks'][ix],minimum=3))
        except ValueError:pass
    return values


def secondary(x,y,q,blocks):
    try:return base.effect(prepare(x,y,q,blocks)),'ESTIMABLE'
    except ValueError as error:return np.nan,str(error)
