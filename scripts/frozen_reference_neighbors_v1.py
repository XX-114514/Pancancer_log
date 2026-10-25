"""Existing approximate neighbor runtime with explicit self exclusion and probes."""
import numpy as np
import pynndescent
from scipy.spatial.distance import cdist


class FrozenReferenceNeighbors:
    def __init__(self,reference,seed,parameters):
        self.reference=np.asarray(reference,dtype=np.float32)
        if self.reference.ndim!=2 or len(self.reference)<2 or not np.isfinite(self.reference).all():raise ValueError('Invalid reference neighbor axes')
        self.parameters=parameters
        self.index=pynndescent.NNDescent(self.reference.copy(),metric='euclidean',n_neighbors=parameters['n_neighbors'],n_trees=parameters['n_trees'],n_search_trees=parameters['n_search_trees'],n_iters=parameters['n_iters'],max_candidates=parameters['max_candidates'],random_state=seed,low_memory=True,n_jobs=4,parallel_batch_queries=True,verbose=False)
    def distances(self,query,excluded_reference,block_rows=4096):
        query=np.asarray(query,dtype=np.float32);excluded=np.asarray(excluded_reference,np.int64)
        if query.shape!=(len(excluded),self.reference.shape[1]) or not np.isfinite(query).all():raise ValueError('Query neighbor axes or finite domain differ')
        out=np.empty(len(query),np.float32)
        for low in range(0,len(query),block_rows):
            z=query[low:low+block_rows];j,d=self.index.query(z,k=2,epsilon=self.parameters['query_epsilon'])
            if np.any(j<0) or np.any(j>=len(self.reference)) or not np.isfinite(d).all():raise ValueError('Approximate neighbor search returned invalid rows')
            for i in range(len(z)):
                take=j[i]!=excluded[low+i]
                if not take.any():raise ValueError('Self exclusion left no other reference neighbor')
                out[low+i]=d[i,take].min()
        return out
    def exact_probe_check(self,query,excluded_reference,approximate):
        exact=np.empty(len(query),float)
        for low in range(0,len(query),16):
            z=np.asarray(query[low:low+16],float);d=cdist(z,self.reference.astype(float))
            for i,exclude in enumerate(excluded_reference[low:low+16]):
                if exclude>=0:d[i,int(exclude)]=np.inf
            exact[low:low+len(z)]=d.min(axis=1)
        approximate=np.asarray(approximate,float);error=approximate-exact
        if np.any(error<-(2e-4+2e-5*exact)):raise ValueError('Approximate distance below independent exact distance')
        matched=np.isclose(approximate,exact,atol=2e-4,rtol=2e-5)
        return {'n_exact_probes':len(query),'distance_exact_match_fraction':float(matched.mean()) if len(query) else None,'maximum_distance_overestimate':float(error.max()) if len(query) else None,'median_distance_overestimate':float(np.median(error)) if len(query) else None}


def actual_runtime_smoke(seed,parameters):
    rng=np.random.default_rng(seed);center=rng.normal(size=(4,12))*5;reference=np.vstack([center[j]+rng.normal(size=(512,12))*.5 for j in range(4)]).astype(np.float32);query=np.vstack([reference[:64],center[0]+rng.normal(size=(64,12))*.5]).astype(np.float32);exclude=np.r_[np.arange(64),np.full(64,-1)]
    search=FrozenReferenceNeighbors(reference,seed,parameters);distance=search.distances(query,exclude);report=search.exact_probe_check(query,exclude,distance)
    if report['distance_exact_match_fraction']<parameters['minimum_exact_probe_fraction']:raise ValueError('Existing neighbor runtime failed frozen synthetic exact-distance gate')
    return report
