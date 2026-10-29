"""Reproducible display sampling preserves observed rare clusters and masks."""
import numpy as np
def balanced_indices(codes,eligible,cap,minimum,seed):
 codes=np.asarray(codes);eligible=np.asarray(eligible,bool)
 if codes.ndim!=1 or eligible.shape!=codes.shape or cap<1 or minimum<1:raise ValueError('Invalid display sampling contract')
 rng=np.random.default_rng(seed);base=np.flatnonzero(eligible);chosen=[]
 for k in np.unique(codes[base]):
  q=base[codes[base]==k];chosen.extend(rng.choice(q,min(minimum,len(q)),replace=False).tolist())
 chosen=np.asarray(chosen,np.int64);remaining=np.setdiff1d(base,chosen,assume_unique=False);budget=min(len(base),max(cap,len(chosen)))-len(chosen)
 if budget:chosen=np.concatenate([chosen,rng.choice(remaining,budget,replace=False)])
 return np.sort(chosen)
