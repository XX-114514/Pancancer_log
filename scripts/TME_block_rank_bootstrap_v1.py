"""Exact Spearman correlations for integer cluster-bootstrap multiplicities."""
import numpy as np

def weighted_midranks(values,weights):
 values=np.asarray(values,float);weights=np.asarray(weights,float)
 assert values.ndim==1 and weights.ndim==2 and weights.shape[1]==len(values)
 assert np.isfinite(values).all() and np.isfinite(weights).all() and (weights>=0).all()
 order=np.argsort(values,kind='stable');s=values[order];starts=np.r_[0,np.flatnonzero(np.diff(s)!=0)+1];sizes=np.diff(np.r_[starts,len(values)])
 group=np.add.reduceat(weights[:,order],starts,axis=1);rank=np.cumsum(group,axis=1)-group/2+.5
 expanded=np.repeat(rank,sizes,axis=1);out=np.empty_like(expanded);out[:,order]=expanded;return out

def weighted_spearman(values_x,values_y,weights):
 w=np.asarray(weights,float);total=w.sum(1);rx=weighted_midranks(values_x,w);ry=weighted_midranks(values_y,w)
 mx=np.divide((rx*w).sum(1),total,out=np.full(len(total),np.nan),where=total>0);my=np.divide((ry*w).sum(1),total,out=np.full(len(total),np.nan),where=total>0)
 dx=rx-mx[:,None];dy=ry-my[:,None];den=np.sqrt((w*dx*dx).sum(1)*(w*dy*dy).sum(1))
 return np.divide((w*dx*dy).sum(1),den,out=np.full(len(total),np.nan),where=den>0)

def block_multiplicities(block_codes,repeats,seed):
 block=np.asarray(block_codes);unique,inverse=np.unique(block,return_inverse=True);nb=len(unique);rng=np.random.default_rng(seed)
 draws=rng.integers(0,nb,(repeats,nb));counts=np.zeros((repeats,nb),np.int32)
 np.add.at(counts,(np.arange(repeats)[:,None],draws),1)
 return counts[:,inverse]
