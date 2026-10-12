"""Null-imposed, block wild bootstrap joint max t. Approximate finite-sample inference.
All receiver TFs share each block weight. Complete design retained per replicate.
"""
import numpy as np

def prepare(Z,x,blocks,cr2):
 U,s,_=np.linalg.svd(Z,full_matrices=False);rank=int((s>s[0]*1e-10).sum());Q=U[:,:rank];Hz=Q@Q.T;groups=[np.flatnonzero(blocks==g) for g in np.unique(blocks)];adjustments=[]
 for ix in groups:
  v,A=np.linalg.eigh(np.eye(len(ix))-Hz[np.ix_(ix,ix)])
  if v.min()<=1e-8:return {'status':'null_block_adjustment_singular'}
  adjustments.append((A/np.sqrt(v))@A.T)
 return dict(status='prepared',Q=Q,groups=groups,adjustments=adjustments,c=cr2['coef_vector'],P=cr2['adjusted_cluster_vectors'],G=len(groups))

def parts(prep,Y):
 Q=prep['Q'];e0=Y-Q@(Q.T@Y);G=prep['G'];K=Y.shape[1];bp=np.zeros((G,K));sp=np.zeros((G,G,K))
 for g,(ix,A) in enumerate(zip(prep['groups'],prep['adjustments'])):
  E=A@e0[ix];bp[g]=prep['c'][ix]@E;sp[g]=prep['P'][ix].T@E
 return bp,sp

def bootstrap_max_statistics(prep,Y,weights):
 bp,sp=parts(prep,Y);b=weights@bp;s=(weights@sp.reshape(prep['G'],-1)).reshape(len(weights),prep['G'],Y.shape[1]);se=np.sqrt(np.sum(s*s,axis=1));t=np.divide(abs(b),se,out=np.zeros_like(b),where=se>1e-14);t[(se<=1e-14)&(abs(b)>1e-12)]=np.inf;return np.max(t,axis=1)

def test(prep,Y,observed_statistics,rng,draws=20000,exact_G_threshold=14,batch=256):
 G=prep['G'];K=Y.shape[1];obs=float(np.nanmax(abs(observed_statistics)));exact=G<=exact_G_threshold;B=2**(G-1) if exact else draws;bp,sp=parts(prep,Y);extreme=0;bad=0
 for start in range(0,B,batch):
  end=min(B,start+batch)
  if exact:
   ix=np.arange(start,end,dtype=np.uint64)[:,None];bits=(ix>>np.arange(G-1,dtype=np.uint64))&1;w=np.column_stack((np.ones(end-start),2*bits.astype(float)-1))
  else:w=rng.choice(np.array([-1.,1.]),size=(end-start,G));w[:,0]=1
  beta=w@bp;scores=(w@sp.reshape(G,-1)).reshape(end-start,G,K);se=np.sqrt(np.sum(scores*scores,axis=1));t=np.divide(abs(beta),se,out=np.zeros_like(beta),where=se>1e-14);t[(se<=1e-14)&(abs(beta)>1e-12)]=np.inf;max_t=np.max(t,axis=1);bad+=int((~np.isfinite(max_t)).sum());extreme+=int((max_t>=obs-1e-10*max(1,obs)).sum())
 p=extreme/B if exact else (extreme+1)/(B+1)
 return dict(p=p,draws=B,extreme=extreme,nonfinite_bootstrap_statistics=bad,enumerated_all_sign_pairs=exact,observed_max_abs_t=obs,finite_sample_method='approximate wild bootstrap; exact enumeration of weights not exact randomization inference')
