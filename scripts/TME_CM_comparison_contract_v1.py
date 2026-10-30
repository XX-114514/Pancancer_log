"""Measured composition and annotation-partition comparisons, no new CM model."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import coo_matrix

def frequency_with_measurement_mask(node_counts,node_parent,parent_counts,minimum_parent_cells):
    """Keep absent parents missing while retaining observed subgroup zeros.

    Parent counts include all cells of that parent, including uncertain fine
    labels and nodes omitted from a particular discovery panel.
    """
    c=np.asarray(node_counts);p=np.asarray(parent_counts);mapping=np.asarray(node_parent)
    if minimum_parent_cells<1:raise ValueError('Positive observed-parent threshold required')
    if c.ndim!=2 or p.ndim!=2 or c.shape[1]!=p.shape[1] or mapping.shape!=(len(c),):raise ValueError('Composition axes differ')
    if mapping.dtype.kind not in 'iu' or (mapping<0).any() or (mapping>=len(p)).any():raise ValueError('Parent mapping invalid')
    if not np.isfinite(c).all() or not np.isfinite(p).all() or (c<0).any() or (p<0).any():raise ValueError('Invalid counts')
    if not np.equal(c,np.floor(c)).all() or not np.equal(p,np.floor(p)).all():raise ValueError('Expected observed counts, not uncalibrated probability counts')
    denominator=p[mapping]
    if (c>denominator).any():raise ValueError('Subset count exceeds complete parent count')
    for parent in np.unique(mapping):
        if (c[mapping==parent].sum(axis=0)>p[parent]).any():raise ValueError('Overlapping or duplicate children exceed parent denominator')
    mask=denominator>=minimum_parent_cells
    result=np.full(c.shape,np.nan,dtype=float)
    np.divide(c,denominator,out=result,where=mask)
    return result,mask

def partition_transport(anchor,node,n_anchor,n_nodes,training_mask):
    """Column-normalized member contingency; no display-string matching.

    Each row is one atom/cell with one anchor and one node. Unknown assignments
    are -1 and excluded explicitly. Training only defines the correspondence.
    """
    a=np.asarray(anchor);b=np.asarray(node);train=np.asarray(training_mask)
    if a.shape!=b.shape or a.ndim!=1 or train.shape!=a.shape or train.dtype.kind!='b':raise ValueError('Partition axes differ')
    if a.dtype.kind not in 'iu' or b.dtype.kind not in 'iu' or (a<-1).any() or (b<-1).any() or (a>=n_anchor).any() or (b>=n_nodes).any():raise ValueError('Partition codes invalid')
    use=train&(a>=0)&(b>=0)
    contingency=coo_matrix((np.ones(use.sum()),(a[use],b[use])),shape=(n_anchor,n_nodes)).toarray()
    total=contingency.sum(axis=0);transport=np.divide(contingency,total[None,:],out=np.zeros_like(contingency),where=total[None,:]>0)
    return transport,total>0,dict(n_training_atoms=int(train.sum()),n_jointly_mapped=int(use.sum()),n_unmapped_training_atoms=int((train&~use).sum()))

def effective_basis(W,S=None):
    W=np.asarray(W,dtype=float)
    if W.ndim!=2 or not np.isfinite(W).all() or (W<0).any():raise ValueError('Invalid nonnegative basis')
    if S is None:return W.copy()
    S=np.asarray(S,dtype=float)
    if S.shape!=(W.shape[1],W.shape[1]) or not np.isfinite(S).all() or (S<0).any():raise ValueError('Invalid nsNMF smoothing matrix')
    return W@S

def map_module_mass(W,transport):
    W=np.asarray(W,float);T=np.asarray(transport,float)
    if W.ndim!=2 or T.ndim!=2 or T.shape[1]!=W.shape[0] or not np.isfinite(W).all() or not np.isfinite(T).all() or (W<0).any() or (T<0).any():raise ValueError('Invalid module transport')
    if np.any((W.sum(axis=1)>0)&~np.isclose(T.sum(axis=0),1)):raise ValueError('Loaded node lacks a complete training correspondence')
    total=W.sum(axis=0);mass=np.divide(W,total[None,:],out=np.zeros_like(W),where=total[None,:]>0)
    return T@mass

def match_modules(A,B,minimum_similarity):
    """One-to-one cosine assignment with explicit unmatched dummy nodes."""
    A=np.asarray(A,float);B=np.asarray(B,float)
    if A.ndim!=2 or B.ndim!=2 or A.shape[0]!=B.shape[0] or not np.isfinite(A).all() or not np.isfinite(B).all() or (A<0).any() or (B<0).any() or not 0<=minimum_similarity<=1:raise ValueError('Invalid aligned module matrices')
    normA=np.linalg.norm(A,axis=0);normB=np.linalg.norm(B,axis=0)
    a=np.divide(A,normA[None,:],out=np.zeros_like(A),where=normA[None,:]>0);b=np.divide(B,normB[None,:],out=np.zeros_like(B),where=normB[None,:]>0);similarity=np.clip(a.T@b,0,1)
    na,nb=similarity.shape;utility=np.zeros((na+nb,na+nb));utility[:na,:nb]=similarity-minimum_similarity;utility[:na,:nb][(normA[:,None]==0)|(normB[None,:]==0)]=-2
    ix,jx=linear_sum_assignment(-utility);pairs=[(int(i),int(j),float(similarity[i,j])) for i,j in zip(ix,jx) if i<na and j<nb and utility[i,j]>0]
    return dict(matches=pairs,unmatched_A=sorted(set(range(na))-{i for i,j,s in pairs}),unmatched_B=sorted(set(range(nb))-{j for i,j,s in pairs}),cosine_matrix=similarity)

def weighted_jaccard(a,b):
    a=np.asarray(a,float);b=np.asarray(b,float)
    if a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all() or (a<0).any() or (b<0).any():raise ValueError('Invalid nonnegative module masses')
    denominator=np.maximum(a,b).sum()
    return float(np.minimum(a,b).sum()/denominator) if denominator else np.nan
