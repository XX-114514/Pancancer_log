"""Exact descriptive partition on a fixed NNLS active face; not causal effects."""
import numpy as np


def projection_parts(Y,B,U,V,available,coefficient):
    k=B.shape[1];L,S,_=Y.shape;active=np.flatnonzero(coefficient>0);lr=np.zeros((L,k));pair=np.zeros((S,S,k))
    if not len(active):return lr,pair
    U=U*available[:,None];V=V*available[:,None];G=(B.T@B)*(U.T@U)*(V.T@V);G=G[np.ix_(active,active)]
    lr_rhs=np.einsum('lsr,lk,sk,rk->lk',Y,B[:,active],U[:,active],V[:,active],optimize=True)
    pair_rhs=np.einsum('lsr,lk,sk,rk->srk',Y,B[:,active],U[:,active],V[:,active],optimize=True)
    lr[:,active]=np.linalg.solve(G,lr_rhs.T).T;pair[:,:,active]=np.linalg.solve(G,pair_rhs.reshape(-1,len(active)).T).T.reshape(S,S,len(active))
    np.testing.assert_allclose(lr.sum(0),coefficient,atol=1e-7,rtol=1e-7);np.testing.assert_allclose(pair.sum((0,1)),coefficient,atol=1e-7,rtol=1e-7)
    return lr,pair


def log_activity_multiplier(coefficient,median,sd):
    out=1/(median*sd);return np.divide(np.log1p(coefficient/median),coefficient*sd,out=out.copy(),where=coefficient>0)
