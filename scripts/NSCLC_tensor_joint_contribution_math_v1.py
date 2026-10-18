"""Two independent operators for descriptive joint partitions, not new inference."""
import numpy as np
from scipy.optimize import nnls


def gram_joint(Y,B,U,V,available,coefficient):
    active=np.flatnonzero(coefficient>0);k=B.shape[1];L,S,_=Y.shape
    out=np.zeros((k,L,S,S))
    if not len(active):return out
    U=U*available[:,None];V=V*available[:,None]
    G=(B.T@B)*(U.T@U)*(V.T@V);G=G[np.ix_(active,active)]
    inverse=np.linalg.solve(G,np.eye(len(active)))
    template=np.einsum('la,sa,ra->lsra',B[:,active],U[:,active],V[:,active])
    out[active]=np.einsum('lsra,ak->klsr',template,inverse,optimize=True)*Y[None]
    np.testing.assert_allclose(out.sum((1,2,3)),coefficient,atol=1e-7,rtol=1e-7)
    return out


def qr_joint(Y,B,U,V,available,coefficient):
    active=np.flatnonzero(coefficient>0);k=B.shape[1];L,S,_=Y.shape
    out=np.zeros((k,L,S,S));ix=np.flatnonzero(available)
    if not len(active):return out
    matrix=np.einsum('la,sa,ra->lsra',B[:,active],U[available][:,active],V[available][:,active]).reshape(-1,len(active))
    Q,R=np.linalg.qr(matrix,mode='reduced');operator=np.linalg.solve(R,Q.T)
    values=operator.reshape(len(active),L,len(ix),len(ix))*Y[:,available][:,:,available][None]
    for s,source in enumerate(ix):
        for r,target in enumerate(ix):out[np.ix_(active,np.arange(L),[source],[target])]=values[:,:,s:s+1,r:r+1]
    np.testing.assert_allclose(out.sum((1,2,3)),coefficient,atol=1e-7,rtol=1e-7)
    return out


def fixture():
    rng=np.random.default_rng(2026101931);L,S,K=7,5,3;B=rng.uniform(.1,1,(L,K));B[-1]=0;U=rng.uniform(.1,1,(S,K));V=rng.uniform(.1,1,(S,K));maximum=0.
    for i in range(12):
        a=np.ones(S,bool);a[i%S]=False;wanted=rng.uniform(.1,1,K)
        if i==0:wanted[:]=0
        if i==1:wanted[-1]=0
        D=np.einsum('lk,sk,rk->lsrk',B,U,V);Y=np.einsum('lsrk,k->lsr',D,wanted);Y[:,~a,:]=0;Y[:,:,~a]=0
        c,_=nnls(D[:,a][:,:,a].reshape(-1,K),Y[:,a][:,:,a].ravel())
        g=gram_joint(Y,B,U,V,a,c);q=qr_joint(Y,B,U,V,a,c)
        np.testing.assert_allclose(g,q,atol=1e-9,rtol=1e-9);maximum=max(maximum,float(np.max(np.abs(g-q))))
        assert not np.any(g[:,:,~a,:]) and not np.any(g[:,:,:,~a])
    return dict(status='PASS_TRUE_ZERO_MISSING_ACTIVE_BOUNDARY_ALL_JOINT_POSITIONS_GRAM_VS_QR',maximum_difference=maximum)
