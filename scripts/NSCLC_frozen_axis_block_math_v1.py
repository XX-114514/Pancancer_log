"""CV3 and transformed-residual WCR-B/WCU-B on fixed non-outcome blocks."""
import hashlib
import numpy as np
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design


def seed(axis,stream):
    return int.from_bytes(hashlib.sha256(('2026101933|'+axis+'|'+stream).encode()).digest()[:8],'little')


def geometry(model,blocks):
    D=model['D'];inverse=np.linalg.pinv(D,rcond=1e-12);H=D@inverse;Q=D[:,:-1];H0=Q@Q.T;n=len(D);labels=sorted(set(blocks));G=len(labels)
    if G<50:raise ValueError('FEWER_THAN50_NONOUTCOME_BLOCKS')
    membership=np.asarray([labels.index(x) for x in blocks]);adjust=np.zeros((n,n));adjust0=np.zeros((n,n));A=np.zeros((G,n));maxblock=0.
    for g in range(G):
        ix=np.flatnonzero(membership==g);h=H[np.ix_(ix,ix)];h0=H0[np.ix_(ix,ix)];eigen=np.linalg.eigvalsh(h);maxblock=max(maxblock,float(eigen[-1]));m=np.eye(len(ix))-h;m0=np.eye(len(ix))-h0
        if np.linalg.eigvalsh(m).min()<=1e-8 or np.linalg.eigvalsh(m0).min()<=1e-8:raise ValueError('DELETE_BLOCK_DESIGN_NOT_IDENTIFIABLE')
        a=np.linalg.solve(m,np.eye(len(ix)));a0=np.linalg.solve(m0,np.eye(len(ix)));adjust[np.ix_(ix,ix)]=a;adjust0[np.ix_(ix,ix)]=a0;A[g,ix]=inverse[-1,ix]@a
    return dict(D=D,inverse=inverse,H=H,H0=H0,adjust=adjust,adjust0=adjust0,score=A@(np.eye(n)-H),row=inverse[-1],groups=G,membership=membership,correction=(G-1)/G,max_block_leverage=maxblock,df=n-D.shape[1])


def evaluate(g,y):
    beta=float(g['row']@y);influence=g['score']@y;variance=float(g['correction']*(influence@influence))
    if not np.isfinite(variance) or variance<=1e-24:raise ValueError('DEGENERATE_CV3_VARIANCE')
    se=np.sqrt(variance)
    return dict(effect=beta,se=float(se),t=float(beta/se),leave_block_effect=beta-influence)


def wild(g,y,draws,random_seed,restricted):
    residual=y-(g['H0']@y if restricted else g['H']@y);adjusted=(g['adjust0'] if restricted else g['adjust'])@residual;rng=np.random.default_rng(random_seed);values=np.empty(draws)
    for start in range(0,draws,1000):
        take=min(1000,draws-start);signs=np.where(rng.random((take,g['groups']))<.5,-1.,1.).T;noise=adjusted[:,None]*signs[g['membership']];scores=g['score']@noise;se=np.sqrt(g['correction']*np.sum(scores**2,axis=0));values[start:start+take]=(g['row']@noise)/se
    if not np.isfinite(values).all():raise ValueError('NONFINITE_BLOCK_WILD_STATISTIC')
    return values


def confidence(fit,full_statistics):
    lo,hi=np.quantile(full_statistics,[.025,.975]);return np.array([fit['effect']-hi*fit['se'],fit['effect']-lo*fit['se']])


def fixture():
    # Direct deletion refits independently verify the faster CV3 operators.
    rng=np.random.default_rng(2026101934);maximum=0.;cases=0
    for pattern in [[1]*75+[2,3],[2]*40,[1]*60+[5]*4]:
        blocks=np.repeat(np.arange(len(pattern)).astype(str),pattern);n=len(blocks);x=np.tile([0.,1.],n//2+1)[:n];C=rng.normal(size=(n,3));C[:,2]=0;C[:5,2]=1;y=.4*x+C@np.array([.2,-.3,.1])+rng.normal(size=n)
        model,status=design(x,y,C,['a','b','rare'],'clinical');assert model is not None;g=geometry(model,blocks);f=evaluate(g,y);D=model['D'];beta=np.linalg.lstsq(D,y,rcond=1e-12)[0][-1];deletes=[]
        for label in sorted(set(blocks)):
            keep=blocks!=label;deletes.append(np.linalg.lstsq(D[keep],y[keep],rcond=1e-12)[0][-1])
        se=np.sqrt((len(deletes)-1)/len(deletes)*np.sum((np.asarray(deletes)-beta)**2));np.testing.assert_allclose(f['leave_block_effect'],deletes,atol=1e-9,rtol=1e-9);np.testing.assert_allclose(f['se'],se,atol=1e-10);maximum=max(maximum,float(np.max(np.abs(f['leave_block_effect']-deletes))));cases+=n
    n=60;x=np.tile([0.,1.],30);rare=np.zeros(n);rare[:2]=1;blocks=np.array(['pair']*2+['s'+str(i) for i in range(n-2)]);y=rng.normal(size=n);m,status=design(x,y,rare[:,None],['rare'],'clinical');assert m is not None
    try:geometry(m,blocks)
    except ValueError as error:assert str(error)=='DELETE_BLOCK_DESIGN_NOT_IDENTIFIABLE'
    else:raise AssertionError('SINGULAR_DELETE_BLOCK_MUST_NOT_FALL_BACK')
    return dict(status='PASS_CV3_EQUALS_ALL_DIRECT_BLOCK_DELETIONS_AND_SINGULAR_BLOCK_REJECTED',synthetic_records=cases,maximum_direct_deletion_difference=maximum)
