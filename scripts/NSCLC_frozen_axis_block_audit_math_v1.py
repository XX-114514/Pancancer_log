"""Independent explicit deletion QR operators and block-bootstrap replay."""
import numpy as np
from NSCLC_immune_RNA_MCP_cohort_audit_math_v1 import prepare


def geometry(model,blocks):
    D=model['D'];Q,R=np.linalg.qr(D,mode='reduced');inverse=np.linalg.solve(R,Q.T);N=D[:,:-1];q0,r0=np.linalg.qr(N,mode='reduced');labels=sorted(set(blocks));G=len(labels);n=len(D)
    if G<50:raise ValueError('FEWER_THAN50_NONOUTCOME_BLOCKS')
    row=inverse[-1];deleted=np.zeros((G,n));adjust=np.zeros((n,n));adjust0=np.zeros((n,n));members=np.asarray([labels.index(v) for v in blocks]);maximum=0.
    for g,label in enumerate(labels):
        ix=np.flatnonzero(members==g);keep=members!=g;h=Q[ix]@Q[ix].T;h0=q0[ix]@q0[ix].T;maximum=max(maximum,float(np.linalg.eigvalsh(h)[-1]))
        if np.linalg.eigvalsh(np.eye(len(ix))-h).min()<=1e-8 or np.linalg.eigvalsh(np.eye(len(ix))-h0).min()<=1e-8:raise ValueError('DELETE_BLOCK_DESIGN_NOT_IDENTIFIABLE')
        q,r=np.linalg.qr(D[keep],mode='reduced');p=np.linalg.solve(r,q.T);deleted[g,keep]=p[-1];adjust[np.ix_(ix,ix)]=np.eye(len(ix));adjust[np.ix_(ix,np.flatnonzero(keep))]=-D[ix]@p
        q,r=np.linalg.qr(N[keep],mode='reduced');p=np.linalg.solve(r,q.T);adjust0[np.ix_(ix,ix)]=np.eye(len(ix));adjust0[np.ix_(ix,np.flatnonzero(keep))]=-N[ix]@p
    return dict(D=D,row=row,deleted=deleted,delta_rows=deleted-row[None],adjust=adjust,adjust0=adjust0,groups=G,membership=members,correction=(G-1)/G,max_block_leverage=maximum,df=n-D.shape[1])


def evaluate(g,y):
    beta=float(g['row']@y);deleted=g['deleted']@y;se=np.sqrt(g['correction']*np.sum((deleted-beta)**2))
    if not np.isfinite(se) or se<=1e-12:raise ValueError('DEGENERATE_CV3_VARIANCE')
    return dict(effect=beta,se=float(se),t=float(beta/se),leave_block_effect=deleted)


def wild(g,y,draws,random_seed,restricted):
    adjusted=(g['adjust0'] if restricted else g['adjust'])@y;rng=np.random.default_rng(random_seed);out=np.empty(draws)
    for start in range(0,draws,733):
        size=min(733,draws-start);signs=np.where(rng.random((size,g['groups']))<.5,-1.,1.).T;noise=adjusted[:,None]*signs[g['membership']];beta=g['row']@noise;deltas=g['delta_rows']@noise;se=np.sqrt(g['correction']*np.sum(deltas**2,axis=0));out[start:start+size]=beta/se
    assert np.isfinite(out).all();return out


def prepare_axis(a,j,mode='main'):
    s=int(a['sender_index'][j]);r=int(a['receiver_index'][j]);C=np.column_stack([a['clinical_controls'],a['global_quality'],np.log1p(a['unit_cells'][:,s]),np.log1p(a['unit_cells'][:,r]),a['unit_quality'][:,s],a['unit_quality'][:,r],a['within_major_fraction'][:,s],a['within_major_fraction'][:,r]])
    names=a['clinical_control_names'].tolist()+['logN_all_immune','all_immune_mean_logRNA','all_immune_mean_loggenes','logN_sender','logN_receiver','sender_meanlogRNA','sender_meanloggenes','receiver_meanlogRNA','receiver_meanloggenes','sender_within_major_fraction','receiver_within_major_fraction'];y=a['standardized_axis_score'][:,j];valid=a['endpoint_measurable'][:,j]&np.isfinite(np.c_[a['response'],y,C]).all(1);x=a['response'][valid];y=y[valid];C=C[valid]
    if mode=='unadjusted':C=np.empty((len(x),0));names=[]
    elif mode=='clinical_RNA_counts_no_fractions':C=C[:,:-2];names=names[:-2]
    elif mode=='main_plus_total_communication':C=np.column_stack([C,a['log1p_mean_observed_normalized_LRscore'][valid],a['captured_author_types'][valid]]);names+=['total_observed_LR_background','captured_author_types']
    elif mode!='main':raise ValueError('UNDECLARED_SENSITIVITY')
    if min(int((x==0).sum()),int((x==1).sum()))<10:return None,None,'FEWER_THAN10_PER_RESPONSE',valid
    m,status=prepare(x,y,C,names,'clinical')
    if m is None:return None,None,status,valid
    if m['df']<20:return None,None,'FEWER_THAN20_RESIDUAL_DF',valid
    try:g=geometry(m,a['profile_blocks'][valid])
    except ValueError as e:return m,None,str(e),valid
    return m,g,'ESTIMABLE',valid


def synthetic_stream_gate():
    from NSCLC_immune_RNA_MCP_cohort_math_v1 import design
    from NSCLC_frozen_axis_block_math_v2 import geometry as production_geometry,evaluate as production_evaluate,wild as production_wild
    rng=np.random.default_rng(2026101935);maximum=0.
    for sizes in [[1]*75+[2,3],[2]*60,[1]*60+[5]*4]:
        blocks=np.repeat(np.arange(len(sizes)).astype(str),sizes);n=len(blocks);x=np.tile([0.,1.],n//2+1)[:n];C=rng.normal(size=(n,3));C[:,2]=0;C[:5,2]=1;y=.2*x+C@np.array([.4,-.2,.3])+rng.normal(size=n)
        p,_=design(x,y,C,['a','b','rare'],'clinical');a,_=prepare(x,y,C,['a','b','rare'],'clinical');g=production_geometry(p,blocks);h=geometry(a,blocks);f=production_evaluate(g,y);v=evaluate(h,y)
        np.testing.assert_allclose([f['effect'],f['se'],f['t']],[v['effect'],v['se'],v['t']],atol=1e-9,rtol=1e-9);np.testing.assert_allclose(f['leave_block_effect'],v['leave_block_effect'],atol=1e-9,rtol=1e-9)
        for restricted in [True,False]:
            left=production_wild(g,y,257,2026101936,restricted);right=wild(h,y,257,2026101936,restricted);np.testing.assert_allclose(left,right,atol=1e-9,rtol=1e-9);maximum=max(maximum,float(np.max(np.abs(left-right))))
    return dict(status='PASS_INDEPENDENT_DELETED_QR_ALL_SYNTHETIC_CV3_WCRB_WCUB_STREAMS',maximum_statistic_difference=maximum)
