"""Fixed continuous marker and conditional spatial randomization kernels.

No threshold-defined phenotype, learned signature, clinical model or LR claim.
All arrays are supplied in RAM; functions return RAM arrays or source summaries.
"""
from collections import Counter
import numpy as np
import pandas as pd

AXES_X=['FOXP3_CD4','ICOS_OX40_GITR','Ki67_reference']
AXES_Y=['Ki67','PD1']
CONTROLS=['Histone H3','DNA1','DNA2','AreaShape_Area']


def rankz(values):
    a=np.asarray(values,dtype=float)
    if not np.isfinite(a).all():raise ValueError('nonfinite_marker')
    r=pd.Series(a).rank(method='average').to_numpy();s=r.std()
    return (r-r.mean())/s if s>1e-12 else np.full(len(a),np.nan)


def residualize(score,controls,min_df=8):
    score=np.asarray(score,float);x=np.column_stack([np.ones(len(score)),controls])
    if not np.isfinite(score).all():return None,'CONSTANT_OR_NONFINITE_AXIS'
    if not np.isfinite(x).all():return None,'CONSTANT_OR_NONFINITE_FIXED_CONTROL'
    singular=np.linalg.svd(x,compute_uv=False)
    if x.shape[0]<x.shape[1] or singular[-1]<=singular[0]*1e-10:return None,'FIXED_DESIGN_NOT_FULL_RANK'
    if singular[0]/singular[-1]>1e8:return None,'FIXED_DESIGN_ILL_CONDITIONED'
    if len(score)-x.shape[1]<min_df:return None,'RESIDUAL_DF_BELOW_FIXED_MINIMUM'
    resid=score-x@np.linalg.lstsq(x,score,rcond=None)[0];sd=resid.std()
    if sd<1e-10:return None,'NO_RESIDUAL_MARKER_VARIATION'
    return resid/sd,'ESTIMABLE'


def marker_residuals(frame,parent):
    bg=np.column_stack([rankz(np.log1p(frame[c]) if c=='AreaShape_Area' else frame[c]) for c in CONTROLS])
    cycle=rankz(frame['Ki-67'])
    if parent=='CD8':
        scores=[cycle,rankz(frame['PD-1'])];names=AXES_Y
    else:
        scores=[np.mean([rankz(frame[c]) for c in ['FOXP3','CD4']],axis=0),
                np.mean([rankz(frame[c]) for c in ['ICOS','OX40','GITR']],axis=0),cycle];names=AXES_X
    result=[];reasons=[]
    for name,score in zip(names,scores):
        controls=bg if name in ('Ki67','Ki67_reference') else np.column_stack([bg,cycle])
        resid,reason=residualize(score,controls)
        result.append(np.full(len(frame),np.nan) if resid is None else resid);reasons.append(reason)
    return np.column_stack(result),reasons


def graph_edges(xy_source,xy_target,radius):
    xs=np.asarray(xy_source,float);yt=np.asarray(xy_target,float)
    if not len(xs) or not len(yt):return None,'EMPTY_PARENT'
    # Chunk source coordinates, avoiding a giant all-cell distance matrix.
    left=[];right=[]
    for start in range(0,len(xs),256):
        delta=xs[start:start+256,None,:]-yt[None,:,:]
        i,j=np.nonzero(np.einsum('ijk,ijk->ij',delta,delta)<=radius**2)
        left.extend((i+start).tolist());right.extend(j.tolist())
    i=np.asarray(left,dtype=int);j=np.asarray(right,dtype=int)
    if len(i)<20:return None,'HETEROTYPIC_EDGES_BELOW_20'
    ds=np.bincount(i,minlength=len(xs));dt=np.bincount(j,minlength=len(yt))
    ns=int((ds>0).sum());nt=int((dt>0).sum())
    if min(ns,nt)<10:return None,'CONNECTED_PARENT_CELLS_BELOW_10'
    weights=.5*(1/(ns*ds[i])+1/(nt*dt[j]))
    if not np.isclose(weights.sum(),1):raise ValueError('edge_weight_normalization')
    return (i,j,weights),'ESTIMABLE'


def weighted_corr(source,target,edges):
    i,j,w=edges;s=np.asarray(source)[i];t=np.asarray(target)[j]
    ms=np.einsum('e,ek->k',w,s);mt=np.einsum('e,el->l',w,t)
    covariance=np.einsum('e,ek,el->kl',w,s,t)-np.outer(ms,mt)
    vs=np.einsum('e,ek->k',w,s*s)-ms*ms;vt=np.einsum('e,el->l',w,t*t)-mt*mt
    scale=np.sqrt(np.maximum(vs[:,None],0)*np.maximum(vt[None,:],0))
    return np.divide(covariance,scale,out=np.full_like(covariance,np.nan),where=scale>1e-12)


def permutation_groups(frame,bounds,grid,interfaces=False):
    lo=np.asarray(bounds[0]);hi=np.asarray(bounds[1]);extent=hi-lo
    if not np.isfinite(extent).all() or (extent<=0).any():return [],0
    xy=frame[['Location_Center_X','Location_Center_Y']].to_numpy()
    bins=np.minimum(grid-1,np.maximum(0,np.floor((xy-lo)/extent*grid))).astype(int)
    key=bins[:,0]*grid+bins[:,1]
    if interfaces:
        # -1 is an explicit unknown group, never false.
        a=frame['is_interface'].to_numpy();b=frame['is_perivascular'].to_numpy()
        key=(key*3+np.where(a<0,2,a))*3+np.where(b<0,2,b)
    groups=[np.flatnonzero(key==k) for k in np.unique(key)]
    return groups,sum(len(x) for x in groups if len(x)>1)


def permutations(groups,n,B,rng):
    out=np.broadcast_to(np.arange(n),(B,n)).copy()
    for group in groups:
        if len(group)>1:
            for b in range(B):out[b,group]=rng.permutation(group)
    return out


def permuted_corr(source,target,edges,source_index,target_index):
    i,j,w=edges
    s=np.asarray(source)[source_index[:,i]];t=np.asarray(target)[target_index[:,j]]
    ms=np.einsum('e,bek->bk',w,s);mt=np.einsum('e,bel->bl',w,t)
    covariance=np.einsum('e,bek,bel->bkl',w,s,t)-ms[:,:,None]*mt[:,None,:]
    vs=np.einsum('e,bek->bk',w,s*s)-ms*ms;vt=np.einsum('e,bel->bl',w,t*t)-mt*mt
    scale=np.sqrt(np.maximum(vs[:,:,None],0)*np.maximum(vt[:,None,:],0))
    return np.divide(covariance,scale,out=np.full_like(covariance,np.nan),where=scale>1e-12)


def bootstrap_mean_ci(values,B,rng):
    a=np.asarray(values,float)
    if len(a)<10:return None
    means=np.mean(a[rng.integers(0,len(a),size=(B,len(a)))],axis=1)
    return np.quantile(means,[.025,.975]).tolist()


def paired_roi_contrast(rows_a,rows_b):
    a={row[0]:row[1:] for row in rows_a};b={row[0]:row[1:] for row in rows_b}
    if len(a)!=len(rows_a) or len(b)!=len(rows_b):raise ValueError('duplicate_ROI_slot')
    common=sorted(set(a)&set(b))
    if not common:return None
    return (float(np.mean([a[k][0]-b[k][0] for k in common])),
            np.mean([a[k][1]-b[k][1] for k in common],axis=0))


def cohort_summary(observed,null,case_keys,minimum=20,inference=True,seed=20261022):
    a=np.asarray(observed,float);n=np.asarray(null,float);N=len(a)
    result=dict(n_cases=N,case_unit='author_patient_key',status='COVERAGE_BELOW_PRIMARY_MINIMUM')
    if not N:return result
    if n.ndim!=2 or n.shape[0]!=N or not np.isfinite(a).all() or not np.isfinite(n).all():
        raise ValueError('invalid_case_randomization_arrays')
    delta=a-n.mean(axis=1);effect=float(delta.mean())
    result.update(observed_r=float(a.mean()),conditional_null_mean_r=float(n.mean()),
                  excess_r=effect,CI95=bootstrap_mean_ci(delta,2000,np.random.default_rng(seed)))
    if N>1:
        loo=(delta.sum()-delta)/(N-1)
        result['largest_leave_one_case_change']=float(np.max(np.abs(loo-effect)))
    if not inference:result['status']='DESCRIPTIVE_NO_P';return result
    if N<minimum:return result
    draws=n.mean(axis=0);actual=float(a.mean());B=len(draws)
    p=min(1.,2*min((1+int((draws>=actual).sum()))/(B+1),(1+int((draws<=actual).sum()))/(B+1)))
    result.update(status='ESTIMABLE',conditional_randomization_P=p)
    return result


def bh_fixed_family(results,family_size):
    # Unestimable prespecified slots remain in the denominator.
    slots=[(k,v['conditional_randomization_P']) for k,v in results.items() if 'conditional_randomization_P' in v]
    slots.sort(key=lambda x:x[1]);q=1.
    for rank in range(len(slots),0,-1):
        key,p=slots[rank-1];q=min(q,p*family_size/rank);results[key]['BH_fixed_family']=float(min(q,1))
