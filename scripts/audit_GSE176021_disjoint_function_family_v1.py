"""Independently replay all 400 functional tests, random arrays and sensitivities."""
import hashlib
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import beta
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import wait


def ranks(x):
    x=np.asarray(x);return 2*(x[:,None]>x[None,:]).sum(axis=1)+(x[:,None]==x[None,:]).sum(axis=1)-len(x)


def rho(x,y):
    a,b=ranks(x),ranks(y);den=float(np.sqrt(sum(a*a)*sum(b*b)));return float(sum(a*b)/den) if den else np.nan


def correct(p,by=False):
    p=np.asarray(p,float);unique=sorted(set(p));h=sum(1/k for k in range(1,len(p)+1)) if by else 1.
    return np.array([min(1.,min(len(p)*h*v/np.count_nonzero(p<=v) for v in unique if v>=x)) for x in p])


def close(a,b):np.testing.assert_allclose(a,b,atol=1e-10,rtol=0,equal_nan=True)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);out=root/'GSE176021_disjoint_function_association_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    src=root/'GSE176021_disjoint_function_associations_20261019_v1';term=wait(root,out,'test_GSE176021_disjoint_function_family_v1',src/'private_terminal.json',{'COMPLETE_ALL400_FUNCTIONAL_ASSOCIATIONS3200_SENSITIVITIES_AND_FULL_STREAMS_AUDIT_PENDING'})
    for name,v in term['files'].items():assert sha(src/name)==v['sha256']
    inp=root/'GSE176021_disjoint_function_inputs_20261019_v1';a=dict(np.load(inp/'private_all400_case_inputs.npz'));reg=pd.read_csv(inp/'private_all400_registered_functional_questions.tsv',sep='\t');result=pd.read_csv(src/'private_all400_primary_functional_results.tsv',sep='\t');sens=pd.read_csv(src/'private_all3200_point_sensitivity_results.tsv',sep='\t');common=dict(np.load(src/'private_common_case_and_bootstrap_arrays.npz'));assert len(result)==400 and len(sens)==3200
    pd.testing.assert_frame_equal(result[reg.columns],reg);assert np.array_equal(common['cases'],a['cases'])
    bsseed=int.from_bytes(hashlib.sha256(b'2026101908|all10_shared_case_pairs').digest()[:8],'little');ix=np.random.default_rng(bsseed).integers(0,10,size=(2000,10));assert np.array_equal(common['bootstrap_indices'],ix)
    ps=[];us=[];cis=[];same=[];random_count=0;max_null_error=0;max_bs_error=0.
    for i,r in enumerate(reg.itertuples(index=False)):
        x,y=a['x'][i],a['y'][i];ar,br=ranks(x),ranks(y);numerator=int(sum(ar*br));den=float(np.sqrt(sum(ar*ar)*sum(br*br)));effect=numerator/den if den else np.nan
        row=result.iloc[i];stream=dict(np.load(src/('private_'+r.model_id+'_complete_stream.npz')));close([row.effect,row.observed_numerator,row.rank_denominator],[effect,numerator,den])
        if den:
            s=int.from_bytes(hashlib.sha256(('2026101908|'+r.model_id+'|full_case_permutation').encode()).digest()[:8],'little');rng=np.random.default_rng(s);expected=[]
            for lo in range(0,99999,717):
                order=np.argsort(rng.random((min(717,99999-lo),10)),axis=1);expected.extend(np.sum(br[order]*ar[None,:],axis=1).tolist())
            expected=np.array(expected,np.int64);assert len(stream['permutation_numerators'])==99999;error=int(np.max(np.abs(expected-stream['permutation_numerators'])));max_null_error=max(max_null_error,error);assert error==0
            k=int(np.count_nonzero(np.abs(expected)>=abs(numerator)));p=(k+1)/100000;upper=1. if k==99999 else float(beta.ppf(1-.01/400,k+1,99999-k));u=min(1.,(1+99999*upper)/100000);close([row.P,row.MC_upper_P,row.exceedances],[p,u,k]);assert row.status=='ESTIMABLE';random_count+=99999
            bx,by=x[ix],y[ix]
            rx=2*(bx[:,:,None]>bx[:,None,:]).sum(axis=2)+(bx[:,:,None]==bx[:,None,:]).sum(axis=2)-10
            ry=2*(by[:,:,None]>by[:,None,:]).sum(axis=2)+(by[:,:,None]==by[:,None,:]).sum(axis=2)-10
            bd=np.sqrt(np.sum(rx*rx,axis=1)*np.sum(ry*ry,axis=1));bs=np.divide(np.sum(rx*ry,axis=1),bd,out=np.full(2000,np.nan),where=bd>0);loo=np.array([rho(np.delete(x,j),np.delete(y,j)) for j in range(10)])
        else:
            p,u=np.nan,np.nan;bs=np.full(2000,np.nan);loo=np.full(10,np.nan);assert len(stream['permutation_numerators'])==0 and row.status=='CONSTANT_INPUT_NO_PRIMARY_ASSOCIATION'
        close(stream['bootstrap_effects'],bs);close(stream['leave_one_effects'],loo)
        if np.isfinite(bs).any():max_bs_error=max(max_bs_error,float(np.nanmax(np.abs(stream['bootstrap_effects']-bs))))
        valid=int(np.isfinite(bs).sum());ci=np.percentile(bs[np.isfinite(bs)],[2.5,97.5]) if valid>=1900 else [np.nan,np.nan];loomin=float(np.min(loo[np.isfinite(loo)])) if np.isfinite(loo).any() else np.nan;loomax=float(np.max(loo[np.isfinite(loo)])) if np.isfinite(loo).any() else np.nan;st=bool(np.isfinite(loo).all() and (loo*effect>0).all())
        close([row.CI_low,row.CI_high,row.valid_bootstrap,row.LOO_min,row.LOO_max],[*ci,valid,loomin,loomax]);assert row.LOO_same_sign==st;ps.append(p);us.append(u);cis.append(bool(ci[0]*ci[1]>0));same.append(st)
        for j,mode in enumerate(a['modes']):
            sr=sens[sens.model_id.eq(r.model_id)&sens['mode'].eq(mode)];assert len(sr)==1;sr=sr.iloc[0]
            if mode=='response_quality_partial_two':
                D=np.column_stack([np.ones(10),ranks(a['response_quality'][i,:,0]),ranks(a['response_quality'][i,:,1])])
                if np.linalg.matrix_rank(D)<3:ee=np.nan;status='COLLINEAR_RESPONSE_QUALITY_NO_COVARIATE_DROPPED'
                else:
                    Q,_=np.linalg.qr(D);xx=ar.astype(float);yy=br.astype(float);xx-=Q@(Q.T@xx);yy-=Q@(Q.T@yy)
                    if min(np.linalg.norm(xx),np.linalg.norm(yy))<1e-10:ee=np.nan;status='ZERO_RESIDUAL_VARIANCE_NO_POINT'
                    else:ee=float(xx@yy/np.sqrt((xx@xx)*(yy@yy)));status='POINT_ONLY_TWO_RESPONSE_QUALITY_CONTROLS'
            else:
                ee=rho(a['sensitivity_x'][i,j],a['sensitivity_y'][i,j]);status='COMPLETE_POINT_ONLY_NO_P' if np.isfinite(ee) else 'CONSTANT_SENSITIVITY_NO_POINT'
            close(sr.effect,ee);assert sr.status==status and sr.case_units==10 and sr.no_new_P
        if (i+1)%25==0:(out/'private_progress.json').write_text(json.dumps({'phase':'FULL400_COMPLETE_RANDOM_AND_SENSITIVITY_REPLAY','positions_replayed':i+1,'positions':400,'observed_UTC':datetime.now(timezone.utc).isoformat()})+'\n')
    good=np.isfinite(ps);bh=correct(np.where(good,ps,1));by=correct(np.where(good,ps,1),True);mc=correct(np.where(good,us,1));expected=(bh<.05)&(mc<.05)&np.array(cis)&np.array(same);bh[~good]=np.nan;by[~good]=np.nan;mc[~good]=np.nan
    close(result.BH,bh);close(result.BY,by);close(result.MC_upper_BH,mc);assert np.array_equal(result.CI_excludes_zero,cis) and np.array_equal(result.priority_functional_association,expected)
    assert random_count==term['permutation_statistics']
    receipt={'status':'ACCEPTED_ALL400_FUNCTIONAL_RESULTS3200_SENSITIVITIES_COMPLETE_RANDOM_STREAM_REPLAY','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(src/'private_terminal.json'),'resource_guard':guard,'primary_positions':400,'point_sensitivities':3200,'permutation_statistics':random_count,'bootstrap_positions':800000,'leave_one_positions':4000,'maximum_integer_null_error':max_null_error,'maximum_bootstrap_error':max_bs_error,'private_scientific_results_not_exported':True}
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':receipt['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
