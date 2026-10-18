"""Independent full exact/permutation/bootstrap/FDR and all descriptive replay."""
import hashlib
import itertools
import json
import os
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binom,beta
from statsmodels.stats.multitest import multipletests
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import wait


def close(a,b,tol=1e-10):assert np.allclose(a,b,rtol=0,atol=tol,equal_nan=True)


def seed(s):return int.from_bytes(hashlib.sha256(('2026101906|'+s).encode()).digest()[:8],'little')


def rho(x,y):
    if min(np.ptp(x),np.ptp(y))==0:return np.nan
    a=pd.Series(x).rank().to_numpy();b=pd.Series(y).rank().to_numpy();a-=a.mean();b-=b.mean();return float(np.sum(a*b)/np.sqrt(np.sum(a*a)*np.sum(b*b)))


def effect(x,y,co):return rho(x,y) if co else float(pd.Series(y).median())


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);src=root/'GSE176021_paired_tissue_inputs_20261019_v1';producer=root/'GSE176021_paired_tissue_associations_20261019_v1';out=root/'GSE176021_paired_tissue_association_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700);term=wait(root,out,'test_GSE176021_paired_tissue_family_v1',producer/'private_terminal.json',{'COMPLETE_ALL6_TISSUE_TESTS54_SENSITIVITIES30_CONTEXT_POINTS_AND_FULL_STREAMS_AUDIT_PENDING'})
    for n,v in term['files'].items():assert sha(producer/n)==v['sha256']
    assert sha(root/'GSE176021_paired_tissue_input_audit_20261019_v1/private_terminal.json')=='1998d626c57e0165ab4b66573170d16a904c600e442823b86dfb37dd39ea5f89';table=pd.read_csv(producer/'private_all6_primary_paired_tissue_results.tsv',sep='\t');registry=pd.read_csv(src/'private_all6_registered_primary_tests.tsv',sep='\t');assert len(table)==6 and table.model_id.tolist()==registry.model_id.tolist();index=table.set_index('model_id');ix=np.random.default_rng(seed('all_primary_paired_cases')).integers(0,10,size=(2000,10));exact_count=0;random_count=0;maxerr=0.
    for r in registry.itertuples(index=False):
        a=np.load(src/('private_'+r.model_id+'_inputs.npz'));s=np.load(producer/('private_'+r.model_id+'_complete_stream.npz'));assert np.array_equal(s['cases'],a['cases']) and np.array_equal(s['bootstrap_indices'],ix);x=a['x'];y=a['y'];co=r.outcome=='RE_tissue_cochange';e=effect(x,y,co);row=index.loc[r.model_id];close(row.effect,e);assert row.n_pairs==10 and row.program==r.program and row.outcome==r.outcome
        if co and not np.isfinite(e):assert row.status=='CONSTANT_COCHANGE_TARGET' and pd.isna(row.P);b=np.full(2000,np.nan);loo=np.full(10,np.nan)
        else:
            if co:
                xr=pd.Series(x).rank().to_numpy();yr=pd.Series(y).rank().to_numpy();xr-=xr.mean();yr-=yr.mean();den=np.sqrt(np.sum(xr*xr)*np.sum(yr*yr));rng=np.random.default_rng(seed(r.model_id+'|cochange'));vv=[]
                for lo in range(0,19999,512):
                    ids=np.argsort(rng.random((min(512,19999-lo),10)),axis=1);vv.extend((np.sum(yr[ids]*xr,axis=1)/den).tolist())
                vals=np.array(vv);close(s['null_statistics'],vals,1e-12);maxerr=max(maxerr,float(np.max(np.abs(s['null_statistics']-vals))));k=int((np.abs(vals)>=abs(e)-1e-12).sum());P=(k+1)/20000;u=1. if k==19999 else float(beta.ppf(1-.01/2,k+1,19999-k));u=min(1.,(1+19999*u)/20000);close([row.P,row.MC_upper_P,row.statistic],[P,u,e],1e-12);assert row.exceedances==k;random_count+=len(vals)
            else:
                d=y[y!=0];ranks=pd.Series(np.abs(d)).rank().to_numpy();observed=float(np.sum(np.sign(d)*ranks));vals=np.array([np.sum(np.array(sign)*ranks) for sign in itertools.product([-1,1],repeat=len(d))]);P=float((np.abs(vals)>=abs(observed)-1e-12).mean());close(np.sort(s['null_statistics']),np.sort(vals));close([row.P,row.MC_upper_P,row.statistic],[P,P,observed],1e-12);assert row.exact and row.nonzero_differences==len(d);positive=int((d>0).sum());sp=min(1.,2*binom.cdf(min(positive,len(d)-positive),len(d),.5)) if len(d) else 1.;close(row.sign_P,sp,1e-14);exact_count+=len(vals)
            b=np.array([effect(x[j],y[j],co) for j in ix]);loo=np.array([effect(np.delete(x,j),np.delete(y,j),co) for j in range(10)])
        close(s['bootstrap_effects'],b);close(s['leave_one_effects'],loo);good=np.isfinite(b);ci=np.quantile(b[good],[.025,.975]) if good.sum()>=1900 else [np.nan,np.nan];close([row.CI_low,row.CI_high],ci);assert row.valid_bootstrap==good.sum();close([row.LOO_min,row.LOO_max],[np.nanmin(loo),np.nanmax(loo)] if np.isfinite(loo).any() else [np.nan,np.nan]);assert row.LOO_same_sign==bool(np.isfinite(loo).all() and (loo*e>0).all())
    for col,P,method in [('BH',table.P.fillna(1),'fdr_bh'),('BY',table.P.fillna(1),'fdr_by'),('MC_upper_BH',table.MC_upper_P.fillna(1),'fdr_bh')]:
        expected=multipletests(P,method=method)[1];mask=table.P.notna().to_numpy();close(table.loc[mask,col],expected[mask]);assert table.loc[~mask,col].isna().all()
    mask=table.sign_P.notna();close(table.loc[mask,'sign_BH'],multipletests(table.loc[mask,'sign_P'],method='fdr_bh')[1]);assert table.loc[~mask,'sign_BH'].isna().all();ci=table.CI_low*table.CI_high>0;assert np.array_equal(ci,table.CI_excludes_zero);priority=table.BH.lt(.05)&table.MC_upper_BH.lt(.05)&ci&table.LOO_same_sign;assert np.array_equal(priority,table.priority_association)
    sens=pd.read_csv(producer/'private_all54_point_sensitivity_results.tsv',sep='\t');raw=pd.read_csv(src/'private_all_prespecified_case_sensitivity_inputs.tsv',sep='\t');assert len(sens)==54
    for r in sens.itertuples(index=False):
        d=raw[raw.model_id.eq(r.model_id)&raw['mode'].eq(r.mode)];good=np.isfinite(d.x)&np.isfinite(d.y);co=index.loc[r.model_id,'outcome']=='RE_tissue_cochange';expected=effect(d.x.to_numpy(),d.y.to_numpy(),co) if good.all() else np.nan;close(r.effect,expected);assert r.usable_case_rows==good.sum() and r.registered_case_rows==len(d) and r.primary_N==10 and r.no_new_P and r.status==('COMPLETE_POINT_ONLY_NO_P' if good.all() else 'INCOMPLETE_PRIMARY_CASE_UNIVERSE_NO_POINT')
    follow=pd.read_csv(producer/'private_all30_coactivity_context_points_only.tsv',sep='\t');cases=pd.read_csv(src/'private_all_scope_case_tissue_values.tsv',sep='\t');cases=cases[cases['mode'].eq('primary')];assert len(follow)==30
    for r in follow.itertuples(index=False):
        use=cases[cases.scope.eq(r.scope)&cases.tissue.eq(r.tissue)&cases.program.eq(r.program)];wide=use.pivot(index='case_block',columns='view',values='value');wide=wide.reindex(columns=['R_FOXP3_capture','E_effector_capture']).dropna();x=wide.iloc[:,0].to_numpy();y=wide.iloc[:,1].to_numpy();n=len(wide);assert n==r.case_units and r.no_P;close(r.rho,rho(x,y) if n>=3 else np.nan);loo=np.array([rho(np.delete(x,j),np.delete(y,j)) for j in range(n)]) if n>=4 else np.array([]);close([r.LOO_min,r.LOO_max],[np.nanmin(loo),np.nanmax(loo)] if np.isfinite(loo).any() else [np.nan,np.nan]);assert r.status==('POINT_ONLY_NO_NEW_INFERENCE' if n>=10 else 'BELOW10_DESCRIPTIVE_ONLY' if n>=3 else 'INSUFFICIENT_CASES_NO_POINT')
    assert exact_count==term['exact_sign_statistics'] and random_count==term['permutation_statistics'];d={'status':'ACCEPTED_ALL6_PRIMARY54_SENSITIVITIES30_CONTEXTS_AND_COMPLETE_RANDOM_STREAMS','completed_UTC':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(__file__),'producer_terminal_sha256':sha(producer/'private_terminal.json'),'exact_sign_statistics':exact_count,'permutation_statistics':random_count,'bootstrap_values':12000,'leave_one_values':60,'maximum_random_statistic_error':maxerr,'no_new_MPR_or_old_source_P':True}
    (out/'private_terminal.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'status':d['status'],'terminal_sha256':sha(out/'private_terminal.json')}))


if __name__=='__main__':main()
