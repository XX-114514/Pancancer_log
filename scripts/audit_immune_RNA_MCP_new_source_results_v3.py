"""Every main/sensitivity estimate and complete random stream independently replayed."""
import json
import os
import hashlib
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import beta as beta_distribution
from statsmodels.stats.multitest import multipletests
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from immune_RNA_MCP_new_source_io_v1 import read,wait
from immune_RNA_MCP_new_source_math_v1 import independent_design,evaluate,independent_wild,bootstrap,leave_one,seed


def close(a,b,tol=1e-8):assert np.allclose(np.asarray(a,float),np.asarray(b,float),rtol=0,atol=tol,equal_nan=True)


def rankref(x,y,binary=False):
    if min(np.ptp(x),np.ptp(y))==0:return np.nan
    if binary:
        one=y[x==1];zero=y[x==0];u=(one[:,None]>zero).sum()+.5*(one[:,None]==zero).sum();return float(2*u/(len(one)*len(zero))-1)
    a=pd.Series(x).rank().to_numpy();b=pd.Series(y).rank().to_numpy();a-=a.mean();b-=b.mean();return float(np.sum(a*b)/np.sqrt(np.sum(a*a)*np.sum(b*b)))


def permutations(x,y,slot,binary=False):
    yr=pd.Series(y).rank().to_numpy();xr=pd.Series(x).rank().to_numpy();xr-=xr.mean();center=yr-yr.mean();rng=np.random.default_rng(seed(slot,'rank_permutation'));values=[]
    for lo in range(0,19999,1000):
        idx=np.argsort(rng.random((min(1000,19999-lo),len(x))),axis=1)
        if binary:
            v=yr[idx];first=v[:,x==1];second=v[:,x==0];u=(first[:,:,None]>second[:,None,:]).sum(axis=(1,2))+.5*(first[:,:,None]==second[:,None,:]).sum(axis=(1,2));values.extend((2*u/(first.shape[1]*second.shape[1])-1).tolist())
        else:values.extend((np.einsum('ij,j->i',center[idx],xr)/np.sqrt(np.sum(center**2)*np.sum(xr**2))).tolist())
    return np.asarray(values)


def indices(row,a):
    key=('clinical|'+row.endpoint if row.family.startswith('clinical') else 'source|'+row.context_id)+'|'+row.pair+'|'+hashlib.sha256(a['cases'].tobytes()).hexdigest();rng=np.random.default_rng(seed(key,'case_bootstrap'));pieces=[]
    for level in sorted(set(a['strata'])):
        ids=np.flatnonzero(a['strata']==level);pieces.append(ids[rng.integers(0,len(ids),size=(2000,len(ids)))])
    return np.column_stack(pieces).astype(np.int16),key


def check_summary(row,b,loo,effect):
    valid=np.isfinite(b);finite=np.isfinite(loo);ci=np.quantile(b[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan];close([row.CI_low,row.CI_high],ci);assert row.valid_bootstrap==valid.sum() and row.LOO_valid==finite.sum();close([row.LOO_min,row.LOO_max],[loo[finite].min(),loo[finite].max()] if finite.any() else [np.nan,np.nan]);assert row.LOO_same_sign==bool(finite.all() and (loo*effect>0).all())


def check_p(row,values,observed,streams):
    count=int((np.abs(values)>=abs(observed)-1e-12).sum());p=(1+count)/20000;u=1. if count==19999 else beta_distribution.ppf(1-.01/streams,count+1,19999-count);u=min(1.,(1+19999*u)/20000);close([row.P,row.MC_upper_P],[p,u],1e-12);assert row.wild_exceedances==count


def ppoint(x,y,C,names,kind,marginal=False):
    if marginal:
        if len(x)<3 or min(np.ptp(x),np.ptp(y))==0:return np.nan,'INSUFFICIENT_VARIATION'
        return rankref(x,y),'DESCRIPTIVE_RANK_POINT'
    m,status=independent_design(x,y,C,names,kind)
    if m is None:return np.nan,status
    if kind=='clinical':e=float(np.linalg.solve(m['r'],m['q'].T@m['y'])[-1])
    else:
        xx=m['x']-m['Q']@(m['Q'].T@m['x']);yy=m['y']-m['Q']@(m['Q'].T@m['y']);e=float(np.sum(xx*yy)/np.sqrt(np.sum(xx**2)*np.sum(yy**2)))
    return e,'ESTIMABLE_POINT'


def check_sensitivities(r,a,table):
    clinical=r.family.startswith('clinical');marginal=r.family=='source_marginal';x=a['x'];y=a['y'];C=a['C'];names=a['names']
    for mode in ['unadjusted_metric_matched','all_captured_own_quality','fixed15_stress_complement','RNA_quality_only','abundance_only','plus_view_lineage_mix','component_A','component_B']:
        e=np.nan;st='NOT_APPLICABLE'
        if not r.metadata_preeligible:st='INPUT_NOT_METADATA_PREELIGIBLE'
        elif mode=='unadjusted_metric_matched':
            if clinical and r.kind=='clinical':e=float(y[x==1].mean()-y[x==0].mean());st='DESCRIPTIVE_GROUP_MEAN_DIFFERENCE'
            else:e,st=ppoint(x,y,np.empty((len(x),0)),np.array([],str),r.kind,True)
        elif mode=='plus_view_lineage_mix' and clinical:e,st=ppoint(x,y,a['fineC'],a['fine_names'],r.kind)
        elif mode in ['component_A','component_B'] and clinical:
            if r.kind=='clinical':
                m,st=independent_design(x,y,C,names,r.kind)
                if m is not None:e=float(np.linalg.solve(m['r'],m['q'].T@a[mode])[-1]);st='COMPONENT_OLS_POINT_NOT_INDEPENDENT_CLINICAL_TEST'
            else:e,st=ppoint(x,a[mode],C,names,r.kind)
        elif mode=='all_captured_own_quality':e,st=ppoint(a['allcapture_x'],a['allcapture_y'],a['allcapture_C'],names,r.kind,marginal)
        elif mode=='fixed15_stress_complement':e,st=ppoint(x if clinical else a['stress_x'],a['stress_y'],C,names,r.kind,marginal)
        elif mode in ['RNA_quality_only','abundance_only']:
            keep=[i for i,n in enumerate(names) if not n.startswith('quality::') or n.split('::')[-1] in (['cellQ','geneQ'] if mode=='RNA_quality_only' else ['fraction','logN'])];e,st=ppoint(x,y,C[:,keep],names[keep],r.kind,marginal)
        z=table.loc[(r.model_id,mode)];assert z.status==st and z.n_units==len(x) and z.no_new_P;close(z.effect,e)


def correction(frame,separate=True):
    groups=frame.groupby('family').groups if separate else {'all':frame.index}
    for _,group in groups.items():
        use=[i for i in group if frame.loc[i,'metadata_preeligible']]
        if not use:continue
        p=frame.loc[use,'P'].fillna(1).to_numpy(float);u=frame.loc[use,'MC_upper_P'].fillna(1).to_numpy(float);bh=multipletests(p,method='fdr_bh')[1];by=multipletests(p,method='fdr_by')[1];mc=multipletests(u,method='fdr_bh')[1];present=frame.loc[use,'P'].notna().to_numpy()
        close(frame.loc[np.asarray(use)[present],'BH'],bh[present],1e-12);close(frame.loc[np.asarray(use)[present],'BY'],by[present],1e-12);close(frame.loc[np.asarray(use)[present],'MC_upper_BH'],mc[present],1e-12)
    assert frame.loc[frame.P.isna(),['BH','BY','MC_upper_BH']].isna().all().all();exclude=frame.CI_low*frame.CI_high>0;priority=frame.BH.lt(.05)&frame.MC_upper_BH.lt(.05)&exclude&frame.LOO_same_sign;assert np.array_equal(frame.CI_excludes_zero,exclude) and np.array_equal(frame.priority_association,priority)


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2);source=root/'immune_RNA_MCP_new_source_associations_20261019_v2';out=root/'immune_RNA_MCP_new_source_association_audit_20261019_v3';assert not out.exists();out.mkdir(mode=0o700)
    terminal=wait(root,out,'test_immune_RNA_MCP_new_source_family_v2',source/'private_terminal.json',{'COMPLETE_ALL420_NEW_SOURCE_AND_CLINICAL_RESULTS_AND_FULL_RANDOM_STREAMS_AUDIT_PENDING'})
    for name,v in terminal['files'].items():assert sha(source/name)==v['sha256']
    _,_,_,registry,arrays=read(root);main=pd.read_csv(source/'private_all420_registered_main_results.tsv',sep='\t');raw=pd.read_csv(source/'private_all48_unadjusted_clinical_rank_results.tsv',sep='\t');sens=pd.read_csv(source/'private_all3360_prespecified_point_sensitivities.tsv',sep='\t').set_index(['model_id','mode']);assert len(main)==420 and len(raw)==48 and len(sens)==3360
    assert main.model_id.tolist()==registry.model_id.tolist();assert raw.model_id.tolist()==registry.loc[registry.family.str.startswith('clinical'),'model_id'].tolist();midx=main.set_index('model_id');ridx=raw.set_index('model_id');cal=pd.read_csv(root/'immune_RNA_MCP_new_source_family_null_check_20261019_v1/private_complete_original_model_design_mapping.tsv',sep='\t').set_index('model_id');streams=int(registry.metadata_preeligible.sum()+registry.loc[registry.family.str.startswith('clinical'),'metadata_preeligible'].sum());main_streams=0;rank_streams=0;bootstrap_values=0;maxerr=0.
    for r in registry.itertuples(index=False):
        a=arrays[r.model_id];z=midx.loc[r.model_id];clinical=r.family.startswith('clinical');status=r.status
        for field in ['family','dataset','context_id','pair','program','contrast','endpoint','n_units','materials','metadata_preeligible','kind','binary_endpoint']:assert z[field]==getattr(r,field)
        if r.metadata_preeligible:
            saved=np.load(source/('private_'+r.model_id+'_complete_stream.npz'));ix,key=indices(r,a);assert np.array_equal(saved['indices'],ix) and saved['cases'].tolist()==a['cases'].tolist() and str(saved['bootstrap_key'])==key
            if r.status=='ESTIMABLE' and r.family=='source_marginal':
                e=rankref(a['x'],a['y']);vals=permutations(a['x'],a['y'],r.model_id);b=np.array([rankref(a['x'][j],a['y'][j]) for j in ix]);loo=np.array([rankref(np.delete(a['x'],j),np.delete(a['y'],j)) for j in range(len(a['x']))]);close(saved['main_statistics'],vals);close(saved['bootstrap_effects'],b);close(saved['leave_one'],loo);assert np.array_equal(saved['bootstrap_reasons'],np.where(np.isfinite(b),'ESTIMABLE','CONSTANT_BOOTSTRAP_TARGET'));check_p(z,vals,e,streams);check_summary(z,b,loo,e);close(z.effect,e);main_streams+=len(vals);bootstrap_values+=len(b)
            elif r.status=='ESTIMABLE':
                if not bool(cal.loc[r.model_id,'pass_all_null_scenarios']):status='NULL_DIAGNOSTIC_FAILED_NO_ADJUSTED_INFERENCE'
                else:
                    m,st=independent_design(a['x'],a['y'],a['C'],a['names'],r.kind);assert st=='ESTIMABLE'
                    try:g,vals=independent_wild(m,19999,seed(r.model_id,'HC3_wild'))
                    except ValueError as err:status=str(err)
                    else:
                        e=g['effect'];b,reasons=bootstrap(a['x'],a['y'],a['C'],a['names'],r.kind,ix,m['nuisance_rank'],independent=True);loo,lr=leave_one(a['x'],a['y'],a['C'],a['names'],r.kind,m['nuisance_rank'],independent=True);close(saved['main_statistics'],vals);close(saved['bootstrap_effects'],b);close(saved['leave_one'],loo);assert np.array_equal(saved['bootstrap_reasons'],reasons) and np.array_equal(saved['leave_one_reasons'],lr);check_p(z,vals,g['t'],streams);check_summary(z,b,loo,e);close([z.effect,z.HC3_t],[e,g['t']]);maxerr=max(maxerr,float(np.max(np.abs(saved['main_statistics']-vals))));main_streams+=len(vals);bootstrap_values+=len(b)
            if clinical:
                rr=ridx.loc[r.model_id]
                if min(np.ptp(a['x']),np.ptp(a['y']))>0:
                    e=rankref(a['x'],a['y'],r.binary_endpoint);vals=permutations(a['x'],a['y'],r.model_id+'::unadjusted',r.binary_endpoint);b=np.array([rankref(a['x'][j],a['y'][j],r.binary_endpoint) for j in ix]);loo=np.array([rankref(np.delete(a['x'],j),np.delete(a['y'],j),r.binary_endpoint) for j in range(len(a['x']))]);assert rr.status=='ESTIMABLE';close(saved['unadjusted_statistics'],vals);close(saved['unadjusted_bootstrap_effects'],b);close(saved['unadjusted_leave_one'],loo);close(rr.effect,e);check_p(rr,vals,e,streams);check_summary(rr,b,loo,e);rank_streams+=len(vals);bootstrap_values+=len(b)
                else:assert rr.status=='CONSTANT_UNADJUSTED_TARGET_OR_READOUT' and pd.isna(rr.P)
        else:
            assert pd.isna(z.P) and pd.isna(z.effect)
            if clinical:assert ridx.loc[r.model_id,'status']=='INPUT_NOT_METADATA_PREELIGIBLE'
        assert z.status==status;check_sensitivities(r,a,sens)
        if r.metadata_preeligible and clinical and r.kind=='clinical' and np.isfinite(z.effect):
            aa=sens.loc[(r.model_id,'component_A'),'effect'];bb=sens.loc[(r.model_id,'component_B'),'effect'];close(z.effect,.5*(aa+(1 if r.contrast=='shared_mean' else -1)*bb))
        if main.index[main.model_id.eq(r.model_id)][0]%24==0:(out/'private_progress.json').write_text(json.dumps(dict(phase='INDEPENDENT_COMPLETE_RESULTS_AND_RANDOM_STREAMS',models_checked=int(main.index[main.model_id.eq(r.model_id)][0])+1,total_models=420,partial_review_allowed=False,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    correction(main);correction(raw,False)
    receipt=dict(status='ACCEPTED_ALL420_MAIN_48_RAW_3360_SENSITIVITIES_AND_COMPLETE_RANDOM_STREAMS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),resource_guard=guard,registered_main_models=420,unadjusted_clinical_slots=48,all_point_sensitivities=3360,main_random_statistics_replayed=main_streams,unadjusted_clinical_random_statistics_replayed=rank_streams,bootstrap_effects_replayed=bootstrap_values,maximum_wild_statistic_error=maxerr,no_new_gene_fit_or_identity_change=True)
    (out/'private_terminal.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(dict(status=receipt['status'],terminal_sha256=sha(out/'private_terminal.json'),maximum_error=maxerr)))


if __name__=='__main__':main()
