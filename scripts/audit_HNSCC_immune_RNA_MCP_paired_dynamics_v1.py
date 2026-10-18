"""All12 true results/48 sensitivities and every random value independently replayed."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta,binom,binomtest
from statsmodels.stats.multitest import multipletests

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from HNSCC_immune_RNA_MCP_paired_io_v1 import wait,read
from HNSCC_immune_RNA_MCP_paired_audit_math_v1 import reference_model,reference_estimate,reference_wild,reference_bootstrap


def equal(a,b):assert np.allclose(a,b,rtol=0,atol=1e-8,equal_nan=True)


def sd(slot,stream):return int.from_bytes(hashlib.sha256(('2026101903|'+str(slot)+'|'+stream).encode()).digest()[:8],'little')


def point(x,y,C,names):
    good=np.isfinite(np.column_stack([x,y,C])).all(axis=1);m,status=reference_model(x[good],y[good],C[good],names)
    if m is None:return np.nan,status,int(good.sum())
    rx=m['x']-m['Q']@m['Q'].T@m['x'];ry=m['y']-m['Q']@m['Q'].T@m['y']
    return float(np.corrcoef(rx,ry)[0,1]),status,int(good.sum())


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'HNSCC_immune_RNA_MCP_paired_dynamics_20261019_v1';out=root/'HNSCC_immune_RNA_MCP_paired_dynamics_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all_true_random_and_results_replay.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,scope='All12 exact sign/median or cochange results,48points, all19999 wild/2000 bootstrap/LOO per eligible cochange, all family adjustments and candidate gates; separate pandas/pivotQR/statsmodels and scipy exact binomial.'),indent=2)+'\n')
    terminal=wait(root,out,'test_HNSCC_immune_RNA_MCP_paired_dynamics_v1',source/'private_terminal.json',{'COMPLETE_ALL12_PAIRED_DYNAMICS_AND48_SENSITIVITIES_FULL_RANDOM_AUDIT_PENDING'})
    for n,v in terminal['files'].items():assert sha(source/n)==v['sha256']
    _,_,plan,_,arrays=read(root);table=pd.read_csv(source/'private_all12_paired_dynamics_associations.tsv',sep='\t').set_index('slot');sensitivity=pd.read_csv(source/'private_all48_point_sensitivities.tsv',sep='\t').set_index(['slot','mode']);assert len(table)==12 and len(sensitivity)==48
    pp=np.full(12,np.nan);upper=np.full(12,np.nan);eff=np.full(12,np.nan);cis=np.full((12,2),np.nan);loostable=np.zeros(12,bool);wildcount=0;bootcount=0;maxwild=0.;maxboot=0.;sens_checked=0
    for item in plan['slots']:
        i=item['slot'];a=arrays[i];r=table.loc[i];kind=item['kind'];assert r.kind==kind and r.program==item['program'] and r.n_patients==len(a['patients'])
        if kind=='paired_median_shift':
            d=a['delta'];used=d[np.abs(d)>1e-10];N=len(d);n=len(used)
            status='FEWER_THAN20_COMPLETE_PAIRED_PATIENTS' if N<20 else 'FEWER_THAN10_NONZERO_PAIRED_DELTAS' if n<10 else 'ESTIMABLE';assert r.status==status
            if status=='ESTIMABLE':
                pos=int(np.sum(used>0));p=float(binomtest(pos,n,.5,alternative='two-sided').pvalue);k=max(k for k in range(1,N//2+1) if 2*binom.cdf(k-1,N,.5)<=.05);ordered=np.sort(d);ci=[ordered[k-1],ordered[N-k]];rho=float(np.median(d));ranks=pd.Series(np.abs(used)).rank().to_numpy();rb=float(np.sum(ranks*np.sign(used))/ranks.sum())
                equal([r.effect,r.P,r.MC_upper_P,r.CI_low,r.CI_high,r.mean_delta,r.paired_rank_biserial,r.exact_interval_minimum_coverage],[rho,p,p,*ci,d.mean(),rb,1-2*binom.cdf(k-1,N,.5)])
                assert r.positive==pos and r.negative==n-pos and r.numerical_zero==N-n and r.nonzero==n
                pp[i]=upper[i]=p;eff[i]=rho;cis[i]=ci
            for mode,key in [('all_captured_mean','all_captured_delta'),('fixed_stress_complement','fixed_stress_delta')]:
                d=a[key];nz=d[np.abs(d)>1e-10];s=sensitivity.loc[(i,mode)];equal(s.point,np.median(d));assert s.positive==int((nz>0).sum()) and s.negative==int((nz<0).sum()) and s.numerical_zero==len(d)-len(nz) and s.n_patients==len(d) and s.status=='POINT_ONLY';sens_checked+=1
        else:
            x,y,C,names=a['x'],a['y'],a['C'],a['names'];m,status=reference_model(x,y,C,names);assert status==r.status
            if m is not None:
                b,t,rho=reference_estimate(m);saved=np.load(source/('private_slot%02d_all_random_streams.npz'%i));wv=reference_wild(m['y'],m['D'],19999,sd(i,'wild'));equal(saved['wild'],wv);maxwild=max(maxwild,float(np.max(np.abs(saved['wild']-wv))))
                count=int((np.abs(wv)>=abs(t)-1e-12).sum());p=(1+count)/20000;up=1. if count==19999 else float(beta.ppf(1-.01/6,count+1,19999-count));up=min(1.,(1+19999*up)/20000)
                indices=np.random.default_rng(sd(i,'patient_bootstrap')).integers(0,len(x),size=(2000,len(x)),dtype=np.int16);assert np.array_equal(indices,saved['bootstrap_indices'])
                bv,bs=reference_bootstrap(x,y,C,names,indices,m['nuisance_rank']);equal(saved['bootstrap_values'],bv);assert np.array_equal(saved['bootstrap_statuses'],bs);valid=np.isfinite(bv);be=float(np.nanmax(np.abs(saved['bootstrap_values']-bv))) if valid.any() else 0.;maxboot=max(maxboot,be);ci=np.quantile(bv[valid],[.025,.975]) if valid.sum()>=1800 else [np.nan,np.nan]
                lv=[];ls=[]
                for j in range(len(x)):
                    keep=np.arange(len(x))!=j;mm,ss=reference_model(x[keep],y[keep],C[keep],names,expected_rank=m['nuisance_rank']);ls.append(ss)
                    if mm is None:lv.append(np.nan)
                    else:
                        rx=mm['x']-mm['Q']@mm['Q'].T@mm['x'];ry=mm['y']-mm['Q']@mm['Q'].T@mm['y'];lv.append(float(np.corrcoef(rx,ry)[0,1]))
                lv=np.asarray(lv);equal(saved['LOO_values'],lv);assert np.array_equal(saved['LOO_statuses'],ls);lg=np.isfinite(lv);stable=bool(lg.all() and (np.sign(lv)==np.sign(rho)).all())
                equal([r.effect,r.P,r.MC_upper_P,r.HC3_t,r.CI_low,r.CI_high],[rho,p,up,t,*ci]);assert r.wild_exceedances==count and r.valid_bootstrap==int(valid.sum()) and r.LOO_valid==int(lg.sum()) and r.LOO_same_sign==stable and r.residual_df==m['residual_df'] and r.nuisance_rank==m['nuisance_rank']
                equal([r.LOO_min,r.LOO_max],[np.nanmin(lv),np.nanmax(lv)] if lg.any() else [np.nan,np.nan]);pp[i]=p;upper[i]=up;eff[i]=rho;cis[i]=ci;loostable[i]=stable;wildcount+=19999;bootcount+=2000
            for mode in plan['cochange_point_sensitivities']:
                if mode=='unadjusted':
                    value=float(pd.Series(x).rank().corr(pd.Series(y).rank())) if min(np.ptp(x),np.ptp(y))>0 else np.nan;status='POINT_ONLY' if np.isfinite(value) else 'CONSTANT_TARGET_OR_READOUT';n=len(x)
                elif mode=='treatment_HPV_only':value,status,n=point(x,y,C[:,:2],names[:2])
                elif mode=='both_RNA_quality_changes_only':
                    ix=[2,3,6,7];value,status,n=point(x,y,C[:,ix],names[ix])
                elif mode=='all_captured_mean':value,status,n=point(a['all_captured_x'],a['all_captured_y'],C,names)
                elif mode=='fixed_stress_complement':value,status,n=point(a['fixed_stress_x'],a['fixed_stress_y'],C,names)
                else:value,status,n=point(x,y,a['C_plus_composition'],np.r_[names,['composition_A','composition_B']])
                s=sensitivity.loc[(i,mode)];equal(s.point,value);assert s.status==status and s.n_patients==n;sens_checked+=1
        (out/'private_progress.json').write_text(json.dumps(dict(phase='INDEPENDENT_ALL_TRUE_EFFECTS_AND_RANDOM_REPLAY',slots_checked=i+1,slots=12,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    bh=multipletests(np.nan_to_num(pp,nan=1.),method='fdr_bh')[1];by=multipletests(np.nan_to_num(pp,nan=1.),method='fdr_by')[1];ub=multipletests(np.nan_to_num(upper,nan=1.),method='fdr_bh')[1]
    for v in [bh,by,ub]:v[~np.isfinite(pp)]=np.nan
    equal(table.BH12.to_numpy(),bh);equal(table.BY12.to_numpy(),by);equal(table.MC_upper_BH12.to_numpy(),ub)
    candidate=(table.kind.to_numpy()=='paired_cochange')&(bh<.05)&(ub<.05)&((cis[:,0]>0)|(cis[:,1]<0))&loostable;assert np.array_equal(table.priority_cochange_candidate,candidate) and sens_checked==48
    t=dict(status='ACCEPTED_ALL12_PAIRED_RESULTS_48_SENSITIVITIES_AND_COMPLETE_RANDOM_STREAMS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),all_slots_replayed=12,all_sensitivities_replayed=sens_checked,all_wild_values_replayed=wildcount,all_bootstrap_values_replayed=bootcount,maximum_wild_error=maxwild,maximum_bootstrap_error=maxboot)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
