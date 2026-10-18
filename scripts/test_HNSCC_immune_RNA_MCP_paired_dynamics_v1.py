"""Entire frozen12 paired family after all input/math/null gates, no selection."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_distribution, spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from HNSCC_immune_RNA_MCP_paired_io_v1 import wait,read
from HNSCC_immune_RNA_MCP_paired_math_v1 import median_shift,design,effect,estimate,wild_values,seed,patient_takes,bootstrap,leave_one,adjust


def point(x,y,C,names):
    good=np.isfinite(np.column_stack([x,y,C])).all(axis=1);m,status=design(x[good],y[good],C[good],names)
    return effect(m) if m is not None else np.nan,status,int(good.sum())


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,24*1024**2)
    out=root/'HNSCC_immune_RNA_MCP_paired_dynamics_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_any_true_paired_effect.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,fixed_family='All6 paired median/sign and6partial-rank cochange; all12 BH/BY/MCupperBH, all48 point sensitivities. No source/endpoint/contrast selection.',prerequisite='Complete native RNA/input/math and independent finite-null pass before true effects.'),indent=2)+'\n')
    audit=root/'HNSCC_immune_RNA_MCP_paired_calibration_audit_20261019_v1/private_terminal.json'
    wait(root,out,'audit_HNSCC_immune_RNA_MCP_paired_calibration_v1',audit,{'ACCEPTED_ALL_PAIRED_NULL_SIMULATIONS_AND_PASS_GATE'})
    assert sha(root/'HNSCC_immune_RNA_MCP_paired_math_gate_20261019_v1/private_terminal.json')=='0c8a80c69549ea41806f4323913baa05135dd908d85b2bfb5ba2b2773a61aa4c'
    _,_,plan,scope,arrays=read(root);results=[];sens=[]
    for item in plan['slots']:
        slot=item['slot'];a=arrays[slot];kind=item['kind'];row=dict(**item,n_patients=len(a['patients']),status='UNAVAILABLE',effect=np.nan,P=np.nan,MC_upper_P=np.nan,CI_low=np.nan,CI_high=np.nan,valid_bootstrap=0,LOO_valid=0,LOO_same_sign=False)
        if kind=='paired_median_shift':
            r,status=median_shift(a['delta']);row['status']=status
            if r is not None:row.update(r,P=r['p'],MC_upper_P=r['p']);row.pop('p',None)
            for mode,key in [('all_captured_mean','all_captured_delta'),('fixed_stress_complement','fixed_stress_delta')]:
                d=a[key];s=d[np.abs(d)>1e-10];sens.append(dict(slot=slot,kind=kind,mode=mode,n_patients=len(d),status='POINT_ONLY',point=float(np.median(d)),positive=int((s>0).sum()),negative=int((s<0).sum()),numerical_zero=len(d)-len(s),new_P=False))
        else:
            x,y,C,names=a['x'],a['y'],a['C'],a['names'];m,status=design(x,y,C,names);row['status']=status
            if m is not None:
                g=estimate(m);rho=effect(m);wv=wild_values(g,19999,seed(slot,'wild'));count=int((np.abs(wv)>=abs(g['t'])-1e-12).sum());p=(1+count)/20000
                upper=1. if count==19999 else float(beta_distribution.ppf(1-.01/6,count+1,19999-count));upper=min(1.,(1+19999*upper)/20000)
                indices=patient_takes(len(x),slot);bv,bs=bootstrap(x,y,C,names,indices,m['nuisance_rank']);valid=np.isfinite(bv);ci=np.quantile(bv[valid],[.025,.975]) if valid.sum()>=1800 else [np.nan,np.nan]
                lv,ls=leave_one(x,y,C,names,m['nuisance_rank']);loo=np.isfinite(lv)
                row.update(effect=rho,P=float(p),MC_upper_P=upper,HC3_t=g['t'],wild_exceedances=count,CI_low=float(ci[0]),CI_high=float(ci[1]),valid_bootstrap=int(valid.sum()),LOO_valid=int(loo.sum()),LOO_same_sign=bool(loo.all() and (np.sign(lv)==np.sign(rho)).all()),LOO_min=float(np.nanmin(lv)) if loo.any() else np.nan,LOO_max=float(np.nanmax(lv)) if loo.any() else np.nan,residual_df=m['residual_df'],nuisance_rank=m['nuisance_rank'])
                np.savez_compressed(out/('private_slot%02d_all_random_streams.npz'%slot),wild=wv,bootstrap_indices=indices,bootstrap_values=bv,bootstrap_statuses=bs,LOO_values=lv,LOO_statuses=ls)
            for mode in plan['cochange_point_sensitivities']:
                if mode=='unadjusted':
                    value=float(spearmanr(x,y).statistic) if min(np.ptp(x),np.ptp(y))>0 else np.nan;status='POINT_ONLY' if np.isfinite(value) else 'CONSTANT_TARGET_OR_READOUT';n=len(x)
                elif mode=='treatment_HPV_only':value,status,n=point(x,y,C[:,:2],names[:2])
                elif mode=='both_RNA_quality_changes_only':
                    ix=[2,3,6,7];value,status,n=point(x,y,C[:,ix],names[ix])
                elif mode=='all_captured_mean':value,status,n=point(a['all_captured_x'],a['all_captured_y'],C,names)
                elif mode=='fixed_stress_complement':value,status,n=point(a['fixed_stress_x'],a['fixed_stress_y'],C,names)
                elif mode=='plus_both_frozen_original_node_composition_changes':value,status,n=point(x,y,a['C_plus_composition'],np.r_[names,['composition_A','composition_B']])
                else:raise ValueError('UNDECLARED_SENSITIVITY')
                sens.append(dict(slot=slot,kind=kind,mode=mode,n_patients=n,status=status,point=value,new_P=False))
        results.append(row)
        (out/'private_progress.json').write_text(json.dumps(dict(phase='ALL_FIXED_TRUE_EFFECTS',slots_completed=len(results),total_slots=12,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    table=pd.DataFrame(results);assert len(table)==12
    table['BH12']=adjust(table.P.fillna(1));table['BY12']=adjust(table.P.fillna(1),by=True);table['MC_upper_BH12']=adjust(table.MC_upper_P.fillna(1))
    table.loc[table.P.isna(),['BH12','BY12','MC_upper_BH12']]=np.nan
    table['priority_cochange_candidate']=table.kind.eq('paired_cochange')&table.BH12.lt(.05)&table.MC_upper_BH12.lt(.05)&((table.CI_low>0)|(table.CI_high<0))&table.LOO_same_sign
    table.to_csv(out/'private_all12_paired_dynamics_associations.tsv',sep='\t',index=False);s=pd.DataFrame(sens);assert len(s)==48;s.to_csv(out/'private_all48_point_sensitivities.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()};assert sum(v['bytes'] for v in files.values())<24*1024**2
    t=dict(status='COMPLETE_ALL12_PAIRED_DYNAMICS_AND48_SENSITIVITIES_FULL_RANDOM_AUDIT_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),null_gate_acceptance_sha256=sha(audit),files=files)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
