"""All six joint programme/clinical bundles after complete input/null gates."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_distribution,rankdata,spearmanr

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait,read
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design,estimate,wild_values,takes,bootstrap,leave_one,adjust,seed


def point(x,y,C,names,kind,response):
    if min(np.count_nonzero(response==0),np.count_nonzero(response==1))<10:return np.nan,'FEWER_THAN10_IN_ONE_PATHOLOGY_GROUP'
    m,status=design(x,y,C,names,kind)
    if m is None:return np.nan,status
    if np.linalg.norm(m['y']-m['Q']@(m['Q'].T@m['y']))<=1e-8*max(1.,np.linalg.norm(m['y'])):return np.nan,'READOUT_EXPLAINED_BY_FULL_CONTROL_SPAN'
    try:_,effect=estimate(m);return effect,status
    except ValueError as e:return np.nan,str(e)


def sensitivities(item,a):
    rows=[];kind=str(a['kind']);names=a['names'];n=len(a['x'])
    for mode in ['unadjusted','clinical_background_only','author_subtype_proportions','all_captured_mean','coordination_plus_response','LUAD_only','LUSC_only']:
        keep=np.ones(n,bool);x=a['x'];y=a['y'];C=a['C'];nn=names;effect=np.nan;status='NOT_ESTIMATED'
        if mode=='unadjusted':
            if kind=='coordination':effect=float(spearmanr(x,y).statistic) if n>=3 and min(np.ptp(x),np.ptp(y))>0 else np.nan
            else:
                n1=int((x==1).sum());n0=int((x==0).sum())
                effect=float(2*(rankdata(y)[x==1].sum()-n1*(n1+1)/2)/(n1*n0)-1) if n1 and n0 else np.nan
            status='DESCRIPTIVE_UNADJUSTED' if np.isfinite(effect) else 'INSUFFICIENT_UNADJUSTED_VARIATION'
        elif mode=='coordination_plus_response' and kind!='coordination':status='NOT_APPLICABLE_TO_CLINICAL_COMPONENT'
        else:
            if mode=='clinical_background_only':columns=np.array([not s.startswith('quality::') for s in names]);C=C[:,columns];nn=names[columns]
            elif mode=='author_subtype_proportions':C=a['fineC'];nn=a['fine_names']
            elif mode=='all_captured_mean':x=a['allcapture_x'];y=a['allcapture_y']
            elif mode=='coordination_plus_response':C=np.c_[C,a['response']];nn=np.concatenate([names,np.array(['response'])])
            elif mode in ['LUAD_only','LUSC_only']:keep=a['cancer']==mode.split('_')[0]
            effect,status=point(x[keep],y[keep],C[keep],nn,kind,a['response'][keep])
        rows.append(dict(slot=item.slot,component=item.component,mode=mode,n_cases=int(keep.sum()),effect=effect,status=status,no_new_P_or_CI=True))
    return rows


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,32*1024**2)
    out=root/'NSCLC_immune_RNA_MCP_cohort_associations_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_any_true_joint_effect.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        scope='All6 frozen pair bundles/18 component tests, same case sets, all126 prespecified point sensitivities. No model or candidate replacement.',
        inference='19999 restricted-null HC3 wild,2000 same-slot histology-agent casebootstrap,>=1900 valid CI,LOO;all6 IUT max(Pcoord,PA,PB) BH/BY and18simultaneous MCupperP.',
        seed=2026101901,residual_readout_gate='Norm after complete nuisance projection >1e-8*max(1,norm(y)); otherwise missing, never near-zero numerical correlation.',
        random_storage='All main numerical wild/bootstrap/LOO values,all bootstrap indices/case order,seeds/NumPy version. No genes or row data public.',
        limits='Author-case observational post-treatment association. No newly reviewed independent patient,baseline prediction,causal cell interaction or receptor-binding claim.'),indent=2)+'\n')
    cal=root/'NSCLC_immune_RNA_MCP_cohort_calibration_audit_20261019_v2/private_terminal.json'
    wait(root,out,'audit_NSCLC_immune_RNA_MCP_cohort_calibration_v2',cal,{'ACCEPTED_ALL_NULL_CALIBRATION_VALUES_AND_PASS_GATE'})
    _,inputaudit,_,scope,arrays=read(root);rows=[];sensitivity=[];indices={}
    for item in scope.itertuples(index=False):
        a=arrays[(item.slot,item.component)];kind=str(a['kind'])
        row=dict(slot=item.slot,program=item.program,view_A=item.view_A,view_B=item.view_B,component=item.component,kind=kind,n_cases=item.n_cases,response0=item.response0,response1=item.response1,status=item.status,
            effect=np.nan,HC3_t=np.nan,P=np.nan,MC_upper_P=np.nan,CI_low=np.nan,CI_high=np.nan,valid_bootstrap=0,LOO_valid=0,LOO_same_sign=False,LOO_min=np.nan,LOO_max=np.nan)
        if item.status=='ESTIMABLE':
            m,status=design(a['x'],a['y'],a['C'],a['names'],kind);assert status=='ESTIMABLE'
            residual=m['y']-m['Q']@(m['Q'].T@m['y'])
            if np.linalg.norm(residual)<=1e-8*max(1.,np.linalg.norm(m['y'])):row['status']='READOUT_EXPLAINED_BY_FULL_CONTROL_SPAN'
            else:
                g,effect=estimate(m);values=wild_values(g,19999,seed(item.slot,item.component,'wild'));count=int((np.abs(values)>=abs(g['t'])-1e-12).sum());p=(1+count)/20000
                upper=1. if count==19999 else float(beta_distribution.ppf(1-.01/18,count+1,19999-count));upper=min(1.,(1+19999*upper)/20000)
                if item.slot not in indices:
                    ix=takes(a['strata'],2000,seed(item.slot,'bundle','bootstrap'));indices[item.slot]=(a['cases'].copy(),ix)
                    np.savez_compressed(out/('private_slot%02d_shared_bootstrap_indices.npz'%item.slot),indices=ix,cases=a['cases'],strata=a['strata'],numpy_version=np.array(np.__version__))
                assert indices[item.slot][0].tolist()==a['cases'].tolist();ix=indices[item.slot][1]
                boot,reasons=bootstrap(a['x'],a['y'],a['C'],a['names'],kind,ix,m['nuisance_rank']);valid=np.isfinite(boot)
                ci=np.quantile(boot[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan]
                loo=leave_one(a['x'],a['y'],a['C'],a['names'],kind,m['nuisance_rank']);finite=np.isfinite(loo)
                row.update(status='ESTIMABLE',effect=effect,HC3_t=g['t'],P=p,MC_upper_P=upper,CI_low=float(ci[0]),CI_high=float(ci[1]),valid_bootstrap=int(valid.sum()),wild_exceedances=count,
                    LOO_valid=int(finite.sum()),LOO_same_sign=bool(finite.all() and (loo*effect>0).all()),LOO_min=float(loo[finite].min()) if finite.any() else np.nan,LOO_max=float(loo[finite].max()) if finite.any() else np.nan)
                np.savez_compressed(out/('private_slot%02d_%s_complete_stream.npz'%(item.slot,item.component)),wild_statistics=values,bootstrap_effects=boot,bootstrap_reasons=reasons,leave_one=loo,wild_seed=np.array(seed(item.slot,item.component,'wild'),np.uint64),cases=a['cases'])
        rows.append(row);sensitivity.extend(sensitivities(item,a))
        (out/'private_progress.json').write_text(json.dumps(dict(stage='ALL_FIXED_MAIN_COMPONENTS',completed_models=len(rows),planned_models=18,partial_effect_review_allowed=False,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    components=pd.DataFrame(rows);assert len(components)==18 and len(sensitivity)==126
    bundles=[]
    for slot,g in components.groupby('slot',sort=True):
        ready=g.status.eq('ESTIMABLE').all() and len(g)==3;effects=g.set_index('component').effect
        sign=bool(np.sign(effects['clinical_A']*effects['clinical_B'])==np.sign(effects['coordination'])) if ready else False
        bundles.append(dict(slot=slot,program=g.program.iloc[0],view_A=g.view_A.iloc[0],view_B=g.view_B.iloc[0],n_cases=int(g.n_cases.iloc[0]),all_three_estimable=ready,
            P=float(g.P.max()) if ready else np.nan,MC_upper_P=float(g.MC_upper_P.max()) if ready else np.nan,
            coherent_direction=sign,all_CI_exclude_zero=bool((g.CI_low*g.CI_high>0).all()) if ready else False,all_LOO_same_sign=bool(g.LOO_same_sign.all()) if ready else False))
    joint=pd.DataFrame(bundles);assert len(joint)==6
    for name,by in [('BH6',False),('BY6',True)]:joint[name]=adjust(joint.P.fillna(1),by);joint.loc[joint.P.isna(),name]=np.nan
    joint['MC_upper_BH6']=adjust(joint.MC_upper_P.fillna(1));joint.loc[joint.P.isna(),'MC_upper_BH6']=np.nan
    joint['priority_joint_candidate']=joint.BH6.lt(.05)&joint.MC_upper_BH6.lt(.05)&joint.coherent_direction&joint.all_CI_exclude_zero&joint.all_LOO_same_sign
    components.to_csv(out/'private_all18_component_associations.tsv',sep='\t',index=False)
    joint.to_csv(out/'private_all6_joint_cohort_clinical_bundles.tsv',sep='\t',index=False)
    pd.DataFrame(sensitivity).to_csv(out/'private_all126_point_sensitivities.tsv',sep='\t',index=False)
    files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()};assert sum(v['bytes'] for v in files.values())<32*1024**2
    t=dict(status='COMPLETE_ALL6_JOINT_BUNDLES_18_COMPONENTS_AND_FULL_RANDOM_STREAMS_ACCEPTANCE_PENDING',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_acceptance_sha256=sha(inputaudit),calibration_acceptance_sha256=sha(cal),files=files,
        independent_patient_or_causal_interaction_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
