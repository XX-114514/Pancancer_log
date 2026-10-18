"""All true effects, sensitivities and every numerical random value via QR."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_distribution

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait,read
from NSCLC_immune_RNA_MCP_cohort_math_v1 import seed
from NSCLC_immune_RNA_MCP_cohort_audit_math_v1 import prepare,evaluate,wild,boot


def equal(a,b):assert np.allclose(a,b,atol=1e-8,rtol=0,equal_nan=True)


def refpoint(x,y,C,names,kind,response):
    if min(np.count_nonzero(response==0),np.count_nonzero(response==1))<10:return np.nan,'FEWER_THAN10_IN_ONE_PATHOLOGY_GROUP'
    m,status=prepare(x,y,C,names,kind)
    if m is None:return np.nan,status
    if np.linalg.norm(m['y']-m['Q']@(m['Q'].T@m['y']))<=1e-8*max(1.,np.linalg.norm(m['y'])):return np.nan,'READOUT_EXPLAINED_BY_FULL_CONTROL_SPAN'
    try:return evaluate(m)['effect'],status
    except ValueError as e:return np.nan,str(e)


def qvalues(values,by=False):
    p=np.asarray(values,float);m=len(p);answer=np.empty(m)
    for i in range(m):
        ordered=np.sort(p,kind='stable');positions=np.flatnonzero(ordered>=p[i]);lower=int(positions[0])+1
        answer[i]=min(1.,min(ordered[j-1]*m/j for j in range(lower,m+1))*(np.sum(1/np.arange(1,m+1)) if by else 1))
    return answer


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'NSCLC_immune_RNA_MCP_cohort_associations_20261019_v1'
    out=root/'NSCLC_immune_RNA_MCP_cohort_associations_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all_true_random_stream_replay.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        scope='All18 component effects,all6 IUT/MC/BH/BY/priority rules,all126 point sensitivities,and every wild/bootstrap/LOO numerical value. Independent pandas ranks and pivotedQR/statistic implementation. No partial scientific review.'),indent=2)+'\n')
    terminal=wait(root,out,'test_NSCLC_immune_RNA_MCP_cohort_v1',source/'private_terminal.json',{'COMPLETE_ALL6_JOINT_BUNDLES_18_COMPONENTS_AND_FULL_RANDOM_STREAMS_ACCEPTANCE_PENDING'})
    for name,item in terminal['files'].items():assert sha(source/name)==item['sha256']
    _,_,_,scope,arrays=read(root)
    components=pd.read_csv(source/'private_all18_component_associations.tsv',sep='\t').set_index(['slot','component'])
    joint=pd.read_csv(source/'private_all6_joint_cohort_clinical_bundles.tsv',sep='\t').set_index('slot')
    sensitivities=pd.read_csv(source/'private_all126_point_sensitivities.tsv',sep='\t').set_index(['slot','component','mode'])
    assert len(components)==18 and len(joint)==6 and len(sensitivities)==126
    wild_count=0;boot_count=0;maxwild=0.;maxboot=0.;completed=0
    for item in scope.itertuples(index=False):
        a=arrays[(item.slot,item.component)];r=components.loc[(item.slot,item.component)];kind=str(a['kind'])
        m,status=prepare(a['x'],a['y'],a['C'],a['names'],kind)
        if item.status!='ESTIMABLE':m=None;status=item.status
        if m is not None and np.linalg.norm(m['y']-m['Q']@(m['Q'].T@m['y']))<=1e-8*max(1.,np.linalg.norm(m['y'])):m=None;status='READOUT_EXPLAINED_BY_FULL_CONTROL_SPAN'
        assert r.status==status and r.n_cases==len(a['x'])
        if m is not None:
            g,values=wild(m,19999,seed(item.slot,item.component,'wild'));equal(r.effect,g['effect']);equal(r.HC3_t,g['t'])
            z=np.load(source/('private_slot%02d_%s_complete_stream.npz'%(item.slot,item.component)),allow_pickle=False)
            assert z['cases'].tolist()==a['cases'].tolist() and int(z['wild_seed'])==seed(item.slot,item.component,'wild')
            error=float(np.max(np.abs(values-z['wild_statistics'])));assert error<1e-8;maxwild=max(maxwild,error)
            count=int((np.abs(values)>=abs(g['t'])-1e-12).sum());assert count==r.wild_exceedances;equal(r.P,(1+count)/20000)
            upper=1. if count==19999 else float(beta_distribution.ppf(1-.01/18,count+1,19999-count));equal(r.MC_upper_P,min(1.,(1+19999*upper)/20000))
            rng=np.random.default_rng(seed(item.slot,'bundle','bootstrap'));indices=np.empty((2000,len(a['x'])),np.int16);offset=0
            for group in sorted(set(a['strata'])):
                ix=np.flatnonzero(a['strata']==group);indices[:,offset:offset+len(ix)]=ix[rng.integers(0,len(ix),size=(2000,len(ix)))];offset+=len(ix)
            saved=np.load(source/('private_slot%02d_shared_bootstrap_indices.npz'%item.slot),allow_pickle=False)
            assert np.array_equal(saved['indices'],indices) and saved['cases'].tolist()==a['cases'].tolist() and str(saved['numpy_version'])==np.__version__
            b,reasons=boot(a['x'],a['y'],a['C'],a['names'],kind,indices,m['nuisance_rank'])
            assert reasons.tolist()==z['bootstrap_reasons'].tolist();equal(b,z['bootstrap_effects']);valid=np.isfinite(b)
            if valid.any():maxboot=max(maxboot,float(np.max(np.abs(b[valid]-z['bootstrap_effects'][valid]))))
            ci=np.quantile(b[valid],[.025,.975]) if valid.sum()>=1900 else [np.nan,np.nan];equal(r.CI_low,ci[0]);equal(r.CI_high,ci[1]);assert r.valid_bootstrap==valid.sum()
            loo=[]
            for i in range(len(a['x'])):
                keep=np.arange(len(a['x']))!=i;mm,_=prepare(a['x'][keep],a['y'][keep],a['C'][keep],a['names'],kind,expected_rank=m['nuisance_rank']);value=np.nan
                if mm:
                    try:value=evaluate(mm)['effect']
                    except ValueError:pass
                loo.append(value)
            loo=np.array(loo);equal(loo,z['leave_one']);finite=np.isfinite(loo)
            assert r.LOO_valid==finite.sum() and bool(r.LOO_same_sign)==bool(finite.all() and (loo*g['effect']>0).all())
            equal(r.LOO_min,float(loo[finite].min()) if finite.any() else np.nan);equal(r.LOO_max,float(loo[finite].max()) if finite.any() else np.nan)
            wild_count+=19999;boot_count+=2000
        else:assert pd.isna(r.P) and pd.isna(r.effect)
        for mode in ['unadjusted','clinical_background_only','author_subtype_proportions','all_captured_mean','coordination_plus_response','LUAD_only','LUSC_only']:
            keep=np.ones(len(a['x']),bool);x=a['x'];y=a['y'];C=a['C'];names=a['names'];effect=np.nan
            if mode=='unadjusted':
                ry=pd.Series(y).rank().to_numpy()
                if kind=='coordination':effect=float(np.corrcoef(pd.Series(x).rank().to_numpy(),ry)[0,1]) if len(x)>=3 and min(np.ptp(x),np.ptp(y))>0 else np.nan
                else:
                    n1=int((x==1).sum());n0=int((x==0).sum());effect=float(2*(np.sum(ry[x==1])-n1*(n1+1)/2)/(n1*n0)-1) if n1 and n0 else np.nan
                reason='DESCRIPTIVE_UNADJUSTED' if np.isfinite(effect) else 'INSUFFICIENT_UNADJUSTED_VARIATION'
            elif mode=='coordination_plus_response' and kind!='coordination':reason='NOT_APPLICABLE_TO_CLINICAL_COMPONENT'
            else:
                if mode=='clinical_background_only':columns=np.array([not s.startswith('quality::') for s in names]);C=C[:,columns];names=names[columns]
                elif mode=='author_subtype_proportions':C=a['fineC'];names=a['fine_names']
                elif mode=='all_captured_mean':x=a['allcapture_x'];y=a['allcapture_y']
                elif mode=='coordination_plus_response':C=np.column_stack([C,a['response']]);names=np.concatenate([names,np.array(['response'])])
                elif mode in ['LUAD_only','LUSC_only']:keep=a['cancer']==mode.split('_')[0]
                effect,reason=refpoint(x[keep],y[keep],C[keep],names,kind,a['response'][keep])
            s=sensitivities.loc[(item.slot,item.component,mode)];assert s.n_cases==keep.sum() and s.status==reason;equal(s.effect,effect)
        completed+=1
        (out/'private_progress.json').write_text(json.dumps(dict(stage='INDEPENDENT_ALL_TRUE_STREAMS',completed_components=completed,wild_values_replayed=wild_count,bootstrap_values_replayed=boot_count,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    p=[];upper=[]
    for slot in range(6):
        g=components.loc[slot];ready=g.status.eq('ESTIMABLE').all();r=joint.loc[slot]
        assert bool(r.all_three_estimable)==ready
        effect=g.effect;coherent=bool(np.sign(effect.loc['clinical_A']*effect.loc['clinical_B'])==np.sign(effect.loc['coordination'])) if ready else False
        ci=bool((g.CI_low*g.CI_high>0).all()) if ready else False;loo=bool(g.LOO_same_sign.all()) if ready else False
        assert bool(r.coherent_direction)==coherent and bool(r.all_CI_exclude_zero)==ci and bool(r.all_LOO_same_sign)==loo
        value=float(g.P.max()) if ready else np.nan;value_upper=float(g.MC_upper_P.max()) if ready else np.nan
        equal(r.P,value);equal(r.MC_upper_P,value_upper);p.append(value if ready else 1.);upper.append(value_upper if ready else 1.)
    bh=qvalues(p);by=qvalues(p,True);mc=qvalues(upper)
    for slot,r in joint.iterrows():
        if r.all_three_estimable:equal(r.BH6,bh[slot]);equal(r.BY6,by[slot]);equal(r.MC_upper_BH6,mc[slot])
        else:assert pd.isna(r.BH6) and pd.isna(r.BY6) and pd.isna(r.MC_upper_BH6)
        expected=bool(bh[slot]<.05 and mc[slot]<.05 and r.coherent_direction and r.all_CI_exclude_zero and r.all_LOO_same_sign)
        assert bool(r.priority_joint_candidate)==expected
    accepted=dict(status='ACCEPTED_ALL6_JOINT_BUNDLES_18_COMPONENTS_126_SENSITIVITIES_AND_COMPLETE_RANDOM_STREAMS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),
        wild_values_replayed=wild_count,bootstrap_values_replayed=boot_count,max_wild_error=maxwild,max_bootstrap_error=maxboot,clinical_or_causal_validation=False)
    (out/'private_terminal.json').write_text(json.dumps(accepted,indent=2)+'\n');print(json.dumps(dict(status=accepted['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
