"""New orchestration/conditional branches on synthetic data only.

All103 assembly paths plus four full-resampling main branches and matched
reference streams are checked; this is not all103 full synthetic resampling.
Actual production must independently replay every fitted real design.
"""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from pathlib import Path
import argparse
import copy
import json
import os
import numpy as np
import NSCLC_author_state_family_statistics_v1 as family
import NSCLC_observed_T_state_statistics_v1 as established
import GSE288199_state_pair_math_v1 as math


def fixture():
    R=['authorTreg_'+g+'_pos' for g in ['TNFRSF9','LAG3','MKI67']]
    E=['authorCD8_state_'+str(i) for i in range(10)]
    plan=dict(R_states=R,E_states=E,all30_pairs=[[r,e] for r in R for e in E],synthetic_only=True)
    labels=['authorTreg_type_'+str(i) for i in range(3)]+E
    rng=np.random.default_rng(147231);rows=[]
    for i in range(243):
        response=i%2;latent=.5*response+rng.normal(0,.3)
        Rn=100+i%31;En=250+i%53;On=170+i%29;N=Rn+En+On
        stateR=rng.multinomial(Rn,[.3,.4,.3]);prob=np.exp(np.linspace(-.4,.4,10)+np.linspace(-.5,.5,10)*latent);prob/=prob.sum()
        stateE=rng.multinomial(En,prob)
        countR=rng.binomial(Rn,np.clip(np.array([.3,.4,.3])+np.array([.1,-.15,.2])*latent,.05,.9)).tolist()
        countO=rng.binomial(On,[.2,.3,.2]).tolist()
        quality={p:[float(rng.uniform(6,8)),float(rng.uniform(4,6)),float(np.log1p(n)),n/N]
                 for p,n in [('Treg',Rn),('CD8T',En),('CD4T_nonTreg',On)]}
        c=[float(rng.integers(45,80)),float(rng.integers(0,2)),float(rng.integers(0,2)),float(i>=40),
           float(rng.integers(1,4)),float(rng.integers(2,5)),float(rng.integers(0,2)),float(rng.integers(0,2))]
        row=dict(case_index=i,QC=family.MAIN_QC,clinical_eligible=i<80,response=response,cancer_type='LUAD' if i<40 else 'LUSC',
                 drug='synthetic_agent(C3)',clinical=c,parent_counts=dict(Treg=Rn,CD8T=En,CD4T_nonTreg=On),
                 quality=quality,author_counts=np.r_[stateR,stateE].astype(int).tolist(),marker_counts=dict(Treg=countR,CD4T_nonTreg=countO))
        rows.append(row)
    payload=dict(rows=rows+[dict(copy.deepcopy(r),QC='all_author_retained') for r in rows],state_names=labels)
    return payload,plan


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(exist_ok=False)
    payload,plan=fixture();rows,coverage=family.select(payload);assert len(rows)==80
    os.environ['PANCANCER_R_BLOCK_REFERENCE']=os.environ['PANCANCER_R_BLOCK_REFERENCE']
    specs=family.specifications(plan);assembled=[]
    for i,spec in enumerate(specs):
        x,y,c=family.vectors(spec,rows,payload,plan);_,fit,status=family.fit(rows,x,y,c)
        assert status=='ESTIMABLE'
        assembled.append(dict(**spec,n=len(rows),status=status,rho=fit['rho'],raw_P=(i+1)/104))
    tables=family.tables(assembled,plan);assert [len(v) for name,v,key in tables]==[13,30,60,30]
    primary=tables[-1][1]
    lookup={(r['kind'],r.get('R_state'),r.get('E_state')):r for r in assembled if r['kind']!='marginal'}
    for r in primary:
        pair=(r['R_state'],r['E_state'])
        assert r['regulatory_addition_IUT_P']==max(lookup[('coordination',*pair)]['raw_P'],lookup[('conditional_R',*pair)]['raw_P'])
    # Complete attempted random streams for each new branch, two fit workers.
    chosen=[next(r for r in specs if r['kind']==kind) for kind in ['marginal','coordination','conditional_R','conditional_E']]
    indices=established.stratified_indices([r['response'] for r in rows],2000,family.NAMESPACE+'|synthetic_common')
    with ThreadPoolExecutor(max_workers=2) as executor:
        checked=list(executor.map(lambda s:family.full_component(s,rows,payload,plan,indices),chosen))
    assert all(r['status']=='ESTIMABLE' and audit['full_wild_statistics_checked']==19999 and audit['full_bootstrap_statistics_checked']==2000
               for r,audit,calibration in checked)
    reference=[];rstate,estate=plan['all30_pairs'][0];c=family.controls(rows,other=True);y=family.author(rows,estate,payload)
    for which in ['R','O']:
        x=family.marker(rows,rstate,plan,other=which=='O');m,f,status=family.fit(rows,x,y,c);assert f is not None
        boot,reasons=math.bootstrap(x,y,c,indices)
        reference.append(established.replay('synthetic_matched_'+which,x,y,c,f,indices=indices,bootstrap=boot))
    result=dict(status='PASS_SYNTHETIC_ALL103_ASSEMBLY_FOUR_FULL_MAIN_BRANCHES_AND_MATCHED_REFERENCES',
        completed_UTC=datetime.now(timezone.utc).isoformat(),human_data_used=False,all103point_designs_and_complete_families_assembled=True,
        four_main_full19999wild2000bootstrap_R_replayed=True,two_matched_reference_full2000bootstrap_R_replayed=True,
        full103synthetic_resampling_not_claimed=True,
        actual_production_all_fitted_real_streams_still_required=True,
        main_audits=[audit for r,audit,calibration in checked],matched_audits=reference,
        synthetic_working_null_diagnostics=[calibration for r,audit,calibration in checked])
    (a.output/'private_synthetic_complete_acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['main_audits','matched_audits','synthetic_working_null_diagnostics']}),flush=True)


if __name__=='__main__':main()
