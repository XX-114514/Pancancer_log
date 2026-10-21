"""Pure synthetic assembly and full R replay of every component kind."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from pathlib import Path
import argparse
import copy
import json
import os
import numpy as np
import NSCLC_state_vs_subtype_statistics_v1 as family
import NSCLC_author_state_family_statistics_v1 as prior
import NSCLC_observed_T_state_statistics_v1 as established


def fixture():
    genes=['TNFRSF9','LAG3','MKI67'];R=['authorTreg_'+g+'_pos' for g in genes]
    rt=['synthetic_R_type_'+str(i) for i in range(3)];ot=['synthetic_O_type_'+str(i) for i in range(6)];E=['synthetic_E_type_'+str(i) for i in range(10)]
    plan=dict(R_markers=genes,R_states=R,E_states=E,all30_pairs=[[r,e] for r in R for e in E],author_lineage_names=dict(Treg=rt,CD8T=E,CD4T_nonTreg=ot))
    rng=np.random.default_rng(374261);rows=[]
    for i in range(243):
        response=i%2;latent=rng.normal()+.5*response
        rc=rng.multinomial(800+i%97,rng.dirichlet([8,15,4]));oc=rng.multinomial(1000+i%53,rng.dirichlet(np.full(6,7.)))
        er=np.exp(np.linspace(-.5,.5,10)*latent);er/=er.sum();ec=rng.multinomial(2000+i%41,er)
        Rn=int(rc.sum());On=int(oc.sum());En=int(ec.sum());total=Rn+On+En;sub={}
        for j,name in enumerate(rt+ot):
            n=int(rc[j] if j<3 else oc[j-3]);prob=np.clip(np.array([.3,.2,.15])+np.array([.01,-.02,.02])*latent+rng.normal(0,.02,3),.01,.8)
            sub[name]=dict(n=n,markers=rng.binomial(n,prob).tolist(),quality=[float(rng.uniform(6,8)),float(rng.uniform(4,6)),float(np.log1p(n)),n/total])
        parents=dict(Treg=Rn,CD8T=En,CD4T_nonTreg=On)
        quality={p:[float(rng.uniform(6,8)),float(rng.uniform(4,6)),float(np.log1p(n)),n/total] for p,n in parents.items()}
        rows.append(dict(case_index=i,QC=prior.MAIN_QC,clinical_eligible=i<160,response=response,cancer_type='LUAD' if i<70 else 'LUSC',drug='synthetic_agent(C3)',
            clinical=[float(rng.integers(45,80)),float(rng.integers(0,2)),float(rng.integers(0,2)),float(i>=70),float(rng.integers(1,4)),float(rng.integers(2,5)),float(rng.integers(0,2)),float(rng.integers(0,2))],
            quality=quality,parent_counts=parents,author_counts=np.r_[rc,ec].tolist(),subtype_measurements=sub,
            marker_counts={p:np.sum([sub[x]['markers'] for x in names],axis=0).tolist() for p,names in [('Treg',rt),('CD4T_nonTreg',ot)]}))
    return dict(rows=rows+[dict(copy.deepcopy(r),QC='all_author_retained') for r in rows],state_names=rt+E),plan


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);a=parser.parse_args();a.output.mkdir(exist_ok=False)
    payload,plan=fixture();main,_=family.select(payload);assert len(main)==160
    ids={r['case_index'] for r in main};specs=family.specifications(plan);assert len(specs)==112
    assembled=[]
    for i,spec in enumerate(specs):
        rows,_=family.select(payload,spec.get('subtype'),ids);x,y,c=family.vectors(spec,rows,payload,plan)
        _,f,why=prior.fit(rows,x,y,c);assert why=='ESTIMABLE'
        assembled.append(dict(**spec,n=len(rows),status=why,rho=f['rho'],raw_P=(i+1)/113))
    tables=prior.tables(assembled[:103],plan);assert [len(t) for n,t,k in tables]==[13,30,60,30]
    c0=family.controls(main,plan,'base');cf=family.controls(main,plan,'full');assert cf.shape[1]==c0.shape[1]+12
    missing=copy.deepcopy(payload)
    for row in missing['rows']:
        if row['case_index']>=30:row['subtype_measurements'][plan['author_lineage_names']['Treg'][2]]['n']=10
    rr,_=family.select(missing,plan['author_lineage_names']['Treg'][2]);assert len(rr)==30
    assert established.qualification(rr)=='FEWER_THAN50_AUTHOR_CASE_UNITS'
    chosen=[next(s for s in specs if s['kind']==kind) for kind in ['marginal','coordination','conditional_R','conditional_E','within_subtype_clinical']]
    def work(spec):
        rows,_=family.select(payload,spec.get('subtype'),ids)
        indices=established.stratified_indices([r['response'] for r in rows],2000,family.NAMESPACE+'|synthetic_QA|'+str(spec.get('subtype')))
        return family.full_component(spec,rows,payload,plan,indices)
    with ThreadPoolExecutor(max_workers=2) as executor:checked=list(executor.map(work,chosen))
    assert all(row['status']=='ESTIMABLE' and audit['full_wild_statistics_checked']==19999 and audit['full_bootstrap_statistics_checked']==2000 for row,audit,cal in checked)
    assert all(s['passed_finite_alarm'] for row,audit,cal in checked for s in cal['scenarios'])
    result=dict(status='PASS_ALL112_POINT_ASSEMBLY_FIVE_COMPLETE_RANDOM_BRANCHES_R_AND_LOW_CAPTURE_GATE',completed_UTC=datetime.now(timezone.utc).isoformat(),human_data_used=False,
        all112_full_synthetic_resampling_not_claimed=True,actual_real_fitted_designs_all_full_R_required=True,
        five_complete_audits=[audit for row,audit,cal in checked],five_working_null_diagnostics=[cal for row,audit,cal in checked])
    (a.output/'private_complete_synthetic_acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['five_complete_audits','five_working_null_diagnostics']}),flush=True)


if __name__=='__main__':main()
