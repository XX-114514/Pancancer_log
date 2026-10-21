"""Finite author-lineage marker/CD8 subtype associations; arrays RAM only.

This source-derived exploratory question is distinct from the closed marker
parent family. Conditional rank association is not predictive clinical gain.
"""
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from pathlib import Path
import hashlib
import json
import os
import subprocess
import numpy as np
import GSE288199_state_pair_math_v1 as math
import NSCLC_observed_T_state_statistics_v1 as established
import NSCLC_author_lineage_marker_capture_RAM_v1 as capture

MAIN_QC='native_QC500_250_MT20'
NAMESPACE='NSCLC_author_lineage_markers_v1'
clean=established.clean


def select(payload,minimum=20,qc=MAIN_QC,ids=None,cancer=None,other=False):
    source=[r for r in payload['rows'] if r['QC']==qc]
    assert len(source)==243 and len({r['case_index'] for r in source})==243
    base=[r for r in source if r['clinical_eligible'] and min(r['parent_counts'][p] for p in ['Treg','CD8T'])>=minimum
          and (ids is None or r['case_index'] in ids) and (cancer is None or r['cancer_type']==cancer)
          and (not other or r['parent_counts']['CD4T_nonTreg']>=20)]
    agents=Counter(established.agent(r['drug']) for r in base)
    rows=sorted([r for r in base if agents[established.agent(r['drug'])]>=10],key=lambda r:r['case_index'])
    return rows,dict(clinical_and_both_author_parent_capture=len(base),retained_after_agent_group10=len(rows),
        removed_for_agent_group_below10=len(base)-len(rows),response1=sum(r['response']==1 for r in rows),
        response0=sum(r['response']==0 for r in rows),reviewed_independent_patients=0)


def controls(rows,mode='full',other=False):
    if not rows:return np.empty((0,0))
    names=sorted({established.agent(r['drug']) for r in rows})
    clinical=np.asarray([r['clinical'] for r in rows],float)
    drugs=np.asarray([[established.agent(r['drug'])==v for v in names] for r in rows],float)
    parents=['Treg','CD8T']+(['CD4T_nonTreg'] if other else [])
    quality=np.asarray([sum([r['quality'][p] for p in parents],[]) for r in rows],float)
    return np.empty((len(rows),0)) if mode=='unadjusted' else np.c_[clinical,drugs] if mode=='clinical_only' else np.c_[clinical,quality,drugs]


def marker(rows,state,plan,other=False):
    index=plan['R_states'].index(state);parent='CD4T_nonTreg' if other else 'Treg'
    return np.asarray([r['marker_counts'][parent][index]/r['parent_counts'][parent] for r in rows],float)


def author(rows,state,payload):
    j=payload['state_names'].index(state);parent='Treg' if j<3 else 'CD8T'
    return np.asarray([r['author_counts'][j]/r['parent_counts'][parent] for r in rows],float)


def specifications(plan):
    specs=[dict(label='marginal|'+s,kind='marginal',state=s) for s in plan['R_states']+plan['E_states']]
    for r,e in plan['all30_pairs']:
        for kind in ['coordination','conditional_R','conditional_E']:
            specs.append(dict(label=kind+'|'+r+'|'+e,kind=kind,R_state=r,E_state=e))
    assert len(specs)==103 and len({r['label'] for r in specs})==103
    return specs


def vectors(spec,rows,payload,plan,mode='full'):
    c=controls(rows,mode);response=np.asarray([r['response'] for r in rows],float)
    if spec['kind']=='marginal':
        state=spec['state'];value=marker(rows,state,plan) if state in plan['R_states'] else author(rows,state,payload)
        return response,value,c
    r=marker(rows,spec['R_state'],plan);e=author(rows,spec['E_state'],payload)
    if spec['kind']=='coordination':return r,e,c
    return (response,r,np.c_[c,e]) if spec['kind']=='conditional_R' else (response,e,np.c_[c,r])


def fit(rows,x,y,c):
    why=established.qualification(rows)
    return (None,None,why) if why!='ELIGIBLE_SOURCE_CASE_DESIGN' else math.association(x,y,c)


def full_component(spec,rows,payload,plan,indices):
    x,y,c=vectors(spec,rows,payload,plan);m,f,status=fit(rows,x,y,c)
    result=dict(**spec,n=len(rows),status=status);audit=None;calibration=None
    label=NAMESPACE+'|'+spec['label']
    if f is None:
        if established.qualification(rows)=='ELIGIBLE_SOURCE_CASE_DESIGN':audit=established.replay(label,x,y,c,None)
        return result,audit,calibration
    result.update({k:f[k] for k in ['beta','se','t','rho']})
    signs=np.where(np.random.default_rng(math.seed(label+'|wild')).random((19999,len(rows)))<.5,-1,1).astype(np.int8)
    wild=math.wild_statistics(m,f,signs)
    result['raw_P']=float((1+np.sum(np.abs(wild)>=abs(f['t'])-1e-12))/20000)
    boot,reasons=math.bootstrap(x,y,c,indices);valid=np.isfinite(boot)
    result.update(CI95=np.quantile(boot[valid],[.025,.975]).tolist() if valid.sum()>=1900 else None,
                  bootstrap_valid=int(valid.sum()),bootstrap_status_counts=dict(Counter(reasons)),
                  nuisance_rank=m['rank'],residual_df=m['df'])
    audit=established.replay(label,x,y,c,f,signs,wild,indices,boot)
    simulations=math.calibration(m,label)
    for s in simulations:s.pop('simulation_pvalues')
    calibration=dict(label=spec['label'],scenarios=simulations)
    loo=[];missing=Counter()
    for i in range(len(rows)):
        take=np.arange(len(rows))!=i;subset=[r for j,r in enumerate(rows) if j!=i]
        _,ff,why=fit(subset,x[take],y[take],c[take])
        if ff is None:missing[why]+=1
        else:loo.append(ff['rho'])
    result['leave_one_author_case_out']=dict(estimable=len(loo),same_direction=sum(np.sign(v)==np.sign(f['rho']) for v in loo),
        rho_range=[min(loo),max(loo)] if loo else None,missing_status_counts=dict(missing))
    return result,audit,calibration


def tables(results,plan):
    marginal=[r for r in results if r['kind']=='marginal'];coord=[r for r in results if r['kind']=='coordination']
    conditional=[r for r in results if r['kind'].startswith('conditional_')]
    assert len(marginal)==13 and len(coord)==30 and len(conditional)==60
    lookup={(r['kind'],r.get('R_state'),r.get('E_state')):r for r in results if r['kind']!='marginal'}
    primary=[]
    for r,e in plan['all30_pairs']:
        c=lookup[('coordination',r,e)];a=lookup[('conditional_R',r,e)];b=lookup[('conditional_E',r,e)]
        good=all(v['status']=='ESTIMABLE' for v in [c,a])
        primary.append(dict(R_state=r,E_state=e,n=c['n'],status='ESTIMABLE' if good else 'FIXED_COMPONENT_NOT_ESTIMABLE',
            regulatory_addition_IUT_P=max(c.get('raw_P',1.),a.get('raw_P',1.)),
            coordination_CI_qualified=c.get('CI95') is not None,regulatory_addition_CI_qualified=a.get('CI95') is not None,
            conditional_R_rho=a.get('rho'),conditional_E_rho=b.get('rho'),
            interpretation='Regulatory-state conditional association beyond CD8 state; not predictive gain'))
    families=[('marginal13',marginal,'raw_P'),('coordination30',coord,'raw_P'),('conditional60',conditional,'raw_P'),
              ('primary30',primary,'regulatory_addition_IUT_P')]
    for name,rows,key in families:
        pv=[r.get(key,1.) for r in rows]
        for row,bh,by in zip(rows,math.adjustments(pv),math.adjustments(pv,by=True)):
            row[name+'_BH']=float(bh);row[name+'_BY']=float(by)
    return families


def run(payload,task,plan):
    code=Path(__file__).parent
    for name,digest in plan['code_SHA256'].items():assert hashlib.sha256((code/name).read_bytes()).hexdigest()==digest
    os.environ['PANCANCER_R_BLOCK_REFERENCE']=plan['Rscript_path']
    if not plan.get('synthetic_only'):payload=capture.augment(payload,task,plan)
    out=task/'state_family_v1';out.mkdir(exist_ok=False)
    save=lambda n,v:(out/n).write_text(json.dumps(clean(v),indent=2,allow_nan=False)+'\n')
    started=datetime.now(timezone.utc).isoformat();rows,coverage=select(payload);save('private_main_coverage.json',coverage)
    specs=specifications(plan);response=np.asarray([r['response'] for r in rows]);ids={r['case_index'] for r in rows}
    indices=established.stratified_indices(response,2000,NAMESPACE+'|common_main') if established.qualification(rows)=='ELIGIBLE_SOURCE_CASE_DESIGN' else None
    results=[];audits=[];calibrations=[]
    with ThreadPoolExecutor(max_workers=2) as executor:
        for result,audit,calibration in executor.map(lambda s:full_component(s,rows,payload,plan,indices),specs):
            results.append(result)
            if audit is not None:audits.append(audit)
            if calibration is not None:calibrations.append(calibration)
            save('private_checkpoint.json',dict(stage='COMPLETE_MAIN_COMPONENT',completed_components=len(results),expected_components=103,
                 last_component=result['label'],observed_UTC=datetime.now(timezone.utc).isoformat()))
    families=tables(results,plan)
    pv={name:[r.get(key,1.) for r in table] for name,table,key in families}
    rcode="suppressPackageStartupMessages(library(jsonlite));a<-fromJSON(file('stdin'));cat(toJSON(lapply(a,function(x)list(BH=p.adjust(x,'BH'),BY=p.adjust(x,'BY'))),digits=16))"
    result=subprocess.run([plan['Rscript_path'],'--vanilla','-e',rcode],input=json.dumps(pv),text=True,capture_output=True,check=True,timeout=60)
    adjusted=json.loads(result.stdout)
    for key,values in pv.items():
        for field,by in [('BH',False),('BY',True)]:np.testing.assert_allclose(adjusted[key][field],math.adjustments(values,by=by),atol=1e-13,rtol=0)
    for name,table,key in families:save('private_'+name+'.json',table)
    settings=[('unadjusted',{},'unadjusted'),('clinical_background_only',{},'clinical_only'),
        ('all_retained_same_main',dict(qc='all_author_retained',ids=ids),'full'),('parent_min10',dict(minimum=10),'full'),
        ('parent_min30',dict(minimum=30),'full'),('LUAD_same_main',dict(ids=ids,cancer='LUAD'),'full'),
        ('LUSC_same_main',dict(ids=ids,cancer='LUSC'),'full')]
    points=[]
    for setting,options,mode in settings:
        rr,_=select(payload,**options)
        for spec in specs:
            x,y,c=vectors(spec,rr,payload,plan,mode);_,f,why=fit(rr,x,y,c)
            points.append(dict(setting=setting,component=spec['label'],n=len(rr),status=why,rho=f['rho'] if f else None))
    for r,e in plan['all30_pairs']:
        _,f,why=fit(rows,marker(rows,r,plan),author(rows,e,payload),np.c_[controls(rows),response])
        points.append(dict(setting='coordination_plus_posttreatment_response_diagnostic',component='coordination|'+r+'|'+e,
                           n=len(rows),status=why,rho=f['rho'] if f else None))
    assert len(points)==751;save('private_all751_point_diagnostics.json',points)
    reference_rows,reference_coverage=select(payload,ids=ids,other=True)
    reference_indices=established.stratified_indices([r['response'] for r in reference_rows],2000,NAMESPACE+'|matched_other') if established.qualification(reference_rows)=='ELIGIBLE_SOURCE_CASE_DESIGN' else None
    references=[]
    for r,e in plan['all30_pairs']:
        c=controls(reference_rows,other=True);y=author(reference_rows,e,payload);fits=[];boots=[]
        record=dict(R_state=r,E_state=e,n=len(reference_rows))
        for which in ['R','O']:
            x=marker(reference_rows,r,plan,other=which=='O');m,f,why=fit(reference_rows,x,y,c);fits.append(f);record[which+'_status']=why
            if f is not None:
                boot,reasons=math.bootstrap(x,y,c,reference_indices);boots.append(boot)
                audits.append(established.replay(NAMESPACE+'|matched|'+which+'|'+r+'|'+e,x,y,c,f,indices=reference_indices,bootstrap=boot))
            else:boots.append(None)
        if all(f is not None for f in fits):
            delta=boots[0]-boots[1];good=np.isfinite(delta)
            record.update(R_rho=fits[0]['rho'],other_CD4T_rho=fits[1]['rho'],rho_difference=fits[0]['rho']-fits[1]['rho'],
                          bootstrap_valid=int(good.sum()),CI95=np.quantile(delta[good],[.025,.975]).tolist() if good.sum()>=1900 else None)
        references.append(record)
        save('private_checkpoint.json',dict(stage='MATCHED_NONREGULATORY_CD4_REFERENCE',completed_pairs=len(references),expected_pairs=30,
             observed_UTC=datetime.now(timezone.utc).isoformat()))
    context=[]
    for r in payload['state_names'][:3]:
        for e in plan['E_states']:
            _,f,why=fit(rows,author(rows,r,payload),author(rows,e,payload),controls(rows))
            context.append(dict(kind='author_subtype_coordination',R_state=r,E_state=e,n=len(rows),status=why,rho=f['rho'] if f else None))
        _,f,why=fit(rows,response,author(rows,r,payload),controls(rows))
        context.append(dict(kind='author_Treg_subtype_clinical',state=r,n=len(rows),status=why,rho=f['rho'] if f else None))
    assert len(context)==33 and len(references)==30
    save('private_all30_matched_other_CD4_references.json',references);save('private_reference_coverage.json',reference_coverage)
    save('private_all33_author_subtype_context_points.json',context);save('private_all_real_random_stream_R_audits.json',audits)
    save('private_all_main_fixed_working_null_diagnostics.json',calibrations)
    accepted=bool(calibrations) and all(s['passed_finite_alarm'] for r in calibrations for s in r['scenarios'])
    terminal=dict(status='COMPLETE_AUTHOR_MARKER_FAMILY_WORKING_NULL_ACCEPTED' if accepted else 'COMPLETE_AUTHOR_MARKER_FAMILY_NULL_UNACCEPTED_OR_NO_MAIN_FITS',
        started_UTC=started,completed_UTC=datetime.now(timezone.utc).isoformat(),all103_main_slots_retained=True,
        all30_primary_and60conditional_and13marginal_and30coordination_family_BH_BY_R_checked=True,
        all_fitted_main_complete19999wild2000bootstrap_and_all_reference2000bootstrap_R_replayed=True,
        individual_inputs_and_resampling_arrays_RAM_only=True,
        source_derived_exploratory_not_independent_validation=True,clinical_rank_increment_not_prediction_or_causal_increment=True,
        original_project_identity_CM_and_closed_families_not_changed=True,
        frozen_contract_SHA256=hashlib.sha256((task/'private_frozen_author_state_family_statistics_v1.json').read_bytes()).hexdigest())
    save('private_terminal.json',terminal);print(json.dumps(dict(status=terminal['status'])),flush=True)
