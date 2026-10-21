"""Conditional source associations beyond author mixture and other-CD4 markers.

This is an exploratory follow-up, not a redefinition of a closed family.
All case, clinical and resampling objects remain in RAM.
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
import NSCLC_author_state_family_statistics_v1 as prior
import NSCLC_observed_T_state_statistics_v1 as established
import GSE288199_state_pair_math_v1 as math

NAMESPACE='NSCLC_state_vs_subtype_v1'


def select(payload,subtype=None,ids=None,cancer=None):
    rows,coverage=prior.select(payload,ids=ids,cancer=cancer,other=True)
    if subtype:
        rows=[r for r in rows if r['subtype_measurements'][subtype]['n']>=20]
        counts=Counter(established.agent(r['drug']) for r in rows)
        rows=[r for r in rows if counts[established.agent(r['drug'])]>=10]
    coverage.update(retained_for_current_question=len(rows),response1=sum(r['response']==1 for r in rows),response0=sum(r['response']==0 for r in rows))
    return rows,coverage


def controls(rows,plan,mode='full',subtype=None):
    if not rows:return np.empty((0,0))
    base=prior.controls(rows,other=True)
    defs=plan['author_lineage_names']
    mixture=lambda names,parent:np.asarray([[r['subtype_measurements'][n]['n']/r['parent_counts'][parent] for n in names] for r in rows])
    r_mix=mixture(defs['Treg'],'Treg')
    o_mix=mixture(defs['CD4T_nonTreg'],'CD4T_nonTreg')
    o_marker=np.asarray([[v/r['parent_counts']['CD4T_nonTreg'] for v in r['marker_counts']['CD4T_nonTreg']] for r in rows])
    arrays=[base]
    if mode in ['R_mixture_only','full']:arrays.append(r_mix)
    if mode in ['other_CD4_only','full']:arrays.extend([o_mix,o_marker])
    assert mode in ['base','R_mixture_only','other_CD4_only','full']
    if subtype:arrays.append(np.asarray([r['subtype_measurements'][subtype]['quality'] for r in rows],float))
    return np.column_stack(arrays)


def specifications(plan):
    specs=prior.specifications(plan)
    for name in plan['author_lineage_names']['Treg']:
        for marker in plan['R_markers']:
            specs.append(dict(kind='within_subtype_clinical',subtype=name,marker=marker,
                              label='within_subtype_clinical|'+name+'|'+marker))
    assert len(specs)==112 and len({r['label'] for r in specs})==112
    return specs


def vectors(spec,rows,payload,plan,mode='full'):
    subtype=spec.get('subtype');c=controls(rows,plan,mode,subtype)
    if subtype:
        j=plan['R_markers'].index(spec['marker'])
        y=np.asarray([r['subtype_measurements'][subtype]['markers'][j]/r['subtype_measurements'][subtype]['n'] for r in rows])
        return np.asarray([r['response'] for r in rows]),y,c
    if spec['kind']=='marginal':
        s=spec['state'];y=prior.marker(rows,s,plan) if s in plan['R_states'] else prior.author(rows,s,payload)
        return np.asarray([r['response'] for r in rows]),y,c
    r=prior.marker(rows,spec['R_state'],plan);e=prior.author(rows,spec['E_state'],payload)
    if spec['kind']=='coordination':return r,e,c
    response=np.asarray([r['response'] for r in rows])
    return (response,r,np.c_[c,e]) if spec['kind']=='conditional_R' else (response,e,np.c_[c,r])


def full_component(spec,rows,payload,plan,indices):
    x,y,c=vectors(spec,rows,payload,plan);m,f,why=prior.fit(rows,x,y,c)
    out=dict(**spec,n=len(rows),status=why);audit=None;calibration=None
    label=NAMESPACE+'|'+spec['label']
    if f is None:
        if established.qualification(rows)=='ELIGIBLE_SOURCE_CASE_DESIGN':audit=established.replay(label,x,y,c,None)
        return out,audit,calibration
    out.update({k:f[k] for k in ['rho','beta','se','t']})
    signs=np.where(np.random.default_rng(math.seed(label+'|wild')).random((19999,len(rows)))<.5,-1,1).astype(np.int8)
    wild=math.wild_statistics(m,f,signs);boot,reasons=math.bootstrap(x,y,c,indices);good=np.isfinite(boot)
    out.update(raw_P=float((1+np.sum(np.abs(wild)>=abs(f['t'])-1e-12))/20000),
        CI95=np.quantile(boot[good],[.025,.975]).tolist() if good.sum()>=1900 else None,
        bootstrap_valid=int(good.sum()),bootstrap_status_counts=dict(Counter(reasons)),
        nuisance_rank=m['rank'],residual_df=m['df'])
    audit=established.replay(label,x,y,c,f,signs,wild,indices,boot)
    trials=math.calibration(m,label)
    for row in trials:row.pop('simulation_pvalues')
    calibration=dict(label=spec['label'],scenarios=trials)
    effects=[];missing=Counter()
    for j in range(len(rows)):
        take=np.arange(len(rows))!=j;sub=[r for i,r in enumerate(rows) if i!=j]
        _,ff,status=prior.fit(sub,x[take],y[take],c[take])
        if ff is None:missing[status]+=1
        else:effects.append(ff['rho'])
    out['leave_one_author_case_out']=dict(estimable=len(effects),same_direction=sum(np.sign(z)==np.sign(f['rho']) for z in effects),
        rho_range=[min(effects),max(effects)] if effects else None,missing_status_counts=dict(missing))
    return out,audit,calibration


def run(payload,task,plan):
    code=Path(__file__).parent
    for name,digest in plan['code_SHA256'].items():assert hashlib.sha256((code/name).read_bytes()).hexdigest()==digest
    os.environ['PANCANCER_R_BLOCK_REFERENCE']=plan['Rscript_path']
    out=task/'statistics_v1';out.mkdir(exist_ok=False)
    save=lambda n,v:(out/n).write_text(json.dumps(established.clean(v),indent=2,allow_nan=False)+'\n')
    started=datetime.now(timezone.utc).isoformat();main,coverage=select(payload);ids={r['case_index'] for r in main}
    cohorts={None:main};coverages={'main':coverage};indices={}
    for s in plan['author_lineage_names']['Treg']:
        cohorts[s],coverages[s]=select(payload,s,ids=ids)
    for key,rows in cohorts.items():
        indices[key]=established.stratified_indices([r['response'] for r in rows],2000,NAMESPACE+'|common|'+str(key)) if established.qualification(rows)=='ELIGIBLE_SOURCE_CASE_DESIGN' else None
    save('private_all_fixed_cohort_coverage.json',coverages)
    specs=specifications(plan);results=[];audits=[];calibrations=[]
    def work(spec):
        key=spec.get('subtype');return full_component(spec,cohorts[key],payload,plan,indices[key])
    with ThreadPoolExecutor(max_workers=2) as executor:
        for row,audit,cal in executor.map(work,specs):
            results.append(row)
            if audit is not None:audits.append(audit)
            if cal is not None:calibrations.append(cal)
            save('private_checkpoint.json',dict(stage='FULL_FIXED_COMPONENT',completed=len(results),expected=112,last_component=row['label'],observed_UTC=datetime.now(timezone.utc).isoformat()))
    families=prior.tables([r for r in results if r['kind']!='within_subtype_clinical'],plan)
    fine=[r for r in results if r['kind']=='within_subtype_clinical']
    for field,by in [('within9_BH',False),('within9_BY',True)]:
        for r,q in zip(fine,math.adjustments([r.get('raw_P',1.) for r in fine],by=by)):r[field]=float(q)
    families.append(('within_subtype_clinical9',fine,'raw_P'))
    pvalues={name:[r.get(key,1.) for r in table] for name,table,key in families}
    rcode="suppressPackageStartupMessages(library(jsonlite));a<-fromJSON(file('stdin'));cat(toJSON(lapply(a,function(x)list(BH=p.adjust(x,'BH'),BY=p.adjust(x,'BY'))),digits=16))"
    r=subprocess.run([plan['Rscript_path'],'--vanilla','-e',rcode],input=json.dumps(pvalues),capture_output=True,text=True,timeout=60,check=True)
    checks=json.loads(r.stdout)
    for key,values in pvalues.items():
        for field,by in [('BH',False),('BY',True)]:np.testing.assert_allclose(checks[key][field],math.adjustments(values,by),atol=1e-13,rtol=0)
    for name,table,key in families:save('private_'+name+'.json',table)
    points=[]
    for setting,mode,cancer in [('base','base',None),('R_mixture_only','R_mixture_only',None),('other_CD4_only','other_CD4_only',None),('full','full',None),('LUAD_full','full','LUAD'),('LUSC_full','full','LUSC')]:
        for spec in specs:
            rows,_=select(payload,spec.get('subtype'),ids,cancer);x,y,c=vectors(spec,rows,payload,plan,mode)
            _,f,why=prior.fit(rows,x,y,c);points.append(dict(setting=setting,component=spec['label'],n=len(rows),status=why,rho=f['rho'] if f else None))
    for spec in specs:
        if spec['kind']=='coordination':
            x,y,c=vectors(spec,main,payload,plan);_,f,why=prior.fit(main,x,y,np.c_[c,[r['response'] for r in main]])
            points.append(dict(setting='full_plus_postresponse_diagnostic',component=spec['label'],n=len(main),status=why,rho=f['rho'] if f else None))
    assert len(points)==702
    save('private_all702_point_diagnostics.json',points);save('private_all112_components.json',results)
    save('private_all_real_random_stream_R_audits.json',audits);save('private_all_fixed_working_null_diagnostics.json',calibrations)
    accepted=bool(calibrations) and all(s['passed_finite_alarm'] for r in calibrations for s in r['scenarios'])
    save('private_terminal.json',dict(status='COMPLETE_SUBTYPE_OTHER_CD4_FAMILY_WORKING_NULL_ACCEPTED' if accepted else 'COMPLETE_FAMILY_NULL_UNACCEPTED_OR_NO_MAIN_FITS',started_UTC=started,completed_UTC=datetime.now(timezone.utc).isoformat(),all112_fixed_slots_retained=True,all_fitted_designs_complete19999wild2000bootstrap_R_replayed=True,all_complete_families_R_adjusted=True,individual_and_resampling_arrays_RAM_only=True,original_CM_identity_and_closed_families_unchanged=True,source_exploratory_not_prediction_or_direct_interaction=True))
    print(json.dumps(dict(status='SUBTYPE_OTHER_CD4_FAMILY_FINISHED',working_null_accepted=accepted)),flush=True)
