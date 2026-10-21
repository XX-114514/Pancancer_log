"""Finite source-case T-state coordination/clinical family, all arrays RAM only.

Selected external pairs come only from the pre-effect private contract. These
are observed marker patterns, not original CM scores or causal interactions.
"""
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import hashlib
import itertools
import json
import os
import re
import subprocess
import numpy as np
import observed_T_marker_evidence_v1 as evidence
import GSE288199_state_pair_math_v1 as math
import GSE288199_state_pair_R_audit_v1 as reference

STATES=['R_TNFRSF9','R_LAG3','R_cycling','E_dysfunction','E_cycling','E_memory']
MAIN_QC='native_QC500_250_MT20'
MAIN_PARENT='dual_regulatory_evidence'


def clean(x):
    if isinstance(x,np.ndarray):return clean(x.tolist())
    if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
    if isinstance(x,np.integer):return int(x)
    if isinstance(x,dict):return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [clean(v) for v in x]
    return x


def agent(value):
    return re.sub(r'\(C\d+(?:-C\d+)?\)','',str(value)).strip()


def select(payload,qc=MAIN_QC,parent=MAIN_PARENT,minimum=20,main_ids=None,cancer=None,require_other=False):
    rows=[r for r in payload['rows'] if r['QC']==qc and r['parent']==parent]
    assert len(rows)==243 and len({r['case_index'] for r in rows})==243
    base=[r for r in rows if r['clinical_eligible'] and min(r['metrics'][1:3])>=minimum and
        (main_ids is None or r['case_index'] in main_ids) and (cancer is None or r['cancer_type']==cancer) and
        (not require_other or r['parent_counts']['O']>=20)]
    count=Counter(agent(r['drug']) for r in base)
    chosen=[r for r in base if count[agent(r['drug'])]>=10]
    chosen=sorted(chosen,key=lambda r:r['case_index'])
    return chosen,dict(clinical_and_both_parent_capture=len(base),retained_after_agent_group10=len(chosen),
        removed_for_agent_group_below10=len(base)-len(chosen),agent_categories_retained=sum(n>=10 for n in count.values()),
        response0=sum(r['response']==0 for r in chosen),response1=sum(r['response']==1 for r in chosen),
        reviewed_independent_patients=0)


def value(rows,state,other=False):
    if other:
        j=['R_TNFRSF9','R_LAG3','R_cycling'].index(state)
        return np.asarray([r['other_T_states'][j]/r['parent_counts']['O'] for r in rows])
    j=evidence.METRICS.index(state);k=1 if state.startswith('R_') else 2
    return np.asarray([r['metrics'][j]/r['metrics'][k] for r in rows])


def controls(rows,parents,mode='full',response=False):
    if not rows:return np.empty((0,0))
    n=len(rows);categories=sorted({agent(r['drug']) for r in rows})
    clinical=np.asarray([r['clinical'] for r in rows],float).reshape(n,8)
    dummy=np.asarray([[agent(r['drug'])==a for a in categories] for r in rows],float)
    quality=np.asarray([sum([r['quality'][p] for p in parents],[]) for r in rows],float).reshape(n,4*len(parents))
    c=np.empty((n,0)) if mode=='unadjusted' else np.c_[clinical,dummy] if mode=='clinical_only' else np.c_[clinical,quality,dummy]
    return np.c_[c,[r['response'] for r in rows]] if response else c


def qualification(rows):
    if len(rows)<50:return 'FEWER_THAN50_AUTHOR_CASE_UNITS'
    n0=sum(r['response']==0 for r in rows);n1=sum(r['response']==1 for r in rows)
    if min(n0,n1)<10:return 'FEWER_THAN10_IN_ONE_PATHOLOGY_GROUP'
    return 'ELIGIBLE_SOURCE_CASE_DESIGN'


def association(rows,x,y,c):
    why=qualification(rows)
    if why!='ELIGIBLE_SOURCE_CASE_DESIGN':return None,None,why
    return math.association(x,y,c)


def stratified_indices(response,n_boot,label):
    response=np.asarray(response);rng=np.random.default_rng(math.seed('NSCLC_state_v1|bootstrap|'+label))
    groups=[np.flatnonzero(response==g) for g in [0,1]]
    assert all(len(g)>=10 for g in groups)
    return np.column_stack([rng.choice(g,size=(n_boot,len(g)),replace=True) for g in groups])


def replay(label,x,y,c,f,signs=None,wild=None,indices=None,bootstrap=None):
    error=0.;wild_count=0;boot_count=0;blocks=0
    expected=[f[k] for k in ['beta','se','t','rho']] if f else None
    def one(**kwargs):
        nonlocal error,blocks
        result=reference.audit([clean(dict(label=label,x=x,y=y,C=c,expected=expected,**kwargs))])
        assert len(result['blocks'])==1
        for k,v in result['blocks'][0].items():
            if 'error' in k:error=max(error,float(v))
        blocks+=1
    if signs is not None:
        for start in range(0,len(signs),1000):
            one(signs=signs[start:start+1000],wild=wild[start:start+1000]);wild_count+=len(signs[start:start+1000])
    if indices is not None:
        one(indices=indices,bootstrap=bootstrap);boot_count=len(indices)
    if signs is None and indices is None:one()
    return dict(label=label,status='PASS_ALL_SUBMITTED_REAL_RANDOM_STREAMS_R_QR',R_blocks=blocks,
        full_wild_statistics_checked=wild_count,full_bootstrap_statistics_checked=boot_count,max_numeric_error=error)


def run(payload,task,plan):
    code=Path(__file__).parent
    rscript=Path(plan['Rscript_path']);assert rscript.is_file() and rscript.name=='Rscript'
    os.environ['PANCANCER_R_BLOCK_REFERENCE']=str(rscript)
    for n,d in plan['code_SHA256'].items():assert hashlib.sha256((code/n).read_bytes()).hexdigest()==d
    pairs=plan['fixed_pairs'];assert len(pairs) in [3,9] and len({(r,e) for r,e,s in pairs})==len(pairs)
    assert all(r in STATES[:3] and e in STATES[3:] and s in [-1,1,None] for r,e,s in pairs)
    if len(pairs)==9:assert {(r,e) for r,e,s in pairs}==set(itertools.product(STATES[:3],STATES[3:]))
    out=task/'state_statistics_v1';out.mkdir(exist_ok=False)
    (task/'native_measurement_v1/private_checkpoint.json').write_text(json.dumps(dict(
        stage='FROZEN_STATISTICS_MODULE_RUNNING_CASE_ARRAYS_RAM_ONLY',observed_UTC=datetime.now(timezone.utc).isoformat()),indent=2)+'\n')
    save=lambda n,v:(out/n).write_text(json.dumps(clean(v),indent=2,allow_nan=False)+'\n')
    started=datetime.now(timezone.utc).isoformat();rows,coverage=select(payload)
    save('private_main_same_case_coverage.json',coverage)
    ids={r['case_index'] for r in rows};response=np.asarray([r['response'] for r in rows])
    indices=stratified_indices(response,2000,'same_main_cases') if qualification(rows)=='ELIGIBLE_SOURCE_CASE_DESIGN' else None
    marginals=[];coordinates=[];audits=[];calibrations=[];sensitivities=[]

    def evaluate(label,x,y,c):
        save('private_checkpoint.json',dict(stage='FIXED_CLINICAL_AND_COORDINATION_FAMILY',label=label,started_UTC=started))
        m,f,status=association(rows,x,y,c);result=dict(label=label,n=len(rows),status=status)
        if f is None:
            if qualification(rows)=='ELIGIBLE_SOURCE_CASE_DESIGN':audits.append(replay(label,x,y,c,None))
            return result
        result.update({k:f[k] for k in ['beta','se','t','rho']})
        rng=np.random.default_rng(math.seed('NSCLC_state_v1|wild|'+label))
        signs=np.where(rng.random((19999,len(x)))<.5,-1,1).astype(np.int8)
        wild=math.wild_statistics(m,f,signs)
        result['raw_P']=float((1+np.sum(np.abs(wild)>=abs(f['t'])-1e-12))/20000)
        boot,reasons=math.bootstrap(x,y,c,indices);valid=np.isfinite(boot)
        result.update(CI95=np.quantile(boot[valid],[.025,.975]).tolist() if valid.sum()>=1900 else None,
            bootstrap_valid=int(valid.sum()),bootstrap_status_counts=dict(Counter(reasons)),
            nuisance_rank=m['rank'],residual_df=m['df'])
        audits.append(replay(label,x,y,c,f,signs,wild,indices,boot))
        cal=math.calibration(m,'NSCLC_state_v1|'+label)
        for z in cal:z.pop('simulation_pvalues')
        calibrations.append(dict(label=label,scenarios=cal))
        loo=[]
        for i in range(len(x)):
            take=np.arange(len(x))!=i;rr=[r for j,r in enumerate(rows) if j!=i]
            _,ff,why=association(rr,x[take],y[take],c[take]);loo.append(ff['rho'] if ff else np.nan)
        good=np.isfinite(loo)
        result['leave_one_author_case_out']=dict(estimable=int(good.sum()),
            same_direction=int(np.sum(np.sign(np.asarray(loo)[good])==np.sign(f['rho']))),
            rho_range=[float(np.nanmin(loo)),float(np.nanmax(loo))] if good.any() else None)
        return result

    for state in STATES:
        c=controls(rows,['R' if state.startswith('R_') else 'E'])
        result=evaluate('clinical|'+state,response,value(rows,state),c);result['state']=state;marginals.append(result)
    for rs,es,direction in pairs:
        c=controls(rows,['R','E']);result=evaluate('coordination|'+rs+'|'+es,value(rows,rs),value(rows,es),c)
        result.update(R_state=rs,E_state=es,external_covariance_direction=direction);coordinates.append(result)
    lookup={r['state']:r for r in marginals};joint=[]
    for coord,(rs,es,direction) in zip(coordinates,pairs):
        a,b=lookup[rs],lookup[es];accepted=all(z['status']=='ESTIMABLE' for z in [a,b,coord])
        joint.append(dict(R_state=rs,E_state=es,n=len(rows),status='ESTIMABLE' if accepted else 'AT_LEAST_ONE_FIXED_COMPONENT_NOT_ESTIMABLE',
            triple_IUT_P=max(a.get('raw_P',1.),b.get('raw_P',1.),coord.get('raw_P',1.)),
            clinical_both_endpoint_IUT_P=max(a.get('raw_P',1.),b.get('raw_P',1.)),
            covariance_direction_matches_external=bool(np.sign(coord['rho'])==direction) if 'rho' in coord and direction is not None else None,
            clinical_direction_product_matches_covariance=bool(np.sign(a['rho']*b['rho'])==np.sign(coord['rho'])) if accepted else None))
    for table,pkey,prefix in [(marginals,'raw_P','auxiliary6'),(coordinates,'raw_P','coordination'+str(len(pairs))),(joint,'triple_IUT_P','primary'+str(len(pairs)))]:
        pv=[r.get(pkey,1.) for r in table]
        for r,bh,by in zip(table,math.adjustments(pv),math.adjustments(pv,by=True)):r.update({prefix+'_BH':float(bh),prefix+'_BY':float(by)})
    # Fixed point-only settings, never alternative significance families.
    settings=[('unadjusted',dict(), 'unadjusted'),('clinical_background_only',dict(),'clinical_only'),
        ('all_retained_same_main_cases',dict(qc='all_retained_native',main_ids=ids),'full'),
        ('costate_parent_same_main_cases',dict(parent='costate_regulatory_sensitivity',main_ids=ids),'full'),
        ('parent_min10',dict(minimum=10),'full'),('parent_min30',dict(minimum=30),'full'),
        ('LUAD_same_main_cases',dict(main_ids=ids,cancer='LUAD'),'full'),('LUSC_same_main_cases',dict(main_ids=ids,cancer='LUSC'),'full')]
    for setting,options,mode in settings:
        rr,_=select(payload,**options);yy=np.asarray([r['response'] for r in rr])
        for state in STATES:
            _,ff,why=association(rr,yy,value(rr,state),controls(rr,['R' if state.startswith('R_') else 'E'],mode))
            sensitivities.append(dict(setting=setting,component='clinical',state=state,n=len(rr),status=why,rho=ff['rho'] if ff else None))
        for rs,es,direction in pairs:
            _,ff,why=association(rr,value(rr,rs),value(rr,es),controls(rr,['R','E'],mode))
            sensitivities.append(dict(setting=setting,component='coordination',R_state=rs,E_state=es,n=len(rr),status=why,rho=ff['rho'] if ff else None))
    for rs,es,direction in pairs:
        _,ff,why=association(rows,value(rows,rs),value(rows,es),controls(rows,['R','E'],response=True))
        sensitivities.append(dict(setting='coordination_plus_posttreatment_response_diagnostic',component='coordination',R_state=rs,E_state=es,n=len(rows),status=why,rho=ff['rho'] if ff else None))
    # Matched other-T comparison: paired bootstrap, both relations use all3-parent controls.
    refrows,refcoverage=select(payload,main_ids=ids,require_other=True);refresults=[]
    refindices=stratified_indices([r['response'] for r in refrows],2000,'matched_other_T') if qualification(refrows)=='ELIGIBLE_SOURCE_CASE_DESIGN' else None
    for rs,es,direction in pairs:
        c=controls(refrows,['R','E','O']);y=value(refrows,es);fits=[];boots=[];record=dict(R_state=rs,E_state=es,n=len(refrows))
        for source in ['R','O']:
            x=value(refrows,rs,other=source=='O');m,f,why=association(refrows,x,y,c);fits.append(f)
            if f:
                b,reasons=math.bootstrap(x,y,c,refindices);boots.append(b)
                audits.append(replay('matched_other_T|'+source+'|'+rs+'|'+es,x,y,c,f,indices=refindices,bootstrap=b))
            else:boots.append(None)
            record[source+'_status']=why
        if all(f is not None for f in fits):
            delta=boots[0]-boots[1];valid=np.isfinite(delta)
            record.update(R_rho=fits[0]['rho'],other_T_rho=fits[1]['rho'],delta_R_minus_other=fits[0]['rho']-fits[1]['rho'],
                paired_bootstrap_valid=int(valid.sum()),delta_CI95=np.quantile(delta[valid],[.025,.975]).tolist() if valid.sum()>=1900 else None)
        refresults.append(record)
    pv=dict(marginal=[r.get('raw_P',1.) for r in marginals],coordination=[r.get('raw_P',1.) for r in coordinates],primary=[r['triple_IUT_P'] for r in joint])
    rcode="suppressPackageStartupMessages(library(jsonlite));a<-fromJSON(file('stdin'));cat(toJSON(lapply(a,function(x)list(BH=p.adjust(x,'BH'),BY=p.adjust(x,'BY'))),digits=16))"
    result=subprocess.run([os.environ['PANCANCER_R_BLOCK_REFERENCE'],'--vanilla','-e',rcode],input=json.dumps(pv),capture_output=True,text=True,check=True,timeout=60)
    check=json.loads(result.stdout)
    for name,values in pv.items():
        for field,by in [('BH',False),('BY',True)]:np.testing.assert_allclose(check[name][field],math.adjustments(values,by=by),rtol=0,atol=1e-13)
    assert len(marginals)==6 and len(coordinates)==len(pairs) and len(joint)==len(pairs)
    assert len(sensitivities)==8*(6+len(pairs))+len(pairs) and len(refresults)==len(pairs)
    accepted=bool(calibrations) and all(s['passed_finite_alarm'] for r in calibrations for s in r['scenarios'])
    save('private_all6_clinical_state_associations.json',marginals);save('private_pair_state_coordination.json',coordinates)
    save('private_joint_coordination_clinical_IUT.json',joint);save('private_point_sensitivities.json',sensitivities)
    save('private_matched_other_T_coverage.json',refcoverage);save('private_matched_other_T_differences.json',refresults)
    save('private_all_real_design_R_random_stream_audit.json',audits);save('private_all_fixed_working_null_diagnostics.json',calibrations)
    status=('COMPLETE_NO_ESTIMABLE_MAIN_COMPONENTS' if not calibrations else
        'COMPLETE_FIXED_FAMILY_WORKING_NULL_ACCEPTED' if accepted else 'COMPLETE_FIXED_FAMILY_NULL_ALARM_INFERENCE_NOT_ACCEPTED')
    save('private_terminal.json',dict(status=status,
        started_UTC=started,completed_UTC=datetime.now(timezone.utc).isoformat(),
        primary_slots=len(pairs),auxiliary_clinical_slots=6,coordination_slots=len(pairs),
        point_sensitivity_slots=len(sensitivities),all_fixed_slots_retained=True,
        all_fitted_main_design_full19999_wild_and2000_bootstrap_R_replayed=True,
        all_fitted_matched_reference_full2000_bootstrap_R_replayed=True,complete_family_BH_BY_R_p_adjust_checked=True,
        individual_rows_scores_and_resampling_arrays_RAM_only=True,not_reviewed_independent_patient_units=True,
        source_posttreatment_association_not_baseline_prediction_or_direct_communication=True,
        clinical_unique_increment_and_formal_state_identities_not_tested=True,
        frozen_statistics_SHA256=hashlib.sha256((task/'private_frozen_NSCLC_state_statistics_v1.json').read_bytes()).hexdigest()))
    (task/'native_measurement_v1/private_checkpoint.json').write_text(json.dumps(dict(
        stage='FROZEN_STATE_STATISTICS_COMPLETE_CASE_ARRAYS_RELEASE_ON_EXIT',observed_UTC=datetime.now(timezone.utc).isoformat()),indent=2)+'\n')
    print(json.dumps(dict(status='FIXED_NSCLC_STATE_FAMILY_FINISHED',working_null_accepted=accepted)),flush=True)
