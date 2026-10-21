"""Finite exploratory clinical family. Individual inputs/random arrays RAM only.

Read-only source transport uses an already configured SSH target and optional
request-scoped ProxyCommand supplied by the caller; no connection changes.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import ast
import hashlib
import json
import os
import shlex
import subprocess
import sys
import time
import numpy as np
import observed_T_marker_evidence_v1 as evidence
import GSE288199_state_pair_math_v1 as math
import GSE288199_state_pair_R_audit_v1 as independent

STATES = ['R_TNFRSF9','R_LAG3','R_cycling','E_dysfunction','E_cycling','E_memory']


def clean(x):
    if isinstance(x,np.ndarray): return clean(x.tolist())
    if isinstance(x,(float,np.floating)): return float(x) if np.isfinite(x) else None
    if isinstance(x,np.integer): return int(x)
    if isinstance(x,dict): return {k:clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [clean(v) for v in x]
    return x


def conjunction(a,b):
    return max(a,b)


def source_payload(code):
    original = (code/'measure_CM02_native_marker_author_correspondence_v3.py').read_text()
    tree = ast.parse(original)
    fn = next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='native_counts')
    source = 'import types\nevidence=types.ModuleType("fixed_evidence")\n'
    source += 'exec('+repr((code/'observed_T_marker_evidence_v1.py').read_text())+',evidence.__dict__)\n'
    source += (code/'HNSCC_observed_marker_clinical_reader_v1.py').read_text()+'\n'
    source += ast.get_source_segment(original,fn)+'\n'
    source += 'print(json.dumps(read_source('+repr(os.environ['PANCANCER_PREFLIGHT_ROOT'])+'),allow_nan=False))\n'
    ssh = ['ssh','-T','-o','ClearAllForwardings=yes','-o','BatchMode=yes','-o','ConnectTimeout=15']
    if os.environ.get('PANCANCER_READONLY_PROXY_COMMAND'):
        ssh += ['-o','ProxyCommand='+os.environ['PANCANCER_READONLY_PROXY_COMMAND']]
    ssh += [os.environ['PANCANCER_READONLY_SSH_TARGET'],shlex.quote(os.environ['PANCANCER_REMOTE_PYTHON'])+' -B -']
    result = subprocess.run(ssh,input=source,capture_output=True,text=True,timeout=3600)
    if result.returncode:
        # Do not leak connection details, paths or patient rows from exceptions.
        raise RuntimeError('READONLY_RAM_SOURCE_FAILED_EXIT_'+str(result.returncode))
    return json.loads(result.stdout)


def select_rows(payload,phase,qc='human_RNA_QC500_250_MT20',minimum=20):
    rows = [r for r in payload['rows'] if r['phase']==phase and r['QC']==qc]
    eligible = [r for r in rows if min(r['metrics'][1:3])>=minimum]
    return rows,eligible


def arrays(rows,state,controls='full'):
    parent = 0 if state.startswith('R_') else 1
    x = np.asarray([r['metrics'][evidence.METRICS.index(state)]/r['metrics'][parent+1] for r in rows])
    y = np.asarray([r['regression'] for r in rows])
    c = np.asarray([r['clinical']+(r['quality'][parent] if controls=='full' else []) for r in rows],float).reshape(len(rows),6 if controls=='full' else 2)
    return x,y,c


def main():
    out = Path(os.environ['PANCANCER_CLINICAL_TASK'])
    contract = json.loads((out/'private_frozen_clinical_implementation_v1.json').read_text())
    pairs_fixed = contract['fixed_external_pairs']
    assert len(pairs_fixed)==3 and len({(r,e) for r,e,d in pairs_fixed})==3
    assert all(r in STATES[:3] and e in STATES[3:] and d in [-1,1] for r,e,d in pairs_fixed)
    code = Path(__file__).parent
    sha = lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    for n,d in contract['code_SHA256'].items(): assert sha(code/n)==d
    production=out/'production';production.mkdir(exist_ok=False)
    save=lambda n,v:(production/n).write_text(json.dumps(clean(v),indent=2,allow_nan=False)+'\n')
    started=datetime.now(timezone.utc).isoformat();begin=time.monotonic()
    save('private_checkpoint.json',dict(stage='READONLY_SOURCE_IN_RAM',started_UTC=started))
    payload=source_payload(code)
    save('private_source_input_proof.json',{k:v for k,v in payload.items() if k!='rows'})
    sourcecoverage=[];results=[];audits=[];calibrations=[];sensitivity=[]
    for phase in ['pre','post']:
        allrows,rows=select_rows(payload,phase)
        sourcecoverage.append(dict(phase=phase,reviewed_patients=len(allrows),both_parents_atleast20=len(rows),
            below_R20=sum(r['metrics'][1]<20 for r in allrows),below_E20=sum(r['metrics'][2]<20 for r in allrows),
            continuous_pathology_distinct_values=len({r['regression'] for r in rows})))
        for state in STATES:
            label=phase+'|'+state
            save('private_checkpoint.json',dict(stage='FINITE_CLINICAL_FAMILY',label=label,started_UTC=started))
            x,y,c=arrays(rows,state)
            m,f,status=math.association(x,y,c)
            row=dict(phase=phase,state=state,n=len(rows),status=status)
            if f is not None:
                rng=np.random.default_rng(math.seed('HNSCC_clinical_v1|'+label))
                signs=np.where(rng.random((19999,len(x)))<.5,-1.,1.)
                wild=math.wild_statistics(m,f,signs)
                p=(1+np.sum(np.abs(wild)>=abs(f['t'])-1e-12))/20000
                indices=rng.integers(len(x),size=(2000,len(x)))
                boots,reasons=math.bootstrap(x,y,c,indices)
                valid=np.isfinite(boots)
                row.update({k:f[k] for k in ['beta','se','t','rho']})
                row.update(raw_P=float(p),CI95=np.quantile(boots[valid],[.025,.975]).tolist() if valid.sum()>=1900 else None,
                    bootstrap_valid=int(valid.sum()),bootstrap_status_counts=dict(Counter(reasons)),
                    nuisance_column_space_rank=m['rank'],residual_df=m['df'])
                replay=dict(label=label,x=x,y=y,C=c,expected=[f[k] for k in ['beta','se','t','rho']],
                    signs=signs,wild=wild,indices=indices,bootstrap=boots)
                audits.append(independent.audit([clean(replay)]))
                cal=math.calibration(m,'HNSCC_clinical_v1|'+label)
                for v in cal:v.pop('simulation_pvalues')
                calibrations.append(dict(label=label,scenarios=cal))
                loo=[]
                for j in range(len(x)):
                    use=np.arange(len(x))!=j
                    _,ff,why=math.association(x[use],y[use],c[use])
                    loo.append(ff['rho'] if ff is not None else np.nan)
                finite=np.isfinite(loo)
                row['leave_one_patient_out']=dict(estimable=int(np.sum(finite)),
                    same_direction=int(np.sum(np.sign(np.asarray(loo)[finite])==np.sign(f['rho']))),
                    rho_range=[float(np.nanmin(loo)),float(np.nanmax(loo))] if finite.any() else None)
            else:
                audits.append(independent.audit([clean(dict(label=label,x=x,y=y,C=c,expected=None))]))
            results.append(row)
            for setting,qc,minimum,controls in [
                ('clinical_background_only','human_RNA_QC500_250_MT20',20,'clinical'),
                ('all_retained_RNA','all_retained_authority_cells',20,'full'),
                ('parent_min10','human_RNA_QC500_250_MT20',10,'full'),
                ('parent_min30','human_RNA_QC500_250_MT20',30,'full')]:
                _,rr=select_rows(payload,phase,qc,minimum)
                xx,yy,cc=arrays(rr,state,controls)
                _,ff,why=math.association(xx,yy,cc)
                sensitivity.append(dict(phase=phase,state=state,setting=setting,n=len(rr),status=why,rho=ff['rho'] if ff else None))
    # Missing/failed designs remain explicit slots with P=1 for the fixed family.
    p=[r.get('raw_P',1.) for r in results]
    for r,bh,by in zip(results,math.adjustments(p),math.adjustments(p,by=True)):
        r.update(auxiliary_BH12=float(bh),auxiliary_BY12=float(by))
    pairs=[]
    lookup={(r['phase'],r['state']):r for r in results}
    for ph in ['pre','post']:
        for rs,es,direction in pairs_fixed:
            a,b=lookup[ph,rs],lookup[ph,es]
            assert a['n']==b['n']
            good=all(v['status']=='ESTIMABLE' for v in [a,b])
            pairs.append(dict(phase=ph,R_state=rs,E_state=es,n=a['n'],
                status='ESTIMABLE' if good else 'AT_LEAST_ONE_ENDPOINT_NOT_ESTIMABLE',
                conjunction_P=conjunction(a.get('raw_P',1.),b.get('raw_P',1.)),
                endpoint_clinical_directions_compatible_with_external_covariance=(bool(np.sign(a['rho']*b['rho'])==direction) if good else None)))
    qp=[r['conjunction_P'] for r in pairs]
    for r,bh,by in zip(pairs,math.adjustments(qp),math.adjustments(qp,by=True)):
        r.update(primary_BH6=float(bh),primary_BY6=float(by))
    # A separate implementation checks the complete fixed families, not only positives.
    rcode="suppressPackageStartupMessages(library(jsonlite));a<-fromJSON(file('stdin'));cat(toJSON(lapply(a,function(x)list(BH=p.adjust(x,'BH'),BY=p.adjust(x,'BY'))),digits=16))"
    check=subprocess.run([os.environ['PANCANCER_R_BLOCK_REFERENCE'],'--vanilla','-e',rcode],
        input=json.dumps(dict(marginal=p,conjunction=qp)),capture_output=True,text=True,timeout=60,check=True)
    rc=json.loads(check.stdout)
    for name,ps in [('marginal',p),('conjunction',qp)]:
        for key,by in [('BH',False),('BY',True)]:
            np.testing.assert_allclose(rc[name][key],math.adjustments(ps,by=by),rtol=0,atol=1e-13)
    passed=all(s['passed_finite_alarm'] for r in calibrations for s in r['scenarios'])
    save('private_phase_clinical_coverage.json',sourcecoverage)
    save('private_all12_marginal_clinical_results.json',results)
    save('private_all6_pair_clinical_conjunctions.json',pairs)
    save('private_all48_descriptive_sensitivities.json',sensitivity)
    save('private_independent_R_all_random_streams.json',audits)
    save('private_all_fixed_design_calibration.json',calibrations)
    save('private_terminal.json',dict(status='COMPLETE_CALIBRATION_ACCEPTED' if passed else 'COMPLETE_FINITE_CALIBRATION_ALARM_INFERENCE_NOT_ACCEPTED',
        started_UTC=started,completed_UTC=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-begin,
        fixed_primary_slots=6,fixed_auxiliary_slots=12,descriptive_sensitivity_slots=48,
        all_original_libraries_checked=True,all_real_primary_endpoint_random_streams_R_replayed=True,
        complete_family_BH_BY_R_p_adjust_checked=True,individual_raw_clinical_expression_scores_and_random_arrays_RAM_only=True,
        original_CM_identity_and_old_families_unchanged=True,
        source_previously_explored_not_untouched_validation=True,clinical_prediction_increment_not_tested=True,
        frozen_implementation_SHA256=sha(out/'private_frozen_clinical_implementation_v1.json')))
    print(json.dumps(dict(status='SOURCE_CLINICAL_FAMILY_FINISHED',finite_calibration_accepted=passed)))


if __name__=='__main__': main()
