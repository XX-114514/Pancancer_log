"""Replay every finite null simulation and all wild exceedance counts via QR."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v1 import wait,read
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design
from NSCLC_immune_RNA_MCP_cohort_audit_math_v1 import prepare,evaluate,wild


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,4*1024**2)
    source=root/'NSCLC_immune_RNA_MCP_cohort_calibration_20261019_v1'
    out=root/'NSCLC_immune_RNA_MCP_cohort_calibration_audit_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_all_null_random_replay.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        scope='All54 fixed null-design scenarios, every simulated data set/HC3beta,t and every complete wild stream exceedance count. Independent QR statistic, regenerate exact errors/signs, no true source association.'),indent=2)+'\n')
    terminal=wait(root,out,'calibrate_NSCLC_immune_RNA_MCP_cohort_v1',source/'private_terminal.json',{'PASS_ALL_FIXED_DESIGN_NULL_DIAGNOSTICS','FAIL_NULL_DIAGNOSTIC_DO_NOT_COMPUTE_TRUE_EFFECTS'})
    for name,item in terminal['files'].items():assert sha(source/name)==item['sha256']
    _,_,_,scope,arrays=read(root);statistics=pd.read_csv(source/'private_all_null_simulation_statistics.tsv',sep='\t').set_index(['slot','component','error_model','simulation'])
    summary=pd.read_csv(source/'private_all54_null_design_scenarios.tsv',sep='\t').set_index(['slot','component','error_model'])
    assert len(summary)==54;rng=np.random.default_rng(2026101902);models=0;values_regenerated=0;maxerror=0.;passed=True
    for index,item in enumerate(scope.itertuples(index=False)):
        if item.status!='ESTIMABLE':
            for error in ['normal','t3_unit_variance','age_heteroskedastic_normal']:assert summary.loc[(item.slot,item.component,error),'status']=='ORIGINAL_DESIGN_UNAVAILABLE'
            continue
        a=arrays[(item.slot,item.component)];original,status=design(a['x'],a['y'],a['C'],a['names'],str(a['kind']));reference,status2=prepare(a['x'],a['y'],a['C'],a['names'],str(a['kind']));assert status==status2=='ESTIMABLE'
        # SVD basis here reproduces the producer's declared synthetic data mean;
        # estimation and all null statistics below use independent pivoted QR.
        D=original['D'];mu=D[:,:-1]@np.linspace(-.2,.2,D.shape[1]-1);n=len(D)
        age=a['C'][:,list(a['names']).index('age_numeric')];age=(age-age.mean())/age.std(ddof=0)
        for error_id,error in enumerate(['normal','t3_unit_variance','age_heteroskedastic_normal']):
            reject=0
            for simulation in range(256):
                e=rng.standard_t(3,n)/np.sqrt(3) if error=='t3_unit_variance' else rng.normal(size=n)
                if error=='age_heteroskedastic_normal':e*=np.exp(.6*age)
                y=mu+e;random_seed=2026101902+index*100000+error_id*10000+simulation
                g,values=wild(reference,999,random_seed,y=y);count=int((np.abs(values)>=abs(g['t'])-1e-12).sum());row=statistics.loc[(item.slot,item.component,error,simulation)]
                assert row.wild_seed==random_seed and row.wild_exceedances==count
                error_numeric=max(abs(row.beta-g['beta']),abs(row.HC3_t-g['t']));assert error_numeric<1e-8;maxerror=max(maxerror,error_numeric)
                reject+=(1+count)/1000<.05;models+=1;values_regenerated+=999
            r=summary.loc[(item.slot,item.component,error)];assert r.rejections==reject and r.simulations==256 and abs(r.rejection_rate-reject/256)<1e-12
            passed=passed and reject/256<=.10
            (out/'private_progress.json').write_text(json.dumps(dict(stage='INDEPENDENT_NULL_ONLY_REPLAY',completed_data_sets=models,wild_statistics_regenerated=values_regenerated,model=index+1,error_model=error,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
    assert models==len(statistics)==terminal['simulated_data_sets']
    assert (terminal['status']=='PASS_ALL_FIXED_DESIGN_NULL_DIAGNOSTICS')==passed
    status='ACCEPTED_ALL_NULL_CALIBRATION_VALUES_AND_PASS_GATE' if passed else 'ACCEPTED_NULL_DIAGNOSTIC_FAILURE_NO_TRUE_EFFECTS'
    t=dict(status=status,completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),producer_terminal_sha256=sha(source/'private_terminal.json'),
        all_null_scenarios=54,simulated_data_sets_checked=models,complete_wild_statistics_regenerated=values_regenerated,max_statistic_error=maxerror,actual_source_association_effects_computed=False)
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
