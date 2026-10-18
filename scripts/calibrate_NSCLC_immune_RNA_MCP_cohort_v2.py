"""Finite null working-model gate on all eligible fixed designs, no true fits."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait,read
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design,wild_values
from HNSCC_author_state_clinical_kernel_v2 import fit


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    out=root/'NSCLC_immune_RNA_MCP_cohort_calibration_20261019_v2';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_real_design_null_calibration.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        scope='All18 fixed model positions, each eligible design3 null scenarios x256 data x999 HC3 wild. No actual association estimates.',
        seed=2026101902,gate='Every available null rejection rate(alpha.05)<=.10; fail stops true effects. Working-model anomaly diagnostic, not guarantee under arbitrary misspecification.',
        errors=['normal','t3_unit_variance','age_heteroskedastic_normal'],
        synthetic_mean='Original nuisance basis times linspace(-.2,.2); target coefficient0. Heteroskedastic scale exp(.6*standardized original age).',
        random_storage='Keep every simulated HC3 beta/t and wild exceedance count plus exact seeds. Independent audit regenerates every wild stream; complete actual-family numerical streams stored separately.'),indent=2)+'\n')
    inputaudit=root/'NSCLC_immune_RNA_MCP_cohort_inputs_audit_20261019_v2/private_terminal.json'
    wait(root,out,'audit_NSCLC_immune_RNA_MCP_cohort_inputs_v2',inputaudit,{'ACCEPTED_ALL6_BUNDLES_18_INPUT_DESIGNS_AND_ALL_EXCLUSIONS'})
    mathgate=root/'NSCLC_immune_RNA_MCP_cohort_math_gate_20261019_v1/private_terminal.json'
    assert sha(mathgate)=='4ee96898bfbd355f90032a34faa9477d072afc93c22c6d330627c1a70c2dadbf'
    _,_,_,scope,arrays=read(root);rng=np.random.default_rng(2026101902);rows=[];summary=[]
    for index,item in enumerate(scope.itertuples(index=False)):
        a=arrays[(item.slot,item.component)]
        if item.status!='ESTIMABLE':
            for error in ['normal','t3_unit_variance','age_heteroskedastic_normal']:summary.append(dict(slot=item.slot,component=item.component,error_model=error,status='ORIGINAL_DESIGN_UNAVAILABLE',simulations=0,rejections=0,rejection_rate=np.nan))
            continue
        m,status=design(a['x'],a['y'],a['C'],a['names'],str(a['kind']));assert status=='ESTIMABLE'
        D=m['D'];n=len(D);mu=D[:,:-1]@np.linspace(-.2,.2,D.shape[1]-1)
        age=a['C'][:,list(a['names']).index('age_numeric')];age=(age-age.mean())/age.std(ddof=0)
        for error_id,error in enumerate(['normal','t3_unit_variance','age_heteroskedastic_normal']):
            rejected=0
            for simulation in range(256):
                noise=rng.standard_t(3,n)/np.sqrt(3) if error=='t3_unit_variance' else rng.normal(size=n)
                if error=='age_heteroskedastic_normal':noise*=np.exp(.6*age)
                y=mu+noise;g=fit(y,D);random_seed=2026101902+index*100000+error_id*10000+simulation
                values=wild_values(g,999,random_seed);count=int((np.abs(values)>=abs(g['t'])-1e-12).sum());p=(1+count)/1000
                rejected+=p<.05
                rows.append(dict(slot=item.slot,component=item.component,error_model=error,simulation=simulation,wild_seed=random_seed,beta=g['beta'],HC3_t=g['t'],wild_exceedances=count))
                if (simulation+1)%64==0:(out/'private_progress.json').write_text(json.dumps(dict(stage='SYNTHETIC_NULL_ONLY',model=index+1,total_models=18,error_model=error,simulations=simulation+1,observed_UTC=datetime.now(timezone.utc).isoformat()))+'\n')
            summary.append(dict(slot=item.slot,component=item.component,error_model=error,status='CALIBRATED',simulations=256,rejections=int(rejected),rejection_rate=float(rejected/256)))
    pd.DataFrame(rows).to_csv(out/'private_all_null_simulation_statistics.tsv',sep='\t',index=False)
    table=pd.DataFrame(summary);assert len(table)==54;table.to_csv(out/'private_all54_null_design_scenarios.tsv',sep='\t',index=False)
    tested=table[table.status.eq('CALIBRATED')];passed=bool(len(tested)>0 and tested.rejection_rate.le(.10).all())
    t=dict(status='PASS_ALL_FIXED_DESIGN_NULL_DIAGNOSTICS' if passed else 'FAIL_NULL_DIAGNOSTIC_DO_NOT_COMPUTE_TRUE_EFFECTS',
        completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_acceptance_sha256=sha(inputaudit),
        simulated_data_sets=len(rows),available_scenarios=len(tested),all_scenarios=54,actual_association_effects_computed=False,
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(out.iterdir()) if p.is_file()})
    assert sum(v['bytes'] for v in t['files'].values())<8*1024**2
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
