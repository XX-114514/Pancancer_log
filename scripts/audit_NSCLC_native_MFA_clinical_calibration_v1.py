"""Replay every synthetic null/power scenario with independent QR numerical math."""
import json
from datetime import datetime,timezone

import numpy as np
import pandas as pd

from NSCLC_native_MFA_clinical_inputs_v1 import read,matched_quota,sha,VIEWS
from NSCLC_native_MFA_clinical_audit_math_v1 import fit,wild


def main():
    root,pre,plan,scope,pairs,arrays=read();src=root/'NSCLC_native_MFA_clinical_calibration_20261016_v1'
    terminal=json.loads((src/'private_terminal.json').read_text())
    for name,item in terminal['files'].items():assert sha(src/name)==item['sha256']
    out=root/'NSCLC_native_MFA_clinical_calibration_audit_20261016_v1';assert not out.exists();matched_quota(root);out.mkdir(mode=0o700)
    observed=pd.read_csv(src/'private_all16_null_and_power_scenarios.tsv',sep='\t');assert len(observed)==16
    rng=np.random.default_rng(plan['calibration_seed']);checked=0
    errors=['normal','t3_unit_variance','heteroskedastic_normal','strong_target_normal']
    for vi,view in enumerate(VIEWS):
        a=arrays[('internal_holdout',view)];d=a['D'];n=len(d);mu=d[:,:-1]@np.linspace(-.2,.2,d.shape[1]-1)
        for ei,error in enumerate(errors):
            count=0
            for k in range(1024):
                noise=rng.standard_t(3,n)/np.sqrt(3) if error=='t3_unit_variance' else rng.normal(size=n)
                if error=='heteroskedastic_normal':noise*=np.exp(.6*a['Zraw'][:,1])
                y=mu+noise+(a['X'] if error=='strong_target_normal' else 0.)
                g=fit(y,d);p=(1+wild(g,1999,plan['calibration_seed']+vi*100000+ei*10000+k))/2000
                count+=p<.05;checked+=1
                if (k+1)%64==0:(out/'private_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),view=view,error_model=error,simulations=k+1,total_checked=checked))+'\n')
            row=observed[(observed['view']==view)&(observed.error_model==error)].iloc[0]
            assert int(row.rejections)==count and row.simulations==1024 and abs(row.rejection_rate-count/1024)<1e-12
        rng.normal(size=n) # production's additional synthetic parity fixture
    passed=bool(observed.loc[observed.error_model!='strong_target_normal','rejection_rate'].le(.075).all())
    assert terminal['status']==('PASS_ALL12_NATIVE_CLINICAL_WORKING_MODEL_NULL_SCENARIOS' if passed else 'FAIL_NULL_CALIBRATION_DO_NOT_COMPUTE_ACTUAL_EFFECTS')
    report=dict(status='PASS_FULL16384_SYNTHETIC_SCENARIOS_AND32751616_WILD_DRAWS_INDEPENDENT_QR_REPLAY',finished_UTC=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),math_sha256=sha(__import__('NSCLC_native_MFA_clinical_audit_math_v1').__file__),
        producer_terminal_sha256=sha(src/'private_terminal.json'),synthetic_datasets_checked=checked,wild_draws_checked=checked*1999,
        original_null_gate_passed=passed,actual_clinical_scores_used=False)
    (out/'private_terminal.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
