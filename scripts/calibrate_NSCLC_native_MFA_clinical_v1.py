"""All4 fixed holdout designs,12 null scenarios and4 power diagnostics."""
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

from NSCLC_native_MFA_clinical_inputs_v1 import read,matched_quota,VIEWS,sha
from HNSCC_author_state_clinical_kernel_v2 import fit,wild_count


def main():
    root,pre,plan,scope,pairs,arrays=read();proof=pre/'private_independent_clinical_input_audit.json'
    assert json.loads(proof.read_text())['status'].startswith('PASS_')
    out=root/'NSCLC_native_MFA_clinical_calibration_20261016_v1';assert not out.exists();matched_quota(root);out.mkdir(mode=0o700)
    contract=dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_audit_sha256=sha(proof),
        original_plan_sha256=sha(pre/'plan_before_any_native_factor_clinical_effect.json'),
        simulations_each=1024,wild_each=1999,null_scenarios=12,power_scenarios=4,
        gate='Every12 null rejection rates at alpha.05<=.075; never infer actual effects if a gate fails.',
        error_models=['normal','t3_unit_variance','heteroskedastic_normal'],
        nuisance_mean='Fixed complete SVD nuisance basis times coefficients linspace(-.2,.2); target response term is0 under null.',
        heteroskedasticity='Normal noise multiplied by exp(.6*original standardized age); same fixed clinical/design covariates.',
        independent_error_seed=plan['calibration_seed'],actual_clinical_scores_used=False,
        caveat='Working-model diagnostic, not exact observational inference or guaranteed typeI error under all misspecification.')
    (out/'plan_before_any_true_native_clinical_effect.json').write_text(json.dumps(contract,indent=2)+'\n')
    rng=np.random.default_rng(plan['calibration_seed']);rows=[]
    for vi,view in enumerate(VIEWS):
        a=arrays[('internal_holdout',view)];d=a['D'];n=len(d);mu=d[:,:-1]@np.linspace(-.2,.2,d.shape[1]-1)
        for ei,error in enumerate(contract['error_models']+['strong_target_normal']):
            count=0
            for k in range(1024):
                noise=rng.standard_t(3,n)/np.sqrt(3) if error=='t3_unit_variance' else rng.normal(size=n)
                if error=='heteroskedastic_normal':noise*=np.exp(.6*a['Zraw'][:,1])
                y=mu+noise+(a['X'] if error=='strong_target_normal' else 0.)
                g=fit(y,d);pv=(1+wild_count(g,1999,plan['calibration_seed']+vi*100000+ei*10000+k))/2000
                count+=pv<.05
                if (k+1)%64==0:(out/'private_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),view=view,error_model=error,simulations=k+1))+'\n')
            rows.append(dict(view=view,error_model=error,simulations=1024,rejections=int(count),rejection_rate=float(count/1024)))
        test=mu+.5*a['X']+rng.normal(size=n);g=fit(test,d);reference=sm.OLS(test,d).fit(cov_type='HC3')
        assert np.allclose([g['beta'],g['t']],[reference.params[-1],reference.tvalues[-1]],atol=1e-9,rtol=1e-9)
    table=pd.DataFrame(rows);table.to_csv(out/'private_all16_null_and_power_scenarios.tsv',sep='\t',index=False)
    passed=bool(table.loc[table.error_model!='strong_target_normal','rejection_rate'].le(.075).all())
    terminal=dict(status='PASS_ALL12_NATIVE_CLINICAL_WORKING_MODEL_NULL_SCENARIOS' if passed else 'FAIL_NULL_CALIBRATION_DO_NOT_COMPUTE_ACTUAL_EFFECTS',
        finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),input_audit_sha256=sha(proof),
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()},actual_clinical_scores_used=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n')
    print(json.dumps(dict(status=terminal['status'],finished_UTC=terminal['finished_UTC'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
