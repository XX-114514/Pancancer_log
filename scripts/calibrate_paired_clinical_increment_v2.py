"""Expanded pre-effect binary-null calibration of conditional and IUT rules."""
import hashlib,json,os
from pathlib import Path
from datetime import datetime,timezone
import numpy as np
from paired_clinical_increment_kernel_v1 import fit,wild_count


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    r=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);v1=r/'paired_clinical_increment_calibration_20261016_v1/private_terminal.json'
    old=json.loads(v1.read_text());assert old['status'].startswith('PASS_')
    kernel=Path(__file__).with_name('paired_clinical_increment_kernel_v1.py');assert sha(kernel)==old['kernel_sha256']
    out=r/'paired_clinical_increment_calibration_20261016_v2';assert not out.exists();out.mkdir(mode=0o700)
    plan=dict(created_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),kernel_sha256=sha(kernel),original_calibration_sha256=sha(v1),
              why='Coarse256 null scenarios included liberal point estimates. Expand all n and both asymmetric nulls BEFORE any real effect, including the actual max-P intersection-union rule. No change to kernel or scientific family.',
              scenarios=['both_null','A_null_B_nonzero','B_null_A_nonzero'],n=[50,100,200],simulations_each=2048,wild_draws_each=1999,
              diagnostic_gate='Every genuinely null conditional test rejection rate AND primary IUT rejection rate <=.075. Diagnostic only, not exact-size guarantee.',seed=2026101612,real_effects_used=False)
    (out/'plan_before_expanded_binary_null_diagnostics.json').write_text(json.dumps(plan,indent=2)+'\n');rng=np.random.default_rng(plan['seed']);diagnostics=[]
    for n in plan['n']:
        for scenario in plan['scenarios']:
            nulls=[0,1] if scenario=='both_null' else [0] if scenario=='A_null_B_nonzero' else [1]
            probabilities=[]
            for k in range(2048):
                c=rng.uniform(-1,1,(n,16));a=rng.uniform(-1,1,n)+.2*c[:,1];b=.3*a+rng.uniform(-1,1,n)
                pr=.5+.1*c[:,0]+(.2*b if scenario=='A_null_B_nonzero' else .2*a if scenario=='B_null_A_nonzero' else 0.)
                assert pr.min()>0 and pr.max()<1;y=rng.binomial(1,pr);d=np.column_stack([np.ones(n),c,a,b]);g=fit(y,d,min_df=10)
                seed=plan['seed']+n*100000+plan['scenarios'].index(scenario)*10000+k*2
                probabilities.append([(1+wild_count(g,-2+j,1999,seed+j))/2000 for j in range(2)])
                if (k+1)%256==0:
                    (out/'private_checkpoint.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(),n=n,scenario=scenario,simulations_done=k+1))+'\n')
            p=np.asarray(probabilities);rates={('A' if j==0 else 'B'):float(np.mean(p[:,j]<.05)) for j in nulls}
            iut=float(np.mean(p.max(1)<.05));diagnostics.append(dict(n=n,scenario=scenario,simulations=2048,null_conditional_rates=rates,IUT_false_rejection_rate=iut))
    passed=all(x['IUT_false_rejection_rate']<=.075 and max(x['null_conditional_rates'].values())<=.075 for x in diagnostics)
    terminal=dict(status='PASS_EXPANDED_BINARY_NULL_CONDITIONAL_AND_IUT_DIAGNOSTICS' if passed else 'FAIL_EXPANDED_BINARY_NULL_DIAGNOSTICS_DO_NOT_RUN_REAL_EFFECTS',
                  finished_UTC=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),kernel_sha256=sha(kernel),diagnostics=diagnostics,real_effects_used=False,exact_small_sample_size_guaranteed=False)
    (out/'private_terminal.json').write_text(json.dumps(terminal,indent=2)+'\n');print(json.dumps(terminal))


if __name__=='__main__':main()
