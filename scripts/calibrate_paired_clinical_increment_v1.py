"""Synthetic binary outcomes only; freeze checks before real clinical effects."""
import hashlib,json,os
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,statsmodels.api as sm
from paired_clinical_increment_kernel_v1 import fit,wild_statistics,wild_count


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);out=root/'paired_clinical_increment_calibration_20261016_v1'
    assert not out.exists();out.mkdir(mode=0o700);rng=np.random.default_rng(2026101610)
    maximum=0.;checks=0
    for n in [50,100,200]:
        for trial in range(16):
            c=rng.uniform(-1,1,(n,16));a=rng.uniform(-1,1,n)+.2*c[:,1];b=.3*a+rng.uniform(-1,1,n)
            probability=.5+.1*c[:,0]+.2*b;assert probability.min()>0 and probability.max()<1
            y=rng.binomial(1,probability);d=np.column_stack([np.ones(n),c,a,b]);g=fit(y,d,min_df=10)
            ref=sm.OLS(y,d).fit(cov_type='HC3')
            assert np.allclose(g['beta'],ref.params,atol=1e-10) and np.allclose(g['t'],ref.tvalues,atol=1e-9)
            maximum=max(maximum,float(np.max(np.abs(g['t']-ref.tvalues))))
            signs=np.where(rng.random((n,16))<.5,-1.,1.)
            for target in [-2,-1]:
                got=wild_statistics(g,target,signs);z=np.delete(d,target,axis=1);null=sm.OLS(y,z).fit()
                h0=sm.OLS(y,z).fit().get_influence().hat_matrix_diag
                for k in range(16):
                    yy=null.fittedvalues+null.resid/np.sqrt(1-h0)*signs[:,k]
                    expected=sm.OLS(yy,d).fit(cov_type='HC3').tvalues[target]
                    assert np.isclose(got[k],expected,atol=1e-8,rtol=1e-8);checks+=1
    diag=[]
    for n in [50,100,200]:
        for scenario in ['both_null','A_null_B_nonzero']:
            ps=[]
            for trial in range(256):
                c=rng.uniform(-1,1,(n,16));a=rng.uniform(-1,1,n)+.2*c[:,1];b=.3*a+rng.uniform(-1,1,n)
                pr=.5+.1*c[:,0]+(.2*b if scenario=='A_null_B_nonzero' else 0.)
                assert pr.min()>0 and pr.max()<1;y=rng.binomial(1,pr);d=np.column_stack([np.ones(n),c,a,b]);g=fit(y,d,min_df=10)
                ps.append((1+wild_count(g,-2,1999,2026101610+n*1000+trial))/2000)
            rate=float(np.mean(np.asarray(ps)<.05));diag.append(dict(n=n,scenario=scenario,simulations=256,wild_draws=1999,null_A_rejection_rate=rate))
            assert rate<=.10,'Gross binary-null diagnostic failure; no real inference'
    try:fit(y,np.column_stack([d,d[:,1]]));raise AssertionError('Duplicate control accepted')
    except ValueError:pass
    result=dict(status='PASS_BINARY_SYNTHETIC_HC3_EXPLICIT_WILD_AND_COARSE_NULL_DIAGNOSTICS_NOT_EXACT_SIZE_GUARANTEE',finished_UTC=datetime.now(timezone.utc).isoformat(),
                code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),kernel_sha256=hashlib.sha256(Path(__file__).with_name('paired_clinical_increment_kernel_v1.py').read_bytes()).hexdigest(),
                explicit_wild_checks=checks,maximum_t_difference=maximum,diagnostics=diag,real_data_effects_used=False)
    (out/'private_terminal.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
