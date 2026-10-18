"""Complete synthetic positive/negative/unidentifiable and random-stream gate."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import statsmodels.api as sm

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design,estimate,wild_values,bootstrap,takes,seed
from NSCLC_immune_RNA_MCP_cohort_audit_math_v1 import prepare,evaluate,wild,boot


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);guard=resource(root,8*1024**2)
    out=root/'NSCLC_immune_RNA_MCP_cohort_math_gate_20261019_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_complete_synthetic_gate.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),resource_guard=guard,
        scope='Pure synthetic positive/negative clinical and rank coordination, unidentifiable/constant/lowN cases. Full19999 wild and2000 block-bootstrap values in each available scenario. SVD vs pandas/pivotedQR and statsmodelsHC3. No real scores/outcomes read.'),indent=2)+'\n')
    rng=np.random.default_rng(2026101903);n=120;age=rng.normal(size=n);male=rng.integers(0,2,n);chemo=rng.integers(0,2,n);agents=np.arange(n)%3
    C=np.c_[age,male,chemo,1-chemo,np.zeros(n),np.eye(3)[agents]]
    names=np.array(['age','male','chemo','targeted','zero','agent::A','agent::B','agent::C'])
    response=rng.integers(0,2,n).astype(float);latent=.3*age+rng.normal(size=n);noise=rng.normal(size=n)
    scenarios=[('positive_clinical',response,.9*response+.3*age+noise,'clinical'),('negative_clinical',response,-.9*response+.3*age+noise,'clinical'),
        ('positive_coordination',latent,.8*latent+.3*age+noise,'coordination'),('negative_coordination',latent,-.8*latent+.3*age+noise,'coordination')]
    results=[];maxwild=0.;maxboot=0.
    for slot,(name,x,y,kind) in enumerate(scenarios):
        a,status=design(x,y,C,names,kind);b,status2=prepare(x,y,C,names,kind);assert status==status2=='ESTIMABLE'
        g,effect=estimate(a);other=evaluate(b);reference=sm.OLS(a['y'],np.c_[a['Zraw'],a['x']]).fit(cov_type='HC3')
        assert np.allclose([g['beta'],g['t']],[other['beta'],other['t']],atol=1e-9,rtol=0)
        assert np.allclose([g['beta'],g['t']],[reference.params[-1],reference.tvalues[-1]],atol=1e-9,rtol=0)
        assert abs(effect-other['effect'])<1e-9 and effect*(1 if name.startswith('positive') else -1)>0
        random_seed=seed(slot,'synthetic','wild');values=wild_values(g,19999,random_seed);_,alternate=wild(b,19999,random_seed)
        error=float(np.max(np.abs(values-alternate)));assert error<1e-8;maxwild=max(maxwild,error)
        strata=agents.astype(str);indices=takes(strata,2000,seed(slot,'synthetic','bootstrap'))
        v,reasons=bootstrap(x,y,C,names,kind,indices,a['nuisance_rank']);w,reasons2=boot(x,y,C,names,kind,indices,b['nuisance_rank'])
        assert reasons.tolist()==reasons2.tolist() and np.allclose(v,w,atol=1e-8,rtol=0,equal_nan=True)
        valid=np.isfinite(v);error=float(np.max(np.abs(v[valid]-w[valid])));maxboot=max(maxboot,error)
        np.savez_compressed(out/('private_synthetic_%d_full_stream.npz'%slot),wild=values,bootstrap=v,indices=indices,reasons=reasons)
        results.append(dict(scenario=name,status='PASS_COMPLETE_INDEPENDENT_VALUES',wild_values=19999,bootstrap_values=2000,valid_bootstrap=int(valid.sum())))
    for label,x,y,cc,nn,expected in [('unidentifiable',male,noise,C,names,'TARGET_NOT_IDENTIFIABLE'),('constant',np.ones(n),noise,C,names,'CONSTANT_EXPOSURE_OR_READOUT')]:
        a,status=design(x,y,cc,nn,'clinical');b,status2=prepare(x,y,cc,nn,'clinical');assert a is b is None and status==status2==expected
    a,status=design(response[:40],noise[:40],C[:40],names,'clinical');assert a is None and status=='FEWER_THAN50_AUTHOR_CASES'
    t=dict(status='ACCEPTED_COMPLETE_SYNTHETIC_COHORT_MATH_AND_RANDOM_STREAMS',completed_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),
        scenarios=results,maximum_wild_error=maxwild,maximum_bootstrap_error=maxboot,real_data_effects_read=False,
        files={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    (out/'private_terminal.json').write_text(json.dumps(t,indent=2)+'\n');print(json.dumps(dict(status=t['status'],terminal_sha256=sha(out/'private_terminal.json'))))


if __name__=='__main__':main()
