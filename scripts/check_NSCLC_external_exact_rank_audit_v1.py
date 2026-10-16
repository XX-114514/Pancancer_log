"""Known failed draw and exact-span certificates; original output unchanged."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np

import NSCLC_external_programme_coordination_math_v1 as original
import NSCLC_external_exact_rank_guards_v1 as guard
import audit_NSCLC_external_programme_coordination_v2 as corrected
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from test_NSCLC_external_programme_coordination_v1 import read_inputs,QCOLS


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);resource(root,1048576)
    out=root/'NSCLC_external_exact_rank_audit_check_20261016_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_exact_rank_audit_checks.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),audit_v2_sha256=sha(corrected.__file__),exact_guard_sha256=sha(guard.__file__),scope='Observed failed draw only plus independent synthetic zero/nonzero and prime-rank fallback checks; no source identity/effect or clinical output.'))+'\n')
    _,_,_,slots,cases=read_inputs(root)
    import pandas as pd
    scope=pd.read_csv(root/'NSCLC_external_programme_coordination_20261016_v1/private_all15600_primary_test_eligibility.tsv',sep='\t')
    item=scope[scope.test_index.eq(62)].iloc[0];meta=slots.set_index('slot').loc[item.slot]
    g=cases.iloc[int(meta.row_start):int(meta.row_stop)]
    if item.family=='reviewed_patients':g=g[g.reviewed_patient]
    p=original.prepare(g.raw_A.to_numpy(),g.raw_B.to_numpy(),g[QCOLS].to_numpy(),g.block_key.to_numpy(str))
    rng=np.random.default_rng(original.seed(item.record,'bootstrap'));groups=[np.flatnonzero(p['blocks']==s) for s in sorted(set(p['blocks']))]
    for _ in range(734):take=np.concatenate([ids[rng.integers(len(ids),size=len(ids))] for ids in groups])
    args=(p['a'][take],p['b'][take],p['quality'][take],p['blocks'][take])
    for fn in [original.prepare,corrected.calculate]:
        try:fn(*args)
        except ValueError as e:assert str(e)=='INDIVIDUAL_PROGRAMME_RESIDUAL_CONSTANT'
        else:raise AssertionError('Known exact-zero bootstrap draw was falsely called variable')
    # Full row rank after duplicate compression proves exact fit of any response.
    z=np.eye(6,dtype=int)[np.r_[np.arange(6),np.arange(6)]]
    y=np.column_stack([np.arange(6),np.arange(6)[::-1]])[np.r_[np.arange(6),np.arange(6)]]
    assert guard.any_exactly_fitted(z,y)
    # A rank drop modulo one prime is not proof of rational dependence.
    prime=1000003;z=np.array([[1,0],[0,prime],[1,prime]],np.int64)
    y=np.array([[0,0],[1,2],[2,4]],np.int64)
    assert guard.modular_rank(z)==1 and guard.rational_rank(z)==2
    assert not guard.any_exactly_fitted(z,y)
    assert guard.any_exactly_fitted(z,np.column_stack([z[:,0],y[:,0]]))
    rng=np.random.default_rng(202610162105);n=45;quality=rng.normal(size=(n,10));blocks=np.repeat(['a','b','c'],15)
    a=rng.normal(size=n)+quality[:,0];b=.5*a+rng.normal(size=n)
    p=original.prepare(a,b,quality,blocks);other=corrected.calculate(a,b,quality,blocks)
    assert abs(original.effect(p)-other['rho'])<1e-12 and p['design_rank']==other['rank']
    proof=out/'private_independent_exact_rank_audit_check.json'
    report=dict(status='PASS_KNOWN_EXACT_ZERO_DRAW_SYNTHETIC_FULL_RANK_NONFIT_MODULAR_FALLBACK_AND_NORMAL_OLS',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),audit_v2_sha256=sha(corrected.__file__),exact_guard_sha256=sha(guard.__file__),original_production_changed=False,scientific_thresholds_or_families_changed=False)
    proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
