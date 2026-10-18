"""Independent OLS checks before any real ligand-programme association."""
import json
import os
from datetime import datetime,timezone
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

import NSCLC_native_interface_association_math_v1 as kernel
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def independent(a,b,quality,blocks):
    q=rankdata(quality,axis=0);sd=q.std(axis=0)
    z=np.column_stack([np.ones(len(a)),*[(blocks==s).astype(float) for s in sorted(set(blocks))],(q-q.mean(axis=0))/np.where(sd>0,sd,1)])
    y=rankdata(np.column_stack([a,b]),axis=0);coefficient=np.linalg.lstsq(z,y,rcond=1e-10)[0]
    residual=y-z@coefficient
    return float(np.corrcoef(residual.T)[0,1]),z,y,residual


def main():
    root=Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);resource(root,1048576)
    out=root/'NSCLC_native_interface_association_kernel_check_20261018_v1';assert not out.exists();out.mkdir(mode=0o700)
    (out/'plan_before_synthetic_kernel_checks.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),math_sha256=sha(kernel.__file__),scope='Synthetic positive/negative/alias/zero/missing/small; full513 FL,41 bootstrap,80 LOO independent OLS. No real effects.'))+'\n')
    rng=np.random.default_rng(202610180739);n=80;blocks=np.repeat(['a','b','c','d'],20)
    quality=rng.normal(size=(n,20));quality[:,19]=quality[:,0]
    latent=rng.normal(size=n);a=latent+quality[:,0]+rng.normal(size=n)*.4;b=.9*latent+quality[:,1]+rng.normal(size=n)*.4
    p=kernel.prepare(a,b,quality,blocks);rho,z,y,r=independent(a,b,quality,blocks)
    assert abs(kernel.effect(p)-rho)<1e-12 and rho>0
    assert kernel.effect(kernel.prepare(a,-b,quality,blocks))<0
    generated=np.random.default_rng(192819);count=0;fit=y[:,1]-r[:,1]
    for start in range(0,513,256):
        m=min(256,513-start);indices=np.broadcast_to(np.arange(n),(m,n)).copy()
        for level in sorted(set(blocks)):
            g=np.flatnonzero(blocks==level);indices[:,g]=g[np.argsort(generated.random((m,len(g))),axis=1)]
        for ix in indices:
            null=fit+r[ix,1];residual=null-z@np.linalg.lstsq(z,null,rcond=1e-10)[0]
            value=residual@r[:,0]/(np.linalg.norm(residual)*np.linalg.norm(r[:,0]));count+=abs(value)>=abs(rho)-1e-12
    assert count==kernel.permutation_count(p,513,192819)
    actual=kernel.bootstrap(p,41,810381);generated=np.random.default_rng(810381);reference=[]
    groups=[np.flatnonzero(blocks==s) for s in sorted(set(blocks))]
    for _ in range(41):
        ix=np.concatenate([g[generated.integers(len(g),size=len(g))] for g in groups]);reference.append(independent(a[ix],b[ix],quality[ix],blocks[ix])[0])
    assert np.allclose(actual,reference,atol=1e-12,rtol=1e-12)
    reference=[]
    for i in range(n):
        keep=np.arange(n)!=i;reference.append(independent(a[keep],b[keep],quality[keep],blocks[keep])[0])
    assert np.allclose(kernel.leave_one(p),reference,atol=1e-12,rtol=1e-12)
    failures=[]
    for name,aa,bb,qq,bs in [('small',a[:49],b[:49],quality[:49],blocks[:49]),('constant',np.ones(n),b,quality,blocks),('confounded',quality[:,0],b,quality,blocks),('saturated',a,b,np.eye(n),blocks),('missing',np.full(n,np.nan),b,quality,blocks)]:
        try:kernel.prepare(aa,bb,qq,bs)
        except ValueError:failures.append(name)
    assert len(failures)==5
    # Exact response in a duplicate-rank nuisance span cannot yield correlation.
    group=np.tile(np.arange(20),4);q=np.eye(20)[group];aa=group.astype(float);bb=(group%7).astype(float)
    try:kernel.prepare(aa,bb,q,blocks)
    except ValueError as error:assert str(error)=='INDIVIDUAL_PROGRAMME_RESIDUAL_CONSTANT'
    else:raise AssertionError('Exact zero residual called variable')
    report=dict(status='PASS_INDEPENDENT_OLS_POSITIVE_NEGATIVE_ALIAS_ZERO_MISSING_513FL_41BOOTSTRAP_80LOO',finished_UTC=datetime.now(timezone.utc).isoformat(),source_sha256=sha(__file__),math_sha256=sha(kernel.__file__),no_real_association_or_clinical_effects_loaded=True,finite_numeric_check_not_universal_typeI_or_exchangeability_guarantee=True)
    proof=out/'private_independent_native_interface_kernel_check.json';proof.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],finished_UTC=report['finished_UTC'],audit_sha256=sha(proof))))


if __name__=='__main__':main()
