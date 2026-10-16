"""Independent numerical fixtures for the fixed conditional rank statistic."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

import NSCLC_external_programme_coordination_math_v1 as k
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def independent(a, b, quality, blocks):
    counts = {s: int(np.sum(blocks == s)) for s in set(blocks)}
    keep = np.array([counts[s] >= 2 for s in blocks])
    a, b, quality, blocks = a[keep], b[keep], quality[keep], blocks[keep]
    qr = rankdata(quality, axis=0)
    std = qr.std(axis=0)
    qr = (qr-qr.mean(axis=0))/np.where(std>0, std, 1)
    z = np.column_stack([np.ones(len(a)), *[(blocks == s).astype(float) for s in sorted(set(blocks))], qr])
    responses = rankdata(np.column_stack([a, b]), axis=0)
    coef, _, rank, _ = np.linalg.lstsq(z, responses, rcond=1e-10)
    residuals = responses-z@coef
    rho = np.corrcoef(residuals.T)[0, 1]
    return float(rho), z, responses, residuals, rank


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    resource(root, 1048576)
    out = root/'NSCLC_external_coordination_kernel_check_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    (out/'plan_before_synthetic_numerical_checks.json').write_text(json.dumps(dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), math_sha256=sha(k.__file__), no_real_effects_loaded=True))+'\n')
    rng = np.random.default_rng(202610162021)
    n = 45
    blocks = np.repeat(['a', 'b', 'c'], 15)
    quality = rng.normal(size=(n, 10))
    quality[:, 9] = quality[:, 0]  # Full redundant quality column is retained.
    latent = rng.normal(size=n)
    a = latent + quality[:, 0]*2 + rng.normal(size=n)*.3
    b = .8*latent + quality[:, 1]*2 + rng.normal(size=n)*.3
    p = k.prepare(a, b, quality, blocks)
    expected, z, response, residual, rank = independent(a, b, quality, blocks)
    assert abs(k.effect(p)-expected) < 1e-12 and rank == p['design_rank']
    assert k.effect(p) > 0 and k.effect(k.prepare(a, -b, quality, blocks)) < 0
    # Full Freedman-Lane fitted response + permuted residual, recomputed OLS.
    random_seed = 918271
    draws = 513
    generated = np.random.default_rng(random_seed)
    independent_count = 0
    fit = response[:, 1]-residual[:, 1]
    for start in range(0, draws, 256):
        m = min(256, draws-start)
        take = np.broadcast_to(np.arange(n), (m, n)).copy()
        for level in sorted(set(blocks)):
            ids = np.flatnonzero(blocks == level)
            take[:, ids] = ids[np.argsort(generated.random((m, len(ids))), axis=1)]
        for row in take:
            null_y = fit + residual[row, 1]
            coefficient = np.linalg.lstsq(z, null_y, rcond=1e-10)[0]
            null_residual = null_y-z@coefficient
            rho = residual[:, 0]@null_residual/(np.linalg.norm(residual[:, 0])*np.linalg.norm(null_residual))
            independent_count += abs(rho) >= abs(expected)-1e-12
    assert independent_count == k.permutation_count(p, draws, random_seed)
    actual_boot = k.bootstrap(p, 41, 81037)
    generated = np.random.default_rng(81037)
    groups = [np.flatnonzero(blocks == level) for level in sorted(set(blocks))]
    reference_boot = []
    for _ in range(41):
        take = np.concatenate([ids[generated.integers(len(ids), size=len(ids))] for ids in groups])
        reference_boot.append(independent(a[take], b[take], quality[take], blocks[take])[0])
    assert np.allclose(actual_boot, reference_boot, atol=1e-12, rtol=1e-12)
    actual_loo = k.leave_one(p)
    reference_loo = []
    for i in range(n):
        take = np.arange(n) != i
        reference_loo.append(independent(a[take], b[take], quality[take], blocks[take])[0])
    assert np.allclose(actual_loo, reference_loo, atol=1e-12, rtol=1e-12)
    failures = []
    for name, aa, bb, qq, block in [
        ('small', a[:19], b[:19], quality[:19], blocks[:19]),
        ('constant', np.ones(n), b, quality, blocks),
        ('exact_quality_confounding', quality[:, 0], b, quality, blocks),
        ('saturated', a, b, np.eye(n), blocks)]:
        try:
            k.prepare(aa, bb, qq, block)
        except ValueError:
            failures.append(name)
    assert len(failures) == 4
    more_blocks = blocks.copy()
    more_blocks[-1] = 'x'
    p2 = k.prepare(a, b, quality, more_blocks)
    assert p2['keep'].sum() == 44 and not p2['keep'][-1]
    proof = dict(status='PASS_SYNTHETIC_QR_VS_OLS_POSITIVE_NEGATIVE_FULL_FL_BOOTSTRAP_LOO_AND_UNESTIMABLE_FIXTURES',
                 finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), math_sha256=sha(k.__file__),
                 permutation_values_checked=draws, bootstrap_values_checked=41, LOO_values_checked=45,
                 numerical_correctness_not_universal_exchangeability_or_typeI_guarantee=True,
                 real_effects_or_P_loaded=False)
    path = out/'private_independent_kernel_check.json'
    path.write_text(json.dumps(proof, indent=2)+'\n')
    print(json.dumps(dict(status=proof['status'], finished_UTC=proof['finished_UTC'], audit_sha256=sha(path))))


if __name__ == '__main__':
    main()
