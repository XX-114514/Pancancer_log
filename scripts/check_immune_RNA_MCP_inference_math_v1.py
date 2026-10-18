"""Synthetic numerical gate: signed effects, failures and independent random replay."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from review_immune_RNA_MCP_seed_stability_v1 import centred_rank_correlation, signed_weight_cosine
from test_immune_RNA_MCP_source_coordination_v1 import fit_rank_partial, resample, adjust_family


def independent_effect(x, y, controls):
    # Independently solve the ranked-outcome nuisance regressions without QR.
    rank_c = np.column_stack([rankdata(controls[:, j]) for j in range(controls.shape[1])])
    rank_c = (rank_c - rank_c.mean(axis=0)) / rank_c.std(axis=0, ddof=1)
    design = np.column_stack([np.ones(len(x)), rank_c])
    assert np.linalg.matrix_rank(design) == design.shape[1]
    rx = rankdata(x) - design @ np.linalg.lstsq(design, rankdata(x), rcond=None)[0]
    ry = rankdata(y) - design @ np.linalg.lstsq(design, rankdata(y), rcond=None)[0]
    return float(np.corrcoef(rx, ry)[0, 1])


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); budget = resource(root, 1024**2)
    output = root / 'immune_RNA_MCP_inference_math_check_20261018_v1'
    assert not output.exists(); output.mkdir(mode=0o700)
    sources = {n: sha(Path(__file__).with_name(n)) for n in [
        'test_immune_RNA_MCP_source_coordination_v1.py', 'review_immune_RNA_MCP_seed_stability_v1.py']}
    (output / 'plan_before_synthetic_math_check.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        production_sources=sources, budget=budget,
        scope='Only generated synthetic values. Positive/negative/unavailable, all19999 independent OLS permutation values, all2000 independent OLS bootstrap effects, BH/BY reference and source-centred matching.'), indent=2) + '\n')
    rng = np.random.default_rng(621018)
    controls = rng.normal(size=(40, 8)); x = rng.normal(size=40)
    y = .8*x + .4*rng.normal(size=40) + .5*controls[:, 0]
    fit, reason = fit_rank_partial(x, y, controls)
    reference = independent_effect(x, y, controls)
    assert reason == 'ESTIMABLE' and reference > .5 and abs(reference-fit['effect']) < 1e-12
    negative, _ = fit_rank_partial(x, -y, controls)
    assert abs(negative['effect'] + reference) < 1e-12
    result, draws = resample(x, y, controls, fit, 0, 2026101803)
    design = np.c_[np.ones(40), np.apply_along_axis(rankdata, 0, controls)]
    xr = rankdata(x); yr = rankdata(y)
    rx = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    fitted_y = design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    ry = yr - fitted_y
    permuted = []
    for take in draws['permutation_indices']:
        yp = fitted_y + ry[take]
        residual = yp - design @ np.linalg.lstsq(design, yp, rcond=None)[0]
        permuted.append(np.corrcoef(rx, residual)[0, 1])
    perm_error = float(np.max(np.abs(np.asarray(permuted)-draws['permutation_effects'])))
    assert perm_error < 1e-10
    p = (1 + np.count_nonzero(np.abs(permuted) >= abs(reference)-1e-12)) / 20000
    assert p == result['permutation_p']
    bootstrap = np.asarray([independent_effect(x[t], y[t], controls[t]) for t in draws['bootstrap_indices']])
    assert result['valid_bootstrap'] == 2000
    boot_error = float(np.max(np.abs(bootstrap-draws['bootstrap_effects'])))
    assert boot_error < 1e-10
    assert np.allclose(np.quantile(bootstrap, [.025, .975]), [result['lower'], result['upper']], atol=1e-10)
    assert fit_rank_partial(x[:9], y[:9], controls[:9])[1] == 'FEWER_THAN20_COMPLETE_UNITS'
    assert fit_rank_partial(x, y, np.c_[controls[:, :7], controls[:, 0]])[1] == 'COLLINEAR_REQUIRED_DESIGN_NO_AUTOMATIC_DELETION'
    assert fit_rank_partial(x, y, np.c_[controls[:, :7], np.zeros(40)])[1] == 'CONSTANT_REQUIRED_COVARIATE_NO_AUTOMATIC_DELETION'
    expected_bh = np.array([.006, .06, .2, 1, 1, 1])
    test_p = [.001, .02, .1, 1, 1, 1]
    assert np.allclose(adjust_family(test_p), expected_bh)
    assert np.allclose(adjust_family(test_p, by=True), np.minimum(expected_bh * np.sum(1/np.arange(1, 7)), 1))
    a = pd.DataFrame(dict(sample=['s%d' % i for i in range(38)],
        dataset=['GSE286827']*28 + ['GSE294109']*10, score=rng.normal(size=38)))
    b = a.copy(); b['score'] = a.score.to_numpy() + np.r_[np.zeros(28), np.ones(10)*100]
    rho, _, n = centred_rank_correlation(a, b)
    assert abs(rho-1) < 1e-12 and n == 38
    left = pd.DataFrame(dict(gene=['A', 'B'], weight=[1., 2.]))
    right = pd.DataFrame(dict(gene=['A', 'B', 'C'], weight=[-1., -2., 0.]))
    assert abs(signed_weight_cosine(left, right)+1) < 1e-12
    record = dict(status='PASS_SYNTHETIC_SIGNED_UNAVAILABLE_ALL_RANDOM_REPLAY',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        production_sources=sources, independent_permutation_values=19999,
        independent_bootstrap_values=2000, maximum_permutation_difference=perm_error,
        maximum_bootstrap_difference=boot_error,
        real_expression_learned_weights_or_clinical_effects_read=False)
    (output / 'private_terminal.json').write_text(json.dumps(record, indent=2) + '\n')
    print('SYNTHETIC_GATE_PASS_NO_REAL_SCIENTIFIC_RESULT')


if __name__ == '__main__':
    main()
