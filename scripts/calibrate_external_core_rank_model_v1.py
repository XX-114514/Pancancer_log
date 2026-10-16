"""Validate complete-control HC3/wild algebra before external effects."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import statsmodels.api as sm
from external_core_rank_model_v1 import point, wild_statistics, wild_count


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'external_core_rank_model_calibration_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    rng = np.random.default_rng(2026101603)
    checked = 0; maximum = 0.
    for duplicate in [False, True]:
        for trial in range(32):
            c = rng.normal(size=(26, 4))
            if duplicate:
                c = np.column_stack([c, c[:, 0], np.zeros((26, 5))])
            x = rng.lognormal(.3 * c[:, 0], 1); x[:5] = 0
            y = .6 * x + c[:, 0] + rng.normal(size=26)
            g = point(x, y, c)
            d = np.column_stack([g['z'], g['x']])
            fit = sm.OLS(g['y'], d).fit(cov_type='HC3')
            assert np.isclose(g['coefficient'], fit.params[-1], rtol=1e-9, atol=1e-9)
            assert np.isclose(g['HC3_t'], fit.tvalues[-1], rtol=1e-8, atol=1e-8)
            assert int(fit.df_resid) == g['residual_df']
            maximum = max(maximum, abs(g['HC3_t'] - fit.tvalues[-1]))
            signs = np.where(rng.random((26, 32)) < .5, -1., 1.)
            got = wild_statistics(g, signs)
            h0 = np.sum(g['q0'] ** 2, axis=1)
            mu = g['q0'] @ (g['q0'].T @ g['y'])
            e0 = (g['y'] - mu) / np.sqrt(1-h0)
            for j in range(32):
                independent = sm.OLS(mu + e0 * signs[:, j], d).fit(cov_type='HC3').tvalues[-1]
                assert np.isclose(got[j], independent, rtol=1e-8, atol=1e-8)
            checked += 1
    c = rng.normal(size=(26, 10))
    for x, reason in [(np.ones(26), 'CONSTANT'), (c[:, 0].copy(), 'NOT_IDENTIFIABLE')]:
        # Exact collinearity refers to the ranked-X scale used by this model.
        if reason == 'NOT_IDENTIFIABLE':
            from scipy.stats import rankdata
            c[:, 0] = rankdata(x)
        try:
            point(x, rng.normal(size=26), c)
            raise AssertionError('Unidentifiable test was accepted')
        except ValueError as error:
            assert reason in str(error)
    calibration = []
    for heterogeneous in [False, True]:
        probabilities = []
        for trial in range(64):
            c = rng.normal(size=(26, 10)); x = rng.normal(size=26)
            error = rng.normal(size=26)
            if heterogeneous:
                error *= .5 + np.abs(c[:, 0])
            y = c[:, 0] * .5 + error
            g = point(x, y, c)
            probability = (1 + wild_count(g, 1999, 2026101604 + trial + 1000*heterogeneous)) / 2000
            probabilities.append(probability)
        rate = float(np.mean(np.asarray(probabilities) < .05))
        calibration.append(dict(heterogeneous=heterogeneous, simulations=64, wild_draws_each=1999, fraction_p_below05=rate))
        assert rate <= .15, 'Gross null-calibration failure; do not run external effects'
    result = dict(status='PASS_SYNTHETIC_FULL_DUPLICATED_CONSTANT_CONTROLS_IDENTIFIABILITY_HC3_ALL_WILD_ALGEBRA_AND_COARSE_NULL_CALIBRATION',
                  finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), kernel_sha256=sha(Path(__file__).with_name('external_core_rank_model_v1.py')),
                  algebra_scenarios=checked, explicit_wild_checks=checked*32, maximum_t_difference=maximum,
                  coarse_null_diagnostics=calibration, small_sample_exact_size_guaranteed=False, external_data_effects_used=False)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
