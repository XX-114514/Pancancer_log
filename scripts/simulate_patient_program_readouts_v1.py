"""Artificial design experiment; no file/data input and no biological conclusions.

Compare mean and distributional programme summaries with patient permutations.
The five-quantile statistic is inspired by Pelka, not an exact reproduction:
we use an absolute-statistic Monte Carlo test with +1 correction throughout.
"""
import argparse
import hashlib
import itertools
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

QUANTILES = np.array([.25, .5, .75, .95, .99])
METHODS = ['mean', 'q75', 'q95', 'q99', 'five_quantile_fisher']


def summaries(cells):
    return np.column_stack([cells.mean(axis=1), np.quantile(cells, QUANTILES, axis=1).T])


def unit_columns(a):
    a = np.asarray(a, dtype=float)
    if a.ndim != 2 or not np.isfinite(a).all():
        raise ValueError('finite two-dimensional input required')
    centered = a - a.mean(axis=0)
    norms = np.linalg.norm(centered, axis=0)
    if np.any(norms < 1e-12):
        raise ValueError('constant summary: do not silently omit a quantile')
    return centered / norms


def statistics(x, y, permutations):
    a, b = unit_columns(x), unit_columns(y)
    r = np.einsum('nk,bnk->bk', a, b[permutations])
    return np.column_stack([r[:, 0], r[:, 3:6],
                            np.arctanh(np.clip(r[:, 1:6], -1+1e-12, 1-1e-12)).mean(axis=1)])


def perm_indices(rng, n, count, blocks=None):
    if blocks is None:
        return np.argsort(rng.random((count, n)), axis=1)
    out = np.tile(np.arange(n), (count, 1))
    for block in np.unique(blocks):
        ix = np.flatnonzero(blocks == block)
        out[:, ix] = ix[np.argsort(rng.random((count, len(ix))), axis=1)]
    return out


def pvalues(x, y, perms):
    obs = statistics(x, y, np.arange(len(x))[None, :])[0]
    null = statistics(x, y, perms)
    return (1 + (np.abs(null) >= np.abs(obs)-1e-12).sum(axis=0)) / (1 + len(perms))


def check_numerics():
    rng = np.random.default_rng(24001)
    x, y = rng.normal(size=(6, 6)), rng.normal(size=(6, 6))
    permutations = np.array(list(itertools.permutations(range(6))))
    got = statistics(x, y, permutations)
    expected = []
    for ix in permutations:
        r = np.array([np.corrcoef(x[:, j], y[ix, j])[0, 1] for j in range(6)])
        expected.append([r[0], *r[3:6], np.mean(np.arctanh(r[1:6]))])
    error = float(np.max(np.abs(got - np.array(expected))))
    assert error < 1e-12
    groups = np.repeat(np.arange(2), 3)
    pp = perm_indices(rng, 6, 100, groups)
    assert np.array_equal(groups[pp], np.tile(groups, (100, 1)))
    np.testing.assert_allclose(pvalues(x, y, pp), pvalues(y, x, np.argsort(pp, axis=1)))
    for bad in [np.ones((6, 6)), np.full((6, 6), np.nan)]:
        try:
            unit_columns(bad)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid summary accepted')
    tied = np.round(x)
    assert np.isfinite(statistics(tied, y, pp)).all()
    return dict(exhaustive_permutations=720, max_statistic_error=error,
                block_preservation=True, swap_inverse_parity=True,
                constant_and_missing_rejected=True, tied_scores_supported=True)


def wilson(k, n):
    z = 1.959963984540054
    p = k / n
    center = (p + z*z/(2*n)) / (1 + z*z/n)
    half = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / (1 + z*z/n)
    return [float(center-half), float(center+half)]


def simulate(repetitions, permutations, seed):
    rng = np.random.default_rng(seed)
    n = 32
    groups = np.repeat(np.arange(2), n//2)
    out = []
    for scenario in ['within_source_null', 'source_shift_only', 'whole_population_signal', 'rare_5percent_signal']:
        for ncell in [50, 200, 1000]:
            rejects = np.zeros((2, len(METHODS)), dtype=int)
            for _ in range(repetitions):
                # Independent patient offsets create realistic within-type heterogeneity.
                x = rng.normal(size=(n, ncell)) + rng.normal(0, .5, (n, 1))
                y = rng.normal(size=(n, ncell)) + rng.normal(0, .5, (n, 1))
                if scenario == 'source_shift_only':
                    x += 3 * (groups[:, None]-.5)
                    y += 3 * (groups[:, None]-.5)
                elif scenario == 'whole_population_signal':
                    signal = rng.normal(0, .55, (n, 1))
                    x += signal
                    y += signal
                elif scenario == 'rare_5percent_signal':
                    signal = np.exp(rng.normal(1.5, .5, (n, 1)))
                    x += (rng.random((n, ncell)) < .05) * signal
                    y += (rng.random((n, ncell)) < .05) * signal
                sx, sy = summaries(x), summaries(y)
                for i, blocks in enumerate([None, groups]):
                    p = pvalues(sx, sy, perm_indices(rng, n, permutations, blocks))
                    rejects[i] += p <= .05
            for i, scheme in enumerate(['unrestricted_patient', 'within_source_patient']):
                for j, method in enumerate(METHODS):
                    k = int(rejects[i, j])
                    out.append(dict(scenario=scenario, patients=n, cells_per_type_patient=ncell,
                                    permutation_scheme=scheme, readout=method,
                                    rejections=k, repetitions=repetitions, rejection_fraction=k/repetitions,
                                    monte_carlo_Wilson95=wilson(k, repetitions)))
            print(json.dumps(dict(completed=scenario, cells_per_type_patient=ncell)), flush=True)
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--repetitions', type=int, default=300)
    parser.add_argument('--permutations', type=int, default=499)
    parser.add_argument('--seed', type=int, default=2026102401)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('preserve existing experiment')
    if args.repetitions < 1 or args.permutations < 19:
        raise ValueError('invalid Monte Carlo size')
    checks = check_numerics()
    start = datetime.now(timezone.utc).isoformat()
    result = dict(kind='ARTIFICIAL_DESIGN_SIMULATION_NOT_REAL_DATA', started_UTC=start,
                  code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  numpy_version=np.__version__, seed=args.seed, permutations=args.permutations,
                  alpha=.05, quantiles=QUANTILES.tolist(), checks=checks,
                  caveats=['Synthetic scores are not sparse RNA or learned programmes.',
                           'Fixed 32 independent patients, two equal sources; not a general power analysis.',
                           'Per-test rejection rates, not BH family power; no method selected by real p-values.',
                           'Source shifts are real marginal associations but null within source.',
                           'Blocked absolute-statistic tests condition on source; no causal conclusion.',
                           'No repeated measures, unequal counts, spatial autocorrelation or missingness modeled.'],
                  results=simulate(args.repetitions, args.permutations, args.seed))
    result['completed_UTC'] = datetime.now(timezone.utc).isoformat()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write('\n')
    print(json.dumps(dict(status='completed', checks=checks, result_rows=len(result['results']))))


if __name__ == '__main__':
    main()
