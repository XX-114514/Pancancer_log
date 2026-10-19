"""Finite two-question rank family; human/random arrays never persisted."""
from collections import Counter
import numpy as np
from scipy.stats import beta as beta_distribution
from statsmodels.stats.multitest import multipletests
from immune_RNA_MCP_new_source_math_v1 import (
    design, independent_design, estimate, evaluate, wild_values, independent_wild,
    bootstrap_indices, bootstrap, leave_one, seed, adjust)
from HNSCC_author_state_clinical_kernel_v2 import fit


def close(a, b, tol=1e-8):
    assert np.allclose(a, b, rtol=0, atol=tol, equal_nan=True)


def infer(questions, C, names, strata, extra, concordant, save):
    models = []
    for slot, x, y in questions:
        m, status = design(x, y, C, names, 'coordination')
        im, istatus = independent_design(x, y, C, names, 'coordination')
        assert status == istatus
        models.append(dict(slot=slot, x=x, y=y, m=m, im=im, status=status))
    nulls = []
    for entry in models:
        if entry['m'] is None:
            continue
        m, im = entry['m'], entry['im']; D = m['D']; n = len(D)
        mu = D[:, :-1] @ np.linspace(-.2, .2, D.shape[1] - 1)
        quality = C[:, 0]
        quality = (quality - quality.mean()) / quality.std(ddof=0) if quality.std(ddof=0) > 0 else np.zeros(n)
        passed = True
        for scenario in ['normal', 't3_unit_variance', 'RNA_quality_heteroskedastic_normal']:
            rng = np.random.default_rng(seed(entry['slot'], 'null_errors_' + scenario))
            rejected = 0; maximum = 0.
            for simulation in range(256):
                err = rng.standard_t(3, n) / np.sqrt(3) if scenario == 't3_unit_variance' else rng.normal(size=n)
                if scenario == 'RNA_quality_heteroskedastic_normal':
                    err *= np.exp(.6 * quality)
                yy = mu + err; g = fit(yy, D)
                random_seed = seed(entry['slot'], scenario + '|' + str(simulation))
                v = wild_values(g, 999, random_seed)
                ig, iv = independent_wild(im, 999, random_seed, y=yy)
                close([g['beta'], g['t']], [ig['beta'], ig['t']]); close(v, iv)
                count = int((np.abs(v) >= abs(g['t']) - 1e-12).sum())
                assert count == int((np.abs(iv) >= abs(ig['t']) - 1e-12).sum())
                maximum = max(maximum, float(np.max(np.abs(v - iv))))
                rejected += (count + 1) / 1000 < .05
            rate = rejected / 256; passed &= rate <= .10
            nulls.append(dict(slot=entry['slot'], scenario=scenario, simulations=256, rejections=rejected,
                              rejection_rate=rate, maximum_independent_error=maximum))
            save('private_progress.json', dict(stage='FIXED_WORKING_NULLS', completed_scenarios=len(nulls)))
        entry['null_passed'] = bool(passed)
    save('private_all_working_null_diagnostics.json', nulls)
    result = []
    for e in models:
        m, im = e['m'], e['im']
        r = dict(slot=e['slot'], n_patients=len(e['x']), status=e['status'], effect=None, P=None,
                 MC_upper_P=None, CI_low=None, CI_high=None)
        if m is not None and not e['null_passed']:
            r['status'] = 'FIXED_NULL_DIAGNOSTIC_FAILED_NO_ADJUSTED_INFERENCE'
        elif m is not None:
            g, effect = estimate(m)
            ig, iv = independent_wild(im, 19999, seed(e['slot'], 'main_wild'))
            v = wild_values(g, 19999, seed(e['slot'], 'main_wild'))
            close([g['beta'], g['t'], effect], [ig['beta'], ig['t'], ig['effect']]); close(v, iv)
            count = int((np.abs(v) >= abs(g['t']) - 1e-12).sum())
            assert count == int((np.abs(iv) >= abs(ig['t']) - 1e-12).sum())
            takes = bootstrap_indices(strata, e['slot'], 2000)
            b, br = bootstrap(e['x'], e['y'], C, names, 'coordination', takes, m['nuisance_rank'])
            ib, ibr = bootstrap(e['x'], e['y'], C, names, 'coordination', takes, im['nuisance_rank'], independent=True)
            loo, lr = leave_one(e['x'], e['y'], C, names, 'coordination', m['nuisance_rank'])
            iloo, ilr = leave_one(e['x'], e['y'], C, names, 'coordination', im['nuisance_rank'], independent=True)
            close(b, ib); close(loo, iloo)
            assert np.array_equal(br, ibr) and np.array_equal(lr, ilr)
            valid = np.isfinite(b); lv = np.isfinite(loo)
            ci = np.quantile(b[valid], [.025, .975]) if valid.sum() >= 1900 else [None, None]
            upper = 1. if count == 19999 else float(beta_distribution.ppf(.995, count + 1, 19999 - count))
            r.update(status='ESTIMABLE_SOURCE_EXPLORATORY_ASSOCIATION', effect=effect, P=(count + 1) / 20000,
                     MC_upper_P=min(1., (1 + 19999 * upper) / 20000), HC3_t=g['t'],
                     CI_low=ci[0], CI_high=ci[1], valid_bootstrap=int(valid.sum()),
                     bootstrap_failure_counts=dict(Counter(br[~valid].tolist())), nuisance_rank=m['nuisance_rank'],
                     residual_df=m['df'], condition=m['condition'], max_leverage=m['max_leverage'],
                     LOO_valid=int(lv.sum()), LOO_min=float(loo[lv].min()) if lv.any() else None,
                     LOO_max=float(loo[lv].max()) if lv.any() else None,
                     LOO_all_same_sign=bool(lv.all() and (loo * effect > 0).all()),
                     maximum_independent_main_statistic_error=float(np.max(np.abs(v - iv))))
        masks = [('unadjusted', np.zeros(len(names), bool)),
                 ('clinical_only', np.array([s.startswith('agent::') for s in names])),
                 ('RNA_quality_only', np.array(['logRNA' in s or 'loggenes' in s for s in names])),
                 ('capture_only', np.array(['fraction' in s or 'logcells' in s for s in names]))]
        points = []
        for mode, take in masks + [('plus_frozen_subtypes', np.ones(len(names), bool)),
                                  ('author_binary_concordant', np.ones(len(names), bool))]:
            cc = C[:, take]; nn = names[take]; choose = np.ones(len(C), bool)
            if mode == 'plus_frozen_subtypes':
                cc = np.c_[cc, extra]
                nn = np.r_[nn, ['candidate_CD8_fraction', 'candidate_Treg_fraction', 'candidate_Macro_fraction']]
            if mode == 'author_binary_concordant':
                choose = concordant
            pm, ps = design(e['x'][choose], e['y'][choose], cc[choose], nn, 'coordination')
            ipm, ips = independent_design(e['x'][choose], e['y'][choose], cc[choose], nn, 'coordination')
            assert ps == ips
            point = None
            if pm is not None:
                _, point = estimate(pm); close(point, evaluate(ipm)['effect'])
            points.append(dict(mode=mode, n_patients=int(choose.sum()), status=ps, effect=point, no_new_P=True))
        r['predeclared_point_sensitivities'] = points
        result.append(r)
    p = np.array([r['P'] if r['P'] is not None else 1. for r in result])
    upper = np.array([r['MC_upper_P'] if r['MC_upper_P'] is not None else 1. for r in result])
    bh = adjust(p); by = adjust(p, by=True); mc = adjust(upper)
    close(bh, multipletests(p, method='fdr_bh')[1], 1e-12)
    close(by, multipletests(p, method='fdr_by')[1], 1e-12)
    close(mc, multipletests(upper, method='fdr_bh')[1], 1e-12)
    for k, r in enumerate(result):
        r.update(BH2=float(bh[k]) if r['P'] is not None else None,
                 BY2=float(by[k]) if r['P'] is not None else None,
                 MC_upper_BH2=float(mc[k]) if r['P'] is not None else None)
    return result
