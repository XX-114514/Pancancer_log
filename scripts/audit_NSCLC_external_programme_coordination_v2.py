"""Independent OLS reconstruction of all real conditional resampling streams."""
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist, rankdata, spearmanr

from NSCLC_external_programme_coordination_math_v1 import seed
import NSCLC_external_exact_rank_guards_v1 as exact_guard
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from test_NSCLC_external_programme_coordination_v1 import read_inputs, QCOLS


def calculate(a, b, quality, blocks, minimum=20):
    blocks = np.asarray(blocks, str)
    counts = {s: int(np.sum(blocks == s)) for s in set(blocks)}
    keep = np.array([counts[s] >= 2 for s in blocks], bool)
    if keep.sum() < minimum:
        raise ValueError('FEWER_THAN20_UNITS_IN_NON_SINGLETON_BLOCKS' if minimum == 20 else 'INSUFFICIENT_INFORMATIVE_UNITS')
    a, b, quality, blocks = np.asarray(a)[keep], np.asarray(b)[keep], np.asarray(quality)[keep], blocks[keep]
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise ValueError('NONFINITE_PROGRAMME')
    if not np.isfinite(quality).all():
        raise ValueError('NONFINITE_QUALITY')
    qr = rankdata(quality, axis=0) if quality.shape[1] else quality
    sd = qr.std(axis=0)
    standard = (qr-qr.mean(axis=0))/np.where(sd>0, sd, 1)
    z = np.column_stack([np.ones(len(a)), *[(blocks == s).astype(float) for s in sorted(set(blocks))], standard])
    y = rankdata(np.column_stack([a, b]), axis=0)
    coefficient, _, rank, _ = np.linalg.lstsq(z, y, rcond=1e-10)
    if len(a)-rank < 8:
        raise ValueError('QUALITY_AND_BLOCK_DESIGN_INSUFFICIENT_RESIDUAL_DF')
    residual = y-z@coefficient
    # Only diagnose near-zero floating residuals; exact span decides zero.
    if np.min(np.linalg.norm(residual, axis=0)) < 1e-6:
        integer_design = np.column_stack([np.ones(len(a)), *[(blocks == s).astype(int) for s in sorted(set(blocks))], 2*qr]).astype(np.int64)
        integer_y = (2*y).astype(np.int64)
        joint = np.unique(np.column_stack([integer_design, integer_y]), axis=0)
        exact_z = joint[:, :integer_design.shape[1]]
        full_row_certificate = len(joint) == rank and exact_guard.modular_rank(exact_z) == rank
        # A full exact span must agree with the unchanged numeric rank.
        # Otherwise it is not evidence about the specified truncated span.
        if full_row_certificate or (exact_guard.rational_rank(exact_z) == rank and exact_guard.any_exactly_fitted(integer_design, integer_y)):
            raise ValueError('INDIVIDUAL_PROGRAMME_RESIDUAL_CONSTANT')
    if np.any(np.linalg.norm(residual, axis=0) < 1e-10):
        raise ValueError('INDIVIDUAL_PROGRAMME_RESIDUAL_CONSTANT')
    rho = float(np.clip(np.corrcoef(residual.T)[0, 1], -1, 1))
    return dict(a=a, b=b, quality=quality, blocks=blocks, keep=keep, z=z, y=y, residual=residual, rank=rank, rho=rho)


def fdr(values, by=False):
    order = np.argsort(values)
    n = len(values)
    factor = np.sum(1/np.arange(1, n+1)) if by else 1.
    adjusted = np.minimum.accumulate((values[order]*n*factor/np.arange(1, n+1))[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(adjusted, 1.)
    return out


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    src = root/'NSCLC_external_programme_coordination_20261016_v1'
    owned = root/'HNSCC_receiver_followup_job_control_20261016_v2/test_NSCLC_external_programme_coordination_v1/launch_receipt.json'
    receipt = json.loads(owned.read_text())
    deadline = time.monotonic()+7200
    while not (src/'private_terminal.json').exists():
        proc = Path('/proc')/str(receipt['pid'])/'stat'
        assert proc.exists(), 'Producer missing without terminal; inspect original log'
        values = proc.read_text().rsplit(')', 1)[1].split()
        assert values[0] != 'Z' and values[19] == str(receipt['proc_start_ticks'])
        assert time.monotonic() < deadline, 'Observation ended; do not infer producer stopped'
        time.sleep(3)
    resource(root, 1048576)
    terminal = json.loads((src/'private_terminal.json').read_text())
    assert terminal['source_sha256'] == receipt['source_sha256']
    for name, item in terminal['files'].items():
        assert sha(src/name) == item['sha256']
    _, _, contract, slots, cases = read_inputs(root)
    plan = json.loads((src/'plan_before_all_external_coordination_effects.json').read_text())
    assert plan['contract'] == contract
    scope = pd.read_csv(src/'private_all15600_primary_test_eligibility.tsv', sep='\t', keep_default_na=False)
    results = pd.read_csv(src/'private_all_estimated_programme_coordination.tsv', sep='\t')
    assert len(scope) == 15600 and scope.record.is_unique
    assert results.record.is_unique and results.test_index.tolist() == list(range(len(results)))
    saved = np.load(src/'private_all_bootstrap_and_LOO.npz', allow_pickle=False)
    assert saved['records'].astype(str).tolist() == results.record.tolist()
    assert saved['bootstrap'].shape == (len(results), 2000)
    scope_index = scope.set_index('record')
    prepared = []
    for row in slots.itertuples():
        g0 = cases.iloc[row.row_start:row.row_stop]
        for family in ['all_units', 'reviewed_patients']:
            record = row.slot+'|'+family
            expected_status = 'ELIGIBLE_FOR_FIXED_PRIMARY_TEST'
            g = g0 if family == 'all_units' else g0[g0.reviewed_patient]
            observed = scope_index.loc[record]
            assert observed.candidate_units == len(g) and observed.candidate_reviewed_units == g.reviewed_patient.sum()
            if row.in_model_source:
                expected_status = 'ORIGINAL_MODEL_SOURCE_NOT_EXTERNAL_FAMILY'
            else:
                try:
                    p = calculate(g.raw_A.to_numpy(), g.raw_B.to_numpy(), g[QCOLS].to_numpy(), g.block_key.to_numpy(str))
                    assert observed.test_index == len(prepared)
                    assert observed.active_units == len(p['a']) and observed.quality_rank == p['rank']
                    assert observed.residual_df == len(p['a'])-p['rank']
                    prepared.append((record, p, g.iloc[np.flatnonzero(p['keep'])]))
                except ValueError as e:
                    expected_status = str(e)
            assert observed.status == expected_status, (record, expected_status, observed.status)
    assert len(prepared) == len(results)
    sizes = {f: int((results.family == f).sum()) for f in scope.family.unique()}
    assert sizes == plan['family_sizes_before_any_effect']
    maximum_boot = maximum_effect = 0.
    perms_checked = boots_checked = loo_checked = 0
    checkpoint = root/'HNSCC_receiver_followup_job_control_20261016_v2'/Path(__file__).stem/'private_audit_checkpoint.json'
    for i, (record, p, frame) in enumerate(prepared):
        row = results.iloc[i]
        assert row.record == record
        assert row.active_units == len(frame) and row.reviewed_units == frame.reviewed_patient.sum()
        maximum_effect = max(maximum_effect, abs(row.partial_rank_rho-p['rho']))
        assert np.isclose(row.partial_rank_rho, p['rho'], atol=1e-10, rtol=1e-10)
        z, residual, blocks = p['z'], p['residual'], p['blocks']
        pinv = np.linalg.pinv(z, rcond=1e-10)
        fit = p['y'][:, 1]-residual[:, 1]
        groups = [np.flatnonzero(blocks == s) for s in sorted(set(blocks))]
        rng = np.random.default_rng(seed(record, 'Freedman_Lane'))
        count = 0
        for start in range(0, 19999, 256):
            m = min(256, 19999-start)
            take = np.broadcast_to(np.arange(len(frame)), (m, len(frame))).copy()
            for ids in groups:
                take[:, ids] = ids[np.argsort(rng.random((m, len(ids))), axis=1)]
            null = fit[None, :] + residual[take, 1]
            rnull = null-(null@pinv.T)@z.T
            rho = rnull@residual[:, 0]/(np.linalg.norm(rnull, axis=1)*np.linalg.norm(residual[:, 0]))
            count += int((abs(rho) >= abs(p['rho'])-1e-12).sum())
        assert count == row.permutation_exceedances
        assert np.isclose(row.P, (1+count)/20000, atol=1e-15)
        upper = 1. if count == 19999 else float(beta_dist.ppf(1-.01/sizes[row.family], count+1, 19999-count))
        assert np.isclose(row.MC_upper_P, min(1., (1+19999*upper)/20000), atol=1e-10)
        perms_checked += 19999
        rng = np.random.default_rng(seed(record, 'bootstrap'))
        bootstrap = np.full(2000, np.nan)
        for b in range(2000):
            take = np.concatenate([ids[rng.integers(len(ids), size=len(ids))] for ids in groups])
            try:
                bootstrap[b] = calculate(p['a'][take], p['b'][take], p['quality'][take], blocks[take])['rho']
            except ValueError:
                pass
        assert np.allclose(bootstrap, saved['bootstrap'][i], equal_nan=True, atol=1e-9, rtol=1e-9)
        valid = np.isfinite(bootstrap)
        maximum_boot = max(maximum_boot, float(np.nanmax(abs(bootstrap-saved['bootstrap'][i]))))
        assert row.valid_bootstraps == valid.sum()
        ci = np.quantile(bootstrap[valid], [.025, .975]) if valid.sum() >= 1900 else [np.nan, np.nan]
        assert np.allclose([row.CI_low, row.CI_high], ci, equal_nan=True, atol=1e-9, rtol=1e-9)
        boots_checked += 2000
        loo = np.full(len(frame), np.nan)
        for j in range(len(frame)):
            take = np.arange(len(frame)) != j
            try:
                loo[j] = calculate(p['a'][take], p['b'][take], p['quality'][take], blocks[take], minimum=3)['rho']
            except ValueError:
                pass
        bounds = saved['LOO_offsets'][i:i+2]
        assert np.allclose(loo, saved['LOO'][bounds[0]:bounds[1]], equal_nan=True, atol=1e-9, rtol=1e-9)
        assert row.LOO_valid == np.isfinite(loo).sum() and row.LOO_all_valid == np.isfinite(loo).all()
        assert row.LOO_same_sign == bool(np.isfinite(loo).all() and (np.sign(loo) == np.sign(p['rho'])).all())
        assert np.allclose([row.LOO_min, row.LOO_max], [np.nanmin(loo), np.nanmax(loo)] if np.isfinite(loo).any() else [np.nan, np.nan], equal_nan=True, atol=1e-9)
        loo_checked += len(frame)
        assert np.isclose(row.raw_Spearman, spearmanr(p['a'], p['b']).statistic, atol=1e-12)
        empty = np.empty((len(frame), 0))
        for target, aa, bb, quality, status_field in [
            ('block_only_rho', p['a'], p['b'], empty, 'block_only_status'),
            ('source_QC_block_rho', frame.source_QC_A.to_numpy(), frame.source_QC_B.to_numpy(), empty, 'source_QC_status'),
            ('abundance_adjusted_rho', p['a'], p['b'], np.column_stack([p['quality'], frame[['fraction_A', 'fraction_B']].to_numpy()]), 'abundance_status')]:
            try:
                expected = calculate(aa, bb, quality, blocks)['rho']
                expected_status = 'ESTIMABLE'
            except ValueError as e:
                expected = np.nan
                expected_status = str(e)
            assert row[status_field] == expected_status
            assert np.isclose(row[target], expected, equal_nan=True, atol=1e-9, rtol=1e-9)
        if (i+1) % 10 == 0 or i+1 == len(results):
            checkpoint.write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), completed_tests=i+1, total_tests=len(results)))+'\n')
    for family, g in results.groupby('family'):
        for target, source, by in [('BH_q', 'P', False), ('BY_q', 'P', True), ('MC_upper_BH_q', 'MC_upper_P', False)]:
            assert np.allclose(g[target], fdr(g[source].to_numpy(), by), atol=1e-12)
    proof = src/'private_independent_coordination_audit.json'
    assert not proof.exists()
    report = dict(status='PASS_ALL15600_ELIGIBILITY_ALL_REAL_EFFECTS_RESAMPLING_SENSITIVITIES_AND_FAMILY_CORRECTIONS',
                  finished_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
                  producer_terminal_sha256=sha(src/'private_terminal.json'),
                  exact_rank_guard_sha256=sha(exact_guard.__file__),
                  implementation_repair='For near-zero OLS residuals only, exact integer-rank span certificates matching the unchanged numerical nuisance rank correct spurious nonzero residuals caused by cancellation. Original production, thresholds, test family, saved values and random streams unchanged.',
                  all_eligibility_slots_checked=len(scope), estimated_tests_checked=len(results),
                  permutations_replayed=perms_checked, bootstraps_replayed=boots_checked, LOO_replayed=loo_checked,
                  maximum_effect_error=maximum_effect, maximum_bootstrap_error=maximum_boot,
                  scientific_clinical_or_direct_LR_validation=False)
    proof.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status=report['status'], finished_UTC=report['finished_UTC'], audit_sha256=sha(proof))))


if __name__ == '__main__':
    main()
