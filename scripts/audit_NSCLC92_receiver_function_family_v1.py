"""Separate complete rerun of every function statistic and random fingerprint."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC92_receiver_function_math_v2 import prepare, slots, seed, SENSITIVITIES
from NSCLC_frozen_axis_block_math_v2 import evaluate, wild
from NSCLC_frozen_axis_block_audit_math_v1 import evaluate as audit_evaluate, wild as audit_wild
from run_NSCLC92_receiver_function_family_v1 import digest


def correction(values, by=False):
    values = np.asarray(values, float); n = len(values)
    order = np.argsort(values, kind='stable')
    q = values[order] * n / np.arange(1, n + 1)
    if by: q *= np.sum(1. / np.arange(1, n + 1))
    output = np.empty(n); output[order] = np.minimum(1., np.minimum.accumulate(q[::-1])[::-1])
    return output


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 2 * 1024**2)
    out = root / 'NSCLC92_receiver_function_family_audit_20261019_v1'
    assert not out.exists(); out.mkdir(mode=0o700)
    src = root / 'NSCLC92_receiver_function_family_20261019_v1'
    terminal = json.loads((src / 'private_terminal.json').read_text())
    assert terminal['status'] == 'COMPLETE_FULL_FUNCTION_FAMILY_ALL_RANDOM_VALUES_PAIRED_QR_STANDALONE_AUDIT_PENDING'
    for name, d in terminal['files'].items(): assert sha(src / name) == d['sha256']
    inp = root / 'NSCLC92_receiver_function_inputs_20261019_v2'
    proof = root / 'NSCLC92_receiver_function_input_audit_20261019_v2/private_terminal.json'
    assert terminal['input_audit_sha256'] == sha(proof)
    assert json.loads(proof.read_text())['producer_terminal_sha256'] == sha(inp / 'private_terminal.json')
    a = dict(np.load(inp / 'private_complete_function_inputs.npz'))
    r = pd.read_csv(src / 'private_complete_function_family_results.tsv', sep='\t')
    fingerprints = pd.read_csv(src / 'private_all_random_stream_fingerprints.tsv', sep='\t').set_index('slot')
    sensitivity = pd.read_csv(src / 'private_complete_function_fixed_sensitivities.tsv', sep='\t').set_index(['slot', 'sensitivity'])
    definitions = list(slots(a)); assert r.slot.tolist() == [s[0] for s in definitions]
    assert len(sensitivity) == len(definitions) * 2
    cal = json.loads((src / 'private_finite_null_terminal.json').read_text())
    assert cal['status'] == 'ACCEPTED_FINITE_NULL_DIAGNOSTICS'
    assert cal['every_random_value_independently_compared']
    calibration = pd.read_csv(src / 'private_all_finite_null_experiments.tsv', sep='\t')
    for summary in cal['summaries']:
        frame = calibration[calibration.slot.eq(summary['slot'])]
        assert len(frame) == 200
        assert float((frame.null_P <= .05).mean()) == summary['rejection_rate'] <= .10
        assert float(frame.covers_zero.mean()) == summary['coverage'] >= .90
    P = []; maximum = 0.; count = 0; accepted_slots = []
    for i, (slot, k, f, partition, question) in enumerate(definitions):
        m, h, status, valid = prepare(a, k, f, partition, question, independent=True)
        row = r.iloc[i]
        if h is not None:
            try: fit = audit_evaluate(h, m['y'])
            except ValueError as error: status = str(error); h = None
        assert row.status == status and row.author_records == int(valid.sum())
        assert row.profile_blocks == len(set(a['profile_blocks'][valid]))
        if h is None:
            assert row.wild_P == 1 and np.isnan(row.effect)
            assert slot not in fingerprints.index
            P.append(1.)
        else:
            pm, g, _, _ = prepare(a, k, f, partition, question)
            pf = evaluate(g, pm['y'])
            null = audit_wild(h, m['y'], 39999, seed(slot, 'WCRB'), True)
            full = audit_wild(h, m['y'], 19999, seed(slot, 'WCUB'), False)
            prod_null = wild(g, pm['y'], 39999, seed(slot, 'WCRB'), True)
            prod_full = wild(g, pm['y'], 19999, seed(slot, 'WCUB'), False)
            saved = fingerprints.loc[slot]
            assert digest(prod_null) == saved.WCRB_SHA256
            assert digest(prod_full) == saved.WCUB_SHA256
            assert digest(null) == saved.independent_WCRB_SHA256
            assert digest(full) == saved.independent_WCUB_SHA256
            np.testing.assert_allclose(null, prod_null, atol=1e-8, rtol=1e-8)
            np.testing.assert_allclose(full, prod_full, atol=1e-8, rtol=1e-8)
            np.testing.assert_allclose([pf['effect'], pf['se'], pf['t']], [fit['effect'], fit['se'], fit['t']], atol=1e-8, rtol=1e-8)
            extreme = int((np.abs(null) >= abs(fit['t']) - 1e-12).sum())
            p = (1 + extreme) / 40000
            assert extreme == row.extreme_count
            lo, hi = np.quantile(full, [.025, .975])
            np.testing.assert_allclose([row.effect, row.CV3_se, row['t'], row.CI_low, row.CI_high, row.wild_P], [fit['effect'], fit['se'], fit['t'], fit['effect'] - hi * fit['se'], fit['effect'] - lo * fit['se'], p], atol=1e-8, rtol=1e-8)
            saved_leave = np.load(src / ('private_leave_blocks_%04d.npz' % i))['leave_profile_block_effect']
            np.testing.assert_allclose(saved_leave, fit['leave_block_effect'], atol=1e-9, rtol=1e-9)
            maximum = max(maximum, float(np.max(np.abs(null - prod_null))), float(np.max(np.abs(full - prod_full))))
            count += len(null) + len(full); P.append(p); accepted_slots.append(slot)
        for mode in SENSITIVITIES:
            sm, sh, why, _ = prepare(a, k, f, partition, question, mode, independent=True)
            ss = sensitivity.loc[(slot, mode)]
            if sh is not None:
                try: sf = audit_evaluate(sh, sm['y'])
                except ValueError as error: why = str(error); sh = None
            assert ss.status == why
            if sh is None: assert np.isnan(ss.effect)
            else:
                critical = student_t.ppf(.975, sh['groups'] - 1)
                np.testing.assert_allclose([ss.effect, ss.CV3_t_CI_low, ss.CV3_t_CI_high], [sf['effect'], sf['effect'] - critical * sf['se'], sf['effect'] + critical * sf['se']], atol=1e-8, rtol=1e-8)
        (out / 'private_progress.json').write_text(json.dumps(dict(stage='FULL_RANDOM_VECTOR_AND_DELETED_QR_REPLAY', completed_slots=i + 1, total_slots=len(definitions), observed_UTC=datetime.now(timezone.utc).isoformat())) + '\n')
    assert set(accepted_slots) == set(fingerprints.index)
    np.testing.assert_allclose(r.wild_P, P, atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(r.BH_complete_family, correction(P), atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(r.BY_complete_family, correction(P, True), atol=1e-12, rtol=1e-12)
    assert count == terminal['all_random_statistics_independently_compared']
    receipt = dict(status='ACCEPTED_COMPLETE_FUNCTION_FAMILY_ALL_MAIN_AND_SENSITIVITY_EFFECTS_RANDOM_VECTORS_FINGERPRINTS_BH_BY_AND_BLOCK_DELETIONS', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), input_audit_sha256=sha(proof), random_statistics_compared=count, maximum_independent_difference=maximum, resource_guard=guard, scientific_outputs_public=False)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
