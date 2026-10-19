"""Finite full function family, all random vectors compared only in RAM."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_math_v1 import adjust
from NSCLC_frozen_axis_block_math_v2 import evaluate, wild, confidence
from NSCLC_frozen_axis_block_audit_math_v1 import evaluate as other_evaluate, wild as other_wild
from NSCLC92_receiver_function_math_v2 import prepare, slots, seed, SENSITIVITIES
from scipy.stats import t as student_t


def digest(vector):
    return hashlib.sha256(np.asarray(vector, dtype='<f8').tobytes()).hexdigest()


def compare(fit, other, null, n2, full, f2):
    np.testing.assert_allclose([fit['effect'], fit['se'], fit['t']], [other['effect'], other['se'], other['t']], atol=1e-8, rtol=1e-8)
    np.testing.assert_allclose(fit['leave_block_effect'], other['leave_block_effect'], atol=1e-8, rtol=1e-8)
    np.testing.assert_allclose(null, n2, atol=1e-8, rtol=1e-8)
    np.testing.assert_allclose(full, f2, atol=1e-8, rtol=1e-8)
    return max(float(np.max(np.abs(null - n2))), float(np.max(np.abs(full - f2))))


def calibrate(a, coverage, out):
    table = coverage[coverage.status.eq('ESTIMABLE')]
    selected = []
    for _, group in table.groupby(['partition', 'question'], sort=True):
        selected.append(group.sort_values('slot', kind='stable').iloc[0].slot)
    if len(table):
        selected.append(table.sort_values(['max_block_leverage', 'slot'], ascending=[False, True], kind='stable').iloc[0].slot)
    selected = list(dict.fromkeys(selected))
    mapping = {v[0]: v[1:] for v in slots(a)}
    plan = dict(selected_designs=selected, selection='first eligible slot per partition/question plus globally largest block leverage; no true effects', repeats_per_design=200, each_stream_draws=999, sigma='.5+.35abs(zage)+.25abs(ztarget); .6block+.8individual Gaussian', gate='rejectrate<=.10 and coverage>=.90 for every design, otherwise no true effects', created_UTC=datetime.now(timezone.utc).isoformat())
    (out / 'private_plan_before_finite_null_calibration.json').write_text(json.dumps(plan, indent=2) + '\n')
    rows = []; maximum = 0.
    for number, slot in enumerate(selected):
        k, f, partition, question = mapping[slot]
        m, g, _, _ = prepare(a, k, f, partition, question)
        _, h, _, _ = prepare(a, k, f, partition, question, independent=True)
        age = m['controls'][:, 0]; zage = (age - age.mean()) / (age.std() or 1.)
        x = m['x']; zx = (x - x.mean()) / (x.std() or 1.)
        sigma = .5 + .35 * np.abs(zage) + .25 * np.abs(zx)
        mu = m['D'][:, :-1] @ np.linspace(.1, .3, m['nuisance_rank'])
        rng = np.random.default_rng(seed(slot, 'calibration_outcomes'))
        for repeat in range(200):
            y = mu + sigma * (.6 * rng.normal(size=g['groups'])[g['membership']] + .8 * rng.normal(size=len(mu)))
            fit = evaluate(g, y); other = other_evaluate(h, y)
            ns = seed(slot, 'calibration_null_' + str(repeat)); fs = seed(slot, 'calibration_full_' + str(repeat))
            null = wild(g, y, 999, ns, True); full = wild(g, y, 999, fs, False)
            n2 = other_wild(h, y, 999, ns, True); f2 = other_wild(h, y, 999, fs, False)
            maximum = max(maximum, compare(fit, other, null, n2, full, f2))
            p = (1 + int((np.abs(null) >= abs(fit['t']) - 1e-12).sum())) / 1000
            ci = confidence(fit, full)
            rows.append(dict(slot=slot, repeat=repeat, null_P=p, covers_zero=bool(ci[0] <= 0 <= ci[1])))
        (out / 'private_progress.json').write_text(json.dumps(dict(stage='FINITE_NULL_CALIBRATION_NO_TRUE_EFFECTS', completed_designs=number + 1, total_designs=len(selected), observed_UTC=datetime.now(timezone.utc).isoformat())) + '\n')
    frame = pd.DataFrame(rows, columns=['slot', 'repeat', 'null_P', 'covers_zero'])
    frame.to_csv(out / 'private_all_finite_null_experiments.tsv', sep='\t', index=False)
    summary = [dict(slot=slot, rejection_rate=float((g.null_P <= .05).mean()), coverage=float(g.covers_zero.mean())) for slot, g in frame.groupby('slot', sort=False)]
    passed = bool(selected) and all(s['rejection_rate'] <= .10 and s['coverage'] >= .90 for s in summary)
    proof = dict(status='ACCEPTED_FINITE_NULL_DIAGNOSTICS' if passed else 'FAILED_FINITE_NULL_DIAGNOSTICS_TRUE_EFFECTS_MUST_NOT_RUN', completed_UTC=datetime.now(timezone.utc).isoformat(), summaries=summary, every_random_value_independently_compared=True, maximum_independent_difference=maximum, synthetic_experiments=len(rows))
    (out / 'private_finite_null_terminal.json').write_text(json.dumps(proof, indent=2) + '\n')
    return passed


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 4 * 1024**2)
    out = root / 'NSCLC92_receiver_function_family_20261019_v1'
    assert not out.exists(); out.mkdir(mode=0o700)
    src = root / 'NSCLC92_receiver_function_inputs_20261019_v2'
    proof = root / 'NSCLC92_receiver_function_input_audit_20261019_v2/private_terminal.json'
    at = json.loads(proof.read_text())
    assert at['status'] == 'ACCEPTED_AMENDED_FUNCTION_INPUTS_COMPLETE_COMPONENT_SCALES_AND_THREE_QUESTION_DESIGNS'
    assert at['producer_terminal_sha256'] == sha(src / 'private_terminal.json')
    t = json.loads((src / 'private_terminal.json').read_text())
    for name, d in t['files'].items(): assert sha(src / name) == d['sha256']
    save = lambda n, d: (out / n).write_text(json.dumps(d, indent=2) + '\n')
    save('plan_before_complete_function_family.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), math_source_sha256=sha(Path(__file__).with_name('NSCLC92_receiver_function_math_v2.py')), input_audit_sha256=sha(proof), resource_guard=guard, main_draws=39999, interval_draws=19999, full_joint_family='All50 x2original partitions x3questions xall originalBH receiver scopes. All failures p1; one BH/BY denominator.', all_scientific_outputs_private=True))
    a = dict(np.load(src / 'private_complete_function_inputs.npz'))
    coverage = pd.read_csv(src / 'private_complete_family_design_coverage.tsv', sep='\t')
    if not calibrate(a, coverage, out):
        save('private_terminal.json', dict(status='STOPPED_FINITE_NULL_GATE_FAILED_NO_TRUE_FUNCTION_EFFECTS', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), scientific_outputs_public=False))
        print(json.dumps(dict(status='STOPPED_FINITE_NULL_GATE_FAILED_NO_TRUE_FUNCTION_EFFECTS', terminal_sha256=sha(out / 'private_terminal.json')))); return
    results = []; fingerprints = []; sensitivity = []; maximum = 0.; random_values = 0
    lookup = coverage.set_index('slot')
    for number, (slot, k, f, partition, question) in enumerate(slots(a)):
        m, g, status, valid = prepare(a, k, f, partition, question)
        row = lookup.loc[slot].to_dict()
        row.update(slot=slot, status=status, effect=np.nan, CV3_se=np.nan, t=np.nan, CI_low=np.nan, CI_high=np.nan, wild_P=1., extreme_count=0)
        if g is not None:
            try: fit = evaluate(g, m['y'])
            except ValueError as error: row['status'] = str(error); g = None
        if g is not None:
            _, h, hs, _ = prepare(a, k, f, partition, question, independent=True)
            assert h is not None and hs == 'ESTIMABLE'
            null = wild(g, m['y'], 39999, seed(slot, 'WCRB'), True)
            full = wild(g, m['y'], 19999, seed(slot, 'WCUB'), False)
            other = other_evaluate(h, m['y'])
            n2 = other_wild(h, m['y'], 39999, seed(slot, 'WCRB'), True)
            f2 = other_wild(h, m['y'], 19999, seed(slot, 'WCUB'), False)
            difference = compare(fit, other, null, n2, full, f2)
            maximum = max(maximum, difference); random_values += len(null) + len(full)
            extreme = int((np.abs(null) >= abs(fit['t']) - 1e-12).sum())
            assert extreme == int((np.abs(n2) >= abs(other['t']) - 1e-12).sum())
            interval = confidence(fit, full)
            row.update(effect=fit['effect'], CV3_se=fit['se'], t=fit['t'], CI_low=interval[0], CI_high=interval[1], wild_P=(1 + extreme) / 40000, extreme_count=extreme)
            fingerprints.append(dict(slot=slot, WCRB_SHA256=digest(null), WCUB_SHA256=digest(full), independent_WCRB_SHA256=digest(n2), independent_WCUB_SHA256=digest(f2), maximum_independent_difference=difference))
            np.savez_compressed(out / ('private_leave_blocks_%04d.npz' % number), leave_profile_block_effect=fit['leave_block_effect'])
        results.append(row)
        for mode in SENSITIVITIES:
            sm, sg, status, _ = prepare(a, k, f, partition, question, mode)
            s = dict(slot=slot, sensitivity=mode, status=status, effect=np.nan, CV3_t_CI_low=np.nan, CV3_t_CI_high=np.nan)
            if sg is not None:
                try: sf = evaluate(sg, sm['y'])
                except ValueError as error: s['status'] = str(error); sg = None
            if sg is not None:
                critical = student_t.ppf(.975, sg['groups'] - 1)
                s.update(effect=sf['effect'], CV3_t_CI_low=sf['effect'] - critical * sf['se'], CV3_t_CI_high=sf['effect'] + critical * sf['se'])
            sensitivity.append(s)
        save('private_progress.json', dict(stage='COMPLETE_FROZEN_FUNCTION_FAMILY', completed_slots=number + 1, total_slots=len(coverage), observed_UTC=datetime.now(timezone.utc).isoformat()))
    r = pd.DataFrame(results)
    assert len(r) == len(coverage)
    r['BH_complete_family'] = adjust(r.wild_P)
    r['BY_complete_family'] = adjust(r.wild_P, by=True)
    r.to_csv(out / 'private_complete_function_family_results.tsv', sep='\t', index=False)
    pd.DataFrame(sensitivity).to_csv(out / 'private_complete_function_fixed_sensitivities.tsv', sep='\t', index=False)
    pd.DataFrame(fingerprints).to_csv(out / 'private_all_random_stream_fingerprints.tsv', sep='\t', index=False)
    files = {p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}
    assert sum(d['bytes'] for d in files.values()) < 4 * 1024**2
    receipt = dict(status='COMPLETE_FULL_FUNCTION_FAMILY_ALL_RANDOM_VALUES_PAIRED_QR_STANDALONE_AUDIT_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), input_audit_sha256=sha(proof), files=files, all_random_statistics_independently_compared=random_values, maximum_independent_statistic_difference=maximum, random_vectors_persisted=False, scientific_outputs_public=False)
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
