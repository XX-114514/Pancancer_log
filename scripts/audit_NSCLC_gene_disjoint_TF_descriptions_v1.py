"""Independent least-squares replay of every fixed TF point and case LOO."""
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def independent_points(x, y, nuisance):
    a = rankdata(x)
    b = np.column_stack([rankdata(y[:, j]) for j in range(y.shape[1])])
    if np.std(a) <= 1e-12:
        return np.full(y.shape[1], np.nan)
    full = np.column_stack([nuisance, (a - np.mean(a)) / np.std(a)])
    if np.linalg.matrix_rank(nuisance) != nuisance.shape[1] or np.linalg.matrix_rank(full) != full.shape[1] or np.linalg.cond(full) >= 1e8:
        return np.full(y.shape[1], np.nan)
    combined = np.column_stack([a, b])
    coef = np.linalg.lstsq(nuisance, combined, rcond=None)[0]
    e = combined - nuisance @ coef
    result = []
    for j in range(y.shape[1]):
        denominator = float(np.sqrt(np.dot(e[:, 0], e[:, 0]) * np.dot(e[:, j + 1], e[:, j + 1])))
        result.append(float(np.dot(e[:, 0], e[:, j + 1]) / denominator) if denominator > 1e-12 else np.nan)
    return np.asarray(result)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'NSCLC_gene_disjoint_TF_response_descriptions_20261016_v1'
    job = root / 'HNSCC_receiver_followup_job_control_20261016_v2/describe_NSCLC_gene_disjoint_TF_response_v1/launch_receipt.json'
    launch = json.loads(job.read_text())
    deadline = time.monotonic() + 14400
    while not (out / 'private_terminal.json').exists():
        proc = Path('/proc') / str(launch['pid']) / 'stat'
        assert proc.exists()
        state = proc.read_text().rsplit(')', 1)[1].split()
        assert state[19] == str(launch['proc_start_ticks']) and state[0] != 'Z' and time.monotonic() < deadline
        time.sleep(15)
    dest = out / 'private_independent_audit.json'
    assert not dest.exists()
    terminal = json.loads((out / 'private_terminal.json').read_text())
    for name, digest in terminal['files'].items():
        assert sha(out / name) == digest
    plan = json.loads((out / 'plan_before_any_TF_association_description.json').read_text())
    assert plan['code_sha256'] == launch['source_sha256']
    inputs = root / 'NSCLC_gene_disjoint_receiver_TF_measurement_20261016_v1'
    assert sha(inputs / 'private_independent_audit.json') == plan['input_audit_sha256']
    scope_path = inputs / 'private_all_TF_measurement_scope.tsv'
    hash_path = inputs / 'private_ready_design_SHA256.tsv'
    assert sha(scope_path) == plan['input_scope_sha256'] and sha(hash_path) == plan['input_hash_table_sha256']
    scope = pd.read_csv(scope_path, sep='\t')
    data = pd.read_csv(out / 'private_all2568_TF_functional_descriptions.tsv', sep='\t')
    for col in scope.columns:
        pd.testing.assert_series_equal(scope[col], data[col], check_dtype=False, check_names=False)
    hashes = pd.read_csv(hash_path, sep='\t').set_index('TF_record')
    loo = pd.read_csv(out / 'private_all_matched_case_LOO.tsv.gz', sep='\t')
    assert not loo.duplicated(['TF_record', 'removed_case_index']).any()
    ready = scope[scope.measurement_status.eq('READY_FOR_INDEPENDENT_MEASUREMENT_AUDIT')]
    assert set(ready.TF_record) == set(hashes.index) == set(loo.TF_record)
    numerical = ['raw_spearman', 'partial_rank_rho', 'sender_programme_adjusted_rho', 'pathology_adjusted_rho', 'TF_to_original_masked_receiver_MCP_rho']
    max_error, checked, loo_checked = 0., 0, 0
    for parent_id, rows in ready.groupby('record', sort=True):
        designs = []
        for row in rows.itertuples():
            p = inputs / 'private_frozen_designs' / (row.TF_record + '.npz')
            assert sha(p) == hashes.loc[row.TF_record, 'sha256']
            designs.append(np.load(p, allow_pickle=False))
        z = designs[0]
        Y = np.column_stack([d['Y'] for d in designs])
        expected = {'raw_spearman': independent_points(z['X'], Y, np.ones((len(Y), 1))), 'partial_rank_rho': independent_points(z['X'], Y, z['Z']), 'TF_to_original_masked_receiver_MCP_rho': independent_points(z['original_masked_receiver_programme'], Y, z['Z'])}
        for name, key in [('sender_programme', 'sender_programme'), ('pathology', 'response')]:
            extra = z[key]
            if np.std(extra) > 1e-12:
                zz = np.column_stack([z['Z'], (extra - np.mean(extra)) / np.std(extra)])
                expected[name + '_adjusted_rho'] = independent_points(z['X'], Y, zz)
            else:
                expected[name + '_adjusted_rho'] = np.full(len(rows), np.nan)
        all_loo = np.full((len(Y), len(rows)), np.nan)
        for k in range(len(Y)):
            take = np.arange(len(Y)) != k
            all_loo[k] = independent_points(z['X'][take], Y[take], z['Z'][take])
        for j, row in enumerate(rows.itertuples()):
            f = data[data.TF_record.eq(row.TF_record)].iloc[0]
            for col in numerical:
                value = float(expected[col][j])
                observed = float(f[col])
                assert np.isclose(observed, value, rtol=1e-10, atol=1e-10, equal_nan=True)
                if np.isfinite(value):
                    max_error = max(max_error, abs(observed - value))
            status = 'DESCRIPTIVE_FULL_CONTROL_POINT_ONLY' if np.isfinite(expected['partial_rank_rho'][j]) else 'UNDEFINED_FULL_CONTROL_RESIDUAL_SCORE'
            assert f.description_status == status and bool(f.no_new_P_BH) and not bool(f.independent_samples)
            ll = loo[loo.TF_record.eq(row.TF_record)].sort_values('removed_case_index')
            assert ll.removed_case_index.tolist() == list(range(len(Y)))
            assert np.allclose(ll.rho, all_loo[:, j], rtol=1e-10, atol=1e-10, equal_nan=True)
            valid = np.isfinite(all_loo[:, j])
            same = np.sign(all_loo[valid, j]) == np.sign(expected['partial_rank_rho'][j])
            assert f.LOO_estimable == valid.sum() and f.LOO_same_sign == same.sum()
            for col, value in [('LOO_rho_minimum', np.min(all_loo[valid, j]) if valid.any() else np.nan), ('LOO_rho_maximum', np.max(all_loo[valid, j]) if valid.any() else np.nan)]:
                assert np.isclose(f[col], value, atol=1e-10, rtol=1e-10, equal_nan=True)
            checked += 1
            loo_checked += len(ll)
        for d in designs:
            d.close()
    missing = data[~data.TF_record.isin(ready.TF_record)]
    assert missing[numerical].isna().all().all()
    assert missing.description_status.tolist() == missing.measurement_status.tolist()
    assert terminal['status_counts'] == data.description_status.value_counts().to_dict()
    assert terminal['all_fixed_slots'] == len(data) and terminal['described_slots'] == checked
    assert terminal['all_case_LOO_rows'] == loo_checked == len(loo)
    result = dict(status='PASS_ALL_FIXED_TF_POINTS_ALL_REQUIRED_CONTROL_SENSITIVITIES_ALL_CASE_LOO_AND_UNAVAILABLE_SLOTS_BY_INDEPENDENT_LEAST_SQUARES', finished_UTC=datetime.now(timezone.utc).isoformat(), audit_code_sha256=sha(__file__), producer_terminal_sha256=sha(out / 'private_terminal.json'), all_fixed_slots=len(data), all_described_TFs=checked, all_case_LOO_records=loo_checked, maximum_point_error=max_error, no_new_P_BH=True, independent_samples_or_protein_validation=False)
    dest.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
