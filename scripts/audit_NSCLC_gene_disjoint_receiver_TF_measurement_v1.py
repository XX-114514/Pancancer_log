"""All fixed TF scores versus pinned official ULM, and every saved case design."""
import ast
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


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'NSCLC_gene_disjoint_receiver_TF_measurement_20261016_v1'
    launch = json.loads((root / 'HNSCC_receiver_followup_job_control_20261016_v2/measure_NSCLC_gene_disjoint_receiver_TF_v1/launch_receipt.json').read_text())
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
    plan = json.loads((out / 'plan_before_any_TF_measurement.json').read_text())
    assert plan['code_sha256'] == launch['source_sha256']
    source = root / 'NSCLC_gene_disjoint_receiver_TF_preflight_20261016_v1'
    assert sha(source / 'private_independent_audit.json') == plan['scope_audit_sha256']
    assert sha(source / 'private_all_deduplicated_TF_slots.tsv') == plan['scope_sha256']
    assert sha(source / 'plan_before_any_TF_activity_or_association.json') == plan['source_scope_plan_sha256']
    scope_plan = json.loads((source / 'plan_before_any_TF_activity_or_association.json').read_text())
    assert sha(root / 'NSCLC_candidate_lineage_cofactor_RNA_20261016_v1/private_independent_audit.json') == plan['native_cofactor_audit_sha256']
    original_scope = pd.read_csv(source / 'private_all_deduplicated_TF_slots.tsv', sep='\t').sort_values(['record', 'TF']).reset_index(drop=True)
    scope = pd.read_csv(out / 'private_all_TF_measurement_scope.tsv', sep='\t')
    pd.testing.assert_frame_equal(original_scope, scope[original_scope.columns], check_dtype=False)
    assert scope.TF_record.tolist() == [f'TFD{i:05d}' for i in range(1, len(scope) + 1)]
    scores = pd.read_csv(out / 'private_all_matched_case_TF_scores.tsv.gz', sep='\t')
    assert not scores.duplicated(['TF_record', 'case_index']).any()
    hashes = pd.read_csv(out / 'private_ready_design_SHA256.tsv', sep='\t').set_index('TF_record')
    ready = scope[scope.measurement_status.eq('READY_FOR_INDEPENDENT_MEASUREMENT_AUDIT')]
    assert set(hashes.index) == set(ready.TF_record)
    net_path = root / 'fullgene_receiver_target_code_v1/collectri_reference_v1/human_collectri_signed_complexes.tsv'
    params_path = root / 'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv'
    official = root / 'fullgene_receiver_target_code_v1/collectri_reference_v1/official_decoupler192_ULM_source.txt'
    assert sha(net_path) == plan['signed_TF_network_sha256'] and sha(params_path) == plan['private_frozen_parameters_sha256']
    assert sha(official) == scope_plan['official_ULM_source_sha256']
    # Execute only the three score algebra functions from the pinned public file.
    tree = ast.parse(official.read_text())
    functions = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in {'mat_cov', 'mat_cor', 't_val'}]
    assert {f.name for f in functions} == {'mat_cov', 'mat_cor', 't_val'}
    ns = {'np': np}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(official), 'exec'), ns)
    net = pd.read_csv(net_path, sep='\t')
    params = pd.read_csv(params_path, sep='\t')
    parents = pd.read_csv(root / 'NSCLC_automatic_ligand_response_results_20261016_v1/private_all_unique_ligand_response_results.tsv', sep='\t').set_index('record')
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    genes = pd.Index(pd.read_csv(native / 'private_native_gene_axis.tsv', sep='\t').gene)
    matrix = np.load(native / 'mean_logCP10K.npy', mmap_mode='r')
    old_root = root / 'NSCLC_automatic_ligand_response_preflight_20261016_v1/private_frozen_designs'
    maximum_error, checked_scores, checked_designs = 0., 0, 0
    for parent_id, all_rows in scope.groupby('record', sort=True):
        active = all_rows[all_rows.status.eq('READY_FOR_INPUT_AUDIT_NO_ACTIVITY_COMPUTED')]
        if not len(active):
            continue
        parent = parents.loc[parent_id]
        pp = params[params.context.eq(parent.context) & params.program.eq(parent.program)]
        forbidden = set(pp.gene) | {parent.ligand}
        forbidden.update(g for receptor in parent.qualified_receptors.split('|') for g in receptor.split('_'))
        bg = pd.Index([g for g in genes if g not in forbidden])
        z = np.load(old_root / (parent_id + '.npz'), allow_pickle=False)
        assert sha(old_root / (parent_id + '.npz')) == active.original_design_sha256.iloc[0]
        expression = np.asarray(matrix[np.ix_(z['receiver_units'].astype(int), genes.get_indexer(bg))], dtype=np.float64)
        # Dense official algebra differs from producer's sparse covariance route.
        dense = np.zeros((len(bg), len(active)), dtype=np.float64)
        for j, tf in enumerate(active.TF):
            f = net[net.source.eq(tf) & net.target.isin(bg)]
            assert len(f) >= 20
            dense[bg.get_indexer(f.target), j] = f.weight.to_numpy(float)
        expected = ns['t_val'](ns['mat_cor'](dense, expression.T), len(bg) - 2)
        assert expected.shape == (len(z['cases']), len(active))
        for j, row in enumerate(active.itertuples()):
            y = expected[:, j]
            f = scores[scores.TF_record.eq(row.TF_record)].sort_values('case_index')
            assert len(f) == len(y) and f.case_index.tolist() == list(range(len(y)))
            assert f.case.astype(str).tolist() == z['cases'].astype(str).tolist()
            assert f.receiver_unit.tolist() == z['receiver_units'].tolist()
            assert f.record.eq(parent_id).all() and f.TF.eq(row.TF).all()
            assert np.allclose(y, f.ULM_score, atol=1e-8, rtol=1e-8, equal_nan=True)
            maximum_error = max(maximum_error, float(np.nanmax(np.abs(y - f.ULM_score.to_numpy()))))
            checked_scores += len(y)
            status = 'READY_FOR_INDEPENDENT_MEASUREMENT_AUDIT'
            if not np.isfinite(y).all():
                status = 'NONFINITE_TF_SCORE_NOT_NEGATIVE'
            elif np.ptp(y) <= 1e-12:
                status = 'NO_TF_SCORE_VARIATION_NOT_NEGATIVE'
            else:
                x = rankdata(z['X'])
                x = (x - x.mean()) / x.std(ddof=1)
                d = np.column_stack([z['Z'], x])
                if np.linalg.matrix_rank(d) != 20 or np.linalg.cond(d) >= 1e8 or len(y) < 50 or len(y) - 20 < 10:
                    status = 'ORIGINAL_REQUIRED_DESIGN_UNIDENTIFIABLE'
                elif np.max(np.sum(d * np.linalg.pinv(d).T, axis=1)) >= .999999:
                    status = 'ORIGINAL_REQUIRED_DESIGN_EXTREME_LEVERAGE'
            assert row.measurement_status == status
            if status == 'READY_FOR_INDEPENDENT_MEASUREMENT_AUDIT':
                path = out / 'private_frozen_designs' / (row.TF_record + '.npz')
                assert sha(path) == hashes.loc[row.TF_record, 'sha256']
                new = np.load(path, allow_pickle=False)
                for key in ['X', 'Z', 'covariates', 'cases', 'sender_units', 'receiver_units', 'sender_programme', 'response']:
                    assert np.array_equal(new[key], z[key])
                assert np.allclose(new['Y'], y, rtol=1e-8, atol=1e-8)
                assert np.array_equal(new['original_masked_receiver_programme'], z['Y'])
                assert np.array_equal(new['original_unmasked_receiver_programme'], z['receiver_full_programme'])
                checked_designs += 1
    unavailable = scope[~scope.status.eq('READY_FOR_INPUT_AUDIT_NO_ACTIVITY_COMPUTED')]
    assert unavailable.measurement_status.tolist() == unavailable.status.tolist()
    assert checked_scores == len(scores) and checked_designs == len(hashes)
    assert terminal['status_counts'] == scope.measurement_status.value_counts().to_dict()
    result = dict(status='PASS_ALL_DISJOINT_TF_SCORES_AGAINST_PINNED_OFFICIAL_ULM_ALL_MASKS_ALL_CASE_DESIGNS_AND_UNAVAILABLE_SLOTS', finished_UTC=datetime.now(timezone.utc).isoformat(), audit_code_sha256=sha(__file__), producer_terminal_sha256=sha(out / 'private_terminal.json'), complete_TF_scope=len(scope), all_measured_score_case_records=checked_scores, all_ready_designs=checked_designs, maximum_ULM_score_error=maximum_error, no_ligand_association_P_BH_computed=True, no_learned_programme_gene_lists_exported=True, independent_samples_or_protein_validation=False)
    dest.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
