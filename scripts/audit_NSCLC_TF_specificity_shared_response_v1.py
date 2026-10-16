"""All source signature scores, RNA diagnostics and target/score redundancies."""
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.stats import rankdata


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def correlation(x, y, Z):
    a, b = rankdata(x).astype(float), rankdata(y).astype(float)
    if np.std(a) <= 1e-12 or np.std(b) <= 1e-12:
        return np.nan
    D = np.column_stack([Z, (a - np.mean(a)) / np.std(a)])
    if np.linalg.matrix_rank(Z) != Z.shape[1] or np.linalg.matrix_rank(D) != D.shape[1] or np.linalg.cond(D) >= 1e8:
        return np.nan
    m = np.column_stack([a, b])
    e = m - Z @ np.linalg.lstsq(Z, m, rcond=None)[0]
    norm = np.linalg.norm(e[:, 0]) * np.linalg.norm(e[:, 1])
    return float(np.dot(e[:, 0], e[:, 1]) / norm) if norm > 1e-12 else np.nan


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'NSCLC_TF_specificity_shared_response_20261016_v1'
    launch = json.loads((root / 'HNSCC_receiver_followup_job_control_20261016_v2/diagnose_NSCLC_TF_specificity_and_shared_response_v1/launch_receipt.json').read_text())
    deadline = time.monotonic() + 14400
    while not (out / 'private_terminal.json').exists():
        p = Path('/proc') / str(launch['pid']) / 'stat'
        assert p.exists()
        s = p.read_text().rsplit(')', 1)[1].split()
        assert s[19] == str(launch['proc_start_ticks']) and s[0] != 'Z' and time.monotonic() < deadline
        time.sleep(15)
    dest = out / 'private_independent_audit.json'
    assert not dest.exists()
    terminal = json.loads((out / 'private_terminal.json').read_text())
    for name, value in terminal['files'].items():
        assert sha(out / name) == value
    plan = json.loads((out / 'plan_before_specificity_diagnostics.json').read_text())
    assert plan['producer_sha256'] == launch['source_sha256']
    desc = root / 'NSCLC_gene_disjoint_TF_response_descriptions_20261016_v1'
    assert sha(desc / 'private_independent_audit.json') == plan['functional_audit_sha256']
    assert sha(desc / 'private_all2568_TF_functional_descriptions.tsv') == plan['functional_table_sha256']
    original = pd.read_csv(desc / 'private_all2568_TF_functional_descriptions.tsv', sep='\t')
    active = original[original.description_status.eq('DESCRIPTIVE_FULL_CONTROL_POINT_ONLY')]
    records = pd.read_csv(out / 'private_all_regulator_specificity_and_shared_response.tsv', sep='\t').set_index('TF_record')
    pairs = pd.read_csv(out / 'private_all_regulator_pair_redundancy.tsv.gz', sep='\t')
    pcs = pd.read_csv(out / 'private_all_parent_residual_PC_descriptions.tsv', sep='\t').set_index('record')
    source_scores = pd.read_csv(out / 'private_all_sender_same_regulator_case_scores.tsv.gz', sep='\t')
    assert not source_scores.duplicated(['TF_record', 'case_index']).any()
    assert set(records.index) == set(active.TF_record)
    assert not pairs.duplicated(['record', 'TF_A', 'TF_B']).any()
    inputs = root / 'NSCLC_gene_disjoint_receiver_TF_measurement_20261016_v1'
    net = pd.read_csv(root / 'fullgene_receiver_target_code_v1/collectri_reference_v1/human_collectri_signed_complexes.tsv', sep='\t')
    params = pd.read_csv(root / 'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv', sep='\t')
    parents = pd.read_csv(root / 'NSCLC_automatic_ligand_response_results_20261016_v1/private_all_unique_ligand_response_results.tsv', sep='\t').set_index('record')
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    genes = pd.Index(pd.read_csv(native / 'private_native_gene_axis.tsv', sep='\t').gene)
    units = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t')
    matrix = np.load(native / 'mean_logCP10K.npy', mmap_mode='r')
    positive = np.load(native / 'positive_cell_counts.npy', mmap_mode='r')
    score_error, point_error, checked_pairs, checked_cases = 0., 0., 0, 0
    for parent_id, f in active.groupby('record', sort=True):
        parent = parents.loc[parent_id]
        pp = params[params.context.eq(parent.context) & params.program.eq(parent.program)]
        forbidden = set(pp.gene) | {parent.ligand} | {g for rec in parent.qualified_receptors.split('|') for g in rec.split('_')}
        bg = pd.Index([g for g in genes if g not in forbidden])
        tf_names = pd.Index(f.TF)
        nn = net[net.source.isin(tf_names) & net.target.isin(bg)]
        W = csr_matrix((nn.weight.to_numpy(float), (tf_names.get_indexer(nn.source), bg.get_indexer(nn.target))), shape=(len(f), len(bg)))
        loaded = [np.load(inputs / 'private_frozen_designs' / (r.TF_record + '.npz'), allow_pickle=False) for r in f.itertuples()]
        z = loaded[0]
        X = np.asarray(matrix[np.ix_(z['sender_units'].astype(int), genes.get_indexer(bg))], float)
        centered = X - X.mean(axis=1, keepdims=True)
        wsum = np.asarray(W.sum(axis=1)).ravel()
        wss = np.asarray(W.multiply(W).sum(axis=1)).ravel()
        sw = np.sqrt((wss - wsum ** 2 / len(bg)) / (len(bg) - 1))
        sx = np.sqrt(np.einsum('ij,ij->i', centered, centered) / (len(bg) - 1))
        r = (W @ centered.T) / ((len(bg) - 1) * sw[:, None] * sx[None, :])
        sender_scores = (r * np.sqrt((len(bg) - 2) / ((1 - r + 1e-16) * (1 + r + 1e-16)))).T
        Y = np.column_stack([a['Y'] for a in loaded])
        R = rankdata(Y, axis=0).astype(float)
        E = R - z['Z'] @ np.linalg.lstsq(z['Z'], R, rcond=None)[0]
        E = E / np.linalg.norm(E, axis=0)
        score_cor = E.T @ E
        weighted = (W @ W.T).toarray() / np.sqrt(wss[:, None] * wss[None, :])
        eigval, eigvec = np.linalg.eigh(score_cor)
        v = eigvec[:, -1]
        if v[np.argmax(np.abs(v))] < 0:
            v = -v
        pc = E @ v
        xr = rankdata(z['X']).astype(float)
        xr -= z['Z'] @ np.linalg.lstsq(z['Z'], xr, rcond=None)[0]
        observed_pc = pcs.loc[parent_id]
        assert np.isclose(observed_pc.first_PC_variance_fraction, eigval[-1] / eigval.sum(), atol=1e-10, rtol=1e-10)
        assert np.isclose(observed_pc.residual_ligand_first_PC_correlation, np.dot(xr, pc) / (np.linalg.norm(xr) * np.linalg.norm(pc)), atol=1e-10, rtol=1e-10)
        pair_rows = pairs[pairs.record.eq(parent_id)].set_index(['TF_A', 'TF_B'])
        assert len(pair_rows) == len(f) * (len(f) - 1) // 2
        for j, row in enumerate(f.itertuples()):
            observed = records.loc[row.TF_record]
            scores = source_scores[source_scores.TF_record.eq(row.TF_record)].sort_values('case_index')
            assert scores.case_index.tolist() == list(range(len(X)))
            assert np.allclose(scores.sender_same_regulator_score, sender_scores[:, j], atol=1e-8, rtol=1e-8)
            score_error = max(score_error, float(np.max(np.abs(scores.sender_same_regulator_score.to_numpy() - sender_scores[:, j]))))
            sx = sender_scores[:, j]
            Z = z['Z']
            full = np.column_stack([Z, (sx - np.mean(sx)) / np.std(sx)]) if np.std(sx) > 1e-12 else None
            values = dict(sender_same_regulator_adjusted_rho=correlation(z['X'], Y[:, j], full) if full is not None else np.nan, sender_regulator_to_ligand_rho=correlation(z['X'], sx, Z), sender_to_receiver_same_regulator_rho=correlation(sx, Y[:, j], Z))
            for name, expected in values.items():
                assert np.isclose(observed[name], expected, atol=1e-8, rtol=1e-8, equal_nan=True)
                if np.isfinite(expected):
                    point_error = max(point_error, abs(observed[name] - expected))
            position = genes.get_indexer([row.TF])[0]
            assert observed.literal_regulator_gene_measured == (position >= 0)
            if position >= 0:
                for role in ['sender', 'receiver']:
                    indices = z[role + '_units'].astype(int)
                    values = positive[indices, position]
                    assert observed[role + '_regulator_RNA_positive_cells'] == int(np.sum(values))
                    assert observed[role + '_regulator_RNA_positive_cases'] == int(np.count_nonzero(values))
                    assert np.isclose(observed[role + '_regulator_RNA_positive_fraction'], np.sum(values) / units.n_cells.to_numpy()[indices].sum(), atol=1e-12, rtol=1e-12)
            for k in range(j + 1, len(f)):
                pair = pair_rows.loc[(tf_names[j], tf_names[k])]
                assert np.isclose(pair.signed_target_cosine, weighted[j, k], atol=1e-12, rtol=1e-12)
                assert np.isclose(pair.receiver_residual_rank_correlation, score_cor[j, k], atol=1e-10, rtol=1e-10)
                checked_pairs += 1
            checked_cases += len(X)
        for a in loaded:
            a.close()
    assert checked_cases == len(source_scores) and checked_pairs == len(pairs)
    result = dict(status='PASS_ALL_SENDER_REGULATOR_SCORES_RNA_COVERAGE_SHARED_CONTROL_POINTS_ALL_TARGET_SCORE_PAIRS_AND_PARENT_PC_DESCRIPTIONS', finished_UTC=datetime.now(timezone.utc).isoformat(), audit_code_sha256=sha(__file__), producer_terminal_sha256=sha(out / 'private_terminal.json'), regulators=len(records), all_sender_case_scores=checked_cases, all_regulator_pairs=checked_pairs, maximum_sender_score_error=score_error, maximum_adjusted_point_error=point_error, no_new_P_BH_or_expression_based_selection=True, independent_validation=False)
    dest.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
