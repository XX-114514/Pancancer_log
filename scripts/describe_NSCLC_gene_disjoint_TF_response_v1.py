"""All fixed TF response points and complete case LOO; no new P or BH.

These same-source postselection descriptions are not independent confirmation.
"""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import qr
from scipy.stats import rankdata


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8388608), b''):
            h.update(b)
    return h.hexdigest()


def correlations(x, Y, Z):
    x = rankdata(x).astype(float)
    y = rankdata(Y, axis=0).astype(float)
    if np.linalg.matrix_rank(Z) != Z.shape[1]:
        return np.full(y.shape[1], np.nan)
    D = np.column_stack([Z, (x - x.mean()) / x.std()]) if x.std() > 1e-12 else None
    if D is None or np.linalg.matrix_rank(D) != D.shape[1] or np.linalg.cond(D) >= 1e8:
        return np.full(y.shape[1], np.nan)
    Q = qr(Z, mode='economic')[0]
    x = x - Q @ (Q.T @ x)
    y = y - Q @ (Q.T @ y)
    denominator = np.linalg.norm(x) * np.linalg.norm(y, axis=0)
    return np.divide(x @ y, denominator, out=np.full(y.shape[1], np.nan), where=denominator > 1e-12)


def extra_design(Z, extra):
    sd = np.std(extra)
    return np.column_stack([Z, (extra - np.mean(extra)) / sd]) if sd > 1e-12 else None


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inputs = root / 'NSCLC_gene_disjoint_receiver_TF_measurement_20261016_v1'
    audit_path = inputs / 'private_independent_audit.json'
    audit = json.loads(audit_path.read_text())
    assert audit['status'].startswith('PASS_ALL_DISJOINT_TF_SCORES_')
    terminal_path = inputs / 'private_terminal.json'
    assert sha(terminal_path) == audit['producer_terminal_sha256']
    terminal = json.loads(terminal_path.read_text())
    scope_path = inputs / 'private_all_TF_measurement_scope.tsv'
    hashes_path = inputs / 'private_ready_design_SHA256.tsv'
    for path in [scope_path, hashes_path]:
        assert sha(path) == terminal['files'][path.name]
    scope = pd.read_csv(scope_path, sep='\t')
    hashes = pd.read_csv(hashes_path, sep='\t').set_index('TF_record')
    ready = scope[scope.measurement_status.eq('READY_FOR_INDEPENDENT_MEASUREMENT_AUDIT')]
    assert set(ready.TF_record) == set(hashes.index)
    parent_root = root / 'NSCLC_automatic_ligand_response_results_20261016_v1'
    parents = pd.read_csv(parent_root / 'private_all_unique_ligand_response_results.tsv', sep='\t').set_index('record')
    out = root / 'NSCLC_gene_disjoint_TF_response_descriptions_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), input_audit_sha256=sha(audit_path), input_scope_sha256=sha(scope_path), input_hash_table_sha256=sha(hashes_path),
        family='ALL2568 frozen parent/TF slots, including all998 original coverage failures; no top-TF, effect or receptor-sign selection.',
        main='Rank original sender ligand and native gene-disjoint receiver TF score; partial rank correlation with ALL original18 raw standardized nuisance variables and intercept. No new P, BH, intervals or null simulation.',
        sensitivities='Add original raw frozen sender programme or pathology response separately, requiring identifiable full extra design. Pathology control may induce collider bias. Descriptive TF-to-original-masked-MCP correlation with the same full Z.',
        LOO='ALL matched author-case deletions, rerank X and every TF Y, retain full original nuisance set/rank/condition gates; failed cases and undefined residual TF scores stay NA.',
        reference_sign='Only a signed receptor-to-TF reference path, NOT proof that the ligand activates this receptor; all cross-receptor sign ambiguity remains.',
        interpretation='Same-source postselection exploration. Gene-disjoint from original programmes/interfaces, not independent samples or all revised biochemical alternatives. No original LR P/q or molecular caveat changed.')
    (out / 'plan_before_any_TF_association_description.json').write_text(json.dumps(plan, indent=2) + '\n')
    rows, loo_rows = [], []
    for parent_id, f in ready.groupby('record', sort=True):
        zs = []
        for row in f.itertuples():
            path = inputs / 'private_frozen_designs' / (row.TF_record + '.npz')
            assert sha(path) == hashes.loc[row.TF_record, 'sha256']
            zs.append(np.load(path, allow_pickle=False))
        z = zs[0]
        x, Z = z['X'], z['Z']
        for zz in zs:
            for name in ['X', 'Z', 'cases', 'sender_units', 'receiver_units', 'sender_programme', 'response', 'original_masked_receiver_programme']:
                assert np.array_equal(zz[name], z[name])
        Y = np.column_stack([zz['Y'] for zz in zs])
        main_rho = correlations(x, Y, Z)
        raw = correlations(x, Y, np.ones((len(x), 1)))
        linkage = correlations(z['original_masked_receiver_programme'], Y, Z)
        extra = {}
        for name, value in [('sender_programme', z['sender_programme']), ('pathology', z['response'])]:
            design = extra_design(Z, value)
            extra[name] = correlations(x, Y, design) if design is not None else np.full(Y.shape[1], np.nan)
        lo = np.full((len(x), Y.shape[1]), np.nan)
        for k in range(len(x)):
            take = np.arange(len(x)) != k
            lo[k] = correlations(x[take], Y[take], Z[take])
        for j, row in enumerate(f.itertuples()):
            valid = np.isfinite(lo[:, j])
            signs = np.sign(lo[valid, j]) == np.sign(main_rho[j])
            parent = parents.loc[parent_id]
            rows.append(dict(TF_record=row.TF_record, status='DESCRIPTIVE_FULL_CONTROL_POINT_ONLY' if np.isfinite(main_rho[j]) else 'UNDEFINED_FULL_CONTROL_RESIDUAL_SCORE', ligand=parent.ligand, cases=len(x), raw_spearman=float(raw[j]), partial_rank_rho=float(main_rho[j]), sender_programme_adjusted_rho=float(extra['sender_programme'][j]), pathology_adjusted_rho=float(extra['pathology'][j]), TF_to_original_masked_receiver_MCP_rho=float(linkage[j]), LOO_estimable=int(valid.sum()), LOO_same_sign=int(signs.sum()), LOO_rho_minimum=float(np.min(lo[valid, j])) if valid.any() else None, LOO_rho_maximum=float(np.max(lo[valid, j])) if valid.any() else None, no_new_P_BH=True, independent_samples=False))
            loo_rows.extend(dict(TF_record=row.TF_record, removed_case_index=k, rho=float(lo[k, j])) for k in range(len(x)))
        for zz in zs:
            zz.close()
        (out / 'private_progress.json').write_text(json.dumps(dict(observed_UTC=datetime.now(timezone.utc).isoformat(), completed_parent=parent_id, completed_TF_rows=len(rows))))
    points = pd.DataFrame(rows)
    table = scope.merge(points, on='TF_record', how='left', validate='one_to_one', suffixes=('', '_effect'))
    assert len(table) == len(scope) and len(points) == len(ready)
    table['description_status'] = table.status_effect.fillna(table.measurement_status)
    table.to_csv(out / 'private_all2568_TF_functional_descriptions.tsv', sep='\t', index=False, na_rep='NA')
    pd.DataFrame(loo_rows).to_csv(out / 'private_all_matched_case_LOO.tsv.gz', sep='\t', index=False, na_rep='NA', compression='gzip')
    files = ['private_all2568_TF_functional_descriptions.tsv', 'private_all_matched_case_LOO.tsv.gz', 'plan_before_any_TF_association_description.json']
    result = dict(status='COMPLETE_ALL_FIXED_TF_CONTROLLED_DESCRIPTIONS_AND_FULL_CASE_LOO_AUDIT_PENDING_NO_NEW_TESTS', finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), all_fixed_slots=len(table), described_slots=len(points), all_case_LOO_rows=len(loo_rows), status_counts=table.description_status.value_counts().to_dict(), files={name: sha(out / name) for name in files}, P_BH_or_independent_validation=False)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
