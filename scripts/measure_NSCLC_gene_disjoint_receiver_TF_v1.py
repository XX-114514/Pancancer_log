"""Measure all frozen disjoint TF inputs, with no ligand association testing."""
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


def wait_audit(root, folder, job):
    receipt = json.loads((root / 'HNSCC_receiver_followup_job_control_20261016_v2' / job / 'launch_receipt.json').read_text())
    deadline = time.monotonic() + 14400
    path = root / folder / 'private_independent_audit.json'
    while not path.exists():
        proc = Path('/proc') / str(receipt['pid']) / 'stat'
        assert proc.exists()
        st = proc.read_text().rsplit(')', 1)[1].split()
        assert st[19] == str(receipt['proc_start_ticks']) and st[0] != 'Z' and time.monotonic() < deadline
        time.sleep(15)
    audit = json.loads(path.read_text())
    assert audit['status'].startswith('PASS_')
    return path


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    scope_root = root / 'NSCLC_gene_disjoint_receiver_TF_preflight_20261016_v1'
    scope_audit = wait_audit(root, scope_root.name, 'audit_NSCLC_gene_disjoint_receiver_TF_scope_v1')
    native_audit = wait_audit(root, 'NSCLC_candidate_lineage_cofactor_RNA_20261016_v1', 'audit_NSCLC_candidate_lineage_cofactor_RNA_v1')
    scope = pd.read_csv(scope_root / 'private_all_deduplicated_TF_slots.tsv', sep='\t').sort_values(['record', 'TF']).reset_index(drop=True)
    parents_root = root / 'NSCLC_automatic_ligand_response_results_20261016_v1'
    parents = pd.read_csv(parents_root / 'private_all_unique_ligand_response_results.tsv', sep='\t').set_index('record')
    scope_plan = json.loads((scope_root / 'plan_before_any_TF_activity_or_association.json').read_text())
    params_path = root / 'NSCLC_context_frozen_parameters_20261016_v1/private_all53_frozen_gene_parameters.tsv'
    net_path = root / 'fullgene_receiver_target_code_v1/collectri_reference_v1/human_collectri_signed_complexes.tsv'
    assert sha(params_path) == scope_plan['private_frozen_parameters_sha256'] and sha(net_path) == scope_plan['signed_TF_network_sha256']
    params = pd.read_csv(params_path, sep='\t')
    net = pd.read_csv(net_path, sep='\t')
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    genes = pd.Index(pd.read_csv(native / 'private_native_gene_axis.tsv', sep='\t').gene)
    matrix = np.load(native / 'mean_logCP10K.npy', mmap_mode='r')
    assert matrix.shape[1] == len(genes)
    out = root / 'NSCLC_gene_disjoint_receiver_TF_measurement_20261016_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    designs = out / 'private_frozen_designs'
    designs.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), scope_audit_sha256=sha(scope_audit), native_cofactor_audit_sha256=sha(native_audit), scope_sha256=sha(scope_root / 'private_all_deduplicated_TF_slots.tsv'), source_scope_plan_sha256=sha(scope_root / 'plan_before_any_TF_activity_or_association.json'), signed_TF_network_sha256=sha(net_path), private_frozen_parameters_sha256=sha(params_path),
        measurement='Official decoupler1.9.2 ULM-equivalent signed gene-pattern score using ALL native background genes outside both frozen programmes and all original interface components. Sparse dot-products, double precision. TF-score t uses gene-background degrees of freedom, not case-association P.',
        case_design='Original matched units/cases, sender ligand X, original full Z and all18 nuisance variables; no auto deletion. Require finite variable Y and full20-column design, condition<1e8, leverage<.999999 and original case/df gates.',
        scope='All original deduplicated parent/TF slots including coverage failures. No original LR score, endpoint, P, BH or exclusion changes.',
        independent_check='Next audit must compare ALL TF scores against official source algorithm and all original case/ligand/nuisance arrays; measurement alone is not acceptance.',
        no_ligand_TF_association_or_P_or_BH=True, same_source_postselection_exploratory=True, independent_samples_or_protein_validation=False)
    (out / 'plan_before_any_TF_measurement.json').write_text(json.dumps(plan, indent=2) + '\n')
    scope['TF_record'] = [f'TFD{i:05d}' for i in range(1, len(scope) + 1)]
    scope['measurement_status'] = scope.status
    score_rows, hashes = [], []
    old_designs = root / 'NSCLC_automatic_ligand_response_preflight_20261016_v1/private_frozen_designs'
    for parent_id, all_slots in scope.groupby('record', sort=True):
        active = all_slots[all_slots.status.eq('READY_FOR_INPUT_AUDIT_NO_ACTIVITY_COMPUTED')]
        if not len(active):
            continue
        parent = parents.loc[parent_id]
        path = old_designs / (parent_id + '.npz')
        assert sha(path) == active.original_design_sha256.iloc[0]
        z = np.load(path, allow_pickle=False)
        pp = params[params.context.eq(parent.context) & params.program.eq(parent.program)]
        forbidden = set(pp.gene) | {parent.ligand} | {g for r in parent.qualified_receptors.split('|') for g in r.split('_')}
        bg = pd.Index([g for g in genes if g not in forbidden])
        assert len(bg) == int(active.native_disjoint_background.iloc[0])
        expression = np.asarray(matrix[np.ix_(z['receiver_units'].astype(int), genes.get_indexer(bg))], dtype=np.float64)
        centered = expression - expression.mean(axis=1, keepdims=True)
        sd_expression = np.sqrt(np.einsum('ij,ij->i', centered, centered) / (len(bg) - 1))
        tfs = pd.Index(active.TF)
        f = net[net.source.isin(tfs) & net.target.isin(bg)]
        weights = csr_matrix((f.weight.to_numpy(float), (tfs.get_indexer(f.source), bg.get_indexer(f.target))), shape=(len(tfs), len(bg)))
        sum_weights = np.asarray(weights.sum(axis=1)).ravel()
        sum_squared = np.asarray(weights.multiply(weights).sum(axis=1)).ravel()
        sd_weights = np.sqrt((sum_squared - sum_weights ** 2 / len(bg)) / (len(bg) - 1))
        correlations = (weights @ centered.T) / ((len(bg) - 1) * sd_weights[:, None] * sd_expression[None, :])
        scores = correlations * np.sqrt((len(bg) - 2) / ((1 - correlations + 1e-16) * (1 + correlations + 1e-16)))
        assert scores.shape == (len(tfs), len(z['cases']))
        for j, row in enumerate(active.itertuples()):
            y = scores[j]
            status = 'READY_FOR_INDEPENDENT_MEASUREMENT_AUDIT'
            if not np.isfinite(y).all():
                status = 'NONFINITE_TF_SCORE_NOT_NEGATIVE'
            elif np.ptp(y) <= 1e-12:
                status = 'NO_TF_SCORE_VARIATION_NOT_NEGATIVE'
            else:
                X = rankdata(z['X'])
                X = (X - X.mean()) / X.std(ddof=1)
                D = np.column_stack([z['Z'], X])
                if np.linalg.matrix_rank(D) != 20 or np.linalg.cond(D) >= 1e8 or len(y) < 50 or len(y) - 20 < 10:
                    status = 'ORIGINAL_REQUIRED_DESIGN_UNIDENTIFIABLE'
                elif np.max(np.sum(D * np.linalg.pinv(D).T, axis=1)) >= .999999:
                    status = 'ORIGINAL_REQUIRED_DESIGN_EXTREME_LEVERAGE'
            scope.loc[scope.TF_record.eq(row.TF_record), 'measurement_status'] = status
            for i, case in enumerate(z['cases']):
                score_rows.append(dict(TF_record=row.TF_record, record=parent_id, TF=row.TF, case_index=i, case=str(case), receiver_unit=int(z['receiver_units'][i]), ULM_score=float(y[i])))
            if status == 'READY_FOR_INDEPENDENT_MEASUREMENT_AUDIT':
                target = designs / (row.TF_record + '.npz')
                np.savez_compressed(target, X=z['X'], Y=y, Z=z['Z'], covariates=z['covariates'], cases=z['cases'], sender_units=z['sender_units'], receiver_units=z['receiver_units'], sender_programme=z['sender_programme'], original_masked_receiver_programme=z['Y'], original_unmasked_receiver_programme=z['receiver_full_programme'], response=z['response'])
                hashes.append(dict(TF_record=row.TF_record, parent_record=parent_id, sha256=sha(target)))
    scope.to_csv(out / 'private_all_TF_measurement_scope.tsv', sep='\t', index=False)
    pd.DataFrame(score_rows).to_csv(out / 'private_all_matched_case_TF_scores.tsv.gz', sep='\t', index=False, compression='gzip')
    pd.DataFrame(hashes).to_csv(out / 'private_ready_design_SHA256.tsv', sep='\t', index=False)
    names = ['private_all_TF_measurement_scope.tsv', 'private_all_matched_case_TF_scores.tsv.gz', 'private_ready_design_SHA256.tsv', 'plan_before_any_TF_measurement.json']
    result = dict(status='COMPLETE_ALL_FIXED_DISJOINT_TF_MEASUREMENTS_INPUT_AUDIT_PENDING_NO_ASSOCIATION_TESTS', finished_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__), complete_TF_scope=len(scope), measured_TF_slots=int(scope.status.eq('READY_FOR_INPUT_AUDIT_NO_ACTIVITY_COMPUTED').sum()), score_case_records=len(score_rows), ready_designs=len(hashes), status_counts=scope.measurement_status.value_counts().to_dict(), files={name: sha(out / name) for name in names}, new_ligand_association_P_BH=False, independent_validation=False)
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
