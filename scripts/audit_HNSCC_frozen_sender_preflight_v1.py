"""Audit all planned designs with independent keys, manual ranks and pivoted QR."""
import hashlib
import itertools
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import qr


def manual_rank(x):
    x = np.asarray(x, float)
    return (x[:, None] > x[None, :]).sum(axis=1) + ((x[:, None] == x[None, :]).sum(axis=1)+1)/2


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root/'HNSCC_frozen_sender_association_preflight_20261016_v1'
    target = out/'private_independent_audit.json'
    assert not target.exists()
    terminal = json.loads((out/'private_terminal.json').read_text())
    contract = json.loads((out/'contract_before_any_association_effects.json').read_text())
    design_path = out/'private_all_planned_patient_designs.tsv'
    result_path = out/'all400_private_design_identifiability.tsv'
    assert hashlib.sha256(design_path.read_bytes()).hexdigest() == terminal['design_sha256']
    assert hashlib.sha256(result_path.read_bytes()).hexdigest() == terminal['preflight_sha256']
    frame = pd.read_csv(design_path, sep='\t')
    result = pd.read_csv(result_path, sep='\t', keep_default_na=False)
    prior = root/'HNSCC_receiver_prior_ligand_ranking_20261016_v1'
    ranks = pd.read_csv(prior/'private_all_ligand_prior_rankings.tsv', sep='\t')
    selected = ranks[ranks.prospective_top20]
    senders = pd.read_csv(prior/'private_sender_capture_summary.tsv', sep='\t')
    senders = senders[senders.prospective_sender_eligible]
    expected = set(itertools.product(selected.ligand, senders.sender_group))
    assert len(expected) == 400 and expected == set(zip(result.ligand, result.sender_group))
    inputs = root/'HNSCC_all_sender_ligand_measurements_20261016_v1'
    rawmeta = pd.read_csv(inputs/'private_all_patient_group_capture_quality.tsv', sep='\t').set_index(['sender_group', 'patient_id'])
    expression = pd.read_csv(inputs/'private_all_patient_group_gene_expression.tsv', sep='\t').set_index(['sender_group', 'patient_id', 'gene'])
    receptors = pd.read_csv(inputs/'private_all_patient_receptor_complex_RNA.tsv', sep='\t')
    old = pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_patient_subtype_scores.tsv', sep='\t')
    old = old[old.subtype_candidate.eq('Macro_C1QC')].set_index('patient_id')
    author = pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv', sep='\t')
    author = author[author.time.eq('Pre')].set_index('patient_registry')
    numeric = contract['numeric_nuisance']; categorical = contract['categorical_nuisance']
    assert len(numeric) == 10 and categorical == ['D_plus_T', 'HPV_positive']
    for row in frame.itertuples():
        s = rawmeta.loc[(row.sender_group, row.patient_id)]; r = rawmeta.loc[('Receiver_C1QC', row.patient_id)]
        e = expression.loc[(row.sender_group, row.patient_id, row.ligand)]
        own = expression.loc[('Receiver_C1QC', row.patient_id, row.ligand)]
        rec = receptors[receptors.patient_id.eq(row.patient_id) & receptors.ligand.eq(row.ligand) & receptors.original_pooled_reference_and_RNA_qualified]
        expected_values = [e.mean_log1p_CPTT, old.loc[row.patient_id, 'score'], s.mean_logRNA, s.mean_loggenes,
            r.mean_logRNA, r.mean_loggenes, np.log1p(s.native_cells), np.log1p(r.native_cells),
            s.retained_TME_fraction, r.receiver_fraction_within_Myeloid, own.mean_log1p_CPTT,
            rec.same_cell_positive_fraction.max() if len(rec) else np.nan]
        observed = [row.exposure, row.response] + [getattr(row, field) for field in numeric]
        assert np.allclose(observed, expected_values, equal_nan=True, atol=1e-11, rtol=1e-10)
        assert row.D_plus_T == int(author.loc[row.patient_id, 'neo_type'] == 'D+T')
        assert row.HPV_positive == int(author.loc[row.patient_id, 'HPV/p16'] == 'positive')
        included = s.native_cells >= 20 and r.native_cells >= 20 and np.isfinite(observed).all()
        assert included == bool(row.included_in_full_design)
    positive = 0; invalid = 0
    for row in result.itertuples():
        sub = frame[frame.ligand.eq(row.ligand) & frame.sender_group.eq(row.sender_group)]
        assert len(sub) == len(old) and sub.patient_id.is_unique
        use = sub[sub.included_in_full_design]
        assert len(use) == int(row.complete_matched_patients)
        df = len(use)-14; assert df == int(row.residual_df)
        constants = []
        full_rank = 0; condition = np.nan
        if len(use):
            z = []
            for field in numeric+['exposure']:
                rr = manual_rank(use[field])
                if np.std(rr) < 1e-12:
                    constants.append(field); z.append(np.zeros(len(use)))
                else:
                    z.append((rr-np.mean(rr))/np.std(rr))
            D = np.column_stack([np.ones(len(use)), use[categorical].to_numpy(float), np.column_stack(z)])
            _, R, _ = qr(D, mode='economic', pivoting=True)
            tol = max(D.shape)*np.finfo(float).eps*abs(R).max()
            full_rank = int(np.count_nonzero(np.abs(np.diag(R)) > tol))
            singular = np.linalg.svd(D, compute_uv=False)
            condition = float(singular[0]/singular[-1])
        assert full_rank == int(row.full_design_rank)
        assert '|'.join(constants) == row.constant_fields
        expected_ok = bool(row.reference_and_RNA_qualified) and len(use) >= 20 and df >= 10 and not constants and full_rank == 14 and np.isfinite(condition) and condition <= 1e8
        assert expected_ok == (row.status == 'full_design_estimable_effects_pending')
        positive += int(expected_ok); invalid += int(not expected_ok)
    assert positive == terminal['estimable_slots'] and positive+invalid == 400
    report = dict(status='PASS_ALL400_KEYS_NATIVE_COVARIATES_MANUAL_RANKS_AND_PIVOTED_QR_IDENTIFIABILITY',
        finished_UTC=datetime.now(timezone.utc).isoformat(), planned_slots=400, estimable_slots=positive,
        design_sha256=terminal['design_sha256'], preflight_sha256=terminal['preflight_sha256'],
        no_effects_or_covariate_deletion=True)
    target.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report))


if __name__ == '__main__':
    main()
