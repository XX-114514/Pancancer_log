"""Freeze full patient-level sender/receiver design before any effect tests."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

NUMERIC = ['sender_logRNA', 'sender_loggenes', 'receiver_logRNA', 'receiver_loggenes',
           'sender_logcells', 'receiver_logcells', 'sender_fraction_TME',
           'receiver_fraction_Myeloid', 'receiver_own_ligand', 'receiver_joint_RNA_max']
CATEGORICAL = ['D_plus_T', 'HPV_positive']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    inputs = root/'HNSCC_all_sender_ligand_measurements_20261016_v1'
    a = json.loads((inputs/'private_independent_audit.json').read_text())
    assert a['status'].startswith('PASS_ALL400_LIGAND_SLOTS')
    for name, key in [('private_all_patient_group_gene_expression.tsv', 'expression_sha256'),
                      ('private_all_patient_group_capture_quality.tsv', 'quality_sha256'),
                      ('all400_private_sender_ligand_coverage.tsv', 'coverage_sha256')]:
        assert sha(inputs/name) == a[key]
    out = root/'HNSCC_frozen_sender_association_preflight_20261016_v1'
    assert not out.exists()
    out.mkdir()
    contract = dict(created_UTC=datetime.now(timezone.utc).isoformat(), producer_sha256=sha(Path(__file__)),
        inputs_audit_sha256=sha(inputs/'private_independent_audit.json'),
        hypothesis='Native sender ligand mean versus fixed original C1QC full macrophage MCP score; no clinical outcome re-screen.',
        primary_family_slots=400, original_sign_slots='Forty prospectively frozen ligands, no replacements, all ten capture-eligible sender groups.',
        numeric_nuisance=NUMERIC, categorical_nuisance=CATEGORICAL,
        numeric_transform='Average-tie rank, zero mean and population SD one within each matched patient set. Rank transformation recomputed in resampling.',
        exposure='Mean log1p(10000*ligand/native full RNA) across all native sender cells; measured zeros retained.',
        response='Original full fixed partial-transfer MCP from original deterministic C1QC cell draws, unchanged; do not switch to a favorable signed component.',
        abundance='Both native cell counts; sender full retained TME fraction, receiver C1QC fraction within Myeloid. Major-level sender within-major fraction would be identically one; full-TME sender fraction is defined prospectively instead.',
        receptor_nuisance='Maximum patient same-cell RNA-positive fraction across that ligand original reference-and-pooled-RNA-qualified receptor entries. Fixed measurement definition, not chosen by outcome or effect; not complex protein abundance or a receptor-specific prior.',
        own_ligand='Mean ligand expression in receiver cells; absence/constancy is an explicit identifiability limitation, never an automatically removed covariate.',
        minimum_sender_and_receiver_cells=20, minimum_matched_patients=20, minimum_residual_df=10,
        constant_or_collinear_required_controls='Unestimable; retain planned slot, do not auto-delete covariates.',
        maximum_scaled_design_condition=1e8,
        reference_gate='At least one original constituent-reference and pooled same-cell receptor RNA-qualified entry; no binding claim.',
        future_inference='If proceeding, all400 primary slots retain internal P1 for unavailable/unestimable tests and one BH family; display missing. Twenty-thousand residual permutations within treatment x HPV, two-thousand patient bootstrap within same strata, fixed seeded streams, full-LOO and interface-disjoint sensitivity before interpreting an axis.',
        permutation_limitation='Conditional residual exchangeability is assumed, not established by this preflight. Singleton strata retained without moving that residual.',
        no_effect_P_BH=True, outcome_selected_receiver_is_exploratory=True)
    (out/'contract_before_any_association_effects.json').write_text(json.dumps(contract, indent=2)+'\n')
    expression = pd.read_csv(inputs/'private_all_patient_group_gene_expression.tsv', sep='\t')
    meta = pd.read_csv(inputs/'private_all_patient_group_capture_quality.tsv', sep='\t').set_index(['sender_group', 'patient_id'])
    coverage = pd.read_csv(inputs/'all400_private_sender_ligand_coverage.tsv', sep='\t')
    receptor = pd.read_csv(inputs/'private_all_patient_receptor_complex_RNA.tsv', sep='\t')
    old = pd.read_csv(root/'HNSCC_within_macro_subtype_diagnostics_20261014_v1/private_patient_subtype_scores.tsv', sep='\t')
    old = old[old.subtype_candidate.eq('Macro_C1QC')].set_index('patient_id').sort_index()
    author = pd.read_csv(root/'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv', sep='\t')
    author = author[author.time.eq('Pre')].set_index('patient_registry')
    assert old.index.is_unique and author.index.is_unique
    assert not expression.duplicated(['sender_group', 'patient_id', 'gene']).any()
    expr = expression.set_index(['sender_group', 'patient_id', 'gene'])
    subjects = old.index.tolist()
    rows = []; all_frames = []
    for slot in coverage.itertuples():
        data = []
        for patient in subjects:
            sender = meta.loc[(slot.sender_group, patient)]; rx = meta.loc[('Receiver_C1QC', patient)]
            e = expr.loc[(slot.sender_group, patient, slot.ligand)]
            own = expr.loc[('Receiver_C1QC', patient, slot.ligand)]
            rr = receptor[receptor.patient_id.eq(patient) & receptor.ligand.eq(slot.ligand) & receptor.original_pooled_reference_and_RNA_qualified]
            cm = author.loc[patient]
            assert int(rx.native_cells) == int(old.loc[patient, 'native_subtype_cells'])
            data.append(dict(patient_id=patient, ligand=slot.ligand, sender_group=slot.sender_group,
                programme_sign=slot.programme_sign, sender_native_cells=int(sender.native_cells), receiver_native_cells=int(rx.native_cells),
                exposure=e.mean_log1p_CPTT, response=old.loc[patient, 'score'],
                sender_logRNA=sender.mean_logRNA, sender_loggenes=sender.mean_loggenes,
                receiver_logRNA=rx.mean_logRNA, receiver_loggenes=rx.mean_loggenes,
                sender_logcells=np.log1p(sender.native_cells), receiver_logcells=np.log1p(rx.native_cells),
                sender_fraction_TME=sender.retained_TME_fraction, receiver_fraction_Myeloid=rx.receiver_fraction_within_Myeloid,
                receiver_own_ligand=own.mean_log1p_CPTT,
                receiver_joint_RNA_max=rr.same_cell_positive_fraction.max() if len(rr) else np.nan,
                D_plus_T=int(cm.neo_type == 'D+T'), HPV_positive=int(cm['HPV/p16'] == 'positive')))
        frame = pd.DataFrame(data)
        count_ok = frame.sender_native_cells.ge(20) & frame.receiver_native_cells.ge(20)
        complete = frame[['exposure', 'response']+NUMERIC+CATEGORICAL].notna().all(axis=1)
        frame['included_in_full_design'] = count_ok & complete
        use = frame[frame.included_in_full_design]
        reasons = []
        if not slot.receiver_prior_qualified:
            reasons.append('receiver_reference_or_pooled_RNA_not_qualified')
        if len(use) < 20:
            reasons.append('fewer_than20_complete_matched_patients')
        df = len(use)-(len(NUMERIC)+len(CATEGORICAL)+2)
        if df < 10:
            reasons.append('fewer_than10_residual_degrees_of_freedom')
        constants = []
        rank = 0; condition = np.nan; design_cols = len(NUMERIC)+len(CATEGORICAL)+2
        if len(use):
            values = []
            for field in NUMERIC+['exposure']:
                r = rankdata(use[field])
                if r.std() < 1e-12:
                    constants.append(field); values.append(np.zeros(len(use)))
                else:
                    values.append((r-r.mean())/r.std())
            values = np.column_stack(values)
            D = np.column_stack((np.ones(len(use)), use[CATEGORICAL].to_numpy(float), values))
            rank = int(np.linalg.matrix_rank(D)); condition = float(np.linalg.cond(D))
            if constants:
                reasons.append('constant_required_continuous_fields')
            if rank != design_cols:
                reasons.append('full_design_rank_deficient')
            if not np.isfinite(condition) or condition > 1e8:
                reasons.append('scaled_design_condition_exceeds1e8')
        status = 'full_design_estimable_effects_pending' if not reasons else 'unestimable_or_ineligible'
        rows.append(dict(programme_sign=slot.programme_sign, ligand=slot.ligand, sender_group=slot.sender_group,
            reference_and_RNA_qualified=bool(slot.receiver_prior_qualified), captured_patients=int(count_ok.sum()),
            complete_matched_patients=len(use), residual_df=df, design_columns=design_cols, full_design_rank=rank,
            scaled_condition=condition, constant_fields='|'.join(constants), reasons='|'.join(reasons), status=status,
            patients_with_any_ligand=slot.patients_with_any_positive_ligand_cells, no_effect_tests=True))
        all_frames.append(frame)
    results = pd.DataFrame(rows)
    assert len(results) == 400
    results.to_csv(out/'all400_private_design_identifiability.tsv', sep='\t', index=False, na_rep='NA')
    pd.concat(all_frames, ignore_index=True).to_csv(out/'private_all_planned_patient_designs.tsv', sep='\t', index=False, na_rep='NA')
    terminal = dict(status='ALL400_PROSPECTIVE_FULL_COVARIATE_DESIGNS_PREFLIGHTED_NO_EFFECTS',
        finished_UTC=datetime.now(timezone.utc).isoformat(), planned_slots=len(results),
        estimable_slots=int(results.status.eq('full_design_estimable_effects_pending').sum()),
        status_counts=results.status.value_counts().to_dict(),
        reason_counts=results.reasons.str.split('|', regex=False).explode().value_counts().to_dict(),
        constant_control_counts=results.constant_fields.str.split('|', regex=False).explode().value_counts().to_dict(),
        design_sha256=sha(out/'private_all_planned_patient_designs.tsv'),
        preflight_sha256=sha(out/'all400_private_design_identifiability.tsv'),
        candidate_family_unchanged=True, no_association_or_P=True)
    (out/'private_terminal.json').write_text(json.dumps(terminal, indent=2)+'\n')
    print(json.dumps(terminal))
    print(coverage.groupby('sender_group').agg(slots=('ligand', 'size'),
        any_in20_patients=('patients_with_any_positive_ligand_cells', lambda s: int(s.ge(20).sum())),
        any_in_any_patient=('patients_with_any_positive_ligand_cells', lambda s: int(s.gt(0).sum()))).to_json(orient='index'))
    print(coverage[coverage.ligand.str.startswith('IFN')].groupby('ligand').agg(
        max_positive_patients=('patients_with_any_positive_ligand_cells', 'max'),
        max_sender_positive_fraction=('pooled_sender_positive_fraction', 'max')).to_json(orient='index'))


if __name__ == '__main__':
    main()
