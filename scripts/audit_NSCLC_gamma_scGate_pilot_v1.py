"""Complete pilot sampling/R-label/aggregate structural audit, not TCR truth."""
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); guard = resource(root, 1024**2)
    src = root / 'NSCLC_gamma_scGate_fixed_pilot_20261019_v1'
    out = root / 'NSCLC_gamma_scGate_fixed_pilot_audit_20261019_v1'; assert not out.exists(); out.mkdir(mode=0o700)
    t = wait(root, out, 'run_NSCLC_gamma_scGate_pilot_v1', src / 'private_terminal.json', {'COMPLETE_REPRESENTATIVE_FIXED_SCGATE_PILOT_ALL_MODEL_SLOTS_STRUCTURAL_REVIEW_PENDING'})
    for name, d in t['files'].items(): assert sha(src / name) == d['sha256']
    plan = json.loads((src / 'plan_before_real_fixed_scGate_pilot_RNA.json').read_text())
    worker = src / 'official_fixed_pilot_worker.R'; assert sha(worker) == plan['R_source_sha256']
    assert sha(root / 'gamma_scGate_official_runtime_qualification_20261019_v1/private_terminal.json') == plan['runtime_acceptance_sha256']
    assert sha(root / 'gamma_scGate_native_model_gene_qualification_20261019_v1/private_complete_fixed_two_model_native_qualification.json') == plan['native_model_scope_sha256']
    choice = dict(np.load(src / 'private_outcome_free_uniform_native_selection.npz', allow_pickle=False))['native_row']
    expected = np.sort(np.random.default_rng(2026101901).choice(1254749, 20000, replace=False))
    np.testing.assert_array_equal(choice, expected)
    assert plan['seed'] == 2026101901 and plan['cells'] == 20000 and plan['native_gene_universe'] == 31831
    assert sha(src / 'private_outcome_free_uniform_native_selection.npz') == plan['selection_sha256']
    proof = json.loads((src / 'private_full_input_and_output_transfer_proofs.json').read_text())
    assert proof['raw_matrix_sha256'] == plan['source_matrix_sha256'] and proof['raw_entries'] == 2010550708
    assert proof['all_native_RNA_and_detection_exact'] and proof['all_selected_old_hard_rule_bits_exact'] and proof['anonymous_R_cells']
    assert not proof['clinical_or_author_labels_supplied_to_R'] and proof['no_full_RNA_or_score_matrix_persisted']
    z = dict(np.load(src / 'private_complete_pilot_candidate_classification_sidecar.npz', allow_pickle=False))
    np.testing.assert_array_equal(z['native_row'], choice)
    names = ['HiTME_full', 'PBMC_exact', 'HiTME_native_subset', 'NK_HiTME_exact']
    assert z['models'].tolist() == names and z['classes'].shape == (20000, 4)
    c = pd.DataFrame(z['classes'], columns=names)
    assert c.HiTME_full.eq('UNMEASURABLE').all()
    for name in names[1:]: assert c[name].isin(['Pure', 'Impure']).all()
    raw = c.copy(); raw.insert(0, 'position', np.arange(20000))
    assert hashlib.sha256(raw.to_csv(sep='\t', index=False).encode()).hexdigest() == proof['classification_TSV_sha256']
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    units = np.load(native / 'private_cell_unit_index.npy', allow_pickle=False)
    u = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t', keep_default_na=False).sort_values('unit_index')
    metadata = u.iloc[units[choice]].reset_index(drop=True)
    d = np.load(native / 'private_detected_genes.npy', allow_pickle=False)
    assert proof['selected_all_gene_nonzeros'] == int(d[choice].sum()) == plan['selected_nonzero_bound']
    recorded = pd.read_csv(src / 'private_complete_author_type_model_classification_summary.tsv', sep='\t', keep_default_na=True).set_index(['model', 'major_cell_type', 'sub_cell_type'])
    checked = 0
    for name in names:
        for (major, subtype), group in metadata.groupby(['major_cell_type', 'sub_cell_type'], sort=True):
            row = recorded.loc[(name, major, subtype)]; assert row.source == 'GSE243013' and row.pilot_cells == len(group)
            if name == 'HiTME_full':
                assert pd.isna(row.pure_cells) and row.model_status == 'NATIVE_GENE_UNMEASURABLE'
            else:
                assert row.pure_cells == c.loc[group.index, name].eq('Pure').sum()
                assert row.model_status == 'FIXED_MODEL_EXECUTED_NOT_FORMAL_IDENTITY'
            checked += 1
    assert checked == len(recorded)
    # Full joint counts retain discordance rather than selecting a favorable model.
    joint = metadata[['major_cell_type', 'sub_cell_type']].copy()
    for name in names[1:]: joint[name] = c[name].eq('Pure')
    summary = joint.groupby(['major_cell_type', 'sub_cell_type'] + names[1:], sort=True).size().reset_index(name='pilot_cells')
    summary.to_csv(out / 'private_complete_joint_author_reference_model_flags.tsv', sep='\t', index=False)
    receipt = dict(status='ACCEPTED_ALL_UNIFORM_PILOT_SAMPLING_ANONYMOUS_R_LABEL_TRANSFER_AND_AUTHOR_AGGREGATES_STRUCTURAL_ONLY', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), resource_guard=guard, pilot_rows_checked=20000, model_slots_checked=4, full_R_classification_table_sha256=proof['classification_TSV_sha256'], biological_TCR_or_external_identity_validation=False, independent_graph_or_full_gene_raw_replay=False, clinical_or_formal_identity_adoption=False, scientific_outputs_public=False)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file()) < 1024**2
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
