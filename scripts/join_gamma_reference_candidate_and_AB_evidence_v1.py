"""Join accepted metadata evidence without adopting new cell identities."""
import io
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 512 * 1024)
    prior = root / 'gamma_candidate_evidence_nomenclature_reconciliation_20261019_v1'
    old_audit = root / 'gamma_candidate_evidence_reconciliation_audit_20261019_v1/private_terminal.json'
    assert sha(old_audit) == 'a5dabd6d1003f0b4576ebe2c9dcdeaa003b875a10fd55bea344a60852d3ab706'
    assert json.loads(old_audit.read_text())['producer_terminal_sha256'] == sha(prior / 'private_terminal.json')
    src = root / 'gamma_alpha_beta_annotation_diagnostic_20261019_v2'
    audit = root / 'gamma_alpha_beta_annotation_diagnostic_audit_20261019_v2/private_terminal.json'
    assert sha(src / 'private_terminal.json') == 'fd1bc131cc4fda0d2ddf55faaac21f67b4b69df24ba493c20732493d6d2644bc'
    accepted = json.loads(audit.read_text())
    assert accepted['status'] == 'ACCEPTED_ALL_NONSEQUENCE_METADATA_EXACT_KEYS_FLAGS_51_GROUPS_3264_PATTERNS_20_CHANNELS_AND_SOURCE_PREFIX8192'
    assert accepted['producer_terminal_sha256'] == sha(src / 'private_terminal.json')
    for folder in [prior, src]:
        for name, info in json.loads((folder / 'private_terminal.json').read_text())['files'].items():
            assert sha(folder / name) == info['sha256']
    out = root / 'gamma_reference_candidate_AB_evidence_join_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_joint_evidence_counts.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), original_union_sha256=sha(prior / 'private_union_reference_and_candidate_identity_evidence.npz'), accepted_annotation_audit_sha256=sha(audit), all_original_union_rows_and_columns_preserved=True, AB_status_codes={'0': 'NO_SOURCE_RECORD_NOT_NEGATIVE', '1': 'SOURCE_RECORD_WITHOUT_SAME_ROW_COMPLETE_ALPHA_BETA', '2': 'SOURCE_RECORD_WITH_SAME_ROW_COMPLETE_ALPHA_BETA'}, no_new_RNA_or_clinical_tests=True, no_final_identity_adoption=True, resource_guard=guard))
    original = dict(np.load(prior / 'private_union_reference_and_candidate_identity_evidence.npz', allow_pickle=False))
    annotation = dict(np.load(src / 'private_native_annotation_evidence_sidecar.npz', allow_pickle=False))
    rows = original['native_row']
    assert len(rows) > 0 and np.all(np.diff(rows) > 0)
    assert np.all(np.diff(annotation['native_row']) > 0)
    pos = np.searchsorted(annotation['native_row'], rows)
    safe = np.minimum(pos, len(annotation['native_row']) - 1)
    present = (pos < len(annotation['native_row'])) & (annotation['native_row'][safe] == rows)
    joined = dict(original)
    for field in ['valid_component_union', 'unknown_component_union', 'row_annotation_flags', 'source_rows', 'distinct_nonsequence_rows']:
        value = np.zeros(len(rows), annotation[field].dtype)
        value[present] = annotation[field][pos[present]]
        joined['AB_' + field] = value
    status = np.zeros(len(rows), np.uint8)
    status[present] = 1
    status[present & ((joined['AB_row_annotation_flags'] & 4) != 0)] = 2
    joined['AB_annotation_status'] = status
    blob = io.BytesIO()
    np.savez_compressed(blob, **joined)
    assert len(blob.getvalue()) <= 256 * 1024
    sidecar = out / 'private_all_original_union_rows_with_AB_evidence.npz'
    sidecar.write_bytes(blob.getvalue())
    got = dict(np.load(sidecar, allow_pickle=False))
    for name, value in original.items():
        np.testing.assert_array_equal(got[name], value)
    # A separately built dictionary checks every selected row and every new field.
    lookup = {int(row): j for j, row in enumerate(annotation['native_row'])}
    for i, row in enumerate(rows):
        j = lookup.get(int(row))
        expected_status = 0 if j is None else (2 if int(annotation['row_annotation_flags'][j]) & 4 else 1)
        assert int(got['AB_annotation_status'][i]) == expected_status
        for field in ['valid_component_union', 'unknown_component_union', 'row_annotation_flags', 'source_rows', 'distinct_nonsequence_rows']:
            assert int(got['AB_' + field][i]) == (0 if j is None else int(annotation[field][j]))
    summary = []
    for reference in range(3):
        for raw_rule in [False, True]:
            for ab_code in range(3):
                mask = (original['author_reference_code'] == reference) & ((original['raw_rule_bits'] != 0) == raw_rule) & (status == ab_code)
                summary.append(dict(source='GSE243013', author_reference_code=reference, any_old_raw_rule=raw_rule, AB_annotation_status=ab_code, cells=int(mask.sum()), scope='ORIGINAL_EVIDENCE_UNION_NOT_ALL_NONREFERENCE_CELLS'))
    table = pd.DataFrame(summary)
    assert len(table) == 18 and table.cells.sum() == len(rows)
    table.to_csv(out / 'private_complete18_joint_evidence_descriptions.tsv', sep='\t', index=False)
    original_coverage = pd.read_csv(src / 'private_all51_RNA_author_annotation_coverage.tsv', sep='\t').set_index('RNA_author_type')
    for code, label in [(1, 'T_gdT_TRDV1'), (2, 'T_gdT_TRDV2')]:
        mask = original['author_reference_code'] == code
        assert mask.sum() == original_coverage.loc[label, 'native_cells']
        assert np.count_nonzero(mask & present) == original_coverage.loc[label, 'TCR_record_present']
        assert np.count_nonzero(mask & (status == 2)) == original_coverage.loc[label, 'both_complete_in_same_row']
    receipt = dict(status='COMPLETE_ALL_ORIGINAL_UNION_METADATA_JOIN_EVERY_ROW_VERIFIED_NO_IDENTITY_ADOPTION', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), annotation_audit_sha256=sha(audit), original_union_rows_preserved=True, original_columns_exact=True, new_metadata_fields_dictionary_checked=True, descriptive_positions=18, no_receptor_productivity_or_gamma_truth=True, no_final_identity_or_state_write=True, no_new_clinical_or_communication_test=True, public_scientific_results=False, resource_guard=guard, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file()) < 512 * 1024
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
