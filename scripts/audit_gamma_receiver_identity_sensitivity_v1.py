"""Independent all-mask metadata reconstruction and bounded receptor-row raw audit."""
import gzip
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


def digest(a):
    return hashlib.sha256(np.asarray(a, dtype='<i4').tobytes(order='C')).hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 512 * 1024)
    out = root / 'gamma_receiver_identity_sensitivity_audit_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    src = root / 'gamma_receiver_identity_sensitivity_20261019_v1'
    t = wait(root, out, 'run_gamma_receiver_identity_sensitivity_v1', src / 'private_terminal.json', {'COMPLETE_FIXED7_IDENTITY_SENSITIVITY_DESCRIPTIONS_FULL_NATIVE_PARITY_AND_ALL_MODEL_QR_RAW_PREFIX_AUDIT_PENDING'})
    for name, info in t['files'].items():
        assert sha(src / name) == info['sha256']
    plan = json.loads((src / 'plan_before_new_receiver_expression_and_effects.json').read_text())
    proof = json.loads((src / 'private_full_native_integrity_and_prefix_digests.json').read_text())
    a = dict(np.load(root / 'NSCLC_frozen92_axis_source_inputs_20261019_v2/private_all129_all92_native_axis_and_design_inputs.npz', allow_pickle=False))
    j = a['axis_ids'].tolist().index(plan['selected_axis_ids'][0])
    receiver = str(a['author_states'][int(a['receiver_index'][j])])
    labels = pd.read_csv(root / 'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz', usecols=['sub_cell_type'], dtype=str, keep_default_na=False).sub_cell_type.to_numpy(str)
    rows = np.where(labels == receiver)[0].astype(np.int32)
    original = dict(np.load(root / 'gamma_reference_candidate_AB_evidence_join_20261019_v1/private_all_original_union_rows_with_AB_evidence.npz', allow_pickle=False))
    lookup = {int(row): (int(raw), int(flag)) for row, raw, flag in zip(original['native_row'], original['raw_rule_bits'], original['AB_row_annotation_flags'])}
    expected = []
    for row in rows:
        raw, ab = lookup[int(row)]
        support = raw != 0
        no_complete_ab = (ab & 4) == 0
        expected.append(1 + 2 * no_complete_ab + 4 * support + 8 * (support and no_complete_ab))
    expected = np.asarray(expected, np.uint8)
    got = dict(np.load(src / 'private_receiver_metadata_mask_sidecar.npz', allow_pickle=False))
    np.testing.assert_array_equal(got['native_row'], rows)
    np.testing.assert_array_equal(got['identity_sensitivity_mask_bits'], expected)
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    unit = np.load(native / 'private_cell_unit_index.npy', allow_pickle=False)
    u = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t', keep_default_na=False)
    coverage = pd.read_csv(src / 'private_complete4_effect_free_receiver_coverage.tsv', sep='\t').set_index('mask')
    target = a['original_unit_indices'][:, int(a['receiver_index'][j])]
    for k, name in enumerate(plan['masks']):
        selected = rows[(expected & (1 << k)) != 0]
        counts = pd.Series(unit[selected]).value_counts().reindex(range(len(u)), fill_value=0).to_numpy(int)
        counts = np.where(target >= 0, counts[np.maximum(target, 0)], 0)
        assert coverage.loc[name, 'total_native_receiver_cells'] == len(selected)
        assert coverage.loc[name, 'receiver_atleast20'] == (counts >= 20).sum()
    effects = pd.read_csv(src / 'private_complete7_identity_sensitivity_effects.tsv', sep='\t')
    exclusions = pd.read_csv(src / 'private_complete516_record_exclusions.tsv.gz', sep='\t')
    assert len(effects) == 7 and len(exclusions) == 516
    assert not any('P' == c or 'BH' in c or 'FDR' in c for c in effects.columns)
    estimable = effects.status == 'ESTIMABLE'
    assert effects.loc[estimable, 'all_deleted_QR_effect_and_SE_checks'].eq(True).all()
    assert proof['all_estimable_models_independent_nuisance_and_deleted_QR_checked'] and proof['original_main_effect_and_CV3_SE_reconstructed']
    prefix = rows[:plan['independent_raw_prefix_reference_cells']]
    assert digest(prefix) == proof['raw_prefix_native_row_sha256']
    assert len(prefix) == 1024 and int(prefix[-1]) == plan['first_prefix_max_native_row']
    x = np.zeros((len(prefix), len(plan['native_columns'])), np.int32)
    column = {int(c): j for j, c in enumerate(plan['native_columns'])}
    row_lookup = pd.Index(prefix)
    download = json.loads((root / 'count_matrix_download_v1/download_state.json').read_text())
    assert download['sha256'] == proof['raw_matrix_sha256']
    reached = False
    examined = 0
    with gzip.open(download['path'], 'rb') as stream:
        assert stream.readline().startswith(b'%%MatrixMarket matrix coordinate real general')
        line = stream.readline()
        while line.startswith(b'%'):
            line = stream.readline()
        assert tuple(map(int, line.split())) == (1254749, 31831, 2010550708)
        reader = pd.read_csv(stream, sep=r'\s+', header=None, chunksize=250000, dtype=np.float64)
        for frame in reader:
            values = frame.to_numpy()
            assert np.array_equal(values, np.floor(values)) and np.isfinite(values).all() and (values[:, 2] > 0).all()
            rr, cc, vv = values[:, 0].astype(np.int64) - 1, values[:, 1].astype(np.int64) - 1, values[:, 2].astype(np.int64)
            ii = row_lookup.get_indexer(rr)
            for c, dest in column.items():
                take = (ii >= 0) & (cc == c)
                x[ii[take], dest] = vv[take]
            examined += len(frame)
            if rr[-1] > prefix[-1]:
                reached = True
                break
    assert reached and digest(x) == proof['raw_prefix_gene_counts_sha256']
    scope = dict(full_metadata_masks_independently_rebuilt=True, complete4_coverage_positions_checked=True, independent_raw_reference_cells=len(prefix), raw_prefix_coordinates_examined=examined, raw_selected_gene_positions=len(prefix) * len(column), all_same_process_models_checked_by_separate_nuisance_and_explicit_deleted_QR_implementation=True, original_complete_raw_file_independently_replayed=False, full_new_expression_or_score_arrays_persisted=False, identity_or_clinical_validation=False)
    p = out / 'private_complete_audit_scope.json'
    p.write_text(json.dumps(scope, indent=2) + '\n')
    receipt = dict(status='ACCEPTED_ALL_IDENTITY_METADATA_MASKS_AND_FIXED1024_REFERENCE_RAW_PREFIX_MODEL_QR_CHECKED_WITHIN_PRODUCER', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), audit_scope_sha256=sha(p), resource_guard=guard, new_P_or_FDR=False, formal_identity_adoption=False, scientific_outputs_public=False)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
