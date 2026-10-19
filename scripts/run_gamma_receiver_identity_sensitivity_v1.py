"""Post-selection identity sensitivity; native RNA and all readouts stay in RAM."""
import gzip
import hashlib
import io
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_full_native_gene_units_v1 import TrackedReader, chunks
from NSCLC_frozen92_axis_io_v2 import controls, eligibility
from NSCLC_immune_RNA_MCP_cohort_math_v1 import design
from NSCLC_immune_RNA_MCP_cohort_audit_math_v1 import prepare as independent_design
from NSCLC_frozen_axis_block_math_v2 import geometry, evaluate
from NSCLC_frozen_axis_block_audit_math_v1 import geometry as independent_geometry, evaluate as independent_evaluate

MASKS = ['AUTHOR_REFERENCE', 'REFERENCE_WITHOUT_COMPLETE_AB_ANNOTATION',
         'REFERENCE_WITH_OLD_RNA_RULE_SUPPORT', 'RNA_SUPPORTED_REFERENCE_WITHOUT_COMPLETE_AB_ANNOTATION']


def digest(a):
    return hashlib.sha256(np.asarray(a, dtype='<i4').tobytes(order='C')).hexdigest()


def fit_pair(x, y, C, names, blocks, valid):
    x, y, C, blocks = x[valid], y[valid], C[valid], blocks[valid]
    base = dict(author_records=len(x), nonMPR=int((x == 0).sum()), pCR_or_MPR=int((x == 1).sum()), profile_blocks=len(set(blocks)), reviewed_independent_patients=0)
    if min(base['nonMPR'], base['pCR_or_MPR']) < 10:
        return dict(base, status='FEWER_THAN10_PER_RESPONSE')
    m, status = design(x, y, C, names, 'clinical')
    if m is None:
        return dict(base, status=status)
    if m['df'] < 20:
        return dict(base, status='FEWER_THAN20_RESIDUAL_DF')
    try:
        g = geometry(m, blocks)
        f = evaluate(g, y)
    except ValueError as e:
        return dict(base, status=str(e))
    # Independent nuisance construction, explicit deleted-design QR, all blocks.
    other, why = independent_design(x, y, C, names, 'clinical')
    assert other is not None, why
    h = independent_geometry(other, blocks)
    z = independent_evaluate(h, y)
    np.testing.assert_allclose([f['effect'], f['se']], [z['effect'], z['se']], atol=1e-9, rtol=1e-9)
    np.testing.assert_allclose(f['leave_block_effect'], z['leave_block_effect'], atol=1e-9, rtol=1e-9)
    critical = student_t.ppf(.975, g['groups'] - 1)
    return dict(base, status='ESTIMABLE', effect=f['effect'], CV3_se=f['se'], CV3_t_CI_low=f['effect'] - critical * f['se'], CV3_t_CI_high=f['effect'] + critical * f['se'], nuisance_rank=m['nuisance_rank'], residual_df=m['df'], max_block_leverage=g['max_block_leverage'], leave_block_effect_min=float(f['leave_block_effect'].min()), leave_block_effect_max=float(f['leave_block_effect'].max()), all_deleted_QR_effect_and_SE_checks=True)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 2 * 1024**2)
    proof = root / 'gamma_reference_candidate_AB_evidence_join_20261019_v1/private_terminal.json'
    assert sha(proof) == '2d19ac4fa9dfdbb0f246d71e5cd9e8a5c72ae863fbc87ea84208ed821d60749f'
    result_audit = root / 'NSCLC_frozen92_axis_source_result_audit_20261019_v1/private_terminal.json'
    assert sha(result_audit) == '79e8aa29aec12db05d437f6e4ebc24dec7e97772bf62ca594434efe36522b9b7'
    prior = root / 'NSCLC_frozen92_axis_source_consistency_20261019_v1'
    assert json.loads(result_audit.read_text())['producer_terminal_sha256'] == sha(prior / 'private_terminal.json')
    for folder in [proof.parent, prior]:
        for name, info in json.loads((folder / 'private_terminal.json').read_text())['files'].items():
            assert sha(folder / name) == info['sha256']
    inp = root / 'NSCLC_frozen92_axis_source_inputs_20261019_v2'
    ip = json.loads((inp / 'private_terminal.json').read_text())
    for name, info in ip['files'].items():
        assert sha(inp / name) == info['sha256']
    prior_terminal = json.loads((prior / 'private_terminal.json').read_text())
    ia = root / 'NSCLC_frozen92_axis_source_input_audit_20261019_v1/private_terminal.json'
    assert sha(ia) == prior_terminal['input_audit_sha256']
    assert json.loads(ia.read_text())['producer_terminal_sha256'] == sha(inp / 'private_terminal.json')
    previous = pd.read_csv(prior / 'private_all92_source_consistency_results.tsv', sep='\t').set_index('axis_id')
    selected = previous.index[previous.BH92 < .05].tolist()
    assert len(selected) == 1
    a = dict(np.load(inp / 'private_all129_all92_native_axis_and_design_inputs.npz', allow_pickle=False))
    axis = selected[0]
    j = a['axis_ids'].tolist().index(axis)
    s, r = int(a['sender_index'][j]), int(a['receiver_index'][j])
    receiver = str(a['author_states'][r])
    assert receiver in ['T_gdT_TRDV1', 'T_gdT_TRDV2']
    vocabpath = root / 'NSCLC_communication_tensor_family_20261019_v1/private_selected_LR_vocabulary.json'
    vt = json.loads((vocabpath.parent / 'private_terminal.json').read_text())
    assert sha(vocabpath) == vt['files'][vocabpath.name]['sha256']
    entry = json.loads(vocabpath.read_text())[int(a['LR_index'][j])]
    cols = sorted(set(entry['native_ligand_columns'] + entry['native_receptor_columns']))
    lc = [cols.index(c) for c in entry['native_ligand_columns']]
    rc = [cols.index(c) for c in entry['native_receptor_columns']]
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    raw_proof = root / 'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json'
    assert sha(raw_proof) == 'bb20365ee83758a7df8478ceadede2596fe749928dab5765e3bc0866d72dd8c1'
    nt = json.loads((native / 'private_terminal.json').read_text())
    used = ['private_cell_unit_index.npy', 'private_full_RNA.npy', 'private_detected_genes.npy', 'private_all_native_units_and_quality.tsv', 'sum_raw_counts.npy', 'positive_cell_counts.npy', 'mean_logCP10K.npy']
    for name in used:
        assert sha(native / name) == nt['outputs'][name]['sha256']
    unit = np.load(native / used[0], allow_pickle=False)
    RNA = np.load(native / used[1], allow_pickle=False)
    detected = np.load(native / used[2], allow_pickle=False)
    u = pd.read_csv(native / used[3], sep='\t', keep_default_na=False).sort_values('unit_index')
    assert len(unit) == 1254749 and u.unit_index.tolist() == list(range(11701))
    z = dict(np.load(proof.parent / 'private_all_original_union_rows_with_AB_evidence.npz', allow_pickle=False))
    reference = u.sub_cell_type.to_numpy(str)[unit] == receiver
    raw = np.zeros(len(unit), bool)
    ab = raw.copy()
    raw[z['native_row']] = z['raw_rule_bits'] != 0
    ab[z['native_row']] = (z['AB_row_annotation_flags'] & 4) != 0
    masks = [reference, reference & ~ab, reference & raw, reference & raw & ~ab]
    rows = np.flatnonzero(reference).astype(np.int32)
    prefix = rows[:min(1024, len(rows))]
    bits = sum((mask[rows].astype(np.uint8) << k for k, mask in enumerate(masks)), np.zeros(len(rows), np.uint8))
    out = root / 'gamma_receiver_identity_sensitivity_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    download_path = root / 'count_matrix_download_v1/download_state.json'
    download = json.loads(download_path.read_text())
    matrix = Path(download['path'])
    assert download['status'] == 'COMPLETE_DOWNLOAD_GZIP_INTEGRITY_AND_SHA256' and matrix.stat().st_size == download['bytes']
    save('plan_before_new_receiver_expression_and_effects.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), resource_guard=guard, previous_turn_classification='progress_complete_annotation_audit_and_evidence_join', selection='All original92 main-BH candidates; follow-up sensitivity only, old family/P/BH/BY unchanged', selected_axis_ids=selected, masks=MASKS, mask_note='No AB record is not a negative receptor assay; excluding annotated cells is a sensitivity, not validated gamma identity', native_columns=cols, ligand_component_positions=lc, receptor_component_positions=rc, independent_raw_prefix_reference_cells=len(prefix), first_prefix_max_native_row=int(prefix[-1]), raw_audit_sha256=sha(raw_proof), identity_proof_sha256=sha(proof), prior_result_audit_sha256=sha(result_audit), input_terminal_sha256=sha(inp / 'private_terminal.json'), download_receipt_sha256=sha(download_path), native_input_sha256={name: nt['outputs'][name]['sha256'] for name in used}, scale='Original frozen score mean/SD; no new within-mask rescaling', denominator='Original complete TNK major and global immune background unchanged; recompute receiver cell count, within-major fraction and per-cell RNA/detected quality', comparisons='Four own-eligible masks plus original baseline on same records as each of three restricted masks; seven fixed descriptions, no P or FDR', model='Original full clinical/RNA/count/fraction nuisance span; >=50 records and blocks, >=10 each response, >=20 residual df; reject nonidentifiable deletions. CV3 t pointwise intervals and every explicit deleted QR comparison in RAM', no_new_download=True, full_expression_or_score_arrays_persisted=False, output_cap_bytes=1024**2, formal_identity_or_state_write=False))
    buf = io.BytesIO()
    np.savez_compressed(buf, native_row=rows, identity_sensitivity_mask_bits=bits)
    assert len(buf.getvalue()) < 256 * 1024
    (out / 'private_receiver_metadata_mask_sidecar.npz').write_bytes(buf.getvalue())
    counts = np.array([np.bincount(unit[mask], minlength=len(u)) for mask in masks])
    qualities = np.array([[np.divide(np.bincount(unit[mask], weights=np.log1p(v[mask]), minlength=len(u)), counts[k], out=np.full(len(u), np.nan), where=counts[k] > 0) for v in [RNA, detected]] for k, mask in enumerate(masks)])
    den = u.groupby(['sampleID', 'major_cell_type']).n_cells.sum()
    majors = u.loc[u.sub_cell_type == receiver, 'major_cell_type'].unique().tolist()
    assert len(majors) == 1
    denom = np.array([den.get((case, majors[0]), 0) for case in a['author_record_ids']])
    uid = a['original_unit_indices'][:, r]
    safe = np.maximum(uid, 0)
    C0, names = controls(a, j)
    original_valid, _, _ = eligibility(a, j)
    record_counts = counts[:, safe].copy()
    record_counts[:, uid < 0] = 0
    np.testing.assert_array_equal(record_counts[0], a['unit_cells'][:, r])
    coverage = []
    for k, name in enumerate(MASKS):
        valid = original_valid & (record_counts[k] >= 20)
        coverage.append(dict(source='GSE243013', mask=name, total_native_receiver_cells=int(masks[k].sum()), original_129_records=129, receiver_atleast20=int((record_counts[k] >= 20).sum()), old_complete_design_and_receiver_atleast20=int(valid.sum()), reviewed_independent_patients=0, clinical_effects_computed=False))
    pd.DataFrame(coverage).to_csv(out / 'private_complete4_effect_free_receiver_coverage.tsv', sep='\t', index=False)
    # Existing native matrix only. All cells provide source integrity/parity checks.
    x = np.zeros((len(unit), len(cols)), np.int32)
    cmap = np.full(31831, -1, np.int16)
    cmap[cols] = np.arange(len(cols))
    totals = np.zeros(len(unit), np.int64)
    ngenes = np.zeros(len(unit), np.int32)
    entries, last, checkpoint = 0, -1, -1
    start = time.monotonic()
    with matrix.open('rb') as binary:
        tracked = TrackedReader(binary)
        with gzip.GzipFile(fileobj=tracked, mode='rb') as stream:
            assert stream.readline().startswith(b'%%MatrixMarket matrix coordinate real general')
            line = stream.readline()
            while line.startswith(b'%'):
                line = stream.readline()
            dims = tuple(map(int, line.split()))
            assert dims == (len(unit), 31831, 2010550708)
            for block in chunks(stream):
                assert np.isfinite(block).all() and np.array_equal(block, np.floor(block)) and (block[:, 2] > 0).all() and block[:, 2].max() < 2**31
                row, col, value = block[:, 0].astype(np.int64) - 1, block[:, 1].astype(np.int64) - 1, block[:, 2].astype(np.int32)
                coord = row * 31831 + col
                assert row.min() >= 0 and row.max() < len(unit) and col.min() >= 0 and col.max() < 31831 and coord[0] > last and (np.diff(coord) > 0).all()
                last = int(coord[-1])
                lo, hi = int(row[0]), int(row[-1]) + 1
                totals[lo:hi] += np.bincount(row - lo, weights=value, minlength=hi - lo).astype(np.int64)
                ngenes[lo:hi] += np.bincount(row - lo, minlength=hi - lo).astype(np.int32)
                keep = cmap[col] >= 0
                x[row[keep], cmap[col[keep]]] = value[keep]
                entries += len(block)
                if entries // 50000000 != checkpoint:
                    checkpoint = entries // 50000000
                    save('private_progress.json', dict(stage='STREAMING_EXISTING_NATIVE_RNA', entries=entries, expected_entries=dims[2], elapsed_seconds=time.monotonic() - start, observed_UTC=datetime.now(timezone.utc).isoformat()))
    assert entries == dims[2] and tracked.bytes == download['bytes'] and tracked.digest.hexdigest() == download['sha256']
    np.testing.assert_array_equal(totals, RNA)
    np.testing.assert_array_equal(ngenes, detected)
    logx = np.log1p(x.astype(float) * 10000 / RNA[:, None])
    mean = np.load(native / 'mean_logCP10K.npy', mmap_mode='r')
    sums = np.load(native / 'sum_raw_counts.npy', mmap_mode='r')
    positive = np.load(native / 'positive_cell_counts.npy', mmap_mode='r')
    for k, col in enumerate(cols):
        np.testing.assert_array_equal(np.bincount(unit, weights=x[:, k], minlength=len(u)).astype(np.int64), sums[:, col])
        np.testing.assert_array_equal(np.bincount(unit[x[:, k] > 0], minlength=len(u)), positive[:, col])
        values = np.bincount(unit, weights=logx[:, k], minlength=len(u)) / u.n_cells.to_numpy()
        np.testing.assert_allclose(values, mean[:, col], atol=1e-10, rtol=1e-9)
    ligand = np.zeros(len(uid))
    sender_uid = a['original_unit_indices'][:, s]
    for i in np.flatnonzero(sender_uid >= 0):
        ligand[i] = np.float32(min(float(mean[sender_uid[i], c]) for c in entry['native_ligand_columns']))
    mu, sd = a['blind_score_mean_SD'][j]
    output, exclusions = [], []
    for k, name in enumerate(MASKS):
        receiver_mean = np.array([np.divide(np.bincount(unit[mask], weights=logx[mask, c], minlength=len(u)), counts[k], out=np.full(len(u), np.nan), where=counts[k] > 0) for c in rc for mask in [masks[k]]]).min(axis=0)
        rec = np.where(uid >= 0, receiver_mean[safe], np.nan).astype(np.float32)
        score = np.log1p(ligand.astype(np.float32).astype(float) * rec.astype(float)).astype(np.float32).astype(float)
        C = C0.copy()
        C[:, names.index('logN_receiver')] = np.log1p(record_counts[k])
        C[:, names.index('receiver_meanlogRNA')] = np.where(uid >= 0, qualities[k, 0, safe], np.nan)
        C[:, names.index('receiver_meanloggenes')] = np.where(uid >= 0, qualities[k, 1, safe], np.nan)
        C[:, names.index('receiver_within_major_fraction')] = np.divide(record_counts[k], denom, out=np.full(len(uid), np.nan), where=denom > 0)
        y = (score - mu) / sd
        valid = original_valid & (record_counts[k] >= 20) & np.isfinite(np.c_[y, C]).all(axis=1)
        if k == 0:
            np.testing.assert_allclose(C[original_valid], C0[original_valid], atol=1e-10, rtol=1e-10)
            np.testing.assert_array_equal(score[original_valid], a['raw_axis_score'][original_valid, j])
        f = fit_pair(a['response'], y, C, names, a['profile_blocks'], valid)
        output.append(dict(source='GSE243013', axis_id=axis, mask=name, comparison='REDEFINED_RECEIVER_OWN_ELIGIBLE_RECORDS', **f))
        if k == 0:
            assert f['status'] == 'ESTIMABLE'
            np.testing.assert_allclose([f['effect'], f['CV3_se']], previous.loc[axis, ['effect', 'CV3_se']].to_numpy(float), atol=1e-9, rtol=1e-9)
        else:
            baseline = fit_pair(a['response'], a['standardized_axis_score'][:, j], C0, names, a['profile_blocks'], valid)
            output.append(dict(source='GSE243013', axis_id=axis, mask=name, comparison='ORIGINAL_RECEIVER_SAME_ELIGIBLE_RECORDS', **baseline))
        for i, case in enumerate(a['author_record_ids']):
            reason = 'INCLUDED' if valid[i] else ('ORIGINAL_MAIN_RECORD_INELIGIBLE' if not original_valid[i] else 'REDEFINED_RECEIVER_FEWER_THAN20_CELLS' if record_counts[k, i] < 20 else 'NONFINITE_REDEFINED_READOUT_OR_CONTROLS')
            exclusions.append(dict(mask=name, author_record_id=case, included=bool(valid[i]), reason=reason))
    assert len(output) == 7 and len(exclusions) == 516
    pd.DataFrame(output).to_csv(out / 'private_complete7_identity_sensitivity_effects.tsv', sep='\t', index=False)
    pd.DataFrame(exclusions).to_csv(out / 'private_complete516_record_exclusions.tsv.gz', sep='\t', index=False, compression='gzip')
    save('private_full_native_integrity_and_prefix_digests.json', dict(raw_matrix_sha256=tracked.digest.hexdigest(), raw_entries=entries, full_native_RNA_and_detected_exact=True, all_unit_selected_gene_sum_detection_and_mean_parity=True, raw_prefix_reference_cells=len(prefix), raw_prefix_native_row_sha256=digest(prefix), raw_prefix_gene_counts_sha256=digest(x[prefix]), all_estimable_models_independent_nuisance_and_deleted_QR_checked=True, original_main_effect_and_CV3_SE_reconstructed=True, score_arrays_persisted=False))
    receipt = dict(status='COMPLETE_FIXED7_IDENTITY_SENSITIVITY_DESCRIPTIONS_FULL_NATIVE_PARITY_AND_ALL_MODEL_QR_RAW_PREFIX_AUDIT_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), resource_guard=guard, identity_mask_positions=4, source_followup_positions=7, record_exclusion_positions=516, new_P_or_FDR=False, original_inference_changed=False, final_identity_or_state_changed=False, scientific_outputs_public=False, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    assert sum(v['bytes'] for v in receipt['files'].values()) < 1024**2
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
