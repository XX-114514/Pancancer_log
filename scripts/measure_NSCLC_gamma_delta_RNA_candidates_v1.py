"""First gamma-delta identity calibration on all existing native NSCLC RNA.

Keep original author identities frozen. All cells count arrays remain in RAM;
only candidate indices/rule bits and private aggregate descriptions persist.
"""
import gzip
import io
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from gamma_delta_RNA_rules_v1 import MARKERS, RULES, flags, digest_array, fixture
from prepare_NSCLC_full_native_gene_units_v1 import TrackedReader, chunks
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_native_MFA_interpretation_inputs_v1 import sha


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); guard = resource(root, 3 * 1024**2)
    proof = root / 'NSCLC_full_native_gene_unit_audit_20261016_v1/private_terminal.json'
    assert sha(proof) == 'bb20365ee83758a7df8478ceadede2596fe749928dab5765e3bc0866d72dd8c1'
    metadata_proof = root / 'NSCLC_original_author_lineage_audit_20261019_v1/private_terminal.json'
    assert sha(metadata_proof) == 'ac887dd08ca0710f763252624ff622137c7b918dc9130f651782dbdb888fd03d'
    out = root / 'NSCLC_gamma_delta_RNA_candidates_20261019_v1'; assert not out.exists(); out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    source = root / 'NSCLC_full_native_gene_units_20261016_v1'
    axis = pd.read_csv(root / 'GSE243013_genes.csv.gz', keep_default_na=False).geneSymbol.str.strip().str.upper().to_numpy(str)
    assert len(axis) == 31831 and len(set(axis)) == len(axis)
    columns = pd.Index(axis).get_indexer(MARKERS); assert (columns >= 0).all()
    unit = np.load(source / 'private_cell_unit_index.npy', allow_pickle=False)
    RNA = np.load(source / 'private_full_RNA.npy', allow_pickle=False)
    detected = np.load(source / 'private_detected_genes.npy', allow_pickle=False)
    u = pd.read_csv(source / 'private_all_native_units_and_quality.tsv', sep='\t', keep_default_na=False).sort_values('unit_index')
    assert len(unit) == len(RNA) == len(detected) == 1254749 and u.unit_index.tolist() == list(range(11701))
    download_path = root / 'count_matrix_download_v1/download_state.json'; download = json.loads(download_path.read_text())
    assert download['status'] == 'COMPLETE_DOWNLOAD_GZIP_INTEGRITY_AND_SHA256'
    matrix = Path(download['path']); assert matrix.stat().st_size == download['bytes']
    save('plan_before_any_gamma_delta_RNA_counts.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), rule_module_sha256=sha(Path(__file__).with_name('gamma_delta_RNA_rules_v1.py')), resource_guard=guard, raw_unit_audit_sha256=sha(proof), author_metadata_audit_sha256=sha(metadata_proof), download_receipt_sha256=sha(download_path), native_gene_axis_sha256=sha(root / 'GSE243013_genes.csv.gz'), markers=MARKERS, rules=RULES, independent_raw_prefix_cells=8192, native_cells=len(unit), native_units=len(u), scope='All original native cells, all author types and243 records. No outcome, current LR endpoint or functional score used. Candidate sidecar only. Author tags are reference, not flow/TCR truth. Counts>=1 and>=2 both reported, no threshold retuning.', original_authority_mutated=False, first_phase_output_cap_bytes=2 * 1024**2))
    save('private_rule_fixture.json', fixture())
    x = np.zeros((len(unit), len(MARKERS)), np.int32)
    totals = np.zeros(len(unit), np.int64); genes = np.zeros(len(unit), np.int32)
    cmap = np.full(len(axis), -1, np.int16); cmap[columns] = np.arange(len(columns))
    entries = 0; last = -1; start = time.monotonic(); checkpoint = -1
    with matrix.open('rb') as binary:
        tracked = TrackedReader(binary)
        with gzip.GzipFile(fileobj=tracked, mode='rb') as stream:
            assert stream.readline().startswith(b'%%MatrixMarket matrix coordinate real general')
            line = stream.readline()
            while line.startswith(b'%'): line = stream.readline()
            dims = tuple(map(int, line.split())); assert dims == (len(unit), len(axis), 2010550708)
            for a in chunks(stream):
                assert np.isfinite(a).all() and np.array_equal(a, np.floor(a)) and (a[:, 2] > 0).all() and a[:, 2].max() < 2**31
                row = a[:, 0].astype(np.int64) - 1; col = a[:, 1].astype(np.int64) - 1; value = a[:, 2].astype(np.int32)
                assert row.min() >= 0 and row.max() < len(unit) and col.min() >= 0 and col.max() < len(axis)
                coord = row * len(axis) + col; assert coord[0] > last and (np.diff(coord) > 0).all(); last = int(coord[-1])
                lo = int(row[0]); hi = int(row[-1]) + 1
                totals[lo:hi] += np.bincount(row - lo, weights=value, minlength=hi - lo).astype(np.int64)
                genes[lo:hi] += np.bincount(row - lo, minlength=hi - lo).astype(np.int32)
                keep = cmap[col] >= 0; x[row[keep], cmap[col[keep]]] = value[keep]
                entries += len(a)
                if entries // 50000000 != checkpoint:
                    checkpoint = entries // 50000000
                    save('private_progress.json', dict(stage='READING_EXISTING_NATIVE_RNA_ONLY', entries=entries, expected_entries=dims[2], elapsed_seconds=time.monotonic() - start, observed_UTC=datetime.now(timezone.utc).isoformat()))
    assert entries == dims[2] and tracked.bytes == download['bytes'] and tracked.digest.hexdigest() == download['sha256']
    np.testing.assert_array_equal(totals, RNA); np.testing.assert_array_equal(genes, detected)
    expected_sums = np.load(source / 'sum_raw_counts.npy', mmap_mode='r')[:, columns]
    expected_positive = np.load(source / 'positive_cell_counts.npy', mmap_mode='r')[:, columns]
    for k in range(len(MARKERS)):
        actual = np.bincount(unit, weights=x[:, k], minlength=len(u)).astype(np.int64)
        positive = np.bincount(unit[x[:, k] > 0], minlength=len(u))
        np.testing.assert_array_equal(actual, expected_sums[:, k]); np.testing.assert_array_equal(positive, expected_positive[:, k])
    bits = flags(x); candidates = np.flatnonzero(bits).astype(np.int32)
    payload = io.BytesIO()
    np.savez_compressed(payload, native_row=candidates, unit_index=unit[candidates], rule_bits=bits[candidates], TRDV1_detected=(x[candidates, 7] > 0), TRDV2_detected=(x[candidates, 8] > 0))
    blob = payload.getvalue(); assert len(blob) < 1024**2
    (out / 'private_candidate_identity_sidecar.npz').write_bytes(blob)
    records = []
    for j, rule in enumerate(RULES):
        found = (bits & (1 << j)) > 0
        count = np.bincount(unit[found], minlength=len(u))
        v1 = np.bincount(unit[found & (x[:, 7] > 0)], minlength=len(u))
        v2 = np.bincount(unit[found & (x[:, 8] > 0)], minlength=len(u))
        for i, r in enumerate(u.itertuples(index=False)):
            records.append(dict(unit_index=i, sampleID=r.sampleID, major_cell_type=r.major_cell_type, sub_cell_type=r.sub_cell_type, rule=rule, native_cells=int(r.n_cells), candidate_cells=int(count[i]), TRDV1_detected_candidate_cells=int(v1[i]), TRDV2_detected_candidate_cells=int(v2[i])))
    table = pd.DataFrame(records); table.to_csv(out / 'private_all_author_unit_rule_counts.tsv.gz', sep='\t', index=False, compression='gzip')
    summary = table.groupby(['major_cell_type', 'sub_cell_type', 'rule'], sort=True)[['native_cells', 'candidate_cells', 'TRDV1_detected_candidate_cells', 'TRDV2_detected_candidate_cells']].sum().reset_index()
    summary.to_csv(out / 'private_complete_author_type_rule_summary.tsv', sep='\t', index=False)
    source_rows = []
    for rule, g in table.groupby('rule', sort=True):
        bysample = g.groupby('sampleID').candidate_cells.sum()
        source_rows.append(dict(source='GSE243013', rule=rule, author_records=len(bysample), records_candidate_atleast20=int(bysample.ge(20).sum()), candidate_cells=int(g.candidate_cells.sum()), all_native_cells=len(unit), reviewed_independent_patients=0, no_clinical_association=True))
    pd.DataFrame(source_rows).to_csv(out / 'private_source_candidate_coverage.tsv', sep='\t', index=False)
    save('private_prefix8192_raw_and_rule_digests.json', dict(cells=8192, marker_columns=columns.tolist(), marker_counts_sha256=digest_array(x[:8192], '<i4'), RNA_sha256=digest_array(RNA[:8192], '<i8'), detected_sha256=digest_array(detected[:8192], '<i4'), rule_bits_sha256=digest_array(bits[:8192], 'u1')))
    save('private_all_cell_rule_digest.json', dict(rule_bits_sha256=digest_array(bits, 'u1'), core2_subset_core1=bool(np.all(((bits & 4) == 0) | ((bits & 1) > 0))), supported2_subset_supported1=bool(np.all(((bits & 8) == 0) | ((bits & 2) > 0))), TRDC_positive_cells=int((x[:, 4] > 0).sum()), TRDC_positive_noncandidate_cells=int(((x[:, 4] > 0) & (bits == 0)).sum()), candidate_membership_sha256=sha(out / 'private_candidate_identity_sidecar.npz')))
    files = {p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()}
    assert sum(v['bytes'] for v in files.values()) < 2 * 1024**2
    receipt = dict(status='COMPLETE_ALL_NATIVE_NSCLC_GAMMA_DELTA_FIXED_RULE_CANDIDATES_RAW_PREFIX_AUDIT_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), files=files, raw_matrix_sha256=tracked.digest.hexdigest(), raw_nonzero_coordinates=entries, native_cells=len(unit), all_native_RNA_and_detection_exact=True, full_unit_marker_sum_and_detection_positions_checked=len(u) * len(MARKERS), full_marker_arrays_persisted=False, clinical_or_formal_identity_adoption=False)
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
