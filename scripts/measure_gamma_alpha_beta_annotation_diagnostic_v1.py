"""Exact RNA-key bridge from processed annotation; sequence columns excluded."""
import gzip
import hashlib
import io
import json
import os
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from gamma_alpha_beta_annotation_rules_v1 import KEYS, GENES, SELECTED, classify, fixture
from gamma_delta_RNA_rules_v1 import RULES
from qualify_public_GSE243013_TCR_header_v1 import URL
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 11 * 1024**2)
    available = int(next(v.split()[1] for v in Path('/proc/meminfo').read_text().splitlines() if v.startswith('MemAvailable:'))) * 1024
    assert available > 4 * 1024**3
    header_root = root / 'public_GSE243013_TCR_header_qualification_20261019_v1'
    assert sha(header_root / 'private_terminal.json') == 'd8dab0d2258cc06419e38be3c51244a149274fec41037a8f4a3f0fc13917069f'
    previous = root / 'gamma_group_lineage_scientific_review_20261019_v1/private_terminal.json'
    assert sha(previous) == 'a9751cec1cf76fcd77e951f9e01789f7a34d65ff6f30a9b3b30464cf7a5d166a'
    for producer, audit_name, expected, filename in [
        ('NSCLC_gamma_delta_RNA_candidates_20261019_v1', 'NSCLC_gamma_delta_RNA_candidate_audit_20261019_v1', 'e4d6c57c1727d33dd4d123178929125787706afa64964520c7c0ab1ff7429502', 'private_candidate_identity_sidecar.npz'),
        ('NSCLC_gamma_scGate_fixed_pilot_20261019_v1', 'NSCLC_gamma_scGate_fixed_pilot_audit_20261019_v1', '519a8671d37cc311d45c799acbd6ce70a2c76fafe956a70f209ee7fa62175cbf', 'private_complete_pilot_candidate_classification_sidecar.npz'),
    ]:
        audit_path = root / audit_name / 'private_terminal.json'
        parent_path = root / producer / 'private_terminal.json'
        assert sha(audit_path) == expected
        assert sha(parent_path) == json.loads(audit_path.read_text())['producer_terminal_sha256']
        assert sha(root / producer / filename) == json.loads(parent_path.read_text())['files'][filename]['sha256']
    out = root / 'gamma_alpha_beta_annotation_diagnostic_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_processed_annotation_rows.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), rule_module_sha256=sha(Path(__file__).with_name('gamma_alpha_beta_annotation_rules_v1.py')), resource_guard=guard, MemAvailable_bytes=available, public_URL=URL, expected_compressed_bytes=15940803, compressed_cap_bytes=32 * 1024**2, selected_nonsequence_fields=SELECTED, exact_key=KEYS, full_source_download_to_disk=False, processed_annotation_stream_in_RAM=True, sequence_fields_used_or_persisted=False, minimum_source_prefix_records_for_independent_check=8192, scientific_scope='All native author types and complete old candidate channels; annotation syntax, no productivity, clinical P or final labels. Missing alpha/beta is not gamma.'))
    save('private_rule_fixture.json', fixture())
    with urllib.request.urlopen(urllib.request.Request(URL, headers={'Accept-Encoding': 'identity', 'Cache-Control': 'no-cache'}), timeout=30) as response:
        assert response.status == 200
        chunks = []
        total = 0
        for chunk in iter(lambda: response.read(1024**2), b''):
            total += len(chunk)
            assert total <= 32 * 1024**2
            chunks.append(chunk)
        blob = b''.join(chunks)
        assert total == len(blob) == 15940803
        declared = response.headers.get('Content-Length')
        assert declared is None or int(declared) == len(blob)
        source_headers = dict(ETag=response.headers.get('ETag'), Last_Modified=response.headers.get('Last-Modified'))
    hp = json.loads((header_root / 'private_header_qualification.json').read_text())
    assert hashlib.sha256(blob[:8192]).hexdigest() == hp['partial_compressed_sha256']
    with gzip.GzipFile(fileobj=io.BytesIO(blob), mode='rb') as raw:
        frame = pd.read_csv(raw, usecols=SELECTED, dtype=str, keep_default_na=False, engine='c')[SELECTED]
        assert raw.read() == b''  # Full stream EOF and gzip CRC validation.
    assert len(frame) > 0 and list(frame) == SELECTED
    canonical = frame.to_csv(index=False, lineterminator='\n').encode()
    reduced = gzip.compress(canonical, compresslevel=9, mtime=0)
    assert len(reduced) <= 8 * 1024**2
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    metadata = root / 'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz'
    assert sha(metadata) == '24155466454df53d01eb6cbbfcd56abeb45adb68f74e721af22a9dce41040b2f'
    rna = pd.read_csv(metadata, usecols=KEYS + ['sub_cell_type'], dtype=str, keep_default_na=False)
    index = pd.MultiIndex.from_frame(rna[KEYS])
    assert len(rna) == 1254749 and index.is_unique
    unit = np.load(native / 'private_cell_unit_index.npy', allow_pickle=False)
    u = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t', keep_default_na=False, dtype={'sampleID': str}).set_index('unit_index').loc[range(11701)]
    np.testing.assert_array_equal(u.sampleID.to_numpy(str)[unit], rna.sampleID)
    np.testing.assert_array_equal(u.sub_cell_type.to_numpy(str)[unit], rna.sub_cell_type)
    source_keys = pd.MultiIndex.from_frame(frame[KEYS])
    mapped = index.get_indexer(source_keys)
    matched = mapped >= 0
    valid, unknown, rowflags = classify(frame)
    label = frame.sub_cell_type.str.strip()
    empty_label = label.str.upper().isin({'', 'NA', 'NAN', 'NONE', 'NULL', '<NA>', 'NO_DATA'}).to_numpy()
    mismatch = np.zeros(len(frame), bool)
    mismatch[matched] = label.to_numpy(str)[matched] != rna.sub_cell_type.to_numpy(str)[mapped[matched]]
    rowflags |= (mismatch.astype(np.uint8) << 5) | (empty_label.astype(np.uint8) << 6)
    n = len(rna)
    count = np.bincount(mapped[matched], minlength=n).astype(np.uint32)
    distinct = ~frame.duplicated(SELECTED).to_numpy()
    unique_count = np.bincount(mapped[matched & distinct], minlength=n).astype(np.uint32)
    vv = np.zeros(n, np.uint8)
    uu = vv.copy()
    ff = vv.copy()
    np.bitwise_or.at(vv, mapped[matched], valid[matched])
    np.bitwise_or.at(uu, mapped[matched], unknown[matched])
    np.bitwise_or.at(ff, mapped[matched], rowflags[matched])
    present = count > 0
    selected = np.flatnonzero(present).astype(np.int32)
    packed = io.BytesIO()
    np.savez_compressed(packed, native_row=selected, valid_component_union=vv[selected], unknown_component_union=uu[selected], row_annotation_flags=ff[selected], source_rows=count[selected], distinct_nonsequence_rows=unique_count[selected])
    assert len(packed.getvalue()) <= 2 * 1024**2
    (out / 'private_complete_selected_nonsequence_annotation.csv.gz').write_bytes(reduced)
    (out / 'private_native_annotation_evidence_sidecar.npz').write_bytes(packed.getvalue())
    summary = []

    def describe(mask):
        return dict(native_cells=int(mask.sum()), TCR_record_present=int((mask & present).sum()), no_TCR_record=int((mask & ~present).sum()), present_all_six_fields_missing=int((mask & present & (vv == 0) & (uu == 0)).sum()), any_expected_component=int((mask & present & (vv > 0)).sum()), nonstandard_annotation_format=int((mask & (uu > 0)).sum()), complete_alpha_in_any_row=int((mask & ((ff & 1) > 0)).sum()), complete_beta_in_any_row=int((mask & ((ff & 2) > 0)).sum()), both_complete_in_same_row=int((mask & ((ff & 4) > 0)).sum()), alpha_JC_in_any_row=int((mask & ((ff & 8) > 0)).sum()), beta_JC_in_any_row=int((mask & ((ff & 16) > 0)).sum()), author_label_disagreement=int((mask & ((ff & 32) > 0)).sum()), TCR_author_label_missing=int((mask & ((ff & 64) > 0)).sum()), repeated_key=int((mask & (count > 1)).sum()), nonidentical_nonsequence_rows=int((mask & (unique_count > 1)).sum()))

    patterns = []
    labels = rna.sub_cell_type.to_numpy(str)
    for name in sorted(set(labels)):
        mask = labels == name
        summary.append(dict(source='GSE243013', RNA_author_type=name, **describe(mask)))
        for pattern in range(64):
            chosen = mask & present & (vv == pattern)
            patterns.append(dict(source='GSE243013', RNA_author_type=name, valid_component_pattern=pattern, TCR_present_cells=int(chosen.sum()), cells_with_nonstandard_component=int((chosen & (uu > 0)).sum())))
    assert len(summary) == 51 and len(patterns) == 3264
    pd.DataFrame(summary).to_csv(out / 'private_all51_RNA_author_annotation_coverage.tsv', sep='\t', index=False)
    pd.DataFrame(patterns).to_csv(out / 'private_all51_by64_component_patterns.tsv', sep='\t', index=False)
    hard = dict(np.load(root / 'NSCLC_gamma_delta_RNA_candidates_20261019_v1/private_candidate_identity_sidecar.npz', allow_pickle=False))
    bits = np.zeros(n, np.uint8)
    bits[hard['native_row']] = hard['rule_bits']
    channels = []
    for j, rule in enumerate(RULES):
        channels.append(dict(source='GSE243013', channel=rule, class_status='RAW_RULE_CANDIDATE', **describe((bits & (1 << j)) > 0)))
    pilot = dict(np.load(root / 'NSCLC_gamma_scGate_fixed_pilot_20261019_v1/private_complete_pilot_candidate_classification_sidecar.npz', allow_pickle=False))
    for j, model in enumerate(pilot['models']):
        for state in ['Pure', 'Impure', 'UNMEASURABLE', 'NOT_ASSESSED']:
            mask = np.zeros(n, bool)
            if state == 'NOT_ASSESSED':
                mask[:] = True
                mask[pilot['native_row']] = False
            else:
                mask[pilot['native_row'][pilot['classes'][:, j] == state]] = True
            channels.append(dict(source='GSE243013', channel=str(model), class_status=state, **describe(mask)))
    assert len(channels) == 20
    pd.DataFrame(channels).to_csv(out / 'private_all20_existing_candidate_channel_coverage.tsv', sep='\t', index=False)
    provenance = dict(full_processed_annotation_sha256=hashlib.sha256(blob).hexdigest(), full_compressed_bytes=len(blob), source_HTTP_headers=source_headers, original_gzip_CRC_validated=True, source_selected_metadata_rows=len(frame), full_selected_metadata_canonical_sha256=hashlib.sha256(canonical).hexdigest(), source_first2MiB_sha256=hashlib.sha256(blob[:2 * 1024**2]).hexdigest(), source_rows_exact_RNA_key_match=int(matched.sum()), source_rows_without_exact_RNA_key=int((~matched).sum()), source_rows_repeated_exact_selected_metadata=int(frame.duplicated(SELECTED).sum()), source_distinct_keys=int(source_keys.nunique()), RNA_keys_unique=True, RNA_metadata_sha256=sha(metadata), native_cell_unit_sha256=sha(native / 'private_cell_unit_index.npy'), native_unit_table_sha256=sha(native / 'private_all_native_units_and_quality.tsv'), matched_native_cells=len(selected), source_nine_fields_no_CDR3_or_clonotype=True, flag_bits={'0': 'alpha_VJC_same_original_row', '1': 'beta_VJC_same_original_row', '2': 'both_VJC_same_original_row', '3': 'alpha_JC_same_original_row', '4': 'beta_JC_same_original_row', '5': 'any_author_label_disagreement', '6': 'any_missing_TCR_author_label'}, no_productivity_or_gamma_truth=True)
    save('private_complete_source_and_bridge_provenance.json', provenance)
    assert sum(p.stat().st_size for p in out.iterdir() if p.is_file()) < 11 * 1024**2
    receipt = dict(status='COMPLETE_FULL_PROCESSED_NONSEQUENCE_ALPHA_BETA_ANNOTATION_EXACT_RNA_BRIDGE_AUDIT_PENDING', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), rule_module_sha256=sha(Path(__file__).with_name('gamma_alpha_beta_annotation_rules_v1.py')), original_author_types=51, pattern_positions=3264, candidate_channel_positions=20, resource_guard=guard, raw_expression_or_raw_sequencing_downloaded=False, sequence_columns_used_or_written=False, original_full_processed_file_only_in_RAM=True, public_scientific_results=False, files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in out.iterdir() if p.is_file()})
    save('private_terminal.json', receipt)
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
