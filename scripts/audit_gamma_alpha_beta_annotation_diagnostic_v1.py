"""Full selected-metadata reparse and bounded independent public-source prefix."""
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
from gamma_alpha_beta_annotation_rules_v1 import SELECTED, GENES, KEYS, PREFIX, MISSING
from gamma_delta_RNA_rules_v1 import RULES
from qualify_public_GSE243013_TCR_header_v1 import URL
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    guard = resource(root, 1024**2)
    out = root / 'gamma_alpha_beta_annotation_diagnostic_audit_20261019_v1'
    assert not out.exists()
    out.mkdir(mode=0o700)
    src = root / 'gamma_alpha_beta_annotation_diagnostic_20261019_v1'
    terminal = wait(root, out, 'measure_gamma_alpha_beta_annotation_diagnostic_v1', src / 'private_terminal.json', {'COMPLETE_FULL_PROCESSED_NONSEQUENCE_ALPHA_BETA_ANNOTATION_EXACT_RNA_BRIDGE_AUDIT_PENDING'})
    for name, info in terminal['files'].items():
        assert sha(src / name) == info['sha256']
    p = json.loads((src / 'private_complete_source_and_bridge_provenance.json').read_text())
    raw = gzip.decompress((src / 'private_complete_selected_nonsequence_annotation.csv.gz').read_bytes())
    assert hashlib.sha256(raw).hexdigest() == p['full_selected_metadata_canonical_sha256']
    d = pd.read_csv(io.BytesIO(raw), dtype=str, keep_default_na=False, engine='python')
    assert list(d) == SELECTED and len(d) == p['source_selected_metadata_rows']
    assert len(d.drop_duplicates(SELECTED)) == len(d) - p['source_rows_repeated_exact_selected_metadata']
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    metadata = root / 'GSE243013_NSCLC_immune_scRNA_metadata.csv.gz'
    assert sha(metadata) == p['RNA_metadata_sha256']
    rna = pd.read_csv(metadata, usecols=KEYS + ['sub_cell_type'], dtype=str, keep_default_na=False)
    lookup = {(a, b): i for i, (a, b) in enumerate(rna[KEYS].itertuples(index=False, name=None))}
    assert len(lookup) == len(rna) == 1254749
    mapped = np.fromiter((lookup.get((a, b), -1) for a, b in d[KEYS].itertuples(index=False, name=None)), dtype=np.int64, count=len(d))
    good = mapped >= 0
    assert good.sum() == p['source_rows_exact_RNA_key_match'] and (~good).sum() == p['source_rows_without_exact_RNA_key']
    assert len(set(d[KEYS].itertuples(index=False, name=None))) == p['source_distinct_keys']
    # Independent scalar syntax classifier and Boolean components.
    valid = np.zeros((len(d), 6), bool)
    other = valid.copy()
    for j, (field, prefix) in enumerate(zip(GENES, PREFIX)):
        for i, value in enumerate(d[field]):
            text = value.strip().upper()
            if text in MISSING:
                continue
            if text.startswith(prefix):
                valid[i, j] = True
            else:
                other[i, j] = True
    alpha = np.logical_and.reduce(valid[:, :3], axis=1)
    beta = np.logical_and.reduce(valid[:, 3:], axis=1)
    missing_label = d.sub_cell_type.map(lambda v: v.strip().upper() in MISSING).to_numpy()
    disagreement = np.zeros(len(d), bool)
    expected_labels = rna.sub_cell_type.to_numpy(str)
    disagreement[good] = d.sub_cell_type.str.strip().to_numpy(str)[good] != expected_labels[mapped[good]]
    flags = np.column_stack([alpha, beta, alpha & beta, valid[:, 1] & valid[:, 2], valid[:, 4] & valid[:, 5], disagreement, missing_label])
    data = {'native_row': mapped}
    for j in range(6):
        data['v' + str(j)] = valid[:, j]
        data['u' + str(j)] = other[:, j]
    for j in range(7):
        data['f' + str(j)] = flags[:, j]
    byrow = pd.DataFrame(data)
    matched = byrow.loc[good]
    group = matched.groupby('native_row', sort=True)
    g = group.max()
    rows = g.index.to_numpy(np.int32)
    count = group.size().to_numpy(np.uint32)
    different = byrow.loc[good & ~d.duplicated(SELECTED).to_numpy()].groupby('native_row').size().reindex(g.index, fill_value=0).to_numpy(np.uint32)
    vv = g[['v' + str(j) for j in range(6)]].to_numpy(np.uint16) @ (2**np.arange(6, dtype=np.uint16))
    uu = g[['u' + str(j) for j in range(6)]].to_numpy(np.uint16) @ (2**np.arange(6, dtype=np.uint16))
    ff = g[['f' + str(j) for j in range(7)]].to_numpy(np.uint16) @ (2**np.arange(7, dtype=np.uint16))
    z = dict(np.load(src / 'private_native_annotation_evidence_sidecar.npz', allow_pickle=False))
    for name, value in [('native_row', rows), ('source_rows', count), ('distinct_nonsequence_rows', different), ('valid_component_union', vv), ('unknown_component_union', uu), ('row_annotation_flags', ff)]:
        np.testing.assert_array_equal(z[name], value)
    assert len(rows) == p['matched_native_cells']
    n = len(rna)
    present = np.zeros(n, bool)
    present[rows] = True
    valid_all = np.zeros(n, int)
    unknown_all = valid_all.copy()
    flags_all = valid_all.copy()
    nrow = valid_all.copy()
    distinct = valid_all.copy()
    valid_all[rows], unknown_all[rows], flags_all[rows], nrow[rows], distinct[rows] = vv, uu, ff, count, different

    def check(row, mask):
        predicates = dict(native_cells=np.ones(n, bool), TCR_record_present=present, no_TCR_record=~present, present_all_six_fields_missing=present & (valid_all == 0) & (unknown_all == 0), any_expected_component=present & (valid_all > 0), nonstandard_annotation_format=unknown_all > 0, complete_alpha_in_any_row=(flags_all & 1) != 0, complete_beta_in_any_row=(flags_all & 2) != 0, both_complete_in_same_row=(flags_all & 4) != 0, alpha_JC_in_any_row=(flags_all & 8) != 0, beta_JC_in_any_row=(flags_all & 16) != 0, author_label_disagreement=(flags_all & 32) != 0, TCR_author_label_missing=(flags_all & 64) != 0, repeated_key=nrow > 1, nonidentical_nonsequence_rows=distinct > 1)
        for name, condition in predicates.items():
            assert row[name] == np.count_nonzero(mask & condition)

    summary = pd.read_csv(src / 'private_all51_RNA_author_annotation_coverage.tsv', sep='\t').set_index('RNA_author_type')
    labels = rna.sub_cell_type.to_numpy(str)
    assert len(summary) == 51 and set(summary.index) == set(labels)
    pattern = pd.read_csv(src / 'private_all51_by64_component_patterns.tsv', sep='\t').set_index(['RNA_author_type', 'valid_component_pattern'])
    assert len(pattern) == 3264
    for label in summary.index:
        check(summary.loc[label], labels == label)
        for code in range(64):
            mask = (labels == label) & present & (valid_all == code)
            row = pattern.loc[(label, code)]
            assert row.TCR_present_cells == mask.sum() and row.cells_with_nonstandard_component == np.count_nonzero(mask & (unknown_all > 0))
    channels = pd.read_csv(src / 'private_all20_existing_candidate_channel_coverage.tsv', sep='\t').set_index(['channel', 'class_status'])
    assert len(channels) == 20
    hard = dict(np.load(root / 'NSCLC_gamma_delta_RNA_candidates_20261019_v1/private_candidate_identity_sidecar.npz', allow_pickle=False))
    for j, name in enumerate(RULES):
        mask = np.zeros(n, bool)
        mask[hard['native_row'][(hard['rule_bits'] & 2**j) != 0]] = True
        check(channels.loc[(name, 'RAW_RULE_CANDIDATE')], mask)
    pilot = dict(np.load(root / 'NSCLC_gamma_scGate_fixed_pilot_20261019_v1/private_complete_pilot_candidate_classification_sidecar.npz', allow_pickle=False))
    for j, name in enumerate(pilot['models']):
        for state in ['Pure', 'Impure', 'UNMEASURABLE', 'NOT_ASSESSED']:
            mask = np.ones(n, bool) if state == 'NOT_ASSESSED' else np.zeros(n, bool)
            if state == 'NOT_ASSESSED':
                mask[pilot['native_row']] = False
            else:
                mask[pilot['native_row'][pilot['classes'][:, j] == state]] = True
            check(channels.loc[(str(name), state)], mask)
    # Bounded independent parsing of original public records, not whole replay.
    cap = 2 * 1024**2
    with urllib.request.urlopen(urllib.request.Request(URL, headers={'Range': 'bytes=0-2097151'}), timeout=30) as response:
        assert response.status == 206 and response.headers.get('Content-Range') == 'bytes 0-2097151/15940803'
        prefix = response.read(cap + 1)
        assert len(prefix) == cap
    assert hashlib.sha256(prefix).hexdigest() == p['source_first2MiB_sha256']
    k = min(8192, len(d))
    with gzip.GzipFile(fileobj=io.BytesIO(prefix), mode='rb') as stream:
        original = pd.read_csv(stream, usecols=SELECTED, nrows=k, dtype=str, keep_default_na=False, engine='python')[SELECTED]
    pd.testing.assert_frame_equal(original.reset_index(drop=True), d.iloc[:k].reset_index(drop=True))
    proof = dict(independent_original_source_prefix_records=k, prefix_network_bytes=cap, prefix_sha256=hashlib.sha256(prefix).hexdigest(), all_minimized_source_rows_reparsed=True, full_exact_RNA_key_bridge_recomputed=True, all_flags_and51_author_profiles_3264_patterns_20_channels_exact=True, original_full_file_independently_replayed=False, sequence_columns_used_or_exported=False, missing_alpha_beta_does_not_define_gamma=True)
    (out / 'private_complete_audit_scope.json').write_text(json.dumps(proof, indent=2) + '\n')
    receipt = dict(status='ACCEPTED_ALL_NONSEQUENCE_METADATA_EXACT_KEYS_FLAGS_51_GROUPS_3264_PATTERNS_20_CHANNELS_AND_SOURCE_PREFIX8192', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), audit_scope_sha256=sha(out / 'private_complete_audit_scope.json'), resource_guard=guard, scientific_outputs_public=False, formal_identity_or_productivity_accepted=False)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
