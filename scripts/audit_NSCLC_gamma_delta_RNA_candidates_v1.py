"""Independent first8192 raw cells plus complete sidecar and native-axis scope.

The raw proof is deliberately bounded; full-cell independent raw replay is not
claimed. All complete candidate rows and original unit aggregates are checked.
"""
import gzip
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from gamma_delta_RNA_rules_v1 import MARKERS, RULES, digest_array
from external_functional_coverage_inputs_v1 import read_inputs
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); guard = resource(root, 1024**2)
    src = root / 'NSCLC_gamma_delta_RNA_candidates_20261019_v1'
    out = root / 'NSCLC_gamma_delta_RNA_candidate_audit_20261019_v1'; assert not out.exists(); out.mkdir(mode=0o700)
    t = wait(root, out, 'measure_NSCLC_gamma_delta_RNA_candidates_v1', src / 'private_terminal.json', {'COMPLETE_ALL_NATIVE_NSCLC_GAMMA_DELTA_FIXED_RULE_CANDIDATES_RAW_PREFIX_AUDIT_PENDING'})
    assert t['status'] == 'COMPLETE_ALL_NATIVE_NSCLC_GAMMA_DELTA_FIXED_RULE_CANDIDATES_RAW_PREFIX_AUDIT_PENDING'
    for name, d in t['files'].items(): assert sha(src / name) == d['sha256']
    plan = json.loads((src / 'plan_before_any_gamma_delta_RNA_counts.json').read_text())
    assert plan['independent_raw_prefix_cells'] == 8192 and plan['markers'] == MARKERS and plan['rules'] == RULES
    n = 8192; raw = np.zeros((n, 9), np.int32); total = np.zeros(n, np.int64); detected = np.zeros(n, np.int32)
    ds = json.loads((root / 'count_matrix_download_v1/download_state.json').read_text()); matrix = Path(ds['path'])
    assert matrix.stat().st_size == ds['bytes'] and ds['sha256'] == t['raw_matrix_sha256']
    gene = pd.read_csv(root / 'GSE243013_genes.csv.gz', keep_default_na=False).geneSymbol.str.upper().str.strip().tolist()
    positions = {gene.index(name): i for i, name in enumerate(MARKERS)}
    with gzip.open(matrix, 'rt', encoding='ascii') as handle:
        assert handle.readline().strip() == '%%MatrixMarket matrix coordinate real general'
        line = handle.readline()
        while line.startswith('%'): line = handle.readline()
        assert tuple(map(int, line.split())) == (1254749, 31831, 2010550708)
        for d in pd.read_csv(handle, sep=r'\s+', names=['row', 'column', 'count'], dtype=np.float64, chunksize=250000, engine='c'):
            stop = bool(d.row.max() > n)
            d = d[d.row.le(n)]
            assert np.isfinite(d.to_numpy()).all() and np.array_equal(d.to_numpy(), np.floor(d.to_numpy())) and d['count'].gt(0).all()
            rr = d.row.to_numpy(np.int64) - 1; cc = d.column.to_numpy(np.int64) - 1; vv = d['count'].to_numpy(np.int64)
            np.add.at(total, rr, vv); np.add.at(detected, rr, 1)
            for col, k in positions.items():
                choose = cc == col; np.add.at(raw[:, k], rr[choose], vv[choose].astype(np.int32))
            if stop: break
    bits = np.zeros(n, np.uint8)
    for i, row in enumerate(raw):
        for k, minimum in enumerate([1, 2]):
            core = all(int(row[j]) >= minimum for j in range(5))
            supported = sum(int(row[j]) >= minimum for j in range(4)) >= 3 and int(row[4]) >= minimum and max(int(row[5]), int(row[6])) >= minimum
            if core: bits[i] += 1 << (2 * k)
            if supported: bits[i] += 1 << (2 * k + 1)
    prefix = json.loads((src / 'private_prefix8192_raw_and_rule_digests.json').read_text())
    for key, array, dtype in [('marker_counts', raw, '<i4'), ('RNA', total, '<i8'), ('detected', detected, '<i4'), ('rule_bits', bits, 'u1')]:
        assert digest_array(array, dtype) == prefix[key + '_sha256']
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    unit = np.load(native / 'private_cell_unit_index.npy', allow_pickle=False)
    u = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t', keep_default_na=False).sort_values('unit_index')
    z = dict(np.load(src / 'private_candidate_identity_sidecar.npz', allow_pickle=False))
    rows = z['native_row']; assert np.array_equal(rows, np.unique(rows)) and (rows >= 0).all() and (rows < len(unit)).all()
    np.testing.assert_array_equal(z['unit_index'], unit[rows]); assert (z['rule_bits'] > 0).all() and (z['rule_bits'] < 16).all()
    fullbits = np.zeros(len(unit), np.uint8); fullbits[rows] = z['rule_bits']
    np.testing.assert_array_equal(fullbits[:n], bits)
    fullproof = json.loads((src / 'private_all_cell_rule_digest.json').read_text()); assert digest_array(fullbits, 'u1') == fullproof['rule_bits_sha256']
    assert fullproof['core2_subset_core1'] and fullproof['supported2_subset_supported1']
    tab = pd.read_csv(src / 'private_all_author_unit_rule_counts.tsv.gz', sep='\t', keep_default_na=False)
    assert len(tab) == 4 * len(u) and not tab.duplicated(['unit_index', 'rule']).any()
    for j, rule in enumerate(RULES):
        g = tab[tab.rule.eq(rule)].sort_values('unit_index'); assert g.unit_index.tolist() == list(range(len(u)))
        for key in ['sampleID', 'major_cell_type', 'sub_cell_type', 'n_cells']:
            np.testing.assert_array_equal(g['native_cells' if key == 'n_cells' else key], u[key])
        chosen = (z['rule_bits'] & (1 << j)) != 0
        for column, mask in [('candidate_cells', chosen), ('TRDV1_detected_candidate_cells', chosen & z['TRDV1_detected']), ('TRDV2_detected_candidate_cells', chosen & z['TRDV2_detected'])]:
            expected = np.bincount(z['unit_index'][mask], minlength=len(u))
            np.testing.assert_array_equal(g[column], expected)
    calculated = tab.groupby(['major_cell_type', 'sub_cell_type', 'rule'], sort=True)[['native_cells', 'candidate_cells', 'TRDV1_detected_candidate_cells', 'TRDV2_detected_candidate_cells']].sum().reset_index()
    recorded = pd.read_csv(src / 'private_complete_author_type_rule_summary.tsv', sep='\t', keep_default_na=False)
    pd.testing.assert_frame_equal(calculated, recorded, check_dtype=False)
    cov = pd.read_csv(src / 'private_source_candidate_coverage.tsv', sep='\t').set_index('rule')
    for rule, g in tab.groupby('rule'):
        s = g.groupby('sampleID').candidate_cells.sum(); r = cov.loc[rule]
        assert r.source == 'GSE243013' and r.author_records == len(s) == 243 and r.records_candidate_atleast20 == s.ge(20).sum() and r.candidate_cells == g.candidate_cells.sum() and r.reviewed_independent_patients == 0
    axes = read_inputs(root); measurements = []
    for (arm, dataset), axis in sorted(axes['source_axes'].items()):
        mapping = [axis['lookup'].get(g, set()) for g in MARKERS]
        measured = [len(v) == 1 for v in mapping]
        columns = [next(iter(v)) if len(v) == 1 else -1 for v in mapping]
        collision = len(set(v for v in columns if v >= 0)) != sum(v >= 0 for v in columns)
        measurements.append(dict(arm=arm, dataset=dataset, core5_uniquely_measurable=all(measured[:5]) and not collision, supported7_uniquely_measurable=all(measured[:7]) and not collision, TRDV1_measurable=measured[7], TRDV2_measurable=measured[8], marker_columns=columns, source_column_collision=collision, expression_not_read=True))
    pd.DataFrame(measurements).to_csv(out / 'private_complete61_prospective_candidate_gene_scope.tsv', sep='\t', index=False)
    receipt = dict(status='ACCEPTED_COMPLETE_NATIVE_CANDIDATE_SIDECAR_AUTHOR_AGGREGATES_AND_INDEPENDENT_FIRST8192_RAW_CELLS', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), resource_guard=guard, independently_raw_cells=n, independently_raw_marker_positions=n * len(MARKERS), complete_candidate_and_unit_aggregates_checked=True, whole_source_independent_raw_replay=False, original61_source_native_gene_scope_only=True, clinical_or_formal_identity_adoption=False, scientific_outputs_public=False, prospective_source_scope_sha256=sha(out / 'private_complete61_prospective_candidate_gene_scope.tsv'))
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
