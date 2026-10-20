"""Independent depth-bin and native-unit aggregation; no raw RNA remeasurement."""
import csv
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def median(values):
    a = sorted(float(v) for v in values)
    if not a:
        return None
    n = len(a)
    return a[n // 2] if n % 2 else (a[n // 2 - 1] + a[n // 2]) / 2


def check_scalar(row, name, expected):
    if expected is None:
        assert row[name] == ''
    else:
        assert abs(float(row[name]) - expected) <= 1e-10 * max(1, abs(expected))


def main():
    parent = Path(os.environ['PANCANCER_DEPTH_PRODUCER_OUTPUT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    assert sha(parent / 'private_terminal.json') == os.environ['PANCANCER_PARENT_TERMINAL_SHA256']
    terminal = json.loads((parent / 'private_terminal.json').read_text())
    assert terminal['source_sha256'] == 'd47c018d709bc34928b82ae77c26fb44e6fac2ab7d376532805c040df934b456'
    for name, info in terminal['files'].items():
        assert sha(parent / name) == info['sha256']
    plan = json.loads((parent / 'plan_before_depth_descriptions.json').read_text())
    native = root / 'NSCLC_full_native_gene_units_20261016_v1'
    for name, expected in plan['input_SHA256'].items():
        assert sha(native / name) == expected
    hpath = root / 'NSCLC_gamma_delta_RNA_candidates_20261019_v1/private_candidate_identity_sidecar.npz'
    assert sha(hpath) == plan['candidate_SHA256']
    unit = np.load(native / 'private_cell_unit_index.npy', allow_pickle=False)
    RNA = np.load(native / 'private_full_RNA.npy', allow_pickle=False)
    detected = np.load(native / 'private_detected_genes.npy', allow_pickle=False)
    u = pd.read_csv(native / 'private_all_native_units_and_quality.tsv', sep='\t', keep_default_na=False).sort_values('unit_index')
    h = dict(np.load(hpath, allow_pickle=False))
    bits = np.zeros(len(unit), np.uint8)
    bits[h['native_row']] = h['rule_bits']
    np.testing.assert_array_equal(unit[h['native_row']], h['unit_index'])
    assert len(unit) == len(RNA) == len(detected)
    assert len(np.unique(h['native_row'])) == len(h['native_row'])
    assert int(np.count_nonzero(bits)) == int(np.count_nonzero(h['rule_bits']))
    labels = u.sub_cell_type.to_numpy(str)
    roles = np.full(len(u), 'Non_T_NK', dtype=object)
    roles[u.major_cell_type.eq('T/NK cell').to_numpy()] = 'Other_T_NK'
    roles[np.char.startswith(labels, 'NK_')] = 'Author_NK'
    roles[labels == 'T_gdT_TRDV1'] = 'Author_Vdelta1'
    roles[labels == 'T_gdT_TRDV2'] = 'Author_Vdelta2'
    reference = u.major_cell_type.eq('T/NK cell').to_numpy()[unit]
    records = pd.Categorical(u.sampleID).codes
    assert records.max() + 1 == 243
    with (parent / 'private_all448_author_group_depth_descriptions.tsv').open() as handle:
        table = list(csv.DictReader(handle, delimiter='\t'))
    lookup = {(r['scope'], r['author_group'], r['axis'], int(r['depth_bin'])): r for r in table}
    assert len(lookup) == len(table) == 448
    contracts = {r['axis']: r for r in json.loads((parent / 'private_depth_bin_contract.json').read_text())}
    rule_names = [f'original_rule_bit_{j}' for j in range(4)] + ['fixed_RNA_candidate_union']
    flags = [(bits & (1 << j)) != 0 for j in range(4)] + [bits != 0]
    checked = 0
    for axis, values in [('native_total_RNA', RNA), ('native_detected_genes', detected)]:
        sorted_reference = np.sort(values[reference])
        edges = []
        for q in [.25, .5, .75]:
            pos = (len(sorted_reference) - 1) * q
            low = int(np.floor(pos)); high = int(np.ceil(pos))
            edges.append(float(sorted_reference[low] + (sorted_reference[high] - sorted_reference[low]) * (pos - low)))
        assert edges == contracts[axis]['boundaries']
        assert int(reference.sum()) == contracts[axis]['reference_cells']
        # Comparison sums retain ties together, independent of searchsorted.
        bins = np.sum(values[:, None] >= np.asarray(edges)[None, :], axis=1)
        key = unit.astype(np.int64) * 4 + bins
        counts = np.zeros((len(u) * 4, 6), np.int64)
        counts[:, 0] = np.bincount(key, minlength=len(u) * 4)
        for j, flag in enumerate(flags):
            counts[:, j + 1] = np.bincount(key[flag], minlength=len(u) * 4)
        counts = counts.reshape(len(u), 4, 6)
        quality = []
        for v in [RNA, detected]:
            assert int(np.sum(v, dtype=np.int64)) < 2**53
            sums = np.bincount(key, weights=v, minlength=len(u) * 4)
            assert np.array_equal(sums, np.floor(sums))
            quality.append(sums.astype(np.int64).reshape(len(u), 4))
        assert counts[:, :, 0].sum() == len(unit)
        for scope, groups in [('all51_author_types', labels), ('five_disjoint_author_roles', roles)]:
            for group in sorted(set(groups)):
                take = groups == group
                per_record = np.zeros((243, 4, 6), np.int64)
                np.add.at(per_record, records[take], counts[take])
                for b in range(4):
                    row = lookup[(scope, group, axis, b + 1)]
                    assert row['source'] == 'GSE243013'
                    n = int(per_record[:, b, 0].sum())
                    assert n == int(row['cells'])
                    assert int(np.count_nonzero(per_record[:, b, 0])) == int(row['author_records_with_cells'])
                    eligible = per_record[:, b, 0] >= 20
                    assert int(eligible.sum()) == int(row['author_records_atleast20'])
                    for j, name in enumerate(rule_names):
                        k = int(per_record[:, b, j + 1].sum())
                        assert k == int(row[name + '_cells'])
                        check_scalar(row, name + '_fraction', k / n if n else None)
                        check_scalar(row, name + '_equal_record_median', median(per_record[eligible, b, j + 1] / per_record[eligible, b, 0]))
                    for name, sums in zip(['RNA_mean', 'detected_genes_mean'], quality):
                        check_scalar(row, name, int(sums[take, b].sum()) / n if n else None)
                    checked += 1
    assert checked == 448
    out.mkdir(mode=0o700)
    result = dict(status='ACCEPTED_ALL_DEPTH_BINS_NATIVE_UNIT_RECONSTRUCTION_AND_SOURCE_AGGREGATES', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), parent_terminal_sha256=sha(parent / 'private_terminal.json'), all_aggregate_positions_checked=checked, independent_sorted_quantiles_and_comparison_bins=True, independent_native_unit_then_author_record_aggregation=True, all_rule_numerators_denominators_quality_means_and_equal_record_medians=True, missing_empty_bins_preserved=True, no_new_raw_RNA_read_or_identity_adoption=True, individual_arrays_persisted=False)
    p = out / 'private_independent_acceptance.json'
    p.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], acceptance_sha256=sha(p))))


if __name__ == '__main__':
    main()
