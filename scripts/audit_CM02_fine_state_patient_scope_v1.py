"""Counter/Fraction replay of every fine-reference coverage and version result.

Shares the pinned, previously independently audited reference join materializer;
recomputes eligibility, fractions, coverage, quantiles and original counts.
"""
from collections import Counter, defaultdict
from datetime import datetime, timezone
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import csv
import hashlib
import importlib.util
import json
import math
import os


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def quantile(values, probability):
    a = sorted(values)
    index = (len(a) - 1) * probability
    low, high = math.floor(index), math.ceil(index)
    return a[low] + (a[high] - a[low]) * (index - low)


def main():
    parent = Path(os.environ['PANCANCER_PATIENT_SCOPE_OUTPUT'])
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    terminal = json.loads((parent / 'private_terminal.json').read_text())
    for name, value in terminal['files'].items():
        assert sha(parent / name) == value['SHA256']
    source = parent.parent / 'qualify_CM02_fine_state_patient_scope_v1.py'
    assert sha(source) == terminal['source_SHA256']
    spec = importlib.util.spec_from_file_location('accepted_input_materializer', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    tab, cur, unique, shared, matrices, table = module.load_inputs()
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    oldpath = root / 'HNSCC_author_state_count_layer_preflight_20261016_v1/private_all57_original_author_counts.tsv'
    with oldpath.open() as f:
        original = {r['library_id']: {k: int(r[k]) for k in module.OLD} for r in csv.DictReader(f, delimiter='\t')}
    for r in tab.to_dict('records'):
        assert all(int(r[k]) == original[r['library_id']][k] for k in module.OLD)
    metadata = {r['library_id']: r for r in tab.to_dict('records')}
    cells = Counter((r.library_id, r.label) for r in unique.itertuples())
    m = matrices['unique_current_reference']
    for library, values in m.iterrows():
        for state, n in values.items():
            assert int(n) == cells[library, state]
    scopes = {r['phase']: r for r in json.loads((parent / 'private_all3_phase_reference_scope.json').read_text())}
    actual = json.loads((parent / 'private_all84_fine_state_opportunities.json').read_text())
    indexed = {(r['phase'], r['CD4_author_state'], r['CD8_author_state']): r for r in actual}
    assert len(indexed) == len(actual) == 84
    for phase in ['pre', 'post', 'post2']:
        candidates = [l for l, r in metadata.items() if r['timepoint_registry'] == phase]
        patients = Counter(metadata[l]['patient_registry'] for l in candidates)
        ready = []
        for l in candidates:
            n4 = sum(cells[l, k] for k in module.CD4)
            n8 = sum(cells[l, k] for k in module.CD8)
            r = metadata[l]
            if n4 >= 20 and n8 >= 20 and patients[r['patient_registry']] == 1 and r['neo_type'] in ['D', 'D+T'] and r['HPV/p16'] in ['positive', 'negative']:
                ready.append((l, n4, n8))
        s = scopes[phase]
        assert s['reviewed_libraries'] == len(candidates) and s['reviewed_patients'] == len(patients)
        assert s['duplicate_patient_phase_libraries'] == sum(v for v in patients.values() if v > 1)
        assert s['eligible_for_future_descriptive_structure'] == len(ready)
        assert s['both_reference_parents_atleast20'] == sum(sum(cells[l, k] for k in module.CD4) >= 20 and sum(cells[l, k] for k in module.CD8) >= 20 for l in candidates)
        for a in s['by_arm']:
            assert a['eligible'] == sum(metadata[l]['neo_type'] == a['arm'] for l, _, _ in ready)
        if phase in ['pre', 'post']:
            for a, b in product(module.CD4, module.CD8):
                r = indexed[phase, a, b]
                v4 = len({Fraction(cells[l, a], n4) for l, n4, _ in ready})
                v8 = len({Fraction(cells[l, b], n8) for l, _, n8 in ready})
                assert r['eligible_author_patients'] == len(ready)
                assert r['CD4_unique_fractions'] == v4 and r['CD8_unique_fractions'] == v8
                assert r['CD4_positive_patients'] == sum(cells[l, a] > 0 for l, _, _ in ready)
                assert r['CD8_positive_patients'] == sum(cells[l, b] > 0 for l, _, _ in ready)
                assert r['status'].startswith('REFERENCE_COUNT_SCOPE_READY') == (len(ready) >= 20 and v4 > 1 and v8 > 1)
    total = Counter((r.library_id, r.node_id) for r in cur.itertuples())
    ref = Counter((r.library_id, r.node_id) for r in unique.itertuples())
    overlap = Counter((r.library_id, r.node_id) for r in cur.itertuples() if r.position in shared)
    groups = defaultdict(list)
    for (l, n), v in total.items():
        groups[metadata[l]['timepoint_registry'], metadata[l]['neo_type'], n].append((l, v))
    coverage = json.loads((parent / 'private_all_node_phase_arm_patient_coverage.json').read_text())
    assert set(groups) == {(r['phase'], r['arm'], r['current_node']) for r in coverage}
    for r in coverage:
        values = [(l, v) for l, v in groups[r['phase'], r['arm'], r['current_node']] if v >= 20]
        assert len(values) == r['libraries_node_atleast20']
        n = r['current_node']
        fs = [ref[l, n] / v for l, v in values]
        if fs:
            expected = [quantile(fs, p) for p in [0, .25, .5, .75, 1]]
            assert all(abs(a-b) < 1e-12 for a, b in zip(expected, r['unique_reference_fraction_min_q25_median_q75_max']))
        else:
            assert r['unique_reference_fraction_min_q25_median_q75_max'] is None
        assert r['libraries_reference_below_half'] == sum(f < .5 for f in fs)
        assert r['all_current_cells'] == sum(v for l, v in values)
        assert r['unique_reference_cells'] == sum(ref[l, n] for l, _ in values)
        assert r['overlap_cells'] == sum(overlap[l, n] for l, _ in values)
        assert r['unannotated_cells'] == r['all_current_cells'] - r['unique_reference_cells'] - r['overlap_cells']
    comparisons = json.loads((parent / 'private_all_count_version_comparisons.json').read_text())
    assert len(comparisons) == 3 * 4 * (10 + 63)
    for r in comparisons:
        libs = [l for l, v in metadata.items() if r['phase'] == 'all' or v['timepoint_registry'] == r['phase']]
        counts = matrices[r['count_scope']]
        pairs = [(sum(int(counts.at[l, k]) for k in r['fine_labels']), original[l][r['original_count_label']]) for l in libs]
        differences = [a-b for a, b in pairs]
        assert r['reviewed_libraries'] == len(libs)
        assert r['exact_count_libraries'] == sum(v == 0 for v in differences)
        assert r['exact_nonzero_libraries'] == sum(a == b and (a > 0 or b > 0) for a, b in pairs)
        assert r['all_libraries_exact'] == all(v == 0 for v in differences)
        assert r['original_total'] == sum(b for _, b in pairs) and r['fine_total'] == sum(a for a, _ in pairs)
        assert r['sum_absolute_difference'] == sum(abs(v) for v in differences)
        assert r['maximum_absolute_difference'] == (max(abs(v) for v in differences) if differences else None)
    result = dict(status='ACCEPTED_ALL84_FINE_REFERENCE_OPPORTUNITIES_ALL_COVERAGE_AND876_COUNT_VERSION_COMPARISONS', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), parent_terminal_SHA256=sha(parent / 'private_terminal.json'), shared_previously_audited_reference_join_materializer=True, all_unique_reference_library_state_counts_replayed_by_Counter=True, original_counts_independently_read_from_prior_count_table=True, original_and_matched_scope_version_checks_share_pinned_input_count_matrices=True, all_patient_scope_fractions_replayed_with_Fraction=True, all_coverage_quantiles_replayed_with_sorted_linear_interpolation=True, no_new_independence_or_functional_truth_claim=True, no_individual_rows_or_arrays_exported=True)
    out.mkdir(mode=0o700)
    (out / 'private_independent_acceptance.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], acceptance_SHA256=sha(out / 'private_independent_acceptance.json'))))


if __name__ == '__main__':
    main()
