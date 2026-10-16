"""Independent native marginal, cell-set and matched-case cofactor audit.

All selected records and all diagnostic component sets are retained. No new
clinical test, programme-weight export, identity change or protein inference.
"""
import gzip
import hashlib
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8388608), b''):
            h.update(block)
    return h.hexdigest()


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    out = root / 'NSCLC_candidate_lineage_cofactor_RNA_20261016_v1'
    receipt_path = root / 'HNSCC_receiver_followup_job_control_20261016_v2/prepare_NSCLC_candidate_lineage_cofactor_RNA_v1/launch_receipt.json'
    launch = json.loads(receipt_path.read_text())
    deadline = time.monotonic() + 14400
    while not (out / 'private_terminal.json').exists():
        proc = Path('/proc') / str(launch['pid']) / 'stat'
        assert proc.exists(), 'Original producer is missing without terminal'
        state = proc.read_text().rsplit(')', 1)[1].split()
        assert state[19] == str(launch['proc_start_ticks']) and state[0] != 'Z'
        assert time.monotonic() < deadline
        time.sleep(15)
    dest = out / 'private_independent_audit.json'
    assert not dest.exists()
    terminal = json.loads((out / 'private_terminal.json').read_text())
    plan = json.loads((out / 'plan_before_native_lineage_cofactor_read.json').read_text())
    assert plan['code_sha256'] == launch['source_sha256']
    for filename, checksum in terminal['files'].items():
        assert sha(out / filename) == checksum
    original = root / 'NSCLC_automatic_ligand_response_results_20261016_v1'
    assert sha(original / 'private_independent_audit.json') == plan['result_audit_sha256']
    result_path = original / 'private_all_unique_ligand_response_results.tsv'
    assert sha(result_path) == plan['result_sha256'] == terminal['result_sha256']
    selected = pd.read_csv(result_path, sep='\t')
    selected = selected[selected.BH_q_complete_unique_family.lt(.05)]
    assert selected.record.tolist() == plan['selected_records']
    rules_path = root / 'NSCLC_automatic_ligand_molecular_review_20261016_v1/private_public_literature_component_rules.json'
    assert sha(rules_path) == plan['molecular_rules_sha256']
    rules = json.loads(rules_path.read_text())['ligand_rules']
    native_root = root / 'NSCLC_full_native_gene_units_20261016_v1'
    axis_native = pd.Index(pd.read_csv(native_root / 'private_native_gene_axis.tsv', sep='\t').gene)
    assert sha(native_root / 'private_native_gene_axis.tsv') == plan['native_axis_sha256']
    assert sha(native_root / 'private_cell_unit_index.npy') == plan['cell_unit_sha256']
    units = pd.read_csv(native_root / 'private_all_native_units_and_quality.tsv', sep='\t')
    unit = np.load(native_root / 'private_cell_unit_index.npy', mmap_mode='r')
    assert units.unit_index.tolist() == list(range(len(units)))
    assert np.array_equal(np.bincount(unit, minlength=len(units)), units.n_cells)
    axis = pd.read_csv(out / 'private_native_panel_axis.tsv', sep='\t')
    assert axis.gene.is_unique
    assert axis.native_column.tolist() == axis_native.get_indexer(axis.gene).tolist()
    assert axis.measured.tolist() == (axis.native_column >= 0).tolist()
    markers = plan['public_marker_reference_sets']
    assert all(axis.set_index('gene').measured.loc[genes].all() for genes in markers.values()), 'A marker reference is unmeasured; do not interpret its zeros'
    panel_genes = set(g for genes in markers.values() for g in genes)
    for row in selected.itertuples():
        panel_genes.add(row.ligand)
        panel_genes.update(g for receptor in row.qualified_receptors.split('|') for g in receptor.split('_'))
        rule = rules[row.ligand]
        panel_genes.update(g for role in ['sender_sets', 'receiver_sets'] for group in rule[role] for g in group)
        panel_genes.update(rule['auxiliary_genes'])
    assert axis.gene.tolist() == sorted(panel_genes)
    counts = np.load(out / 'private_native_capped2_panel.npy', mmap_mode='r')
    assert counts.dtype == np.uint8 and counts.shape == (len(unit), len(axis))
    assert int(counts.max()) <= 2
    positive = np.load(native_root / 'positive_cell_counts.npy', mmap_mode='r')
    for j, row in enumerate(axis.itertuples()):
        observed = pd.Series(unit[counts[:, j] > 0]).value_counts().reindex(range(len(units)), fill_value=0).to_numpy(np.int64)
        if row.measured:
            assert np.array_equal(observed, positive[:, int(row.native_column)])
        else:
            assert not observed.any()
    cells = np.load(out / 'private_prespecified_prefix64_cells.npy')
    assert cells.tolist() == np.sort(np.random.default_rng(2026101613).choice(min(4096, len(unit)), 64, replace=False)).tolist()
    lookup = {int(cell): j for j, cell in enumerate(cells)}
    columns = {int(row.native_column): j for j, row in enumerate(axis.itertuples()) if row.measured}
    expected = np.zeros((len(cells), len(axis)), np.uint8)
    download = json.loads((root / 'count_matrix_download_v1/download_state.json').read_text())
    assert terminal['source_compressed_sha256'] == download['sha256'] == plan['matrix_expected_sha256']
    with gzip.open(download['path'], 'rt') as stream:
        assert stream.readline().startswith('%%MatrixMarket matrix coordinate')
        line = stream.readline()
        while line.startswith('%'):
            line = stream.readline()
        nr, ng, nnz = map(int, line.split())
        assert (nr, ng, nnz) == (len(unit), len(axis_native), terminal['source_coordinates'])
        reader = pd.read_csv(stream, sep=r'\s+', header=None, names=['row', 'col', 'count'], dtype=np.int64, engine='c', chunksize=1000000)
        for frame in reader:
            sub = frame[frame.row.isin(cells + 1) & frame.col.isin(np.asarray(list(columns)) + 1)]
            for row in sub.itertuples(index=False):
                expected[lookup[int(row.row) - 1], columns[int(row.col) - 1]] = min(int(row.count), 2)
            if int(frame.row.iloc[-1]) > int(cells.max()) + 1:
                break
        reader.close()
    assert np.array_equal(expected, counts[cells])
    sets = {(gene, level): set(np.flatnonzero(counts[:, j] >= level).tolist()) for j, gene in enumerate(axis.gene) for level in [1, 2]}
    marker_sets = {}
    for name, genes in markers.items():
        for level in [1, 2]:
            overlaps = set()
            for j, first in enumerate(genes):
                for second in genes[j + 1:]:
                    overlaps.update(sets[first, level].intersection(sets[second, level]))
            marker_sets[name, level] = overlaps
    diagnostics, summaries, cache = [], [], {}
    designs = root / 'NSCLC_automatic_ligand_response_preflight_20261016_v1/private_frozen_designs'
    measured_axis = axis.set_index('gene').measured
    for row in selected.itertuples():
        z = np.load(designs / (row.record + '.npz'), allow_pickle=False)
        for role in ['sender', 'receiver']:
            indices = z[role + '_units'].astype(int)
            assert len(set(indices.tolist())) == len(indices)
            subtype = getattr(row, role + '_author_subtype')
            matched = units.iloc[indices]
            assert matched.sub_cell_type.eq(subtype).all()
            assert matched.sampleID.astype(str).tolist() == z['cases'].astype(str).tolist()
            n = matched.n_cells.to_numpy(int)
            core = 'myeloid_lineage' if subtype.startswith('M\u03c6_') else 'T_lineage'
            groups = [(row.ligand,)] if role == 'sender' else [tuple(r.split('_')) for r in row.qualified_receptors.split('|')]
            groups += list(map(tuple, rules[row.ligand][role + '_sets']))
            groups = sorted(set(tuple(sorted(g)) for g in groups))
            for group in groups:
                measured = all(g in measured_axis.index and measured_axis.loc[g] for g in group)
                for level in [1, 2]:
                    common = dict(record=row.record, role=role, author_subtype=subtype, group='|'.join(group), threshold=level)
                    if not measured:
                        summaries.append(dict(**common, measured=False, status='REQUIRED_NATIVE_GENE_NOT_MEASURED', native_cells=int(n.sum()), cases=len(indices)))
                        continue
                    key = group, level
                    if key not in cache:
                        active = set.intersection(*(sets[g, level] for g in group))
                        groups_of_cells = {'GOI_positive_cells': active}
                        groups_of_cells.update({name + '_co_positive_cells': active.intersection(marker_sets[name, level]) for name in markers})
                        cache[key] = {}
                        for name, cellset in groups_of_cells.items():
                            rows = np.fromiter(cellset, np.int64)
                            cache[key][name] = pd.Series(unit[rows]).value_counts().reindex(range(len(units)), fill_value=0).to_numpy(np.int64)
                    values = {name: value[indices] for name, value in cache[key].items()}
                    for j, u in enumerate(indices):
                        diagnostics.append(dict(**common, unit_index=int(u), native_cells=int(n[j]), measured=True, **{name: int(v[j]) for name, v in values.items()}))
                    totals = {name: int(v.sum()) for name, v in values.items()}
                    positive_cells = totals['GOI_positive_cells']
                    summaries.append(dict(**common, measured=True, native_cells=int(n.sum()), cases=len(indices), positive_cases=int(np.count_nonzero(values['GOI_positive_cells'])), GOI_positive_fraction=positive_cells / int(n.sum()), core_lineage_marker_set=core, core_lineage_among_GOI_positive=totals[core + '_co_positive_cells'] / positive_cells if positive_cells else None, **totals))
    for expected_frame, name, keys in [
        (pd.DataFrame(diagnostics), 'private_all_matched_unit_lineage_cofactor_counts.tsv.gz', ['record', 'role', 'group', 'threshold', 'unit_index']),
        (pd.DataFrame(summaries), 'private_all_record_lineage_cofactor_summaries.tsv', ['record', 'role', 'group', 'threshold']),
    ]:
        observed = pd.read_csv(out / name, sep='\t')
        assert not observed.duplicated(keys).any()
        if 'status' not in expected_frame and 'status' in observed:
            assert observed.status.isna().all()
            observed = observed.drop(columns='status')
        assert set(expected_frame.columns) == set(observed.columns)
        a = expected_frame.set_index(keys).sort_index()
        b = observed.set_index(keys).sort_index().reindex(columns=a.columns)
        pd.testing.assert_frame_equal(a, b, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-12)
    result = dict(status='PASS_ALL_NATIVE_UNIT_MARGINALS_ALL_MATCHED_CELL_SET_INTERSECTIONS_ALL_SUMMARIES_AND_FIXED64_RAW_CELL_REPLAY', finished_UTC=datetime.now(timezone.utc).isoformat(), audit_code_sha256=sha(__file__), producer_terminal_sha256=sha(out / 'private_terminal.json'), records=len(selected), native_cells=len(unit), native_units=len(units), panel_columns=len(axis), matched_unit_diagnostic_rows=len(diagnostics), summary_rows=len(summaries), raw_replay_cells=len(cells), all_native_cells_independently_raw_replayed=False, all_marker_reference_genes_measured=True, learned_programme_gene_lists_exported=False, new_P_BH_or_identity_changes=False, signaling_or_protein_interaction_proved=False)
    dest.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
