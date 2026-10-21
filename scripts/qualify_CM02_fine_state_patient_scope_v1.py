"""Outcome-free patient coverage and count-version checks of author references.

Cell/patient arrays stay in RAM. Source-level aggregates only are written.
Reference-subset proportions are not whole-population abundance estimates.
"""
from collections import Counter
from datetime import datetime, timezone
from itertools import combinations, product
from pathlib import Path
import hashlib
import importlib.util
import json
import os

import numpy as np
import pandas as pd
import zarr
from anndata._io.specs import read_elem

CD4 = ['4-1BB+ Treg', 'HLA-DR+ Treg', 'LAG3+ Treg', 'Prolif Treg', 'CXCL13+ CD4T', 'IFN-stim CD4T', 'Tn/Tmem']
CD8 = ['CXCL13+Tex', 'Tdys', 'GZMK+T', 'Tmem', 'Tprolif', 'Tstr']
OLD = CD4[:4] + ['Tn/Tmem', 'IFN-stim CD4T', 'GZMK+T', 'Tmem', 'Teff', 'Tstr', 'Tprolif']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    parent = Path(os.environ['PANCANCER_PARENT_OUTPUT'])
    terminal = json.loads((parent / 'private_terminal.json').read_text())
    source = parent.parent / 'probe_CM02_author_library_barcode_bridge_v1.py'
    assert sha(source) == terminal['source_SHA256']
    for name, v in terminal['files'].items():
        assert sha(parent / name) == v['SHA256']
    spec = importlib.util.spec_from_file_location('pinned_reference_reader', source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    table = root / 'HNSCC_author_clinical_table_user_20261014_v1/private_key_and_response_audit_v1/all57_private_author_table_join.tsv'
    cols = list(dict.fromkeys(['library_id', 'patient_registry', 'timepoint_registry', 'neo_type', 'HPV/p16', 'Cells analyzed'] + module.PATIENT_FIELDS + module.TIME_FIELDS + OLD))
    tab = pd.read_csv(table, sep='\t', usecols=cols, dtype=str, keep_default_na=False)
    assert tab.library_id.is_unique
    atlas = root.parents[1] / 'runs/20261013_immunotherapy_atlas_v1'
    z = zarr.open_group(str(atlas / 'outputs/immune_validation_authority_v1.zarr'), mode='r')
    mask = np.asarray(read_elem(z['obs/dataset'])).astype(str) == 'GSE286827'
    cur = pd.DataFrame({k: np.asarray(read_elem(z['obs'][k])).astype(str)[mask] for k in ['library_id', 'timepoint', 'major', 'node_id', 'original_cell_id']})
    cur['position'] = np.arange(len(cur))
    cur['barcode'] = cur.original_cell_id.map(lambda v: module.barcode(v, 'within_reviewed_library_16nt'))
    assert cur.barcode.notna().all() and not cur.duplicated(['library_id', 'barcode']).any()
    assert set(cur.library_id) == set(tab.library_id)
    opportunities = json.loads((parent / 'private_all18_library_barcode_opportunities.json').read_text())
    cross = json.loads((parent / 'private_all_source_node_state_crosscounts.json').read_text())
    raw, matched = [], []
    for name, field, expected in module.SOURCES[1:]:
        o = next(r for r in opportunities if r['annotation_source'] == name)
        a = module.read_author(name, field, expected)
        lookup = {}
        for r in tab.to_dict('records'):
            key = (module.norm(r[o['patient_alias_field']]), module.norm(r[o['time_alias_field']]))
            assert key not in lookup
            lookup[key] = r['library_id']
        a['library_id'] = [lookup.get((module.norm(p), module.norm(t))) for p, t in zip(a.patient, a.phase)]
        a['barcode'] = a.cell_key.map(lambda v: module.barcode(v, 'within_reviewed_library_16nt'))
        assert a.library_id.notna().all() and a.barcode.notna().all() and not a.duplicated(['library_id', 'barcode']).any()
        kind = 'CD4_reference' if field == 'cd4t_subcluster' else 'CD8_reference'
        assert set(a.label) == set(CD4 if kind == 'CD4_reference' else CD8)
        a['reference'] = kind
        m = cur.merge(a[['library_id', 'barcode', 'label', 'reference']], on=['library_id', 'barcode'], validate='one_to_one')
        got = {(p, g, n, l): int(v) for (p, g, n, l), v in m.groupby(['timepoint', 'major', 'node_id', 'label']).size().items()}
        old = {(r['phase'], r['current_major'], r['current_node'], r['author_label']): r['matched_cells'] for r in cross if r['annotation_source'] == name and r['barcode_mode'] == 'within_reviewed_library_16nt' and r['current_cell_key_field'] == 'original_cell_id'}
        assert got == old
        raw.append(a)
        matched.append(m)
    shared = set(matched[0].position) & set(matched[1].position)
    unique = pd.concat([m[~m.position.isin(shared)] for m in matched], ignore_index=True)
    assert unique.position.is_unique
    scopes = {'original_author_reference': pd.concat(raw, ignore_index=True), 'matched_current_reference': pd.concat(matched, ignore_index=True), 'unique_current_reference': unique}
    matrices = {k: f.groupby(['library_id', 'label']).size().unstack(fill_value=0).reindex(index=tab.library_id, columns=CD4 + CD8, fill_value=0).astype(int) for k, f in scopes.items()}
    return tab, cur, unique, shared, matrices, table


def main():
    out = Path(os.environ['PANCANCER_RAM_OUTPUT'])
    assert out.parent.parent == Path('/dev/shm') and not out.exists()
    out.mkdir(mode=0o700)
    save = lambda name, value: (out / name).write_text(json.dumps(value, indent=2) + '\n')
    save('plan_before_patient_coverage_and_versions.json', dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), CD4_reference_labels=CD4, CD8_reference_labels=CD8, scope='All reviewed libraries; all 7x6 author reference combinations in pre/post. Post2 coverage only.', reference_denominators='within each fine reference file, overlap excluded; not all tissue/CD4/CD8 truth', parent_minimum=20, minimum_patients_for_future_structure=20, whole_node_coverage='all current nodes, phase and treatment, with no reference assigned to missing or overlap cells', version_tests='All shared exact labels; every nonempty subset of six CD8 labels against original wide Teff for all three count scopes; no automatic renaming.', outcome_expression_TCR_values_not_read=True, no_new_clinical_or_pair_effects=True, no_cell_or_patient_arrays_saved=True))
    tab, cur, unique, shared, matrices, table = load_inputs()
    phase_rows, opportunities, cover, versions = [], [], [], []
    matrix = matrices['unique_current_reference']
    metadata = tab.set_index('library_id')
    for phase in ['pre', 'post', 'post2']:
        g = tab[tab.timepoint_registry.eq(phase)].copy()
        x = matrix.loc[g.library_id].to_numpy(float)
        n4, n8 = x[:, :7].sum(axis=1), x[:, 7:].sum(axis=1)
        duplicates = g.patient_registry.duplicated(keep=False).to_numpy()
        eligible = (n4 >= 20) & (n8 >= 20) & ~duplicates & g.neo_type.isin(['D', 'D+T']).to_numpy() & g['HPV/p16'].isin(['positive', 'negative']).to_numpy()
        phase_rows.append(dict(phase=phase, reviewed_libraries=len(g), reviewed_patients=g.patient_registry.nunique(), duplicate_patient_phase_libraries=int(duplicates.sum()), both_reference_parents_atleast20=int(((n4 >= 20) & (n8 >= 20)).sum()), eligible_for_future_descriptive_structure=int(eligible.sum()), by_arm=[dict(arm=a, eligible=int((eligible & g.neo_type.eq(a).to_numpy()).sum())) for a in ['D', 'D+T']]))
        if phase in ['pre', 'post']:
            y = x[eligible] / np.column_stack([np.repeat(n4[eligible, None], 7, axis=1), np.repeat(n8[eligible, None], 6, axis=1)])
            for i, j in product(range(7), range(6)):
                v4, v8 = len(np.unique(y[:, i])), len(np.unique(y[:, 7+j]))
                opportunities.append(dict(phase=phase, CD4_author_state=CD4[i], CD8_author_state=CD8[j], eligible_author_patients=int(eligible.sum()), CD4_positive_patients=int((x[eligible, i] > 0).sum()), CD8_positive_patients=int((x[eligible, 7+j] > 0).sum()), CD4_unique_fractions=v4, CD8_unique_fractions=v8, status='REFERENCE_COUNT_SCOPE_READY_NOT_POPULATION_OR_CLINICAL_TRUTH' if eligible.sum() >= 20 and v4 > 1 and v8 > 1 else 'INSUFFICIENT_REFERENCE_COUNT_SCOPE'))
    ref_counts = unique.groupby(['library_id', 'node_id']).size()
    overlap_counts = cur[cur.position.isin(shared)].groupby(['library_id', 'node_id']).size()
    node_counts = cur.groupby(['library_id', 'node_id']).size().rename('all_cells').reset_index()
    node_counts['unique_reference_cells'] = [int(ref_counts.get((l, n), 0)) for l, n in zip(node_counts.library_id, node_counts.node_id)]
    node_counts['overlap_cells'] = [int(overlap_counts.get((l, n), 0)) for l, n in zip(node_counts.library_id, node_counts.node_id)]
    node_counts['unannotated_cells'] = node_counts.all_cells - node_counts.unique_reference_cells - node_counts.overlap_cells
    assert (node_counts.unannotated_cells >= 0).all()
    node_counts['phase'] = node_counts.library_id.map(metadata.timepoint_registry)
    node_counts['arm'] = node_counts.library_id.map(metadata.neo_type)
    for (phase, arm, node), g in node_counts.groupby(['phase', 'arm', 'node_id']):
        g = g[g.all_cells >= 20]
        fractions = g.unique_reference_cells.to_numpy(float) / g.all_cells.to_numpy(float)
        cover.append(dict(phase=phase, arm=arm, current_node=node, libraries_node_atleast20=len(g), unique_reference_fraction_min_q25_median_q75_max=np.quantile(fractions, [0, .25, .5, .75, 1]).tolist() if len(g) else None, libraries_reference_below_half=int((fractions < .5).sum()), all_current_cells=int(g.all_cells.sum()), unique_reference_cells=int(g.unique_reference_cells.sum()), overlap_cells=int(g.overlap_cells.sum()), unannotated_cells=int(g.unannotated_cells.sum())))
    for scope, counts in matrices.items():
        for phase in ['all', 'pre', 'post', 'post2']:
            g = tab if phase == 'all' else tab[tab.timepoint_registry.eq(phase)]
            for labels, original_label in [([l], l) for l in CD4 + CD8 if l in OLD] + [(list(s), 'Teff') for k in range(1, 7) for s in combinations(CD8, k)]:
                actual = counts.loc[g.library_id, labels].sum(axis=1).to_numpy(int)
                original = g[original_label].to_numpy(int)
                difference = actual - original
                versions.append(dict(count_scope=scope, phase=phase, fine_labels=labels, original_count_label=original_label, reviewed_libraries=len(g), exact_count_libraries=int((difference == 0).sum()), exact_nonzero_libraries=int(((difference == 0) & ((actual > 0) | (original > 0))).sum()), all_libraries_exact=bool((difference == 0).all()), original_total=int(original.sum()), fine_total=int(actual.sum()), maximum_absolute_difference=int(np.abs(difference).max()) if len(g) else None, sum_absolute_difference=int(np.abs(difference).sum()), no_automatic_label_equivalence=True))
    assert len(opportunities) == 84
    save('private_all3_phase_reference_scope.json', phase_rows)
    save('private_all84_fine_state_opportunities.json', opportunities)
    save('private_all_node_phase_arm_patient_coverage.json', cover)
    save('private_all_count_version_comparisons.json', versions)
    result = dict(status='COMPLETE_ALL_FINE_REFERENCE_PATIENT_COVERAGE_AND_COUNT_VERSION_OPPORTUNITIES', completed_UTC=datetime.now(timezone.utc).isoformat(), source_SHA256=sha(__file__), original_author_join_SHA256=sha(table), all84_opportunities_preserved=True, both_reference_joins_match_accepted_source_crosscounts=True, all_node_library_denominators_partitioned=True, no_cell_or_patient_rows_saved=True, no_new_clinical_expression_TCR_or_pair_effects=True, files={p.name: dict(SHA256=sha(p), bytes=p.stat().st_size) for p in out.iterdir()})
    save('private_terminal.json', result)
    print(json.dumps(dict(status=result['status'], terminal_SHA256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
