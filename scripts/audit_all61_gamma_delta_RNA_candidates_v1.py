"""UNEXECUTED hard-rule extension audit; science gate currently not accepted.

Complete candidate/metadata/count audit plus bounded independent CSR rows.
"""
import json
import os
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
import zarr
from gamma_delta_RNA_rules_v1 import MARKERS, RULES, digest_array
from NSCLC_external_lineage_inputs_v2 import load, cell_arms
from external_functional_coverage_inputs_v1 import read_inputs
from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from NSCLC_immune_RNA_MCP_cohort_io_v2 import wait


def independent_bits(x, known):
    result = np.zeros(len(x), np.uint8)
    for i, row in enumerate(x):
        for j, minimum in enumerate([1, 2]):
            if all(known[:5]) and all(int(v) >= minimum for v in row[:5]): result[i] |= 1 << (2 * j)
            if all(known[:7]) and sum(int(v) >= minimum for v in row[:4]) >= 3 and row[4] >= minimum and max(row[5:7]) >= minimum: result[i] |= 1 << (2 * j + 1)
    return result


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']); guard = resource(root, 2 * 1024**2)
    src = root / 'all61_gamma_delta_RNA_candidates_20261019_v1'
    out = root / 'all61_gamma_delta_RNA_candidate_audit_20261019_v1'; assert not out.exists(); out.mkdir(mode=0o700)
    t = wait(root, out, 'measure_all61_gamma_delta_RNA_candidates_v1', src / 'private_terminal.json', {'COMPLETE_ALL61_FIXED_GAMMA_DELTA_RNA_CANDIDATES_BOUNDED_RAW_AUDIT_PENDING'})
    for name, d in t['files'].items(): assert sha(src / name) == d['sha256']
    plan = json.loads((src / 'plan_before_all61_native_candidate_RNA.json').read_text())
    assert plan['markers'] == MARKERS and plan['rules'] == RULES
    data = load(); axes = read_inputs(root); lib = data['lib']
    assert data['hashes'] == plan['input_metadata_hashes'] and axes['input_sha256'] == plan['native_axis_input_hashes']
    scope = pd.read_csv(src / 'private_complete61_native_measurement_scope.tsv', sep='\t', keep_default_na=False).set_index(['arm', 'dataset'])
    proof = {(v['arm'], v['dataset']): v for v in json.loads((src / 'private_all_source_count_and_candidate_proofs.json').read_text())}
    selected = pd.read_csv(src / 'private_prespecified_source_raw_audit_rows.tsv.gz', sep='\t', keep_default_na=False)
    totals = np.zeros((len(lib), 4, 4), np.int64); measurable = np.zeros((len(lib), 4), bool)
    raw_cells = 0; source_entries = 0
    for chunk in cell_arms(data):
        arm = chunk['arm']; original = arm == 'frozen_original_annotation'
        z = zarr.open_group(str(data['old_zarr'] if original else data['atlas']), mode='r'); raw = z['X'] if original else z['layers/counts']
        if original:
            with np.load(data['old'] / 'private_canonical_row_selection.npz', allow_pickle=False) as d: rows = d['global_rows']
        else: rows = np.arange(chunk['cell_count'])
        inverse = np.full(int(raw.attrs['shape'][0]), -1, np.int32); inverse[rows] = np.arange(len(rows))
        path = src / ('private_' + arm + '_candidate_identity_sidecar.npz')
        record = dict(np.load(path, allow_pickle=False)) if path.exists() else {k: np.array([], int) for k in ['global_row', 'library_position', 'frozen_node_code', 'parent_group', 'rule_bits', 'TRDV1_detection', 'TRDV2_detection']}
        gr = record['global_row']; local = inverse[gr]
        assert len(gr) == len(np.unique(gr)) and (local >= 0).all() and (record['rule_bits'] > 0).all() and (record['rule_bits'] < 16).all()
        assert np.all(((record['rule_bits'] & 4) == 0) | ((record['rule_bits'] & 1) > 0))
        assert np.all(((record['rule_bits'] & 8) == 0) | ((record['rule_bits'] & 2) > 0))
        np.testing.assert_array_equal(record['library_position'], chunk['library_index'][local])
        np.testing.assert_array_equal(record['frozen_node_code'], chunk['node'].codes[local])
        assert json.loads((src / ('private_' + arm + '_frozen_node_categories.json')).read_text()) == [str(v) for v in chunk['node'].categories]
        majors = chunk['major'][chunk['node'].codes[local]]
        parents = np.array([0 if m in {'T', 'T cell'} else 1 if m in {'NK', 'T_NK_boundary', 'T/NK lymphoid', 'Lymphoid cell'} else 2 for m in majors])
        np.testing.assert_array_equal(record['parent_group'], parents)
        for j in range(4):
            choose = (record['rule_bits'] & (1 << j)) > 0
            totals[:, j, 0] += np.bincount(record['library_position'][choose], minlength=len(lib))
            for p in range(3): totals[:, j, p + 1] += np.bincount(record['library_position'][choose & (parents == p)], minlength=len(lib))
        for dataset in sorted(str(v) for v in chunk['source'].categories):
            source_entries += 1; s = scope.loc[(arm, dataset)]
            li = lib.index[lib.arm.eq(arm) & lib.dataset.eq(dataset)].to_numpy(int)
            axis = axes['source_axes'][(arm, dataset)]
            if original:
                g = pd.read_csv(root.parents[2] / 'Project_v3/annotation_rebuild/20261010_uniform_expression_annotation_v1/by_gse' / dataset / 'genes.tsv', sep='\t', keep_default_na=False).set_index('gene')
                cols = [int(g.loc[m, 'core_column']) if m in g.index else -1 for m in MARKERS]
            else:
                cols = [next(iter(axis['lookup'][m])) if len(axis['lookup'].get(m, set())) == 1 else -1 for m in MARKERS]
            known = np.array([c >= 0 for c in cols]); assert s.physical_columns == str(cols)
            assert s.core5_measurable == all(known[:5]) and s.supported7_measurable == all(known[:7])
            gates = [all(known[:5]), all(known[:7]), all(known[:5]), all(known[:7])]; measurable[li] = gates
            assert s.libraries == len(li) and s.canonical_cells == int(lib.loc[li, 'annotated_scope_cells'].sum())
            members = np.isin(record['library_position'], li)
            if not gates[0]:
                assert not members.any() and not s.RNA_read and (arm, dataset) not in proof; continue
            p = proof[(arm, dataset)]; assert p['physical_columns'] == cols and p['measured_markers'] == known.tolist()
            rr = selected[selected.arm.eq(arm) & selected.dataset.eq(dataset)]
            sample = rr.global_row.to_numpy(np.int64); assert len(sample) == p['raw_audit_cells'] and np.array_equal(sample, np.unique(sample))
            assert digest_array(sample, '<i8') == p['raw_audit_rows_sha256']
            assert (inverse[sample] >= 0).all()
            np.testing.assert_array_equal(rr.library_position, chunk['library_index'][inverse[sample]])
            x = np.zeros((len(sample), 9), np.int32)
            for i, row in enumerate(sample):
                bounds = np.asarray(raw['indptr'][int(row):int(row) + 2], np.int64)
                index = np.asarray(raw['indices'][int(bounds[0]):int(bounds[1])]); values = np.asarray(raw['data'][int(bounds[0]):int(bounds[1])])
                for k, column in enumerate(cols):
                    if column >= 0:
                        v = values[index == column]
                        assert np.isfinite(v).all() and (v >= 0).all() and np.array_equal(v, np.floor(v))
                        x[i, k] = int(v.sum())
            assert digest_array(x, '<i4') == p['raw_audit_count_sha256']
            bits = independent_bits(x, known); assert digest_array(bits, 'u1') == p['raw_audit_rule_bits_sha256']
            found = pd.Series(record['rule_bits'], index=gr).reindex(sample, fill_value=0).to_numpy(np.uint8)
            np.testing.assert_array_equal(bits, found)
            assert digest_array(gr[members], '<i8') == p['all_candidate_global_rows_sha256'] and digest_array(record['rule_bits'][members], 'u1') == p['all_candidate_rule_bits_sha256']
            for k, key in [(7, 'TRDV1_detection'), (8, 'TRDV2_detection')]:
                if not known[k]: assert (record[key][members] == -1).all()
                else: assert np.isin(record[key][members], [0, 1]).all()
            raw_cells += len(sample)
            (out / 'private_progress.json').write_text(json.dumps(dict(stage='BOUNDED_INDEPENDENT_SOURCE_CSR_ROWS', checked_source_entries=source_entries, observed_UTC=datetime.now(timezone.utc).isoformat())) + '\n')
    assert source_entries == 61
    tab = pd.read_csv(src / 'private_all1690_library_four_rule_counts.tsv.gz', sep='\t', keep_default_na=True)
    assert len(tab) == 6760 and not tab.duplicated(['arm', 'dataset', 'library_id', 'rule']).any()
    columns = ['candidate_cells', 'frozen_T_parent_candidate_cells', 'frozen_boundary_or_NK_candidate_cells', 'frozen_other_parent_candidate_cells']
    for j, rule in enumerate(RULES):
        g = tab[tab.rule.eq(rule)].reset_index(drop=True)
        pd.testing.assert_frame_equal(g[['arm', 'dataset', 'library_id']], lib[['arm', 'dataset', 'library_id']], check_dtype=False)
        np.testing.assert_array_equal(g.RNA_rule_measurable, measurable[:, j])
        assert g.loc[~measurable[:, j], columns].isna().all().all()
        np.testing.assert_array_equal(g.loc[measurable[:, j], columns].to_numpy(np.int64), totals[measurable[:, j], j])
        np.testing.assert_array_equal(g.annotated_scope_cells, lib.annotated_scope_cells)
    srcsum = pd.read_csv(src / 'private_complete61_source_rule_summary.tsv', sep='\t').set_index(['arm', 'dataset', 'rule'])
    for key, g in tab.groupby(['arm', 'dataset', 'rule'], sort=True):
        r = srcsum.loc[key]; good = g.RNA_rule_measurable
        assert r.libraries == len(g) and r.measurable_libraries == good.sum() and r.libraries_candidate_atleast20 == g.loc[good, 'candidate_cells'].ge(20).sum()
        for column in columns:
            if good.any(): assert r[column] == g.loc[good, column].sum()
            else: assert pd.isna(r[column])
    receipt = dict(status='ACCEPTED_ALL61_CANDIDATE_KEYS_FOUR_RULE_COUNTS_MISSINGNESS_PARENT_GROUPS_AND_BOUNDED_SOURCE_RAW_ROWS', completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), producer_terminal_sha256=sha(src / 'private_terminal.json'), resource_guard=guard, source_arms_checked=61, library_rule_positions_checked=6760, independently_raw_cells=raw_cells, whole_source_independent_raw_replay=False, clinical_or_formal_identity_adoption=False, scientific_outputs_public=False)
    (out / 'private_terminal.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(dict(status=receipt['status'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
