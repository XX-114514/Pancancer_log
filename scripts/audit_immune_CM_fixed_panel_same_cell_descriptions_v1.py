"""Independent full derived count, overlap benchmark and RNA attribution audit."""
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def numeric(value, expected, tolerance=1e-9):
    if expected is None:
        assert pd.isna(value)
    elif isinstance(expected, (int, np.integer)):
        assert float(value) == expected
    else:
        assert np.isfinite(float(value)) and np.isclose(float(value), expected, rtol=tolerance, atol=1e-7)


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT']);budget = resource(root, 8 * 1024**2)
    source = root / 'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    derived = root / 'immune_CM_fixed_panel_same_cell_descriptions_20261018_v1'
    proof_path = root / 'immune_CM_fixed_panel_same_cell_RNA_audit_20261018_v1/private_terminal.json'
    proof = json.loads(proof_path.read_text())
    assert proof['status'].startswith('ACCEPTED_ALL_ORIGINAL_T_CELLS_')
    assert sha(source / 'private_terminal.json') == proof['producer_terminal_sha256']
    for directory in [source, derived]:
        for n, v in json.loads((directory / 'private_terminal.json').read_text())['files'].items():
            assert sha(directory / n) == v['sha256']
    plan = json.loads((derived / 'plan_before_all_same_cell_descriptions.json').read_text())
    assert sha(proof_path) == plan['independent_raw_audit_sha256']
    out = root / 'immune_CM_fixed_panel_same_cell_descriptions_audit_20261018_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    (out / 'plan_before_complete_derived_numeric_audit.json').write_text(json.dumps(dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), budget=budget,
        derived_terminal_sha256=sha(derived / 'private_terminal.json'),
        scope='All unit/panel/pair/support arrays, every source/scope/node table row, all gene attribution rows. CSR grouping instead of bincount; original-cell direct sums instead of unit RNA sums. No new P/clinical tests.'), indent=2) + '\n')
    u = pd.read_csv(source / 'private_original_T_unit_scope.tsv', sep='\t', keep_default_na=False)
    z = np.load(source / 'private_fixed34_gene_T_counts_and_depth.npz', allow_pickle=False)
    saved = np.load(derived / 'private_all_unit_fixed_panel_descriptions.npz', allow_pickle=False)
    raw = z['raw_gene_counts'];unit = z['original_T_unit_index'];total = z['full_native_total_counts'];genes = z['genes'].tolist()
    names = list(plan['panels']);panels = plan['panels'];pairs = [tuple(p) for p in saved['pairs']]
    assert list(saved['panels']) == names and len(pairs) == 21
    nu, nc = len(u), len(unit)
    contract = json.loads((source / 'private_source_gene_measurement_contract.json').read_text())
    by_source = {x['dataset']: x for x in contract}
    availability = np.array([[by_source[s]['complete_panels'][name] for name in names] for s in u.dataset])
    ga = np.array([[i >= 0 for i in by_source[s]['physical_columns']] for s in u.dataset])
    assert np.array_equal(availability, saved['unit_panel_available']) and np.array_equal(ga, saved['unit_gene_available'])
    depth = np.zeros(nc, np.int8)
    groups = pd.Series(np.arange(nc)).groupby(unit, sort=True).groups
    for i, ix in groups.items():
        rows = np.asarray(ix);cuts = np.unique(np.percentile(total[rows], list(range(10, 100, 10))))
        depth[rows] = (total[rows, None] > cuts[None, :]).sum(axis=1)
    assert np.array_equal(depth, saved['depth_decile'])
    detect = {th: np.column_stack([np.count_nonzero(raw[:, [genes.index(g) for g in panels[name]]] >= th, axis=1) for name in names]) for th in [1, 2]}
    pan, reg, cyt = [names.index(k) for k in ['pan_T', 'regulatory_associated', 'cytotoxic']]
    gates = [np.ones(nc, bool), (detect[1][:, pan] >= 2) & availability[unit, pan]]
    scope_masks = [np.ones(nu, bool), u.original_CM_scope_library.to_numpy(bool), u.original_CM_scope_pre_library.to_numpy(bool)]
    selectors = {}
    for sk, scope in zip(plan['scopes'], scope_masks):
        for source_name in sorted(u.dataset.unique()):
            for node in ['ALL_ORIGINAL_T_NODES', *sorted(u.node_id.unique())]:
                select = scope & u.dataset.eq(source_name).to_numpy()
                if node != 'ALL_ORIGINAL_T_NODES': select &= u.node_id.eq(node).to_numpy()
                selectors[sk, source_name, node] = select
    numeric_arrays = {};positive_states = {};supported = {}
    for th in [1, 2]:
        pos = (detect[th] >= 2) & availability[unit]
        fox = pos[:, reg] & (raw[:, genes.index('FOXP3')] >= th)
        eff = (np.count_nonzero(raw[:, [genes.index(g) for g in ['GZMB', 'PRF1', 'CTSW', 'NKG7']]] >= th, axis=1) >= 2) & availability[unit, cyt]
        cd8 = (raw[:, genes.index('CD8A')] >= th) & (raw[:, genes.index('CD8B')] >= th)
        eff_zero = np.count_nonzero(raw[:, [genes.index(g) for g in ['GZMB', 'PRF1', 'CTSW', 'NKG7']]] >= th, axis=1) == 0
        special = np.column_stack([fox, eff, cd8 & eff_zero & availability[unit, cyt], fox & eff])
        positive_states[th] = pos;supported[th] = (fox, eff)
        for gate_i, gate in enumerate(gates):
            rows = np.flatnonzero(gate)
            group = csr_matrix((np.ones(len(rows), np.int64), (unit[rows], rows)), shape=(nu, nc))
            gd = csr_matrix((np.ones(len(rows), np.int64), (unit[rows] * 10 + depth[rows], rows)), shape=(nu * 10, nc))
            n = np.asarray(group.sum(axis=1)).ravel();pc = np.asarray(group @ pos.astype(np.int8))
            anyc = np.asarray(group @ (((detect[th] > 0) & availability[unit]).astype(np.int8)))
            joint = np.asarray(group @ np.column_stack([pos[:, a] & pos[:, b] for a, b in pairs]).astype(np.int8))
            ndepth = np.asarray(gd.sum(axis=1)).ravel();cdepth = np.asarray(gd @ pos.astype(np.int8))
            exp = np.column_stack([np.divide(cdepth[:, a] * cdepth[:, b], ndepth, out=np.zeros(nu * 10), where=ndepth > 0).reshape(nu, 10).sum(axis=1) for a, b in pairs])
            sc = np.asarray(group @ special.astype(np.int8))
            data = dict(n=n, panel_ge2=pc, panel_any=anyc, pair_joint=joint, pair_expected_RNA_depth=exp, supported_capture=sc)
            for k, value in data.items():
                observed = saved['threshold%d_gate%d_%s' % (th, gate_i, k)]
                assert np.allclose(value, observed, rtol=1e-12, atol=1e-10)
            numeric_arrays[th, plan['gates'][gate_i]] = data
    table = pd.read_csv(derived / 'private_all_scope_source_node_panel_and_pair_descriptions.tsv', sep='\t')
    support = pd.read_csv(derived / 'private_all_scope_source_node_supported_capture.tsv', sep='\t')
    expected_rows = len(selectors) * 2 * 2 * (7 + 21)
    assert len(table) == expected_rows and len(support) == len(selectors) * 2 * 2
    assert not table.duplicated(['dataset', 'scope', 'original_node', 'threshold', 'identity_gate', 'kind', 'panelA', 'panelB']).any()
    pair_lookup = {tuple(p): j for j, p in enumerate(pairs)}
    def base(row):
        select = selectors[row.scope, row.dataset, row.original_node]
        values = numeric_arrays[row.threshold, row.identity_gate];n = values['n']
        ia = np.ones(nu, bool) if row.identity_gate == plan['gates'][0] else availability[:, pan]
        for k, expected in dict(original_units=int(select.sum()), original_cells=int(u.loc[select, 'n_cells'].sum()),
            captured_cells=int(n[select].sum()), measurable_identity_units=int((select & ia).sum()), units_ge20=int((select & (n >= 20)).sum())).items(): numeric(getattr(row, k), expected)
        return select, ia, values
    for row in table.itertuples():
        select, ia, v = base(row);a = names.index(row.panelA)
        if row.kind == 'panel':
            valid = select & ia & availability[:, a]
            numeric(row.positive_cells, int(v['panel_ge2'][valid, a].sum()) if valid.any() else None)
            numeric(row.any_gene_cells, int(v['panel_any'][valid, a].sum()) if valid.any() else None)
        else:
            b = names.index(row.panelB);j = pair_lookup[a, b];valid = select & ia & availability[:, a] & availability[:, b]
            exp = float(v['pair_expected_RNA_depth'][valid, j].sum())
            numeric(row.positive_cells, int(v['pair_joint'][valid, j].sum()) if valid.any() else None)
            numeric(row.joint_RNA_depth_expected, exp if valid.any() else None)
            numeric(row.ratio_observed_to_depth_expected, float(v['pair_joint'][valid, j].sum() / exp) if exp > 0 else None)
            shared = ';'.join(sorted(set(panels[names[a]]) & set(panels[names[b]])))
            assert ('' if pd.isna(row.shared_genes) else row.shared_genes) == shared
            assert (str(row.one_panel_defines_identity_gate).lower() == 'true') == (row.identity_gate == plan['gates'][1] and pan in [a, b])
        numeric(row.available_units, int(valid.sum()));numeric(row.valid_denominator, int(v['n'][valid].sum()))
    for row in support.itertuples():
        select, ia, v = base(row);valid = select & ia & availability[:, reg] & availability[:, cyt]
        for k, expected in dict(available_both_panel_units=int(valid.sum()), valid_complete_both_panels_cells=int(v['n'][valid].sum()),
            regulatory_ge2=int(v['panel_ge2'][valid, reg].sum()), cytotoxic_ge2=int(v['panel_ge2'][valid, cyt].sum()),
            FOXP3_supported=int(v['supported_capture'][valid, 0].sum()), effector_supported=int(v['supported_capture'][valid, 1].sum()),
            CD8_pair_without_any_effector_gene=int(v['supported_capture'][valid, 2].sum()), FOXP3_and_effector_supported=int(v['supported_capture'][valid, 3].sum())).items(): numeric(getattr(row, k), expected)
    attribution = pd.read_csv(derived / 'private_all_scope_source_gene_capture_attribution.tsv', sep='\t')
    assert len(attribution) == 3 * u.dataset.nunique() * 2 * 2 * 2 * 4 * 34
    group_keys = ['dataset', 'scope', 'threshold', 'identity_gate', 'axis']
    for (source_name, scope_name, threshold, identity_gate, axis), frame in attribution.groupby(group_keys, sort=True):
        select = selectors[scope_name, source_name, 'ALL_ORIGINAL_T_NODES'].copy()
        ia = np.ones(nu, bool) if identity_gate == plan['gates'][0] else availability[:, pan]
        select &= ia & availability[:, reg] & availability[:, cyt]
        mask = select[unit] & gates[plan['gates'].index(identity_gate)]
        ix = np.flatnonzero(mask);pos = positive_states[threshold]
        left, right = (pos[:, reg], pos[:, cyt]) if axis == 'original_panels' else supported[threshold]
        state = left[ix].astype(np.int8) + 2 * right[ix].astype(np.int8)
        counts = np.maximum(raw[ix], 0).astype(np.int64)
        cp = counts * (10000 / total[ix])[:, None]
        rna = np.array([counts[state == k].sum(axis=0, dtype=np.int64) for k in range(4)])
        norm = np.array([cp[state == k].sum(axis=0) for k in range(4)])
        ns = np.array([(state == k).sum() for k in range(4)])
        for row in frame.itertuples():
            g = genes.index(row.gene);pattern = row.capture_pattern;available_units = int((select & ga[:, g]).sum())
            numeric(row.measured_units, available_units);numeric(row.cells_in_pattern, int(ns[pattern]))
            numeric(row.raw_transcripts, int(rna[pattern, g]) if available_units else None)
            numeric(row.full_native_CP10K_sum, float(norm[pattern, g]) if available_units else None)
            numeric(row.raw_gene_RNA_share, float(rna[pattern, g] / rna[:, g].sum()) if available_units and rna[:, g].sum() > 0 else None)
            numeric(row.CP10K_gene_RNA_share, float(norm[pattern, g] / norm[:, g].sum()) if available_units and norm[:, g].sum() > 0 else None)
    result = dict(status='ACCEPTED_ALL_UNIT_PANEL_PAIR_SCOPE_NODE_AND_GENE_ATTRIBUTION_DESCRIPTIONS',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        derived_terminal_sha256=sha(derived / 'private_terminal.json'),
        full_original_raw_audit_sha256=sha(proof_path), table_rows=len(table), support_rows=len(support), gene_attribution_rows=len(attribution),
        all_unit_arrays_recomputed_by_independent_grouping=True, every_table_row_and_every_gene_attribution_direct_cell_sum_checked=True,
        raw_cache_reused_after_full_original_replay_not_new_biological_replicate=True, no_new_P_clinical_or_state_claim=True)
    (out / 'private_independent_acceptance.json').write_text(json.dumps(result, indent=2) + '\n')
    files = [p for p in out.iterdir() if p.is_file()]
    result['files'] = {p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in files}
    (out / 'private_terminal.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], completed_UTC=result['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
