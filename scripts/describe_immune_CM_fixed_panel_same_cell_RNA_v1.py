"""Complete fixed-panel co-detection/source attribution; no new P or labels."""
import json
import os
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from NSCLC_native_MFA_interpretation_inputs_v1 import sha
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource


def main():
    root = Path(os.environ['PANCANCER_PREFLIGHT_ROOT'])
    budget = resource(root, 64 * 1024**2)
    source = root / 'immune_CM_fixed_panel_same_cell_RNA_20261018_v1'
    audit = root / 'immune_CM_fixed_panel_same_cell_RNA_audit_20261018_v1'
    proof = json.loads((audit / 'private_terminal.json').read_text())
    assert proof['status'].startswith('ACCEPTED_ALL_ORIGINAL_T_CELLS_')
    assert sha(source / 'private_terminal.json') == proof['producer_terminal_sha256']
    for folder in [source, audit]:
        for n, v in json.loads((folder / 'private_terminal.json').read_text())['files'].items():
            assert sha(folder / n) == v['sha256']
    fixed_path = root.parents[1] / 'runs/20261018_priority_immune_CM_bulk_v1/analysis_v2/frozen_plan_before_effects.json'
    fixed = json.loads(fixed_path.read_text());panels = fixed['panels'];names = list(panels)
    u = pd.read_csv(source / 'private_original_T_unit_scope.tsv', sep='\t', keep_default_na=False)
    z = np.load(source / 'private_fixed34_gene_T_counts_and_depth.npz', allow_pickle=False)
    raw = z['raw_gene_counts'];genes = z['genes'].tolist();unit = z['original_T_unit_index'];total = z['full_native_total_counts']
    nf = len(genes);nu = len(u);nc = len(unit)
    source_contract = json.loads((source / 'private_source_gene_measurement_contract.json').read_text())
    available = np.array([[next(s for s in source_contract if s['dataset'] == dataset)['complete_panels'][k] for k in names] for dataset in u.dataset])
    gene_available = np.array([[v >= 0 for v in next(s for s in source_contract if s['dataset'] == dataset)['physical_columns']] for dataset in u.dataset])
    assert nf == 34 and len(names) == 7 and (total > 0).all()
    pairs = list(combinations(range(len(names)), 2));assert len(pairs) == 21
    out = root / 'immune_CM_fixed_panel_same_cell_descriptions_20261018_v1'
    assert not out.exists();out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__), budget=budget,
        raw_input_terminal_sha256=sha(source / 'private_terminal.json'), independent_raw_audit_sha256=sha(audit / 'private_terminal.json'),
        fixed_panels_sha256=sha(fixed_path), scopes=['all_native_T', 'original_CM_scope_libraries', 'original_CM_scope_pre_libraries'],
        gates=['original_T_annotation', 'pan_T_RNA_ge2'], thresholds=[1, 2],
        panels=panels, pairs=[[names[a], names[b]] for a, b in pairs],
        support_definitions='Regulatory support=FOXP3 captured plus>=1 other ORIGINAL regulatory gene. Effector support=>=2 of ORIGINAL GZMB/PRF1/CTSW/NKG7. CD8A/CD8B alone cannot prove cytotoxic activity. These are descriptive decompositions, not changed clinical signatures or cell labels.',
        primary_descriptions='Every unit/panel/pair; source and all16-node counts and all corresponding scope memberships. Capture-positive requires>=2 distinct genes; panT gate uses threshold1 in both sensitivity views.',
        depth_benchmark='Expected overlap from product of panel detection fractions within original library-node x full-RNA deciles, summed to units. Decile boundaries preserve ties. No P, CI, patient independence or inference that this removes confounding.',
        contribution='All34 native raw transcript sums and full-native CP10K sums in four original R/C capture patterns and four FOXP3/effector-supported patterns. Gene-specific within-source contribution shares, not bulk percentile/min score decomposition or learned weights.',
        limitations='Source capture, dropout, original T annotation, ex-vivo/material/timepoint, overlapping signatures and doublets remain. FOXP3-not-detected is not absence of Treg. Same-cell RNA is not functional regulation, secretion, direct communication or clinical prediction. Missing panels stay unavailable; empty groups retained.',
        no_clinical_refit_or_new_statistical_test=True)
    (out / 'plan_before_all_same_cell_descriptions.json').write_text(json.dumps(plan, indent=2) + '\n')
    order = np.argsort(unit, kind='stable');ends = np.cumsum(u.n_cells.to_numpy());starts = np.r_[0, ends[:-1]]
    depth = np.zeros(nc, np.int8)
    for i, (start, stop) in enumerate(zip(starts, ends)):
        ix = order[start:stop];assert (unit[ix] == i).all()
        cuts = np.unique(np.quantile(total[ix], np.arange(1, 10) / 10))
        depth[ix] = np.searchsorted(cuts, total[ix], side='left')
    depth_group = unit * 10 + depth
    detection = {threshold: np.column_stack([(raw[:, [genes.index(g) for g in panels[k]]] >= threshold).sum(axis=1) for k in names]).astype(np.int8) for threshold in [1, 2]}
    pan = names.index('pan_T');reg = names.index('regulatory_associated');cyto = names.index('cytotoxic')
    identity_gates = [np.ones(nc, bool), (detection[1][:, pan] >= 2) & available[unit, pan]]
    scopes = [np.ones(nu, bool), u.original_CM_scope_library.to_numpy(bool), u.original_CM_scope_pre_library.to_numpy(bool)]
    source_names = sorted(u.dataset.unique());node_names = sorted(u.node_id.unique())
    summary = [];attribution = [];special_summary = []
    saved = dict(depth_decile=depth, unit_panel_available=available, unit_gene_available=gene_available,
        panels=np.asarray(names), pairs=np.asarray(pairs), original_unit_n_cells=u.n_cells.to_numpy())
    raw_nonnegative = np.maximum(raw, 0)
    cp = raw_nonnegative.astype(np.float64) * (10000 / total)[:, None]
    for threshold in [1, 2]:
        positive = (detection[threshold] >= 2) & available[unit]
        fox = positive[:, reg] & (raw[:, genes.index('FOXP3')] >= threshold)
        effector = ((raw[:, [genes.index(g) for g in ['GZMB', 'PRF1', 'CTSW', 'NKG7']]] >= threshold).sum(axis=1) >= 2) & available[unit, cyto]
        cd8_only = ((raw[:, [genes.index(g) for g in ['CD8A', 'CD8B']]] >= threshold).all(axis=1) &
                    (raw[:, [genes.index(g) for g in ['GZMB', 'PRF1', 'CTSW', 'NKG7']]] < threshold).all(axis=1) & available[unit, cyto])
        special = np.column_stack([fox, effector, cd8_only, fox & effector])
        for gate_number, gate in enumerate(identity_gates):
            gate_name = plan['gates'][gate_number]
            n = np.bincount(unit[gate], minlength=nu)
            pc = np.column_stack([np.bincount(unit[gate & positive[:, k]], minlength=nu) for k in range(len(names))])
            any_counts = np.column_stack([np.bincount(unit[gate & (detection[threshold][:, k] > 0) & available[unit, k]], minlength=nu) for k in range(len(names))])
            joint = np.column_stack([np.bincount(unit[gate & positive[:, a] & positive[:, b]], minlength=nu) for a, b in pairs])
            ng = np.bincount(depth_group[gate], minlength=nu * 10).astype(float)
            pg = np.column_stack([np.bincount(depth_group[gate & positive[:, k]], minlength=nu * 10) for k in range(len(names))])
            expected = np.column_stack([np.divide(pg[:, a] * pg[:, b], ng, out=np.zeros_like(ng), where=ng > 0).reshape(nu, 10).sum(axis=1) for a, b in pairs])
            sc = np.column_stack([np.bincount(unit[gate & special[:, k]], minlength=nu) for k in range(special.shape[1])])
            prefix = 'threshold%d_gate%d_' % (threshold, gate_number)
            saved.update({prefix + k: v for k, v in dict(n=n, panel_ge2=pc, panel_any=any_counts, pair_joint=joint, pair_expected_RNA_depth=expected, supported_capture=sc).items()})
            identity_available = np.ones(nu, bool) if gate_number == 0 else available[:, pan]
            for scope_number, scope_mask in enumerate(scopes):
                for source_name in source_names:
                    source_mask = scope_mask & u.dataset.eq(source_name).to_numpy()
                    for node in ['ALL_ORIGINAL_T_NODES', *node_names]:
                        select = source_mask if node == 'ALL_ORIGINAL_T_NODES' else source_mask & u.node_id.eq(node).to_numpy()
                        base = dict(dataset=source_name, scope=plan['scopes'][scope_number], original_node=node,
                                    threshold=threshold, identity_gate=gate_name, original_units=int(select.sum()),
                                    original_cells=int(u.loc[select, 'n_cells'].sum()), captured_cells=int(n[select].sum()),
                                    measurable_identity_units=int((select & identity_available).sum()), units_ge20=int((select & (n >= 20)).sum()))
                        for k, name in enumerate(names):
                            valid = select & identity_available & available[:, k]
                            summary.append(dict(**base, kind='panel', panelA=name, panelB='', shared_genes='',
                                available_units=int(valid.sum()), valid_denominator=int(n[valid].sum()),
                                positive_cells=int(pc[valid, k].sum()) if valid.any() else None,
                                any_gene_cells=int(any_counts[valid, k].sum()) if valid.any() else None,
                                joint_RNA_depth_expected=None, ratio_observed_to_depth_expected=None))
                        for j, (a, b) in enumerate(pairs):
                            valid = select & identity_available & available[:, a] & available[:, b]
                            exp = float(expected[valid, j].sum())
                            summary.append(dict(**base, kind='panel_pair', panelA=names[a], panelB=names[b],
                                shared_genes=';'.join(sorted(set(panels[names[a]]) & set(panels[names[b]]))),
                                available_units=int(valid.sum()), valid_denominator=int(n[valid].sum()),
                                positive_cells=int(joint[valid, j].sum()) if valid.any() else None, any_gene_cells=None,
                                joint_RNA_depth_expected=exp if valid.any() else None,
                                ratio_observed_to_depth_expected=float(joint[valid, j].sum() / exp) if exp > 0 else None,
                                one_panel_defines_identity_gate=bool(gate_number == 1 and pan in [a, b])))
                        valid = select & identity_available & available[:, reg] & available[:, cyto]
                        special_summary.append(dict(**base, available_both_panel_units=int(valid.sum()), valid_complete_both_panels_cells=int(n[valid].sum()),
                            regulatory_ge2=int(pc[valid, reg].sum()), cytotoxic_ge2=int(pc[valid, cyto].sum()),
                            FOXP3_supported=int(sc[valid, 0].sum()), effector_supported=int(sc[valid, 1].sum()),
                            CD8_pair_without_any_effector_gene=int(sc[valid, 2].sum()), FOXP3_and_effector_supported=int(sc[valid, 3].sum())))
            # Aggregate signal attribution to units in RAM, persist source summaries only.
            for axis, state in [('original_panels', positive[:, reg].astype(np.int8) + 2 * positive[:, cyto].astype(np.int8)),
                                ('FOXP3_effector_supported', fox.astype(np.int8) + 2 * effector.astype(np.int8))]:
                eligible = gate & available[unit, reg] & available[unit, cyto]
                columns = np.flatnonzero(eligible)
                grouping = csr_matrix((np.ones(len(columns), np.int64), (unit[columns] * 4 + state[columns], columns)), shape=(nu * 4, nc))
                raw_mass = np.asarray(grouping @ raw_nonnegative).reshape(nu, 4, nf)
                cp_mass = np.asarray(grouping @ cp).reshape(nu, 4, nf)
                state_n = np.asarray(grouping.sum(axis=1)).ravel().reshape(nu, 4)
                for scope_number, scope_mask in enumerate(scopes):
                    for source_name in source_names:
                        select = scope_mask & u.dataset.eq(source_name).to_numpy() & identity_available & available[:, reg] & available[:, cyto]
                        nstate = state_n[select].sum(axis=0)
                        raw_sum = raw_mass[select].sum(axis=0);cp_sum = cp_mass[select].sum(axis=0)
                        for g, gene in enumerate(genes):
                            measured_units = int((select & gene_available[:, g]).sum())
                            for pattern in range(4):
                                attribution.append(dict(dataset=source_name, scope=plan['scopes'][scope_number], threshold=threshold,
                                    identity_gate=gate_name, axis=axis, capture_pattern=pattern, gene=gene,
                                    measured_units=measured_units, cells_in_pattern=int(nstate[pattern]),
                                    raw_transcripts=int(raw_sum[pattern, g]) if measured_units else None,
                                    full_native_CP10K_sum=float(cp_sum[pattern, g]) if measured_units else None,
                                    raw_gene_RNA_share=float(raw_sum[pattern, g] / raw_sum[:, g].sum()) if measured_units and raw_sum[:, g].sum() > 0 else None,
                                    CP10K_gene_RNA_share=float(cp_sum[pattern, g] / cp_sum[:, g].sum()) if measured_units and cp_sum[:, g].sum() > 0 else None))
    np.savez_compressed(out / 'private_all_unit_fixed_panel_descriptions.npz', **saved)
    pd.DataFrame(summary).to_csv(out / 'private_all_scope_source_node_panel_and_pair_descriptions.tsv', sep='\t', index=False)
    pd.DataFrame(special_summary).to_csv(out / 'private_all_scope_source_node_supported_capture.tsv', sep='\t', index=False)
    pd.DataFrame(attribution).to_csv(out / 'private_all_scope_source_gene_capture_attribution.tsv', sep='\t', index=False)
    files = [p for p in out.iterdir() if p.is_file()]
    assert sum(p.stat().st_size for p in files) < 64 * 1024**2
    terminal = dict(status='COMPLETE_ALL_FIXED_PANEL_SOURCE_DESCRIPTIONS_DERIVED_AUDIT_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(), source_sha256=sha(__file__),
        files={p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in files},
        no_new_P_or_clinical_refit_or_formal_cell_states=True)
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
    print(json.dumps(dict(status=terminal['status'], completed_UTC=terminal['completed_UTC'], terminal_sha256=sha(out / 'private_terminal.json'))))


if __name__ == '__main__':
    main()
