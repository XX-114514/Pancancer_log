"""All original internal candidates: fixed, post-selection mechanistic diagnostics.

The ligand/reference scope is server-private. No new P, endpoint replacement,
clinical outcome, threshold rescue or independent-validation claim.
"""
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd

from NSCLC_native_MFA_interpretation_inputs_v1 import load, pool_native, sha, VIEWS
from test_NSCLC_native_interface_association_v1 import inputs, context
from prepare_NSCLC_external_fixed_gene_programme_assay_v1 import resource
from prepare_NSCLC_fixed_programme_continuous_priors_v1 import read_rds, PRIORS
from disjoint_hallmark_ULM_kernel_v1 import score
import NSCLC_native_interface_association_math_v1 as kernel

SCOPE_SHA = 'b87fa57fb8246056d36c619c0edaa86dcd0b0d3d830d670524c1d8c6d3144550'
OUT = 'NSCLC_native_interface_mechanistic_followup_20261018_v1'


def read_scope(data):
    root = data['root']
    path = root / 'NSCLC_native_interface_mechanistic_scope_20261018_v1.json'
    assert sha(path) == SCOPE_SHA
    scope = json.loads(path.read_text())
    src = root / 'NSCLC_native_interface_association_review_20261018_v1'
    assert sha(src / 'private_independent_full_review_audit.json') == scope['full_review_proof_sha256']
    path = src / 'private_all_original_internal_candidates.tsv.gz'
    assert sha(path) == scope['candidate_table_sha256']
    candidates = pd.read_csv(path, sep='\t', keep_default_na=False)
    assert candidates.record.tolist() == scope['records'] and candidates.record.is_unique
    assert set(candidates.ligand) == set(scope['ligand_specs'])
    return scope, candidates


def exclusions(root, scope):
    prior = root / 'immune_same_name_family_source_review_v1/fullgene_receiver_target_readiness_v1/verified_prior_inputs_v2'
    path = prior / 'lr_network_human_21122021.rds'
    assert sha(path) == PRIORS[path.name]
    flat = read_rds(path)
    path = root / 'LR_binding_resource_review_v1/all4624_reference_evidence_status.tsv'
    assert sha(path) == 'df574533d866438e646423f29bd555cb283173dbceadfe9f8264e1622a86fd06'
    refs = pd.read_csv(path, sep='\t', keep_default_na=False)
    result = {}
    for ligand, spec in scope['ligand_specs'].items():
        strings = {ligand, *flat.loc[flat['from'].eq(ligand), 'to'], *refs.loc[refs.ligand.eq(ligand), 'receptor']}
        parts = {str(g) for s in strings for g in str(s).split('_')} | set(strings)
        parts.update(g for complex_ in spec['complexes'] for g in complex_)
        result[ligand] = parts
    return result


def sender_without_ligand(data, vi, factor, ligand):
    p = data['prepared'][vi]
    w, tau = data['model']['weights_' + str(vi)], data['model']['precision_' + str(vi)]
    a = np.linalg.solve(np.eye(10) + w.T @ (tau[:, None] * w), w.T * tau)
    assert np.allclose(p['X'] @ a.T, data['model']['own_view_scores_' + str(vi)], atol=1e-10, rtol=1e-10)
    keep = data['genes'][p['selected_native_columns']] != ligand
    return p['X'][:, keep] @ a[factor, keep], int((~keep).sum())


def main():
    data = load(); root = data['root']; budget = resource(root, 16 * 1024**2)
    available = int(next(s for s in open('/proc/meminfo') if s.startswith('MemAvailable:')).split()[1]) * 1024
    assert available > 2 * 1024**3
    scope, candidates = read_scope(data)
    out = root / OUT; assert not out.exists(); out.mkdir(mode=0o700)
    plan = dict(created_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__),
                scope_sha256=SCOPE_SHA, budget=budget, RAM_available_bytes=available,
                all_original_candidate_records_retained=True,
                score='Binary-membership ULM and target/non-target means; remove union of all learned genes and all original/curated LR components from both sets; no RNA renormalization.',
                support='At least20 remaining target genes and500 non-target genes; median positive targets<10 or target mean<=0 flagged, no gene replacement.',
                effects='Original primary and discovery cases with original complete nuisance span. Same rank/SVD gates min50 df10; no P, intervals or extra selection.',
                controls='Base; recipient own same-ligand RNA; original sender same-factor programme; sender programme with ligand feature removed without refit; both original immune fractions.',
                missing_complex='Component not packed: true joint is missing, not zero; report marginal minimum upper bound separately.',
                post_selection_descriptive_not_independent_confirmation=True,
                new_P_BH_clinical_or_primary_changes=False)
    (out / 'plan_before_direct_functional_descriptions.json').write_text(json.dumps(plan, indent=2) + '\n')
    _, _, _, z, clinical, _ = inputs(root); state = context(data, z, clinical)
    native = root / 'NSCLC_continuous_prior_native_interfaces_20261016_v2'
    manifest = json.loads((native / 'private_terminal.json').read_text())
    for name in ['private_composite_panel_axis.tsv', 'private_positive_panel_packed.npy']:
        assert sha(native / name) == manifest['files'][name]['sha256']
    panel = pd.read_csv(native / 'private_composite_panel_axis.tsv', sep='\t').gene.tolist()
    packed = np.load(native / 'private_positive_panel_packed.npy', mmap_mode='r')
    cell_unit_path = data['native'] / 'private_cell_unit_index.npy'
    original = json.loads((data['native'] / 'private_terminal.json').read_text())
    for name in ['private_cell_unit_index.npy', 'positive_cell_counts.npy']:
        assert sha(data['native'] / name) == original['outputs'][name]['sha256']
    cell_units = np.load(cell_unit_path, mmap_mode='r')
    cell_pool = z['native_unit_to_case_pool'][cell_units]
    hits = np.load(data['native'] / 'positive_cell_counts.npy', mmap_mode='r')
    genes = pd.Index(data['genes']); vi_by = {v: i for i, v in enumerate(VIEWS)}
    pooled, _ = pool_native(data); omitted = exclusions(root, scope)
    rows, receptors, senders, controls, audit_inputs = [], [], [], [], {}
    for ri, item in enumerate(candidates.itertuples()):
        sv, rv = vi_by[item.sender_view], vi_by[item.receiver_view]
        factor = int(item.factor.rsplit('F', 1)[1]) - 1; li = int(item.ligand_axis_index)
        spec = scope['ligand_specs'][item.ligand]
        without, removed = sender_without_ligand(data, sv, factor, item.ligand)
        remain = data['remaining'][~np.isin(data['genes'][data['remaining']], list(omitted[item.ligand]))]
        for part in ['internal_holdout', 'discovery']:
            ix, quality, blocks = state['parts'][(item.sender_view, item.receiver_view, part)]
            x, y = z['ligand_mean'][sv, ix, li], item.direction * z['response'][item.slot_index, ix]
            own = z['ligand_mean'][rv, ix, li]
            alternatives = dict(base=quality, recipient_own_ligand=np.column_stack([quality, own]),
                sender_original_programme=np.column_stack([quality, z['original_scores'][sv, ix, factor]]),
                sender_ligand_removed_programme=np.column_stack([quality, without[ix]]),
                native_immune_fractions=np.column_stack([quality, state['fractions'][[sv, rv]][:, ix].T]))
            audit_inputs[f'{ri}_{part}_x'] = x; audit_inputs[f'{ri}_{part}_blocks'] = blocks
            audit_inputs[f'{ri}_{part}_programme'] = y
            for label, q in alternatives.items():
                audit_inputs[f'{ri}_{part}_q_{label}'] = q
                value, status = kernel.secondary(x, y, q, blocks)
                if part == 'internal_holdout' and label == 'base':
                    assert status == 'ESTIMABLE' and abs(value - item.partial_rank_rho) < 1e-10
                controls.append(dict(record=item.record, partition=part, control=label, rho=value, status=status, n_cases=len(ix), sender_ligand_features_removed=removed))
            # Unit marginal positivity and sender subtype contributions, no case rows exported.
            selected = np.isin(z['native_unit_to_case_pool'], sv * 208 + ix)
            for label, g in data['units'][selected].groupby('sub_cell_type'):
                positives = int(hits[g.unit_index.to_numpy(int), genes.get_loc(item.ligand)].sum())
                senders.append(dict(record=item.record, partition=part, author_subtype=label,
                                    cells=int(g.n_cells.sum()), ligand_positive_cells=positives))
            for complex_ in spec['complexes']:
                rec = '_'.join(complex_)
                unit_take = np.isin(z['native_unit_to_case_pool'], rv * 208 + ix)
                unit_rows = np.flatnonzero(unit_take)
                marginal = [int(hits[unit_rows, genes.get_loc(g)].sum()) for g in complex_]
                denom = int(z['cells'][rv, ix].sum())
                entry = dict(record=item.record, partition=part, receptor=rec,
                    receiver_cells=denom, marginal_counts_json=json.dumps(dict(zip(complex_, marginal))),
                    joint_upper_bound=min(marginal), joint_positive_cells=np.nan, positive_cases=np.nan,
                    joint_status='COMPONENT_NOT_IN_ACCEPTED_PACKED_PANEL')
                if set(complex_).issubset(panel):
                    joint = np.ones(len(cell_units), bool)
                    for gene in complex_:
                        j = panel.index(gene)
                        joint &= ((packed[:, j // 8] >> (j % 8)) & 1).astype(bool)
                    valid = joint & (cell_pool >= 0)
                    counts = np.bincount(cell_pool[valid], minlength=4*208).reshape(4, 208)[rv, ix]
                    assert counts.sum() <= min(marginal)
                    entry.update(joint_positive_cells=int(counts.sum()), positive_cases=int((counts > 0).sum()), joint_status='ACTUAL_SAME_CELL_RNA')
                    if rec in z['receptors']:
                        assert np.array_equal(counts, z['receptor_positive'][rv, ix, list(z['receptors']).index(rec)])
                receptors.append(entry)
        for function in spec['functions']:
            assert function in data['refs']
            targets = np.flatnonzero(np.isin(data['genes'][remain], list(data['refs'][function])))
            qualified = len(targets) >= 20 and len(remain) - len(targets) >= 500
            measurements = score(pooled[rv][:, remain], [targets]) if qualified else None
            for part in ['internal_holdout', 'discovery']:
                ix, quality, blocks = state['parts'][(item.sender_view, item.receiver_view, part)]
                for readout in ['ULM_score', 'target_mean', 'nontarget_mean']:
                    values = measurements[readout][:, 0] if qualified else np.full(208, np.nan)
                    low = bool(qualified and (np.median(measurements['expressed_target_genes'][ix, 0]) < 10 or np.median(measurements['target_mean'][ix, 0]) <= 0))
                    for label, q in [('base', quality),
                                     ('recipient_own_ligand', np.column_stack([quality, z['ligand_mean'][rv, ix, li]])),
                                     ('native_immune_fractions', np.column_stack([quality, state['fractions'][[sv, rv]][:, ix].T]))]:
                        value, status = kernel.secondary(z['ligand_mean'][sv, ix, li], values[ix], q, blocks)
                        rows.append(dict(record=item.record, partition=part, function=function, readout=readout, control=label,
                            remaining_targets=len(targets), remaining_background=len(remain)-len(targets), qualified=qualified,
                            low_foreground_support=low, n_cases=len(ix), rho=value, status=status))
    artifacts = {'private_all_functional_descriptions.tsv': pd.DataFrame(rows),
                 'private_all_canonical_receptor_support.tsv': pd.DataFrame(receptors),
                 'private_all_sender_subtype_contributions.tsv': pd.DataFrame(senders),
                 'private_all_original_programme_control_descriptions.tsv': pd.DataFrame(controls)}
    for name, table in artifacts.items(): table.to_csv(out / name, sep='\t', index=False)
    np.savez_compressed(out / 'private_control_design_replay.npz', **audit_inputs)
    files = [*artifacts, 'private_control_design_replay.npz']
    assert sum((out / name).stat().st_size for name in files) < 16*1024**2
    terminal = dict(status='COMPLETE_FIXED_POST_SELECTION_MECHANISTIC_DESCRIPTIONS_AUDIT_PENDING',
        completed_UTC=datetime.now(timezone.utc).isoformat(), code_sha256=sha(__file__),
        scope_sha256=SCOPE_SHA, files={name: dict(sha256=sha(out/name), bytes=(out/name).stat().st_size) for name in files},
        new_P_BH_or_primary_changes=False, independent_biological_validation=False)
    (out / 'private_terminal.json').write_text(json.dumps(terminal, indent=2) + '\n')
    print(json.dumps(dict(status=terminal['status'], completed_UTC=terminal['completed_UTC'], terminal_sha256=sha(out/'private_terminal.json'))))


if __name__ == '__main__':
    main()
