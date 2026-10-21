"""Prespecified follow-up of two source-level T-state relationships.

Eight exploratory slots: two parents x binary whole/near/far and nuclear
residual whole. No outcomes, independent-patient or original CM02 projection.
Coordinates, source types and local positivity remain conditioned upon.
Only source aggregates leave RAM. Histone/DNA are not pure technical controls.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.util
import io
import numpy as np
import pandas as pd

_parent_path = Path(__file__).with_name('analyze_moldoveanu_fixed_Ki67_source_organization_v3.py')
assert hashlib.sha256(_parent_path.read_bytes()).hexdigest() == '1024838da17c105987e7b44908839a0945ca17e74cc28889cc8e3059c119ef9c'
_spec = importlib.util.spec_from_file_location('frozen_parent_context_core', _parent_path)
_parent = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(_parent)
load, audit_final_author_phenotypes = _parent.load, _parent.audit_final_author_phenotypes
covariance, means_by_group, shuffled, bh = _parent.covariance, _parent.means_by_group, _parent.shuffled, _parent.bh

SEED = 2026102201
B = 2000
MIN_PARENT = 20
MIN_EDGES = 20
MIN_MATERIALS = 5
BACKGROUND = ['Histone.H3', 'X190BCKG', 'X191Ir', 'X193Ir']
CONTEXTS = ['binary_all', 'binary_both_tumor_neighbor',
            'binary_both_not_tumor_neighbor', 'nuclear_residual_all']
PARENTS = [['Treg'], ['Th.ae', 'Th.naive'], ['Tc.ae', 'Tc.naive']]


def exact_edges(coords, neighbors):
    _, count = np.unique(coords, axis=0, return_counts=True)
    assert count.max() <= 50
    tree = neighbors.ExactEuclideanNeighbors(coords)
    distances, _ = tree.query(coords, k=min(11, len(coords)))
    pairs = set()
    for i, found in enumerate(tree.query_ball_point(coords, np.nextafter(distances[:, -1], np.inf))):
        for j in found:
            if i != j:
                pairs.add((min(i, int(j)), max(i, int(j))))
    return np.asarray(sorted(pairs), dtype=int).reshape(-1, 2)


def tumor_neighbors(edges, tumor):
    near = np.zeros(len(tumor), dtype=bool)
    a, b = edges.T
    near[a[tumor[b]]] = True
    near[b[tumor[a]]] = True
    return near


def nuclear_residual(part, rank):
    if len(part) < MIN_PARENT:
        return None, 'background_parent_below_minimum'
    if not set(BACKGROUND) <= set(part):
        return None, 'background_fields_missing'
    raw = part[['Ki67'] + BACKGROUND].to_numpy(dtype=float)
    if not np.isfinite(raw).all():
        return None, 'nonfinite_background_or_Ki67'
    scores, varied = rank(raw)
    if not bool(np.all(varied)):
        return None, 'constant_background_or_Ki67_no_covariate_dropped'
    design = np.column_stack([np.ones(len(part)), scores[:, 1:]])
    if len(part) <= design.shape[1] + 2 or np.linalg.matrix_rank(design) != design.shape[1]:
        return None, 'background_design_not_identifiable'
    if np.linalg.cond(design) > 1e8:
        return None, 'background_design_ill_conditioned'
    coefficient = np.linalg.lstsq(design, scores[:, 0], rcond=None)[0]
    residual = scores[:, 0] - design @ coefficient
    # Verify the least-squares normal equations by an independent scalar sum.
    for column in design.T:
        assert abs(sum(float(a)*float(b) for a, b in zip(column, residual))) <= 1e-7 * len(part)
    scale = float(np.sqrt(np.mean(residual**2)))
    if scale <= 1e-10:
        return None, 'no_residual_variation'
    return residual/scale, None


def pair_indices(edges, left, right, n):
    inverse_l = np.full(n, -1, dtype=int); inverse_r = inverse_l.copy()
    inverse_l[left] = np.arange(len(left)); inverse_r[right] = np.arange(len(right))
    a, b = edges.T
    forward = (inverse_l[a] >= 0) & (inverse_r[b] >= 0)
    reverse = (inverse_l[b] >= 0) & (inverse_r[a] >= 0)
    return (inverse_l[np.r_[a[forward], b[reverse]]],
            inverse_r[np.r_[b[forward], a[reverse]]])


def quadrant_contributions(x, y, a, b):
    left = x[a] - x[a].mean(); right = y[b] - y[b].mean()
    masks = [(left >= 0) & (right >= 0), (left < 0) & (right < 0),
             (left >= 0) & (right < 0), (left < 0) & (right >= 0)]
    pieces = [float(np.sum((left*right)[take])/len(a)) for take in masks]
    assert abs(sum(pieces) - covariance(x, y, a, b)) < 1e-9
    return dict(zip(['above_above', 'below_below', 'above_below', 'below_above'], pieces))


def aggregate(records, reasons, diagnostics):
    rows = []
    for slot in range(8):
        selected = records[slot]
        row = dict(parent_relation=['regulatory_CD8', 'conventional_CD4_CD8'][slot//4],
                   context=CONTEXTS[slot % 4], source_materials=len(selected),
                   missing_reason_counts=dict(reasons[slot]), independent_patients_reviewed=False,
                   P=None, BH=None, effect=None, patient_CI=None)
        if len(selected) < MIN_MATERIALS:
            row['status'] = 'INSUFFICIENT_MATERIALS_FOR_CONDITIONAL_TEST'
        else:
            effects = [r['observed']-r['expected'] for r in selected]
            effect = float(np.mean(effects))
            draws = np.mean(np.stack([r['draws']-r['expected'] for r in selected]), axis=0)
            for j in range(B):
                replay = sum(float(r['draws'][j]-r['expected']) for r in selected)/len(selected)
                assert abs(replay-draws[j]) < 1e-9
            row.update(status='EXPLORATORY_MATERIAL_CONDITIONAL_TEST', effect=effect,
                       P=(1+int(np.sum(np.abs(draws) >= abs(effect)-1e-12)))/(B+1),
                       positive_materials=sum(x > 0 for x in effects),
                       negative_materials=sum(x < 0 for x in effects),
                       leave_one_material_effect_range=[float((sum(effects)-max(effects))/(len(effects)-1)),
                                                        float((sum(effects)-min(effects))/(len(effects)-1))])
        if selected:
            row['material_equal_observed'] = float(np.mean([r['observed'] for r in selected]))
            row['material_equal_null_expected'] = float(np.mean([r['expected'] for r in selected]))
            row['source_pooled_cross_edges'] = sum(r['edges'] for r in selected)
        rows.append(row)
    for row, q in zip(rows, bh([r['P'] if r['P'] is not None else 1. for r in rows])):
        if row['P'] is not None:
            row['BH'] = float(q)
    return dict(family=rows, descriptive_continuous_quadrants=[
        dict(parent_relation=['regulatory_CD8', 'conventional_CD4_CD8'][i],
             source_materials=len(items), contribution_to_observed_covariance={
                 key: float(np.mean([x[key] for x in items])) for key in items[0]} if items else {})
        for i, items in enumerate(diagnostics)])


def statistics(frame, knn, rank, grid):
    records = [[] for _ in range(8)]
    reasons = [Counter() for _ in range(8)]
    diagnostics = [[], []]
    ties_at_threshold = int(frame['Ki67.prob'].eq(.5).sum())
    rng = np.random.default_rng(SEED)
    for _, part in frame.groupby('sample.id', sort=True):
        part = part.reset_index(drop=True)
        coords = part[['coord.x', 'coord.y']].to_numpy(dtype=float)
        edges = knn(coords)
        near = tumor_neighbors(edges, part['Cluster'].eq('melano').to_numpy())
        positions = [np.flatnonzero(part['Cluster'].isin(labels)) for labels in PARENTS]
        if len(positions[2]) < MIN_PARENT or max(len(positions[0]), len(positions[1])) < MIN_PARENT:
            for reason in reasons:
                reason['whole_parent_below_minimum'] += 1
            continue
        group = grid(coords, (coords.min(0), coords.max(0))) * 2 + near.astype(int)
        binary = []; residual = []; residual_reasons = []; rank_values = []
        for pos in positions:
            p = part.iloc[pos]
            prob = p['Ki67.prob'].to_numpy(dtype=float)
            binary.append((prob >= .5).astype(float) if np.isfinite(prob).all() else None)
            value, reason = nuclear_residual(p, rank)
            residual.append(value); residual_reasons.append(reason)
            ranked, _ = rank(p[['Ki67']].to_numpy()) if len(p) else (np.empty((0, 1)), None)
            rank_values.append(ranked[:, 0])
        for parent in range(2):
            a0, b0 = pair_indices(edges, positions[parent], positions[2], len(part))
            if len(positions[parent]) >= MIN_PARENT and len(a0) >= MIN_EDGES:
                diagnostics[parent].append(quadrant_contributions(rank_values[parent], rank_values[2], a0, b0))
            for context in range(4):
                slot = parent*4+context
                left_pos, right_pos = positions[parent], positions[2]
                if len(left_pos) < MIN_PARENT:
                    reasons[slot]['whole_parent_below_minimum'] += 1; continue
                a, b = a0, b0
                if context in (1, 2):
                    wanted = context == 1
                    if sum(near[left_pos] == wanted) < MIN_PARENT or sum(near[right_pos] == wanted) < MIN_PARENT:
                        reasons[slot]['context_parent_below_minimum'] += 1; continue
                    keep = (near[left_pos[a]] == wanted) & (near[right_pos[b]] == wanted)
                    a, b = a[keep], b[keep]
                if len(a) < MIN_EDGES:
                    reasons[slot]['context_cross_edges_below_minimum'] += 1; continue
                vectors = residual if context == 3 else binary
                x, y = vectors[parent], vectors[2]
                if x is None or y is None:
                    reason = (residual_reasons[parent] or residual_reasons[2]) if context == 3 else 'missing_author_Ki67_probability'
                    reasons[slot][reason] += 1; continue
                selected_l = (near[left_pos] == (context == 1)) if context in (1, 2) else np.ones(len(x), dtype=bool)
                selected_r = (near[right_pos] == (context == 1)) if context in (1, 2) else np.ones(len(y), dtype=bool)
                if np.ptp(x[selected_l]) <= 1e-10 or np.ptp(y[selected_r]) <= 1e-10:
                    reasons[slot]['no_state_variation_in_context_parent'] += 1; continue
                gx, gy = group[left_pos], group[right_pos]
                mx, my = means_by_group(x, gx), means_by_group(y, gy)
                metric = covariance if context == 3 else lambda u, v, i, j: float(np.mean(u[i]*v[j]))
                observed = metric(x, y, a, b); expected = metric(mx, my, a, b)
                draws = np.empty(B)
                for j in range(B):
                    draws[j] = metric(shuffled(x, gx, rng), shuffled(y, gy, rng), a, b)
                records[slot].append(dict(observed=observed, expected=expected, draws=draws, edges=len(a)))
    result = aggregate(records, reasons, diagnostics)
    result['source_threshold_equal_cells'] = ties_at_threshold
    result['source_materials_considered'] = int(frame['sample.id'].nunique())
    return result


def analyze(bundle, parameters):
    assert parameters == dict(scope='fixed_followup_context_background_no_clinical', permutations=B, seed=SEED)
    folder = Path(__file__).parent
    qualifier = load(folder/'qualify_moldoveanu_ICI_cells_in_RAM_v1.py',
                     'a498faeb15053036e72003585ba82d5c9b44b2f8f9d3c16d392c3fcf8e693125', 'qualified_followup_cells')
    qualification = qualifier.analyze(bundle, dict(coordinate_tolerance=1e-6, minimum_parent_cells=20))
    archive = bundle['archive']; frame = archive._own_ICI_qualified_cells
    body = qualifier.decoded(qualifier.table_bytes(archive, 'ICI_meanIntensity.tsv.gz'))
    columns = ['Cell_ID', 'Sample_ID', 'Cluster.v2', 'prob.FoxP3', 'prob.CD45RO'] + BACKGROUND
    raw = pd.read_csv(io.BytesIO(body), sep='\t', usecols=columns, dtype={'Cell_ID': str, 'Sample_ID': str})
    # Avoid duplicated Histone columns already provided by the frozen cell reader.
    raw = raw.drop(columns=[x for x in BACKGROUND if x in frame])
    checked = frame.merge(raw, left_on=['sample.id', 'obj.id'], right_on=['Sample_ID', 'Cell_ID'],
                          how='left', validate='one_to_one', suffixes=('', '_raw'), indicator='_followup_join')
    assert len(checked) == len(frame) and checked['_followup_join'].eq('both').all()
    code = qualifier.table_bytes(archive, 'cell_type_clustering.R')
    assert hashlib.sha256(code).hexdigest() == '2d38c485cfdc5e5292d92d5f4ce4f11a08c95767a2eaa382aa92a077c840487f'
    audit = audit_final_author_phenotypes(checked)
    assert audit['author_probability_recipes_all_passed']
    helpers = load(folder/'analyze_melanoma_IMC_continuous_state_neighbors_v2.py',
                   '51b8a4323d806ae011ecb7d9419d56548eba63ca947a6756d3ab1dc5f9baf7bf', 'fixed_rank_grid_followup')
    neighbors = load(folder/'melanoma_IMC_numpy_neighbors_v1.py',
                     'acbde69dae387ee79d6d8397899681bedf00553dd3da08878abd16b392858471', 'fixed_exact_graph_followup')
    result = statistics(checked, lambda coords: exact_edges(coords, neighbors), helpers.rank_z, helpers.tile_groups)
    return dict(status='COMPLETE_EXPLORATORY_SOURCE_CONTEXT_SENSITIVITY_NOT_CM02_OR_CLINICAL_VALIDATION',
                completed_UTC=datetime.now(timezone.utc).isoformat(), source_qualification=qualification,
                final_phenotype_audit=audit, results=result, family_size=8, permutations=B, seed=SEED,
                binary_rule='Published ST4 Ki67.prob >= 0.5, threshold ties reported',
                null='Within material x parent x grid4 x one-hop-tumor-neighbor status',
                tumor_neighbor='One edge to author melano in native k10 union; not micron contact or pathology compartment',
                nuclear_background_fields=BACKGROUND,
                all_background_covariates_required_no_automatic_dropping=True,
                background_markers_can_carry_real_cell_cycle_biology=True,
                binary_HH_enrichment_not_proof_of_unique_HH_mechanism=True,
                patient_identity_and_phase_unreviewed_no_patient_CI_or_clinical_model=True,
                residual_and_binary_effect_units_are_not_comparable=True,
                original_closed_families_unchanged=True, arrays_RAM_only=True,
                producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
