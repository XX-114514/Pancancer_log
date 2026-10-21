"""Fixed matched-material representation/null comparison, source aggregates.

Twelve positions: two parent relationships x continuous rank covariance,
binary standardized covariance, binary HH edge fraction x two spatial nulls.
All six positions per relationship use exactly the same qualified materials.
Comparisons of effects or P values across positions are descriptive, not a
test of differences or a causal adjustment. No clinical endpoints are read.
"""
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import importlib.util
import io
import json
import numpy as np
import pandas as pd

BASE = Path(__file__).with_name('analyze_moldoveanu_Ki67_context_sensitivity_v1.py')
assert hashlib.sha256(BASE.read_bytes()).hexdigest() == '2c0d67729a64e9a9a1831640572f41a458135c21c49ec572f956067795a90352'
SPEC = importlib.util.spec_from_file_location('frozen_context_sensitivity_parent', BASE)
base = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(base)

SEED = 2026102202
B = 2000
MIN_PARENT = 20
MIN_EDGES = 20
MIN_MATERIALS = 5
METRICS = ['continuous_rank_covariance', 'binary_standardized_covariance',
           'binary_double_positive_edge_fraction']
NULLS = ['grid4', 'grid4_and_tumor_neighbor']


def profile_background(frame, rank):
    """Outcome-blind variation/design qualification; no background effect fit."""
    profiles = []
    for role, labels in zip(['regulatory', 'conventional_CD4', 'CD8'], base.PARENTS):
        field_counts = {field: Counter() for field in base.BACKGROUND}
        joint = Counter(); probabilities = Counter()
        for _, part in frame.loc[frame['Cluster'].isin(labels)].groupby('sample.id', sort=True):
            if len(part) < MIN_PARENT:
                joint['parent_below_minimum'] += 1
                continue
            varied = []
            for field in base.BACKGROUND:
                values = part[field].to_numpy(dtype=float)
                state = 'nonfinite' if not np.isfinite(values).all() else 'constant' if np.ptp(values) == 0 else 'variable'
                field_counts[field][state] += 1
                varied.append(state == 'variable')
            prob = part['Ki67.prob'].to_numpy(dtype=float)
            probabilities['source_cells_equal_threshold'] += int(np.sum(prob == .5))
            for name, positive in [('inclusive_ge_half', prob >= .5), ('strict_gt_half', prob > .5)]:
                probabilities[name+'_parents_with_binary_variation'] += int(np.any(positive) and not np.all(positive))
            if not all(varied):
                joint['constant_or_nonfinite_required_background_no_covariate_dropped'] += 1
            else:
                scores, _ = rank(part[base.BACKGROUND].to_numpy(dtype=float))
                design = np.column_stack([np.ones(len(part)), scores])
                if np.linalg.matrix_rank(design) != design.shape[1]:
                    joint['required_background_design_rank_deficient'] += 1
                elif np.linalg.cond(design) > 1e8:
                    joint['required_background_design_ill_conditioned'] += 1
                else:
                    joint['required_background_design_identifiable'] += 1
        profiles.append(dict(role=role, per_channel_source_material_counts=[
            dict(field=field, counts=dict(counts)) for field, counts in field_counts.items()],
            joint_design_counts=dict(joint), threshold_qualification_counts=dict(probabilities)))
    return dict(status='COMPLETE_OUTCOME_BLIND_BACKGROUND_AND_THRESHOLD_INPUT_QUALIFICATION',
                roles=profiles, no_background_effects_or_clinical_outcomes_computed=True,
                no_automatic_control_dropping_or_selection_by_effect=True)


def permute_indices(groups, rng):
    indices = np.arange(len(groups))
    for group in np.unique(groups):
        take = np.flatnonzero(groups == group)
        indices[take] = rng.permutation(take)
    return indices


def three_statistics(scores_left, scores_right, binary_left, binary_right, a, b):
    hh = float(np.mean(binary_left[a]*binary_right[b]))
    ordinary = sum(float(binary_left[i])*float(binary_right[j]) for i, j in zip(a, b))/len(a)
    assert abs(hh-ordinary) < 1e-9
    return np.array([base.covariance(scores_left[:, 0], scores_right[:, 0], a, b),
                     base.covariance(scores_left[:, 1], scores_right[:, 1], a, b), hh])


def summarize(records, coverage):
    rows = []
    for parent in range(2):
        selected = records[parent]
        for mode in range(2):
            for metric in range(3):
                row = dict(parent_relation=['regulatory_CD8', 'conventional_CD4_CD8'][parent],
                           null= NULLS[mode], metric=METRICS[metric],
                           common_source_materials=len(selected), P=None, BH=None, effect=None,
                           patient_CI=None, independent_patient_identity_reviewed=False,
                           missing_reason_counts=dict(coverage[parent]))
                if len(selected) < MIN_MATERIALS:
                    row['status'] = 'INSUFFICIENT_MATCHED_MATERIALS'
                else:
                    observed = [r['observed'][metric] for r in selected]
                    expected = [r['expected'][mode, metric] for r in selected]
                    effects = [float(a-b) for a, b in zip(observed, expected)]
                    point = sum(effects)/len(effects)
                    draws = np.mean(np.stack([r['draws'][mode, :, metric]-r['expected'][mode, metric]
                                              for r in selected]), axis=0)
                    for j in range(B):
                        replay = sum(float(r['draws'][mode, j, metric]-r['expected'][mode, metric])
                                     for r in selected)/len(selected)
                        assert abs(replay-draws[j]) < 1e-9
                    row.update(status='EXPLORATORY_MATCHED_MATERIAL_CONDITIONAL_TEST', effect=float(point),
                               P=(1+int(np.sum(np.abs(draws) >= abs(point)-1e-12)))/(B+1),
                               material_equal_observed=float(np.mean(observed)),
                               material_equal_null_expected=float(np.mean(expected)),
                               positive_materials=sum(x > 0 for x in effects),
                               negative_materials=sum(x < 0 for x in effects),
                               zero_materials=sum(x == 0 for x in effects),
                               leave_one_material_effect_range=[float((sum(effects)-max(effects))/(len(effects)-1)),
                                                                float((sum(effects)-min(effects))/(len(effects)-1))])
                rows.append(row)
    for row, q in zip(rows, base.bh([r['P'] if r['P'] is not None else 1. for r in rows])):
        if row['P'] is not None:
            row['BH'] = float(q)
    return rows


def statistics(frame, knn, rank, grid):
    records = [[], []]; coverage = [Counter(), Counter()]
    source_totals = int(frame['sample.id'].nunique())
    # Separate fixed random streams for nulls. Every within-parent shuffle moves
    # continuous and binary scores together, and CD8 indices are shared by roles.
    rngs = [np.random.default_rng(SEED), np.random.default_rng(SEED+100000)]
    for _, part in frame.groupby('sample.id', sort=True):
        part = part.reset_index(drop=True)
        positions = [np.flatnonzero(part['Cluster'].isin(labels)) for labels in base.PARENTS]
        if len(positions[2]) < MIN_PARENT or max(len(positions[0]), len(positions[1])) < MIN_PARENT:
            for reason in coverage:reason['whole_parent_below_minimum'] += 1
            continue
        coords = part[['coord.x', 'coord.y']].to_numpy(dtype=float)
        edges = knn(coords); near = base.tumor_neighbors(edges, part['Cluster'].eq('melano').to_numpy())
        tiles = grid(coords, (coords.min(0), coords.max(0)))
        pairs = [base.pair_indices(edges, positions[i], positions[2], len(part)) for i in range(2)]
        values = []; binary = []; reasons = []
        for pos in positions:
            if len(pos) < MIN_PARENT:
                values.append(None); binary.append(None); reasons.append('whole_parent_below_minimum');continue
            selected = part.iloc[pos]
            prob = selected['Ki67.prob'].to_numpy(dtype=float)
            if not np.isfinite(prob).all():
                values.append(None); binary.append(None); reasons.append('missing_author_probability');continue
            positive = (prob >= .5).astype(float)
            scores, varied = rank(np.column_stack([selected['Ki67'].to_numpy(dtype=float), positive]))
            reasons.append('no_continuous_variation' if not varied[0] else 'no_binary_variation' if not varied[1] else None)
            values.append(scores); binary.append(positive)
        valid = []
        for parent in range(2):
            reason = reasons[parent] or reasons[2]
            if reason is None and len(pairs[parent][0]) < MIN_EDGES:reason = 'cross_edges_below_minimum'
            valid.append(reason is None)
            if reason is not None:coverage[parent][reason] += 1
        if not any(valid):continue
        observed = [three_statistics(values[i], values[2], binary[i], binary[2], *pairs[i])
                    if valid[i] else None for i in range(2)]
        expected = [np.zeros((2, 3)) for _ in range(2)]
        draws = [np.zeros((2, B, 3)) for _ in range(2)]
        for mode in range(2):
            all_groups = tiles if mode == 0 else tiles*2+near.astype(int)
            groups = [all_groups[p] for p in positions]
            mean_scores = [np.column_stack([base.means_by_group(v[:, j], g) for j in range(2)])
                           if v is not None else None for v, g in zip(values, groups)]
            mean_binary = [base.means_by_group(v, g) if v is not None else None for v, g in zip(binary, groups)]
            for parent in range(2):
                if valid[parent]:expected[parent][mode] = three_statistics(mean_scores[parent], mean_scores[2],
                    mean_binary[parent], mean_binary[2], *pairs[parent])
            for iteration in range(B):
                permutations = [permute_indices(g, rngs[mode]) if values[i] is not None else None
                                for i, g in enumerate(groups)]
                for parent in range(2):
                    if valid[parent]:
                        left, right = permutations[parent], permutations[2]
                        draws[parent][mode, iteration] = three_statistics(values[parent][left], values[2][right],
                            binary[parent][left], binary[2][right], *pairs[parent])
        for parent in range(2):
            if valid[parent]:records[parent].append(dict(observed=observed[parent], expected=expected[parent], draws=draws[parent]))
    return dict(family=summarize(records, coverage), source_materials_considered=source_totals,
                every_relation_uses_identical_materials_edges_across_its_six_positions=True,
                jointly_permuted_representations_and_shared_CD8_indices=True,
                all_covariances_and_HH_scalar_checked=True, permutation_source_means_scalar_replayed=True)


def analyze(bundle, parameters):
    assert set(parameters) == {'scope', 'permutations', 'seed', 'checkpoint_directory'}
    assert parameters['scope'] == 'matched_representation_null_no_clinical'
    assert parameters['permutations'] == B and parameters['seed'] == SEED
    folder = Path(__file__).parent
    qualifier = base.load(folder/'qualify_moldoveanu_ICI_cells_in_RAM_v1.py',
                          'a498faeb15053036e72003585ba82d5c9b44b2f8f9d3c16d392c3fcf8e693125', 'matched_cells_qualified')
    source = qualifier.analyze(bundle, dict(coordinate_tolerance=1e-6, minimum_parent_cells=20))
    archive = bundle['archive']; frame = archive._own_ICI_qualified_cells
    raw_bytes = qualifier.decoded(qualifier.table_bytes(archive, 'ICI_meanIntensity.tsv.gz'))
    fields = ['Cell_ID', 'Sample_ID', 'Cluster.v2', 'prob.FoxP3', 'prob.CD45RO'] + base.BACKGROUND
    raw = pd.read_csv(io.BytesIO(raw_bytes), sep='\t', usecols=fields, dtype={'Cell_ID': str, 'Sample_ID': str})
    raw = raw.drop(columns=[x for x in base.BACKGROUND if x in frame])
    checked = frame.merge(raw, left_on=['sample.id', 'obj.id'], right_on=['Sample_ID', 'Cell_ID'],
                          how='left', validate='one_to_one', suffixes=('', '_raw'), indicator='_matched_join')
    assert len(checked) == len(frame) and checked['_matched_join'].eq('both').all()
    assert hashlib.sha256(qualifier.table_bytes(archive, 'cell_type_clustering.R')).hexdigest() == '2d38c485cfdc5e5292d92d5f4ce4f11a08c95767a2eaa382aa92a077c840487f'
    phenotype = base.audit_final_author_phenotypes(checked)
    assert phenotype['author_probability_recipes_all_passed']
    helper = base.load(folder/'analyze_melanoma_IMC_continuous_state_neighbors_v2.py',
                       '51b8a4323d806ae011ecb7d9419d56548eba63ca947a6756d3ab1dc5f9baf7bf', 'matched_rank_grid')
    neighbors = base.load(folder/'melanoma_IMC_numpy_neighbors_v1.py',
                          'acbde69dae387ee79d6d8397899681bedf00553dd3da08878abd16b392858471', 'matched_graph')
    profile = profile_background(checked, helper.rank_z)
    profile['completed_UTC'] = datetime.now(timezone.utc).isoformat()
    out = Path(parameters['checkpoint_directory'])
    assert out.is_dir() and not out.is_symlink()
    target = out/'private_input_background_qualification_v1.json'
    with target.open('x') as handle:json.dump(profile, handle, indent=2, allow_nan=False)
    result = statistics(checked, lambda c:base.exact_edges(c, neighbors), helper.rank_z, helper.tile_groups)
    return dict(status='COMPLETE_MATCHED_REPRESENTATION_NULL_EXPLORATION_NOT_CLINICAL_VALIDATION',
                completed_UTC=datetime.now(timezone.utc).isoformat(), source_qualification=source,
                final_type_audit=phenotype, preeffect_background_qualification=profile, results=result,
                family_size=12, permutations=B, seed=SEED,
                covariance_inputs='Within-material parent mean-zero unit-RMS average ranks of Ki67 or binary state',
                binary_rule='Published ST4 Ki67.prob >= 0.5; strict > threshold availability only, no adaptive threshold test',
                different_metric_units_not_pooled=True, across_position_effect_changes_descriptive_not_formal_difference_tests=True,
                no_background_effect_fit_or_covariate_dropping=True,
                no_patient_CI_or_clinical_model=True, original_closed_families_unchanged=True,
                no_original_CM02_weight_projection_or_direct_communication_claim=True,
                arrays_RAM_only=True, producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
