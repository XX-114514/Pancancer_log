"""Edge-centered continuous-protein state organization in retained RAM.

Author adjacency variant, fixed before any real state effect.
Prepared algorithm, not a frozen real-data endpoint declaration. Exact author
identity labels must be qualified outcome-blind before a request is created.
No clinical outcomes are read. Relative ranks are not protein positivity,
new formal cell identities, direct communication, or functional assays.
"""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

R_MARKERS = ['Ki-67', 'LAG3', 'OX40']
E_MARKERS = ['Ki-67', 'PD1', 'CD127']
PERMUTATIONS = 2000
BOOTSTRAPS = 2000
SEED = 2026102102
PARENT_MIN = 20
ROI_PARENT_MIN = 5
ROI_EDGE_MIN = 20
SOURCE_MIN = 5


def rank_z(values):
    """Average ranks, mean zero, unit RMS; constant columns stay invalid."""
    values = np.asarray(values, dtype=float)
    assert values.ndim == 2 and np.isfinite(values).all() and len(values) > 0
    ranks = pd.DataFrame(values).rank(method='average').to_numpy()
    centered = ranks-ranks.mean(axis=0)
    rms = np.sqrt(np.mean(centered**2, axis=0))
    valid = rms > 0
    return np.divide(centered, rms, out=np.zeros_like(centered), where=valid), valid


def tile_groups(coords, extent, bins=4):
    """Common ROI grid; invariant to translation and positive axis scaling."""
    lower, upper = extent
    span = upper-lower
    scaled = np.divide(coords-lower, span, out=np.zeros_like(coords, dtype=float), where=span > 0)
    indices = np.minimum(np.maximum(np.floor(scaled*bins).astype(int), 0), bins-1)
    return indices[:, 0]*bins+indices[:, 1]


def shuffled(values, groups, rng):
    """Permute whole co-marker vectors together within the given strata."""
    out = values.copy()
    for group in np.unique(groups):
        selected = np.flatnonzero(groups == group)
        out[selected] = values[rng.permutation(selected)]
    return out


def group_means(values, groups):
    out = np.empty_like(values)
    for group in np.unique(groups):
        selected = groups == group
        out[selected] = values[selected].mean(axis=0)
    return out


def edge_covariance(left, right, a, b):
    assert len(a) == len(b) and len(a) > 0
    return left[a].T @ right[b]/len(a)-np.outer(left[a].mean(axis=0), right[b].mean(axis=0))


def roi_permutations(left, right, a, b, left_groups, right_groups, count, rng):
    """Analytic null center plus independent within-parent vector shuffles."""
    observed = edge_covariance(left, right, a, b)
    expected = edge_covariance(group_means(left, left_groups), group_means(right, right_groups), a, b)
    draws = np.empty((count, left.shape[1], right.shape[1]), dtype=float)
    for index in range(count):
        draws[index] = edge_covariance(shuffled(left, left_groups, rng),
            shuffled(right, right_groups, rng), a, b)
    return observed, expected, draws


def oriented_edges(edges, labels, left_label, right_label):
    a, b = edges[:, 0], edges[:, 1]
    forward = (labels[a] == left_label) & (labels[b] == right_label)
    reverse = (labels[a] == right_label) & (labels[b] == left_label)
    return np.r_[a[forward], b[reverse]], np.r_[b[forward], a[reverse]]


def bh(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind='stable')
    q = np.empty(len(values))
    q[order] = np.minimum(1, np.minimum.accumulate((values[order]*len(values)/np.arange(1, len(values)+1))[::-1])[::-1])
    return q


def finish_family(observed, expected, draws, roi_counts, family, seed):
    """ROI-equal within source case, then equal case weights, no ROI pseudo-n."""
    rows = []
    p_values = []
    for index, candidate in enumerate(family):
        take = roi_counts[:, index] > 0
        n = int(take.sum())
        row = dict(candidate, source_case_units=n,
            contributing_ROIs=int(roi_counts[take, index].sum()),
            clinical_twenty_case_coverage_necessary_only=n >= 20)
        if n < SOURCE_MIN:
            row.update(status='INSUFFICIENT_SOURCE_CASES', effect=None, P=None, BH=None)
            rows.append(row)
            p_values.append(1.0)
            continue
        denominator = roi_counts[take, index]
        effects = (observed[take, index]-expected[take, index])/denominator
        null = (draws[take, :, index]-expected[take, index, None])/denominator[:, None]
        point = float(effects.mean())
        null_source = null.mean(axis=0)
        p = float((1+np.count_nonzero(np.abs(null_source) >= abs(point)-1e-12))/(len(null_source)+1))
        ci = None
        loo = None
        if n >= 10:
            rng = np.random.default_rng(seed+index)
            samples = rng.integers(0, n, size=(BOOTSTRAPS, n))
            ci = np.quantile(effects[samples].mean(axis=1), [.025, .975]).tolist()
            loo_values = (effects.sum()-effects)/(n-1)
            loo = [float(loo_values.min()), float(loo_values.max())]
        row.update(status='ESTIMABLE_CONDITIONAL_ORGANIZATION', effect=point, P=p, BH=None,
            case_bootstrap_95_CI=ci, leave_one_case_effect_range=loo,
            positive_case_count=int((effects > 0).sum()), negative_case_count=int((effects < 0).sum()),
            small_case_interval_suppressed=n < 10,
            P_scope='Conditional state randomization in these observed source cases; not a population clinical test')
        p_values.append(p)
        rows.append(row)
    for row, q in zip(rows, bh(p_values)):
        if row['P'] is not None:
            row['BH'] = float(q)
    return rows


def analyze(bundle, parameters):
    assert set(parameters) == {'regulatory_author_label', 'conventional_CD4_author_label', 'effector_author_label'}
    role_labels = [parameters[k] for k in ['regulatory_author_label', 'conventional_CD4_author_label', 'effector_author_label']]
    assert len(set(role_labels)) == 3 and all(isinstance(x, str) and x for x in role_labels)
    root = Path(__file__).parent.parent/'private_output'
    measurement = json.loads((root/'private_frozen_step_author_node_native_geometry_v1_source_summary.json').read_text())
    assert measurement['status'] == 'COMPLETE_AUTHOR_NODE_BIJECTION_AND_EDGE_PREFIX_NATIVE_COORDINATE_REVIEW'
    assert measurement['independent_CSV_Pandas_node_replay'] and measurement['complete_node_member_CRC_checked']
    frame = bundle['frame']
    assert set(role_labels) <= set(frame['celltype20'])
    coords = bundle['native_coordinates']
    mapping = bundle['ROI_to_source_patient']
    cases = sorted(set(mapping.values()))
    case_index = {value: index for index, value in enumerate(cases)}
    family = [dict(regulatory_context=context, left_author_label=label, right_author_label=role_labels[2],
        left_marker=a, right_marker=b) for context, label in zip(['author_regulatory', 'author_conventional_CD4'], role_labels[:2])
        for a in R_MARKERS for b in E_MARKERS]
    selected = frame['sample_id'].isin(mapping)
    counts = frame.loc[selected, ['sample_id', 'celltype20']].copy()
    counts['source_case'] = counts['sample_id'].map(mapping)
    counts = counts.groupby(['source_case', 'celltype20']).size().unstack(fill_value=0).reindex(index=cases, columns=role_labels, fill_value=0)
    storage = {}
    rngs = {}
    for mode, offset in [('ROI_parent_shuffle', 0), ('ROI_parent_grid4_shuffle', 100000)]:
        storage[mode] = dict(observed=np.zeros((len(cases), 18)), expected=np.zeros((len(cases), 18)),
            draws=np.zeros((len(cases), PERMUTATIONS, 18)), roi_counts=np.zeros((len(cases), 18), dtype=int))
        rngs[mode] = np.random.default_rng(SEED+offset)
    opportunities = dict(explicit_margin_ROIs=len(mapping), role_pair_ROI_opportunities=[0, 0],
        low_cell_ROIs=[0, 0], low_cross_edge_ROIs=[0, 0], no_outcome_values_read=True)
    for roi, indices in frame.groupby('sample_id', sort=True).indices.items():
        if roi not in mapping:
            continue
        case = mapping[roi]
        case_slot = case_index[case]
        roi_frame = frame.iloc[indices]
        labels = roi_frame['celltype20'].to_numpy()
        roi_coords = coords[indices]
        edges = None
        for context, left_label in enumerate(role_labels[:2]):
            if counts.loc[case, left_label] < PARENT_MIN or counts.loc[case, role_labels[2]] < PARENT_MIN:
                continue
            left_positions = np.flatnonzero(labels == left_label)
            right_positions = np.flatnonzero(labels == role_labels[2])
            if min(len(left_positions), len(right_positions)) < ROI_PARENT_MIN:
                opportunities['low_cell_ROIs'][context] += 1
                continue
            if edges is None:
                edges = bundle['author_neighbor_edges'][roi]
            left_edges, right_edges = oriented_edges(edges, labels, left_label, role_labels[2])
            if len(left_edges) < ROI_EDGE_MIN:
                opportunities['low_cross_edge_ROIs'][context] += 1
                continue
            opportunities['role_pair_ROI_opportunities'][context] += 1
            a_map = np.full(len(labels), -1, dtype=int)
            b_map = np.full(len(labels), -1, dtype=int)
            a_map[left_positions] = np.arange(len(left_positions))
            b_map[right_positions] = np.arange(len(right_positions))
            a, b = a_map[left_edges], b_map[right_edges]
            assert (a >= 0).all() and (b >= 0).all()
            left, valid_left = rank_z(roi_frame.iloc[left_positions][R_MARKERS].to_numpy())
            right, valid_right = rank_z(roi_frame.iloc[right_positions][E_MARKERS].to_numpy())
            valid = (valid_left[:, None] & valid_right[None, :]).reshape(-1)
            destinations = np.arange(context*9, (context+1)*9)[valid]
            if not len(destinations):
                continue
            extent = (roi_coords.min(axis=0), roi_coords.max(axis=0))
            for mode, state in storage.items():
                left_groups = tile_groups(roi_coords[left_positions], extent) if 'grid4' in mode else np.zeros(len(left), dtype=int)
                right_groups = tile_groups(roi_coords[right_positions], extent) if 'grid4' in mode else np.zeros(len(right), dtype=int)
                observed, expected, draws = roi_permutations(left, right, a, b,
                    left_groups, right_groups, PERMUTATIONS, rngs[mode])
                state['observed'][case_slot, destinations] += observed.reshape(-1)[valid]
                state['expected'][case_slot, destinations] += expected.reshape(-1)[valid]
                state['draws'][case_slot][:, destinations] += draws.reshape(PERMUTATIONS, 9)[:, valid]
                state['roi_counts'][case_slot, destinations] += 1
    results = {mode: finish_family(**state, family=family, seed=SEED+200000) for mode, state in storage.items()}
    return dict(status='COMPLETE_CONDITIONAL_CONTINUOUS_STATE_ORGANIZATION_NOT_CLINICAL_VALIDATION',
        completed_UTC=datetime.now(timezone.utc).isoformat(),
        producer_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        measurement_bridge_SHA256=hashlib.sha256((root/'private_frozen_step_author_node_native_geometry_v1_source_summary.json').read_bytes()).hexdigest(),
        family_size=18, permutations=PERMUTATIONS, case_bootstraps=BOOTSTRAPS, seed=SEED,
        graph='Complete author-supplied undirected within-ROI graph; paper defines20um centroid neighborhood, threshold not independently reconstructed',
        rank_semantics='Within-author-identity/ROI continuous marker rank with edge-endpoint centering; positive covariance includes high-high AND low-low, not high-high alone',
        null_semantics='Author adjacency and identity counts fixed; whole co-marker rows permuted within each parent and ROI. Native-coordinate grid4 is sensitivity only; exact Position provenance unresolved',
        family=results, opportunities=opportunities,
        direct_comparator_difference_not_tested=True, no_Treg_specificity_from_one_significant_one_nonsignificant=True,
        no_positive_intensity_threshold_or_new_formal_identity=True,
        endpoint_marginal_state_affinity_removed_by_edge_centering=True,
        no_communication_causality_function_or_clinical_increment_claim=True,
        no_clinical_values_read_or_individual_arrays_written=True)
