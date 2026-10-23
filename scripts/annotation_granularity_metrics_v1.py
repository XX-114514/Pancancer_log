"""Source/lineage granularity and expression-unit-out RNA diagnostics.

Expression units are not assumed to be patients. Graph discovery is transductive;
this diagnostic does not replace the separate source-heldout reference analysis.
"""
import numpy as np


def granularity(labels):
    labels = np.asarray(labels)
    known = labels >= 0
    _, size = np.unique(labels[known], return_counts=True)
    total = int(size.sum())
    if not total:return {'status': 'unpartitioned', 'n_cells': len(labels), 'n_partitioned': 0}
    p = size / total
    return {'status': 'estimated', 'n_cells': len(labels), 'n_partitioned': total,
            'n_occupied_groups': len(size), 'effective_group_count': float(np.exp(-np.sum(p * np.log(p)))),
            'largest_group_fraction': float(p.max()), 'fraction_in_groups_below20': float(size[size < 20].sum() / total),
            'median_group_size': float(np.median(size)), 'minimum_group_size': int(size.min()),
            'partitioned_fraction': float(known.mean())}


def unit_out_state(labels, units, depth, values, minimum_other_cells=20, return_predictions=False):
    """Predict a cell's RNA readout from its group in other expression units.

    For each held-out unit, depth regression is fitted on other units only.
    Group residual means also exclude the entire held-out unit. Report the
    improvement over that same training-only depth baseline on identical cells.
    No group/threshold is selected from outcomes here.
    """
    labels, units = np.asarray(labels), np.asarray(units)
    depth, values = np.asarray(depth, float), np.asarray(values, float)
    if any(len(x) != len(labels) for x in [units,depth,values]):raise ValueError('State diagnostic axes differ')
    good = (labels >= 0) & np.isfinite(depth) & np.isfinite(values)
    pred, baseline = np.full(len(labels), np.nan), np.full(len(labels), np.nan)
    if len(np.unique(units[good])) < 2:
        out = {'status': 'fewer_than_two_expression_units', 'n_evaluable': 0, 'n_input_cells': len(labels)}
        if return_predictions:out.update(prediction=pred, baseline=baseline)
        return out
    codes = np.unique(labels[good]);group_map = {int(k): j for j,k in enumerate(codes)}
    groups = np.array([group_map[int(k)] for k in labels[good]])
    valid_rows = np.flatnonzero(good);u = units[good]
    x = np.column_stack([np.ones(good.sum()), depth[good]])
    y = values[good];k = len(codes)
    count = np.bincount(groups, minlength=k)
    ysum = np.bincount(groups, weights=y, minlength=k)
    xsum = np.zeros((k,2));np.add.at(xsum,groups,x)
    fitted_units = 0
    for unit in np.unique(u):
        query = u == unit;train = ~query
        if train.sum() < 30:continue
        beta = np.linalg.lstsq(x[train], y[train], rcond=None)[0]
        query_count = np.bincount(groups[query], minlength=k)
        query_ysum = np.bincount(groups[query], weights=y[query], minlength=k)
        query_xsum = np.zeros((k,2));np.add.at(query_xsum, groups[query], x[query])
        available = count - query_count
        residual_sum = ysum - query_ysum - (xsum - query_xsum) @ beta
        retained = query & (available[groups] >= minimum_other_cells)
        if not retained.any():continue
        rows = valid_rows[retained];base = x[retained] @ beta
        baseline[rows] = base
        pred[rows] = base + residual_sum[groups[retained]] / available[groups[retained]]
        fitted_units += 1
    ok = np.isfinite(pred);n = int(ok.sum())
    out = {'status': 'insufficient_comparable_cells', 'n_input_cells': len(labels), 'n_evaluable': n,
           'evaluable_fraction': float(ok.mean()) if len(ok) else 0., 'n_expression_units_evaluable_not_patients': fitted_units,
           'all_query_unit_rows_excluded_from_depth_and_group_prediction': True,
           'graph_discovery_included_query_cells': True, 'independent_biological_validation': False}
    if n >= 100:
        base_sse = float(np.sum((values[ok] - baseline[ok]) ** 2))
        sse = float(np.sum((values[ok] - pred[ok]) ** 2))
        out.update(status='no_variable_depth_residual' if base_sse < 1e-10 else 'estimated',
                   depth_baseline_SSE=base_sse, group_plus_depth_SSE=sse,
                   RNA_prediction_improvement_over_depth=1 - sse / base_sse if base_sse >= 1e-10 else None)
    if return_predictions:out.update(prediction=pred, baseline=baseline)
    return out


def gamma_group_retention(labels, values):
    labels, values = np.asarray(labels), np.asarray(values,float)
    keep = (labels >= 0) & np.isfinite(values)
    rare = values[keep] > 0;lab = labels[keep]
    n = int(rare.sum())
    if not len(lab):return {'status':'no_measured_partitioned_cells'}
    _, codes = np.unique(lab, return_inverse=True)
    counts = np.bincount(codes);positives = np.bincount(codes,weights=rare)
    eligible = rare & (counts[codes] > 1)
    return {'status':'RNA_candidate_group_diagnostic_not_verified_identity',
            'n_measured_cells':len(lab),'n_candidates':n,'candidate_prevalence':float(rare.mean()),
            'candidate_same_group_fraction_excluding_self':float(np.mean((positives[codes[eligible]]-1)/(counts[codes[eligible]]-1))) if eligible.any() else None,
            'n_singleton_candidates':int((rare & (counts[codes]==1)).sum())}
