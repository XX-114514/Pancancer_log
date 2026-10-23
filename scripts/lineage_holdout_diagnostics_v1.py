"""Training-only feature selection and paired, conditional RNA transfer probes.

Frozen lineage/direct-candidate labels are conditioning inputs, not test truth.
Neither query expression nor clinical outcomes choose model features or anchors.
"""
import numpy as np

from annotation_representation_pilot_v1 import select_features
from run_annotation_source_holdout_v1 import conditional_scores


VARIANTS = ['shared_training_features', 'lineage_training_features']


def fold_members(source, fold, positive, heldout):
    source, fold, positive = np.asarray(source), np.asarray(fold), np.asarray(positive, bool)
    if source.ndim != 1 or fold.shape != source.shape or positive.shape != source.shape:
        raise ValueError('Holdout axes differ')
    train = np.flatnonzero(positive & (fold != heldout))
    query = np.flatnonzero(positive & (fold == heldout))
    if not len(train) or not len(query):
        raise ValueError('Empty training or heldout scope')
    if set(source[train]) & set(source[query]):
        raise ValueError('Source split across training and query')
    return train, query


def training_feature_sets(x, source, route, train, wanted, number, shared=None):
    train = np.asarray(train, np.int64)
    if len(np.unique(train)) != len(train):
        raise ValueError('Training rows repeated')
    local = train[np.asarray(route)[train] == wanted]
    if not len(local):
        raise ValueError('Training lineage unavailable')
    global_features = select_features(x[train], np.asarray(source)[train], number) if shared is None else shared
    lineage_features = select_features(x[local], np.asarray(source)[local], number)
    return dict(zip(VARIANTS, [np.asarray(global_features), np.asarray(lineage_features)]))


def paired_profiles(source, direct, truth, totals, predictions, supported, programs, minimum):
    """Every contrast uses the same evaluable query cells for each readout.

    Cell-level correlations are descriptive diagnostics. No independent-patient
    P values, biological identity accuracy, or clinical validation are claimed.
    """
    source, direct, truth, totals = map(np.asarray, [source, direct, truth, totals])
    if truth.shape != (len(source), len(programs)) or direct.shape != source.shape or totals.shape != source.shape:
        raise ValueError('Paired profile axes differ')
    rows = []
    scopes = ['all_route_cells_not_pure_type'] + sorted(set(direct) - {'unresolved'})
    for scope in scopes:
        if scope not in predictions:
            continue
        wanted = np.ones(len(source), bool) if scope == 'all_route_cells_not_pure_type' else direct == scope
        for ds in np.unique(source[wanted]):
            subset = wanted & (source == ds)
            for j, program in enumerate(programs):
                common = subset & np.isfinite(truth[:, j]) & np.isfinite(totals)
                for variant in VARIANTS:
                    p = predictions[scope][variant]
                    if p.shape != truth.shape or np.asarray(supported[variant]).shape != source.shape:
                        raise ValueError('Prediction/support axes differ')
                    common &= supported[variant] & np.isfinite(p[:, j])
                metrics = {variant: conditional_scores(truth[common, j], predictions[scope][variant][common, j],
                           totals[common], minimum) for variant in VARIANTS}
                first, second = [metrics[v] for v in VARIANTS]
                record = {'source': str(ds), 'identity_scope': str(scope), 'program': str(program),
                    'n_query_scope_cells': int(subset.sum()), 'n_common_query_cells': int(common.sum()),
                    'common_query_fraction': float(common.sum() / subset.sum()), 'metrics': metrics,
                    'status': 'estimated' if all(m['status'] == 'estimated' for m in metrics.values()) else 'insufficient_coverage'}
                if all(m.get('depth_adjusted_structure_status') == 'estimated' for m in metrics.values()):
                    record['delta_lineage_minus_shared_depth_adjusted_Pearson'] = (
                        second['depth_adjusted_Pearson'] - first['depth_adjusted_Pearson'])
                rows.append(record)
    return rows
