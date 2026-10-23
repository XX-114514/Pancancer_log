"""Provisional discovery fit pools; neither archive deletion nor release gates."""
import numpy as np


def boolean_evidence(values, n):
    if values is None:
        return np.zeros(n, bool), np.zeros(n, bool)
    values = np.asarray(values)
    if values.shape != (n,): raise ValueError('Carried boolean evidence axis differs')
    normalized = np.char.lower(np.char.strip(values.astype(str)))
    true = np.isin(normalized, ['true', '1', '1.0', 'yes'])
    false = np.isin(normalized, ['false', '0', '0.0', 'no'])
    return true, true | false


def fit_pool(positive, context_member, expression_units, doublet=None, excluded=None):
    positive, context = np.asarray(positive, bool), np.asarray(context_member, bool)
    units = np.asarray(expression_units)
    if positive.ndim != 1 or context.shape != positive.shape or units.shape != positive.shape:
        raise ValueError('Fit-pool input axes differ')
    n = len(positive)
    doublet_positive, doublet_known = boolean_evidence(doublet, n)
    exclusion_positive, exclusion_known = boolean_evidence(excluded, n)
    valid_unit = ~np.isin(np.char.lower(np.char.strip(units.astype(str))),
        ['', 'nan', 'none', 'na', 'n/a', 'unknown', '__missing__', '__source_value_missing__'])
    return {'fit_pool_candidate': positive & context & valid_unit & ~doublet_positive & ~exclusion_positive,
        'positive_native_RNA': positive, 'in_frozen_context': context, 'declared_expression_unit_available': valid_unit,
        'carried_doublet_positive': doublet_positive, 'carried_doublet_assessed_value': doublet_known,
        'carried_primary_exclusion_positive': exclusion_positive, 'carried_primary_exclusion_assessed_value': exclusion_known}
