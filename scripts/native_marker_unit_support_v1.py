"""Descriptive within-source expression-unit support; never patient-level DE."""
import numpy as np


def summarize_unit_direction(own_n, own_sum, other_n, other_sum, minimum_cells=20):
    own_n, own_sum, other_n, other_sum = [np.asarray(x) for x in [own_n, own_sum, other_n, other_sum]]
    if own_n.ndim != 1 or any(x.shape != own_n.shape for x in [own_sum, other_n, other_sum]):
        raise ValueError('Marker-unit axes differ')
    if minimum_cells < 1 or any(not np.isfinite(x).all() or (x < 0).any() for x in [own_n, own_sum, other_n, other_sum]):
        raise ValueError('Invalid marker-unit denominator or value')
    valid = (own_n >= minimum_cells) & (other_n >= minimum_cells)
    ratio = np.log2((own_sum[valid] / own_n[valid] + .1) / (other_sum[valid] / other_n[valid] + .1))
    return {'n_units_with_group': int((own_n > 0).sum()), 'n_comparable_expression_units': int(valid.sum()),
            'minimum_cells_in_group_and_comparator': minimum_cells,
            'n_units_positive_mean_ratio': int((ratio > 0).sum()), 'n_units_negative_mean_ratio': int((ratio < 0).sum()),
            'median_log2_mean_ratio': float(np.median(ratio)) if len(ratio) else None,
            'descriptive_not_replicate_aware_DEG': True, 'expression_units_not_assumed_independent_patients': True}
