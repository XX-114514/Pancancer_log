"""Require informative repeated estimates before state-feature candidacy."""
import math

def candidate_status(row,n_splits,n_positive_units,minimum_units=10,minimum_nonzero_units=4,
                     minimum_finite_fraction=.8,minimum_median=.6,minimum_positive_fraction=.8):
    if row['status']!='estimable':return False,row['status']
    if row['n_units']<minimum_units:return False,'insufficient_within_source_context_units'
    if row['n_finite_splits']<math.ceil(n_splits*minimum_finite_fraction):return False,'insufficient_finite_split_estimates'
    if n_positive_units<minimum_nonzero_units:return False,'too_few_units_with_observed_positive_programme_score'
    median=row['split_half_Spearman_median'];positive=row['positive_split_fraction']
    if not math.isfinite(median) or not math.isfinite(positive):return False,'nonfinite_reliability_estimate'
    if median<minimum_median or positive<minimum_positive_fraction:return False,'split_reliability_below_project_measurement_rule'
    return True,'measurement_candidate_scope_and_source_specific'
