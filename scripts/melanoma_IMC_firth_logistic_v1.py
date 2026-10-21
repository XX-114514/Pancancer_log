"""Small-sample logistic utility; no data access or real clinical execution.

Jeffreys penalty handles separation. Comparing models of different dimension
requires a separately specified bootstrap statistic, not a naive penalized
likelihood difference interpreted as chi-square. Returned arrays are RAM-only
internal values, never a source-summary interface.
"""
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit


def standardized_design(values, minimum_cases=20, minimum_residual_df=10):
    values = np.asarray(values, dtype=float)
    assert values.ndim == 2 and np.isfinite(values).all()
    n, p = values.shape
    assert n >= minimum_cases and n-p-1 >= minimum_residual_df, 'insufficient_cases_or_residual_df'
    center = values.mean(axis=0)
    scale = values.std(axis=0, ddof=0)
    assert (scale > 0).all(), 'constant_covariate_not_dropped'
    design = np.column_stack([np.ones(n), (values-center)/scale])
    assert np.linalg.matrix_rank(design) == p+1, 'rank_deficient_design_not_reduced'
    assert np.linalg.cond(design) < 1e6, 'ill_conditioned_design_not_reduced'
    return design, center, scale


def objective(beta, design, outcome):
    eta = design @ beta
    probabilities = expit(eta)
    weight = probabilities*expit(-eta)
    information = design.T @ (weight[:, None]*design)
    sign, logdet = np.linalg.slogdet(information)
    if sign <= 0 or not np.isfinite(logdet):
        return np.inf, np.zeros(design.shape[1])
    log_likelihood = float(np.sum(outcome*eta-np.logaddexp(0, eta)))
    leverage = weight*np.sum(design*np.linalg.solve(information, design.T).T, axis=1)
    adjusted_score = design.T @ (outcome-probabilities+(.5-probabilities)*leverage)
    return -(log_likelihood+.5*logdet), -adjusted_score


def fit_firth(design, outcome):
    design = np.asarray(design, dtype=float)
    outcome = np.asarray(outcome, dtype=float)
    assert design.ndim == 2 and outcome.shape == (len(design),)
    assert np.isfinite(design).all() and np.isfinite(outcome).all()
    assert np.isin(outcome, [0., 1.]).all()
    assert np.linalg.matrix_rank(design) == design.shape[1]
    result = minimize(objective, np.zeros(design.shape[1]), args=(design, outcome),
        jac=True, method='BFGS', options=dict(gtol=1e-8, maxiter=500))
    value, gradient = objective(result.x, design, outcome)
    assert np.isfinite(value) and np.isfinite(result.x).all()
    assert np.max(np.abs(gradient)) < 1e-6, 'firth_fit_gradient_not_converged'
    eta = design @ result.x
    ll = float(np.sum(outcome*eta-np.logaddexp(0, eta)))
    return dict(coefficients=result.x, probabilities=expit(eta), log_likelihood=ll,
        penalized_log_likelihood=-float(value), maximum_score_error=float(np.max(np.abs(gradient))),
        optimizer_success=bool(result.success), accepted_by_explicit_score_check=True)
