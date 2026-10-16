"""Complete-control-space rank regression and fixed wild bootstrap kernel."""
import numpy as np
from scipy.stats import rankdata, beta as beta_distribution


def geometry(x, controls, min_n=20, min_df=10, require_leverage=True):
    x = np.asarray(x, float); controls = np.asarray(controls, float)
    n = len(x)
    if n < min_n:
        raise ValueError('INSUFFICIENT_LIBRARY_UNITS')
    if controls.ndim != 2 or len(controls) != n or not np.isfinite(controls).all() or not np.isfinite(x).all():
        raise ValueError('MISSING_INPUT_NOT_ZERO')
    sd = controls.std(axis=0)
    z = np.column_stack([np.ones(n), (controls - controls.mean(axis=0)) / np.where(sd > 1e-12, sd, 1)])
    u, s, _ = np.linalg.svd(z, full_matrices=False)
    rz = int(np.sum(s > s[0] * 1e-12)); q0 = u[:, :rz]
    xx = rankdata(x).astype(float)
    if xx.std() <= 1e-12:
        raise ValueError('CONSTANT_LIGAND')
    xx = (xx - xx.mean()) / xx.std()
    rx = xx - q0 @ (q0.T @ xx)
    ssx = float(rx @ rx)
    singular = np.linalg.svd(np.column_stack([z, xx]), compute_uv=False)
    rank = int(np.sum(singular > singular[0] * 1e-12))
    if rank != rz + 1 or ssx <= 1e-16 * float(xx @ xx):
        raise ValueError('LIGAND_NOT_IDENTIFIABLE_GIVEN_COMPLETE_CONTROLS')
    condition = float(singular[0] / singular[rank - 1])
    if condition >= 1e8 or n - rank < min_df:
        raise ValueError('LOW_EFFECTIVE_DF_OR_ILL_CONDITIONED')
    q = np.column_stack([q0, rx / np.sqrt(ssx)])
    h = np.einsum('ij,ij->i', q, q)
    if require_leverage and h.max() >= .999999:
        raise ValueError('EXTREME_LEVERAGE')
    return dict(n=n, z=z, x=xx, q0=q0, q=q, rx=rx, ssx=ssx, h=h,
                rank=rank, nuisance_rank=rz, residual_df=n-rank, condition=condition,
                strict_full_column_rank=(rank == z.shape[1] + 1), constant_controls=sd <= 1e-12)


def point(x, y, controls, min_n=20, min_df=10, require_leverage=True):
    g = geometry(x, controls, min_n, min_df, require_leverage)
    y = np.asarray(y, float)
    if not np.isfinite(y).all() or np.ptp(y) <= 1e-12:
        raise ValueError('MISSING_OR_CONSTANT_RECEIVER_READOUT')
    yy = rankdata(y).astype(float); yy = (yy - yy.mean()) / yy.std()
    ry = yy - g['q0'] @ (g['q0'].T @ yy)
    sy = float(ry @ ry)
    if sy <= 1e-16 * float(yy @ yy):
        raise ValueError('NO_RECEIVER_VARIATION_AFTER_COMPLETE_CONTROLS')
    coefficient = float(g['rx'] @ yy / g['ssx'])
    residual = yy - g['q'] @ (g['q'].T @ yy)
    variance = np.nan; statistic = np.nan
    if require_leverage:
        variance = float(np.sum((g['rx'] * residual / (1-g['h'])) ** 2) / g['ssx'] ** 2)
        if not np.isfinite(variance) or variance <= 0:
            raise ValueError('DEGENERATE_HC3_VARIANCE')
        statistic = coefficient / np.sqrt(variance)
    rho = float(g['rx'] @ ry / np.sqrt(g['ssx'] * sy))
    g.update(y=yy, residual=residual, coefficient=coefficient, HC3_variance=variance,
             HC3_t=statistic, partial_rank_rho=rho)
    return g


def wild_statistics(g, signs):
    h0 = np.einsum('ij,ij->i', g['q0'], g['q0'])
    fitted0 = g['q0'] @ (g['q0'].T @ g['y'])
    residual0 = (g['y'] - fitted0) / np.sqrt(1-h0)
    y = fitted0[:, None] + residual0[:, None] * signs
    b = g['rx'] @ y / g['ssx']
    e = y - g['q'] @ (g['q'].T @ y)
    variance = np.sum((g['rx'][:, None] * e / (1-g['h'][:, None])) ** 2, axis=0) / g['ssx'] ** 2
    t = b / np.sqrt(variance)
    assert np.isfinite(t).all()
    return t


def wild_count(g, draws, seed, batch=1000):
    rng = np.random.default_rng(seed); exceed = 0
    for start in range(0, draws, batch):
        signs = np.where(rng.random((min(batch, draws-start), g['n'])) < .5, -1., 1.).T
        statistics = wild_statistics(g, signs)
        exceed += int(np.sum(np.abs(statistics) >= abs(g['HC3_t']) - 1e-12))
    return exceed


def monte_carlo_upper(exceed, draws, hypotheses):
    upper = 1. if exceed == draws else float(beta_distribution.ppf(1-.01/hypotheses, exceed+1, draws-exceed))
    return min(1., (1 + draws * upper) / (draws + 1))


def bootstrap_rhos(x, y, controls, draws, seed):
    rng = np.random.default_rng(seed); result = np.full(draws, np.nan)
    for i in range(draws):
        take = rng.integers(0, len(x), size=len(x))
        try:
            result[i] = point(x[take], y[take], controls[take], min_n=1, min_df=1, require_leverage=False)['partial_rank_rho']
        except ValueError:
            pass
    return result


def leave_one_out(x, y, controls):
    result = np.full(len(x), np.nan)
    for i in range(len(x)):
        use = np.arange(len(x)) != i
        try:
            result[i] = point(x[use], y[use], controls[use], min_n=1, min_df=1, require_leverage=False)['partial_rank_rho']
        except ValueError:
            pass
    return result
