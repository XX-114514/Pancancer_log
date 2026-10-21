"""Fixed source-unit partial-rank associations; no outcome or identity fitting.

Inputs and resampling arrays must stay in RAM. Zero/constant readouts and
unidentified designs remain missing. The complete nuisance column space is
preserved; individual collinear nuisance coefficients are not interpreted.
"""
import hashlib
import numpy as np


def seed(label):
    return int.from_bytes(hashlib.sha256(('2026102201|' + label).encode()).digest()[:8], 'little')


def ranks(a):
    a = np.asarray(a, float)
    if a.ndim == 2:
        if a.shape[1] == 0:
            return a.copy()
        return np.column_stack([ranks(v) for v in a.T])
    order = np.argsort(a, kind='stable')
    values = a[order]
    ends = np.r_[0, np.flatnonzero(values[1:] != values[:-1]) + 1, len(a)]
    result = np.empty(len(a))
    for lo, hi in zip(ends[:-1], ends[1:]):
        result[order[lo:hi]] = (lo + hi - 1) / 2 + 1
    return result


def standard(a):
    a = np.asarray(a, float)
    sd = np.std(a, axis=0)
    return (a - np.mean(a, axis=0)) / np.where(sd > 0, sd, 1)


def design(x, c, min_n=20, min_df=10):
    x, c = np.asarray(x, float), np.asarray(c, float)
    if len(x) < min_n:
        return None, 'FEWER_THAN20_SOURCE_UNITS'
    if not np.isfinite(np.c_[x, c]).all():
        return None, 'MISSING_INPUT'
    z = standard(ranks(x))
    raw = np.c_[np.ones(len(x)), standard(ranks(c))]
    u, s, _ = np.linalg.svd(raw, full_matrices=False)
    rank = int(np.sum(s > s[0] * 1e-12))
    if rank != np.linalg.matrix_rank(raw) or s[rank - 1] / s[0] < 1e-8:
        return None, 'UNSTABLE_NUISANCE_COLUMN_SPACE'
    q = u[:, :rank]
    assert np.max(np.abs(raw - q @ (q.T @ raw))) < 1e-9
    rx = z - q @ (q.T @ z)
    if np.linalg.norm(rx) <= 1e-8 * max(1., np.linalg.norm(z)):
        return None, 'EXPOSURE_CONSTANT_OR_EXPLAINED_BY_CONTROLS'
    d = np.c_[q, z]
    if len(x) - rank - 1 < min_df:
        return None, 'RESIDUAL_DF_BELOW10'
    if np.linalg.cond(d) >= 1e8:
        return None, 'ILL_CONDITIONED_DESIGN'
    inverse = np.linalg.pinv(d, rcond=1e-12)
    h = np.sum(d * inverse.T, axis=1)
    if h.max() >= .999999:
        return None, 'UNIT_LEVERAGE_TOO_HIGH'
    return dict(q=q, z=z, rx=rx, d=d, inverse=inverse, h=h,
                h0=np.sum(q*q, axis=1), df=len(x)-rank-1, rank=rank), 'ESTIMABLE'


def fit(m, y, transform=True):
    y = standard(ranks(y)) if transform else np.asarray(y, float)
    ry = y - m['q'] @ (m['q'].T @ y)
    if np.linalg.norm(ry) <= 1e-8 * max(1., np.linalg.norm(y)):
        return None, 'READOUT_CONSTANT_OR_EXPLAINED_BY_CONTROLS'
    operator = m['inverse'][-1]
    residual = y - m['d'] @ (m['inverse'] @ y)
    variance = np.sum((operator * residual / (1-m['h']))**2)
    if not variance > 1e-20:
        return None, 'ZERO_ROBUST_VARIANCE'
    beta = float(operator @ y)
    rho = float(m['rx'] @ ry / (np.linalg.norm(m['rx']) * np.linalg.norm(ry)))
    return dict(beta=beta, se=float(np.sqrt(variance)), t=beta/np.sqrt(variance),
                rho=float(np.clip(rho, -1, 1)), y=y), 'ESTIMABLE'


def association(x, y, c):
    m, status = design(x, c)
    if m is None:
        return None, None, status
    result, status = fit(m, y)
    return m, result, status


def wild_statistics(m, f, signs):
    null_residual = f['y'] - m['q'] @ (m['q'].T @ f['y'])
    adjusted = null_residual / (1-m['h0'])
    operator = m['inverse'][-1]
    student = operator / (1-m['h'])
    answer = np.empty(signs.shape[0])
    for start in range(0, len(answer), 1000):
        noise = adjusted[:, None] * signs[start:start+1000].T
        errors = noise - m['d'] @ (m['inverse'] @ noise)
        se = np.sqrt(np.sum((student[:, None] * errors)**2, axis=0))
        answer[start:start+noise.shape[1]] = (operator @ noise) / se
    assert np.isfinite(answer).all()
    return answer


def adjustments(values, by=False):
    p = np.asarray(values, float)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order]*len(p)/np.arange(1, len(p)+1))[::-1])[::-1]
    if by:
        q *= np.sum(1/np.arange(1, len(p)+1))
    out = np.empty(len(p)); out[order] = np.minimum(1, q)
    return out


def bootstrap(x, y, c, indices):
    values = np.full(len(indices), np.nan)
    reasons = []
    for j, take in enumerate(indices):
        _, f, why = association(x[take], y[take], c[take])
        reasons.append(why)
        if f is not None:
            values[j] = f['rho']
    return values, reasons


def calibration(m, label):
    """All fixed trials; finite working-model alarm, not a validity guarantee."""
    rows = []
    for scenario in ['normal', 'heteroskedastic_normal']:
        rng = np.random.default_rng(seed('calibration|' + label + '|' + scenario))
        pvalues = []
        for _ in range(256):
            scale = 1 if scenario == 'normal' else np.exp(.5*m['z'])
            y = m['q'] @ np.linspace(-.3, .3, m['rank']) + rng.normal(size=len(m['z']))*scale
            f, status = fit(m, y, transform=False)
            assert status == 'ESTIMABLE'
            signs = np.where(rng.random((999, len(y))) < .5, -1., 1.)
            null = wild_statistics(m, f, signs)
            pvalues.append(float((1+np.sum(np.abs(null) >= abs(f['t'])-1e-12))/1000))
        rate = float(np.mean(np.asarray(pvalues) < .05))
        rows.append(dict(scenario=scenario, trials=256, alpha=.05, rejection_rate=rate,
                         passed_finite_alarm=rate <= .10, simulation_pvalues=pvalues))
    return rows
