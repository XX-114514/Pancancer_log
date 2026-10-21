"""Patient-level partial rank descriptions and independent SVD replay.

No records or learned arrays are written by this module. Intervals remain
pointwise descriptions; there are no P values or clinical claims.
"""
import hashlib
import numpy as np
import pandas as pd
from scipy.stats import rankdata

MIN_N = 20
B = 2000
CONTINUOUS = [1, 2, 5, 6, 7, 8]


def estimate(y, covariates, reference=False):
    n, k = covariates.shape
    empty = np.full((3, 3), np.nan)
    if n < MIN_N:
        return empty, 'INSUFFICIENT20_REVIEWED_CASES'
    if n - k - 1 < 10:
        return empty, 'INSUFFICIENT10_RESIDUAL_DF'
    c = covariates.copy()
    ranking = (lambda x: pd.DataFrame(x).rank(method='average').to_numpy()) if reference else (lambda x: rankdata(x, axis=0))
    for j in CONTINUOUS:
        if j < k:
            v = ranking(c[:, j:j+1]).ravel()
            scale = v.std()
            if scale == 0:
                return empty, 'CONSTANT_CONTROL'
            c[:, j] = (v - v.mean()) / scale
    u, s, vt = np.linalg.svd(c, full_matrices=False)
    if np.linalg.matrix_rank(c) != k:
        return empty, 'DESIGN_NOT_FULL_RANK'
    if s[0] / s[-1] >= 1e6:
        return empty, 'DESIGN_CONDITION_ATLEAST1E6'
    a = ranking(y)
    if reference:
        e = a - u @ (u.T @ a)
    else:
        e = a - c @ np.linalg.lstsq(c, a, rcond=1e-12)[0]
        q, _ = np.linalg.qr(c, mode='reduced')
        np.testing.assert_allclose(e, a-q@(q.T@a), atol=1e-9, rtol=1e-9)
    norm = np.linalg.norm(e, axis=0)
    rho = e[:, :3].T @ e[:, 3:] / np.maximum(norm[:3, None]*norm[None, 3:], 1e-30)
    rho[(norm[:3, None] < 1e-10) | (norm[None, 3:] < 1e-10)] = np.nan
    return np.clip(rho, -1, 1), 'DESCRIPTIVE_PARTIAL_RANK_ASSOCIATION'


def checked(y, c):
    first, status = estimate(y, c)
    second, other = estimate(y, c, reference=True)
    assert status == other
    np.testing.assert_allclose(first, second, atol=1e-9, rtol=1e-9, equal_nan=True)
    return first, status


def finite(value):
    return float(value) if np.isfinite(value) else None


def profile(y, c, strata, denominators, label):
    # Every representation uses the exact same eligible patients.
    main, main_status = checked(y, c)
    base, base_status = checked(y, c[:, :5])
    raw, raw_status = checked(y, c[:, :1])
    sensitivities = {name: checked(z, c)[0] for name, z in denominators.items()}
    draws = {name: [] for name in ['capture', 'base']}
    loos = {name: [] for name in ['capture', 'base']}
    rng = np.random.default_rng(int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], 'little'))
    blocks = [np.flatnonzero(strata == value) for value in np.unique(strata)]
    if np.isfinite(main).any() or np.isfinite(base).any():
        for _ in range(B):
            ix = np.concatenate([rng.choice(block, len(block), replace=True) for block in blocks])
            for name, cc in [('capture', c), ('base', c[:, :5])]:
                if np.isfinite(main if name == 'capture' else base).any():
                    draws[name].append(checked(y[ix], cc[ix])[0])
    if len(y) - 1 >= MIN_N:
        for i in range(len(y)):
            for name, cc in [('capture', c), ('base', c[:, :5])]:
                loos[name].append(checked(np.delete(y, i, 0), np.delete(cc, i, 0))[0])
    output = []
    for i in range(3):
        for j in range(3):
            def status_for(global_status, point):
                return 'STATE_OR_RESIDUAL_INVARIANT' if global_status == 'DESCRIPTIVE_PARTIAL_RANK_ASSOCIATION' and not np.isfinite(point) else global_status
            record = dict(i=i, j=j, reviewed_cases=len(y), capture_status=status_for(main_status,main[i,j]),
                          base_status=status_for(base_status,base[i,j]), unadjusted_status=status_for(raw_status,raw[i,j]),
                          capture_adjusted_partial_Spearman=finite(main[i,j]),
                          abundance_treatment_HPV_partial_Spearman=finite(base[i,j]),
                          unadjusted_Spearman=finite(raw[i,j]),
                          sensitivity_partial_Spearman={name: finite(z[i,j]) for name,z in sensitivities.items()},
                          no_P_FDR_or_confirmatory_promotion=True)
            for name in ['capture', 'base']:
                v = np.array([x[i,j] for x in draws[name]]);v=v[np.isfinite(v)]
                lv = np.array([x[i,j] for x in loos[name]]);lv=lv[np.isfinite(lv)]
                point = main[i,j] if name == 'capture' else base[i,j]
                record[name+'_valid_bootstraps'] = len(v)
                record[name+'_pointwise_CI'] = np.quantile(v, [.025,.975]).tolist() if len(v) >= 1800 else None
                record[name+'_valid_leave_one'] = len(lv)
                record[name+'_leave_one_range'] = [float(lv.min()),float(lv.max())] if len(lv) else None
                record[name+'_leave_one_same_direction'] = int((lv*point>0).sum()) if np.isfinite(point) else None
            output.append(record)
    return output


def synthetic():
    rng = np.random.default_rng(2026102121)
    c = np.column_stack([np.ones(32), rng.integers(1,10,(32,2)), np.tile([0,0,1,1],8), np.tile([0,1],16), rng.integers(1,12,(32,4))]).astype(float)
    y = rng.integers(0,8,(32,6)).astype(float)
    z, status = checked(y,c)
    assert np.isfinite(z).all()
    constant=y.copy();constant[:,0]=1
    assert np.isnan(checked(constant,c)[0][0]).all()
    assert checked(y[:19],c[:19])[1]=='INSUFFICIENT20_REVIEWED_CASES'
    singular=c.copy();singular[:,4]=singular[:,3]
    assert checked(y,singular)[1]=='DESIGN_NOT_FULL_RANK'
    # Unadjusted positive and negative tied-rank fixtures.
    v=np.tile([0.,0.,1.,2.],8)
    f=np.column_stack([v,v,v,v,-v,v])
    a,_=checked(f,np.ones((32,1)))
    np.testing.assert_allclose(a[:,0],1,atol=1e-12)
    np.testing.assert_allclose(a[:,1],-1,atol=1e-12)
    strata = 2*c[:,3]+c[:,4]
    full = profile(y,c,strata,{'same_representation':y},'synthetic_complete_resampling')
    assert len(full)==9
    assert all(r['capture_valid_bootstraps']>=1800 and r['base_valid_bootstraps']>=1800 for r in full)
    return dict(tied_ranks_positive_negative_constant_and_patient_gate=True,
                full_covariate_rank_deficiency_rejected_without_dropping=True,
                lstsq_QR_vs_Pandas_rank_SVD_exact=True,
                complete2000_stratified_bootstrap_and_LOO_synthetic_replayed=True)
