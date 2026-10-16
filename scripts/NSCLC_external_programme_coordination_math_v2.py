"""Fixed statistic with exact-zero certificates for cancellation-prone draws.

Version1 remains immutable. Scientific cutoffs and random streams unchanged.
"""
import hashlib

import numpy as np
from scipy.stats import rankdata
import NSCLC_external_exact_rank_guards_v1 as exact_guard


def seed(slot, stream):
    return int.from_bytes(hashlib.sha256((slot+'|'+stream+'|external_programme_coordination_v1').encode()).digest()[:8], 'little')


def design(quality, blocks):
    quality = np.asarray(quality, float)
    if quality.ndim == 1:
        quality = quality[:, None]
    if not np.isfinite(quality).all():
        raise ValueError('NONFINITE_QUALITY')
    r = rankdata(quality, axis=0, method='average') if quality.shape[1] else quality
    mean = r.mean(axis=0)
    sd = r.std(axis=0)
    z = (r - mean) / np.where(sd > 0, sd, 1)
    levels = sorted(set(blocks))
    categorical = np.column_stack([np.asarray(blocks) == s for s in levels]).astype(float)
    return np.column_stack([np.ones(len(blocks)), categorical, z])


def basis(z):
    u, s, _ = np.linalg.svd(z, full_matrices=False)
    rank = int(np.sum(s > s[0] * 1e-10))
    return u[:, :rank], rank


def prepare(a, b, quality, blocks, min_n=20, min_df=8):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    quality = np.asarray(quality, float)
    blocks = np.asarray(blocks, str)
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise ValueError('NONFINITE_PROGRAMME')
    keep = np.array([np.sum(blocks == s) >= 2 for s in blocks])
    if keep.sum() < min_n:
        raise ValueError('FEWER_THAN20_UNITS_IN_NON_SINGLETON_BLOCKS' if min_n == 20 else 'INSUFFICIENT_INFORMATIVE_UNITS')
    a, b, quality, blocks = a[keep], b[keep], quality[keep], blocks[keep]
    z = design(quality, blocks)
    q, rank = basis(z)
    if len(a) - rank < min_df:
        raise ValueError('QUALITY_AND_BLOCK_DESIGN_INSUFFICIENT_RESIDUAL_DF')
    ar, br = rankdata(a, method='average'), rankdata(b, method='average')
    u, v = ar - q @ (q.T @ ar), br - q @ (q.T @ br)
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    # A trigger for exact diagnosis, not a new residual cutoff.
    if min(nu, nv) < 1e-6:
        qr = rankdata(quality, axis=0) if quality.shape[1] else quality
        zi = np.column_stack([np.ones(len(a)), *[(blocks == s).astype(int) for s in sorted(set(blocks))], 2*qr]).astype(np.int64)
        yi = np.column_stack([2*ar, 2*br]).astype(np.int64)
        joint = np.unique(np.column_stack([zi, yi]), axis=0)
        zz = joint[:, :zi.shape[1]]
        certified = len(joint) == rank and exact_guard.modular_rank(zz) == rank
        if certified or (exact_guard.rational_rank(zz) == rank and exact_guard.any_exactly_fitted(zi, yi)):
            raise ValueError('INDIVIDUAL_PROGRAMME_RESIDUAL_CONSTANT')
    if nu < 1e-10 or nv < 1e-10:
        raise ValueError('INDIVIDUAL_PROGRAMME_RESIDUAL_CONSTANT')
    return dict(a=a, b=b, quality=quality, blocks=blocks, keep=keep, design=z, Q=q,
                design_rank=rank, residual_df=len(a)-rank, u=u, v=v,
                norm_u=float(nu), norm_v=float(nv))


def effect(prepared):
    return float(np.clip(prepared['u'] @ prepared['v'] / (prepared['norm_u']*prepared['norm_v']), -1, 1))


def permutation_count(p, draws, random_seed):
    rng = np.random.default_rng(random_seed)
    q, u, v, blocks = p['Q'], p['u'], p['v'], p['blocks']
    obs = abs(effect(p))
    count = 0
    groups = [np.flatnonzero(blocks == s) for s in sorted(set(blocks))]
    for start in range(0, draws, 256):
        m = min(256, draws-start)
        take = np.broadcast_to(np.arange(len(v)), (m, len(v))).copy()
        for ids in groups:
            take[:, ids] = ids[np.argsort(rng.random((m, len(ids))), axis=1)]
        shuffled = v[take]
        projected = shuffled - (shuffled @ q) @ q.T
        den = np.linalg.norm(projected, axis=1)*p['norm_u']
        # Degenerate shuffled outcomes are explicitly undefined, never P=negative.
        if np.any(den < 1e-12):
            raise ValueError('DEGENERATE_CONDITIONAL_PERMUTATION')
        values = projected @ u / den
        count += int((np.abs(values) >= obs-1e-12).sum())
    return count


def bootstrap(p, draws, random_seed):
    rng = np.random.default_rng(random_seed)
    groups = [np.flatnonzero(p['blocks'] == s) for s in sorted(set(p['blocks']))]
    out = np.full(draws, np.nan)
    for i in range(draws):
        take = np.concatenate([ids[rng.integers(len(ids), size=len(ids))] for ids in groups])
        try:
            out[i] = effect(prepare(p['a'][take], p['b'][take], p['quality'][take], p['blocks'][take]))
        except ValueError:
            pass
    return out


def leave_one(p):
    out = np.full(len(p['a']), np.nan)
    for i in range(len(out)):
        take = np.arange(len(out)) != i
        try:
            out[i] = effect(prepare(p['a'][take], p['b'][take], p['quality'][take], p['blocks'][take], min_n=3))
        except ValueError:
            pass
    return out
