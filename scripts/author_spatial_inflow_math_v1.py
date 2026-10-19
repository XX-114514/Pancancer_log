"""Transparent LIANA-defined inflow with fixed-field conditional mark nulls."""
import hashlib
import itertools
import numpy as np
from scipy.sparse import csr_matrix, diags
from scipy.spatial import cKDTree


def seed(*parts):
    return int.from_bytes(hashlib.sha256(('20261019_spatial_inflow|' + '|'.join(map(str, parts))).encode()).digest()[:8], 'little')


def graph(xy, k, radius_multiplier):
    assert len(xy) > 2 and len(np.unique(xy, axis=0)) == len(xy)
    distance, neighbor = cKDTree(xy).query(xy, k=min(k + 1, len(xy)))
    spacing = float(np.median(distance[:, 1]))
    assert spacing > 0
    row = np.broadcast_to(np.arange(len(xy))[:, None], neighbor[:, 1:].shape)
    keep = distance[:, 1:] <= radius_multiplier * spacing
    A = csr_matrix((np.ones(keep.sum()), (row[keep], neighbor[:, 1:][keep])), shape=(len(xy), len(xy)))
    A = A.maximum(A.T)
    A.setdiag(0)
    A.eliminate_zeros()
    degree = np.asarray(A.sum(1)).ravel()
    W = diags(np.divide(1., degree, out=np.zeros_like(degree), where=degree > 0)) @ A
    np.testing.assert_allclose(np.asarray(W.sum(1)).ravel(), degree > 0, atol=1e-12)
    return W.tocsr(), dict(undirected_edges=A.nnz // 2, isolated_nodes=int((degree == 0).sum()), median_nearest_neighbor_distance=spacing, radius_in_source_units=radius_multiplier * spacing)


def bins(total, area, xy, mode):
    depth = np.searchsorted(np.unique(np.quantile(total, [.25, .5, .75])), total, side='right')
    morphology = area > np.median(area)
    code = depth + 4 * morphology
    if mode == 'DEPTH_AREA_XY':
        code = code + 8 * (xy[:, 0] > np.median(xy[:, 0])) + 16 * (xy[:, 1] > np.median(xy[:, 1]))
    elif mode != 'DEPTH_AREA':
        raise ValueError('UNDECLARED_NULL')
    return [np.flatnonzero(code == v) for v in sorted(set(code))]


def group_mean(values, groups):
    out = np.empty_like(values)
    for ids in groups:
        out[ids] = values[ids].mean()
    return out


def permutation(values, groups, count, rng):
    out = np.repeat(values[:, None], count, axis=1)
    for ids in groups:
        if len(ids) > 1:
            order = np.argsort(rng.random((len(ids), count)), axis=0, kind='stable')
            out[ids] = values[ids[order]]
    return out


def describe(W, L, R):
    near = W @ L
    opportunity = np.asarray(W.sum(1)).ravel()
    observed = float(np.mean(R * near))
    lm, rm = float(L.mean()), float(R.mean())
    centered = W @ (L - lm)
    terms = dict(structural_opportunity=lm * rm * float(opportunity.mean()), sender_mark=rm * float(centered.mean()), receiver_mark=lm * float(np.mean((R - rm) * opportunity)), joint_centered_marks=float(np.mean((R - rm) * centered)))
    np.testing.assert_allclose(sum(terms.values()), observed, atol=1e-11, rtol=1e-11)
    coo = W.tocoo()
    explicit = float(np.sum(coo.data * R[coo.row] * L[coo.col]) / len(R))
    np.testing.assert_allclose(explicit, observed, atol=1e-11, rtol=1e-11)
    take = opportunity > 0
    return dict(inflow=observed, mean_MPS_neighbor_fraction=float(opportunity.mean()), receiver_cells_with_MPS_neighbors=int(take.sum()), MPS_composition_normalized_inflow=float(np.mean(R[take] * near[take] / opportunity[take])) if take.any() else np.nan, **terms)


def null(W, L, R, lg, rg, draws, random_seed):
    expected = float(np.mean(group_mean(R, rg) * (W @ group_mean(L, lg))))
    rng = np.random.default_rng(random_seed)
    values = np.empty(draws)
    for start in range(0, draws, 128):
        count = min(128, draws - start)
        left = permutation(L, lg, count, rng)
        right = permutation(R, rg, count, rng)
        values[start:start + count] = np.mean(right * (W @ left), axis=0)
    return values, expected


def pvalue(observed, values):
    upper = (1 + np.count_nonzero(values >= observed - 1e-12)) / (len(values) + 1)
    lower = (1 + np.count_nonzero(values <= observed + 1e-12)) / (len(values) + 1)
    return min(1., 2 * min(upper, lower))


def fixture():
    W = csr_matrix([[.4, .1], [.2, .5]])
    L, R = np.array([1., 3.]), np.array([2., 4.])
    d = describe(W, L, R)
    exact = np.array([np.mean(np.asarray(rr) * (W @ np.asarray(ll))) for ll in itertools.permutations(L) for rr in itertools.permutations(R)])
    expected = float(np.mean(R.mean() * (W @ np.repeat(L.mean(), len(L)))))
    np.testing.assert_allclose(exact.mean(), expected, atol=1e-12)
    assert describe(W, L, np.zeros(2))['inflow'] == 0 and pvalue(0., np.zeros(19)) == 1.
    xy = np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.], [9., 9.]])
    g, meta = graph(xy, 2, 2.)
    assert g.diagonal().sum() == 0 and meta['isolated_nodes'] == 1
    return dict(status='PASS_EXPLICIT_EDGE_FLOW_FOUR_TERM_PARTITION_ENUMERATED_MARK_NULL_AND_ZERO_GRAPH_FIXTURE', exact_null_mean=expected, observed=d['inflow'])
