"""State-label spatial randomization on fixed, disjoint parent cell positions.

No data access, clinical outcomes, learned labels, or plotting in this module.
All overlapping child states move as one row inside a parent/stratum.
"""
from itertools import product, permutations

import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree


def quantile_bins(values):
    """Value-based tertiles; ties stay together, never cell-order tie breaking."""
    v = np.asarray(values, float)
    assert len(v) and np.isfinite(v).all()
    cuts = np.unique(np.quantile(v, [1/3, 2/3]))
    return np.searchsorted(cuts, v, side='left')


def groups_from_keys(keys):
    _, inverse = np.unique(np.asarray(keys), axis=0, return_inverse=True)
    return [np.flatnonzero(inverse == i) for i in range(inverse.max()+1)]


def label_probabilities(labels, groups):
    p = np.empty(labels.shape, float)
    covered = np.zeros(len(labels), int)
    for g in groups:
        assert len(g)
        p[g] = labels[g].mean(axis=0)
        covered[g] += 1
    assert (covered == 1).all()
    return p


def permute_joint(labels, groups, rng):
    result = labels.copy()
    for g in groups:
        result[g] = labels[rng.permutation(g)]
    return result


def graph_for_radius(left_xy, right_xy, radius):
    """One undirected cross-parent pair per edge; no same-parent/self edges."""
    assert radius > 0 and len(left_xy) and len(right_xy)
    neighbours = cKDTree(left_xy).query_ball_tree(cKDTree(right_xy), radius)
    lengths = np.fromiter((len(v) for v in neighbours), dtype=np.int64)
    rows = np.repeat(np.arange(len(neighbours)), lengths)
    cols = np.fromiter((j for v in neighbours for j in v), dtype=np.int64,
                       count=int(lengths.sum()))
    g = csr_matrix((np.ones(len(rows), np.int64), (rows, cols)),
                   shape=(len(left_xy), len(right_xy)))
    assert g.nnz == len(rows)
    return g, rows, cols


def edge_crosscounts(left, right, rows, cols):
    """Independent edge-list arithmetic; never materializes pairwise dense axes."""
    return np.asarray([[np.dot(left[rows, i], right[cols, j])
                        for j in range(right.shape[1])]
                       for i in range(left.shape[1])])


def measure_graphs(left, right, left_groups, right_groups, graphs,
                   permutations_count, seed):
    left = np.asarray(left, np.int64); right = np.asarray(right, np.int64)
    assert ((left == 0) | (left == 1)).all()
    assert ((right == 0) | (right == 1)).all()
    lp = label_probabilities(left, left_groups)
    rp = label_probabilities(right, right_groups)
    observed = []; expected = []
    for graph, rows, cols in graphs:
        obs = left.T @ (graph @ right)
        exp = lp.T @ (graph @ rp)
        np.testing.assert_array_equal(obs, edge_crosscounts(left, right, rows, cols))
        np.testing.assert_allclose(exp, edge_crosscounts(lp, rp, rows, cols),
                                   rtol=1e-12, atol=1e-8)
        observed.append(obs); expected.append(exp)
    draws = np.empty((permutations_count, len(graphs), left.shape[1],
                      right.shape[1]), np.int64)
    rng = np.random.default_rng(seed)
    # Structural checks are independent of whether a scientific contrast is strong.
    for b in range(permutations_count):
        l = permute_joint(left, left_groups, rng)
        r = permute_joint(right, right_groups, rng)
        if b < 8:
            for before, after, groups in [(left, l, left_groups), (right, r, right_groups)]:
                for g in groups:
                    np.testing.assert_array_equal(
                        np.sort(before[g] @ (1 << np.arange(before.shape[1]))),
                        np.sort(after[g] @ (1 << np.arange(after.shape[1]))))
        for k, (graph, rows, cols) in enumerate(graphs):
            draws[b, k] = l.T @ (graph @ r)
            if b < 8:
                np.testing.assert_array_equal(draws[b, k],
                    edge_crosscounts(l, r, rows, cols))
    observed = np.asarray(observed); expected = np.asarray(expected)
    assert (expected >= 0).all() and (observed >= 0).all()
    return dict(observed=observed, expected=expected, draws=draws,
                log2_enrichment=np.log2((observed+.5)/(expected+.5)),
                null_log2_enrichment=np.log2((draws+.5)/(expected[None, ...]+.5)),
                empirical_null_varies=np.ptp(draws, axis=0)>0,
                positive_expected=expected>0)


def synthetic_fixture():
    """Exhaustive small null and direct geometry check, including fixed strata."""
    lxy = np.array([[0,0], [2,0], [7,0], [9,0]], float)
    rxy = np.array([[1,0], [3,0], [8,0], [10,0]], float)
    left = np.array([[1,0,1], [0,1,1], [1,1,0], [0,0,0]], np.int64)
    right = np.array([[1,1,0], [0,1,1], [1,0,1], [0,0,0]], np.int64)
    lg = [np.array([0,1]), np.array([2,3])]
    rg = [np.array([0,1]), np.array([2,3])]
    graph, rows, cols = graph_for_radius(lxy, rxy, 2)
    direct = np.linalg.norm(lxy[:,None]-rxy[None,:], axis=2) <= 2
    np.testing.assert_array_equal(graph.toarray(), direct)
    draws = []
    orderings = list(product(list(permutations([0,1])), list(permutations([2,3]))))
    for lo, ro in product(orderings, orderings):
        l = left[np.concatenate(lo)]; r = right[np.concatenate(ro)]
        draws.append(edge_crosscounts(l,r,rows,cols))
    analytic = label_probabilities(left,lg).T @ (graph @ label_probabilities(right,rg))
    np.testing.assert_allclose(np.mean(draws,axis=0), analytic, rtol=0, atol=0)
    results = measure_graphs(left, right, lg, rg, [(graph,rows,cols)], 32, [20261021,0])
    # All-constant and all-singleton strata are permitted, but cannot identify
    # any state arrangement; this must not become a significant zero effect.
    fixed = measure_graphs(left, right, [np.array([i]) for i in range(4)],
        [np.array([i]) for i in range(4)], [(graph,rows,cols)], 16, [20261021,1])
    assert not fixed['empirical_null_varies'].any()
    assert quantile_bins([1,1,1,1]).tolist() == [0,0,0,0]
    assert len(results['draws']) == 32
    return dict(KDTree_edges_match_bruteforce=True,
                exact_enumeration_matches_analytic_null_expectation=True,
                joint_state_patterns_preserved=True,
                dense_and_edge_list_crosscounts_match=True,
                singleton_null_not_estimable=True,
                tied_quantile_values_kept_together=True)
