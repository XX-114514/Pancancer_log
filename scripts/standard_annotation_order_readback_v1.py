"""Independent CSR-window reconstruction and conditional statistic readback.

No biological identity, calibrated accuracy or independent patient evidence is
created by these numeric checks. Fractional adopted deltas remain fractional.
"""
import numpy as np
import scipy.sparse as sp
from scipy.stats import rankdata
from sklearn.neighbors import NearestNeighbors


def gather_windows(csr, rows, width=3072):
    requested = np.asarray(rows, dtype=np.int64)
    shape = tuple(map(int, csr.attrs['shape']))
    if requested.ndim != 1 or (len(requested) and (requested.min() < 0 or requested.max() >= shape[0])):
        raise ValueError('Native row axis outside sparse matrix')
    if not len(requested):
        return sp.csr_matrix((0, shape[1]), dtype=np.float64)
    ordered, inverse = np.unique(requested, return_inverse=True)
    tiles = ordered // width
    pieces = []
    for tile in np.unique(tiles):
        chosen = ordered[tiles == tile]
        low, high = int(chosen[0]), int(chosen[-1]) + 1
        pointer = np.asarray(csr['indptr'][low:high+1], dtype=np.int64)
        start, end = int(pointer[0]), int(pointer[-1])
        values = np.asarray(csr['data'][start:end], dtype=np.float64)
        columns = np.asarray(csr['indices'][start:end], dtype=np.int64)
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError('Nonfinite or negative stored counts/delta')
        matrix = sp.csr_matrix((values, columns, pointer-start), shape=(high-low, shape[1]))
        matrix.sum_duplicates()
        pieces.append(matrix[chosen-low])
    result = sp.vstack(pieces, format='csr')[inverse]
    result.eliminate_zeros()
    return result


def correlation(a, b):
    a = np.asarray(a, np.float64); b = np.asarray(b, np.float64)
    a = a-a.mean(); b = b-b.mean()
    return float(np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b)))


def independent_conditional_metrics(embedding, states, totals, candidate, indices, neighbors):
    indices = np.asarray(indices, np.int64)
    if len(indices) < max(30, neighbors+2):
        return dict(status='insufficient_conditional_cells')
    coordinates = np.asarray(embedding)[indices]
    _, near = NearestNeighbors(n_neighbors=neighbors+1, n_jobs=4).fit(coordinates).kneighbors(coordinates)
    near = np.asarray([r[r != i][:neighbors] for i,r in enumerate(near)])
    assert near.shape == (len(indices), neighbors)
    depth = np.log1p(np.asarray(totals)[indices])
    design = np.column_stack([np.ones(len(indices)), depth-depth.mean()])
    out = dict(status='conditional_reserved_RNA_diagnostic_not_orthogonal_truth', n_cells=len(indices), states={})
    for name, values in states.items():
        y = np.asarray(values)[indices]
        if not np.isfinite(y).all():
            out['states'][name] = dict(status='measurement_or_unused_gene_coverage_incomplete')
            continue
        residual = y-design@(np.linalg.pinv(design)@y)
        prediction = residual[near].mean(axis=1)
        if residual.std() < 1e-8 or prediction.std() < 1e-8:
            out['states'][name] = dict(status='no_variable_signal')
            continue
        centered = residual-residual.mean()
        # Explicit binary directed kNN Moran formula; not Scanpy's evaluator.
        denominator = np.dot(centered, centered)
        moran = np.sum(centered[:,None]*centered[near])/(neighbors*denominator)
        out['states'][name] = dict(status='estimated',
            depth_adjusted_neighbor_Pearson=correlation(residual,prediction),
            depth_adjusted_neighbor_Spearman=correlation(rankdata(residual),rankdata(prediction)),
            depth_adjusted_official_Scanpy_Morans_I=float(moran))
    rare = np.asarray(candidate)[indices]
    if not np.isfinite(rare).all():
        out['gamma_delta_RNA_candidate'] = dict(status='unmeasured_not_negative', independent_identity_accuracy_not_claimed=True)
    else:
        rare = rare.astype(bool); n = int(rare.sum())
        out['gamma_delta_RNA_candidate'] = dict(status='conditional_RNA_candidate_diagnostic',
            n_candidates_not_verified_identity=n,
            conditional_neighbor_candidate_fraction=float(rare[near[rare]].mean()) if n else None,
            conditional_candidate_prevalence=float(rare.mean()),
            candidate_RNA_genes_overlap_coarse_routing_not_independent_identity_validation=True)
    return out


def compare_record(actual, expected, atol=2e-5):
    if isinstance(actual,dict):
        if set(actual) != set(expected):raise ValueError('Different scientific result fields')
        for key in actual:compare_record(actual[key],expected[key],atol)
    elif isinstance(actual,float):
        if not np.isclose(actual,expected,atol=atol,rtol=atol):raise ValueError('Numeric conditional readback differs')
    elif actual != expected:
        raise ValueError('Status/coverage/axis readback differs')
