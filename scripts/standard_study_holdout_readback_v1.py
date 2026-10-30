"""Preserved independent frozen-reference, distance and native-score checks."""
import json
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.spatial.distance import cdist
from scipy.stats import rankdata

def restore_reference(path, expected_genes, expected_indices):
    import anndata as ad
    with np.load(path, allow_pickle=False) as f:
        values = {k: f[k].copy() for k in f.files}
    if (values['genes'].tolist() != list(expected_genes) or
            not np.array_equal(values['feature_candidate_indices'], expected_indices)):
        raise ValueError('Frozen reference feature axis changed')
    genes, PCs, centers, nr = values['genes'], values['loadings'], values['C'], values['Nr']
    # symphonypy compresses centers as soft clusters × PCA coordinates.
    # This is not Harmony's alternate transposed representation.
    if PCs.ndim != 2 or PCs.shape[0] != len(genes) or centers.shape != (len(nr), PCs.shape[1]):
        raise ValueError('Reference PCA/Harmony axes differ')
    if values['sigma'].shape != nr.shape or (nr < 0).any() or not np.all(values['std'] == 1):
        raise ValueError('Invalid compressed reference parameters')
    for k in ['mean', 'std', 'loadings', 'C', 'Nr', 'sigma']:
        if not np.isfinite(values[k]).all():raise ValueError('Nonfinite reference parameter')
    if not np.allclose(PCs.T @ PCs, np.eye(PCs.shape[1]), atol=2e-4, rtol=2e-4):
        raise ValueError('Reference PCA loadings not orthonormal')
    reference = ad.AnnData(sp.csr_matrix((0, len(genes)), dtype=np.float32),
                          var=pd.DataFrame({'mean': values['mean'], 'std': values['std']}, index=genes))
    reference.varm['PCs'] = PCs
    reference.uns['harmony'] = {'C': centers, 'Nr': nr, 'K': len(nr), 'sigma': values['sigma'],
                                'ref_basis_loadings': 'PCs', 'ref_basis_adjusted': 'X_pca_harmony', 'converged': None}
    return reference, json.loads(str(values['metadata_json']))

def prototype_prediction(centers, outcomes, query, neighbors):
    """Independent Euclidean distances and inverse-distance reconstruction."""
    if len(centers) != len(outcomes) or centers.shape[1] != query.shape[1] or not len(centers):
        raise ValueError('Prototype/target/query axes differ')
    chunks = []
    for low in range(0, len(query), 1024):
        distances = cdist(query[low:low + 1024], centers, metric='euclidean')
        indices = np.argsort(distances, axis=1, kind='stable')[:, :min(neighbors, len(centers))]
        near = np.take_along_axis(distances, indices, axis=1)
        weights = 1 / np.maximum(near, 1e-6)
        weights /= weights.sum(axis=1, keepdims=True)
        chunks.append((outcomes[indices] * weights[..., None]).sum(axis=1))
    return np.vstack(chunks) if chunks else np.empty((0, outcomes.shape[1]))

def independent_scores(observed, predicted, totals, minimum):
    keep = np.isfinite(observed) & np.isfinite(predicted) & np.isfinite(totals)
    n = int(keep.sum())
    if n < minimum:return {'status': 'insufficient_evaluable_query_cells', 'n_cells': n}
    y, pred, depth = observed[keep].astype(float), predicted[keep].astype(float), np.log1p(totals[keep])
    depth -= depth.mean()
    def residual(v):
        centered = v - v.mean()
        norm = depth @ depth
        return centered - depth * (depth @ centered / norm) if norm > 0 else centered
    ry, rp = residual(y), residual(pred)
    out = {'status': 'estimated', 'n_cells': n,
           'MAE_native_logCP10K': float(np.abs(y - pred).mean()),
           'RMSE_native_logCP10K': float(np.sqrt(np.square(y - pred).mean()))}
    if ry.std() < 1e-8 or rp.std() < 1e-8:
        out['depth_adjusted_structure_status'] = 'no_variable_observed_or_predicted_signal'
    else:
        out.update(depth_adjusted_structure_status='estimated',
                   depth_adjusted_Pearson=float(np.corrcoef(ry, rp)[0, 1]),
                   depth_adjusted_Spearman=float(np.corrcoef(rankdata(ry), rankdata(rp))[0, 1]))
    return out

def check_score_bounds(row):
    if row['n_cells'] < 0 or row['n_cells'] > row['n_query_route_cells']:
        raise ValueError('Query score denominator outside scope')
    if row['status'] not in {'estimated', 'insufficient_evaluable_query_cells'}:
        raise ValueError('Unexpected score status')
    for field in ['depth_adjusted_Pearson', 'depth_adjusted_Spearman']:
        if field in row and (not np.isfinite(row[field]) or not -1 <= row[field] <= 1):
            raise ValueError('Invalid structure correlation')
    for field in ['MAE_native_logCP10K', 'RMSE_native_logCP10K']:
        if field in row and (not np.isfinite(row[field]) or row[field] < 0):
            raise ValueError('Invalid reconstruction error')
