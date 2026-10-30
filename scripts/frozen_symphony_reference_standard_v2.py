"""Train-only standard centered-ARPACK Harmony/Symphony references.

The independently hashed optional wheel is loaded only for this process. No
shared environment is changed. Minimal references contain no expression matrix
or cell metadata. Mapping adapts query batches but never refits the reference.
"""
import hashlib
import importlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

from standard_annotation_order_primitives_v1 import pca_fit


def load_symphonypy(wheel, expected_sha256):
    path = Path(wheel)
    if hashlib.sha256(path.read_bytes()).hexdigest() != expected_sha256:
        raise ValueError('Reference mapping dependency changed')
    if str(path) not in sys.path:sys.path.insert(0, str(path))
    module = importlib.import_module('symphonypy')
    if not str(module.__file__).startswith(str(path) + '/'):
        raise ValueError('Unexpected reference mapping dependency origin')
    return module


def valid_log_input(x, genes):
    if x.ndim != 2 or x.shape[1] != len(genes) or len(set(genes)) != len(genes):
        raise ValueError('Invalid reference feature axis')
    values = x.data if sp.issparse(x) else np.asarray(x)
    if (not np.isfinite(values).all() or (values < 0).any() or
            (values > np.log1p(10000) + 1e-5).any()):
        raise ValueError('Expected native log1p(CP10K), not corrected or scaled expression')


def reference_fingerprint(reference):
    digest = hashlib.sha256()
    digest.update('\0'.join(reference.var_names).encode())
    for x in [reference.var['mean'].to_numpy(), reference.var['std'].to_numpy(), reference.varm['PCs'],
              reference.uns['harmony']['C'], reference.uns['harmony']['Nr'], reference.uns['harmony']['sigma']]:
        a = np.asarray(x)
        digest.update(str(a.dtype).encode());digest.update(str(a.shape).encode());digest.update(a.tobytes())
    return digest.hexdigest()


def build_reference(x_train, source_train, genes, seed, parameters, symphony):
    """Caller must select HVGs on training sources before entering this function."""
    import anndata as ad
    valid_log_input(x_train, genes)
    if len(x_train.shape) != 2 or len(source_train) != x_train.shape[0]:
        raise ValueError('Reference/source axis differs')
    raw, pca = pca_fit(x_train, seed, parameters['components'], parameters['pca_iterations'])
    reference = ad.AnnData(x_train.copy(),
        obs=pd.DataFrame({'source': np.asarray(source_train).astype(str)},
                         index=['training_' + str(i) for i in range(len(source_train))]),
        var=pd.DataFrame({'mean': pca.mean_, 'std': np.ones(len(genes), dtype=np.float32)}, index=genes))
    reference.varm['PCs'] = pca.components_.T.copy()
    reference.obsm['X_pca'] = raw
    # The installed backend leaves scalar sigma as a Python float when K=1,
    # then calls .astype() on it. Explicit vectors preserve the same intended
    # defaults for every K and qualify the small-reference boundary case.
    soft_clusters = max(1, min(int(np.round(len(source_train) / 30.)), 100))
    symphony.pp.harmony_integrate(reference, key='source', random_seed=seed,
        theta=parameters['theta'], lamb=parameters['lambda'],
        nclust=soft_clusters, sigma=np.full(soft_clusters, .1),
        max_iter_harmony=parameters['max_iter_harmony'], ncores=parameters['threads'], verbose=False)
    corrected = np.asarray(reference.obsm['X_pca_harmony'])
    if corrected.shape != raw.shape or not np.isfinite(corrected).all():
        raise ValueError('Invalid reference correction')
    # All centered logCP10K values with std=1 lie inside (-10,10). Therefore
    # symphonypy's projection clip does not change this declared PCA model.
    small = ad.AnnData(sp.csr_matrix((0, len(genes)), dtype=np.float32), var=reference.var.copy())
    small.varm['PCs'] = reference.varm['PCs'].copy()
    h = reference.uns['harmony']
    small.uns['harmony'] = {k: np.array(h[k], copy=True) if isinstance(h[k], np.ndarray) else h[k]
                           for k in ['Nr', 'C', 'K', 'sigma', 'ref_basis_loadings', 'ref_basis_adjusted', 'converged']}
    return small, raw, corrected, pca


def map_query(reference, x_query, genes, source_query, symphony):
    import anndata as ad
    valid_log_input(x_query, genes)
    if list(genes) != reference.var_names.tolist() or len(source_query) != x_query.shape[0]:
        raise ValueError('Query genes/order differ; source-missing genes must not be imputed as zero')
    before = reference_fingerprint(reference)
    projected = np.asarray(x_query @ reference.varm['PCs']) - reference.var['mean'].to_numpy() @ reference.varm['PCs']
    # The upstream routine normalizes by each row's maximum before its L2 norm.
    # Preserve numerically unsupported rows instead of passing zero/negative-max
    # vectors through that normalization. These flags are not biological calls.
    supported = (np.isfinite(projected).all(axis=1) & (np.linalg.norm(projected, axis=1) > 1e-10) &
                 (np.max(projected, axis=1) > 0))
    mapped = np.full(projected.shape, np.nan, dtype=np.float32)
    if supported.any():
        query = ad.AnnData(x_query[supported].copy(),
            obs=pd.DataFrame({'source': np.asarray(source_query)[supported].astype(str)},
                             index=['query_' + str(i) for i in np.flatnonzero(supported)]),
            var=pd.DataFrame(index=genes))
        query.uns['log1p'] = {'base': None}
        symphony.tl.map_embedding(query, reference, key='source', lamb=1., use_genes_column=None)
        np.testing.assert_allclose(query.obsm['X_pca_reference'], projected[supported], rtol=2e-5, atol=2e-4)
        mapped[supported] = query.obsm['X_pca_harmony']
        if not np.isfinite(mapped[supported]).all():raise ValueError('Non-finite query mapping')
    if reference_fingerprint(reference) != before:
        raise ValueError('Query mapping changed the frozen reference')
    return projected, mapped, supported


def save_minimal_reference(path, reference, feature_candidate_indices, metadata):
    """Write private model parameters only; no counts or individual records."""
    h = reference.uns['harmony']
    with Path(path).open('xb') as f:
        np.savez_compressed(f, genes=np.asarray(reference.var_names, dtype=str), mean=reference.var['mean'].to_numpy(),
                            std=reference.var['std'].to_numpy(), loadings=reference.varm['PCs'],
                            C=h['C'], Nr=h['Nr'], sigma=h['sigma'],
                            feature_candidate_indices=np.asarray(feature_candidate_indices),
                            metadata_json=np.asarray(json.dumps(metadata, allow_nan=False)))
