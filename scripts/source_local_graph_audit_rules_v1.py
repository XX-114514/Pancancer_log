"""Array invariants for source-local discovery, independent of its producer."""
import numpy as np


def validate_partition(native_rows, groups, fine, eligible, families, n_groups):
    rows, groups, fine, eligible, families = map(np.asarray, (native_rows, groups, fine, eligible, families))
    if rows.ndim != 1 or any(a.shape != rows.shape for a in [groups, fine, eligible, families]):
        raise ValueError('Source graph row axes differ')
    if len(np.unique(rows)) != len(rows) or not np.issubdtype(rows.dtype, np.integer):
        raise ValueError('Native rows must be unique integer indices')
    if not np.issubdtype(groups.dtype, np.integer) or (groups < -1).any():
        raise ValueError('Invalid local group codes')
    if fine.dtype != bool or eligible.dtype != bool or np.any(fine & ~eligible):
        raise ValueError('Fine graph flag exceeds Boolean fit domain')
    if not np.array_equal(groups >= 0, eligible):
        raise ValueError('Fit exclusions or archive retention changed')
    if not np.array_equal(np.unique(groups[eligible]), np.arange(n_groups)):
        raise ValueError('Source group codes are not contiguous')
    for code in range(n_groups):
        mask = groups == code
        if len(np.unique(families[mask])) != 1 or len(np.unique(fine[mask])) != 1:
            raise ValueError('Local group crosses source branch or fine status')
    return True


def validate_model_axes(model, positions, native_rows, group_codes, columns):
    for name, expected in [('source_row_positions', positions), ('native_rows', native_rows),
                           ('source_local_group_code', group_codes), ('native_columns', columns)]:
        if not np.array_equal(model[name], expected):
            raise ValueError('Saved local model axis differs: ' + name)
    n, p = len(positions), len(columns)
    c, e, m, v = (model[k] for k in ['components', 'fit_embedding', 'mean', 'explained_variance'])
    if c.ndim != 2 or c.shape[1] != p or e.shape != (n, c.shape[0]) or m.shape != (p,) or v.shape != (c.shape[0],):
        raise ValueError('Saved PCA dimensions differ')
    if any(not np.isfinite(x).all() for x in [c, e, m, v]) or np.any(v < 0):
        raise ValueError('Saved PCA has invalid numeric values')
    np.testing.assert_allclose(c @ c.T, np.eye(len(c)), rtol=0, atol=2e-5)
    np.testing.assert_allclose(np.square(e.astype(np.float64)).sum(axis=0) / (n - 1), v, rtol=2e-5, atol=1e-6)
    return True
