"""Exact membership partition for candidate-subset expression evidence."""
import numpy as np


def partition_guarded_subsets(parent, celltype_code, eligible):
    parent, celltype_code, eligible = map(np.asarray, (parent, celltype_code, eligible))
    if parent.ndim != 1 or celltype_code.shape != parent.shape or eligible.shape != parent.shape:
        raise ValueError('Guarded subset axes differ')
    if eligible.dtype != bool or not np.issubdtype(parent.dtype,np.integer) or not np.issubdtype(celltype_code.dtype,np.integer):
        raise ValueError('Guarded subset domain/code dtypes differ')
    if np.any(parent < -1) or np.any(celltype_code < 0) or not np.array_equal(parent >= 0,eligible):
        raise ValueError('Guarded subset archive exclusions or candidate codes invalid')
    pairs, inverse = np.unique(np.column_stack([parent[eligible],celltype_code[eligible]]),axis=0,return_inverse=True)
    code = np.full(len(parent),-1,np.int32)
    code[eligible] = inverse
    return code,pairs
