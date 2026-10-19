"""Fixed literature-informed RNA candidates, no fitted annotation or outcomes."""
import hashlib
import numpy as np

MARKERS = ['CD3D', 'CD3E', 'CD3G', 'CD247', 'TRDC', 'TRGC1', 'TRGC2', 'TRDV1', 'TRDV2']
RULES = ['core5_count1', 'supported3_count1', 'core5_count2', 'supported3_count2']


def flags(counts, measured=None):
    assert counts.ndim == 2 and counts.shape[1] == len(MARKERS)
    measured = np.ones(len(MARKERS), bool) if measured is None else np.asarray(measured, bool)
    assert measured.shape == (len(MARKERS),)
    bits = np.zeros(len(counts), np.uint8)
    for k, minimum in enumerate([1, 2]):
        present = counts >= minimum
        core = present[:, :5].all(1) & measured[:5].all()
        supported = (present[:, :4].sum(1) >= 3) & present[:, 4] & present[:, 5:7].any(1) & measured[:7].all()
        bits |= core.astype(np.uint8) << (2 * k)
        bits |= supported.astype(np.uint8) << (2 * k + 1)
    return bits


def digest_array(array, dtype):
    return hashlib.sha256(np.ascontiguousarray(array, dtype=dtype).tobytes()).hexdigest()


def fixture():
    a = np.array([[1, 1, 1, 1, 1, 0, 0, 0, 0],
                  [1, 1, 1, 0, 1, 1, 0, 0, 0],
                  [0, 0, 0, 1, 1, 1, 1, 0, 0],
                  [2, 2, 2, 2, 2, 0, 0, 0, 0],
                  [2, 2, 2, 0, 2, 2, 0, 0, 0],
                  [1, 1, 1, 1, 0, 1, 1, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0]], np.int32)
    np.testing.assert_array_equal(flags(a), [1, 2, 0, 5, 10, 0, 0])
    missing = np.ones(9, bool); missing[4] = False
    assert not flags(a, missing).any()
    missing = np.ones(9, bool); missing[6] = False
    np.testing.assert_array_equal(flags(a, missing), [1, 0, 0, 5, 0, 0, 0])
    assert flags(a[:0]).shape == (0,)
    return dict(status='PASS_FIXED_CORE_DROPOUT_NKLIKE_ZERO_MISSING_AND_COUNT_THRESHOLDS', examples=len(a))
