"""Bounded source-local mapping diagnostics, never final cell identities."""
import numpy as np


def weighted_vote(distances, neighbor_labels):
    d = np.asarray(distances, float)
    labels = np.asarray(neighbor_labels)
    if (d.ndim != 2 or labels.shape != d.shape or not d.shape[1] or
            not np.isfinite(d).all() or (d < 0).any() or
            not np.issubdtype(labels.dtype, np.integer) or (labels < 0).any()):
        raise ValueError('Invalid neighbor distances or labels')
    prediction = np.empty(len(d), np.int32)
    confidence = np.empty(len(d), float)
    for i in range(len(d)):
        zero = d[i] <= 1e-12
        weights = zero.astype(float) if zero.any() else 1. / (d[i] + 1e-8)
        names, codes = np.unique(labels[i], return_inverse=True)
        totals = np.bincount(codes, weights=weights, minlength=len(names))
        winner = int(np.argmax(totals))  # smallest sorted group code on exact ties
        prediction[i] = names[winner]
        confidence[i] = totals[winner] / totals.sum()
    return prediction, confidence


def without_self(distances, indices, self_indices, k):
    d, ix, own = np.asarray(distances, float), np.asarray(indices), np.asarray(self_indices)
    if d.ndim != 2 or ix.shape != d.shape or own.shape != (len(d),) or k < 1:
        raise ValueError('Invalid self-exclusion axes')
    if not np.isfinite(d).all() or (d < 0).any():
        raise ValueError('Invalid self-exclusion distances')
    result = np.empty((len(d), k), float)
    for i in range(len(d)):
        other = d[i][ix[i] != own[i]]
        if len(other) < k:
            raise ValueError('Insufficient non-self neighbors')
        result[i] = other[:k]
    return result


def support_mask(mean_distance, vote_fraction, cutoff, minimum_vote):
    d, v = np.asarray(mean_distance, float), np.asarray(vote_fraction, float)
    if (d.shape != v.shape or not np.isfinite(cutoff) or cutoff < 0 or
            not 0 <= minimum_vote <= 1 or not np.isfinite(d).all() or
            not np.isfinite(v).all() or (d < 0).any() or ((v < 0) | (v > 1)).any()):
        raise ValueError('Invalid support inputs')
    return (d <= cutoff) & (v >= minimum_vote)
