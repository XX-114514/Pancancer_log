"""Exact blocked Euclidean neighbors without an optional SciPy dependency.

This computes full distances in bounded working blocks. It has quadratic
compute cost, and supplies only the query methods used by the fixed native
kNN definition. It is not a KD tree or approximate neighbor implementation.
"""
import resource
import sys
import numpy as np


def peak_rss_bytes():
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == 'darwin' else value*1024


class ExactEuclideanNeighbors:
    def __init__(self, points):
        self.points = np.asarray(points, dtype=float)
        assert self.points.ndim == 2 and self.points.shape[1] == 2
        assert len(self.points) > 0 and np.isfinite(self.points).all()
        self.block = min(128, max(1, 64*1024**2//(32*len(self.points))))

    def _distances(self, queries):
        with np.errstate(over='raise', invalid='raise'):
            x = queries[:, None, 0]-self.points[None, :, 0]
            y = queries[:, None, 1]-self.points[None, :, 1]
            np.square(x, out=x)
            np.square(y, out=y)
            x += y
            np.sqrt(x, out=x)
        assert np.isfinite(x).all()
        return x

    def _queries(self, values):
        values = np.asarray(values, dtype=float)
        assert values.ndim == 2 and values.shape[1] == 2 and np.isfinite(values).all()
        return values

    def query(self, values, k, workers=1):
        assert workers == 1 and isinstance(k, int) and 1 <= k <= len(self.points)
        values = self._queries(values)
        distances = np.empty((len(values), k), dtype=float)
        indices = np.empty((len(values), k), dtype=np.int64)
        for start in range(0, len(values), self.block):
            stop = min(start+self.block, len(values))
            full = self._distances(values[start:stop])
            selected = np.argpartition(full, k-1, axis=1)[:, :k]
            order = np.argsort(np.take_along_axis(full, selected, axis=1), axis=1, kind='stable')
            selected = np.take_along_axis(selected, order, axis=1)
            distances[start:stop] = np.take_along_axis(full, selected, axis=1)
            indices[start:stop] = selected
        return (distances[:, 0], indices[:, 0]) if k == 1 else (distances, indices)

    def query_ball_point(self, values, radius, workers=1):
        assert workers == 1
        values = self._queries(values)
        radius = np.broadcast_to(np.asarray(radius, dtype=float), (len(values),))
        assert np.isfinite(radius).all() and (radius >= 0).all()
        result = []
        for start in range(0, len(values), self.block):
            stop = min(start+self.block, len(values))
            full = self._distances(values[start:stop])
            result.extend(np.flatnonzero(row <= limit).tolist()
                for row, limit in zip(full, radius[start:stop]))
        return result
