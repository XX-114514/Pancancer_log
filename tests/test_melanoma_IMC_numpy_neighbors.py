"""Independent scalar-distance, boundary-tie, block and memory-unit checks."""
import importlib.util
import math
from pathlib import Path
import sys
import unittest
import numpy as np


class ExactNeighbors(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[1]/'scripts'/'melanoma_IMC_numpy_neighbors_v1.py'
        spec = importlib.util.spec_from_file_location('numpy_neighbors', path)
        cls.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.m)

    def test_full_scalar_distance_replay_across_blocks(self):
        rng = np.random.default_rng(20261021)
        points = rng.integers(-100, 101, (173, 2)).astype(float)
        queries = rng.integers(-100, 101, (137, 2)).astype(float)
        tree = self.m.ExactEuclideanNeighbors(points)
        tree.block = 11
        actual, indices = tree.query(queries, 11)
        for query, distance, selected in zip(queries, actual, indices):
            scalar = np.array([math.sqrt((query[0]-x)**2+(query[1]-y)**2) for x, y in points])
            np.testing.assert_allclose(distance, np.sort(scalar)[:11], rtol=0, atol=0)
            np.testing.assert_allclose(distance, scalar[selected], rtol=0, atol=0)
        radii = np.nextafter(actual[:, -1], np.inf)
        balls = tree.query_ball_point(queries, radii)
        for query, radius, neighbors in zip(queries, radii, balls):
            independent = [i for i, point in enumerate(points) if math.dist(query, point) <= radius]
            self.assertEqual(neighbors, independent)

    def test_all_boundary_ties_and_coincident_points_retained(self):
        points = np.array([[0., 0.], [0., 0.], [1., 0.], [-1., 0.], [0., 1.], [0., -1.]])
        tree = self.m.ExactEuclideanNeighbors(points)
        distance, _ = tree.query(points, 3)
        balls = tree.query_ball_point(points, np.nextafter(distance[:, -1], np.inf))
        self.assertEqual(balls[0], list(range(6)))
        self.assertEqual(balls[1], list(range(6)))

    def test_queries_do_not_stitch_or_mutate_input(self):
        points = np.array([[0., 0.], [1., 0.], [10., 0.]])
        before = points.copy()
        tree = self.m.ExactEuclideanNeighbors(points)
        distance, _ = tree.query(points, 2)
        np.testing.assert_array_equal(points, before)
        scaled = self.m.ExactEuclideanNeighbors(points*7)
        changed, _ = scaled.query(points*7, 2)
        np.testing.assert_array_equal(changed, distance*7)

    def test_platform_memory_units_and_failure_gates(self):
        raw = int(self.m.resource.getrusage(self.m.resource.RUSAGE_SELF).ru_maxrss)
        expected = raw if sys.platform == 'darwin' else raw*1024
        self.assertEqual(self.m.peak_rss_bytes(), expected)
        for values in [np.array([[np.nan, 0.]]), np.zeros((2, 3))]:
            with self.assertRaises(AssertionError):self.m.ExactEuclideanNeighbors(values)
        tree = self.m.ExactEuclideanNeighbors(np.zeros((2, 2)))
        with self.assertRaises(AssertionError):tree.query(np.zeros((1, 2)), 3)


if __name__ == '__main__':
    unittest.main()
