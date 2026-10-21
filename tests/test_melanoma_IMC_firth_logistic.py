"""Finite differences, exact intercept correction, separation and gates."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np


class FirthLogistic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[1]/'scripts'/'melanoma_IMC_firth_logistic_v1.py'
        spec = importlib.util.spec_from_file_location('firth_under_test', path)
        cls.m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.m)

    def test_closed_form_intercept_all_binary_outcomes(self):
        for count in [0, 3, 9, 10]:
            outcome = np.r_[np.ones(count), np.zeros(10-count)]
            fit = self.m.fit_firth(np.ones((10, 1)), outcome)
            np.testing.assert_allclose(fit['probabilities'], (count+.5)/11, atol=1e-7)

    def test_penalty_gradient_independent_central_difference(self):
        rng = np.random.default_rng(4)
        x = np.column_stack([np.ones(24), rng.normal(size=(24, 3))])
        y = rng.integers(0, 2, size=24)
        beta = np.array([.3, -.4, .2, -.1])
        _, gradient = self.m.objective(beta, x, y)
        numerical = []
        for j in range(len(beta)):
            delta = np.zeros(len(beta)); delta[j] = 1e-5
            numerical.append((self.m.objective(beta+delta, x, y)[0]-self.m.objective(beta-delta, x, y)[0])/2e-5)
        np.testing.assert_allclose(gradient, numerical, rtol=1e-6, atol=1e-7)

    def test_complete_separation_finite(self):
        x = np.column_stack([np.ones(24), np.linspace(-2, 2, 24)])
        y = (x[:, 1] > 0).astype(int)
        fit = self.m.fit_firth(x, y)
        self.assertGreater(fit['coefficients'][1], 0)
        self.assertTrue(((fit['probabilities'] > 0) & (fit['probabilities'] < 1)).all())

    def test_covariate_reordering_preserves_predictions(self):
        rng = np.random.default_rng(13)
        x = np.column_stack([np.ones(24), rng.normal(size=(24, 3))])
        y = rng.integers(0, 2, 24)
        first = self.m.fit_firth(x, y)
        second = self.m.fit_firth(x[:, [0, 3, 1, 2]], y)
        np.testing.assert_allclose(first['probabilities'], second['probabilities'], atol=1e-7)

    def test_no_silent_covariate_removal(self):
        values = np.column_stack([np.arange(24), np.ones(24)])
        with self.assertRaisesRegex(AssertionError, 'constant'):
            self.m.standardized_design(values)
        values[:, 1] = values[:, 0]
        with self.assertRaisesRegex(AssertionError, 'rank_deficient'):
            self.m.standardized_design(values)
        with self.assertRaisesRegex(AssertionError, 'insufficient'):
            self.m.standardized_design(np.arange(19)[:, None])


if __name__ == '__main__':
    unittest.main()
