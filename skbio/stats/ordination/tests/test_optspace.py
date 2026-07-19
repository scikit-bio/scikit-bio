# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest
from unittest.mock import patch

import numpy as np
import numpy.testing as npt

from skbio.stats.ordination._optspace import optspace, line_search


class TestOptSpace(unittest.TestCase):
    """Tests for OptSpace matrix completion algorithm."""

    def setUp(self):
        """Set up test fixtures."""
        rng = np.random.default_rng(0)
        self.tol = 1e-5

        # Create a low-rank matrix
        self.m, self.n, self.r = 100, 100, 3
        self.U_true = rng.normal(size=(self.m, self.r))
        self.V_true = rng.normal(size=(self.n, self.r))
        self.M_true = self.U_true @ self.V_true.T

        # Introduce noise
        noise_level = 0.01
        Z = noise_level * rng.normal(size=self.M_true.shape)
        M_noisy = self.M_true + Z

        # Mask some entries
        p = 0.4
        self.M_obs = M_noisy.copy()
        self.mask = rng.random((self.m, self.n)) < p
        self.M_obs[~self.mask] = np.nan

    def compute_rmse(self, M_hat):
        return np.sqrt(np.sum((M_hat - self.M_true)**2) / (self.m * self.n))

    def test_basic_completion(self):
        """Test basic matrix completion on low-rank matrix."""
        M_hat = optspace(self.M_obs, self.r, tol=self.tol)
        rmse = self.compute_rmse(M_hat)
        true_rmse = 0.0041953807728951865
        rel_error = (rmse - true_rmse) / true_rmse
        self.assertLessEqual(rel_error, self.tol)

    def test_gauss_newton(self):
        """Test low-rank matrix completion with Gauss-Newton"""
        M_hat = optspace(self.M_obs, self.r, method='GN', tol=self.tol)
        rmse = self.compute_rmse(M_hat)
        true_rmse = 0.004193496955751318
        rel_error = (rmse - true_rmse) / true_rmse
        self.assertLessEqual(rel_error, self.tol)

    def test_gradient_descent(self):
        """Test low-rank matrix completion with gradient descent"""
        M_hat = optspace(self.M_obs, self.r, method='LS', tol=self.tol)
        rmse = self.compute_rmse(M_hat)
        true_rmse = 0.0041953807728951865
        rel_error = (rmse - true_rmse) / true_rmse
        self.assertLessEqual(rel_error, self.tol)

    def test_positive_dimension(self):
        """Test error for nonpositive dimension."""
        dim = 0 # Not positive
        with self.assertRaises(ValueError) as context:
            _ = optspace(self.M_obs, dimensions=dim)
        self.assertIn("positive", str(context.exception))

    def test_exceed_dimension(self):
        """Test error for dimension exceeding array size."""
        dim = max(self.m, self.n) + 1 # Too large
        with self.assertRaises(ValueError) as context:
            _ = optspace(self.M_obs, dimensions=dim)
        self.assertIn("less than", str(context.exception))

    def test_integer_dimension(self):
        """Test error for noninteger dimension."""
        dim = 1.1 # Not an integer
        with self.assertRaises(ValueError) as context:
            _ = optspace(self.M_obs, dimensions=dim)
        self.assertIn("integer", str(context.exception))

    def test_non_2d_error(self):
        """Test error on non-2D input."""
        with self.assertRaises(ValueError) as context:
            _ = optspace(np.random.randn(5, 5, 5))
        self.assertIn("2D", str(context.exception))

    def test_invalid_method(self):
        """Test invalid method parameter."""
        with self.assertRaises(ValueError) as context:
            _ = optspace(self.M_obs, dimensions=self.r, method=1)
        self.assertIn("Method", str(context.exception))

    def test_svds_exception(self):
        with patch("skbio.stats.ordination._optspace.svds",
            side_effect=RuntimeError("PROPACK failed")
        ):
            with self.assertWarns(RuntimeWarning) as cm:
                _ = optspace(self.M_obs, self.r)

        self.assertIn("Sparse SVD failed", str(cm.warning))

    def test_max_iter(self):
        with self.assertWarns(RuntimeWarning) as cm:
            _ = optspace(self.M_obs, self.r, max_iter=1)

        self.assertIn("did not converge", str(cm.warning))

    def test_max_ls(self):
        with patch(
            "skbio.stats.ordination._optspace.retract_grassmann",
            side_effect=lambda X, dX: X
        ):
            with self.assertWarns(RuntimeWarning) as cm:
                _ = optspace(self.M_obs, self.r, method='LS')

        self.assertIn("depth limit was reached", str(cm.warning))

# TODO: More varied sizes, shapes, and missingness among matrices
# Vary total size (m*n), ratio (m/n), and portion p of observed entries


if __name__ == '__main__':
    unittest.main()