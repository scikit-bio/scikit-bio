# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np
import numpy.testing as npt

from skbio.stats.ordination._optspace import OptSpace, _svd_sort


class TestSVDSort(unittest.TestCase):
    """Tests for SVD sorting utility."""

    def test_sort_by_singular_values(self):
        """Test that SVD components are sorted by descending singular values."""
        # U and V need same number of columns as length of S
        U = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).T  # 3x3
        S = np.array([1, 3, 2])
        V = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]]).T  # 3x3

        U_sorted, S_sorted, V_sorted = _svd_sort(U, S, V)

        # S should be sorted in descending order
        npt.assert_array_equal(S_sorted, [3, 2, 1])

    def test_sort_preserves_correspondence(self):
        """Test that U, S, V correspondence is maintained after sorting."""
        np.random.seed(42)
        U = np.random.randn(5, 3)
        S = np.array([1, 3, 2])
        V = np.random.randn(4, 3)

        # Original reconstruction
        M_orig = U.dot(np.diag(S)).dot(V.T)

        U_sorted, S_sorted, V_sorted = _svd_sort(U, S, V)

        # Sorted reconstruction should match (up to sign)
        M_sorted = U_sorted.dot(np.diag(S_sorted)).dot(V_sorted.T)

        npt.assert_almost_equal(np.abs(M_orig), np.abs(M_sorted))


class TestOptSpace(unittest.TestCase):
    """Tests for OptSpace matrix completion algorithm."""

    def setUp(self):
        """Set up test fixtures."""
        np.random.seed(42)

        # Create a low-rank matrix
        self.n, self.m, self.r = 20, 15, 3
        self.U_true = np.random.randn(self.n, self.r)
        self.V_true = np.random.randn(self.m, self.r)
        self.M_true = self.U_true.dot(self.V_true.T)

        # Create partially observed matrix
        self.M_observed = self.M_true.copy()
        # Randomly mask 30% of entries
        mask = np.random.rand(self.n, self.m) < 0.3
        self.M_observed[mask] = np.nan

    def test_basic_completion(self):
        """Test basic matrix completion on low-rank matrix."""
        opt = OptSpace(n_components=self.r, max_iterations=10)
        M_recovered = opt.fit_transform(self.M_observed)

        # Check that observed entries are well recovered
        observed_mask = ~np.isnan(self.M_observed)
        rmse_observed = np.sqrt(np.mean(
            (M_recovered[observed_mask] - self.M_true[observed_mask]) ** 2
        ))

        # RMSE on observed entries should be small
        self.assertLess(rmse_observed, 1.0)

    def test_fit_then_transform(self):
        """Test fit and transform methods separately."""
        opt = OptSpace(n_components=self.r, max_iterations=10)

        # Fit should return self
        result = opt.fit(self.M_observed)
        self.assertIs(result, opt)

        # Transform should return matrix
        M_recovered = opt.transform()
        self.assertEqual(M_recovered.shape, (self.n, self.m))

    def test_get_loadings(self):
        """Test that loadings have correct shapes."""
        opt = OptSpace(n_components=self.r, max_iterations=10)
        opt.fit(self.M_observed)

        sample_loadings, feature_loadings = opt.get_loadings()

        self.assertEqual(sample_loadings.shape, (self.n, self.r))
        self.assertEqual(feature_loadings.shape, (self.m, self.r))

    def test_not_fitted_error(self):
        """Test that methods raise error when not fitted."""
        opt = OptSpace()

        with self.assertRaises(ValueError) as context:
            opt.transform()
        self.assertIn("fitted", str(context.exception))

        with self.assertRaises(ValueError) as context:
            opt.get_loadings()
        self.assertIn("fitted", str(context.exception))

    def test_invalid_dimensions(self):
        """Test error on invalid n_components."""
        opt = OptSpace(n_components=100)  # Too large

        with self.assertRaises(ValueError) as context:
            opt.fit(self.M_observed)
        self.assertIn("cannot exceed", str(context.exception))

    def test_non_2d_error(self):
        """Test error on non-2D input."""
        opt = OptSpace()

        with self.assertRaises(ValueError) as context:
            opt.fit(np.random.randn(5, 5, 5))
        self.assertIn("2D", str(context.exception))

    def test_attributes_after_fit(self):
        """Test that attributes are set after fitting."""
        opt = OptSpace(n_components=self.r, max_iterations=10)
        opt.fit(self.M_observed)

        self.assertIsNotNone(opt.U)
        self.assertIsNotNone(opt.S)
        self.assertIsNotNone(opt.V)

        self.assertEqual(opt.U.shape, (self.n, self.r))
        self.assertEqual(opt.S.shape, (self.r, self.r))
        self.assertEqual(opt.V.shape, (self.m, self.r))

    def test_reproducibility(self):
        """Test that results are reproducible with same seed."""
        np.random.seed(123)
        opt1 = OptSpace(n_components=self.r, max_iterations=5)
        M1 = opt1.fit_transform(self.M_observed)

        np.random.seed(123)
        opt2 = OptSpace(n_components=self.r, max_iterations=5)
        M2 = opt2.fit_transform(self.M_observed)

        npt.assert_almost_equal(M1, M2)

    def test_sparse_matrix(self):
        """Test completion on very sparse matrix."""
        # Create a matrix with 70% missing
        M_sparse = self.M_true.copy()
        mask = np.random.rand(self.n, self.m) < 0.7
        M_sparse[mask] = np.nan

        opt = OptSpace(n_components=self.r, max_iterations=20)
        M_recovered = opt.fit_transform(M_sparse)

        # Check that it returns something reasonable
        self.assertEqual(M_recovered.shape, (self.n, self.m))
        self.assertFalse(np.any(np.isnan(M_recovered)))


class TestOptSpaceNoisyData(unittest.TestCase):
    """Tests for OptSpace with noisy data."""

    def test_noisy_completion(self):
        """Test completion on noisy low-rank matrix."""
        np.random.seed(42)

        n, m, r = 30, 25, 2
        U_true = np.random.randn(n, r)
        V_true = np.random.randn(m, r)
        M_true = U_true.dot(V_true.T)

        # Add noise
        noise_level = 0.1
        M_noisy = M_true + noise_level * np.random.randn(n, m)

        # Mask some entries
        M_observed = M_noisy.copy()
        mask = np.random.rand(n, m) < 0.2
        M_observed[mask] = np.nan

        opt = OptSpace(n_components=r, max_iterations=15)
        M_recovered = opt.fit_transform(M_observed)

        # The reconstruction should be closer to the noisy observed
        # values than to random
        observed_mask = ~np.isnan(M_observed)
        rmse = np.sqrt(np.mean(
            (M_recovered[observed_mask] - M_noisy[observed_mask]) ** 2
        ))

        # RMSE should be reasonable (not perfect due to noise)
        self.assertLess(rmse, 1.0)


if __name__ == '__main__':
    unittest.main()
