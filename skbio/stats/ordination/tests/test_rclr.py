# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np
import numpy.testing as npt

from skbio.stats.composition import clr, closure
from skbio.stats.ordination._rclr import matrix_rclr, tensor_rclr, _matrix_closure


class TestMatrixClosure(unittest.TestCase):
    """Tests for the matrix closure function."""

    def test_basic_closure(self):
        """Test basic normalization to sum 1."""
        mat = np.array([[1, 2, 3], [4, 5, 6]])
        result = _matrix_closure(mat)

        # Each row should sum to 1
        npt.assert_almost_equal(result.sum(axis=1), [1.0, 1.0])

    def test_zero_row(self):
        """Test handling of zero rows."""
        mat = np.array([[0, 0, 0], [1, 2, 3]])
        result = _matrix_closure(mat)

        # First row should remain zeros
        npt.assert_almost_equal(result[0], [0.0, 0.0, 0.0])
        # Second row should sum to 1
        npt.assert_almost_equal(result[1].sum(), 1.0)

    def test_1d_input(self):
        """Test that 1D input is handled correctly."""
        vec = np.array([1, 2, 3])
        result = _matrix_closure(vec)

        self.assertEqual(result.ndim, 2)
        npt.assert_almost_equal(result.sum(), 1.0)


class TestMatrixRclr(unittest.TestCase):
    """Tests for the robust centered log-ratio transformation."""

    def test_basic_rclr(self):
        """Test basic rclr transformation."""
        # Simple matrix with no zeros
        mat = np.array([[1, 2, 3], [4, 5, 6]])
        result = matrix_rclr(mat)

        # For each row, the mean of transformed values should be ~0
        # (only over observed values)
        for i in range(mat.shape[0]):
            observed = ~np.isnan(result[i])
            npt.assert_almost_equal(result[i, observed].mean(), 0.0)

    def test_rclr_with_zeros(self):
        """Test rclr handles zeros by producing NaN."""
        mat = np.array([[1, 0, 3], [4, 5, 0]])
        result = matrix_rclr(mat)

        # Zeros should become NaN
        self.assertTrue(np.isnan(result[0, 1]))
        self.assertTrue(np.isnan(result[1, 2]))

        # Non-zero positions should not be NaN
        self.assertFalse(np.isnan(result[0, 0]))
        self.assertFalse(np.isnan(result[0, 2]))

    def test_rclr_centering(self):
        """Test that rclr centers each row correctly."""
        mat = np.array([[1, 2, 0, 4], [0, 3, 3, 0], [2, 2, 2, 2]])
        result = matrix_rclr(mat)

        # For rows with observed values, mean should be 0
        for i in range(mat.shape[0]):
            observed = ~np.isnan(result[i])
            if np.any(observed):
                npt.assert_almost_equal(result[i, observed].mean(), 0.0)

    def test_rclr_negative_values_error(self):
        """Test that negative values raise an error."""
        mat = np.array([[1, -2, 3], [4, 5, 6]])

        with self.assertRaises(ValueError) as context:
            matrix_rclr(mat)

        self.assertIn("negative", str(context.exception))

    def test_rclr_inf_error(self):
        """Test that infinite values raise an error."""
        mat = np.array([[1, np.inf, 3], [4, 5, 6]])

        with self.assertRaises(ValueError) as context:
            matrix_rclr(mat)

        self.assertIn("infinite", str(context.exception))

    def test_rclr_handles_nan(self):
        """Test that NaN values (missing entries) are handled."""
        mat = np.array([[1, np.nan, 3], [4, 5, 6]])
        result = matrix_rclr(mat)

        # NaN should remain NaN
        self.assertTrue(np.isnan(result[0, 1]))

        # Non-NaN values should be transformed
        self.assertFalse(np.isnan(result[0, 0]))
        self.assertFalse(np.isnan(result[0, 2]))

    def test_rclr_1d_input(self):
        """Test that 1D input is promoted to 2D."""
        vec = np.array([1, 2, 3])
        result = matrix_rclr(vec)

        self.assertEqual(result.ndim, 2)

    def test_rclr_uniform_row(self):
        """Test rclr on uniform row (all same values)."""
        mat = np.array([[2, 2, 2, 2]])
        result = matrix_rclr(mat)

        # Uniform row should have all zeros (log-ratio of equal values)
        npt.assert_almost_equal(result[0], [0.0, 0.0, 0.0, 0.0])

    def test_rclr_preserves_ratios(self):
        """Test that rclr preserves log-ratios between features."""
        mat = np.array([[1, 2, 4]])
        result = matrix_rclr(mat)

        # log(2) - log(1) = log(2)
        expected_ratio = np.log(2)
        observed_ratio = result[0, 1] - result[0, 0]
        npt.assert_almost_equal(observed_ratio, expected_ratio)


class TestTensorRclr(unittest.TestCase):
    """Tests for tensor rclr transformation."""

    def test_basic_tensor_rclr(self):
        """Test basic 3D tensor rclr."""
        tensor = np.array([
            [[1, 2, 3], [4, 5, 6]],
            [[7, 8, 9], [10, 11, 12]]
        ])
        result = tensor_rclr(tensor)

        # Shape should be preserved
        self.assertEqual(result.shape, tensor.shape)

        # Each sample's mean should be approximately 0
        result_2d = result.reshape(-1, tensor.shape[-1])
        for i in range(result_2d.shape[0]):
            observed = ~np.isnan(result_2d[i])
            if np.any(observed):
                npt.assert_almost_equal(result_2d[i, observed].mean(), 0.0,
                                         decimal=5)

    def test_tensor_rclr_with_zeros(self):
        """Test tensor rclr handles zeros."""
        tensor = np.array([
            [[1, 0, 3], [0, 5, 6]],
            [[7, 8, 0], [10, 0, 12]]
        ])
        result = tensor_rclr(tensor)

        # Zeros should become NaN
        self.assertTrue(np.isnan(result[0, 0, 1]))
        self.assertTrue(np.isnan(result[0, 1, 0]))
        self.assertTrue(np.isnan(result[1, 0, 2]))
        self.assertTrue(np.isnan(result[1, 1, 1]))

    def test_tensor_rclr_preserves_shape(self):
        """Test that various tensor shapes are preserved."""
        shapes = [(2, 3, 4), (5, 2, 6), (3, 3, 3)]

        for shape in shapes:
            tensor = np.random.rand(*shape) + 0.1  # Avoid zeros
            result = tensor_rclr(tensor)
            self.assertEqual(result.shape, shape)

    def test_tensor_rclr_negative_error(self):
        """Test that negative values raise an error."""
        tensor = np.array([[[1, -2, 3]]])

        with self.assertRaises(ValueError) as context:
            tensor_rclr(tensor)

        self.assertIn("negative", str(context.exception))


class TestRclrGemelliBehavior(unittest.TestCase):
    """Tests verifying gemelli-compatible behavior for rclr.

    These tests use specific numeric values from the gemelli test suite
    to ensure the scikit-bio implementation produces identical results.
    """

    def test_rclr_dense_equals_clr(self):
        """Test that matrix_rclr equals clr on dense data without zeros.

        Based on gemelli's test_rclr_dense test which validates that
        the rclr transformation produces identical results to scikit-bio's
        CLR transformation when there are no zeros (i.e., when no robust
        handling is needed).
        """
        # Dense count data with no zeros
        count_data = np.array([[2, 2, 6], [4, 4, 2]])

        # Apply matrix_rclr
        rclr_result = matrix_rclr(count_data)

        # Apply closure then clr from skbio.stats.composition
        closed_data = closure(count_data)
        clr_result = clr(closed_data)

        # Results should be identical for dense data
        npt.assert_allclose(rclr_result, clr_result, rtol=1e-10)

    def test_rclr_sparse_with_expected_values(self):
        """Test rclr on sparse data with known expected values.

        Based on gemelli's test_rclr_sparse test with specific expected
        numeric values.
        """
        # Sparse count data with zeros
        count_data = np.array([[3, 3, 0], [0, 4, 2]])

        # Expected values from gemelli (zeros become NaN, observed values centered)
        expected = np.array([[0.0, 0.0, np.nan],
                            [np.nan, 0.34657359, -0.34657359]])

        result = matrix_rclr(count_data)

        # Check non-NaN values match expected
        npt.assert_allclose(result[~np.isnan(result)],
                           expected[~np.isnan(expected)],
                           rtol=1e-5)

        # Check NaN positions match
        npt.assert_array_equal(np.isnan(result), np.isnan(expected))

    def test_closure_with_zeros(self):
        """Test closure with rows containing all zeros.

        Based on gemelli's test_closure_missing test.
        """
        # Table with one row of all zeros
        table = np.array([[2, 2, 6], [0, 0, 0]])

        # Expected: first row normalized, second row becomes NaN
        # (our implementation keeps zeros for zero rows)
        result = _matrix_closure(table)

        # First row should sum to 1
        npt.assert_almost_equal(result[0].sum(), 1.0)
        npt.assert_almost_equal(result[0], [0.2, 0.2, 0.6])

        # Second row should remain zeros (our implementation)
        # Note: gemelli returns NaN, but our impl keeps zeros for zero rows
        npt.assert_almost_equal(result[1], [0.0, 0.0, 0.0])

    def test_closure_without_zeros(self):
        """Test closure normalization without zeros.

        Based on gemelli's test_closure test.
        """
        table = np.array([[2, 2, 6], [2, 2, 6]])
        expected = np.array([[0.2, 0.2, 0.6], [0.2, 0.2, 0.6]])

        result = _matrix_closure(table)

        npt.assert_allclose(result, expected)

    def test_matrix_tensor_rclr_consistency(self):
        """Test that tensor_rclr produces same results as matrix_rclr.

        Based on gemelli's test_matrix_tensor_rclr test which validates
        that the tensor and matrix implementations are consistent.
        """
        # Create a 3D tensor that can be compared element-wise with 2D
        tensor_3d = np.array([
            [[1., 2., 3.], [4., 5., 6.]],
            [[7., 8., 9.], [10., 11., 12.]],
            [[13., 14., 15.], [16., 17., 18.]]
        ])

        # Apply tensor_rclr to 3D
        tensor_result = tensor_rclr(tensor_3d)

        # Apply matrix_rclr to equivalent 2D (reshape and apply)
        tensor_2d = tensor_3d.reshape(-1, tensor_3d.shape[-1])
        matrix_result = matrix_rclr(tensor_2d)
        matrix_result_reshaped = matrix_result.reshape(tensor_3d.shape)

        # Results should be identical
        npt.assert_allclose(tensor_result, matrix_result_reshaped, rtol=1e-10)

    def test_rclr_nan_input_raises(self):
        """Test that NaN values in input are handled appropriately.

        Note: Our implementation allows NaN as missing entries,
        but gemelli raises an error. We test our behavior.
        """
        # Input with NaN (representing missing data)
        data_with_nan = np.array([[1, np.nan, 3], [4, 5, 6]])

        # Our implementation should handle this gracefully
        result = matrix_rclr(data_with_nan)

        # NaN should remain NaN
        self.assertTrue(np.isnan(result[0, 1]))

        # Non-NaN values should be transformed
        self.assertFalse(np.isnan(result[0, 0]))
        self.assertFalse(np.isnan(result[0, 2]))


class TestRclrNumericalAccuracy(unittest.TestCase):
    """Tests for numerical accuracy of rclr transformation."""

    def test_rclr_log_ratio_preservation(self):
        """Test that log-ratios are preserved exactly."""
        # Data where we can compute expected values analytically
        data = np.array([[1, 2, 4, 8]])

        result = matrix_rclr(data)

        # Log-ratio between adjacent pairs should be log(2)
        expected_ratio = np.log(2)
        for i in range(3):
            observed_ratio = result[0, i + 1] - result[0, i]
            npt.assert_almost_equal(observed_ratio, expected_ratio)

    def test_rclr_geometric_mean_centering(self):
        """Test that each row is centered at geometric mean."""
        np.random.seed(42)
        data = np.random.rand(5, 10) * 100 + 1  # No zeros

        result = matrix_rclr(data)

        # Each row should have mean of 0 (centered at geometric mean)
        for i in range(data.shape[0]):
            row_mean = np.nanmean(result[i])
            npt.assert_almost_equal(row_mean, 0.0, decimal=10)


if __name__ == '__main__':
    unittest.main()
