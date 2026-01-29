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
import pandas as pd

from skbio import OrdinationResults
from skbio.stats.ordination._rpca import rpca, _filter_table


class TestFilterTable(unittest.TestCase):
    """Tests for table filtering function."""

    def setUp(self):
        """Set up test fixtures."""
        self.table = pd.DataFrame(
            [[10, 0, 5, 2],
             [20, 3, 0, 8],
             [5, 1, 1, 1],
             [0, 0, 0, 1],
             [15, 2, 3, 4]],
            index=['s1', 's2', 's3', 's4', 's5'],
            columns=['f1', 'f2', 'f3', 'f4']
        )

    def test_no_filtering(self):
        """Test that no filtering preserves table."""
        result = _filter_table(self.table)
        pd.testing.assert_frame_equal(result, self.table)

    def test_filter_by_sample_count(self):
        """Test filtering by minimum sample count."""
        result = _filter_table(self.table, min_sample_count=10)

        # s4 has sum=1, s3 has sum=8 - both should be removed
        self.assertNotIn('s4', result.index)
        self.assertNotIn('s3', result.index)
        self.assertEqual(len(result), 3)

    def test_filter_by_feature_count(self):
        """Test filtering by minimum feature count."""
        result = _filter_table(self.table, min_feature_count=10)

        # Features with sum >= 10: f1=50, f4=16
        self.assertIn('f1', result.columns)
        self.assertIn('f4', result.columns)
        self.assertEqual(len(result.columns), 2)

    def test_filter_by_feature_frequency(self):
        """Test filtering by minimum feature frequency."""
        # f1 appears in 4/5 = 0.8, f2 in 3/5 = 0.6, f3 in 3/5 = 0.6, f4 in 5/5 = 1.0
        result = _filter_table(self.table, min_feature_frequency=0.7)

        # Only f1 (4/5) and f4 (5/5) meet the 0.7 threshold
        self.assertIn('f1', result.columns)
        self.assertIn('f4', result.columns)
        self.assertNotIn('f2', result.columns)
        self.assertNotIn('f3', result.columns)
        self.assertEqual(len(result.columns), 2)


class TestRPCA(unittest.TestCase):
    """Tests for Robust PCA ordination."""

    def setUp(self):
        """Set up test fixtures."""
        np.random.seed(42)

        # Create a simple count table
        n_samples, n_features = 15, 20

        # Generate counts with some zeros
        counts = np.random.poisson(5, size=(n_samples, n_features))
        counts[counts < 2] = 0

        self.table = pd.DataFrame(
            counts,
            index=['sample_%d' % i for i in range(n_samples)],
            columns=['feature_%d' % i for i in range(n_features)]
        )

    def test_basic_rpca(self):
        """Test basic RPCA analysis."""
        ordination = rpca(self.table, n_components=3)

        # Check ordination results
        self.assertIsInstance(ordination, OrdinationResults)
        self.assertEqual(ordination.short_method_name, 'RPCA')

        # Check samples shape
        self.assertEqual(ordination.samples.shape[0], self.table.shape[0])
        self.assertEqual(ordination.samples.shape[1], 3)

        # Check features shape
        self.assertEqual(ordination.features.shape[0], self.table.shape[1])
        self.assertEqual(ordination.features.shape[1], 3)

    def test_rpca_preserves_sample_ids(self):
        """Test that sample IDs are preserved."""
        ordination = rpca(self.table, n_components=2)

        self.assertListEqual(
            list(ordination.samples.index),
            list(self.table.index)
        )

    def test_rpca_preserves_feature_ids(self):
        """Test that feature IDs are preserved."""
        ordination = rpca(self.table, n_components=2)

        self.assertListEqual(
            list(ordination.features.index),
            list(self.table.columns)
        )

    def test_rpca_proportion_explained(self):
        """Test that proportion explained sums to <= 1."""
        ordination = rpca(self.table, n_components=3)

        # Proportion should sum to approximately 1 or less
        total_prop = ordination.proportion_explained.sum()
        self.assertLessEqual(total_prop, 1.01)  # Allow small numerical error

        # Each proportion should be non-negative
        self.assertTrue(all(ordination.proportion_explained >= 0))

    def test_rpca_eigvals_decreasing(self):
        """Test that eigenvalues are in decreasing order."""
        ordination = rpca(self.table, n_components=3)

        eigvals = ordination.eigvals.values
        for i in range(len(eigvals) - 1):
            self.assertGreaterEqual(eigvals[i], eigvals[i + 1])

    def test_rpca_with_filtering(self):
        """Test RPCA with filtering parameters."""
        ordination = rpca(
            self.table,
            n_components=2,
            min_sample_count=5,
            min_feature_count=5,
            min_feature_frequency=0.1
        )

        # Results should still be valid
        self.assertIsInstance(ordination, OrdinationResults)

    def test_rpca_non_dataframe_error(self):
        """Test error on non-DataFrame input."""
        with self.assertRaises(ValueError) as context:
            rpca(self.table.values, n_components=2)

        self.assertIn("DataFrame", str(context.exception))

    def test_rpca_negative_values_error(self):
        """Test error on negative values."""
        table_neg = self.table.copy()
        table_neg.iloc[0, 0] = -5

        with self.assertRaises(ValueError) as context:
            rpca(table_neg, n_components=2)

        self.assertIn("negative", str(context.exception))

    def test_rpca_insufficient_samples_error(self):
        """Test error when too few samples after filtering."""
        small_table = self.table.iloc[:2, :]

        with self.assertRaises(ValueError) as context:
            rpca(small_table, n_components=2)

        self.assertIn("samples", str(context.exception))

    def test_rpca_insufficient_features_error(self):
        """Test error when n_components exceeds features."""
        small_table = self.table.iloc[:, :2]

        with self.assertRaises(ValueError) as context:
            rpca(small_table, n_components=5)

        self.assertIn("features", str(context.exception))


class TestRPCAReproducibility(unittest.TestCase):
    """Tests for RPCA reproducibility."""

    def test_reproducible_with_seed(self):
        """Test that results are reproducible with same random seed."""
        np.random.seed(42)
        counts = np.random.poisson(5, size=(10, 15))
        counts[counts < 2] = 0

        table = pd.DataFrame(
            counts,
            index=['s%d' % i for i in range(10)],
            columns=['f%d' % i for i in range(15)]
        )

        # Run twice with same seed
        np.random.seed(123)
        ord1 = rpca(table, n_components=2)

        np.random.seed(123)
        ord2 = rpca(table, n_components=2)

        # Results should be identical
        npt.assert_almost_equal(
            ord1.samples.values, ord2.samples.values
        )


class TestRPCAAxisLabels(unittest.TestCase):
    """Tests for RPCA axis labeling."""

    def test_axis_labels(self):
        """Test that axis labels are correct."""
        np.random.seed(42)
        counts = np.random.poisson(5, size=(10, 15))
        counts[counts < 2] = 0

        table = pd.DataFrame(
            counts,
            index=['s%d' % i for i in range(10)],
            columns=['f%d' % i for i in range(15)]
        )

        ordination = rpca(table, n_components=3)

        expected_labels = ['PC1', 'PC2', 'PC3']
        self.assertListEqual(list(ordination.eigvals.index), expected_labels)
        self.assertListEqual(list(ordination.samples.columns), expected_labels)
        self.assertListEqual(list(ordination.features.columns), expected_labels)
        self.assertListEqual(
            list(ordination.proportion_explained.index), expected_labels
        )


class TestRPCAIntegration(unittest.TestCase):
    """Integration tests for RPCA based on gemelli tutorials.

    These tests verify that RPCA can reproduce the expected behavior
    demonstrated in the gemelli tutorials, ensuring the implementation
    is functionally equivalent.
    """

    def test_rpca_separates_groups(self):
        """Test that RPCA can separate known sample groups.

        Based on the gemelli RPCA tutorial where samples from different
        body sites should separate in the ordination.
        """
        np.random.seed(42)

        # Simulate microbiome data with two distinct groups
        n_samples_per_group = 15
        n_features = 50

        # Group 1: dominated by features 0-24
        group1_counts = np.random.poisson(8, size=(n_samples_per_group, n_features))
        group1_counts[:, :25] = np.random.poisson(20, size=(n_samples_per_group, 25))
        group1_counts[group1_counts < 2] = 0

        # Group 2: dominated by features 25-49
        group2_counts = np.random.poisson(8, size=(n_samples_per_group, n_features))
        group2_counts[:, 25:] = np.random.poisson(20, size=(n_samples_per_group, 25))
        group2_counts[group2_counts < 2] = 0

        # Combine into single table
        all_counts = np.vstack([group1_counts, group2_counts])
        sample_ids = ['g1_s%d' % i for i in range(n_samples_per_group)] + \
                     ['g2_s%d' % i for i in range(n_samples_per_group)]
        feature_ids = ['f%d' % i for i in range(n_features)]

        table = pd.DataFrame(all_counts, index=sample_ids, columns=feature_ids)

        # Run RPCA
        ordination = rpca(table, n_components=3)

        # Extract PC1 coordinates for each group
        g1_pc1 = ordination.samples.loc[
            [s for s in sample_ids if s.startswith('g1')], 'PC1'
        ].values
        g2_pc1 = ordination.samples.loc[
            [s for s in sample_ids if s.startswith('g2')], 'PC1'
        ].values

        # Groups should be separated - their means should be significantly different
        g1_mean = np.mean(g1_pc1)
        g2_mean = np.mean(g2_pc1)

        # Check separation exists
        self.assertNotAlmostEqual(g1_mean, g2_mean, places=1)

    def test_rpca_biplot_loadings(self):
        """Test that feature loadings in biplot are informative.

        Based on gemelli tutorial showing that feature loadings
        indicate which taxa drive sample separation.
        """
        np.random.seed(123)

        # Create data where specific features drive separation
        n_samples = 20
        n_features = 30

        # Features 0-4 are high in first half of samples
        # Features 25-29 are high in second half
        counts = np.random.poisson(5, size=(n_samples, n_features))

        # Make first 5 features dominant in first half
        counts[:10, :5] = np.random.poisson(50, size=(10, 5))
        # Make last 5 features dominant in second half
        counts[10:, 25:] = np.random.poisson(50, size=(10, 5))

        counts[counts < 2] = 0

        table = pd.DataFrame(
            counts,
            index=['s%d' % i for i in range(n_samples)],
            columns=['f%d' % i for i in range(n_features)]
        )

        ordination = rpca(table, n_components=2)

        # Feature loadings should indicate which features drive separation
        # Features 0-4 and 25-29 should have larger loadings
        feature_loadings = ordination.features

        # Check that feature loadings exist and have correct shape
        self.assertEqual(feature_loadings.shape[0], n_features)
        self.assertEqual(feature_loadings.shape[1], 2)

        # Dominant features should have larger absolute loadings on PC1
        dominant_features = list(range(5)) + list(range(25, 30))
        minor_features = list(range(10, 20))

        pc1_loadings = np.abs(feature_loadings['PC1'].values)

        dominant_mean = np.mean(pc1_loadings[dominant_features])
        minor_mean = np.mean(pc1_loadings[minor_features])

        # Dominant features should have higher loadings
        self.assertGreater(dominant_mean, minor_mean)

    def test_rpca_sparse_data_handling(self):
        """Test that RPCA handles highly sparse microbiome data.

        Microbiome data typically has 70-90% zeros. RPCA should
        handle this gracefully.
        """
        np.random.seed(456)

        n_samples, n_features = 25, 100

        # Create very sparse data (80% zeros)
        counts = np.random.poisson(2, size=(n_samples, n_features))
        counts[np.random.rand(n_samples, n_features) < 0.8] = 0

        # Ensure minimum observations per sample/feature
        counts[:, :10] = np.random.poisson(5, size=(n_samples, 10))
        counts[:10, :] = counts[:10, :] + np.random.poisson(3, size=(10, n_features))

        table = pd.DataFrame(
            counts,
            index=['s%d' % i for i in range(n_samples)],
            columns=['f%d' % i for i in range(n_features)]
        )

        # Verify data is sparse
        sparsity = (table.values == 0).sum() / table.size
        self.assertGreater(sparsity, 0.3)  # At least 30% sparse

        # RPCA should complete successfully
        ordination = rpca(table, n_components=3)

        # Results should be valid
        self.assertFalse(np.any(np.isnan(ordination.samples.values)))
        self.assertFalse(np.any(np.isinf(ordination.samples.values)))

    def test_rpca_distance_computation(self):
        """Test that RPCA coordinates can be used for distance matrices.

        Based on gemelli workflow where RPCA coordinates are used
        with PERMANOVA statistical testing.
        """
        np.random.seed(789)

        n_samples, n_features = 20, 40
        counts = np.random.poisson(10, size=(n_samples, n_features))
        counts[counts < 2] = 0

        table = pd.DataFrame(
            counts,
            index=['s%d' % i for i in range(n_samples)],
            columns=['f%d' % i for i in range(n_features)]
        )

        ordination = rpca(table, n_components=3)

        # Compute Euclidean distances in PC space
        from scipy.spatial.distance import pdist, squareform

        coords = ordination.samples.values
        distances = pdist(coords, metric='euclidean')
        dm = squareform(distances)

        # Distance matrix should be symmetric
        npt.assert_almost_equal(dm, dm.T)

        # Diagonal should be zero
        npt.assert_almost_equal(np.diag(dm), np.zeros(n_samples))

        # Off-diagonal should be positive
        self.assertTrue(np.all(dm[np.triu_indices(n_samples, k=1)] > 0))

    def test_rpca_min_sample_count_filtering(self):
        """Test that min_sample_count filtering works as in gemelli.

        The gemelli tutorial uses min_sample_count=500 to filter
        low-depth samples.
        """
        np.random.seed(101)

        n_samples, n_features = 15, 30

        # Create samples with varying depths
        counts = np.zeros((n_samples, n_features))
        depths = [1000, 2000, 500, 100, 50, 1500, 800, 200, 1200, 300,
                  600, 1100, 150, 900, 400]

        for i, depth in enumerate(depths):
            probs = np.random.dirichlet(np.ones(n_features))
            counts[i, :] = np.random.multinomial(depth, probs)

        counts = counts.astype(int)

        table = pd.DataFrame(
            counts,
            index=['s%d' % i for i in range(n_samples)],
            columns=['f%d' % i for i in range(n_features)]
        )

        # Run with min_sample_count=500
        ordination = rpca(table, n_components=2, min_sample_count=500)

        # Check that only samples with >= 500 reads are included
        n_expected = sum(1 for d in depths if d >= 500)
        self.assertEqual(ordination.samples.shape[0], n_expected)


class TestRPCAGemelliBehavior(unittest.TestCase):
    """Tests verifying gemelli-compatible behavior."""

    def test_rclr_transformation(self):
        """Test that RPCA applies rclr transformation correctly.

        The rclr transformation is central to gemelli's approach.
        """
        from skbio.stats.ordination._rclr import matrix_rclr

        # Simple test case from gemelli
        data = np.array([[1, 2, 3], [4, 5, 6], [0, 8, 9]])
        result = matrix_rclr(data)

        # Zeros should become NaN
        self.assertTrue(np.isnan(result[2, 0]))

        # Non-zero rows should be centered (mean of observed values = 0)
        for i in range(3):
            observed = ~np.isnan(result[i])
            if np.any(observed):
                npt.assert_almost_equal(result[i, observed].mean(), 0.0)

    def test_proportion_explained_calculation(self):
        """Test proportion explained calculation matches expected behavior."""
        np.random.seed(42)

        n_samples, n_features = 20, 30
        counts = np.random.poisson(10, size=(n_samples, n_features))
        counts[counts < 2] = 0

        table = pd.DataFrame(
            counts,
            index=['s%d' % i for i in range(n_samples)],
            columns=['f%d' % i for i in range(n_features)]
        )

        ordination = rpca(table, n_components=5)

        # Proportion explained should sum to <= 1
        self.assertLessEqual(ordination.proportion_explained.sum(), 1.01)

        # Each proportion should be non-negative
        self.assertTrue(all(ordination.proportion_explained >= 0))

        # Proportions should be in decreasing order
        props = ordination.proportion_explained.values
        for i in range(len(props) - 1):
            self.assertGreaterEqual(props[i], props[i + 1])


if __name__ == '__main__':
    unittest.main()
