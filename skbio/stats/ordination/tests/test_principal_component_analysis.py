# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np

from skbio import OrdinationResults
from skbio.stats.ordination import pca
from skbio.util import assert_ordination_results_equal
from skbio.table._tabular import _create_table, _create_table_1d


class TestPCA(TestCase):
    def setUp(self):
        # Simple dataset with 4 samples and 4 features
        self.X = np.array([
            [1, 4, 2, 6],
            [2, 5, 1, 7],
            [1, 5, 2, 4],
            [3, 4, 3, 6]
        ])

        # Full dimensionality expected results

        components = np.array([
            [0.48260815, -0.12451597, -0.06637759, 0.86439526],
            [0.51789148, -0.33553222, 0.73499837, -0.28104082],
            [0.63156772, 0.6869599, -0.23531131, -0.27172945],
            [-0.31622777, 0.63245553, 0.63245553, 0.31622777]
            ])
        
        variances = np.array([1.95862455, 1.18116419, 0.36021126, 0])

        variance_proportion = np.array([0.55960701, 0.33747548, 0.1029175, 0])

        projected_samples = np.array([
            [-0.08359932, -0.29091271, -0.8850881, 0],
            [1.20526571, -1.12459263, 0.39702138, 0],
            [-1.9369058, -0.06436329, 0.3453307, 0],
            [0.8152394, 1.47986863, 0.14273602, 0]
            ])

        row_ids = ["Sample 1", "Sample 2", "Sample 3", "Sample 4"]
        column_ids = ["Feature 1", "Feature 2", "Feature 3", "Feature 4"]
        pc_ids = ["PC1", "PC2", "PC3", "PC4"]
        
        eigvals = _create_table_1d(
            variances, index=pc_ids, backend=None
            )
        samples = _create_table(
            projected_samples, index=row_ids, columns=pc_ids, backend=None
            )
        features = _create_table(
            components, index=pc_ids, columns=column_ids, backend=None
            )
        proportion_explained = _create_table_1d(
            variance_proportion, index=pc_ids, backend=None
            )

        self.full_dim_expected_results = OrdinationResults(
            short_method_name="PCA",
            long_method_name="Principal Component Analysis",
            eigvals=eigvals,
            samples=samples,
            sample_ids=row_ids,
            features=features,
            feature_ids=pc_ids,
            proportion_explained=proportion_explained,
        )

        # Reduced dimensionality expected results (2D)

        components = np.array([
            [0.48260815, -0.12451597, -0.06637759, 0.86439526],
            [0.51789148, -0.33553222, 0.73499837, -0.28104082]
            ])
        
        variances = np.array([1.95862455, 1.18116419])

        variance_proportion = np.array([0.55960701, 0.33747548])

        projected_samples = np.array([
            [-0.08359932, -0.29091271],
            [1.20526571, -1.12459263],
            [-1.9369058, -0.06436329],
            [0.8152394, 1.47986863]
            ])

        row_ids = ["Sample 1", "Sample 2", "Sample 3", "Sample 4"]
        column_ids = ["Feature 1", "Feature 2", "Feature 3", "Feature 4"]
        pc_ids = ["PC1", "PC2"]
        
        eigvals = _create_table_1d(
            variances, index=pc_ids, backend=None
            )
        samples = _create_table(
            projected_samples, index=row_ids, columns=pc_ids, backend=None
            )
        features = _create_table(
            components, index=pc_ids, columns=column_ids, backend=None
            )
        proportion_explained = _create_table_1d(
            variance_proportion, index=pc_ids, backend=None
            )

        self.reduced_dim_expected_results = OrdinationResults(
            short_method_name="PCA",
            long_method_name="Principal Component Analysis",
            eigvals=eigvals,
            samples=samples,
            sample_ids=row_ids,
            features=features,
            feature_ids=pc_ids,
            proportion_explained=proportion_explained,
        )

    def test_simple(self):
        results = pca(self.X, 
                      sample_ids=["Sample 1", "Sample 2", "Sample 3", "Sample 4"],
                      feature_ids=["Feature 1", "Feature 2", "Feature 3", "Feature 4"])

        assert_ordination_results_equal(results, self.full_dim_expected_results,
                                        ignore_directionality=True)
        
    def test_dimensionality_reduction(self):
        results = pca(self.X, 
                      dimensions=2,
                      sample_ids=["Sample 1", "Sample 2", "Sample 3", "Sample 4"],
                      feature_ids=["Feature 1", "Feature 2", "Feature 3", "Feature 4"])

        assert_ordination_results_equal(results, self.reduced_dim_expected_results,
                                        ignore_directionality=True)
        
    def test_svd_full_dim(self):
        results = pca(self.X, 
                      method="svd",
                      sample_ids=["Sample 1", "Sample 2", "Sample 3", "Sample 4"],
                      feature_ids=["Feature 1", "Feature 2", "Feature 3", "Feature 4"])

        assert_ordination_results_equal(results, self.full_dim_expected_results,
                                        ignore_directionality=True)
        
    def test_svd_reduced_dim(self):
        results = pca(self.X, 
                      method="svd",
                      dimensions=2,
                      sample_ids=["Sample 1", "Sample 2", "Sample 3", "Sample 4"],
                      feature_ids=["Feature 1", "Feature 2", "Feature 3", "Feature 4"])

        assert_ordination_results_equal(results, self.reduced_dim_expected_results,
                                        ignore_directionality=True)
        
    def test_svd_full_dim_iterative(self):
        results = pca(self.X, 
                      method="svd",
                      iterative=True,
                      sample_ids=["Sample 1", "Sample 2", "Sample 3", "Sample 4"],
                      feature_ids=["Feature 1", "Feature 2", "Feature 3", "Feature 4"])

        assert_ordination_results_equal(results, self.full_dim_expected_results,
                                        ignore_directionality=True)
        
    def test_svd_reduced_dim_iterative(self):
        results = pca(self.X, 
                      method="svd",
                      dimensions=2,
                      iterative=True,
                      sample_ids=["Sample 1", "Sample 2", "Sample 3", "Sample 4"],
                      feature_ids=["Feature 1", "Feature 2", "Feature 3", "Feature 4"])

        assert_ordination_results_equal(results, self.reduced_dim_expected_results,
                                        ignore_directionality=True)
        
    def test_eigh_full_dim(self):
        results = pca(self.X, 
                      method="eigh",
                      sample_ids=["Sample 1", "Sample 2", "Sample 3", "Sample 4"],
                      feature_ids=["Feature 1", "Feature 2", "Feature 3", "Feature 4"])

        assert_ordination_results_equal(results, self.full_dim_expected_results,
                                        ignore_directionality=True)
        
    def test_eigh_reduced_dim(self):
        results = pca(self.X, 
                      method="eigh",
                      dimensions=2,
                      sample_ids=["Sample 1", "Sample 2", "Sample 3", "Sample 4"],
                      feature_ids=["Feature 1", "Feature 2", "Feature 3", "Feature 4"])

        assert_ordination_results_equal(results, self.reduced_dim_expected_results,
                                        ignore_directionality=True)
        
    def test_eigh_full_dim(self):
        results = pca(self.X, 
                      method="eigh",
                      iterative=True,
                      sample_ids=["Sample 1", "Sample 2", "Sample 3", "Sample 4"],
                      feature_ids=["Feature 1", "Feature 2", "Feature 3", "Feature 4"])

        assert_ordination_results_equal(results, self.full_dim_expected_results,
                                        ignore_directionality=True)
        
    def test_eigh_reduced_dim(self):
        results = pca(self.X, 
                      method="eigh",
                      dimensions=2,
                      iterative=True,
                      sample_ids=["Sample 1", "Sample 2", "Sample 3", "Sample 4"],
                      feature_ids=["Feature 1", "Feature 2", "Feature 3", "Feature 4"])

        assert_ordination_results_equal(results, self.reduced_dim_expected_results,
                                        ignore_directionality=True)
        
    def test_invalid_method(self):
        with self.assertRaises(ValueError):
            pca(self.X, method="invalid_method")

    def test_invalid_dimensions(self):
        with self.assertRaises(ValueError):
            pca(self.X, dimensions=0)
        with self.assertRaises(ValueError):
            pca(self.X, dimensions=0.1)
        with self.assertRaises(ValueError):
            pca(self.X, dimensions=5)

if __name__ == "__main__":
    main()
