# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import pandas as pd

from skbio import DistanceMatrix
from skbio.stats.distance import bioenv
from skbio.stats.distance._bioenv import _scale
from skbio.util import get_data_path, assert_data_frame_almost_equal
from skbio.util._testing import _data_frame_to_default_int_type


class BIOENVTests(TestCase):
    """Results were verified with R 3.0.2 and vegan 2.0-10 (vegan::bioenv)."""

    def setUp(self):
        # The test dataset used here is a subset of the Lauber et al. 2009
        # "88 Soils" dataset. It has been altered to exercise various aspects
        # of the code, including (but not limited to):
        #
        # - order of distance matrix IDs and IDs in data frame (metadata) are
        #   not exactly the same
        # - data frame has an extra sample that is not in the distance matrix
        # - this extra sample has non-numeric and missing values in some of its
        #   cells
        #
        # Additional variations of the distance matrix and data frame are used
        # to test different orderings of rows/columns, extra non-numeric data
        # frame columns, etc.
        #
        # This dataset is also useful because it is non-trivial in size (6
        # samples, 11 environment variables) and it includes positive/negative
        # floats and integers in the data frame.
        self.dm = DistanceMatrix.read(get_data_path('dm.txt'))

        # Reordered rows and columns (i.e., different ID order). Still
        # conceptually the same distance matrix.
        self.dm_reordered = DistanceMatrix.read(
            get_data_path('dm_reordered.txt'))

        self.df = pd.read_csv(get_data_path('df.txt'), sep='\t', index_col=0)

        # Similar to the above data frame, except that it has an extra
        # non-numeric column, and some of the other rows and columns have been
        # reordered.
        self.df_extra_column = pd.read_csv(
            get_data_path('df_extra_column.txt'), sep='\t', index_col=0)

        # All columns in the original data frame (these are all numeric
        # columns).
        self.cols = self.df.columns.tolist()

        # This second dataset is derived from vegan::bioenv's example dataset
        # (varespec and varechem). The original dataset includes a site x
        # species table (e.g., OTU table) and a data frame of environmental
        # variables. Since the bioenv function defined here accepts a distance
        # matrix, we use a Bray-Curtis distance matrix that is derived from the
        # site x species table (this matches what is done by vegan::bioenv when
        # provided an OTU table, using their default distance measure). The
        # data frame only includes the numeric environmental variables we're
        # interested in for these tests: log(N), P, K, Ca, pH, Al
        self.dm_vegan = DistanceMatrix.read(
            get_data_path('bioenv_dm_vegan.txt'))
        self.df_vegan = pd.read_csv(
            get_data_path('bioenv_df_vegan.txt'), sep='\t',
            converters={0: str})
        self.df_vegan.set_index('#SampleID', inplace=True)

        # Load expected results.
        self.exp_results = pd.read_csv(get_data_path('exp_results.txt'),
                                       sep='\t',
                                       index_col=0)
        _data_frame_to_default_int_type(self.exp_results)

        self.exp_results_single_column = pd.read_csv(
            get_data_path('exp_results_single_column.txt'),
            sep='\t',
            index_col=0
        )
        _data_frame_to_default_int_type(self.exp_results_single_column)

        self.exp_results_different_column_order = pd.read_csv(
            get_data_path('exp_results_different_column_order.txt'),
            sep='\t',
            index_col=0
        )
        _data_frame_to_default_int_type(self.exp_results_different_column_order)

        self.exp_results_vegan = pd.read_csv(
            get_data_path('bioenv_exp_results_vegan.txt'),
            sep='\t',
            index_col=0
        )
        _data_frame_to_default_int_type(self.exp_results_vegan)

    def test_bioenv_all_columns_implicit(self):
        # Test with all columns in data frame (implicitly).
        obs = bioenv(self.dm, self.df)
        assert_data_frame_almost_equal(obs, self.exp_results)

        # Should get the same results if order of rows/cols in distance matrix
        # is changed.
        obs = bioenv(self.dm_reordered, self.df)
        assert_data_frame_almost_equal(obs, self.exp_results)

    def test_bioenv_all_columns_explicit(self):
        # Test with all columns being specified.
        obs = bioenv(self.dm, self.df, columns=self.cols)
        assert_data_frame_almost_equal(obs, self.exp_results)

        # Test against a data frame that has an extra non-numeric column and
        # some of the rows and columns reordered (we should get the same
        # result since we're specifying the same columns in the same order).
        obs = bioenv(self.dm, self.df_extra_column, columns=self.cols)
        assert_data_frame_almost_equal(obs, self.exp_results)

    def test_bioenv_single_column(self):
        obs = bioenv(self.dm, self.df, columns=['PH'])
        assert_data_frame_almost_equal(obs, self.exp_results_single_column)

    def test_bioenv_different_column_order(self):
        # Specifying columns in a different order will change the row labels in
        # the results data frame as the column subsets will be reordered, but
        # the actual results (e.g., correlation coefficients) shouldn't change.
        obs = bioenv(self.dm, self.df, columns=self.cols[::-1])
        assert_data_frame_almost_equal(
            obs,
            self.exp_results_different_column_order)

    def test_bioenv_no_side_effects(self):
        # Deep copies of both primary inputs.
        dm_copy = self.dm.copy()
        df_copy = self.df.copy(deep=True)

        bioenv(self.dm, self.df)

        # Make sure we haven't modified the primary input in some way (e.g.,
        # with scaling, type conversions, etc.).
        self.assertEqual(self.dm, dm_copy)
        assert_data_frame_almost_equal(self.df, df_copy)

    def test_bioenv_vegan_example(self):
        # The correlation coefficient in the first row of the
        # results (rho=0.2516) is different from the correlation coefficient
        # computed by vegan (rho=0.2513). This seems to occur due to
        # differences in numerical precision when calculating the Euclidean
        # distances, which affects the rank calculations in Spearman
        # (specifically, dealing with ties). The ranked distances end up being
        # slightly different between vegan and our implementation because some
        # distances are treated as ties in vegan but treated as distinct values
        # in our implementation. This explains the difference in rho values. I
        # verified that using Pearson correlation instead of Spearman on the
        # same distances yields *very* similar results. Thus, the discrepancy
        # seems to stem from differences when computing ranks/ties.
        obs = bioenv(self.dm_vegan, self.df_vegan)
        assert_data_frame_almost_equal(
            obs,
            self.exp_results_vegan,
            rtol=1e-3
            )

    def test_bioenv_no_distance_matrix(self):
        with self.assertRaises(TypeError):
            bioenv('breh', self.df)

    def test_bioenv_no_data_frame(self):
        with self.assertRaises(TypeError):
            bioenv(self.dm, None)

    def test_bioenv_duplicate_columns(self):
        with self.assertRaises(ValueError):
            bioenv(self.dm, self.df, columns=self.cols + ['PH'])

    def test_bioenv_no_columns(self):
        with self.assertRaises(ValueError):
            bioenv(self.dm, self.df, columns=[])

    def test_bioenv_missing_columns(self):
        with self.assertRaises(ValueError):
            bioenv(self.dm, self.df, columns=self.cols + ['brofist'])

    def test_bioenv_missing_distance_matrix_ids(self):
        df = self.df[1:]
        with self.assertRaises(ValueError):
            bioenv(self.dm, df)

    def test_bioenv_nans(self):
        df = self.df.replace(53.9, np.nan)
        with self.assertRaises(ValueError):
            bioenv(self.dm, df)

    def test_bioenv_nonnumeric_columns(self):
        df = self.df.replace(2400, 'no cog yay')
        with self.assertRaises(TypeError):
            bioenv(self.dm, df)

        with self.assertRaises(TypeError):
            bioenv(self.dm, self.df_extra_column)

    def test_scale_single_column(self):
        df = pd.DataFrame([[1], [0], [2]], index=['A', 'B', 'C'],
                          columns=['foo'])
        exp = pd.DataFrame([[0.0], [-1.0], [1.0]], index=['A', 'B', 'C'],
                           columns=['foo'])
        obs = _scale(df)
        assert_data_frame_almost_equal(obs, exp)

    def test_scale_multiple_columns(self):
        # Floats and ints, including positives and negatives.
        df = pd.DataFrame([[7.0, 400, -1],
                           [8.0, 530, -5],
                           [7.5, 450, 1],
                           [8.5, 810, -4]],
                          index=['A', 'B', 'C', 'D'],
                          columns=['pH', 'Elevation', 'negatives'])
        exp = pd.DataFrame([[-1.161895, -0.805979, 0.453921],
                            [0.387298, -0.095625, -0.998625],
                            [-0.387298, -0.532766, 1.180194],
                            [1.161895, 1.434369, -0.635489]],
                           index=['A', 'B', 'C', 'D'],
                           columns=['pH', 'Elevation', 'negatives'])
        obs = _scale(df)
        assert_data_frame_almost_equal(obs, exp)

    def test_scale_no_variance(self):
        df = pd.DataFrame([[-7.0, -1.2], [6.2, -1.2], [2.9, -1.2]],
                          index=['A', 'B', 'C'], columns=['foo', 'bar'])
        with self.assertRaises(ValueError):
            _scale(df)


if __name__ == '__main__':
    main()
