#! /usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils.six import StringIO
from unittest import TestCase, main

import numpy as np
import pandas as pd
from pandas.util.testing import assert_frame_equal


from skbio.core.distance import DissimilarityMatrix, DistanceMatrix
from skbio.math.stats.distance.base import (CategoricalStats,
                                            CategoricalStatsResults, bioenv,
                                            _scale)
from skbio.util.testing import get_data_path


class CategoricalStatsTests(TestCase):
    def setUp(self):
        self.dm = DistanceMatrix([[0.0, 1.0, 2.0], [1.0, 0.0, 3.0],
                                  [2.0, 3.0, 0.0]], ['a', 'b', 'c'])
        self.grouping = [1, 2, 1]
        # Ordering of IDs shouldn't matter, nor should extra IDs.
        self.df = pd.read_csv(
            StringIO('ID,Group\nb,Group1\na,Group2\nc,Group1\nd,Group3'),
            index_col=0)
        self.df_missing_id = pd.read_csv(
            StringIO('ID,Group\nb,Group1\nc,Group1'), index_col=0)
        self.categorical_stats = CategoricalStats(self.dm, self.grouping)
        self.categorical_stats_from_df = CategoricalStats(self.dm, self.df,
                                                          column='Group')

    def test_init_invalid_input(self):
        # Requires a DistanceMatrix.
        with self.assertRaises(TypeError):
            CategoricalStats(DissimilarityMatrix([[0, 2], [3, 0]], ['a', 'b']),
                             [1, 2])

        # Requires column if DataFrame.
        with self.assertRaises(ValueError):
            CategoricalStats(self.dm, self.df)

        # Cannot provide column if not data frame.
        with self.assertRaises(ValueError):
            CategoricalStats(self.dm, self.grouping, column='Group')

        # Column must exist in data frame.
        with self.assertRaises(ValueError):
            CategoricalStats(self.dm, self.df, column='foo')

        # All distance matrix IDs must be in data frame.
        with self.assertRaises(ValueError):
            CategoricalStats(self.dm, self.df_missing_id, column='Group')

        # Grouping vector length must match number of objects in dm.
        with self.assertRaises(ValueError):
            CategoricalStats(self.dm, [1, 2])

        # Grouping vector cannot have only unique values.
        with self.assertRaises(ValueError):
            CategoricalStats(self.dm, [1, 2, 3])

        # Grouping vector cannot have only a single group.
        with self.assertRaises(ValueError):
            CategoricalStats(self.dm, [1, 1, 1])

    def test_call(self):
        with self.assertRaises(NotImplementedError):
            self.categorical_stats()

    def test_call_invalid_permutations(self):
        with self.assertRaises(ValueError):
            self.categorical_stats(-1)


class CategoricalStatsResultsTests(TestCase):

    def setUp(self):
        self.results = CategoricalStatsResults('foo', 'Foo', 'my stat', 42,
                                               ['a', 'b', 'c', 'd'],
                                               0.01234567890, 0.1151111, 99)
        self.p_value = 0.119123123123

    def test_str(self):
        exp = ('Method name  Sample size  Number of groups       my stat  '
               'p-value  Number of permutations\n        foo           42'
               '                 4  0.0123456789     0.12'
               '                      99\n')
        obs = str(self.results)
        self.assertEqual(obs, exp)

    def test_repr_html(self):
        # Not going to test against exact HTML that we expect, as this could
        # easily break and be annoying to constantly update. Do some light
        # sanity-checking to ensure there are some of the expected HTML tags.
        obs = self.results._repr_html_()
        self.assertTrue('<table' in obs)
        self.assertTrue('<thead' in obs)
        self.assertTrue('<tr' in obs)
        self.assertTrue('<th' in obs)
        self.assertTrue('<tbody' in obs)
        self.assertTrue('<td' in obs)
        self.assertTrue('1 rows' in obs)
        self.assertTrue('5 columns' in obs)

    def test_summary(self):
        exp = ('Method name\tSample size\tNumber of groups\tmy stat\tp-value\t'
               'Number of permutations\nfoo\t42\t4\t0.0123456789\t0.12\t99\n')
        obs = self.results.summary()
        self.assertEqual(obs, exp)

    def test_format_p_value(self):
        obs = self.results._format_p_value(self.p_value, 100)
        self.assertEqual(obs, '0.12')

        obs = self.results._format_p_value(self.p_value, 250)
        self.assertEqual(obs, '0.12')

        obs = self.results._format_p_value(self.p_value, 1000)
        self.assertEqual(obs, '0.119')

    def test_format_p_value_few_perms(self):
        obs = self.results._format_p_value(self.p_value, 9)
        self.assertEqual(obs, 'Too few permutations to compute p-value '
                              '(permutations = 9)')

        obs = self.results._format_p_value(self.p_value, 1)
        self.assertEqual(obs, 'Too few permutations to compute p-value '
                              '(permutations = 1)')

        obs = self.results._format_p_value(self.p_value, 0)
        self.assertEqual(obs, 'Too few permutations to compute p-value '
                              '(permutations = 0)')

    def test_format_p_value_none(self):
        obs = self.results._format_p_value(None, 0)
        self.assertEqual(obs, 'N/A')


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

        self.dm = DistanceMatrix.from_file(get_data_path('dm.txt'))

        # Reordered rows and columns (i.e., different ID order). Still
        # conceptually the same distance matrix.
        self.dm_reordered = DistanceMatrix.from_file(
            get_data_path('dm_reordered.txt'))

        self.df = pd.read_csv(get_data_path('df.txt'), sep='\t',
                                       index_col=0)

        # Similar to the above data frame, except that it has an extra
        # non-numeric column, and some of the other rows and columns have been
        # reordered.
        self.df_extra_column = pd.read_csv(
            get_data_path('df_extra_column.txt'), sep='\t', index_col=0)

        # All columns in the original data frame (these are all numeric
        # columns).
        self.cols = self.df.columns.tolist()

        self.exp_results = pd.read_csv(get_data_path('exp_results.txt'),
                                       sep='\t', index_col=0)
        self.exp_results_single_column = pd.read_csv(
            get_data_path('exp_results_single_column.txt'), sep='\t',
            index_col=0)
        self.exp_results_different_column_order = pd.read_csv(
            get_data_path('exp_results_different_column_order.txt'), sep='\t',
            index_col=0)

    def test_bioenv_all_columns_implicit(self):
        # Test with all columns in data frame (implicitly).
        obs = bioenv(self.dm, self.df)
        assert_frame_equal(obs, self.exp_results)

        # Should get the same results if order of rows/cols in distance matrix
        # is changed.
        obs = bioenv(self.dm_reordered, self.df)
        assert_frame_equal(obs, self.exp_results)

    def test_bioenv_all_columns_explicit(self):
        # Test with all columns being specified.
        obs = bioenv(self.dm, self.df, columns=self.cols)
        assert_frame_equal(obs, self.exp_results)

        # Test against a data frame that has an extra non-numeric column and
        # some of the rows and columns reordered (we should get the same
        # result since we're specifying the same columns in the same order).
        obs = bioenv(self.dm, self.df_extra_column, columns=self.cols)
        assert_frame_equal(obs, self.exp_results)

    def test_bioenv_single_column(self):
        obs = bioenv(self.dm, self.df, columns=['PH'])
        assert_frame_equal(obs, self.exp_results_single_column)

    def test_bioenv_different_column_order(self):
        # Specifying columns in a different order will change the row labels in
        # the results data frame as the column subsets will be reordered, but
        # the actual results (e.g., correlation coefficients) shouldn't change.
        obs = bioenv(self.dm, self.df, columns=self.cols[::-1])
        assert_frame_equal(obs, self.exp_results_different_column_order)

    def test_bioenv_no_side_effects(self):
        # Deep copies of both primary inputs.
        dm_copy = self.dm.copy()
        df_copy = self.df.copy(deep=True)

        bioenv(self.dm, self.df)

        # Make sure we haven't modified the primary input in some way (e.g.,
        # with scaling, type conversions, etc.).
        self.assertEqual(self.dm, dm_copy)
        assert_frame_equal(self.df, df_copy)

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

    # TODO add a few more tests for scale
    def test_scale(self):
        df = pd.DataFrame([[7.0, 400], [8.0, 530], [7.5, 450], [8.5, 810]],
                          index=['A','B','C','D'], columns=['pH', 'Elevation'])
        exp = pd.DataFrame([[-1.161895, -0.805979], [0.387298, -0.095625],
                            [-0.387298, -0.532766], [1.161895, 1.434369]],
                           index=['A','B','C','D'],
                           columns=['pH', 'Elevation'])
        obs = _scale(df)
        assert_frame_equal(obs, exp)

    def test_scale_single_column(self):
        df = pd.DataFrame([[1], [0], [2]], index=['A','B','C'],
                          columns=['foo'])
        exp = pd.DataFrame([[0.0], [-1.0], [1.0]], index=['A','B','C'],
                           columns=['foo'])
        obs = _scale(df)
        assert_frame_equal(obs, exp)


if __name__ == '__main__':
    main()
