# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from functools import partial
from unittest import TestCase, main

import numpy as np
import pandas as pd
from pandas.testing import assert_series_equal
from scipy.spatial.distance import squareform

from skbio import DistanceMatrix
from skbio.stats.distance import permanova
from skbio.util import get_data_path
from skbio.stats.distance._base import _preprocess_input_sng


class PERMANOVATests(TestCase):
    """All results were verified with R (vegan::adonis)."""

    def setUp(self):
        # Distance matrices with and without ties in the ranks, with 2 groups
        # of equal size.
        dm_ids = ['s1', 's2', 's3', 's4']
        self.grouping_equal = ['Control', 'Control', 'Fast', 'Fast']
        self.df = pd.read_csv(
            io.StringIO('ID,Group\ns2,Control\ns3,Fast\ns4,Fast\ns5,Control\n'
                        's1,Control'), index_col=0)

        self.dm_ties = DistanceMatrix([[0, 1, 1, 4],
                                       [1, 0, 3, 2],
                                       [1, 3, 0, 3],
                                       [4, 2, 3, 0]], dm_ids)

        self.dm_no_ties = DistanceMatrix([[0, 1, 5, 4],
                                          [1, 0, 3, 2],
                                          [5, 3, 0, 3],
                                          [4, 2, 3, 0]], dm_ids)

        # Test with 3 groups of unequal size.
        self.grouping_unequal = ['Control', 'Treatment1', 'Treatment2',
                                 'Treatment1', 'Control', 'Control']

        # Equivalent grouping but with different labels -- groups should be
        # assigned different integer labels but results should be the same.
        self.grouping_unequal_relabeled = ['z', 42, 'abc', 42, 'z', 'z']

        self.dm_unequal = DistanceMatrix(
            [[0.0, 1.0, 0.1, 0.5678, 1.0, 1.0],
             [1.0, 0.0, 0.002, 0.42, 0.998, 0.0],
             [0.1, 0.002, 0.0, 1.0, 0.123, 1.0],
             [0.5678, 0.42, 1.0, 0.0, 0.123, 0.43],
             [1.0, 0.998, 0.123, 0.123, 0.0, 0.5],
             [1.0, 0.0, 1.0, 0.43, 0.5, 0.0]],
            ['s1', 's2', 's3', 's4', 's5', 's6'])

        # Expected series index is the same across all tests.
        self.exp_index = ['method name', 'test statistic name', 'sample size',
                          'number of groups', 'test statistic', 'p-value',
                          'number of permutations']

        # Stricter series equality testing than the default.
        self.assert_series_equal = partial(assert_series_equal,
                                           check_index_type=True,
                                           check_series_type=True)

    def test_call_ties(self):
        # Ensure we get the same results if we rerun the method using the same
        # inputs. Also ensure we get the same results if we run the method
        # using a grouping vector or a data frame with equivalent groupings.
        exp = pd.Series(index=self.exp_index,
                        data=['PERMANOVA', 'pseudo-F', 4, 2, 2.0, 0.68, 999],
                        name='PERMANOVA results')

        obs = permanova(self.dm_ties, self.grouping_equal, seed=42)
        # pstat can change slightly, depending on random source used
        self.assertAlmostEqual(obs.array[5],exp.array[5],delta=0.05)
        # update exp accordingly for further tests
        exp.array[5] = obs.array[5]

        for _ in range(2):
            obs = permanova(self.dm_ties, self.grouping_equal, seed=42)
            self.assert_series_equal(obs, exp)

        for _ in range(2):
            obs = permanova(self.dm_ties, self.df, column='Group', seed=42)
            self.assert_series_equal(obs, exp)

    def test_call_no_ties(self):
        exp = pd.Series(index=self.exp_index,
                        data=['PERMANOVA', 'pseudo-F', 4, 2, 4.4, 0.345, 999],
                        name='PERMANOVA results')
        obs = permanova(self.dm_no_ties, self.grouping_equal, seed=42)
        # pstat can change slightly, depending on random source used
        self.assertAlmostEqual(obs.array[5],exp.array[5],delta=0.05)
        # update exp accordingly for further tests
        exp.array[5] = obs.array[5]
        self.assert_series_equal(obs, exp)

    def test_call_no_permutations(self):
        exp = pd.Series(index=self.exp_index,
                        data=['PERMANOVA', 'pseudo-F', 4, 2, 4.4, np.nan, 0],
                        name='PERMANOVA results')
        obs = permanova(self.dm_no_ties, self.grouping_equal, permutations=0)
        self.assert_series_equal(obs, exp)

    def test_call_unequal_group_sizes(self):
        exp = pd.Series(
            index=self.exp_index,
            data=['PERMANOVA', 'pseudo-F', 6, 3, 0.578848, 0.655, 999],
            name='PERMANOVA results')

        obs = permanova(self.dm_unequal, self.grouping_unequal, seed=42)
        # pstat can change slightly, depending on random source used
        self.assertAlmostEqual(obs.array[5],exp.array[5],delta=0.05)
        # update exp accordingly for further tests
        exp.array[5] = obs.array[5]
        self.assert_series_equal(obs, exp)

        obs = permanova(self.dm_unequal, self.grouping_unequal_relabeled, seed=42)
        # pstat can change slightly, depending on random source used
        self.assertAlmostEqual(obs.array[5],exp.array[5],delta=0.05)
        # update exp accordingly for further tests
        exp.array[5] = obs.array[5]
        self.assert_series_equal(obs, exp)

    def test_call_via_series(self):
        # test https://github.com/scikit-bio/scikit-bio/issues/1877
        # permanova gives different results if grouping is either
        # a pd.DataFrame or a pd.Series
        dm = DistanceMatrix.read(get_data_path('frameSeries_dm.tsv'))
        grouping = pd.read_csv(get_data_path("frameSeries_grouping.tsv"),
                               sep="\t", index_col=0)

        obs_frame = permanova(dm, grouping, column='tumor', seed=42)

        obs_series = permanova(dm, grouping['tumor'], seed=42)

        # in principle, both tests - if seed is the same - should return the
        # exact same results. However, they don't for the current example ...
        self.assert_series_equal(obs_frame, obs_series)

        # ... which is due to different result in computing "unique" values for
        # the grouping, which is illustrated with the following test
        grp_frame = _preprocess_input_sng(
            dm.ids, dm.shape[0],
            grouping, 'tumor'  # grouping as a pd.DataFrame
            )[-1]
        grp_series = _preprocess_input_sng(
            dm.ids, dm.shape[0],
            grouping['tumor'], column=None  # grouping as a pd.Series, note
                                            # that user of permanova do not
                                            # have to explicitly set
                                            # column=None
            )[-1]
        # convert np.array to tuple to ease comparison for equality
        self.assertEqual(tuple(grp_frame), tuple(grp_series))

        # to better illustrate what is going wrong, we compare the computed
        # grouping labels (0 or 1) with the original user provided data, here
        # "no-tumor mice" and "tumor-bearing mice". We expect a one-to-one
        # correspondens, i.e. if we group on both columns at the same time, we
        # expect exactly two groups, like
        # tumor               series
        # no-tumor mice       0          5
        # tumor-bearing mice  1         37
        # dtype: int64
        # which is not the case of the pd.Series case
        g = pd.DataFrame(data={'series': list(grp_series),
                               'dataframe': list(grp_frame),
                               'tumor': grouping.loc[list(dm.ids), 'tumor']},
                         index=dm.ids)
        self.assertEqual(g.groupby(['tumor', 'dataframe']).size().shape[0], 2)
        self.assertEqual(g.groupby(['tumor', 'series']).size().shape[0], 2)

        # test that ValueError is raised, if use provided column does not match
        # the provided pd.Series name for grouping
        with self.assertRaises(ValueError):
            _preprocess_input_sng(dm.ids, dm.shape[0],
                                  grouping['tumor'], 'foo')

    def test_invalid_input(self):
        with self.assertRaises(TypeError):
            permanova(self.dm_ties.data, self.grouping_equal, seed=42)


class PERMANOVACondensedTests(TestCase):
    """Tests for PERMANOVA with condensed distance matrices.
    
    These tests verify that condensed and redundant forms of the same
    distance matrix produce identical results.
    """

    def setUp(self):
        # Distance matrices with and without ties in the ranks, with 2 groups
        # of equal size.
        dm_ids = ['s1', 's2', 's3', 's4']
        self.grouping_equal = ['Control', 'Control', 'Fast', 'Fast']
        self.df = pd.read_csv(
            io.StringIO('ID,Group\ns2,Control\ns3,Fast\ns4,Fast\ns5,Control\n'
                        's1,Control'), index_col=0)

        # Create both redundant and condensed versions
        self.dm_ties_redundant = DistanceMatrix([[0, 1, 1, 4],
                                                  [1, 0, 3, 2],
                                                  [1, 3, 0, 3],
                                                  [4, 2, 3, 0]], dm_ids)
        
        self.dm_ties_condensed = DistanceMatrix(
            squareform(self.dm_ties_redundant.data), #, checks=False),
            dm_ids,
            condensed=True
        )

        self.dm_no_ties_redundant = DistanceMatrix([[0, 1, 5, 4],
                                                     [1, 0, 3, 2],
                                                     [5, 3, 0, 3],
                                                     [4, 2, 3, 0]], dm_ids)
        
        self.dm_no_ties_condensed = DistanceMatrix(
            squareform(self.dm_no_ties_redundant.data), #, checks=False),
            dm_ids,
            condensed=True
        )

        # Test with 3 groups of unequal size.
        self.grouping_unequal = ['Control', 'Treatment1', 'Treatment2',
                                 'Treatment1', 'Control', 'Control']

        self.dm_unequal_redundant = DistanceMatrix(
            [[0.0, 1.0, 0.1, 0.5678, 1.0, 1.0],
             [1.0, 0.0, 0.002, 0.42, 0.998, 0.0],
             [0.1, 0.002, 0.0, 1.0, 0.123, 1.0],
             [0.5678, 0.42, 1.0, 0.0, 0.123, 0.43],
             [1.0, 0.998, 0.123, 0.123, 0.0, 0.5],
             [1.0, 0.0, 1.0, 0.43, 0.5, 0.0]],
            ['s1', 's2', 's3', 's4', 's5', 's6'])
        
        self.dm_unequal_condensed = DistanceMatrix(
            squareform(self.dm_unequal_redundant.data), #, checks=False),
            ['s1', 's2', 's3', 's4', 's5', 's6'],
            condensed=True
        )

        # Expected series index is the same across all tests.
        self.exp_index = ['method name', 'test statistic name', 'sample size',
                          'number of groups', 'test statistic', 'p-value',
                          'number of permutations']

        # Stricter series equality testing than the default.
        self.assert_series_equal = partial(assert_series_equal,
                                           check_index_type=True,
                                           check_series_type=True)

    def test_condensed_vs_redundant_ties(self):
        """Test that condensed and redundant forms give identical results with ties."""
        # Run with same seed for both
        obs_redundant = permanova(self.dm_ties_redundant, self.grouping_equal, seed=42)
        obs_condensed = permanova(self.dm_ties_condensed, self.grouping_equal, seed=42)
        
        self.assert_series_equal(obs_redundant, obs_condensed)

    def test_condensed_vs_redundant_no_ties(self):
        """Test that condensed and redundant forms give identical results without ties."""
        obs_redundant = permanova(self.dm_no_ties_redundant, self.grouping_equal, seed=42)
        obs_condensed = permanova(self.dm_no_ties_condensed, self.grouping_equal, seed=42)
        
        self.assert_series_equal(obs_redundant, obs_condensed)

    def test_condensed_vs_redundant_unequal_groups(self):
        """Test that condensed and redundant forms give identical results with unequal groups."""
        obs_redundant = permanova(self.dm_unequal_redundant, self.grouping_unequal, seed=42)
        obs_condensed = permanova(self.dm_unequal_condensed, self.grouping_unequal, seed=42)
        
        self.assert_series_equal(obs_redundant, obs_condensed)

    def test_condensed_with_dataframe_grouping(self):
        """Test that condensed form works with DataFrame grouping."""
        obs_redundant = permanova(self.dm_ties_redundant, self.df, column='Group', seed=42)
        obs_condensed = permanova(self.dm_ties_condensed, self.df, column='Group', seed=42)
        
        self.assert_series_equal(obs_redundant, obs_condensed)

    def test_condensed_no_permutations(self):
        """Test that condensed form works with no permutations."""
        obs_redundant = permanova(self.dm_no_ties_redundant, self.grouping_equal, 
                                  permutations=0)
        obs_condensed = permanova(self.dm_no_ties_condensed, self.grouping_equal, 
                                  permutations=0)
        
        self.assert_series_equal(obs_redundant, obs_condensed)

    def test_condensed_different_permutation_counts(self):
        """Test that condensed form gives consistent results across different permutation counts."""
        # Test with various permutation counts
        for n_perms in [99, 199, 999]:
            obs_redundant = permanova(self.dm_ties_redundant, self.grouping_equal, 
                                      permutations=n_perms, seed=42)
            obs_condensed = permanova(self.dm_ties_condensed, self.grouping_equal, 
                                      permutations=n_perms, seed=42)
            
            self.assert_series_equal(obs_redundant, obs_condensed)

    def test_condensed_multiple_runs_same_seed(self):
        """Test that condensed form gives identical results across multiple runs with same seed."""
        obs1 = permanova(self.dm_ties_condensed, self.grouping_equal, seed=42)
        obs2 = permanova(self.dm_ties_condensed, self.grouping_equal, seed=42)
        obs3 = permanova(self.dm_ties_condensed, self.grouping_equal, seed=42)
        
        self.assert_series_equal(obs1, obs2)
        self.assert_series_equal(obs2, obs3)

    def test_condensed_test_statistic_positive(self):
        """Test that condensed form produces positive F-statistics."""
        # This is a regression test for the bug where condensed matrices
        # produced negative F-statistics due to incorrect s_T calculation
        obs = permanova(self.dm_ties_condensed, self.grouping_equal, seed=42)
        
        # F-statistic should always be non-negative
        self.assertGreaterEqual(obs['test statistic'], 0.0)
        
    def test_condensed_test_statistic_matches_expected(self):
        """Test that condensed form produces expected F-statistic values."""
        # These expected values come from the redundant form tests
        obs_ties = permanova(self.dm_ties_condensed, self.grouping_equal, 
                            permutations=0)
        self.assertAlmostEqual(obs_ties['test statistic'], 2.0, places=5)
        
        obs_no_ties = permanova(self.dm_no_ties_condensed, self.grouping_equal, 
                               permutations=0)
        self.assertAlmostEqual(obs_no_ties['test statistic'], 4.4, places=5)

    def test_condensed_large_matrix(self):
        """Test condensed form with a larger matrix."""
        # Create a 10x10 distance matrix
        np.random.seed(42)
        n = 10
        redundant_data = np.random.rand(n, n)
        redundant_data = (redundant_data + redundant_data.T) / 2  # Make symmetric
        np.fill_diagonal(redundant_data, 0)  # Zero diagonal
        
        ids = [f's{i}' for i in range(n)]
        dm_redundant = DistanceMatrix(redundant_data, ids)
        dm_condensed = DistanceMatrix(squareform(redundant_data, checks=False), ids)
        
        grouping = ['A'] * 5 + ['B'] * 5
        
        obs_redundant = permanova(dm_redundant, grouping, seed=42, permutations=99)
        obs_condensed = permanova(dm_condensed, grouping, seed=42, permutations=99)
        
        self.assert_series_equal(obs_redundant, obs_condensed)

    def test_condensed_3x3_matrix(self):
        """Test condensed form with 3x3 matrix."""
        dm_redundant = DistanceMatrix([[0, 1, 2],
                                       [1, 0, 3],
                                       [2, 3, 0]], 
                                      ['s1', 's2', 's3'])
        dm_condensed = DistanceMatrix([1, 2, 3], ['s1', 's2', 's3'])
        
        grouping = ['A', 'A', 'B']
        
        obs_redundant = permanova(dm_redundant, grouping, seed=42, permutations=99)
        obs_condensed = permanova(dm_condensed, grouping, seed=42, permutations=99)
        
        self.assert_series_equal(obs_redundant, obs_condensed)

    def test_condensed_via_series(self):
        """Test condensed form with pd.Series grouping (regression test for #1877)."""
        dm_redundant = DistanceMatrix.read(get_data_path('frameSeries_dm.tsv'))
        dm_condensed = DistanceMatrix(
            squareform(dm_redundant.data, checks=False),
            dm_redundant.ids
        )
        grouping = pd.read_csv(get_data_path("frameSeries_grouping.tsv"),
                               sep="\t", index_col=0)

        # Test with DataFrame
        obs_redundant_frame = permanova(dm_redundant, grouping, column='tumor', seed=42)
        obs_condensed_frame = permanova(dm_condensed, grouping, column='tumor', seed=42)
        self.assert_series_equal(obs_redundant_frame, obs_condensed_frame)

        # Test with Series
        obs_redundant_series = permanova(dm_redundant, grouping['tumor'], seed=42)
        obs_condensed_series = permanova(dm_condensed, grouping['tumor'], seed=42)
        self.assert_series_equal(obs_redundant_series, obs_condensed_series)

        # Verify both grouping types give same results for condensed
        self.assert_series_equal(obs_condensed_frame, obs_condensed_series)

    def test_condensed_flags_set_correctly(self):
        """Test that condensed matrices have correct flags set."""
        self.assertTrue(self.dm_ties_condensed._flags["CONDENSED"])
        self.assertFalse(self.dm_ties_redundant._flags["CONDENSED"])

    def test_condensed_data_shape(self):
        """Test that condensed matrix data has correct shape."""
        n = 4  # 4x4 matrix
        expected_condensed_length = n * (n - 1) // 2

        self.assertEqual(len(self.dm_ties_condensed.data), expected_condensed_length)
        self.assertEqual(self.dm_ties_redundant.data.shape, (n, n))


if __name__ == '__main__':
    main()
