# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from functools import partial
from unittest import TestCase, main

import numpy as np
import pandas as pd
from pandas.util.testing import assert_series_equal

from skbio import DistanceMatrix
from skbio.stats.distance import permanova


class TestPERMANOVA(TestCase):
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
                        data=['PERMANOVA', 'pseudo-F', 4, 2, 2.0, 0.671, 999],
                        name='PERMANOVA results')

        for _ in range(2):
            np.random.seed(0)
            obs = permanova(self.dm_ties, self.grouping_equal)
            self.assert_series_equal(obs, exp)

        for _ in range(2):
            np.random.seed(0)
            obs = permanova(self.dm_ties, self.df, column='Group')
            self.assert_series_equal(obs, exp)

    def test_call_no_ties(self):
        exp = pd.Series(index=self.exp_index,
                        data=['PERMANOVA', 'pseudo-F', 4, 2, 4.4, 0.332, 999],
                        name='PERMANOVA results')
        np.random.seed(0)
        obs = permanova(self.dm_no_ties, self.grouping_equal)
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
            data=['PERMANOVA', 'pseudo-F', 6, 3, 0.578848, 0.645, 999],
            name='PERMANOVA results')

        np.random.seed(0)
        obs = permanova(self.dm_unequal, self.grouping_unequal)
        self.assert_series_equal(obs, exp)

        np.random.seed(0)
        obs = permanova(self.dm_unequal, self.grouping_unequal_relabeled)
        self.assert_series_equal(obs, exp)


if __name__ == '__main__':
    main()
