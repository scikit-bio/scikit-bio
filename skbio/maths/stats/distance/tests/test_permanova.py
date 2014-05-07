#! /usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division
from future.utils.six import StringIO
from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd

from skbio.core.distance import DistanceMatrix
from skbio.maths.stats.distance.permanova import PERMANOVA


class PERMANOVATests(TestCase):
    """All results were verified with R (vegan::adonis)."""

    def setUp(self):
        # Distance matrices with and without ties in the ranks, with 2 groups
        # of equal size.
        dm_ids = ['s1', 's2', 's3', 's4']
        grouping_equal = ['Control', 'Control', 'Fast', 'Fast']
        df = pd.read_csv(
            StringIO('ID,Group\ns2,Control\ns3,Fast\ns4,Fast\ns5,Control\n'
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
        grouping_unequal = ['Control', 'Treatment1', 'Treatment2',
                            'Treatment1', 'Control', 'Control']

        self.dm_unequal = DistanceMatrix(
            [[0.0, 1.0, 0.1, 0.5678, 1.0, 1.0],
             [1.0, 0.0, 0.002, 0.42, 0.998, 0.0],
             [0.1, 0.002, 0.0, 1.0, 0.123, 1.0],
             [0.5678, 0.42, 1.0, 0.0, 0.123, 0.43],
             [1.0, 0.998, 0.123, 0.123, 0.0, 0.5],
             [1.0, 0.0, 1.0, 0.43, 0.5, 0.0]],
            ['s1', 's2', 's3', 's4', 's5', 's6'])

        self.permanova_ties = PERMANOVA(self.dm_ties, grouping_equal)
        self.permanova_no_ties = PERMANOVA(self.dm_no_ties, grouping_equal)
        self.permanova_ties_df = PERMANOVA(self.dm_ties, df, column='Group')
        self.permanova_unequal = PERMANOVA(self.dm_unequal, grouping_unequal)

    def test_call_ties(self):
        # Ensure we get the same results if we rerun the method on the same
        # object. Also ensure we get the same results if we run the method
        # using a grouping vector or a data frame with equivalent groupings.
        for inst in self.permanova_ties, self.permanova_ties_df:
            for trial in range(2):
                np.random.seed(0)
                obs = inst()
                self.assertEqual(obs.sample_size, 4)
                npt.assert_array_equal(obs.groups, ['Control', 'Fast'])
                self.assertAlmostEqual(obs.statistic, 2.0)
                self.assertAlmostEqual(obs.p_value, 0.671)
                self.assertEqual(obs.permutations, 999)

    def test_call_no_ties(self):
        np.random.seed(0)
        obs = self.permanova_no_ties()
        self.assertEqual(obs.sample_size, 4)
        npt.assert_array_equal(obs.groups, ['Control', 'Fast'])
        self.assertAlmostEqual(obs.statistic, 4.4)
        self.assertAlmostEqual(obs.p_value, 0.332)
        self.assertEqual(obs.permutations, 999)

    def test_call_no_permutations(self):
        obs = self.permanova_no_ties(0)
        self.assertEqual(obs.sample_size, 4)
        npt.assert_array_equal(obs.groups, ['Control', 'Fast'])
        self.assertAlmostEqual(obs.statistic, 4.4)
        self.assertEqual(obs.p_value, None)
        self.assertEqual(obs.permutations, 0)

    def test_call_unequal_group_sizes(self):
        np.random.seed(0)
        obs = self.permanova_unequal()
        self.assertEqual(obs.sample_size, 6)
        npt.assert_array_equal(obs.groups,
                               ['Control', 'Treatment1', 'Treatment2'])
        self.assertAlmostEqual(obs.statistic, 0.578848, 6)
        self.assertAlmostEqual(obs.p_value, 0.645)
        self.assertEqual(obs.permutations, 999)


if __name__ == '__main__':
    main()
