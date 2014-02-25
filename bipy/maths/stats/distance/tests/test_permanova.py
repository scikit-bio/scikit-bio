#! /usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The bipy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division

import numpy as np

from bipy.core.distance import DistanceMatrix, SymmetricDistanceMatrix
from bipy.maths.stats.distance.permanova import PERMANOVA
from bipy.util.unit_test import TestCase, main


class PERMANOVATests(TestCase):
    """All results were verified with R (vegan::adonis)."""

    def setUp(self):
        # Distance matrices with and without ties in the ranks, with 2 groups
        # of equal size.
        dm_sids = ['s1', 's2', 's3', 's4']
        grouping_equal = ['Control', 'Control', 'Fast', 'Fast']

        self.dm_ties = SymmetricDistanceMatrix([[0, 1, 1, 4],
                                                [1, 0, 3, 2],
                                                [1, 3, 0, 3],
                                                [4, 2, 3, 0]], dm_sids)

        self.dm_no_ties = SymmetricDistanceMatrix([[0, 1, 5, 4],
                                                   [1, 0, 3, 2],
                                                   [5, 3, 0, 3],
                                                   [4, 2, 3, 0]], dm_sids)

        # Test with 3 groups of unequal size.
        grouping_unequal = ['Control', 'Treatment1', 'Treatment2',
                            'Treatment1', 'Control', 'Control']

        self.dm_unequal = SymmetricDistanceMatrix(
            [[0.0, 1.0, 0.1, 0.5678, 1.0, 1.0],
             [1.0, 0.0, 0.002, 0.42, 0.998, 0.0],
             [0.1, 0.002, 0.0, 1.0, 0.123, 1.0],
             [0.5678, 0.42, 1.0, 0.0, 0.123, 0.43],
             [1.0, 0.998, 0.123, 0.123, 0.0, 0.5],
             [1.0, 0.0, 1.0, 0.43, 0.5, 0.0]],
            ['s1', 's2', 's3', 's4', 's5', 's6'])

        self.permanova_ties = PERMANOVA(self.dm_ties, grouping_equal)
        self.permanova_no_ties = PERMANOVA(self.dm_no_ties, grouping_equal)
        self.permanova_unequal = PERMANOVA(self.dm_unequal, grouping_unequal)

    def test_init_invalid_input(self):
        # Requires a SymmetricDistanceMatrix.
        with self.assertRaises(TypeError):
            _ = PERMANOVA(DistanceMatrix([[0, 2], [3, 0]], ['a', 'b']), [1, 2])

        # Grouping vector length must match number of samples in dm.
        with self.assertRaises(ValueError):
            _ = PERMANOVA(self.dm_ties, [1, 2])

        # Grouping vector cannot have only unique values.
        with self.assertRaises(ValueError):
            _ = PERMANOVA(self.dm_ties, [1, 2, 3, 4])

        # Grouping vector cannot have only a single group.
        with self.assertRaises(ValueError):
            _ = PERMANOVA(self.dm_ties, [1, 1, 1, 1])

    def test_call_ties(self):
        # TODO: update tests to use CategoricalStatsResults objects to store
        # expected results
        exp_f_stat = 2.0
        exp_p_val = 0.671

        np.random.seed(0)
        obs = self.permanova_ties()
        self.assertFloatEqual(obs.statistic, exp_f_stat)
        self.assertFloatEqual(obs.p_value, exp_p_val)

        # Ensure we get the same results if we rerun the method on the same
        # object.
        np.random.seed(0)
        obs = self.permanova_ties()
        self.assertFloatEqual(obs.statistic, exp_f_stat)
        self.assertFloatEqual(obs.p_value, exp_p_val)

    def test_call_no_ties(self):
        np.random.seed(0)
        obs = self.permanova_no_ties()
        self.assertFloatEqual(obs.statistic, 4.4)
        self.assertFloatEqual(obs.p_value, 0.332)

    def test_call_no_permutations(self):
        obs = self.permanova_no_ties(0)
        self.assertFloatEqual(obs.statistic, 4.4)
        self.assertEqual(obs.p_value, None)

    def test_call_invalid_permutations(self):
        with self.assertRaises(ValueError):
            _ = self.permanova_ties(-1)

    def test_call_unequal_group_sizes(self):
        np.random.seed(0)
        obs = self.permanova_unequal()
        self.assertFloatEqual(obs.statistic, 0.578848)
        self.assertFloatEqual(obs.p_value, 0.645)


if __name__ == '__main__':
    main()
