#! /usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division, print_function

import numpy as np

from bipy.core.distance import SymmetricDistanceMatrix
from bipy.maths.stats.distance.anosim import ANOSIM
from bipy.util.unit_test import TestCase, main


class ANOSIMTests(TestCase):

    def setUp(self):
        # Define two small dms for easy testing. One has ties in the ranks.
        dm_sids = ['s1', 's2', 's3', 's4']
        self.grouping = ['Control', 'Control', 'Fast', 'Fast']

        self.dm_ties = SymmetricDistanceMatrix([[0, 1, 1, 4],
                                                [1, 0, 3, 2],
                                                [1, 3, 0, 3],
                                                [4, 2, 3, 0]], dm_sids)

        self.dm_no_ties = SymmetricDistanceMatrix([[0, 1, 5, 4],
                                                   [1, 0, 3, 2],
                                                   [5, 3, 0, 3],
                                                   [4, 2, 3, 0]], dm_sids)

        self.anosim_ties = ANOSIM(self.dm_ties, self.grouping)
        self.anosim_no_ties = ANOSIM(self.dm_no_ties, self.grouping)

    def test_call_ties(self):
        # These results were verified with R.
        exp_r_stat = 0.25
        exp_p_val = 0.671

        np.random.seed(0)
        obs = self.anosim_ties()
        self.assertFloatEqual(obs.r_statistic, exp_r_stat)
        self.assertFloatEqual(obs.p_value, exp_p_val)

        # Ensure we get the same results if we rerun the method on the same
        # object.
        np.random.seed(0)
        obs = self.anosim_ties()
        self.assertFloatEqual(obs.r_statistic, exp_r_stat)
        self.assertFloatEqual(obs.p_value, exp_p_val)

    def test_call_no_ties(self):
        # These results were verified with R.
        np.random.seed(0)
        obs = self.anosim_no_ties()
        self.assertFloatEqual(obs.r_statistic, 0.625)
        self.assertFloatEqual(obs.p_value, 0.332)

    def test_call_no_ties(self):
        # These results were verified with R.
        np.random.seed(0)
        obs = self.anosim_no_ties()
        self.assertFloatEqual(obs.r_statistic, 0.625)
        self.assertFloatEqual(obs.p_value, 0.332)

    def test_call_no_permutations(self):
        # These results were verified with R.
        obs = self.anosim_no_ties(0)
        self.assertFloatEqual(obs.r_statistic, 0.625)
        self.assertEqual(obs.p_value, None)


if __name__ == '__main__':
    main()
