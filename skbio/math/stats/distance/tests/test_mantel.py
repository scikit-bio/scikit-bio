#! /usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd
from pandas.util.testing import assert_frame_equal

from skbio.core.distance import DistanceMatrix
from skbio.core.exception import DissimilarityMatrixError, DistanceMatrixError
from skbio.math.stats.distance import mantel
from skbio.util.testing import get_data_path


class MantelTests(TestCase):
    """Results were verified with R 3.1.0 and vegan 2.0-10 (vegan::mantel).
    
    vegan::mantel performs a one-sided (greater) test and does not have the
    option to specify different alternative hypotheses. In order to test the
    other alternative hypotheses, I modified vegan::mantel to perform the
    appropriate test, source()'d the file and verified the output.

    """

    # TODO: add test from Legendre (where r is undefined, no variation)

    # TODO: add test to ensure inputs aren't modified

    def setUp(self):
        self.methods = ('pearson', 'spearman')

        # Small dataset of minimal size (3x3). Mix of floats and ints in a
        # native Python nested list structure.
        self.minx = [[0, 1, 2], [1, 0, 3], [2, 3, 0]]
        self.miny = [[0, 2, 7], [2, 0, 6], [7, 6, 0]]
        self.minz = [[0, 0.5, 0.25], [0.5, 0, 0.1], [0.25, 0.1, 0]]

        # This second dataset is derived from vegan::mantel's example dataset.
        # The "veg" distance matrix contains Bray-Curtis distances derived from
        # the varespec data (named "veg.dist" in the example). The "env"
        # distance matrix contains Euclidean distances derived from scaled
        # varechem data (named "env.dist" in the example).
        self.veg_dm_vegan = np.loadtxt(
            get_data_path('mantel_veg_dm_vegan.txt'))
        self.env_dm_vegan = np.loadtxt(
            get_data_path('mantel_env_dm_vegan.txt'))

        # Expected test statistic when comparing x and y with method='pearson'.
        self.exp_x_vs_y = 0.7559289

        # Expected test statistic when comparing x and z with method='pearson'.
        self.exp_x_vs_z = -0.9897433

    def test_statistic_same_across_alternatives_and_permutations(self):
        # Varying permutations and alternative hypotheses shouldn't affect the
        # computed test statistics.
        for n in (0, 99, 999):
            for alt in ('twosided', 'greater', 'less'):
                for method, exp in (('pearson', self.exp_x_vs_y),
                                    ('spearman', 0.5)):
                    obs = mantel(self.minx, self.miny, method=method,
                                 permutations=n, alternative=alt)[0]
                    self.assertAlmostEqual(obs, exp)

    def test_comparing_same_matrices(self):
        for method in self.methods:
            obs = mantel(self.minx, self.minx, method=method)[0]
            self.assertAlmostEqual(obs, 1)

            obs = mantel(self.miny, self.miny, method=method)[0]
            self.assertAlmostEqual(obs, 1)

    def test_negative_correlation(self):
        for method, exp in (('pearson', self.exp_x_vs_z), ('spearman', -1)):
            obs = mantel(self.minx, self.minz, method=method)[0]
            self.assertAlmostEqual(obs, exp)

    def test_zero_permutations(self):
        for alt in ('twosided', 'greater', 'less'):
            for method, exp in (('pearson', self.exp_x_vs_y),
                                ('spearman', 0.5)):
                obs = mantel(self.minx, self.miny, permutations=0,
                             method=method, alternative=alt)
                self.assertAlmostEqual(obs[0], exp)
                npt.assert_equal(obs[1], np.nan)

                # TODO test swapping order of matrices -- should get same
                # result

    def test_one_sided_greater(self):
        np.random.seed(0)

        obs = mantel(self.minx, self.miny, alternative='greater')
        self.assertAlmostEqual(obs[0], self.exp_x_vs_y)
        self.assertAlmostEqual(obs[1], 0.324)

        obs = mantel(self.minx, self.minx, alternative='greater')
        self.assertAlmostEqual(obs[0], 1)
        self.assertAlmostEqual(obs[1], 0.172)

    def test_one_sided_less(self):
        # no need to seed here as permuted test statistics will all be less
        # than or equal to the observed test statistic (1.0)
        for method in self.methods:
            obs = mantel(self.minx, self.minx, method=method,
                         alternative='less')
            self.assertEqual(obs, (1, 1))

        np.random.seed(0)

        obs = mantel(self.minx, self.miny, alternative='less')
        self.assertAlmostEqual(obs[0], self.exp_x_vs_y)
        self.assertAlmostEqual(obs[1], 0.843)

        obs = mantel(self.minx, self.minz, alternative='less')
        self.assertAlmostEqual(obs[0], self.exp_x_vs_z)
        self.assertAlmostEqual(obs[1], 0.172)

    def test_two_sided(self):
        np.random.seed(0)

        obs = mantel(self.minx, self.minx, method='spearman',
                     alternative='twosided')
        self.assertEqual(obs[0], 1)
        self.assertAlmostEqual(obs[1], 0.328)

        obs = mantel(self.minx, self.miny, method='spearman',
                     alternative='twosided')
        self.assertAlmostEqual(obs[0], 0.5)
        self.assertAlmostEqual(obs[1], 1.0)

        obs = mantel(self.minx, self.minz, method='spearman',
                     alternative='twosided')
        self.assertAlmostEqual(obs[0], -1)
        self.assertAlmostEqual(obs[1], 0.322)

    def test_vegan_example(self):
        np.random.seed(0)

        # pearson
        obs = mantel(self.veg_dm_vegan, self.env_dm_vegan,
                     alternative='greater')
        self.assertAlmostEqual(obs[0], 0.3047454)
        self.assertAlmostEqual(obs[1], 0.002)

        # spearman
        obs = mantel(self.veg_dm_vegan, self.env_dm_vegan,
                     alternative='greater', method='spearman')
        self.assertAlmostEqual(obs[0], 0.283791)
        self.assertAlmostEqual(obs[1], 0.003)

    def test_mantel_invalid_distance_matrix(self):
        # Single asymmetric, non-hollow distance matrix.
        with self.assertRaises(DissimilarityMatrixError):
            mantel([[1, 2], [3, 4]], [[0, 0], [0, 0]])

        # Two asymmetric distance matrices.
        with self.assertRaises(DistanceMatrixError):
            mantel([[0, 2], [3, 0]], [[0, 1], [0, 0]])

    def test_mantel_invalid_input(self):
        # invalid correlation method
        with self.assertRaises(ValueError):
            mantel([[1]], [[1]], method='brofist')

        # invalid permutations
        with self.assertRaises(ValueError):
            mantel([[1]], [[1]], permutations=-1)

        # invalid alternative
        with self.assertRaises(ValueError):
            mantel([[1]], [[1]], alternative='no cog yay')

        # mismatched shape
        with self.assertRaises(ValueError):
            mantel(self.minx, [[0, 2], [2, 0]])


if __name__ == '__main__':
    main()
