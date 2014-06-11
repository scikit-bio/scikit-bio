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
from skbio.math.stats.distance import mantel, pwmantel
from skbio.util.testing import get_data_path


class MantelTests(TestCase):
    """Results were verified with R 3.1.0 and vegan 2.0-10 (vegan::mantel).

    vegan::mantel performs a one-sided (greater) test and does not have the
    option to specify different alternative hypotheses. In order to test the
    other alternative hypotheses, I modified vegan::mantel to perform the
    appropriate test, source()'d the file and verified the output.

    """

    def setUp(self):
        self.methods = ('pearson', 'spearman')
        self.alternatives = ('two-sided', 'greater', 'less')

        # Small dataset of minimal size (3x3). Mix of floats and ints in a
        # native Python nested list structure.
        self.minx = [[0, 1, 2], [1, 0, 3], [2, 3, 0]]
        self.miny = [[0, 2, 7], [2, 0, 6], [7, 6, 0]]
        self.minz = [[0, 0.5, 0.25], [0.5, 0, 0.1], [0.25, 0.1, 0]]

        # No variation in distances. Taken from Figure 10.20(b), pg. 603 in L&L
        # 3rd edition. Their example is 4x4 but using 3x3 here for easy
        # comparison to the minimal dataset above.
        self.no_variation = [[0, 0.667, 0.667],
                             [0.667, 0, 0.667],
                             [0.667, 0.667, 0]]

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
            for alt in self.alternatives:
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
        for alt in self.alternatives:
            for method, exp in (('pearson', self.exp_x_vs_y),
                                ('spearman', 0.5)):
                obs = mantel(self.minx, self.miny, permutations=0,
                             method=method, alternative=alt)
                self.assertAlmostEqual(obs[0], exp)
                npt.assert_equal(obs[1], np.nan)

                # swapping order of matrices should give same result
                obs = mantel(self.miny, self.minx, permutations=0,
                             method=method, alternative=alt)
                self.assertAlmostEqual(obs[0], exp)
                npt.assert_equal(obs[1], np.nan)

    def test_distance_matrix_instances_as_input(self):
        # IDs shouldn't matter -- the function should only care about the
        # matrix data
        dmx = DistanceMatrix(self.minx)
        dmy = DistanceMatrix(self.miny, ['no', 'cog', 'yay'])

        np.random.seed(0)

        obs = mantel(dmx, dmy, alternative='less')

        self.assertAlmostEqual(obs[0], self.exp_x_vs_y)
        self.assertAlmostEqual(obs[1], 0.843)

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
                     alternative='two-sided')
        self.assertEqual(obs[0], 1)
        self.assertAlmostEqual(obs[1], 0.328)

        obs = mantel(self.minx, self.miny, method='spearman',
                     alternative='two-sided')
        self.assertAlmostEqual(obs[0], 0.5)
        self.assertAlmostEqual(obs[1], 1.0)

        obs = mantel(self.minx, self.minz, method='spearman',
                     alternative='two-sided')
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

    def test_no_variation_pearson(self):
        # Output doesn't match vegan::mantel with method='pearson'. Consider
        # revising output and this test depending on outcome of
        # https://github.com/scipy/scipy/issues/3728
        for alt in self.alternatives:
            # test one or both inputs having no variation in their
            # distances
            obs = mantel(self.miny, self.no_variation, method='pearson',
                         alternative=alt)
            npt.assert_equal(obs, (0.0, 1.0))

            obs = mantel(self.no_variation, self.miny, method='pearson',
                         alternative=alt)
            npt.assert_equal(obs, (0.0, 1.0))

            obs = mantel(self.no_variation, self.no_variation,
                         method='pearson', alternative=alt)
            npt.assert_equal(obs, (1.0, 1.0))

    def test_no_variation_spearman(self):
        exp = (np.nan, np.nan)
        for alt in self.alternatives:
            obs = mantel(self.miny, self.no_variation, method='spearman',
                         alternative=alt)
            npt.assert_equal(obs, exp)

            obs = mantel(self.no_variation, self.miny, method='spearman',
                         alternative=alt)
            npt.assert_equal(obs, exp)

            obs = mantel(self.no_variation, self.no_variation,
                         method='spearman', alternative=alt)
            npt.assert_equal(obs, exp)

    def test_no_side_effects(self):
        minx = np.asarray(self.minx, dtype='float')
        miny = np.asarray(self.miny, dtype='float')

        minx_copy = np.copy(minx)
        miny_copy = np.copy(miny)

        mantel(minx, miny)

        # Make sure we haven't modified the input.
        npt.assert_equal(minx, minx_copy)
        npt.assert_equal(miny, miny_copy)

    def test_invalid_distance_matrix(self):
        # Single asymmetric, non-hollow distance matrix.
        with self.assertRaises(DissimilarityMatrixError):
            mantel([[1, 2], [3, 4]], [[0, 0], [0, 0]])

        # Two asymmetric distance matrices.
        with self.assertRaises(DistanceMatrixError):
            mantel([[0, 2], [3, 0]], [[0, 1], [0, 0]])

    def test_invalid_input(self):
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

        # too small dms
        with self.assertRaises(ValueError):
            mantel([[0, 3], [3, 0]], [[0, 2], [2, 0]])


class PairwiseMantelTests(TestCase):

    # TODO add test with duplicate dms

    def setUp(self):
        self.minx = DistanceMatrix([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        self.miny = DistanceMatrix([[0, 2, 7], [2, 0, 6], [7, 6, 0]])
        self.minz = DistanceMatrix([[0, 0.5, 0.25],
                                    [0.5, 0, 0.1],
                                    [0.25, 0.1, 0]])
        self.min_dms = (self.minx, self.miny, self.minz)

        # Load expected results.
        self.exp_results_minimal = pd.read_csv(
            get_data_path('pwmantel_exp_results_minimal.txt'), sep='\t',
            index_col=(0, 1))
        self.exp_results_minimal_with_labels = pd.read_csv(
            get_data_path('pwmantel_exp_results_minimal_with_labels.txt'),
            sep='\t', index_col=(0, 1))

    def test_minimal_compatible_input(self):
        #obs.to_csv('tests/data/pwmantel_exp_results_minimal.txt', sep='\t')
        np.random.seed(0)

        obs = pwmantel(self.min_dms, alternative='greater')
        assert_frame_equal(obs, self.exp_results_minimal)

    def test_minimal_compatible_input_with_labels(self):
        np.random.seed(0)

        obs = pwmantel(self.min_dms, alternative='greater',
                       labels=('minx', 'miny', 'minz'))
        assert_frame_equal(obs, self.exp_results_minimal_with_labels)


if __name__ == '__main__':
    main()
