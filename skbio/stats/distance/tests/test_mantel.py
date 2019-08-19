# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd

from skbio import DistanceMatrix
from skbio.stats.distance import (DissimilarityMatrixError,
                                  DistanceMatrixError, mantel, pwmantel)
from skbio.stats.distance._mantel import _order_dms
from skbio.util import get_data_path, assert_data_frame_almost_equal


class MantelTestData(TestCase):
    def setUp(self):
        # Small dataset of minimal size (3x3). Mix of floats and ints in a
        # native Python nested list structure.
        self.minx = [[0, 1, 2], [1, 0, 3], [2, 3, 0]]
        self.miny = [[0, 2, 7], [2, 0, 6], [7, 6, 0]]
        self.minz = [[0, 0.5, 0.25], [0.5, 0, 0.1], [0.25, 0.1, 0]]

        # Version of the above dataset stored as DistanceMatrix instances.
        self.minx_dm = DistanceMatrix(self.minx)
        self.miny_dm = DistanceMatrix(self.miny)
        self.minz_dm = DistanceMatrix(self.minz)

        # Versions of self.minx_dm and self.minz_dm that each have an extra ID
        # on the end.
        self.minx_dm_extra = DistanceMatrix([[0, 1, 2, 7],
                                             [1, 0, 3, 2],
                                             [2, 3, 0, 4],
                                             [7, 2, 4, 0]],
                                            ['0', '1', '2', 'foo'])
        self.minz_dm_extra = DistanceMatrix([[0, 0.5, 0.25, 3],
                                             [0.5, 0, 0.1, 24],
                                             [0.25, 0.1, 0, 5],
                                             [3, 24, 5, 0]],
                                            ['0', '1', '2', 'bar'])


class MantelTests(MantelTestData):
    """Results were verified with R 3.1.0 and vegan 2.0-10 (vegan::mantel).

    vegan::mantel performs a one-sided (greater) test and does not have the
    option to specify different alternative hypotheses. In order to test the
    other alternative hypotheses, I modified vegan::mantel to perform the
    appropriate test, source()'d the file and verified the output.

    """

    def setUp(self):
        super(MantelTests, self).setUp()

        self.methods = ('pearson', 'spearman')
        self.alternatives = ('two-sided', 'greater', 'less')

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
                self.assertEqual(obs[2], 3)

                # swapping order of matrices should give same result
                obs = mantel(self.miny, self.minx, permutations=0,
                             method=method, alternative=alt)
                self.assertAlmostEqual(obs[0], exp)
                npt.assert_equal(obs[1], np.nan)
                self.assertEqual(obs[2], 3)

    def test_distance_matrix_instances_as_input(self):
        # Matrices with all matching IDs in the same order.
        np.random.seed(0)

        obs = mantel(self.minx_dm, self.miny_dm, alternative='less')

        self.assertAlmostEqual(obs[0], self.exp_x_vs_y)
        self.assertAlmostEqual(obs[1], 0.843)
        self.assertEqual(obs[2], 3)

    def test_distance_matrix_instances_with_reordering_and_nonmatching(self):
        x = self.minx_dm_extra.filter(['1', '0', 'foo', '2'])
        y = self.miny_dm.filter(['0', '2', '1'])

        # strict=True should disallow IDs that aren't found in both matrices
        with self.assertRaises(ValueError):
            mantel(x, y, alternative='less', strict=True)

        np.random.seed(0)

        # strict=False should ignore IDs that aren't found in both matrices
        obs = mantel(x, y, alternative='less', strict=False)

        self.assertAlmostEqual(obs[0], self.exp_x_vs_y)
        self.assertAlmostEqual(obs[1], 0.843)
        self.assertEqual(obs[2], 3)

    def test_distance_matrix_instances_with_lookup(self):
        self.minx_dm.ids = ('a', 'b', 'c')
        self.miny_dm.ids = ('d', 'e', 'f')
        lookup = {'a': 'A', 'b': 'B', 'c': 'C',
                  'd': 'A', 'e': 'B', 'f': 'C'}

        np.random.seed(0)

        obs = mantel(self.minx_dm, self.miny_dm, alternative='less',
                     lookup=lookup)
        self.assertAlmostEqual(obs[0], self.exp_x_vs_y)
        self.assertAlmostEqual(obs[1], 0.843)
        self.assertEqual(obs[2], 3)

    def test_one_sided_greater(self):
        np.random.seed(0)

        obs = mantel(self.minx, self.miny, alternative='greater')
        self.assertAlmostEqual(obs[0], self.exp_x_vs_y)
        self.assertAlmostEqual(obs[1], 0.324)
        self.assertEqual(obs[2], 3)

        obs = mantel(self.minx, self.minx, alternative='greater')
        self.assertAlmostEqual(obs[0], 1)
        self.assertAlmostEqual(obs[1], 0.172)
        self.assertEqual(obs[2], 3)

    def test_one_sided_less(self):
        # no need to seed here as permuted test statistics will all be less
        # than or equal to the observed test statistic (1.0)
        for method in self.methods:
            obs = mantel(self.minx, self.minx, method=method,
                         alternative='less')
            npt.assert_almost_equal(obs, (1, 1, 3))

        np.random.seed(0)

        obs = mantel(self.minx, self.miny, alternative='less')
        self.assertAlmostEqual(obs[0], self.exp_x_vs_y)
        self.assertAlmostEqual(obs[1], 0.843)
        self.assertEqual(obs[2], 3)

        obs = mantel(self.minx, self.minz, alternative='less')
        self.assertAlmostEqual(obs[0], self.exp_x_vs_z)
        self.assertAlmostEqual(obs[1], 0.172)
        self.assertEqual(obs[2], 3)

    def test_two_sided(self):
        np.random.seed(0)

        obs = mantel(self.minx, self.minx, method='spearman',
                     alternative='two-sided')
        self.assertEqual(obs[0], 1)
        self.assertAlmostEqual(obs[1], 0.328)
        self.assertEqual(obs[2], 3)

        obs = mantel(self.minx, self.miny, method='spearman',
                     alternative='two-sided')
        self.assertAlmostEqual(obs[0], 0.5)
        self.assertAlmostEqual(obs[1], 1.0)
        self.assertEqual(obs[2], 3)

        obs = mantel(self.minx, self.minz, method='spearman',
                     alternative='two-sided')
        self.assertAlmostEqual(obs[0], -1)
        self.assertAlmostEqual(obs[1], 0.322)
        self.assertEqual(obs[2], 3)

    def test_vegan_example(self):
        np.random.seed(0)

        # pearson
        obs = mantel(self.veg_dm_vegan, self.env_dm_vegan,
                     alternative='greater')
        self.assertAlmostEqual(obs[0], 0.3047454)
        self.assertAlmostEqual(obs[1], 0.002)
        self.assertEqual(obs[2], 24)

        # spearman
        obs = mantel(self.veg_dm_vegan, self.env_dm_vegan,
                     alternative='greater', method='spearman')
        self.assertAlmostEqual(obs[0], 0.283791)
        self.assertAlmostEqual(obs[1], 0.003)
        self.assertEqual(obs[2], 24)

    def test_no_variation_pearson(self):
        # Output doesn't match vegan::mantel with method='pearson'. Consider
        # revising output and this test depending on outcome of
        # https://github.com/scipy/scipy/issues/3728
        for alt in self.alternatives:
            # test one or both inputs having no variation in their
            # distances
            obs = mantel(self.miny, self.no_variation, method='pearson',
                         alternative=alt)
            npt.assert_equal(obs, (np.nan, np.nan, 3))

            obs = mantel(self.no_variation, self.miny, method='pearson',
                         alternative=alt)
            npt.assert_equal(obs, (np.nan, np.nan, 3))

            obs = mantel(self.no_variation, self.no_variation,
                         method='pearson', alternative=alt)
            npt.assert_equal(obs, (np.nan, np.nan, 3))

    def test_no_variation_spearman(self):
        exp = (np.nan, np.nan, 3)
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

        # too small dms
        with self.assertRaises(ValueError):
            mantel([[0, 3], [3, 0]], [[0, 2], [2, 0]])


class PairwiseMantelTests(MantelTestData):
    def setUp(self):
        super(PairwiseMantelTests, self).setUp()

        self.min_dms = (self.minx_dm, self.miny_dm, self.minz_dm)

        self.exp_results_minimal = pd.read_csv(
            get_data_path('pwmantel_exp_results_minimal.txt'), sep='\t',
            index_col=(0, 1))

        self.exp_results_minimal_with_labels = pd.read_csv(
            get_data_path('pwmantel_exp_results_minimal_with_labels.txt'),
            sep='\t', index_col=(0, 1))

        self.exp_results_duplicate_dms = pd.read_csv(
            get_data_path('pwmantel_exp_results_duplicate_dms.txt'),
            sep='\t', index_col=(0, 1))

        self.exp_results_na_p_value = pd.read_csv(
            get_data_path('pwmantel_exp_results_na_p_value.txt'),
            sep='\t', index_col=(0, 1))

        self.exp_results_reordered_distance_matrices = pd.read_csv(
            get_data_path('pwmantel_exp_results_reordered_distance_matrices'
                          '.txt'),
            sep='\t', index_col=(0, 1))

        self.exp_results_dm_dm2 = pd.read_csv(
            get_data_path('pwmantel_exp_results_dm_dm2.txt'),
            sep='\t', index_col=(0, 1))

        self.exp_results_all_dms = pd.read_csv(
            get_data_path('pwmantel_exp_results_all_dms.txt'),
            sep='\t', index_col=(0, 1))

    def test_minimal_compatible_input(self):
        # Matrices are already in the correct order and have matching IDs.
        np.random.seed(0)

        # input as DistanceMatrix instances
        obs = pwmantel(self.min_dms, alternative='greater')
        assert_data_frame_almost_equal(obs, self.exp_results_minimal)

        np.random.seed(0)

        # input as array_like
        obs = pwmantel((self.minx, self.miny, self.minz),
                       alternative='greater')
        assert_data_frame_almost_equal(obs, self.exp_results_minimal)

    def test_minimal_compatible_input_with_labels(self):
        np.random.seed(0)

        obs = pwmantel(self.min_dms, alternative='greater',
                       labels=('minx', 'miny', 'minz'))
        assert_data_frame_almost_equal(
            obs,
            self.exp_results_minimal_with_labels)

    def test_duplicate_dms(self):
        obs = pwmantel((self.minx_dm, self.minx_dm, self.minx_dm),
                       alternative='less')
        assert_data_frame_almost_equal(obs, self.exp_results_duplicate_dms)

    def test_na_p_value(self):
        obs = pwmantel((self.miny_dm, self.minx_dm), method='spearman',
                       permutations=0)
        assert_data_frame_almost_equal(obs, self.exp_results_na_p_value)

    def test_reordered_distance_matrices(self):
        # Matrices have matching IDs but they all have different ordering.
        x = self.minx_dm.filter(['1', '0', '2'])
        y = self.miny_dm.filter(['0', '2', '1'])
        z = self.minz_dm.filter(['1', '2', '0'])

        np.random.seed(0)

        obs = pwmantel((x, y, z), alternative='greater')
        assert_data_frame_almost_equal(
            obs,
            self.exp_results_reordered_distance_matrices)

    def test_strict(self):
        # Matrices have some matching and nonmatching IDs, with different
        # ordering.
        x = self.minx_dm_extra.filter(['1', '0', 'foo', '2'])
        y = self.miny_dm.filter(['0', '2', '1'])
        z = self.minz_dm_extra.filter(['bar', '1', '2', '0'])

        np.random.seed(0)

        # strict=False should discard IDs that aren't found in both matrices
        obs = pwmantel((x, y, z), alternative='greater', strict=False)
        assert_data_frame_almost_equal(
            obs,
            self.exp_results_reordered_distance_matrices)

    def test_id_lookup(self):
        # Matrices have mismatched IDs but a lookup is provided.
        self.minx_dm_extra.ids = ['a', 'b', 'c', 'foo']
        self.minz_dm_extra.ids = ['d', 'e', 'f', 'bar']
        lookup = {'a': '0', 'b': '1', 'c': '2', 'foo': 'foo',
                  'd': '0', 'e': '1', 'f': '2', 'bar': 'bar',
                  '0': '0', '1': '1', '2': '2'}

        x = self.minx_dm_extra.filter(['b', 'a', 'foo', 'c'])
        y = self.miny_dm.filter(['0', '2', '1'])
        z = self.minz_dm_extra.filter(['bar', 'e', 'f', 'd'])

        x_copy = x.copy()
        y_copy = y.copy()
        z_copy = z.copy()

        np.random.seed(0)

        obs = pwmantel((x, y, z), alternative='greater', strict=False,
                       lookup=lookup)
        assert_data_frame_almost_equal(
            obs,
            self.exp_results_reordered_distance_matrices)

        # Make sure the inputs aren't modified.
        self.assertEqual(x, x_copy)
        self.assertEqual(y, y_copy)
        self.assertEqual(z, z_copy)

    def test_too_few_dms(self):
        with self.assertRaises(ValueError):
            pwmantel([self.miny_dm])

    def test_wrong_number_of_labels(self):
        with self.assertRaises(ValueError):
            pwmantel(self.min_dms, labels=['foo', 'bar'])

    def test_duplicate_labels(self):
        with self.assertRaises(ValueError):
            pwmantel(self.min_dms, labels=['foo', 'bar', 'foo'])

    def test_mixed_input_types(self):
        # DistanceMatrix, DistanceMatrix, array_like
        with self.assertRaises(TypeError):
            pwmantel((self.miny_dm, self.minx_dm, self.minz))

    def test_filepaths_as_input(self):
        dms = [
            get_data_path('dm.txt'),
            get_data_path('dm2.txt'),
        ]
        np.random.seed(0)

        obs = pwmantel(dms)
        assert_data_frame_almost_equal(obs, self.exp_results_dm_dm2)

    def test_many_filepaths_as_input(self):
        dms = [
            get_data_path('dm2.txt'),
            get_data_path('dm.txt'),
            get_data_path('dm4.txt'),
            get_data_path('dm3.txt')
        ]
        np.random.seed(0)

        obs = pwmantel(dms)
        assert_data_frame_almost_equal(obs, self.exp_results_all_dms)


class OrderDistanceMatricesTests(MantelTestData):
    def setUp(self):
        super(OrderDistanceMatricesTests, self).setUp()

    def test_array_like_input(self):
        obs = _order_dms(self.minx, self.miny)
        self.assertEqual(obs, (self.minx_dm, self.miny_dm))

    def test_reordered_distance_matrices(self):
        # All matching IDs but with different orderings.
        x = self.minx_dm.filter(['1', '0', '2'])
        y = self.miny_dm.filter(['0', '2', '1'])

        exp = (x, y.filter(['1', '0', '2']))
        obs = _order_dms(x, y)
        self.assertEqual(obs, exp)

    def test_reordered_and_nonmatching_distance_matrices(self):
        # Some matching and nonmatching IDs, with different ordering.
        x = self.minx_dm_extra.filter(['1', '0', 'foo', '2'])
        z = self.minz_dm_extra.filter(['bar', '0', '2', '1'])

        exp = (x.filter(['1', '0', '2']), z.filter(['1', '0', '2']))
        obs = _order_dms(x, z, strict=False)
        self.assertEqual(obs, exp)

    def test_id_lookup(self):
        # Matrices have mismatched IDs but a lookup is provided.
        self.minx_dm_extra.ids = ['a', 'b', 'c', 'foo']
        self.minz_dm_extra.ids = ['d', 'e', 'f', 'bar']
        lookup = {'a': '0', 'b': '1', 'c': '2', 'foo': 'foo',
                  'd': '0', 'e': '1', 'f': '2', 'bar': 'bar'}

        x = self.minx_dm_extra.filter(['b', 'a', 'foo', 'c'])
        z = self.minz_dm_extra.filter(['bar', 'e', 'f', 'd'])

        x_copy = x.copy()
        z_copy = z.copy()

        exp = (self.minx_dm.filter(['1', '0', '2']),
               self.minz_dm.filter(['1', '0', '2']))
        obs = _order_dms(x, z, strict=False, lookup=lookup)
        self.assertEqual(obs, exp)

        # Make sure the inputs aren't modified.
        self.assertEqual(x, x_copy)
        self.assertEqual(z, z_copy)

    def test_lookup_with_array_like(self):
        lookup = {'0': 'a', '1': 'b', '2': 'c'}
        with self.assertRaises(ValueError):
            _order_dms(self.minx, self.miny, lookup=lookup)

    def test_shape_mismatch(self):
        with self.assertRaises(ValueError):
            _order_dms(self.minx, [[0, 2], [2, 0]])

    def test_missing_ids_in_lookup(self):
        # Mapping for '1' is missing. Should get an error while remapping IDs
        # for the first distance matrix.
        lookup = {'0': 'a', '2': 'c'}

        with self.assertRaisesRegex(KeyError, r"first.*(x).*'1'\"$"):
            _order_dms(self.minx_dm, self.miny_dm, lookup=lookup)

        # Mapping for 'bar' is missing. Should get an error while remapping IDs
        # for the second distance matrix.
        lookup = {'0': 'a', '1': 'b', '2': 'c',
                  'foo': 'a', 'baz': 'c'}
        self.miny_dm.ids = ('foo', 'bar', 'baz')

        with self.assertRaisesRegex(KeyError, r"second.*(y).*'bar'\"$"):
            _order_dms(self.minx_dm, self.miny_dm, lookup=lookup)

    def test_nonmatching_ids_strict_true(self):
        with self.assertRaises(ValueError):
            _order_dms(self.minx_dm, self.minz_dm_extra, strict=True)

    def test_no_matching_ids(self):
        self.minx_dm.ids = ['foo', 'bar', 'baz']
        self.miny_dm.ids = ['a', 'b', 'c']

        with self.assertRaises(ValueError):
            _order_dms(self.minx_dm, self.miny_dm, strict=False)

    def test_mixed_input_types(self):
        with self.assertRaises(TypeError):
            _order_dms(self.minx, self.minz_dm)

        with self.assertRaises(TypeError):
            _order_dms(self.minz_dm, self.minx)


if __name__ == '__main__':
    main()
