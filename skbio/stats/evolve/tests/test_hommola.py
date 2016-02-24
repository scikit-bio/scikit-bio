# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np
import numpy.testing as npt

from skbio.stats.distance import mantel
from skbio.stats.evolve import hommola_cospeciation
from skbio.stats.evolve._hommola import _get_dist, _gen_lists


class HommolaCospeciationTests(unittest.TestCase):
    def setUp(self):
        # Test matrices, as presented in original paper by Hommola et al.
        self.hdist = np.array([[0, 3, 8, 8, 9], [3, 0, 7, 7, 8], [
            8, 7, 0, 6, 7], [8, 7, 6, 0, 3], [9, 8, 7, 3, 0]])
        self.pdist = np.array([[0, 5, 8, 8, 8], [5, 0, 7, 7, 7], [
            8, 7, 0, 4, 4], [8, 7, 4, 0, 2], [8, 7, 4, 2, 0]])
        self.interact = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [
            0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 1, 1]])

        # Reduced-size host matrix for testing asymmetric interaction matrix
        self.hdist_4x4 = np.array([[0, 3, 8, 8], [3, 0, 7, 7], [8, 7, 0, 6],
                                  [8, 7, 6, 0]])
        self.interact_5x4 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0],
                                      [0, 0, 0, 1], [0, 0, 0, 1]])

        # One to one interaction matrix for comparing against Mantel output
        self.interact_1to1 = np.array([[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [
            0, 0, 1, 0, 0], [0, 0, 0, 1, 0], [0, 0, 0, 0, 1]])

        # interaction matrix yielding non-significant results.
        # this matrix was picked because it will generate an r value that's
        # less than a standard deviation away from the mean of the normal
        # distribution of r vals
        self.interact_ns = np.array(
            [[0, 0, 0, 1, 0], [0, 0, 0, 0, 1], [1, 0, 0, 0, 0],
             [1, 0, 0, 0, 0], [0, 0, 0, 0, 1]])

        # minimal size matrices for sanity checks of inputs
        self.h_dist_3x3 = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        self.h_dist_2x2 = np.array([[0, 3], [3, 0]])
        self.p_dist_3x3 = np.array([[0, 3, 2], [3, 0, 1], [2, 1, 0]])
        self.interact_3x3 = np.array([[0, 1, 1], [1, 0, 1], [0, 0, 1]])
        self.interact_3x2 = np.array([[0, 1], [1, 0], [1, 1]])
        self.interact_2x3 = np.array([[0, 1, 1], [1, 0, 1]])
        self.interact_zero = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

    def test_hommola_cospeciation_sig(self):
        np.random.seed(1)

        obs_r, obs_p, obs_perm_stats = hommola_cospeciation(
            self.hdist, self.pdist, self.interact, 9)
        exp_p = .1
        exp_r = 0.83170965463247915
        exp_perm_stats = np.array([-0.14928122, 0.26299538, -0.21125858,
                                   0.24143838, 0.61557855, -0.24258293,
                                   0.09885203, 0.02858, 0.42742399])
        self.assertAlmostEqual(obs_p, exp_p)
        self.assertAlmostEqual(obs_r, exp_r)

        npt.assert_allclose(obs_perm_stats, exp_perm_stats)

    def test_hommola_cospeciation_asymmetric(self):
        np.random.seed(1)

        obs_r, obs_p, obs_perm_stats = hommola_cospeciation(
            self.hdist_4x4, self.pdist, self.interact_5x4, 9)
        exp_p = 0.2
        exp_r = 0.85732140997411233
        exp_perm_stats = np.array([-0.315244162496, -0.039405520312,
                                   0.093429386594, -0.387835875941,
                                   0.183711730709,  0.056057631956,
                                   0.945732487487,  0.056057631956,
                                   -0.020412414523])
        self.assertAlmostEqual(obs_p, exp_p)
        self.assertAlmostEqual(obs_r, exp_r)

        npt.assert_allclose(obs_perm_stats, exp_perm_stats)

    def test_hommola_cospeciation_no_sig(self):
        np.random.seed(1)

        obs_r, obs_p, obs_perm_stats = hommola_cospeciation(
            self.hdist, self.pdist, self.interact_ns, 9)
        exp_p = .6
        exp_r = -0.013679391379114569
        exp_perm_stats = np.array([-0.22216543, -0.14836061, -0.04434843,
                                   0.1478281, -0.29105645, 0.56395839,
                                   0.47304992, 0.79125657, 0.06804138])
        self.assertAlmostEqual(obs_p, exp_p)
        self.assertAlmostEqual(obs_r, exp_r)
        npt.assert_allclose(obs_perm_stats, exp_perm_stats)

    def test_hommola_vs_mantel(self):
        # we don't compare p-values because the two methods use different
        # permutation strategies
        r_mantel, p_mantel, _ = mantel(
            self.hdist, self.pdist, method='pearson', permutations=0,
            alternative='greater'
        )
        r_hommola, p_hommola, _ = hommola_cospeciation(
            self.hdist, self.pdist, self.interact_1to1, permutations=0
        )

        self.assertAlmostEqual(r_hommola, r_mantel)
        npt.assert_equal(p_hommola, p_mantel)

    def test_zero_permutations(self):
        obs_r, obs_p, obs_perm_stats = hommola_cospeciation(
            self.hdist, self.pdist, self.interact, 0)

        exp_p = np.nan
        exp_r = 0.83170965463247915
        exp_perm_stats = np.array([])

        npt.assert_equal(obs_p, exp_p)
        self.assertAlmostEqual(obs_r, exp_r)
        npt.assert_equal(obs_perm_stats, exp_perm_stats)

    def test_get_dist(self):
        labels = np.array([0, 1, 1, 2, 3])
        k_labels, t_labels = _gen_lists(labels)
        dists = np.array([[0, 2, 6, 3], [2, 0, 5, 4], [6, 5, 0, 7],
                          [3, 4, 7, 0]])
        index = np.array([2, 3, 1, 0])

        expected_vec = np.array([7, 7, 5, 6, 0, 4, 3, 4, 3, 2])
        actual_vec = _get_dist(k_labels, t_labels, dists, index)

        npt.assert_allclose(actual_vec, expected_vec)

    def test_gen_lists(self):
        exp_pars_k_labels = np.array([0, 0, 0, 0, 0, 1, 1, 1,
                                      1, 2, 2, 2, 3, 3, 4])
        exp_pars_t_labels = np.array([1, 2, 3, 4, 4, 2, 3, 4,
                                      4, 3, 4, 4, 4, 4, 4])
        exp_host_k_labels = np.array([0, 0, 0, 0, 0, 1, 1, 1,
                                      1, 2, 2, 2, 3, 3, 3])
        exp_host_t_labels = np.array([1, 2, 3, 3, 4, 2, 3, 3,
                                      4, 3, 3, 4, 3, 4, 4])

        pars, hosts = np.nonzero(self.interact)

        obs_pars_k_labels, obs_pars_t_labels = _gen_lists(pars)
        obs_hosts_k_labels, obs_hosts_t_labels = _gen_lists(hosts)

        npt.assert_allclose(exp_pars_k_labels, obs_pars_k_labels)
        npt.assert_allclose(exp_pars_t_labels, obs_pars_t_labels)
        npt.assert_allclose(exp_host_k_labels, obs_hosts_k_labels)
        npt.assert_allclose(exp_host_t_labels, obs_hosts_t_labels)

    def test_dm_too_small(self):
        with self.assertRaises(ValueError):
            hommola_cospeciation(self.h_dist_2x2, self.p_dist_3x3,
                                 self.interact_3x3)

    def test_host_interaction_not_equal(self):
        with self.assertRaises(ValueError):
            hommola_cospeciation(self.h_dist_3x3, self.p_dist_3x3,
                                 self.interact_2x3)

    def test_par_interaction_not_equal(self):
        with self.assertRaises(ValueError):
            hommola_cospeciation(self.h_dist_3x3, self.p_dist_3x3,
                                 self.interact_3x2)

    def test_interaction_too_few(self):
        with self.assertRaises(ValueError):
            hommola_cospeciation(self.h_dist_3x3, self.p_dist_3x3,
                                 self.interact_zero)

    def test_permutations_too_few(self):
        with self.assertRaises(ValueError):
            hommola_cospeciation(self.h_dist_3x3, self.p_dist_3x3,
                                 self.interact_3x3, -1)


if __name__ == '__main__':
    unittest.main()
