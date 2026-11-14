# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
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
        obs_r, obs_p, obs_perm_stats = hommola_cospeciation(
            self.hdist, self.pdist, self.interact, 9, seed=123)
        exp_p = 0.1
        exp_r = 0.8317096546324791
        exp_perm_stats = np.array([-0.234008859, 0.490495437, -0.116936751,
                                   -0.175562859, 0.189925198, 0.824737151,
                                   0.509368992, 0.016431258, -0.260609898])
        self.assertAlmostEqual(obs_p, exp_p)
        self.assertAlmostEqual(obs_r, exp_r)

        npt.assert_allclose(obs_perm_stats, exp_perm_stats)

    def test_hommola_cospeciation_asymmetric(self):
        obs_r, obs_p, obs_perm_stats = hommola_cospeciation(
            self.hdist_4x4, self.pdist, self.interact_5x4, 9, seed=123)
        exp_p = 0.1
        exp_r = 0.8573214099741122
        exp_perm_stats = np.array([-0.0985138, 0.05605763, -0.38783588,
                                   0.09342939, -0.28028816, -0.13080114,
                                   0.05910828, 0.4286607, 0.04082483])
        self.assertAlmostEqual(obs_p, exp_p)
        self.assertAlmostEqual(obs_r, exp_r)

        npt.assert_allclose(obs_perm_stats, exp_perm_stats)

    def test_hommola_cospeciation_no_sig(self):
        obs_r, obs_p, obs_perm_stats = hommola_cospeciation(
            self.hdist, self.pdist, self.interact_ns, 9, seed=123)
        exp_p = 0.5
        exp_r = -0.01367939137911453
        exp_perm_stats = np.array([0.49493401, -0.15686395, -0.40824829,
                                   -0.17739372, 0.23529593, -0.13187609,
                                   0.23529593, -0.03137279, 0.07664242])
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

    def test_equivalent_pairings_count_as_ties(self):
        # The one-to-one example in #810 has eight exact ties in these draws.
        # Integer covariance ranking gives no strictly greater correlations:
        # n*sum(x*y)-sum(x)*sum(y) is maximal for the observed pairing. Both
        # norms are invariant under these one-to-one permutations.
        _, p_value, _ = hommola_cospeciation(
            self.hdist, self.pdist, self.interact_1to1, 999, seed=0)
        self.assertEqual(p_value, 9 / 1000)

    def test_complete_interaction_permutations_are_all_tied(self):
        # Relabeling a complete bipartite interaction preserves the entire
        # multiset of paired distances. Every draw is therefore an exact tie,
        # regardless of last-bit differences in Pearson's correlation.
        for dtype in (np.float32, np.float64):
            with self.subTest(dtype=dtype):
                corr, p_value, perm_stats = hommola_cospeciation(
                    self.hdist.astype(dtype), self.pdist.astype(dtype),
                    np.ones((5, 5)), 99, seed=0)
                self.assertEqual(p_value, 1.0)
                npt.assert_allclose(perm_stats, corr, rtol=0,
                                    atol=np.finfo(dtype).eps)

    def test_zero_correlation_ties(self):
        # x=[1, 2, 3], y=[2, 1, 2] have exactly zero covariance. A
        # tolerance relative only to the near-zero observed r loses ties.
        host = [[0, 1, 2], [1, 0, 3], [2, 3, 0]]
        parasite = [[0, 2, 1], [2, 0, 2], [1, 2, 0]]
        _, p_value, _ = hommola_cospeciation(
            host, parasite, np.eye(3), 999, seed=0)
        self.assertEqual(p_value, 657 / 1000)

    def test_distinct_nearby_correlations_are_not_ties(self):
        # Break one symmetry with a small but resolved distance change.
        # Formerly tied correlations are separated by about 2.55e-10,
        # far above floating-point roundoff but below np.isclose's defaults.
        parasite = self.pdist * 100000000
        parasite[0, 4] += 1
        parasite[4, 0] += 1
        _, p_value, _ = hommola_cospeciation(
            self.hdist, parasite, self.interact_1to1, 999, seed=0)
        self.assertEqual(p_value, 6 / 1000)

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
