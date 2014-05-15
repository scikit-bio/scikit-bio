#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------
from __future__ import division

from warnings import filterwarnings
from unittest import TestCase, main

import numpy as np

from skbio.math.stats.test import (G_2_by_2, g_fit, t_paired, t_one_sample,
                                   t_two_sample, mc_t_two_sample,
                                   _permute_observations, t_one_observation,
                                   correlation_t, ZeroExpectedError, fisher,
                                   safe_sum_p_log_p, permute_2d, mantel,
                                   mantel_t, _flatten_lower_triangle, pearson,
                                   spearman, ANOVA_one_way, mw_t,
                                   mw_boot, is_symmetric_and_hollow,
                                   reverse_tails, tail, fdr_correction,
                                   benjamini_hochberg_step_down,
                                   bonferroni_correction, fisher_z_transform,
                                   fisher_population_correlation,
                                   inverse_fisher_z_transform,
                                   z_transform_pval, kruskal_wallis, kendall,
                                   kendall_pval, assign_correlation_pval,
                                   cscore, williams_correction, g_stat)

from skbio.math.stats.distribution import chi_high, tprob


class TestsHelper(TestCase):

    """Class with utility methods useful for other tests."""

    # How many times a p-value should be tested to fall in a given range
    # before failing the test.
    p_val_tests = 20

    def assertCorrectPValue(self, exp_min, exp_max, fn, args=None,
                            kwargs=None, p_val_idx=0):
        """Tests that the stochastic p-value falls in the specified range.

        Performs the test self.p_val_tests times and fails if the observed
        p-value does not fall into the specified range at least once. Each
        p-value is also tested that it falls in the range 0.0 to 1.0.

        This method assumes that fn is callable, and will unpack and pass args
        and kwargs to fn if they are provided. It also assumes that fn returns
        a single value (the p-value to be tested) or a tuple of results (any
        length greater than or equal to 1), with the p-value at position
        p_val_idx.

        This is primarily used for testing the Mantel and correlation_test
        functions.
        """
        found_match = False
        for i in range(self.p_val_tests):
            if args is not None and kwargs is not None:
                obs = fn(*args, **kwargs)
            elif args is not None:
                obs = fn(*args)
            elif kwargs is not None:
                obs = fn(**kwargs)
            else:
                obs = fn()

            try:
                p_val = float(obs)
            except TypeError:
                p_val = obs[p_val_idx]

            self.assertTrue(0.0 <= p_val <= 1.0)
            if p_val >= exp_min and p_val <= exp_max:
                found_match = True
                break
        self.assertTrue(found_match)


class TestsTests(TestCase):

    """Tests miscellaneous functions."""

    def test_tail(self):
        """tail should return x/2 if test is true; 1-(x/2) otherwise"""
        np.testing.assert_allclose(tail(0.25, 'a' == 'a'), 0.25 / 2)
        np.testing.assert_allclose(tail(0.25, 'a' != 'a'), 1 - (0.25 / 2))

    def test_fisher(self):
        """fisher results should match p 795 Sokal and Rohlf"""
        np.testing.assert_allclose(fisher([0.073, 0.086, 0.10, 0.080, 0.060]),
                                   0.0045957946540917905, atol=10e-7)

    def test_permute_2d(self):
        """permute_2d permutes rows and cols of a matrix."""
        a = np.reshape(np.arange(9), (3, 3))
        np.testing.assert_allclose(permute_2d(a, [0, 1, 2]), a)
        np.testing.assert_allclose(permute_2d(a, [2, 1, 0]),
                                   np.array([[8, 7, 6], [5, 4, 3],
                                             [2, 1, 0]]))
        np.testing.assert_allclose(permute_2d(a, [1, 2, 0]),
                                   np.array([[4, 5, 3], [7, 8, 6],
                                             [1, 2, 0]]))


class GTests(TestCase):

    """Tests implementation of the G tests for fit and independence."""

    def test_G_2_by_2_2tailed_equal(self):
        """G_2_by_2 should return 0 if all cell counts are equal"""
        np.testing.assert_allclose(0, G_2_by_2(1, 1, 1, 1, False, False)[0])
        np.testing.assert_allclose(0, G_2_by_2(100, 100, 100, 100, False,
                                               False)[0])
        np.testing.assert_allclose(0, G_2_by_2(100, 100, 100, 100, True,
                                               False)[0])

    def test_G_2_by_2_bad_data(self):
        """G_2_by_2 should raise ValueError if any counts are negative"""
        self.assertRaises(ValueError, G_2_by_2, 1, -1, 1, 1)

    def test_G_2_by_2_2tailed_examples(self):
        """G_2_by_2 values should match examples in Sokal & Rohlf"""
        # example from p 731, Sokal and Rohlf (1995)
        # without correction
        np.testing.assert_allclose(G_2_by_2(12, 22, 16, 50, False, False)[0],
                                   1.33249, 0.0001)
        np.testing.assert_allclose(G_2_by_2(12, 22, 16, 50, False, False)[1],
                                   0.24836, 0.0001)
        # with correction
        np.testing.assert_allclose(G_2_by_2(12, 22, 16, 50, True, False)[0],
                                   1.30277, 0.0001)
        np.testing.assert_allclose(G_2_by_2(12, 22, 16, 50, True, False)[1],
                                   0.25371, 0.0001)

    def test_G_2_by_2_1tailed_examples(self):
        """G_2_by_2 values should match values from codon_binding program"""
        # first up...the famous arginine case
        np.testing.assert_allclose(G_2_by_2(36, 16, 38, 106), (29.111609, 0),
                                   atol=10e-7)
        # then some other miscellaneous positive and negative values
        np.testing.assert_allclose(
            G_2_by_2(0, 52, 12, 132), (-7.259930, 0.996474), atol=10e-7)
        np.testing.assert_allclose(
            G_2_by_2(5, 47, 14, 130), (-0.000481, 0.508751), atol=10e-7)
        np.testing.assert_allclose(
            G_2_by_2(5, 47, 36, 108), (-6.065167, 0.993106), atol=10e-7)

    def test_g_fit(self):
        """Test G fit is correct with and without Williams correction."""
        # test with williams correction
        data = [np.array(i) for i in [63, 31, 28, 12, 39, 16, 40, 12]]
        exp_G = 69.030858949133162 / 1.00622406639
        exp_p = 2.8277381487281706e-12
        obs_G, obs_p = g_fit(data, williams=True)
        np.testing.assert_allclose(obs_G, exp_G)
        np.testing.assert_allclose(obs_p, exp_p)
        # test with hand computed example and williams correction
        data = [np.array([75, 65, 48]), np.array([200]), np.array([10, 250, 13,
                                                                   85])]
        exp_G = 85.90859811005285 / 1.0018930430667
        exp_p = 2.4012235241479195e-19
        obs_G, obs_p = g_fit(data, williams=True)
        np.testing.assert_allclose(obs_G, exp_G)
        np.testing.assert_allclose(obs_p, exp_p)
        # test without williams correction on another hand computed example
        data = [np.array([10, 12, 15, 7]), np.array([15, 12, 17, 18]),
                np.array([6, 9, 13])]
        exp_G = 1.6610421781232
        exp_p = 0.43582212499949591
        obs_G, obs_p = g_fit(data, williams=False)
        np.testing.assert_allclose(obs_G, exp_G)
        np.testing.assert_allclose(obs_p, exp_p)

    def test_williams_correction(self):
        """Test that the Williams correction is correctly computed."""
        n = 100
        a = 10
        G = 10.5783
        exp = 10.387855973813421
        np.testing.assert_allclose(williams_correction(n, a, G), exp,
                                   rtol=1e-5)
        # test with an example from Sokal and Rohlf pg 699
        n = 241
        a = 8
        G = 8.82396
        exp = 8.76938
        np.testing.assert_allclose(williams_correction(n, a, G), exp,
                                   rtol=1e-5)

    def test_g_stat(self):
        """Test G-stat is correct when extrinsic hypothesis is equal freqs."""
        # test with equal len=1 vectors
        data = [np.array(i) for i in [63, 31, 28, 12, 39, 16, 40, 12]]
        exp = 69.030858949133162
        np.testing.assert_allclose(g_stat(data), exp, rtol=1e-5)
        # test with a hand computed example
        data = [np.array([75, 65, 48]), np.array([200]), np.array([10, 250, 13,
                                                                   85])]
        exp = 85.908598110052
        np.testing.assert_allclose(g_stat(data), exp, rtol=1e-5)

    def test_safe_sum_p_log_p(self):
        """safe_sum_p_log_p should ignore zero elements, not raise error"""
        m = np.array([2, 4, 0, 8])
        self.assertEqual(safe_sum_p_log_p(m, 2), 2 * 1 + 4 * 2 + 8 * 3)


class StatTests(TestsHelper):

    """Tests that the t and z tests are implemented correctly"""

    def setUp(self):
        super(StatTests, self).setUp()

        self.x = [
            7.33, 7.49, 7.27, 7.93, 7.56,
            7.81, 7.46, 6.94, 7.49, 7.44,
            7.95, 7.47, 7.04, 7.10, 7.64,
        ]

        self.y = [
            7.53, 7.70, 7.46, 8.21, 7.81,
            8.01, 7.72, 7.13, 7.68, 7.66,
            8.11, 7.66, 7.20, 7.25, 7.79,
        ]

        # silence warnings that will be raised by t_one_sample
        filterwarnings('ignore', category=RuntimeWarning)

    def test_t_paired_2tailed(self):
        """t_paired should match values from Sokal & Rohlf p 353"""
        x, y = self.x, self.y
        # check value of t and the probability for 2-tailed
        np.testing.assert_allclose(t_paired(y, x)[0], 19.7203, 1e-4)
        np.testing.assert_allclose(t_paired(y, x)[1], 1.301439e-11, 1e-4)

    def test_t_paired_no_variance(self):
        """t_paired should return None if lists are invariant"""
        x = [1, 1, 1]
        y = [0, 0, 0]
        self.assertEqual(t_paired(x, x), (None, None))
        self.assertEqual(t_paired(x, y), (None, None))

    def test_t_paired_1tailed(self):
        """t_paired should match pre-calculated 1-tailed values"""
        x, y = self.x, self.y
        # check probability for 1-tailed low and high
        np.testing.assert_allclose(
            t_paired(y, x, "low")[1], 1 - (1.301439e-11 / 2), 1e-4)
        np.testing.assert_allclose(
            t_paired(x, y, "high")[1], 1 - (1.301439e-11 / 2), 1e-4)
        np.testing.assert_allclose(
            t_paired(y, x, "high")[1], 1.301439e-11 / 2, 1e-4)
        np.testing.assert_allclose(
            t_paired(x, y, "low")[1], 1.301439e-11 / 2, 1e-4)

    def test_t_paired_specific_difference(self):
        """t_paired should allow a specific difference to be passed"""
        x, y = self.x, self.y
        # difference is 0.2, so test should be non-significant if 0.2 passed
        self.assertFalse(t_paired(y, x, exp_diff=0.2)[0] > 1e-10)
        # same, except that reversing list order reverses sign of difference
        self.assertFalse(t_paired(x, y, exp_diff=-0.2)[0] > 1e-10)
        # check that there's no significant difference from the true mean
        np.testing.assert_allclose(
            t_paired(y, x, exp_diff=0.2)[1], 1, 1e-4)

    def test_t_paired_bad_data(self):
        """t_paired should raise ValueError on lists of different lengths"""
        self.assertRaises(ValueError, t_paired, self.y, [1, 2, 3])

    def test_t_two_sample(self):
        """t_two_sample should match example on p.225 of Sokal and Rohlf"""
        I = np.array([7.2, 7.1, 9.1, 7.2, 7.3, 7.2, 7.5])
        II = np.array([8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2])
        np.testing.assert_allclose(t_two_sample(I, II),
                                   (-0.1184, 0.45385 * 2),
                                   atol=10e-3)

    def test_t_two_sample_no_variance(self):
        """t_two_sample should properly handle lists that are invariant"""
        # By default should return (None, None) to mimic R's t.test.
        x = np.array([1, 1., 1])
        y = np.array([0, 0, 0.0])
        self.assertEqual(t_two_sample(x, x), (None, None))
        self.assertEqual(t_two_sample(x, y), (None, None))

        # Test none_on_zero_variance=False on various tail types. We use
        # self.assertEqual instead of assert_allclose because the latter
        # sees inf and -inf as being equal.

        # Two tailed: a < b
        self.assertEqual(t_two_sample(y, x, none_on_zero_variance=False),
                         (float('-inf'), 0.0))

        # Two tailed: a > b
        self.assertEqual(t_two_sample(x, y, none_on_zero_variance=False),
                         (float('inf'), 0.0))

        # One-tailed 'high': a < b
        self.assertEqual(t_two_sample(y, x, tails='high',
                                      none_on_zero_variance=False),
                         (float('-inf'), 1.0))

        # One-tailed 'high': a > b
        self.assertEqual(t_two_sample(x, y, tails='high',
                                      none_on_zero_variance=False),
                         (float('inf'), 0.0))

        # One-tailed 'low': a < b
        self.assertEqual(t_two_sample(y, x, tails='low',
                                      none_on_zero_variance=False),
                         (float('-inf'), 0.0))

        # One-tailed 'low': a > b
        self.assertEqual(t_two_sample(x, y, tails='low',
                                      none_on_zero_variance=False),
                         (float('inf'), 1.0))

        # Should still receive (None, None) if the lists have no variance and
        # have the same single value.
        self.assertEqual(t_two_sample(x, x, none_on_zero_variance=False),
                         (None, None))
        self.assertEqual(t_two_sample(x, [1, 1], none_on_zero_variance=False),
                         (None, None))

    def test_t_two_sample_invalid_input(self):
        """t_two_sample should raise an error on invalid input."""
        self.assertRaises(ValueError, t_two_sample, [1, 2, 3], [4, 5, 6],
                          tails='foo')

    def test_t_one_sample(self):
        """t_one_sample results should match those from R"""
        x = np.array(range(-5, 5))
        y = np.array(range(-1, 10))
        np.testing.assert_allclose(t_one_sample(x), (-0.5222, 0.6141),
                                   atol=10e-3)
        np.testing.assert_allclose(t_one_sample(y), (4, 0.002518), atol=10e-3)
        # do some one-tailed tests as well
        np.testing.assert_allclose(
            t_one_sample(y, tails='low'), (4, 0.9987), atol=10e-3)
        np.testing.assert_allclose(
            t_one_sample(y, tails='high'), (4, 0.001259), atol=10e-3)

    def test_t_two_sample_switch(self):
        """t_two_sample should call t_one_observation if 1 item in sample."""
        sample = np.array([4.02, 3.88, 3.34, 3.87, 3.18])
        x = np.array([3.02])
        np.testing.assert_allclose(t_two_sample(x, sample),
                                   (-1.5637254, 0.1929248))
        np.testing.assert_allclose(t_two_sample(sample, x),
                                   (1.5637254, 0.1929248))

        # can't do the test if both samples have single item
        self.assertEqual(t_two_sample(x, x), (None, None))

        # Test special case if t=0.
        np.testing.assert_allclose(t_two_sample([2], [1, 2, 3]), (0.0, 1.0))
        np.testing.assert_allclose(t_two_sample([1, 2, 3], [2]), (0.0, 1.0))

    def test_t_one_observation(self):
        """t_one_observation should match p. 228 of Sokal and Rohlf"""
        sample = np.array([4.02, 3.88, 3.34, 3.87, 3.18])
        x = 3.02
        # note that this differs after the 3rd decimal place from what's in
        # the book, because Sokal and Rohlf round their intermediate steps...
        np.testing.assert_allclose(t_one_observation(x, sample),
                                   (-1.5637254, 0.1929248))

    def test_t_one_observation_no_variance(self):
        """t_one_observation should correctly handle an invariant list."""
        sample = np.array([1.0, 1.0, 1.0])

        # Can't perform test if invariant list's single value matches x,
        # regardless of none_on_zero_variance.
        self.assertEqual(t_one_observation(1, sample), (None, None))
        self.assertEqual(t_one_observation(1, sample,
                                           none_on_zero_variance=False),
                         (None, None))

        # Test correct handling of none_on_zero_variance.
        self.assertEqual(t_one_observation(2, sample), (None, None))
        self.assertEqual(t_one_observation(2, sample,
                                           none_on_zero_variance=False),
                         (float('inf'), 0.0))
        self.assertEqual(t_one_observation(2, sample,
                                           none_on_zero_variance=False,
                                           tails='low'),
                         (float('inf'), 1.0))

    def test_mc_t_two_sample(self):
        """Test gives correct results with valid input data."""
        # Verified against R's t.test() and Deducer::perm.t.test().

        # With numpy array as input.
        exp = (-0.11858541225631833, 0.90756579317867436)
        I = np.array([7.2, 7.1, 9.1, 7.2, 7.3, 7.2, 7.5])
        II = np.array([8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2])
        obs = mc_t_two_sample(I, II)
        np.testing.assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 999)
        self.assertCorrectPValue(0.8, 0.9, mc_t_two_sample, [I, II],
                                 p_val_idx=3)

        # With python list as input.
        exp = (-0.11858541225631833, 0.90756579317867436)
        I = [7.2, 7.1, 9.1, 7.2, 7.3, 7.2, 7.5]
        II = [8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2]
        obs = mc_t_two_sample(I, II)
        np.testing.assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 999)
        self.assertCorrectPValue(0.8, 0.9, mc_t_two_sample, [I, II],
                                 p_val_idx=3)

        exp = (-0.11858541225631833, 0.45378289658933718)
        obs = mc_t_two_sample(I, II, tails='low')
        np.testing.assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 999)
        self.assertCorrectPValue(0.4, 0.47, mc_t_two_sample, [I, II],
                                 {'tails': 'low'}, p_val_idx=3)

        exp = (-0.11858541225631833, 0.54621710341066287)
        obs = mc_t_two_sample(I, II, tails='high', permutations=99)
        np.testing.assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 99)
        self.assertCorrectPValue(0.4, 0.62, mc_t_two_sample, [I, II],
                                 {'tails': 'high', 'permutations': 99},
                                 p_val_idx=3)

        exp = (-2.8855783649036986, 0.99315596652421401)
        obs = mc_t_two_sample(I, II, tails='high',
                              permutations=99, exp_diff=1)
        np.testing.assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 99)
        self.assertCorrectPValue(0.55, 0.99, mc_t_two_sample, [I, II],
                                 {'tails': 'high', 'permutations': 99,
                                  'exp_diff': 1}, p_val_idx=3)

    def test_mc_t_two_sample_unbalanced_obs(self):
        """Test gives correct results with unequal number of obs per sample."""
        # Verified against R's t.test() and Deducer::perm.t.test().
        exp = (-0.10302479888889175, 0.91979753020527177)
        I = np.array([7.2, 7.1, 9.1, 7.2, 7.3, 7.2])
        II = np.array([8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2])
        obs = mc_t_two_sample(I, II)
        np.testing.assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 999)
        self.assertCorrectPValue(0.8, 0.9, mc_t_two_sample, [I, II],
                                 p_val_idx=3)

    def test_mc_t_two_sample_single_obs_sample(self):
        """Test works correctly with one sample having a single observation."""
        sample = np.array([4.02, 3.88, 3.34, 3.87, 3.18])
        x = np.array([3.02])
        exp = (-1.5637254, 0.1929248)
        obs = mc_t_two_sample(x, sample)
        np.testing.assert_allclose(obs[:2], exp)
        np.testing.assert_allclose(len(obs[2]), 999)
        self.assertTrue(0.0 <= obs[3] <= 1.0)

        exp = (1.5637254, 0.1929248)
        obs = mc_t_two_sample(sample, x)
        np.testing.assert_allclose(obs[:2], exp)
        np.testing.assert_allclose(len(obs[2]), 999)
        self.assertTrue(0.0 <= obs[3] <= 1.0)

        # Test the case where we can have no variance in the permuted lists.
        x = np.array([1, 1, 2])
        y = np.array([1])
        exp = (0.5, 0.666666666667)
        obs = mc_t_two_sample(x, y)
        np.testing.assert_allclose(obs[:2], exp)
        np.testing.assert_allclose(len(obs[2]), 999)
        self.assertTrue(0.0 <= obs[3] <= 1.0)

    def test_mc_t_two_sample_no_perms(self):
        """Test gives empty permutation results if no perms are given."""
        exp = (-0.11858541225631833, 0.90756579317867436, [], np.nan)
        I = np.array([7.2, 7.1, 9.1, 7.2, 7.3, 7.2, 7.5])
        II = np.array([8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2])
        obs = mc_t_two_sample(I, II, permutations=0)
        np.testing.assert_allclose(obs[0], exp[0])
        np.testing.assert_allclose(obs[1], exp[1])
        self.assertEqual(obs[2], exp[2])
        np.testing.assert_allclose(obs[3], exp[3])

    def test_mc_t_two_sample_no_mc(self):
        """Test no MC stats if initial t-test is bad."""
        x = np.array([1, 1, 1])
        y = np.array([0, 0, 0])
        self.assertEqual(mc_t_two_sample(x, x), (None, None, [], np.nan))

    def test_mc_t_two_sample_no_variance(self):
        """Test input with no variance. Should match Deducer::perm.t.test."""
        x = np.array([1, 1, 1])
        y = np.array([2, 2, 2])

        exp = (float('-inf'), 0.0)
        obs = mc_t_two_sample(x, y, permutations=1000)

        self.assertEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 1000)
        self.assertCorrectPValue(0.09, 0.11, mc_t_two_sample, [x, y],
                                 {'permutations': 1000}, p_val_idx=3)

        exp = (float('inf'), 0.0)
        obs = mc_t_two_sample(y, x, permutations=1000)

        self.assertEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 1000)
        self.assertCorrectPValue(0.09, 0.11, mc_t_two_sample, [y, x],
                                 {'permutations': 1000}, p_val_idx=3)

        exp = (float('-inf'), 1.0)
        obs = mc_t_two_sample(x, y, permutations=1000, tails='high')

        self.assertEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 1000)
        self.assertCorrectPValue(0.9999, 1.0, mc_t_two_sample, [x, y],
                                 {'permutations': 1000, 'tails': 'high'},
                                 p_val_idx=3)

        exp = (float('-inf'), 0.0)
        obs = mc_t_two_sample(x, y, permutations=1000, tails='low')

        self.assertEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 1000)
        self.assertCorrectPValue(0.04, 0.051, mc_t_two_sample, [x, y],
                                 {'permutations': 1000, 'tails': 'low'},
                                 p_val_idx=3)

    def test_mc_t_two_sample_no_permuted_variance(self):
        """Test with chance of getting no variance with some perms."""
        # Verified against R's t.test() and Deducer::perm.t.test().
        x = np.array([1, 1, 2])
        y = np.array([2, 2, 1])

        exp = (-0.70710678118654791, 0.51851851851851838)
        obs = mc_t_two_sample(x, y, permutations=1000)

        np.testing.assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 1000)
        self.assertCorrectPValue(0.97, 1.0, mc_t_two_sample, [x, y],
                                 {'permutations': 1000}, p_val_idx=3)

    def test_mc_t_two_sample_invalid_input(self):
        """Test fails on various invalid input."""
        self.assertRaises(ValueError, mc_t_two_sample, [1, 2, 3],
                          [4., 5., 4.], tails='foo')
        self.assertRaises(ValueError, mc_t_two_sample, [1, 2, 3],
                          [4., 5., 4.], permutations=-1)
        self.assertRaises(ValueError, mc_t_two_sample, [1], [4.])
        self.assertRaises(ValueError, mc_t_two_sample, [1, 2], [])

    def test_permute_observations(self):
        """Test works correctly on small input dataset."""
        I = [10, 20., 1]
        II = [2, 4, 5, 7]
        obs = _permute_observations(I, II, 1)
        self.assertEqual(len(obs[0]), 1)
        self.assertEqual(len(obs[1]), 1)
        self.assertEqual(len(obs[0][0]), len(I))
        self.assertEqual(len(obs[1][0]), len(II))
        np.testing.assert_allclose(sorted(np.concatenate((obs[0][0],
                                                          obs[1][0]))),
                                   sorted(I + II))

    def test_reverse_tails(self):
        """reverse_tails should return 'high' if tails was 'low' or vice versa
        """
        self.assertEqual(reverse_tails('high'), 'low')
        self.assertEqual(reverse_tails('low'), 'high')
        self.assertEqual(reverse_tails(None), None)
        self.assertEqual(reverse_tails(3), 3)

    def test_tail(self):
        """tail should return prob/2 if test is true, or 1-(prob/2) if false
        """
        np.testing.assert_allclose(tail(0.25, True), 0.125)
        np.testing.assert_allclose(tail(0.25, False), 0.875)
        np.testing.assert_allclose(tail(1, True), 0.5)
        np.testing.assert_allclose(tail(1, False), 0.5)
        np.testing.assert_allclose(tail(0, True), 0)
        np.testing.assert_allclose(tail(0, False), 1)


class CorrelationTests(TestsHelper):

    """Tests of correlation coefficients and Mantel test."""

    def setUp(self):
        """Sets up variables used in the tests."""
        super(CorrelationTests, self).setUp()

        # For testing spearman and correlation_test using method='spearman'.
        # Taken from the Spearman wikipedia article. Also used for testing
        # Pearson (verified with R).
        self.data1 = [106, 86, 100, 101, 99, 103, 97, 113, 112, 110]
        self.data2 = [7, 0, 27, 50, 28, 29, 20, 12, 6, 17]

        # For testing spearman.
        self.a = [1, 2, 4, 3, 1, 6, 7, 8, 10, 4]
        self.b = [2, 10, 20, 1, 3, 7, 5, 11, 6, 13]
        self.c = [7, 1, 20, 13, 3, 57, 5, 121, 2, 9]
        self.r = (1.7, 10, 20, 1.7, 3, 7, 5, 11, 6.5, 13)
        self.x = (1, 2, 4, 3, 1, 6, 7, 8, 10, 4, 100, 2, 3, 77)

        # Ranked copies for testing spearman.
        self.b_ranked = [2, 7, 10, 1, 3, 6, 4, 8, 5, 9]
        self.c_ranked = [5, 1, 8, 7, 3, 9, 4, 10, 2, 6]

        # silence the warnings that will tests for correlation_test
        filterwarnings('ignore', category=RuntimeWarning)

    def test_mantel(self):
        """mantel should be significant for same matrix, not for random"""
        a = np.reshape(np.arange(25), (5, 5))
        a = np.tril(a) + np.tril(a).T
        np.fill_diagonal(a, 0)
        b = a.copy()
        # closely related -- should be significant
        self.assertCorrectPValue(0.0, 0.049, mantel, (a, b, 1000))

        c = np.reshape(np.ones(25), (5, 5))
        c[0, 1] = 3.0
        c[1, 0] = 3.0
        np.fill_diagonal(c, 0)
        # not related -- should not be significant
        self.assertCorrectPValue(0.06, 1.0, mantel, (a, c, 1000))

    def test_mantel_test_one_sided_greater(self):
        """Test one-sided mantel test (greater)."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test).
        m1 = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        m2 = np.array([[0, 2, 7], [2, 0, 6], [7, 6, 0]])
        p, stat, perms = mantel_t(m1, m1, 999, alt='greater')
        np.testing.assert_allclose(stat, 1.0)
        self.assertEqual(len(perms), 999)

        self.assertCorrectPValue(0.09, 0.25, mantel_t, (m1, m1, 999),
                                 {'alt': 'greater'})

        p, stat, perms = mantel_t(m1, m2, 999, alt='greater')
        np.testing.assert_allclose(stat, 0.755928946018)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.2, 0.5, mantel_t, (m1, m2, 999),
                                 {'alt': 'greater'})

    def test_mantel_test_one_sided_less(self):
        """Test one-sided mantel test (less)."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test, but I modified their output to do a one-sided
        # less test).
        m1 = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        m2 = np.array([[0, 2, 7], [2, 0, 6], [7, 6, 0]])
        m3 = np.array([[0, 0.5, 0.25], [0.5, 0, 0.1], [0.25, 0.1, 0]])
        p, stat, perms = mantel_t(m1, m1, 999, alt='less')
        np.testing.assert_allclose(p, 1.0)
        np.testing.assert_allclose(stat, 1.0)
        self.assertEqual(len(perms), 999)

        p, stat, perms = mantel_t(m1, m2, 999, alt='less')
        np.testing.assert_allclose(stat, 0.755928946018)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.6, 1.0, mantel_t, (m1, m2, 999),
                                 {'alt': 'less'})

        p, stat, perms = mantel_t(m1, m3, 999, alt='less')
        np.testing.assert_allclose(stat, -0.989743318611)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.1, 0.25, mantel_t, (m1, m3, 999),
                                 {'alt': 'less'})

    def test_mantel_test_two_sided(self):
        """Test two-sided mantel test."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test, but I modified their output to do a two-sided
        # test).
        m1 = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        m2 = np.array([[0, 2, 7], [2, 0, 6], [7, 6, 0]])
        m3 = np.array([[0, 0.5, 0.25], [0.5, 0, 0.1], [0.25, 0.1, 0]])
        p, stat, perms = mantel_t(m1, m1, 999, alt='two sided')
        np.testing.assert_allclose(stat, 1.0)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.20, 0.45, mantel_t, (m1, m1, 999),
                                 {'alt': 'two sided'})

        p, stat, perms = mantel_t(m1, m2, 999, alt='two sided')
        np.testing.assert_allclose(stat, 0.755928946018)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.6, 0.75, mantel_t, (m1, m2, 999),
                                 {'alt': 'two sided'})

        p, stat, perms = mantel_t(m1, m3, 999, alt='two sided')
        np.testing.assert_allclose(stat, -0.989743318611)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.2, 0.45, mantel_t, (m1, m3, 999),
                                 {'alt': 'two sided'})

    def test_mantel_test_invalid_distance_matrix(self):
        """Test mantel test with invalid distance matrix."""
        # Single asymmetric, non-hollow distance matrix.
        self.assertRaises(ValueError, mantel_t, np.array([[1, 2], [3, 4]]),
                          np.array([[0, 0], [0, 0]]), 999)

        # Two asymmetric distance matrices.
        self.assertRaises(ValueError, mantel_t, np.array([[0, 2], [3, 0]]),
                          np.array([[0, 1], [0, 0]]), 999)

    def test_mantel_test_invalid_input(self):
        """Test mantel test with invalid input."""
        self.assertRaises(ValueError, mantel_t, np.array([[1]]),
                          np.array([[1]]), 999, alt='foo')
        self.assertRaises(ValueError, mantel_t, np.array([[1]]),
                          np.array([[1, 2], [3, 4]]), 999)
        self.assertRaises(ValueError, mantel_t, np.array([[1]]),
                          np.array([[1]]), 0)
        self.assertRaises(ValueError, mantel_t, np.array([[1]]),
                          np.array([[1]]), -1)

    def test_is_symmetric_and_hollow(self):
        """Should correctly test for symmetry and hollowness of dist mats."""
        self.assertTrue(is_symmetric_and_hollow(np.array([[0, 1], [1, 0]])))
        self.assertTrue(is_symmetric_and_hollow(np.matrix([[0, 1], [1, 0]])))
        self.assertTrue(is_symmetric_and_hollow(np.matrix([[0.0, 0],
                                                           [0.0, 0]])))
        self.assertTrue(not is_symmetric_and_hollow(
            np.array([[0.001, 1], [1, 0]])))
        self.assertTrue(not is_symmetric_and_hollow(
            np.array([[0, 1.1], [1, 0]])))
        self.assertTrue(not is_symmetric_and_hollow(
            np.array([[0.5, 1.1], [1, 0]])))

    def test_flatten_lower_triangle(self):
        """Test flattening various dms' lower triangulars."""
        self.assertEqual(_flatten_lower_triangle(
            np.array([[8]])), [])
        self.assertEqual(_flatten_lower_triangle(
            np.array([[1, 2], [3, 4]])), [3])
        self.assertEqual(_flatten_lower_triangle(
            np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])), [4, 7, 8])

    def test_pearson(self):
        """Test pearson correlation method on valid data."""
        # This test output was verified by R.
        np.testing.assert_allclose(pearson([1, 2], [1, 2]), 1.0)
        np.testing.assert_allclose(pearson([1, 2, 3], [1, 2, 3]), 1.0)
        np.testing.assert_allclose(pearson([1, 2, 3], [1, 2, 4]), 0.9819805)

    def test_pearson_invalid_input(self):
        """Test running pearson on bad input."""
        self.assertRaises(ValueError, pearson, [1.4, 2.5], [5.6, 8.8, 9.0])
        self.assertRaises(ValueError, pearson, [1.4], [5.6])

    def test_spearman(self):
        """Test the spearman function with valid input."""
        # One vector has no ties.
        exp = 0.3719581
        obs = spearman(self.a, self.b)
        np.testing.assert_allclose(obs, exp)

        # Both vectors have no ties.
        exp = 0.2969697
        obs = spearman(self.b, self.c)
        np.testing.assert_allclose(obs, exp)

        # Both vectors have ties.
        exp = 0.388381
        obs = spearman(self.a, self.r)
        np.testing.assert_allclose(obs, exp)

        exp = -0.17575757575757578
        obs = spearman(self.data1, self.data2)
        np.testing.assert_allclose(obs, exp)

    def test_spearman_no_variation(self):
        """Test the spearman function with a vector having no variation."""
        exp = np.nan
        obs = spearman([1, 1, 1], [1, 2, 3])
        np.testing.assert_allclose(obs, exp)

    def test_spearman_ranked(self):
        """Test the spearman function with a vector that is already ranked."""
        exp = 0.2969697
        obs = spearman(self.b_ranked, self.c_ranked)
        np.testing.assert_allclose(obs, exp)

    def test_spearman_one_obs(self):
        """Test running spearman on a single observation."""
        self.assertRaises(ValueError, spearman, [1.0], [5.0])

    def test_spearman_invalid_input(self):
        """Test the spearman function with invalid input."""
        self.assertRaises(ValueError, spearman, [], [])
        self.assertRaises(ValueError, spearman, self.a, [])

    def test_correlation_test_pearson(self):
        """Test correlation_t using pearson on valid input."""
        # These results were verified with R.
        # Test with non-default confidence level and permutations.
        obs = correlation_t(self.data1, self.data2, method='pearson',
                            confidence_level=0.90, permutations=990)
        np.testing.assert_allclose(obs[:2], (-0.03760147,
                                             0.91786297277172868), atol=10e-7)
        self.assertEqual(len(obs[2]), 990)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.9, 0.93, correlation_t,
                                 (self.data1, self.data2),
                                 {'method': 'pearson',
                                  'confidence_level': 0.90,
                                  'permutations': 990},
                                 p_val_idx=3)
        np.testing.assert_allclose(obs[4], (-0.5779077, 0.5256224))

        # Test with non-default tail type.
        obs = correlation_t(self.data1, self.data2, method='pearson',
                            confidence_level=0.90, permutations=990,
                            tails='low')
        np.testing.assert_allclose(obs[:2], (-0.03760147,
                                             0.45893148638586434), atol=10e-7)
        self.assertEqual(len(obs[2]), 990)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.41, 0.46, correlation_t,
                                 (self.data1, self.data2),
                                 {'method': 'pearson',
                                  'confidence_level': 0.90,
                                  'permutations': 990,
                                  'tails': 'low'},
                                 p_val_idx=3)
        np.testing.assert_allclose(obs[4], (-0.5779077, 0.5256224))

    def test_correlation_test_spearman(self):
        """Test correlation_t using spearman on valid input."""
        # This example taken from Wikipedia page:
        # http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient
        obs = correlation_t(self.data1, self.data2, method='spearman',
                            tails='high')
        np.testing.assert_allclose(obs[:2], (-0.17575757575757578,
                                             0.686405827612))
        self.assertEqual(len(obs[2]), 999)
        for rho in obs[2]:
            self.assertTrue(rho >= -1.0 and rho <= 1.0)
        self.assertCorrectPValue(0.67, 0.7, correlation_t,
                                 (self.data1, self.data2),
                                 {'method': 'spearman',
                                  'tails': 'high'},
                                 p_val_idx=3)
        np.testing.assert_allclose(obs[4], (-0.7251388558041697,
                                            0.51034422964834503))

        # The p-value is off because the example uses a one-tailed test, while
        # we use a two-tailed test. Someone confirms the answer that we get
        # here for a two-tailed test:
        # http://stats.stackexchange.com/questions/22816/calculating-p-value-
        #     for-spearmans-rank-correlation-coefficient-example-on-wikip
        obs = correlation_t(self.data1, self.data2, method='spearman',
                            tails=None)
        np.testing.assert_allclose(obs[:2], (-0.17575757575757578,
                                             0.62718834477648433))
        self.assertEqual(len(obs[2]), 999)
        for rho in obs[2]:
            self.assertTrue(rho >= -1.0 and rho <= 1.0)
        self.assertCorrectPValue(0.60, 0.64, correlation_t,
                                 (self.data1, self.data2),
                                 {'method': 'spearman', 'tails': None},
                                 p_val_idx=3)
        np.testing.assert_allclose(obs[4], (-0.7251388558041697,
                                            0.51034422964834503))

    def test_correlation_test_invalid_input(self):
        """Test correlation_t using invalid input."""
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          method='foo')
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          tails='foo')
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          permutations=-1)
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          confidence_level=-1)
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          confidence_level=1.1)
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          confidence_level=0)
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          confidence_level=0.0)
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          confidence_level=1)
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          confidence_level=1.0)

    def test_correlation_test_no_permutations(self):
        """Test correlation_t with no permutations."""
        # These results were verified with R.
        exp = (-0.2581988897471611, 0.7418011102528389, [], None,
               (-0.97687328610475876, 0.93488023560400879))
        obs = correlation_t([1, 2, 3, 4], [1, 2, 1, 1], permutations=0)
        np.testing.assert_allclose(obs[0], exp[0])
        np.testing.assert_allclose(obs[1], exp[1])
        np.testing.assert_allclose(obs[2], exp[2])
        self.assertEqual(obs[3], exp[3])
        np.testing.assert_allclose(obs[4], exp[4])

    def test_correlation_test_perfect_correlation(self):
        """Test correlation_t with perfectly-correlated input vectors."""
        # These results were verified with R.
        obs = correlation_t([1, 2, 3, 4], [1, 2, 3, 4])
        np.testing.assert_allclose(obs[:2], (1.0, 0.0))
        self.assertEqual(len(obs[2]), 999)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.06, 0.09, correlation_t,
                                 ([1, 2, 3, 4], [1, 2, 3, 4]),
                                 p_val_idx=3)
        np.testing.assert_allclose(obs[4], (0.99999999999998879, 1.0))

    def test_correlation_test_small_obs(self):
        """Test correlation_t with a small number of observations."""
        # These results were verified with R.
        obs = correlation_t([1, 2, 3], [1, 2, 3])
        np.testing.assert_allclose(obs[:2], (1.0, 0))
        self.assertEqual(len(obs[2]), 999)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.3, 0.4, correlation_t,
                                 ([1, 2, 3], [1, 2, 3]),
                                 p_val_idx=3)
        self.assertEqual(obs[4], (None, None))

        obs = correlation_t([1, 2, 3], [1, 2, 3], method='spearman')
        np.testing.assert_allclose(obs[:2], (1.0, 0))
        self.assertEqual(len(obs[2]), 999)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.3, 0.4, correlation_t,
                                 ([1, 2, 3], [1, 2, 3]),
                                 {'method': 'spearman'}, p_val_idx=3)
        self.assertEqual(obs[4], (None, None))

    def test_mw_test(self):
        """mann-whitney test results should match Sokal & Rohlf"""
        # using Sokal and Rolhf and R wilcox.test
        # x <- c(104, 109, 112, 114, 116, 118, 118, 119, 121, 123, 125, 126,
        #        126, 128, 128, 128)
        # y <- c(100, 105, 107, 107, 108, 111, 116, 120, 121, 123)
        # wilcox.test(x,y)
        # W = 123.5, p-value = 0.0232
        x = [104, 109, 112, 114, 116, 118, 118, 119, 121, 123, 125, 126, 126,
             128, 128, 128]
        y = [100, 105, 107, 107, 108, 111, 116, 120, 121, 123]
        u, p = mw_t(x, y, continuity=True, two_sided=True)
        # a return of 123.5 would also be okay, there is a consensus to use the
        # smaller U statistic, but the probability calculated from each is the
        # same
        self.assertTrue(u == 36.5 or u == 123.5)
        np.testing.assert_allclose(p, .0232, rtol=1e-3)

    def test_mw_boot(self):
        """excercising the Monte-carlo variant of mann-whitney"""
        x = [104, 109, 112, 114, 116, 118, 118, 119, 121, 123, 125, 126, 126,
             128, 128, 128]
        y = [100, 105, 107, 107, 108, 111, 116, 120, 121, 123]
        u, p = mw_boot(x, y, 10)
        self.assertTrue(u == 36.5 or u == 123.5)
        self.assertTrue(0 <= p <= 0.5)

    def test_kendall(self):
        """tests new kendall tau implamentation, returns tau, prob"""
        # test from pg. 594 Sokal and Rohlf, Box 15.7
        v1 = [8.7, 8.5, 9.4, 10, 6.3, 7.8, 11.9, 6.5, 6.6, 10.6, 10.2, 7.2,
              8.6, 11.1, 11.6]
        v2 = [5.95, 5.65, 6.00, 5.70, 4.70, 5.53, 6.40, 4.18, 6.15, 5.93, 5.70,
              5.68, 6.13, 6.30, 6.03]
        obs_tau = kendall(v1, v2)
        obs_prob = kendall_pval(obs_tau, len(v1))
        exp_tau = 0.49761335152811925
        exp_prob = 0.0097188572446995618
        np.testing.assert_allclose(obs_tau, exp_tau)
        np.testing.assert_allclose(obs_prob, exp_prob)
        # random vectors checked against scipy. v1 has 33 ties, v2 32
        v1 = np.array(
            [1.2, 9.7, 8.8, 1.7, 8.6, 9.9, 6.8, 7.3, 5.5, 5.4, 8.3,
             3.6, 7.5, 2., 9.3, 5.1, 8.4, 0.3, 8.2, 2.4, 9.8, 8.5,
             2.1, 6., 1.8, 3.7, 1.4, 4.6, 7.6, 5.2, 0.9, 5.2, 4.7,
             2.9, 5., 6.9, 1.3, 6.7, 5.2, 2.4, 6.9, 2., 7.4, 0.4,
             8.2, 9.5, 2.9, 5.7, 2.4, 8.8, 1.6, 3.5, 5.1, 3.6, 3.3,
             7.5, 0.9, 9.3, 5.4, 6.9, 9.3, 2.3, 1.9, 8.1, 3.2, 4.2,
             8.7, 3., 9.8, 5.3, 6.2, 4.8, 9., 2.8, 5.5, 8.4, 4.1,
             5.6, 5.4, 6.9, 3.8, 2.7, 0.3, 3.9, 8.2, 6.6, 1.9, 3.9,
             2., 4.4, 0.8, 6.5, 4.8, 1.5, 9.9, 9.1, 9.9, 6.2, 2.9,
             2.])
        v2 = np.array([6.6, 8.6, 3.9, 6.1, 0.9, 8.4, 10., 3.3, 0.4,
                       3.9, 7.6, 8.2, 8.6, 3., 6.9, 0.6, 8.4, 8.1,
                       6.3, 0.5, 5.2, 6.4, 8., 9.9, 1.2, 6.7, 8.4,
                       2.7, 8.4, 4.1, 4.6, 5.1, 5.2, 5.3, 2.2, 2.2,
                       4.3, 7.1, 1.4, 6.6, 7.6, 4.5, 7.8, 3.5, 7.1,
                       0.6, 4.6, 3.2, 2.2, 0.2, 3.9, 5.9, 7.7, 8.8,
                       1.3, 5.1, 5.6, 8.3, 8.8, 1.7, 5.2, 6.9, 1.3,
                       1.4, 4.9, 9.4, 2.3, 3.7, 9.1, 3.4, 1.6, 4.1,
                       9.7, 2.8, 9.9, 0.5, 2., 2.7, 3.3, 2.4, 3.6,
                       7.9, 6.5, 7., 4.2, 1.8, 1.6, 1.9, 5.5, 0.,
                       1.4, 2.2, 7.2, 8.2, 1.1, 2.5, 5.3, 0.2, 9., 0.2])
        exp_tau, exp_prob = (0.024867511238807951, 0.71392573687923555)
        obs_tau = kendall(v1, v2)
        obs_prob = kendall_pval(obs_tau, len(v1))
        np.testing.assert_allclose(obs_tau, exp_tau)
        np.testing.assert_allclose(obs_prob, exp_prob)


class TestDistMatrixPermutationTest(TestCase):

    """Tests of distance_matrix_permutation_test"""

    def setUp(self):
        """sets up variables for testing"""
        self.matrix = np.array(
            [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
        self.cells = [(0, 1), (1, 3)]
        self.cells2 = [(0, 2), (2, 3)]

    def test_ANOVA_one_way(self):
        """ANOVA one way returns same values as ANOVA on a stats package
        """
        g1 = np.array([10.0, 11.0, 10.0, 5.0, 6.0])
        g2 = np.array([1.0, 2.0, 3.0, 4.0, 1.0, 2.0])
        g3 = np.array([6.0, 7.0, 5.0, 6.0, 7.0])
        i = [g1, g2, g3]
        F, pval = ANOVA_one_way(i)

        np.testing.assert_allclose(F, 18.565450643776831)
        np.testing.assert_allclose(pval, 0.00015486238993089464)

    def test_kruskal_wallis(self):
        """Test kruskal_wallis on Sokal & Rohlf Box 13.6 dataset"""
        d_control = [75, 67, 70, 75, 65, 71, 67, 67, 76, 68]
        d_2_gluc = [57, 58, 60, 59, 62, 60, 60, 57, 59, 61]
        d_2_fruc = [58, 61, 56, 58, 57, 56, 61, 60, 57, 58]
        d_1_1 = [58, 59, 58, 61, 57, 56, 58, 57, 57, 59]
        d_2_sucr = [62, 66, 65, 63, 64, 62, 65, 65, 62, 67]
        data = [d_control, d_2_gluc, d_2_fruc, d_1_1, d_2_sucr]
        kw_stat, pval = kruskal_wallis(data)
        np.testing.assert_allclose(kw_stat, 38.436807439)
        np.testing.assert_allclose(pval, 9.105424085598766e-08)
        # test using a random data set against scipy
        x_0 = np.array([0, 0, 0, 31, 12, 0, 25, 26, 775, 13])
        x_1 = np.array([14, 15, 0, 15, 12, 13])
        x_2 = np.array([0, 0, 0, 55, 92, 11, 11, 11, 555])
        # kruskal(x_0, x_1, x_2) = (0.10761259465923653, 0.94761564440615031)
        exp = (0.10761259465923653, 0.94761564440615031)
        obs = kruskal_wallis([x_0, x_1, x_2])
        np.testing.assert_allclose(obs, exp)


class PvalueTests(TestCase):

    '''Test that the methods for handling Pvalues return the results we expect.

    Note: eps is being set lower on some of these because Sokal and Rohlf
    provide only ~5 sig figs and our integrals diverge by that much or more.
    '''

    def setUp(self):
        '''Nothing needed for all tests.'''
        pass

    def test_fdr_correction(self):
        """Test that the fdr_correction works as anticipated."""
        pvals = np.array([.1, .7, .5, .3, .9])
        exp = np.array([.5, .7 * 5 / 4., .5 * 5 / 3., .3 * 5 / 2., .9])
        obs = fdr_correction(pvals)
        np.testing.assert_allclose(obs, exp)

    def test_benjamini_hochberg_step_down(self):
        """Test that the BH step down procedure behaves as it does in R."""
        # r values
        # q = c(0.64771481,  0.93517796,  0.7169902 ,  0.18223457,  0.26918556,
        #  0.1450153 ,  0.22448242,  0.74723508,  0.89061034,  0.74007906)
        # p.adjust(q, method='BH')
        #  [1] 0.9340439 0.9351780 0.9340439 0.6729639 0.6729639 0.6729639
        #      0.6729639
        #  [8] 0.9340439 0.9351780 0.9340439
        pvals = np.array([0.64771481, 0.93517796, 0.7169902, 0.18223457,
                          0.26918556, 0.1450153, 0.22448242, 0.74723508,
                          0.89061034, 0.74007906])
        exp = np.array([0.9340439, 0.9351780, 0.9340439, 0.6729639, 0.6729639,
                        0.6729639, 0.6729639, 0.9340439, 0.9351780, 0.9340439])
        obs = benjamini_hochberg_step_down(pvals)
        np.testing.assert_allclose(obs, exp)
        # example 2
        pvals = np.array([1.32305426, 1.9345059, 0.87129877, 1.89957702,
                          1.85712616, 0.68757988, 0.41248969, 0.20751712,
                          1.97658599, 1.06209437])
        exp = np.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])
        obs = benjamini_hochberg_step_down(pvals)
        np.testing.assert_allclose(obs, exp)

    def test_bonferroni_correction(self):
        """Test that Bonferroni correction behaves correctly."""
        pvals = np.array([.1, .7, .5, .3, .9])
        exp = pvals * 5.
        obs = bonferroni_correction(pvals)
        np.testing.assert_allclose(obs, exp)

    def test_fisher_z_transform(self):
        '''Test Fisher Z transform is correct.'''
        r = .657
        exp = .5 * np.log(1.657 / .343)
        obs = fisher_z_transform(r)
        np.testing.assert_allclose(exp, obs)
        r = 1
        obs = fisher_z_transform(r)
        np.testing.assert_allclose(obs, np.nan)
        r = -1
        obs = fisher_z_transform(r)
        np.testing.assert_allclose(obs, np.nan)
        r = -5.6
        obs = fisher_z_transform(r)
        np.testing.assert_allclose(obs, np.nan)
        # from sokal and rohlf pg 575
        r = .972
        obs = fisher_z_transform(r)
        exp = 2.12730
        np.testing.assert_allclose(exp, obs, rtol=1e-4)

    def test_z_transform_pval(self):
        '''Test that pval associated with Fisher Z is correct.'''
        r = .6
        n = 100
        obs = z_transform_pval(r, n)
        exp = 3.4353390341723208e-09
        np.testing.assert_allclose(exp, obs)
        r = .5
        n = 3
        obs = z_transform_pval(r, n)
        np.testing.assert_allclose(obs, np.nan)

    def test_inverse_fisher_z_transform(self):
        '''Test that Fisher's Z transform is computed correctly.'''
        z = .65
        exp = 0.5716699660851171
        obs = inverse_fisher_z_transform(z)
        np.testing.assert_allclose(exp, obs)

    def test_fisher_population_correlation(self):
        '''Test that the population rho and homogeneity coeff are correct.'''
        # note: the error tolerances are lower than they would normally be
        # because sokal and rolhf don't give many significant figures
        # example from Sokal and Rohlf Biometry pg. 580 - 582
        rs = np.array([.29, .7, .58, .56, .55, .67, .65, .61, .64, .56])
        ns = np.array([100, 46, 28, 74, 33, 27, 52, 26, 20, 17])
        zbar = .615268
        X2 = 15.26352
        pop_r = .547825
        hval = chi_high(X2, len(ns) - 1)
        obs_p_rho, obs_hval = fisher_population_correlation(rs, ns)
        np.testing.assert_allclose(obs_p_rho, pop_r, rtol=1e-5)
        np.testing.assert_allclose(obs_hval, hval, rtol=1e-5)
        # test with np.nans
        rs = np.array(
            [.29, .7, np.nan, .58, .56, .55, .67, .65, .61, .64, .56])
        ns = np.array([100, 46, 400, 28, 74, 33, 27, 52, 26, 20, 17])
        obs_p_rho, obs_hval = fisher_population_correlation(rs, ns)
        np.testing.assert_allclose(obs_p_rho, pop_r, rtol=1e-5)
        np.testing.assert_allclose(obs_hval, hval, rtol=1e-5)
        # test with short vectors
        rs = [.6, .5, .4, .6, .7]
        ns = [10, 12, 42, 11, 3]
        obs_p_rho, obs_hval = fisher_population_correlation(rs, ns)
        np.testing.assert_allclose(obs_p_rho, np.nan)
        np.testing.assert_allclose(obs_hval, np.nan)
        # test with data with rs >1
        rs = [.6, .5, .4, 1.4]
        ns = [10, 50, 100, 10]
        self.assertRaises(ValueError, fisher_population_correlation, rs, ns)

    def test_assign_correlation_pval(self):
        '''Test that correlation pvalues are assigned correctly with each meth.
        '''
        # test with parametric t distribution, use example from Sokal and Rohlf
        # Biometry pg 576.
        r = .86519
        n = 12
        ts = 5.45618  # only 5 sig figs in sokal and rohlf
        exp = tprob(ts, n - 2)
        obs = assign_correlation_pval(r, n, 'parametric_t_distribution')
        np.testing.assert_allclose(exp, obs, rtol=1e-5)
        # test with too few samples
        n = 3
        self.assertRaises(ValueError, assign_correlation_pval, r, n,
                          'parametric_t_distribution')
        # test with fisher_z_transform
        r = .29
        n = 100
        z = 0.29856626366017841  # .2981 in biometry
        exp = z_transform_pval(z, n)
        obs = assign_correlation_pval(r, n, 'fisher_z_transform')
        np.testing.assert_allclose(exp, obs, rtol=1e-5)
        r = .61
        n = 26
        z = 0.70892135942740819  # .7089 in biometry
        exp = z_transform_pval(z, n)
        obs = assign_correlation_pval(r, n, 'fisher_z_transform')
        np.testing.assert_allclose(exp, obs, rtol=1e-5)
        # prove that we can have specify the other options, and as long as we
        # dont have bootstrapped selected we are fine.
        v1 = np.array([10, 11, 12])
        v2 = np.array([10, 14, 15])
        obs = assign_correlation_pval(r, n, 'fisher_z_transform',
                                      permutations=1000, perm_test_fn=pearson,
                                      v1=v1, v2=v2)
        np.testing.assert_allclose(exp, obs)
        # test with bootstrapping, seed for reproducibility.
        np.random.seed(0)
        v1 = np.array([54, 71, 60, 54, 42, 64, 43, 89, 96, 38])
        v2 = np.array([79, 52, 56, 92, 7, 8, 2, 83, 77, 87])
        # c = corrcoef(v1,v2)[0][1]
        exp = .357
        obs = assign_correlation_pval(0.33112494, 20000, 'bootstrapped',
                                      permutations=1000, perm_test_fn=pearson,
                                      v1=v1, v2=v2)
        np.testing.assert_allclose(exp, obs)
        # make sure it throws an error
        self.assertRaises(ValueError, assign_correlation_pval, 7, 20000,
                          'bootstrapped', perm_test_fn=pearson, v1=None, v2=v2)
        # test that it does properly with kendall
        exp = kendall_pval(r, n)
        obs = assign_correlation_pval(r, n, 'kendall')
        np.testing.assert_allclose(exp, obs)

    def test_cscore(self):
        '''Test cscore is calculated correctly.'''
        # test using example from Stone and Roberts pg 75
        v1 = np.array([1, 0, 0, 0, 1, 1, 0, 1, 0, 1])
        v2 = np.array([1, 1, 1, 0, 1, 0, 1, 1, 1, 0])
        obs = cscore(v1, v2)
        exp = 8
        self.assertEqual(obs, exp)
        # test using examples verified in ecosim
        v1 = np.array([4, 6, 12, 13, 14, 0, 0, 0, 14, 11, 9, 6, 0, 1, 1, 0, 0,
                       4])
        v2 = np.array([4, 0, 0, 113, 1, 2, 20, 0, 1, 0, 19, 16, 0, 13, 6, 0, 5,
                       4])
        # from R
        # library(vegan)
        # library(bipartite)
        # m = matrix(c(4,6,12,13,14,0,0,0,14,11,9,6,0,1,1,0,0,4,4,0,0,113,1,2,
        #              20,0,1,0,19,16,0,13,6,0,5,4), 18,2)
        # C.score(m, normalise=FALSE)
        exp = 9
        obs = cscore(v1, v2)
        self.assertEqual(obs, exp)

# execute tests if called from command line
if __name__ == '__main__':
    main()
