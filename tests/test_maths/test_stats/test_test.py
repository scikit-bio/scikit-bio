#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------
from __future__ import division
from bipy.util.unit_test import TestCase, main
from bipy.maths.stats.test import (G_2_by_2, G_fit, t_paired, t_one_sample,
                                   t_two_sample, mc_t_two_sample,
                                   _permute_observations, t_one_observation,
                                   correlation_test, ZeroExpectedError, fisher,
                                   safe_sum_p_log_p, permute_2d, mantel,
                                   mantel_test, _flatten_lower_triangle,
                                   pearson, spearman, _get_rank,
                                   ANOVA_one_way, mw_test, mw_boot,
                                   is_symmetric_and_hollow, reverse_tails, tail)

from numpy import (array, concatenate, fill_diagonal, reshape, arange, matrix,
                   ones, testing, tril, cov, sqrt)
import math


class TestsHelper(TestCase):

    """Class with utility methods useful for other tests."""

    def setUp(self):
        """Sets up variables used in the tests."""
        # How many times a p-value should be tested to fall in a given range
        # before failing the test.
        self.p_val_tests = 10

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

            self.assertIsProb(p_val)
            if p_val >= exp_min and p_val <= exp_max:
                found_match = True
                break
        self.assertTrue(found_match)


class TestsTests(TestCase):

    """Tests miscellaneous functions."""

    def test_tail(self):
        """tail should return x/2 if test is true; 1-(x/2) otherwise"""
        self.assertFloatEqual(tail(0.25, 'a' == 'a'), 0.25 / 2)
        self.assertFloatEqual(tail(0.25, 'a' != 'a'), 1 - (0.25 / 2))

    def test_fisher(self):
        """fisher results should match p 795 Sokal and Rohlf"""
        self.assertFloatEqual(fisher([0.073, 0.086, 0.10, 0.080, 0.060]),
                              0.0045957946540917905)

    def test_permute_2d(self):
        """permute_2d permutes rows and cols of a matrix."""
        a = reshape(arange(9), (3, 3))
        self.assertFloatEqual(permute_2d(a, [0, 1, 2]), a)
        self.assertFloatEqual(permute_2d(a, [2, 1, 0]),
                              array([[8, 7, 6], [5, 4, 3], [2, 1, 0]]))
        self.assertFloatEqual(permute_2d(a, [1, 2, 0]),
                              array([[4, 5, 3], [7, 8, 6], [1, 2, 0]]))


class GTests(TestCase):

    """Tests implementation of the G tests for fit and independence."""

    def test_G_2_by_2_2tailed_equal(self):
        """G_2_by_2 should return 0 if all cell counts are equal"""
        self.assertFloatEqual(0, G_2_by_2(1, 1, 1, 1, False, False)[0])
        self.assertFloatEqual(0, G_2_by_2(100, 100, 100, 100, False, False)[0])
        self.assertFloatEqual(0, G_2_by_2(100, 100, 100, 100, True, False)[0])

    def test_G_2_by_2_bad_data(self):
        """G_2_by_2 should raise ValueError if any counts are negative"""
        self.assertRaises(ValueError, G_2_by_2, 1, -1, 1, 1)

    def test_G_2_by_2_2tailed_examples(self):
        """G_2_by_2 values should match examples in Sokal & Rohlf"""
        # example from p 731, Sokal and Rohlf (1995)
        # without correction
        self.assertFloatEqual(G_2_by_2(12, 22, 16, 50, False, False)[0],
                              1.33249, 0.0001)
        self.assertFloatEqual(G_2_by_2(12, 22, 16, 50, False, False)[1],
                              0.24836, 0.0001)
        # with correction
        self.assertFloatEqual(G_2_by_2(12, 22, 16, 50, True, False)[0],
                              1.30277, 0.0001)
        self.assertFloatEqual(G_2_by_2(12, 22, 16, 50, True, False)[1],
                              0.25371, 0.0001)

    def test_G_2_by_2_1tailed_examples(self):
        """G_2_by_2 values should match values from codon_binding program"""
        # first up...the famous arginine case
        self.assertFloatEqualAbs(G_2_by_2(36, 16, 38, 106), (29.111609, 0),
                                 0.00001)
        # then some other miscellaneous positive and negative values
        self.assertFloatEqualAbs(
            G_2_by_2(0, 52, 12, 132), (-7.259930, 0.996474),
            0.00001)
        self.assertFloatEqualAbs(
            G_2_by_2(5, 47, 14, 130), (-0.000481, 0.508751),
            0.00001)
        self.assertFloatEqualAbs(
            G_2_by_2(5, 47, 36, 108), (-6.065167, 0.993106),
            0.00001)

    def test_Gfit_unequal_lists(self):
        """Gfit should raise errors if lists unequal"""
        # lists must be equal
        self.assertRaises(ValueError, G_fit, [1, 2, 3], [1, 2])

    def test_Gfit_negative_observeds(self):
        """Gfit should raise ValueError if any observeds are negative."""
        self.assertRaises(ValueError, G_fit, [-1, 2, 3], [1, 2, 3])

    def test_Gfit_nonpositive_expecteds(self):
        """Gfit should raise ZeroExpectedError if expecteds are zero/negative"""
        self.assertRaises(ZeroExpectedError, G_fit, [1, 2, 3], [0, 1, 2])
        self.assertRaises(ZeroExpectedError, G_fit, [1, 2, 3], [-1, 1, 2])

    def test_Gfit_good_data(self):
        """Gfit tests for fit should match examples in Sokal and Rohlf"""
        # example from p. 699, Sokal and Rohlf (1995)
        obs = [63, 31, 28, 12, 39, 16, 40, 12]
        exp = [67.78125, 22.59375, 22.59375, 7.53125, 45.18750,
               15.06250, 45.18750, 15.06250]
        # without correction
        self.assertFloatEqualAbs(G_fit(obs, exp, False)[0], 8.82397, 0.00002)
        self.assertFloatEqualAbs(G_fit(obs, exp, False)[1], 0.26554, 0.00002)
        # with correction
        self.assertFloatEqualAbs(G_fit(obs, exp)[0], 8.76938, 0.00002)
        self.assertFloatEqualAbs(G_fit(obs, exp)[1], 0.26964, 0.00002)

        # example from p. 700, Sokal and Rohlf (1995)
        obs = [130, 46]
        exp = [132, 44]
        # without correction
        self.assertFloatEqualAbs(G_fit(obs, exp, False)[0], 0.12002, 0.00002)
        self.assertFloatEqualAbs(G_fit(obs, exp, False)[1], 0.72901, 0.00002)
        # with correction
        self.assertFloatEqualAbs(G_fit(obs, exp)[0], 0.11968, 0.00002)
        self.assertFloatEqualAbs(G_fit(obs, exp)[1], 0.72938, 0.00002)

    def test_safe_sum_p_log_p(self):
        """safe_sum_p_log_p should ignore zero elements, not raise error"""
        m = array([2, 4, 0, 8])
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

    def test_t_paired_2tailed(self):
        """t_paired should match values from Sokal & Rohlf p 353"""
        x, y = self.x, self.y
        # check value of t and the probability for 2-tailed
        self.assertFloatEqual(t_paired(y, x)[0], 19.7203, 1e-4)
        self.assertFloatEqual(t_paired(y, x)[1], 1.301439e-11, 1e-4)

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
        self.assertFloatEqual(
            t_paired(y, x, "low")[1], 1 - (1.301439e-11 / 2), 1e-4)
        self.assertFloatEqual(
            t_paired(x, y, "high")[1], 1 - (1.301439e-11 / 2), 1e-4)
        self.assertFloatEqual(
            t_paired(y, x, "high")[1], 1.301439e-11 / 2, 1e-4)
        self.assertFloatEqual(
            t_paired(x, y, "low")[1], 1.301439e-11 / 2, 1e-4)

    def test_t_paired_specific_difference(self):
        """t_paired should allow a specific difference to be passed"""
        x, y = self.x, self.y
        # difference is 0.2, so test should be non-significant if 0.2 passed
        self.failIf(t_paired(y, x, exp_diff=0.2)[0] > 1e-10)
        # same, except that reversing list order reverses sign of difference
        self.failIf(t_paired(x, y, exp_diff=-0.2)[0] > 1e-10)
        # check that there's no significant difference from the true mean
        self.assertFloatEqual(
            t_paired(y, x, exp_diff=0.2)[1], 1, 1e-4)

    def test_t_paired_bad_data(self):
        """t_paired should raise ValueError on lists of different lengths"""
        self.assertRaises(ValueError, t_paired, self.y, [1, 2, 3])

    def test_t_two_sample(self):
        """t_two_sample should match example on p.225 of Sokal and Rohlf"""
        I = array([7.2, 7.1, 9.1, 7.2, 7.3, 7.2, 7.5])
        II = array([8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2])
        self.assertFloatEqual(t_two_sample(I, II), (-0.1184, 0.45385 * 2),
                              0.001)

    def test_t_two_sample_no_variance(self):
        """t_two_sample should properly handle lists that are invariant"""
        # By default should return (None, None) to mimic R's t.test.
        x = array([1, 1., 1])
        y = array([0, 0, 0.0])
        self.assertEqual(t_two_sample(x, x), (None, None))
        self.assertEqual(t_two_sample(x, y), (None, None))

        # Test none_on_zero_variance=False on various tail types. We use
        # self.assertEqual instead of self.assertFloatEqual because the latter
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
        x = array(range(-5, 5))
        y = array(range(-1, 10))
        self.assertFloatEqualAbs(t_one_sample(x), (-0.5222, 0.6141), 1e-4)
        self.assertFloatEqualAbs(t_one_sample(y), (4, 0.002518), 1e-4)
        # do some one-tailed tests as well
        self.assertFloatEqualAbs(
            t_one_sample(y, tails='low'), (4, 0.9987), 1e-4)
        self.assertFloatEqualAbs(
            t_one_sample(y, tails='high'), (4, 0.001259), 1e-4)

    def test_t_two_sample_switch(self):
        """t_two_sample should call t_one_observation if 1 item in sample."""
        sample = array([4.02, 3.88, 3.34, 3.87, 3.18])
        x = array([3.02])
        self.assertFloatEqual(t_two_sample(x, sample), (-1.5637254, 0.1929248))
        self.assertFloatEqual(t_two_sample(sample, x), (1.5637254, 0.1929248))

        # can't do the test if both samples have single item
        self.assertEqual(t_two_sample(x, x), (None, None))

        # Test special case if t=0.
        self.assertFloatEqual(t_two_sample([2], [1, 2, 3]), (0.0, 1.0))
        self.assertFloatEqual(t_two_sample([1, 2, 3], [2]), (0.0, 1.0))

    def test_t_one_observation(self):
        """t_one_observation should match p. 228 of Sokal and Rohlf"""
        sample = array([4.02, 3.88, 3.34, 3.87, 3.18])
        x = 3.02
        # note that this differs after the 3rd decimal place from what's in the
        # book, because Sokal and Rohlf round their intermediate steps...
        self.assertFloatEqual(t_one_observation(x, sample),
                              (-1.5637254, 0.1929248))

    def test_t_one_observation_no_variance(self):
        """t_one_observation should correctly handle an invariant list."""
        sample = array([1.0, 1.0, 1.0])

        # Can't perform test if invariant list's single value matches x,
        # regardless of none_on_zero_variance.
        self.assertEqual(t_one_observation(1, sample), (None, None))
        self.assertEqual(t_one_observation(1, sample,
                                           none_on_zero_variance=False), (None, None))

        # Test correct handling of none_on_zero_variance.
        self.assertEqual(t_one_observation(2, sample), (None, None))
        self.assertEqual(t_one_observation(2, sample,
                                           none_on_zero_variance=False), (float('inf'), 0.0))
        self.assertEqual(t_one_observation(2, sample,
                                           none_on_zero_variance=False, tails='low'), (float('inf'), 1.0))

    def test_mc_t_two_sample(self):
        """Test gives correct results with valid input data."""
        # Verified against R's t.test() and Deducer::perm.t.test().

        # With numpy array as input.
        exp = (-0.11858541225631833, 0.90756579317867436)
        I = array([7.2, 7.1, 9.1, 7.2, 7.3, 7.2, 7.5])
        II = array([8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2])
        obs = mc_t_two_sample(I, II)
        self.assertFloatEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 999)
        self.assertCorrectPValue(0.8, 0.9, mc_t_two_sample, [I, II],
                                 p_val_idx=3)

        # With python list as input.
        exp = (-0.11858541225631833, 0.90756579317867436)
        I = [7.2, 7.1, 9.1, 7.2, 7.3, 7.2, 7.5]
        II = [8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2]
        obs = mc_t_two_sample(I, II)
        self.assertFloatEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 999)
        self.assertCorrectPValue(0.8, 0.9, mc_t_two_sample, [I, II],
                                 p_val_idx=3)

        exp = (-0.11858541225631833, 0.45378289658933718)
        obs = mc_t_two_sample(I, II, tails='low')
        self.assertFloatEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 999)
        self.assertCorrectPValue(0.4, 0.47, mc_t_two_sample, [I, II],
                                 {'tails': 'low'}, p_val_idx=3)

        exp = (-0.11858541225631833, 0.54621710341066287)
        obs = mc_t_two_sample(I, II, tails='high', permutations=99)
        self.assertFloatEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 99)
        self.assertCorrectPValue(0.4, 0.62, mc_t_two_sample, [I, II],
                                 {'tails': 'high', 'permutations': 99}, p_val_idx=3)

        exp = (-2.8855783649036986, 0.99315596652421401)
        obs = mc_t_two_sample(I, II, tails='high', permutations=99, exp_diff=1)
        self.assertFloatEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 99)
        self.assertCorrectPValue(0.55, 0.99, mc_t_two_sample, [I, II],
                                 {'tails': 'high', 'permutations': 99, 'exp_diff': 1}, p_val_idx=3)

    def test_mc_t_two_sample_unbalanced_obs(self):
        """Test gives correct results with unequal number of obs per sample."""
        # Verified against R's t.test() and Deducer::perm.t.test().
        exp = (-0.10302479888889175, 0.91979753020527177)
        I = array([7.2, 7.1, 9.1, 7.2, 7.3, 7.2])
        II = array([8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2])
        obs = mc_t_two_sample(I, II)
        self.assertFloatEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 999)
        self.assertCorrectPValue(0.8, 0.9, mc_t_two_sample, [I, II],
                                 p_val_idx=3)

    def test_mc_t_two_sample_single_obs_sample(self):
        """Test works correctly with one sample having a single observation."""
        sample = array([4.02, 3.88, 3.34, 3.87, 3.18])
        x = array([3.02])
        exp = (-1.5637254, 0.1929248)
        obs = mc_t_two_sample(x, sample)
        self.assertFloatEqual(obs[:2], exp)
        self.assertFloatEqual(len(obs[2]), 999)
        self.assertIsProb(obs[3])

        exp = (1.5637254, 0.1929248)
        obs = mc_t_two_sample(sample, x)
        self.assertFloatEqual(obs[:2], exp)
        self.assertFloatEqual(len(obs[2]), 999)
        self.assertIsProb(obs[3])

        # Test the case where we can have no variance in the permuted lists.
        x = array([1, 1, 2])
        y = array([1])
        exp = (0.5, 0.666666666667)
        obs = mc_t_two_sample(x, y)
        self.assertFloatEqual(obs[:2], exp)
        self.assertFloatEqual(len(obs[2]), 999)
        self.assertIsProb(obs[3])

    def test_mc_t_two_sample_no_perms(self):
        """Test gives empty permutation results if no perms are given."""
        exp = (-0.11858541225631833, 0.90756579317867436, [], None)
        I = array([7.2, 7.1, 9.1, 7.2, 7.3, 7.2, 7.5])
        II = array([8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2])
        obs = mc_t_two_sample(I, II, permutations=0)
        self.assertFloatEqual(obs[0], exp[0])
        self.assertFloatEqual(obs[1], exp[1])
        self.assertEqual(obs[2], exp[2])
        self.assertEqual(obs[3], exp[3])

    def test_mc_t_two_sample_no_mc(self):
        """Test no MC stats if initial t-test is bad."""
        x = array([1, 1, 1])
        y = array([0, 0, 0])
        self.assertEqual(mc_t_two_sample(x, x), (None, None, [], None))

    def test_mc_t_two_sample_no_variance(self):
        """Test input with no variance. Should match Deducer::perm.t.test."""
        x = array([1, 1, 1])
        y = array([2, 2, 2])

        exp = (float('-inf'), 0.0)
        obs = mc_t_two_sample(x, y, permutations=10000)

        self.assertEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 10000)
        self.assertCorrectPValue(0.09, 0.11, mc_t_two_sample, [x, y],
                                 {'permutations': 10000}, p_val_idx=3)

        exp = (float('inf'), 0.0)
        obs = mc_t_two_sample(y, x, permutations=10000)

        self.assertEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 10000)
        self.assertCorrectPValue(0.09, 0.11, mc_t_two_sample, [y, x],
                                 {'permutations': 10000}, p_val_idx=3)

        exp = (float('-inf'), 1.0)
        obs = mc_t_two_sample(x, y, permutations=10000, tails='high')

        self.assertEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 10000)
        self.assertCorrectPValue(0.9999, 1.0, mc_t_two_sample, [x, y],
                                 {'permutations': 10000, 'tails': 'high'},
                                 p_val_idx=3)

        exp = (float('-inf'), 0.0)
        obs = mc_t_two_sample(x, y, permutations=10000, tails='low')

        self.assertEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 10000)
        self.assertCorrectPValue(0.04, 0.051, mc_t_two_sample, [x, y],
                                 {'permutations': 10000, 'tails': 'low'},
                                 p_val_idx=3)

    def test_mc_t_two_sample_no_permuted_variance(self):
        """Test with chance of getting no variance with some perms."""
        # Verified against R's t.test() and Deducer::perm.t.test().
        x = array([1, 1, 2])
        y = array([2, 2, 1])

        exp = (-0.70710678118654791, 0.51851851851851838)
        obs = mc_t_two_sample(x, y, permutations=10000)

        self.assertFloatEqual(obs[:2], exp)
        self.assertEqual(len(obs[2]), 10000)
        self.assertCorrectPValue(0.97, 1.0, mc_t_two_sample, [x, y],
                                 {'permutations': 10000}, p_val_idx=3)

    def test_mc_t_two_sample_invalid_input(self):
        """Test fails on various invalid input."""
        self.assertRaises(ValueError, mc_t_two_sample, [1, 2, 3], [4., 5., 4.],
                          tails='foo')
        self.assertRaises(ValueError, mc_t_two_sample, [1, 2, 3], [4., 5., 4.],
                          permutations=-1)
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
        self.assertFloatEqual(sorted(concatenate((obs[0][0], obs[1][0]))),
                              sorted(I + II))

    def test_reverse_tails(self):
        """reverse_tails should return 'high' if tails was 'low' or vice versa"""
        self.assertEqual(reverse_tails('high'), 'low')
        self.assertEqual(reverse_tails('low'), 'high')
        self.assertEqual(reverse_tails(None), None)
        self.assertEqual(reverse_tails(3), 3)

    def test_tail(self):
        """tail should return prob/2 if test is true, or 1-(prob/2) if false"""
        self.assertFloatEqual(tail(0.25, True), 0.125)
        self.assertFloatEqual(tail(0.25, False), 0.875)
        self.assertFloatEqual(tail(1, True), 0.5)
        self.assertFloatEqual(tail(1, False), 0.5)
        self.assertFloatEqual(tail(0, True), 0)
        self.assertFloatEqual(tail(0, False), 1)


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

    def test_mantel(self):
        """mantel should be significant for same matrix, not for random"""
        a = reshape(arange(25), (5, 5))
        a = tril(a) + tril(a).T
        fill_diagonal(a, 0)
        b = a.copy()
        # closely related -- should be significant
        self.assertCorrectPValue(0.0, 0.049, mantel, (a, b, 1000))

        c = reshape(ones(25), (5, 5))
        c[0, 1] = 3.0
        c[1, 0] = 3.0
        fill_diagonal(c, 0)
        # not related -- should not be significant
        self.assertCorrectPValue(0.06, 1.0, mantel, (a, c, 1000))

    def test_mantel_test_one_sided_greater(self):
        """Test one-sided mantel test (greater)."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test).
        m1 = array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        m2 = array([[0, 2, 7], [2, 0, 6], [7, 6, 0]])
        p, stat, perms = mantel_test(m1, m1, 999, alt='greater')
        self.assertFloatEqual(stat, 1.0)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.09, 0.25, mantel_test, (m1, m1, 999),
                                 {'alt': 'greater'})

        p, stat, perms = mantel_test(m1, m2, 999, alt='greater')
        self.assertFloatEqual(stat, 0.755928946018)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.2, 0.5, mantel_test, (m1, m2, 999),
                                 {'alt': 'greater'})

    def test_mantel_test_one_sided_less(self):
        """Test one-sided mantel test (less)."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test, but I modified their output to do a one-sided
        # less test).
        m1 = array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        m2 = array([[0, 2, 7], [2, 0, 6], [7, 6, 0]])
        m3 = array([[0, 0.5, 0.25], [0.5, 0, 0.1], [0.25, 0.1, 0]])
        p, stat, perms = mantel_test(m1, m1, 999, alt='less')
        self.assertFloatEqual(p, 1.0)
        self.assertFloatEqual(stat, 1.0)
        self.assertEqual(len(perms), 999)

        p, stat, perms = mantel_test(m1, m2, 999, alt='less')
        self.assertFloatEqual(stat, 0.755928946018)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.6, 1.0, mantel_test, (m1, m2, 999),
                                 {'alt': 'less'})

        p, stat, perms = mantel_test(m1, m3, 999, alt='less')
        self.assertFloatEqual(stat, -0.989743318611)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.1, 0.25, mantel_test, (m1, m3, 999),
                                 {'alt': 'less'})

    def test_mantel_test_two_sided(self):
        """Test two-sided mantel test."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test, but I modified their output to do a two-sided
        # test).
        m1 = array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        m2 = array([[0, 2, 7], [2, 0, 6], [7, 6, 0]])
        m3 = array([[0, 0.5, 0.25], [0.5, 0, 0.1], [0.25, 0.1, 0]])
        p, stat, perms = mantel_test(m1, m1, 999, alt='two sided')
        self.assertFloatEqual(stat, 1.0)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.20, 0.45, mantel_test, (m1, m1, 999),
                                 {'alt': 'two sided'})

        p, stat, perms = mantel_test(m1, m2, 999, alt='two sided')
        self.assertFloatEqual(stat, 0.755928946018)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.6, 0.75, mantel_test, (m1, m2, 999),
                                 {'alt': 'two sided'})

        p, stat, perms = mantel_test(m1, m3, 999, alt='two sided')
        self.assertFloatEqual(stat, -0.989743318611)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.2, 0.45, mantel_test, (m1, m3, 999),
                                 {'alt': 'two sided'})

    def test_mantel_test_invalid_distance_matrix(self):
        """Test mantel test with invalid distance matrix."""
        # Single asymmetric, non-hollow distance matrix.
        self.assertRaises(ValueError, mantel_test, array([[1, 2], [3, 4]]),
                          array([[0, 0], [0, 0]]), 999)

        # Two asymmetric distance matrices.
        self.assertRaises(ValueError, mantel_test, array([[0, 2], [3, 0]]),
                          array([[0, 1], [0, 0]]), 999)

    def test_mantel_test_invalid_input(self):
        """Test mantel test with invalid input."""
        self.assertRaises(ValueError, mantel_test, array([[1]]), array([[1]]),
                          999, alt='foo')
        self.assertRaises(ValueError, mantel_test, array([[1]]),
                          array([[1, 2], [3, 4]]), 999)
        self.assertRaises(ValueError, mantel_test, array([[1]]),
                          array([[1]]), 0)
        self.assertRaises(ValueError, mantel_test, array([[1]]),
                          array([[1]]), -1)

    def test_is_symmetric_and_hollow(self):
        """Should correctly test for symmetry and hollowness of dist mats."""
        self.assertTrue(is_symmetric_and_hollow(array([[0, 1], [1, 0]])))
        self.assertTrue(is_symmetric_and_hollow(matrix([[0, 1], [1, 0]])))
        self.assertTrue(is_symmetric_and_hollow(matrix([[0.0, 0], [0.0, 0]])))
        self.assertTrue(not is_symmetric_and_hollow(
            array([[0.001, 1], [1, 0]])))
        self.assertTrue(not is_symmetric_and_hollow(
            array([[0, 1.1], [1, 0]])))
        self.assertTrue(not is_symmetric_and_hollow(
            array([[0.5, 1.1], [1, 0]])))

    def test_flatten_lower_triangle(self):
        """Test flattening various dms' lower triangulars."""
        self.assertEqual(_flatten_lower_triangle(array([[8]])), [])
        self.assertEqual(_flatten_lower_triangle(array([[1, 2], [3, 4]])), [3])
        self.assertEqual(_flatten_lower_triangle(array([[1, 2, 3], [4, 5, 6],
                                                        [7, 8, 9]])), [4, 7, 8])

    def test_pearson(self):
        """Test pearson correlation method on valid data."""
        # This test output was verified by R.
        self.assertFloatEqual(pearson([1, 2], [1, 2]), 1.0)
        self.assertFloatEqual(pearson([1, 2, 3], [1, 2, 3]), 1.0)
        self.assertFloatEqual(pearson([1, 2, 3], [1, 2, 4]), 0.9819805)

    def test_pearson_invalid_input(self):
        """Test running pearson on bad input."""
        self.assertRaises(ValueError, pearson, [1.4, 2.5], [5.6, 8.8, 9.0])
        self.assertRaises(ValueError, pearson, [1.4], [5.6])

    def test_spearman(self):
        """Test the spearman function with valid input."""
        # One vector has no ties.
        exp = 0.3719581
        obs = spearman(self.a, self.b)
        self.assertFloatEqual(obs, exp)

        # Both vectors have no ties.
        exp = 0.2969697
        obs = spearman(self.b, self.c)
        self.assertFloatEqual(obs, exp)

        # Both vectors have ties.
        exp = 0.388381
        obs = spearman(self.a, self.r)
        self.assertFloatEqual(obs, exp)

        exp = -0.17575757575757578
        obs = spearman(self.data1, self.data2)
        self.assertFloatEqual(obs, exp)

    def test_spearman_no_variation(self):
        """Test the spearman function with a vector having no variation."""
        exp = 0.0
        obs = spearman([1, 1, 1], [1, 2, 3])
        self.assertFloatEqual(obs, exp)

    def test_spearman_ranked(self):
        """Test the spearman function with a vector that is already ranked."""
        exp = 0.2969697
        obs = spearman(self.b_ranked, self.c_ranked)
        self.assertFloatEqual(obs, exp)

    def test_spearman_one_obs(self):
        """Test running spearman on a single observation."""
        self.assertRaises(ValueError, spearman, [1.0], [5.0])

    def test_spearman_invalid_input(self):
        """Test the spearman function with invalid input."""
        self.assertRaises(ValueError, spearman, [], [])
        self.assertRaises(ValueError, spearman, self.a, [])
        self.assertRaises(TypeError, spearman, {0: 2}, [1, 2, 3])

    def test_get_rank(self):
        """Test the _get_rank function with valid input."""
        exp = (
            [1.5,
             3.5,
             7.5,
             5.5,
             1.5,
             9.0,
             10.0,
             11.0,
             12.0,
             7.5,
             14.0,
             3.5,
             5.5,
             13.0],
            4)
        obs = _get_rank(self.x)
        self.assertFloatEqual(exp[0], obs[0])
        self.assertFloatEqual(exp[1], obs[1])

        exp = ([1.5, 3.0, 5.5, 4.0, 1.5, 7.0, 8.0, 9.0, 10.0, 5.5], 2)
        obs = _get_rank(self.a)
        self.assertFloatEqual(exp[0], obs[0])
        self.assertFloatEqual(exp[1], obs[1])

        exp = ([2, 7, 10, 1, 3, 6, 4, 8, 5, 9], 0)
        obs = _get_rank(self.b)
        self.assertFloatEqual(exp[0], obs[0])
        self.assertFloatEqual(exp[1], obs[1])

        exp = ([1.5, 7.0, 10.0, 1.5, 3.0, 6.0, 4.0, 8.0, 5.0, 9.0], 1)
        obs = _get_rank(self.r)
        self.assertFloatEqual(exp[0], obs[0])
        self.assertFloatEqual(exp[1], obs[1])

        exp = ([], 0)
        obs = _get_rank([])
        self.assertEqual(exp, obs)

    def test_get_rank_invalid_input(self):
        """Test the _get_rank function with invalid input."""
        vec = [1, 'a', 3, 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

        vec = [1, 2, {1: 2}, 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

        vec = [1, 2, [23, 1], 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

        vec = [1, 2, (1,), 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

    def test_correlation_test_pearson(self):
        """Test correlation_test using pearson on valid input."""
        # These results were verified with R.

        # Test with non-default confidence level and permutations.
        obs = correlation_test(self.data1, self.data2, method='pearson',
                               confidence_level=0.90, permutations=990)
        self.assertFloatEqual(obs[:2], (-0.03760147, 0.91786297277172868))
        self.assertEqual(len(obs[2]), 990)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.9, 0.93, correlation_test,
                                 (self.data1, self.data2),
                                 {'method': 'pearson', 'confidence_level': 0.90,
                                  'permutations': 990}, p_val_idx=3)
        self.assertFloatEqual(obs[4], (-0.5779077, 0.5256224))

        # Test with non-default tail type.
        obs = correlation_test(self.data1, self.data2, method='pearson',
                               confidence_level=0.90, permutations=990,
                               tails='low')
        self.assertFloatEqual(obs[:2], (-0.03760147, 0.45893148638586434))
        self.assertEqual(len(obs[2]), 990)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.41, 0.46, correlation_test,
                                 (self.data1, self.data2),
                                 {'method': 'pearson', 'confidence_level': 0.90,
                                  'permutations': 990, 'tails': 'low'}, p_val_idx=3)
        self.assertFloatEqual(obs[4], (-0.5779077, 0.5256224))

    def test_correlation_test_spearman(self):
        """Test correlation_test using spearman on valid input."""
        # This example taken from Wikipedia page:
        # http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient
        obs = correlation_test(self.data1, self.data2, method='spearman',
                               tails='high')
        self.assertFloatEqual(obs[:2], (-0.17575757575757578, 0.686405827612))
        self.assertEqual(len(obs[2]), 999)
        for rho in obs[2]:
            self.assertTrue(rho >= -1.0 and rho <= 1.0)
        self.assertCorrectPValue(0.67, 0.7, correlation_test,
                                 (self.data1, self.data2),
                                 {'method': 'spearman', 'tails': 'high'}, p_val_idx=3)
        self.assertFloatEqual(obs[4],
                              (-0.7251388558041697, 0.51034422964834503))

        # The p-value is off because the example uses a one-tailed test, while
        # we use a two-tailed test. Someone confirms the answer that we get
        # here for a two-tailed test:
        # http://stats.stackexchange.com/questions/22816/calculating-p-value-
        #     for-spearmans-rank-correlation-coefficient-example-on-wikip
        obs = correlation_test(self.data1, self.data2, method='spearman',
                               tails=None)
        self.assertFloatEqual(obs[:2],
                              (-0.17575757575757578, 0.62718834477648433))
        self.assertEqual(len(obs[2]), 999)
        for rho in obs[2]:
            self.assertTrue(rho >= -1.0 and rho <= 1.0)
        self.assertCorrectPValue(0.60, 0.64, correlation_test,
                                 (self.data1, self.data2),
                                 {'method': 'spearman', 'tails': None}, p_val_idx=3)
        self.assertFloatEqual(obs[4],
                              (-0.7251388558041697, 0.51034422964834503))

    def test_correlation_test_invalid_input(self):
        """Test correlation_test using invalid input."""
        self.assertRaises(ValueError, correlation_test, self.data1, self.data2,
                          method='foo')
        self.assertRaises(ValueError, correlation_test, self.data1, self.data2,
                          tails='foo')
        self.assertRaises(ValueError, correlation_test, self.data1, self.data2,
                          permutations=-1)
        self.assertRaises(ValueError, correlation_test, self.data1, self.data2,
                          confidence_level=-1)
        self.assertRaises(ValueError, correlation_test, self.data1, self.data2,
                          confidence_level=1.1)
        self.assertRaises(ValueError, correlation_test, self.data1, self.data2,
                          confidence_level=0)
        self.assertRaises(ValueError, correlation_test, self.data1, self.data2,
                          confidence_level=0.0)
        self.assertRaises(ValueError, correlation_test, self.data1, self.data2,
                          confidence_level=1)
        self.assertRaises(ValueError, correlation_test, self.data1, self.data2,
                          confidence_level=1.0)

    def test_correlation_test_no_permutations(self):
        """Test correlation_test with no permutations."""
        # These results were verified with R.
        exp = (-0.2581988897471611, 0.7418011102528389, [], None,
               (-0.97687328610475876, 0.93488023560400879))
        obs = correlation_test([1, 2, 3, 4], [1, 2, 1, 1], permutations=0)
        self.assertFloatEqual(obs[0], exp[0])
        self.assertFloatEqual(obs[1], exp[1])
        self.assertFloatEqual(obs[2], exp[2])
        self.assertEqual(obs[3], exp[3])
        self.assertFloatEqual(obs[4], exp[4])

    def test_correlation_test_perfect_correlation(self):
        """Test correlation_test with perfectly-correlated input vectors."""
        # These results were verified with R.
        obs = correlation_test([1, 2, 3, 4], [1, 2, 3, 4])
        self.assertFloatEqual(obs[:2],
                              (0.99999999999999978, 2.2204460492503131e-16))
        self.assertEqual(len(obs[2]), 999)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.06, 0.09, correlation_test,
                                 ([1, 2, 3, 4], [1, 2, 3, 4]), p_val_idx=3)
        self.assertFloatEqual(obs[4], (0.99999999999998879, 1.0))

    def test_correlation_test_small_obs(self):
        """Test correlation_test with a small number of observations."""
        # These results were verified with R.
        obs = correlation_test([1, 2, 3], [1, 2, 3])
        self.assertFloatEqual(obs[:2], (1.0, 0))
        self.assertEqual(len(obs[2]), 999)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.3, 0.4, correlation_test,
                                 ([1, 2, 3], [1, 2, 3]), p_val_idx=3)
        self.assertEqual(obs[4], (None, None))

        obs = correlation_test([1, 2, 3], [1, 2, 3], method='spearman')
        self.assertFloatEqual(obs[:2], (1.0, 0))
        self.assertEqual(len(obs[2]), 999)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.3, 0.4, correlation_test,
                                 ([1, 2, 3], [1, 2, 3]), {'method': 'spearman'}, p_val_idx=3)
        self.assertEqual(obs[4], (None, None))


class MannWhitneyTests(TestCase):

    """check accuracy of Mann-Whitney implementation"""
    x = map(int, "104 109 112 114 116 118 118 119 121 123 125 126"
            " 126 128 128 128".split())
    y = map(int, "100 105 107 107 108 111 116 120 121 123".split())

    def test_mw_test(self):
        """mann-whitney test results should match Sokal & Rohlf"""
        U, p = mw_test(self.x, self.y)
        self.assertFloatEqual(U, 123.5)
        self.assertTrue(0.02 <= p <= 0.05)

    def test_mw_boot(self):
        """excercising the Monte-carlo variant of mann-whitney"""
        U, p = mw_boot(self.x, self.y, 10)
        self.assertFloatEqual(U, 123.5)
        self.assertTrue(0 <= p <= 0.5)


class TestDistMatrixPermutationTest(TestCase):

    """Tests of distance_matrix_permutation_test"""

    def setUp(self):
        """sets up variables for testing"""
        self.matrix = array(
            [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
        self.cells = [(0, 1), (1, 3)]
        self.cells2 = [(0, 2), (2, 3)]

    def test_ANOVA_one_way(self):
        """ANOVA one way returns same values as ANOVA on a stats package
        """
        g1 = array([10.0, 11.0, 10.0, 5.0, 6.0])
        g2 = array([1.0, 2.0, 3.0, 4.0, 1.0, 2.0])
        g3 = array([6.0, 7.0, 5.0, 6.0, 7.0])
        i = [g1, g2, g3]
        F, pval = ANOVA_one_way(i)

        self.assertFloatEqual(F, 18.565450643776831)
        self.assertFloatEqual(pval, 0.00015486238993089464)

# execute tests if called from command line
if __name__ == '__main__':
    main()
