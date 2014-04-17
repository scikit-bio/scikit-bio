#!/usr/bin/env python
from __future__ import division

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from itertools import izip
from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio.maths.diversity.alpha import (ace, berger_parker_d, brillouin_d,
                                         chao1, chao1_confidence, _chao1_var,
                                         dominance, doubles, enspie,
                                         equitability, esty_ci, fisher_alpha,
                                         gini_index, goods_coverage, heip_e,
                                         kempton_taylor_q, margalef,
                                         mcintosh_d, mcintosh_e, menhinick,
                                         observed_species, osd, robbins,
                                         shannon, simpson, simpson_e,
                                         simpson_reciprocal, singles, strong,
                                         _indices_to_counts, _lorenz_curve,
                                         _lorenz_curve_integrator,
                                         _expand_counts,
                                         lladser_point_estimates,
                                         get_interval_for_r_new_species,
                                         lladser_ci_series, lladser_ci_from_r)


class AlphaDiversityTests(TestCase):
    def setUp(self):
        self.counts = np.array([0, 1, 1, 4, 2, 5, 2, 4, 1, 2])
        self.no_singles = np.array([0, 2, 2, 4, 5, 0, 0, 0, 0, 0])
        self.no_doubles = np.array([0, 1, 1, 4, 5, 0, 0, 0, 0, 0])

        # For Gini index and related tests.
        self.gini_data = np.array([4.5, 6.7, 3.4, 15., 18., 3.5, 6.7, 14.1])
        self.gini_lorenz_curve_points = [(0.125, 0.047287899860917935),
                                         (0.25, 0.095966620305980521),
                                         (0.375, 0.15855354659248957),
                                         (0.5, 0.2517385257301808),
                                         (0.625, 0.34492350486787204),
                                         (0.75, 0.541029207232267),
                                         (0.875, 0.74965229485396379),
                                         (1.0, 1.0)]

    def diversity(self, indices, f, step=1, start=None):
        """Calculate diversity index for each window of size step.

        indices: vector of indices of species
        f: f(counts) -> diversity measure
        start: first index to sum up to (default: step)
        step: step size (default:1)

        """
        result = []
        if start is None:
            start = step
        freqs = np.zeros(max(indices) + 1)
        i = 0
        for j in range(start, len(indices) + 1, step):
            freqs = _indices_to_counts(indices[i:j], freqs)
            try:
                curr = f(freqs)
            except (ZeroDivisionError, FloatingPointError):
                curr = 0
            result.append(curr)
            i = j
        return np.array(result)

    def test_ace(self):
        self.assertAlmostEqual(ace(np.array([2, 0])), 1.0)
        self.assertAlmostEqual(ace(np.array([12, 0, 9])), 2.0)
        self.assertAlmostEqual(ace(np.array([12, 2, 8])), 3.0)
        self.assertAlmostEqual(ace(np.array([12, 2, 1])), 4.0)
        self.assertAlmostEqual(ace(np.array([12, 1, 2, 1])), 7.0)
        self.assertAlmostEqual(ace(np.array([12, 3, 2, 1])), 4.6)
        self.assertAlmostEqual(ace(np.array([12, 3, 6, 1, 10])), 5.62749672)

        # Just returns the number of species when all are abundant.
        self.assertAlmostEqual(ace(np.array([12, 12, 13, 14])), 4.0)

    def test_berger_parker_d(self):
        self.assertEqual(berger_parker_d(np.array([5])), 1)
        self.assertEqual(berger_parker_d(np.array([5, 5])), 0.5)
        self.assertEqual(berger_parker_d(np.array([1, 1, 1, 1, 0])), 0.25)
        self.assertEqual(berger_parker_d(self.counts), 5 / 22)

    def test_brillouin_d(self):
        self.assertAlmostEqual(brillouin_d(np.array([1, 2, 0, 0, 3, 1])),
                               0.86289353018248782)

    def test_chao1(self):
        self.assertEqual(chao1(self.counts), 9.75)
        self.assertEqual(chao1(self.counts, bias_corrected=False), 10.5)

        self.assertEqual(chao1(self.no_singles), 4)
        self.assertEqual(chao1(self.no_singles, bias_corrected=False), 4)

        self.assertEqual(chao1(self.no_doubles), 5)
        self.assertEqual(chao1(self.no_doubles, bias_corrected=False), 5)

    def test_chao1_confidence(self):
        # Should match observed results from EstimateS. NOTE: EstimateS rounds
        # to 2 dp.
        obs = chao1_confidence(self.counts)
        npt.assert_allclose(obs, (9.07, 17.45), rtol=0.01)

        obs = chao1_confidence(self.counts, bias_corrected=False)
        npt.assert_allclose(obs, (9.17, 21.89), rtol=0.01)

        obs = chao1_confidence(self.no_singles)
        npt.assert_array_almost_equal(obs, (4, 4.95), decimal=2)

        obs = chao1_confidence(self.no_singles, bias_corrected=False)
        npt.assert_array_almost_equal(obs, (4, 4.95), decimal=2)

        obs = chao1_confidence(self.no_doubles)
        npt.assert_array_almost_equal(obs, (4.08, 17.27), decimal=2)

        obs = chao1_confidence(self.no_doubles, bias_corrected=False)
        npt.assert_array_almost_equal(obs, (4.08, 17.27), decimal=2)

    def test_chao1_var(self):
        # Should match observed results from EstimateS.NOTE: EstimateS reports
        # sd, not var, and rounds to 2 dp.
        obs = _chao1_var(self.counts)
        npt.assert_allclose(obs, 1.42 ** 2, rtol=0.01)

        obs = _chao1_var(self.counts, bias_corrected=False)
        npt.assert_allclose(obs, 2.29 ** 2, rtol=0.01)

        obs = _chao1_var(self.no_singles)
        self.assertAlmostEqual(obs, 0.39 ** 2, delta=0.01)

        obs = _chao1_var(self.no_singles, bias_corrected=False)
        self.assertAlmostEqual(obs, 0.39 ** 2, delta=0.01)

        obs = _chao1_var(self.no_doubles)
        self.assertAlmostEqual(obs, 2.17 ** 2, delta=0.01)

        obs = _chao1_var(self.no_doubles, bias_corrected=False)
        self.assertAlmostEqual(obs, 2.17 ** 2, delta=0.01)

    def test_dominance(self):
        self.assertEqual(dominance(np.array([5])), 1)
        self.assertAlmostEqual(dominance(np.array([1, 0, 2, 5, 2])), 0.34)

    def test_doubles(self):
        self.assertEqual(doubles(self.counts), 3)
        self.assertEqual(doubles(np.array([0, 3, 4])), 0)
        self.assertEqual(doubles(np.array([2])), 1)

    def test_enspie(self):
        # Totally even community should have ENS_pie = number of species.
        self.assertAlmostEqual(enspie(np.array([1, 1, 1, 1, 1, 1])), 6)
        self.assertAlmostEqual(enspie(np.array([13, 13, 13, 13])), 4)

        # Hand calculated.
        arr = np.array([1, 41, 0, 0, 12, 13])
        exp = 1 / ((arr / arr.sum()) ** 2).sum()
        self.assertAlmostEqual(enspie(arr), exp)

        # Using dominance.
        exp = 1 / dominance(arr)
        self.assertAlmostEqual(enspie(arr), exp)

    def test_simpson_reciprocal(self):
        arr = np.array([1, 0, 2, 5, 2])
        self.assertAlmostEqual(simpson_reciprocal(arr), 1 / dominance(arr))

    def test_equitability(self):
        self.assertAlmostEqual(equitability(np.array([5, 5])), 1)
        self.assertAlmostEqual(equitability(np.array([1, 1, 1, 1, 0])), 1)

    def test_esty_ci(self):
        data = [1, 1, 2, 1, 1, 3, 2, 1, 3, 4]

        (observed_upper, observed_lower) = zip(
            *self.diversity(data, esty_ci, step=1))

        expected_upper = np.array([1, 1.38590382, 1.40020259, 0.67434465,
                                   0.55060902, 0.71052858, 0.61613483,
                                   0.54041008, 0.43554755, 0.53385652])
        expected_lower = np.array([1, -1.38590382, -0.73353593, -0.17434465,
                                   -0.15060902, -0.04386191, -0.33042054,
                                   -0.29041008, -0.43554755, -0.33385652])

        npt.assert_array_almost_equal(observed_upper, expected_upper)
        npt.assert_array_almost_equal(observed_lower, expected_lower)

    def test_fisher_alpha(self):
        arr = np.array([4, 3, 4, 0, 1, 0, 2])
        obs = fisher_alpha(arr)
        self.assertAlmostEqual(obs, 2.7823795367398798)

    def test_gini_index(self):
        exp = 0.32771210013908214
        obs = gini_index(self.gini_data, 'trapezoids')
        self.assertAlmostEqual(obs, exp)

        exp = 0.20271210013908214
        obs = gini_index(self.gini_data, 'rectangles')
        self.assertAlmostEqual(obs, exp)

    def test_goods_coverage(self):
        counts = [1] * 75 + [2, 2, 2, 2, 2, 2, 3, 4, 4]
        obs = goods_coverage(counts)
        self.assertAlmostEqual(obs, 0.23469387755)

    def test_heip_e(self):
        arr = np.array([1, 2, 3, 1])
        h = shannon(arr, base=np.e)
        expected = np.exp(h - 1) / 3
        self.assertEqual(heip_e(arr), expected)

    def test_kempton_taylor_q(self):
        # Approximate Magurran 1998 calculation p143.
        arr = np.array([2, 3, 3, 3, 3, 3, 4, 4, 4, 6, 6, 7, 7, 9, 9, 11, 14,
                        15, 15, 20, 29, 33, 34, 36, 37, 53, 57, 138, 146, 170])
        exp = 14 / np.log(34 / 4)
        self.assertAlmostEqual(kempton_taylor_q(arr), exp)

        # Should get same answer regardless of input order.
        np.random.shuffle(arr)
        self.assertAlmostEqual(kempton_taylor_q(arr), exp)

    def test_margalef(self):
        self.assertEqual(margalef(self.counts), 8 / np.log(22))

    def test_mcintosh_d(self):
        self.assertAlmostEqual(mcintosh_d(np.array([1, 2, 3])),
                               0.636061424871458)

    def test_mcintosh_e(self):
        num = np.sqrt(15)
        den = np.sqrt(19)
        exp = num / den
        self.assertEqual(mcintosh_e(np.array([1, 2, 3, 1])), exp)

    def test_menhinick(self):
        self.assertEqual(menhinick(self.counts), 9 / np.sqrt(22))

    def test_observed_species(self):
        obs = observed_species(np.array([4, 3, 4, 0, 1, 0, 2]))
        self.assertEqual(obs, 5)

        obs = observed_species(np.array([0, 0, 0]))
        self.assertEqual(obs, 0)

        obs = observed_species(self.counts)
        self.assertEqual(obs, 9)

    def test_osd(self):
        self.assertEqual(osd(self.counts), (9, 3, 3))

    def test_robbins(self):
        self.assertEqual(robbins(np.array([1, 2, 3, 0, 1])), 2 / 7)

    def test_shannon(self):
        self.assertEqual(shannon(np.array([5])), 0)
        self.assertEqual(shannon(np.array([5, 5])), 1)
        self.assertEqual(shannon(np.array([1, 1, 1, 1, 0])), 2)

    def test_simpson(self):
        self.assertAlmostEqual(simpson(np.array([1, 0, 2, 5, 2])), 0.66)
        self.assertAlmostEqual(simpson(np.array([5])), 0)

    def test_simpson_e(self):
        # A totally even community should have simpson_e = 1.
        self.assertEqual(simpson_e(np.array([1, 1, 1, 1, 1, 1, 1])), 1)

        arr = np.array([0, 30, 25, 40, 0, 0, 5])
        freq_arr = arr / arr.sum()
        D = (freq_arr ** 2).sum()
        exp = 1 / (D * 4)
        obs = simpson_e(arr)
        self.assertEqual(obs, exp)

        # From:
        # https://groups.nceas.ucsb.edu/sun/meetings/calculating-evenness-
        #   of-habitat-distributions
        arr = np.array([500, 400, 600, 500])
        D = 0.0625 + 0.04 + 0.09 + 0.0625
        exp = 1 / (D * 4)
        self.assertEqual(simpson_e(arr), exp)

    def test_singles(self):
        self.assertEqual(singles(self.counts), 3)
        self.assertEqual(singles(np.array([0, 3, 4])), 0)
        self.assertEqual(singles(np.array([1])), 1)

    def test_strong(self):
        self.assertAlmostEqual(strong(np.array([1, 2, 3, 1])), 0.214285714)

    def test_indices_to_counts(self):
        exp = np.array([1, 2, 0, 0, 0, 3])
        obs = _indices_to_counts(np.array([5, 0, 1, 1, 5, 5]))
        npt.assert_array_equal(obs, exp)

        exp = np.array([2, 3, 2, 0, 0, 3])
        obs = _indices_to_counts(np.array([2, 2, 1, 0]), obs)
        npt.assert_array_equal(obs, exp)

    def test_lorenz_curve(self):
        self.assertEqual(_lorenz_curve(self.gini_data),
                         self.gini_lorenz_curve_points)

        # Raises error on negative data.
        with self.assertRaises(ValueError):
            _lorenz_curve([1.0, -3.1, 4.5])

    def test_lorenz_curve_integrator(self):
        exp = 0.33614394993045893
        obs = _lorenz_curve_integrator(self.gini_lorenz_curve_points,
                                       'trapezoids')
        self.assertAlmostEqual(obs, exp)

        exp = 0.39864394993045893
        obs = _lorenz_curve_integrator(self.gini_lorenz_curve_points,
                                       'rectangles')
        self.assertAlmostEqual(obs, exp)

        # Raises error on invalid method.
        with self.assertRaises(ValueError):
            _lorenz_curve_integrator(self.gini_lorenz_curve_points, 'brofist')

    def test_expand_counts(self):
        arr = np.array([2, 0, 1, 2])
        npt.assert_array_equal(_expand_counts(arr), np.array([0, 0, 2, 3, 3]))

    def test_lladser_point_estimates(self):
        s = [5, 1, 5, 1, 2, 3, 1, 5, 3, 2, 5, 3]
        r = 3
        observed = list(lladser_point_estimates(s, r))
        self.assertEqual(len(observed), 3)

        for k in (range(3)):
            x = observed[k]
            t = x[2]
            self.assertEqual(x[0], (r - 1) / t)

        # Estimator has variance of (1-p)^2/(r-2),
        # which for r=7 and p=0.5 is 0.05
        seq = "WBWBWBWBWBWBWBWBWBWBWBWBWBWBWBWBWBW"
        reps = 1000
        sum = 0
        for i in range(reps):
            p, _, _ = list(lladser_point_estimates(seq, r=7))[0]
            sum += p
        self.assertTrue(0.45 < sum / reps and sum / reps < 0.55)

    def test_get_interval_for_r_new_species(self):
        s = [5, 1, 5, 1, 2, 3, 1, 5, 3, 2, 5]
        expected = [(3, set([5]), 4, 0),
                    (4, set([5, 1]), 6, 1),
                    (4, set([5, 1, 2]), 9, 4)]
        for x, y in izip(get_interval_for_r_new_species(s, 2), expected):
            self.assertEqual(y, x)

        s = [5, 5, 5, 5, 5]
        # never saw new one
        self.assertEqual(list(get_interval_for_r_new_species(s, 2)), [])

    def test_lladser_ci_series_exact(self):
        # Values are from Manuel's email of 9/11/09
        # have seen RWB
        urn_1 = 'RWBWWBWRRWRYWRPPZ'
        results = list(lladser_ci_series(urn_1, r=4))
        self.assertEqual(len(results), 3)

    def test_lladser_ci_series_random(self):
        seq = "WBWBWBWBWBWB"
        observations = []
        alpha = 0.95
        reps = 1000
        for i in range(reps):
            obs = list(lladser_ci_series(seq, r=4, alpha=alpha))[0]
            observations.append(obs)
        tps = filter(lambda a_b: a_b[0] < 0.5 and 0.5 < a_b[1], observations)
        self.assertTrue(len(tps) >= alpha * reps)  # 100%-95%

    def test_lladser_ci_from_r(self):
        f = 10
        t = 10
        r = 4
        obs_low, obs_high = lladser_ci_from_r(r=r, t=t, f=f)
        self.assertAlmostEqual(obs_low, 0.0806026244)
        self.assertAlmostEqual(obs_high, 0.806026244)

        r = 20
        t = 100
        obs_low, obs_high = lladser_ci_from_r(r=r, t=t, f=f)
        self.assertAlmostEqual(obs_low, 0.02787923964)
        self.assertAlmostEqual(obs_high, 0.2787923964)

        # make sure we test with each possible alpha
        alpha = 0.99
        obs_low, obs_high = lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha)
        self.assertAlmostEqual(obs_low, 0.03184536992)
        self.assertAlmostEqual(obs_high, 0.3184536992)

        alpha = 0.9
        r = 3
        obs_low, obs_high = lladser_ci_from_r(r=r, t=t, f=f, alpha=alpha)
        self.assertAlmostEqual(obs_low, 0.005635941995)
        self.assertAlmostEqual(obs_high, 0.05635941995)

        # test other ci_types
        ci_type = 'ULCU'
        obs_low, obs_high = lladser_ci_from_r(
            r=r, t=t, f=f, alpha=alpha, ci_type=ci_type)
        self.assertAlmostEqual(obs_low, 0.01095834700)
        self.assertAlmostEqual(obs_high, 0.1095834700)

        alpha = 0.95
        t = 10
        ci_type = 'U'
        obs_low, obs_high = lladser_ci_from_r(
            r=r, t=t, f=f, alpha=alpha, ci_type=ci_type)
        self.assertAlmostEqual(obs_low, 0)
        self.assertAlmostEqual(obs_high, 0.6295793622)

        ci_type = 'L'
        obs_low, obs_high = lladser_ci_from_r(
            r=r, t=t, f=f, alpha=alpha, ci_type=ci_type)
        self.assertAlmostEqual(obs_low, 0.0817691447)
        self.assertAlmostEqual(obs_high, 1)

        # Requesting CI for not precomputed values raises Exception
        r = 500
        self.assertRaises(ValueError, lladser_ci_from_r, r=r, t=t, f=f,
                          alpha=alpha, ci_type=ci_type)


if __name__ == '__main__':
    main()
