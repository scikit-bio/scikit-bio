#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The bipy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

"""Tests of statistical probability distribution integrals.

Currently using tests against calculations in R, spreadsheets being unreliable.
"""

from bipy.util.unit_test import TestCase, main
from bipy.maths.stats.distribution import (chi_high, z_high, zprob, f_high,
                                           binomial_high, bdtrc, stdtr)


class DistributionsTests(TestCase):

    """Tests of particular statistical distributions."""

    def setUp(self):
        self.values = [0, 0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 30, 50, 200]
        self.negvalues = [-i for i in self.values]
        self.df = [1, 10, 100]

    def test_z_high(self):
        """z_high should match R's pnorm(lower.tail=FALSE) function"""

        negprobs = [
            0.5000000, 0.5039894, 0.5398278, 0.6914625, 0.8413447,
            0.9772499, 0.9999997, 1.0000000, 1.0000000, 1.0000000,
            1.0000000, 1.0000000,
        ]
        probs = [
            5.000000e-01, 4.960106e-01, 4.601722e-01, 3.085375e-01,
            1.586553e-01, 2.275013e-02, 2.866516e-07, 7.619853e-24,
            2.753624e-89, 4.906714e-198, 0.000000e+00, 0.000000e+00]

        for z, p in zip(self.values, probs):
            self.assertFloatEqual(z_high(z), p)
        for z, p in zip(self.negvalues, negprobs):
            self.assertFloatEqual(z_high(z), p)

    def test_zprob(self):
        """zprob should match twice the z_high probability for abs(z)"""

        probs = [2 * i for i in [
            5.000000e-01, 4.960106e-01, 4.601722e-01, 3.085375e-01,
            1.586553e-01, 2.275013e-02, 2.866516e-07, 7.619853e-24,
            2.753624e-89, 4.906714e-198, 0.000000e+00, 0.000000e+00]]

        for z, p in zip(self.values, probs):
            self.assertFloatEqual(zprob(z), p)
        for z, p in zip(self.negvalues, probs):
            self.assertFloatEqual(zprob(z), p)

    def test_chi_high(self):
        """chi_high should match R's pchisq(lower.tail=FALSE) function"""
        probs = {
            1: [1.000000e+00, 9.203443e-01, 7.518296e-01, 4.795001e-01,
                3.173105e-01, 1.572992e-01, 2.534732e-02, 1.565402e-03,
                7.744216e-06, 4.320463e-08, 1.537460e-12, 2.088488e-45,
                ],
            10: [1.000000e+00, 1.000000e-00, 1.000000e-00, 9.999934e-01,
                 9.998279e-01, 9.963402e-01, 8.911780e-01, 4.404933e-01,
                 2.925269e-02, 8.566412e-04, 2.669083e-07, 1.613931e-37,
                 ],
            100: [1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
                  1.00000e+00, 1.00000e+00, 1.00000e+00, 1.00000e+00,
                  1.00000e+00, 1.00000e+00, 9.99993e-01, 1.17845e-08,
                  ],
        }

        for df in self.df:
            for x, p in zip(self.values, probs[df]):
                self.assertFloatEqual(chi_high(x, df), p)

    def test_binomial_high(self):
        """Binomial high should match values from R for integer successes"""
        expected = {
            (0, 1, 0.5): 0.5,
            (1, 1, 0.5): 0,
            (1, 1, 0.0000001): 0,
            (1, 1, 0.9999999): 0,
            (3, 5, 0.75): 0.6328125,
            (0, 60, 0.5): 1,
            (129, 130, 0.5): 7.34684e-40,
            (299, 300, 0.099): 4.904089e-302,
            (9, 27, 0.0003): 4.958496e-29,
            (1032, 2050, 0.5): 0.3702155,
            (-1, 3, 0.1): 1,  # if successes less than 0, return 1
            (-0.5, 3, 0.1): 1,
        }
        for (key, value) in expected.items():
            self.assertFloatEqualRel(binomial_high(*key), value, 1e-4)
        # should reject if successes > trials or successes < -1
        self.assertRaises(ValueError, binomial_high, 7, 5, 0.5)

    def test_f_high(self):
        """F high should match values from R for integer successes"""
        expected = {
            (1, 1, 0): 1,
            (1, 1, 1): 0.5,
            (1, 1, 20): 0.1400487,
            (1, 1, 1000000): 0.0006366196,
            (1, 10, 0): 1,
            (1, 10, 5): 0.0493322,
            (1, 10, 20): 0.001193467,
            (10, 1, 0): 1,
            (10, 10, 14.7): 0.0001062585,
            # test non-integer degrees of freedom
            (13.7, 11.9, 3.8): 0.01340347,
            # used following series to track down a bug after a failed test
            # case
            (28, 29, 2): 0.03424088,
            (28, 29, 10): 1.053019e-08,
            (28, 29, 20): 1.628245e-12,
            (28, 29, 300): 5.038791e-29,
            (28, 35, 1): 0.4946777,
            (28, 37, 1): 0.4934486,
            (28, 38, 1): 0.4928721,
            (28, 38.001, 1): 0.4928716,
            (28, 38.5, 1): 0.4925927,
            (28, 39, 1): 0.492319,
            (28, 39, 10): 1.431901e-10,
            (28, 39, 20): 1.432014e-15,
            (28, 39, 30): 1.059964e-18,
            (28, 39, 50): 8.846678e-23,
            (28, 39, 10): 1.431901e-10,
            (28, 39, 300): 1.226935e-37,
            (28, 39, 50): 8.846678e-23,
            (28, 39, 304.7): 9.08154e-38,
            (28.4, 39.2, 304.7): 5.573927e-38,
            (1032, 2050, 0): 1,
            (1032, 2050, 4.15): 1.23535e-165,
            (1032, 2050, 0.5): 1,
            (1032, 2050, 0.1): 1,
        }
        e = sorted(expected.items())
        for (key, value) in e:
            self.assertFloatEqualRel(f_high(*key), value)

    def test_bdtrc(self):
        """bdtrc should give same results as cephes"""
        k_s = [0, 1, 2, 3, 5]
        n_s = [5, 10, 1000]
        p_s = [1e-10, .1, .5, .9, .999999]

        exp = [
            4.999999999e-10,
            0.40951,
            0.96875,
            0.99999,
            1.0,
            9.9999999955e-10,
            0.6513215599,
            0.9990234375,
            0.9999999999,
            1.0,
            9.9999995005e-08,
            1.0,
            1.0,
            1.0,
            1.0,
            9.999999998e-20,
            0.08146,
            0.8125,
            0.99954,
            1.0,
            4.4999999976e-19,
            0.2639010709,
            0.9892578125,
            0.9999999909,
            1.0,
            4.99499966766e-15,
            1.0,
            1.0,
            1.0,
            1.0,
            9.9999999985e-30,
            0.00856,
            0.5,
            0.99144,
            1.0,
            1.19999999937e-28,
            0.0701908264,
            0.9453125,
            0.9999996264,
            1.0,
            1.66166987575e-22,
            1.0,
            1.0,
            1.0,
            1.0,
            4.9999999996e-40,
            0.00046,
            0.1875,
            0.91854,
            0.99999999999,
            2.09999999899e-38,
            0.0127951984,
            0.828125,
            0.9999908784,
            1.0,
            4.14171214499e-30,
            1.0,
            1.0,
            1.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            2.09999999928e-58,
            0.0001469026,
            0.376953125,
            0.9983650626,
            1.0,
            1.36817318242e-45,
            1.0,
            1.0,
            1.0,
            1.0,
        ]
        index = 0
        for k in k_s:
            for n in n_s:
                for p in p_s:
                    self.assertFloatEqual(bdtrc(k, n, p), exp[index])
                    index += 1

    def test_stdtr(self):
        """stdtr should match cephes results"""
        t = [-10, -3.1, -0.5, -0.01, 0, 1, 0.5, 10]
        k = [2, 10, 100]
        exp = [
            0.00492622851166,
            7.94776587798e-07,
            4.9508444923e-17,
            0.0451003650651,
            0.00562532860804,
            0.00125696358826,
            0.333333333333,
            0.313946802871,
            0.309086782915,
            0.496464554479,
            0.496108987495,
            0.496020605117,
            0.5,
            0.5,
            0.5,
            0.788675134595,
            0.829553433849,
            0.840137922108,
            0.666666666667,
            0.686053197129,
            0.690913217085,
            0.995073771488,
            0.999999205223,
            1.0,
        ]
        index = 0
        for i in t:
            for j in k:
                self.assertFloatEqual(stdtr(j, i), exp[index])
                index += 1


if __name__ == "__main__":
    main()
