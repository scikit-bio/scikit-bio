#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

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

from skbio.diversity.alpha import (
    berger_parker_d, brillouin_d, dominance, doubles, enspie, equitability,
    esty_ci, fisher_alpha, goods_coverage, heip_e, kempton_taylor_q, margalef,
    mcintosh_d, mcintosh_e, menhinick, michaelis_menten_fit, observed_otus,
    osd, robbins, shannon, simpson, simpson_e, singles, strong)
from skbio.diversity.alpha._base import _validate


class BaseTests(TestCase):
    def setUp(self):
        self.counts = np.array([0, 1, 1, 4, 2, 5, 2, 4, 1, 2])

    def test_validate(self):
        # python list
        obs = _validate([0, 2, 1, 3])
        npt.assert_array_equal(obs, np.array([0, 2, 1, 3]))
        self.assertEqual(obs.dtype, int)

        # numpy array (no copy made)
        data = np.array([0, 2, 1, 3])
        obs = _validate(data)
        npt.assert_array_equal(obs, data)
        self.assertEqual(obs.dtype, int)
        self.assertTrue(obs is data)

        # single element
        obs = _validate([42])
        npt.assert_array_equal(obs, np.array([42]))
        self.assertEqual(obs.dtype, int)
        self.assertEqual(obs.shape, (1,))

        # suppress casting to int
        obs = _validate([42.2, 42.1, 0], suppress_cast=True)
        npt.assert_array_equal(obs, np.array([42.2, 42.1, 0]))
        self.assertEqual(obs.dtype, float)

        # all zeros
        obs = _validate([0, 0, 0])
        npt.assert_array_equal(obs, np.array([0, 0, 0]))
        self.assertEqual(obs.dtype, int)

        # all zeros (single value)
        obs = _validate([0])
        npt.assert_array_equal(obs, np.array([0]))
        self.assertEqual(obs.dtype, int)

    def test_validate_invalid_input(self):
        # wrong dtype
        with self.assertRaises(TypeError):
            _validate([0, 2, 1.2, 3])

        # wrong number of dimensions (2-D)
        with self.assertRaises(ValueError):
            _validate([[0, 2, 1, 3], [4, 5, 6, 7]])

        # wrong number of dimensions (scalar)
        with self.assertRaises(ValueError):
            _validate(1)

        # negative values
        with self.assertRaises(ValueError):
            _validate([0, 0, 2, -1, 3])

    def test_berger_parker_d(self):
        self.assertEqual(berger_parker_d(np.array([5])), 1)
        self.assertEqual(berger_parker_d(np.array([5, 5])), 0.5)
        self.assertEqual(berger_parker_d(np.array([1, 1, 1, 1, 0])), 0.25)
        self.assertEqual(berger_parker_d(self.counts), 5 / 22)

    def test_brillouin_d(self):
        self.assertAlmostEqual(brillouin_d(np.array([1, 2, 0, 0, 3, 1])),
                               0.86289353018248782)

    def test_dominance(self):
        self.assertEqual(dominance(np.array([5])), 1)
        self.assertAlmostEqual(dominance(np.array([1, 0, 2, 5, 2])), 0.34)

    def test_doubles(self):
        self.assertEqual(doubles(self.counts), 3)
        self.assertEqual(doubles(np.array([0, 3, 4])), 0)
        self.assertEqual(doubles(np.array([2])), 1)
        self.assertEqual(doubles(np.array([0, 0])), 0)

    def test_enspie(self):
        # Totally even community should have ENS_pie = number of OTUs.
        self.assertAlmostEqual(enspie(np.array([1, 1, 1, 1, 1, 1])), 6)
        self.assertAlmostEqual(enspie(np.array([13, 13, 13, 13])), 4)

        # Hand calculated.
        arr = np.array([1, 41, 0, 0, 12, 13])
        exp = 1 / ((arr / arr.sum()) ** 2).sum()
        self.assertAlmostEqual(enspie(arr), exp)

        # Using dominance.
        exp = 1 / dominance(arr)
        self.assertAlmostEqual(enspie(arr), exp)

        arr = np.array([1, 0, 2, 5, 2])
        exp = 1 / dominance(arr)
        self.assertAlmostEqual(enspie(arr), exp)

    def test_equitability(self):
        self.assertAlmostEqual(equitability(np.array([5, 5])), 1)
        self.assertAlmostEqual(equitability(np.array([1, 1, 1, 1, 0])), 1)

    def test_esty_ci(self):
        def _diversity(indices, f):
            """Calculate diversity index for each window of size 1.

            indices: vector of indices of OTUs
            f: f(counts) -> diversity measure

            """
            result = []
            max_size = max(indices) + 1
            freqs = np.zeros(max_size, dtype=int)
            for i in range(len(indices)):
                freqs += np.bincount(indices[i:i + 1], minlength=max_size)
                try:
                    curr = f(freqs)
                except (ZeroDivisionError, FloatingPointError):
                    curr = 0
                result.append(curr)
            return np.array(result)

        data = [1, 1, 2, 1, 1, 3, 2, 1, 3, 4]

        observed_lower, observed_upper = zip(*_diversity(data, esty_ci))

        expected_lower = np.array([1, -1.38590382, -0.73353593, -0.17434465,
                                   -0.15060902, -0.04386191, -0.33042054,
                                   -0.29041008, -0.43554755, -0.33385652])
        expected_upper = np.array([1, 1.38590382, 1.40020259, 0.67434465,
                                   0.55060902, 0.71052858, 0.61613483,
                                   0.54041008, 0.43554755, 0.53385652])

        npt.assert_array_almost_equal(observed_lower, expected_lower)
        npt.assert_array_almost_equal(observed_upper, expected_upper)

    def test_fisher_alpha(self):
        exp = 2.7823795367398798
        arr = np.array([4, 3, 4, 0, 1, 0, 2])
        obs = fisher_alpha(arr)
        self.assertAlmostEqual(obs, exp)

        # Should depend only on S and N (number of OTUs, number of
        # individuals / seqs), so we should obtain the same output as above.
        obs = fisher_alpha([1, 6, 1, 0, 1, 0, 5])
        self.assertAlmostEqual(obs, exp)

        # Should match another by hand:
        # 2 OTUs, 62 seqs, alpha is 0.39509
        obs = fisher_alpha([61, 0, 0, 1])
        self.assertAlmostEqual(obs, 0.39509, delta=0.0001)

        # Test case where we have >1000 individuals (SDR-IV makes note of this
        # case). Verified against R's vegan::fisher.alpha.
        obs = fisher_alpha([999, 0, 10])
        self.assertAlmostEqual(obs, 0.2396492)

    def test_goods_coverage(self):
        counts = [1] * 75 + [2, 2, 2, 2, 2, 2, 3, 4, 4]
        obs = goods_coverage(counts)
        self.assertAlmostEqual(obs, 0.23469387755)

    def test_heip_e(self):
        # Calculate "by hand".
        arr = np.array([1, 2, 3, 1])
        h = shannon(arr, base=np.e)
        expected = (np.exp(h) - 1) / 3
        self.assertEqual(heip_e(arr), expected)

        # From Statistical Ecology: A Primer in Methods and Computing, page 94,
        # table 8.1.
        self.assertAlmostEqual(heip_e([500, 300, 200]), 0.90, places=2)
        self.assertAlmostEqual(heip_e([500, 299, 200, 1]), 0.61, places=2)

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
        # observed_otus = 9, total # of individuals = 22
        self.assertEqual(menhinick(self.counts), 9 / np.sqrt(22))

    def test_michaelis_menten_fit(self):
        obs = michaelis_menten_fit([22])
        self.assertAlmostEqual(obs, 1.0)

        obs = michaelis_menten_fit([42])
        self.assertAlmostEqual(obs, 1.0)

        obs = michaelis_menten_fit([34], num_repeats=3, params_guess=(13, 13))
        self.assertAlmostEqual(obs, 1.0)

        obs = michaelis_menten_fit([70, 70], num_repeats=5)
        self.assertAlmostEqual(obs, 2.0, places=1)

        obs_few = michaelis_menten_fit(np.arange(4) * 2, num_repeats=10)
        obs_many = michaelis_menten_fit(np.arange(4) * 100, num_repeats=10)
        # [0,100,200,300] looks like only 3 OTUs.
        self.assertAlmostEqual(obs_many, 3.0, places=1)
        # [0,2,4,6] looks like 3 OTUs with maybe more to be found.
        self.assertTrue(obs_few > obs_many)

    def test_observed_otus(self):
        obs = observed_otus(np.array([4, 3, 4, 0, 1, 0, 2]))
        self.assertEqual(obs, 5)

        obs = observed_otus(np.array([0, 0, 0]))
        self.assertEqual(obs, 0)

        obs = observed_otus(self.counts)
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
        self.assertEqual(singles(np.array([0, 0])), 0)

    def test_strong(self):
        self.assertAlmostEqual(strong(np.array([1, 2, 3, 1])), 0.214285714)


if __name__ == '__main__':
    main()
