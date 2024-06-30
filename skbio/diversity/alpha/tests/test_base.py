# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from io import StringIO
import warnings

import numpy as np
import numpy.testing as npt

from skbio import TreeNode
from skbio.diversity.alpha import (
    berger_parker_d, brillouin_d, dominance, doubles, enspie, esty_ci, fisher_alpha,
    goods_coverage, heip_e, hill, inv_simpson, kempton_taylor_q, margalef, mcintosh_d,
    mcintosh_e, menhinick, michaelis_menten_fit, observed_features, observed_otus, osd,
    pielou_e, renyi, robbins, shannon, simpson, simpson_d, simpson_e, singles, sobs,
    strong, tsallis)


class BaseTests(TestCase):
    def setUp(self):
        self.counts = np.array([0, 1, 1, 4, 2, 5, 2, 4, 1, 2])
        self.sids1 = list('ABCD')
        self.oids1 = ['OTU%d' % i for i in range(1, 6)]
        self.t1 = TreeNode.read(StringIO(
            '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):'
            '0.0,(OTU4:0.75,OTU5:0.75):1.25):0.0)root;'))
        self.t1_w_extra_tips = TreeNode.read(
           StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                    '0.75,(OTU5:0.25,(OTU6:0.5,OTU7:0.5):0.5):0.5):1.25):0.0'
                    ')root;'))

    def test_berger_parker_d(self):
        self.assertEqual(berger_parker_d(np.array([5, 5])), 0.5)
        self.assertEqual(berger_parker_d(np.array([1, 1, 1, 1, 0])), 0.25)
        self.assertEqual(berger_parker_d(self.counts), 5 / 22)
        self.assertEqual(berger_parker_d(np.array([5])), 1)
        self.assertTrue(np.isnan(berger_parker_d([0, 0])))

    def test_brillouin_d(self):
        self.assertAlmostEqual(brillouin_d(np.array([1, 2, 0, 0, 3, 1])),
                               0.86289353018248782)
        self.assertTrue(np.isnan(brillouin_d([0, 0])))

    def test_dominance(self):
        self.assertEqual(dominance(np.array([5])), 1)
        self.assertAlmostEqual(dominance(np.array([1, 0, 2, 5, 2])), 0.34)
        self.assertTrue(np.isnan(dominance([0, 0])))

        # finite sample correction
        self.assertEqual(dominance(np.array([5]), finite=True), 1)
        self.assertAlmostEqual(dominance(
            np.array([1, 0, 2, 5, 2]), finite=True), 0.8 / 3)
        self.assertTrue(np.isnan(dominance([0, 0], finite=True)))

    def test_doubles(self):
        self.assertEqual(doubles(self.counts), 3)
        self.assertEqual(doubles(np.array([0, 3, 4])), 0)
        self.assertEqual(doubles(np.array([2])), 1)
        self.assertEqual(doubles([0, 0]), 0)

    def test_enspie(self):
        for vec in (
            np.array([1, 1, 1, 1, 1, 1]),
            np.array([1, 41, 0, 0, 12, 13]),
            np.array([1, 0, 2, 5, 2])
        ):
            self.assertEqual(enspie(vec), inv_simpson(vec))
        vec = np.array([1, 2, 3, 4])
        self.assertEqual(enspie(vec, finite=True),
                         inv_simpson(vec, finite=True))

    def test_esty_ci(self):
        def _diversity(indices, f):
            """Calculate diversity index for each window of size 1.

            indices: vector of indices of taxa
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

        self.assertTrue(np.isnan(esty_ci([0, 0])))

    def test_fisher_alpha(self):
        exp = 2.7823796
        arr = np.array([4, 3, 4, 0, 1, 0, 2])
        obs = fisher_alpha(arr)
        self.assertAlmostEqual(obs, exp, places=6)

        # Should depend only on S and N (number of taxa, number of
        # individuals / seqs), so we should obtain the same output as above.
        obs = fisher_alpha([1, 6, 1, 0, 1, 0, 5])
        self.assertAlmostEqual(obs, exp, places=6)

        # Should match another by hand:
        # 2 taxa, 62 seqs, alpha is 0.39509
        obs = fisher_alpha([61, 0, 0, 1])
        self.assertAlmostEqual(obs, 0.3950909, places=6)

        # Test case where we have >1000 individuals (SDR-IV makes note of this
        # case). Verified against R's vegan::fisher.alpha.
        obs = fisher_alpha([999, 0, 10])
        self.assertAlmostEqual(obs, 0.2396492, places=6)

        # Should be infinite when all species are singletons
        obs = fisher_alpha([1, 1, 1, 1, 1])
        self.assertEqual(obs, np.inf)

        # Should be large when most species are singletons
        obs = fisher_alpha([1] * 99 + [2])
        self.assertAlmostEqual(obs, 5033.278, places=3)

        # Similar but even larger
        obs = fisher_alpha([1] * 999 + [2])
        TestCase().assertAlmostEqual(obs, 500333.3, places=1)

        self.assertTrue(np.isnan(fisher_alpha([0, 0])))

    def test_goods_coverage(self):
        counts = [1] * 75 + [2, 2, 2, 2, 2, 2, 3, 4, 4]
        obs = goods_coverage(counts)
        self.assertAlmostEqual(obs, 0.23469387755)
        self.assertTrue(np.isnan(goods_coverage([0, 0])))

    def test_heip_e(self):
        # Calculate "by hand".
        arr = np.array([1, 2, 3, 1])
        H = shannon(arr)
        expected = (np.exp(H) - 1) / (arr.size - 1)
        self.assertEqual(heip_e(arr), expected)

        # From Statistical Ecology: A Primer in Methods and Computing, page 94,
        # table 8.1.
        self.assertAlmostEqual(heip_e([500, 300, 200]), 0.90, places=2)
        self.assertAlmostEqual(heip_e([500, 299, 200, 1]), 0.61, places=2)

        # Edge cases
        self.assertEqual(heip_e([5]), 1)
        self.assertTrue(np.isnan(heip_e([0, 0])))

    def test_hill(self):
        orders = [0, 0.5, 1, 2, 10, np.inf]

        # a regular case
        arr = np.array([1, 2, 3, 4, 5])
        obs = [hill(arr, order=x) for x in orders]
        exp = [5, 4.68423304, 4.43598780, 4.09090909, 3.34923645, 3]
        npt.assert_almost_equal(obs, exp)

        # equivalent to observed species richness when q = 0
        self.assertAlmostEqual(hill(arr, order=0), sobs(arr))

        # equivalent to the exponential of Shannon index when q = 1
        self.assertAlmostEqual(hill(arr, order=1), shannon(arr, exp=True))

        # equivalent to inverse Simpson index when q = 2 (default)
        self.assertAlmostEqual(hill(arr), inv_simpson(arr))

        # equivalent to the inverse of Berger-Parker dominance index when q = inf
        self.assertAlmostEqual(hill(arr, order=np.inf), 1 / berger_parker_d(arr))

        # equally abundant taxa: qD = S
        arr = np.array([5, 5, 5])
        obs = [hill(arr, order=x) for x in orders]
        exp = [arr.size] * 6
        npt.assert_almost_equal(obs, exp)

        # single taxon: qD = 1
        self.assertEqual(hill([1]), 1)

        # empty community
        self.assertTrue(np.isnan(hill([0, 0])))

    def test_inv_simpson(self):
        # Totally even community should have 1 / D = number of taxa.
        self.assertAlmostEqual(inv_simpson(np.array([1, 1, 1, 1, 1, 1])), 6)
        self.assertAlmostEqual(inv_simpson(np.array([13, 13, 13, 13])), 4)

        # Hand calculated.
        arr = np.array([1, 41, 0, 0, 12, 13])
        exp = 1 / ((arr / arr.sum()) ** 2).sum()
        self.assertAlmostEqual(inv_simpson(arr), exp)

        # Using dominance.
        exp = 1 / dominance(arr)
        self.assertAlmostEqual(inv_simpson(arr), exp)

        arr = np.array([1, 0, 2, 5, 2])
        exp = 1 / dominance(arr)
        self.assertAlmostEqual(inv_simpson(arr), exp)

        # Finite sample correction.
        self.assertEqual(inv_simpson(
            np.array([1, 0, 2, 5, 2]), finite=True), 3.75)
        self.assertEqual(inv_simpson(np.array([3, 3, 3]), finite=True), 4)

        self.assertTrue(np.isnan(inv_simpson([0, 0])))

    def test_kempton_taylor_q(self):
        # Approximate Magurran 1998 calculation p143.
        arr = np.array([2, 3, 3, 3, 3, 3, 4, 4, 4, 6, 6, 7, 7, 9, 9, 11, 14,
                        15, 15, 20, 29, 33, 34, 36, 37, 53, 57, 138, 146, 170])
        exp = 14 / np.log(34 / 4)
        self.assertAlmostEqual(kempton_taylor_q(arr), exp)

        # Should get same answer regardless of input order.
        np.random.shuffle(arr)
        self.assertAlmostEqual(kempton_taylor_q(arr), exp)

        self.assertTrue(np.isnan(kempton_taylor_q([0, 0])))

    def test_margalef(self):
        self.assertEqual(margalef(self.counts), 8 / np.log(22))
        self.assertTrue(np.isnan(margalef([1])))
        self.assertTrue(np.isnan(margalef([0, 0])))

    def test_mcintosh_d(self):
        self.assertAlmostEqual(mcintosh_d(np.array([1, 2, 3])),
                               0.636061424871458)
        self.assertTrue(np.isnan(mcintosh_d([1])))
        self.assertTrue(np.isnan(mcintosh_d([0, 0])))

    def test_mcintosh_e(self):
        num = np.sqrt(15)
        den = np.sqrt(19)
        exp = num / den
        self.assertEqual(mcintosh_e(np.array([1, 2, 3, 1])), exp)
        self.assertTrue(np.isnan(mcintosh_e([0, 0])))

    def test_menhinick(self):
        # observed species richness = 9, total # of individuals = 22
        self.assertEqual(menhinick(self.counts), 9 / np.sqrt(22))
        self.assertTrue(np.isnan(menhinick([0, 0])))

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
        # [0,100,200,300] looks like only 3 taxa.
        self.assertAlmostEqual(obs_many, 3.0, places=1)
        # [0,2,4,6] looks like 3 taxa with maybe more to be found.
        self.assertTrue(obs_few > obs_many)

        self.assertTrue(np.isnan(michaelis_menten_fit([0, 0])))

    def test_observed_features(self):
        for vec in (np.array([4, 3, 4, 0, 1, 0, 2]), self.counts):
            self.assertEqual(observed_features(vec), sobs(vec))

    def test_observed_otus(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for vec in (np.array([4, 3, 4, 0, 1, 0, 2]), self.counts):
                self.assertEqual(observed_otus(vec), sobs(vec))

    def test_osd(self):
        self.assertEqual(osd(self.counts), (9, 3, 3))

    def test_pielou_e(self):
        # Calculate "by hand".
        arr = np.array([1, 2, 3, 1])
        H = shannon(arr)
        S = arr.size
        expected = H / np.log(S)
        self.assertAlmostEqual(pielou_e(arr), expected)

        # alternative logarithm base
        expected = shannon(arr, base=2) / np.log2(S)
        self.assertAlmostEqual(pielou_e(arr, base=2), expected)

        self.assertAlmostEqual(pielou_e(self.counts), 0.92485490560)

        self.assertAlmostEqual(pielou_e([1, 1]), 1.0)
        self.assertAlmostEqual(pielou_e([1, 1, 1, 1]), 1.0)
        self.assertAlmostEqual(pielou_e([1, 1, 1, 1, 0, 0]), 1.0)

        # Examples from
        # http://ww2.mdsg.umd.edu/interactive_lessons/biofilm/diverse.htm#3
        self.assertAlmostEqual(pielou_e([1, 1, 196, 1, 1]), 0.078, 3)

        # Edge cases
        self.assertEqual(pielou_e([5]), 1)
        self.assertTrue(np.isnan(pielou_e([0, 0])))

    def test_renyi(self):
        orders = [0, 0.5, 1, 2, 10, np.inf]

        # a regular case
        arr = np.array([1, 2, 3, 4, 5])
        obs = [renyi(arr, order=x) for x in orders]
        exp = [1.60943791, 1.54420220, 1.48975032,
               1.40876722, 1.20873239, 1.09861229]
        npt.assert_almost_equal(obs, exp)

        # equivalent to Shannon index when q = 1
        self.assertAlmostEqual(renyi(arr, order=1), shannon(arr))

        # equivalent to log(inverse Simpson index) when q = 2 (default)
        self.assertAlmostEqual(renyi(arr), np.log(inv_simpson(arr)))

        # default q, custom logarithm base
        self.assertAlmostEqual(renyi(arr, base=2), 2.03242148)

        # equally abundant taxa: qH = log(S)
        arr = np.array([5, 5, 5])
        obs = [renyi(arr, order=x) for x in orders]
        exp = [np.log(arr.size)] * 6
        npt.assert_almost_equal(obs, exp)

        # single taxon: qH = 0
        self.assertEqual(renyi([1]), 0)

        # empty community
        self.assertTrue(np.isnan(renyi([0, 0])))

    def test_robbins(self):
        self.assertEqual(robbins(np.array([1, 2, 3, 0, 1])), 2 / 7)
        self.assertTrue(np.isnan(robbins([0, 0])))

    def test_shannon(self):
        self.assertAlmostEqual(shannon([5, 5]), 0.693147181)
        self.assertEqual(shannon([5, 5], base=2), 1)
        self.assertAlmostEqual(shannon([5, 5], base=10), 0.301029996)

        # taxa with 0 counts are excluded from calculation
        self.assertAlmostEqual(shannon([1, 2, 3, 4]), 1.279854226)
        self.assertAlmostEqual(shannon([0, 1, 2, 3, 4]), 1.279854226)

        # Shannon index of a single-taxon community is always 0
        self.assertEqual(shannon(np.array([5])), 0)

        # Shannon index cannot be calculated for an empty community
        self.assertTrue(np.isnan(shannon([0, 0])))

        # NaN still holds if input is empty (instead of 0's), this behavior is
        # different from scipy.stats.entropy, which would return 0.0.
        self.assertTrue(np.isnan(shannon([])))

        # Exponential of Shannon index
        self.assertAlmostEqual(shannon([1, 2, 3, 4], exp=True), 3.596115467)

        # Equally abundant taxa, exp(H) = # taxa
        self.assertAlmostEqual(shannon([5, 5, 5], exp=True), 3.0)

    def test_simpson(self):
        self.assertAlmostEqual(simpson(np.array([1, 0, 2, 5, 2])), 0.66)
        self.assertEqual(simpson(np.array([5])), 0)
        self.assertEqual(simpson(np.array([5]), finite=True), 0)
        self.assertTrue(np.isnan(simpson([0, 0])))

    def  test_simpson_d(self):
        for vec in (np.array([5]), np.array([1, 0, 2, 5, 2])):
            self.assertEqual(dominance(vec), simpson_d(vec))
            self.assertEqual(dominance(vec, finite=True),
                             simpson_d(vec, finite=True))

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

        self.assertTrue(np.isnan(simpson_e([0, 0])))

    def test_singles(self):
        self.assertEqual(singles(self.counts), 3)
        self.assertEqual(singles(np.array([0, 3, 4])), 0)
        self.assertEqual(singles(np.array([1])), 1)
        self.assertEqual(singles([0, 0]), 0)

    def test_sobs(self):
        obs = sobs(np.array([4, 3, 4, 0, 1, 0, 2]))
        self.assertEqual(obs, 5)

        obs = sobs(np.array([0, 0, 0]))
        self.assertEqual(obs, 0)

        obs = sobs(self.counts)
        self.assertEqual(obs, 9)

    def test_strong(self):
        self.assertAlmostEqual(strong(np.array([1, 2, 3, 1])), 0.214285714)
        self.assertTrue(np.isnan(strong([0, 0])))

    def test_tsallis(self):
        orders = [0, 0.5, 1, 2, 10, np.inf]

        # a regular case
        arr = np.array([1, 2, 3, 4, 5])
        obs = [tsallis(arr, order=x) for x in orders]
        exp = [4, 2.32861781, 1.48975032, 0.75555556, 0.11110902, 0]
        npt.assert_almost_equal(obs, exp)

        # equivalent to richess - 1 when q = 0
        self.assertAlmostEqual(tsallis(arr, order=0), sobs(arr) - 1)

        # equivalent to Shannon index when q = 1
        self.assertAlmostEqual(tsallis(arr, order=1), shannon(arr))

        # equivalent to Simpson's diversity index) when q = 2 (default)
        self.assertAlmostEqual(tsallis(arr), simpson(arr))

        # 0 when order is infinity
        self.assertAlmostEqual(tsallis(arr, order=np.inf), 0)

        # 0 when there is a single taxon
        self.assertEqual(tsallis([1]), 0)

        # empty community
        self.assertTrue(np.isnan(tsallis([0, 0])))


if __name__ == '__main__':
    main()
