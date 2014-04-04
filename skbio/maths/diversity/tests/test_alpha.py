#!/usr/bin/env python
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np

from skbio.maths.diversity.alpha import (berger_parker_d, brillouin_d,
                                         dominance, doubles, enspie,
                                         equitability, mcintosh_d, mcintosh_e,
                                         menhinick, observed_species, osd,
                                         shannon, simpson, simpson_e,
                                         simpson_reciprocal, singles, strong)


class AlphaDiversityTests(TestCase):
    def setUp(self):
        self.counts = np.array([0, 1, 1, 4, 2, 5, 2, 4, 1, 2])

    def test_berger_parker_d(self):
        """Should match hand-calculated values."""
        self.assertEqual(berger_parker_d(np.array([5])), 1)
        self.assertEqual(berger_parker_d(np.array([5, 5])), 0.5)
        self.assertEqual(berger_parker_d(np.array([1, 1, 1, 1, 0])), 0.25)
        self.assertEqual(berger_parker_d(self.counts), 5 / 22)

    def test_brillouin_d(self):
        """Should match hand-calculated values."""
        self.assertAlmostEqual(brillouin_d(np.array([1, 2, 0, 0, 3, 1])),
                               0.86289353018248782)

    def test_dominance(self):
        """Should match hand-calculated values."""
        self.assertEqual(dominance(np.array([5])), 1)
        self.assertAlmostEqual(dominance(np.array([1, 0, 2, 5, 2])), 0.34)

    def test_doubles(self):
        """Should return correct number of doubles."""
        self.assertEqual(doubles(self.counts), 3)
        self.assertEqual(doubles(np.array([0, 3, 4])), 0)
        self.assertEqual(doubles(np.array([2])), 1)

    def test_enspie(self):
        """Test that ENS_pie is correctly calculated."""
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
        """Should match (1 / dominance) results."""
        arr = np.array([1, 0, 2, 5, 2])
        self.assertAlmostEqual(simpson_reciprocal(arr), 1 / dominance(arr))

    def test_equitability(self):
        """Should match hand-calculated values."""
        self.assertAlmostEqual(equitability(np.array([5, 5])), 1)
        self.assertAlmostEqual(equitability(np.array([1, 1, 1, 1, 0])), 1)

    def test_mcintosh_d(self):
        """Should match hand-calculated values."""
        self.assertAlmostEqual(mcintosh_d(np.array([1, 2, 3])),
                               0.636061424871458)

    def test_mcintosh_e(self):
        """Should match hand-calculated results."""
        num = np.sqrt(15)
        den = np.sqrt(19)
        exp = num / den
        self.assertEqual(mcintosh_e(np.array([1, 2, 3, 1])), exp)

    def test_menhinick(self):
        """Should match hand-calculated values"""
        self.assertEqual(menhinick(self.counts), 9 / np.sqrt(22))

    def test_observed_species(self):
        """Should return number of observed species."""
        obs = observed_species(np.array([4, 3, 4, 0, 1, 0, 2]))
        self.assertEqual(obs, 5)

        obs = observed_species(np.array([0, 0, 0]))
        self.assertEqual(obs, 0)

        obs = observed_species(self.counts)
        self.assertEqual(obs, 9)

    def test_osd(self):
        """Should return correct number of observed, singles, and doubles."""
        self.assertEqual(osd(self.counts), (9, 3, 3))

    def test_shannon(self):
        """Should match hand-calculated values."""
        self.assertEqual(shannon(np.array([5])), 0)
        self.assertEqual(shannon(np.array([5, 5])), 1)
        self.assertEqual(shannon(np.array([1, 1, 1, 1, 0])), 2)

    def test_simpson(self):
        """Should match hand-calculated values."""
        self.assertAlmostEqual(simpson(np.array([1, 0, 2, 5, 2])), 0.66)
        self.assertAlmostEqual(simpson(np.array([5])), 0)

    def test_simpson_e(self):
        """Should match hand-calculated values."""
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
        """Should return correct number of singles."""
        self.assertEqual(singles(self.counts), 3)
        self.assertEqual(singles(np.array([0, 3, 4])), 0)
        self.assertEqual(singles(np.array([1])), 1)

    def test_strong(self):
        """Strong's dominance index should match hand-calculated values."""
        self.assertAlmostEqual(strong(np.array([1, 2, 3, 1])), 0.214285714)


if __name__ == '__main__':
    main()
