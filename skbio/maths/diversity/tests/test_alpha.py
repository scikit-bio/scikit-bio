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
                                         observed_species, osd)


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


if __name__ == '__main__':
    main()
