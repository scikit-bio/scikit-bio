#!/usr/bin/env python
from __future__ import division

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

from skbio.maths.diversity.alpha.gini import (gini_index, _lorenz_curve,
                                              _lorenz_curve_integrator)


class GiniTests(TestCase):
    def setUp(self):
        self.data = np.array([4.5, 6.7, 3.4, 15., 18., 3.5, 6.7, 14.1])
        self.lorenz_curve_points = [(0.125, 0.047287899860917935),
                                    (0.25, 0.095966620305980521),
                                    (0.375, 0.15855354659248957),
                                    (0.5, 0.2517385257301808),
                                    (0.625, 0.34492350486787204),
                                    (0.75, 0.541029207232267),
                                    (0.875, 0.74965229485396379),
                                    (1.0, 1.0)]

    def test_gini_index(self):
        exp = 0.32771210013908214
        obs = gini_index(self.data, 'trapezoids')
        self.assertAlmostEqual(obs, exp)

        exp = 0.20271210013908214
        obs = gini_index(self.data, 'rectangles')
        self.assertAlmostEqual(obs, exp)

        # Raises error on negative data.
        with self.assertRaises(ValueError):
            gini_index([1.0, -3.1, 4.5])

    def test_lorenz_curve(self):
        self.assertEqual(_lorenz_curve(self.data), self.lorenz_curve_points)

    def test_lorenz_curve_integrator(self):
        exp = 0.33614394993045893
        obs = _lorenz_curve_integrator(self.lorenz_curve_points, 'trapezoids')
        self.assertAlmostEqual(obs, exp)

        exp = 0.39864394993045893
        obs = _lorenz_curve_integrator(self.lorenz_curve_points, 'rectangles')
        self.assertAlmostEqual(obs, exp)

        # Raises error on invalid method.
        with self.assertRaises(ValueError):
            _lorenz_curve_integrator(self.lorenz_curve_points, 'brofist')


if __name__ == '__main__':
    main()
