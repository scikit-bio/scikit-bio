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

from skbio.maths.subsample import subsample


class SubsampleTests(TestCase):

    def test_subsample(self):
        """Should return a random subsample of a vector."""
        a = np.array([0, 5, 0])
        np.testing.assert_equal(subsample(a, 6), a)
        np.testing.assert_equal(subsample(a, 5), a)
        np.testing.assert_equal(subsample(a, 2), np.array([0, 2, 0]))

        # Selecting 2 counts from the vector 1000 times yields each of the two
        # possible results at least once each.
        b = np.array([2, 0, 1])
        actual = {}
        for i in range(1000):
            e = subsample(b, 2)
            actual[tuple(e)] = None
        self.assertEqual(actual, {(1, 0, 1): None, (2, 0, 0): None})

        obs = subsample(b, 2)
        self.assertTrue(np.array_equal(obs, np.array([1, 0, 1])) or
                        np.array_equal(obs, np.array([2, 0, 0])))

    def test_subsample_invalid_input(self):
        """Should raise an error on invalid input."""
        # Wrong number of dimensions.
        with self.assertRaises(ValueError):
            _ = subsample([[1, 2, 3], [4, 5, 6]], 2)

        # Floats.
        with self.assertRaises(TypeError):
            _ = subsample([1, 2.3, 3], 2)


if __name__ == '__main__':
    main()
