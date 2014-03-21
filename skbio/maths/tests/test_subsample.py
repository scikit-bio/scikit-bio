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

    def test_subsample_nonrandom(self):
        """Should function correctly for nonrandom cases."""
        a = np.array([0, 5, 0])

        # Input has too few counts.
        np.testing.assert_equal(subsample(a, 6), a)
        np.testing.assert_equal(subsample(a, 5), a)

        # Can only choose from one bin.
        np.testing.assert_equal(subsample(a, 2), np.array([0, 2, 0]))

    def test_subsample_without_replacement(self):
        """Should return a random subsample (without replacement)."""
        # Selecting 2 counts from the vector 1000 times yields each of the two
        # possible results at least once each.
        a = np.array([2, 0, 1])
        actual = set()
        for i in range(1000):
            obs = subsample(a, 2)
            actual.add(tuple(obs))
        self.assertEqual(actual, {(1, 0, 1), (2, 0, 0)})

        obs = subsample(a, 2)
        self.assertTrue(np.array_equal(obs, np.array([1, 0, 1])) or
                        np.array_equal(obs, np.array([2, 0, 0])))

    def test_subsample_with_replacement(self):
        """Should return a random subsample (with replacement)."""
        # Can choose from all in first bin, all in last bin (since we're
        # sampling with replacement), or split across bins.
        a = np.array([2, 0, 1])
        actual = set()
        for i in range(1000):
            obs = subsample(a, 2, replace=True)
            actual.add(tuple(obs))
        self.assertEqual(actual, {(1, 0, 1), (2, 0, 0), (0, 0, 2)})

        # Test that selecting 35 counts from a 36-count vector 1000 times
        # yields more than 10 different subsamples. If we were subsampling
        # *without* replacement, there would be only 10 possible subsamples
        # because there are 10 nonzero bins in array a. However, there are more
        # than 10 possibilities when sampling *with* replacement.
        a = np.array([2, 0, 1, 2, 1, 8, 6, 0, 3, 3, 5, 0, 0, 0, 5])
        actual = set()
        for i in range(1000):
            obs = subsample(a, 35, replace=True)
            self.assertEqual(obs.sum(), 35)
            actual.add(tuple(obs))
        self.assertTrue(len(actual) > 10)

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
