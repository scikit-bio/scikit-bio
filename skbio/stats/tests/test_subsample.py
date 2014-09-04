# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
try:
    # future >= 0.12
    from future.backports.test.support import import_fresh_module
except ImportError:
    from future.standard_library.test.support import import_fresh_module

import unittest
import warnings

import numpy as np


cy_subsample = import_fresh_module('skbio.stats._subsample',
                                   fresh=['skbio.stats.__subsample'])
py_subsample = import_fresh_module('skbio.stats._subsample',
                                   blocked=[
                                       'skbio.stats.__subsample'
                                   ])


def setup():
    """Ignore warnings during tests."""
    warnings.simplefilter("ignore")


def teardown():
    """Clear the list of warning filters, so that no filters are active."""
    warnings.resetwarnings()


class SubsampleTests(object):

    def test_subsample_nonrandom(self):
        """Should function correctly for nonrandom cases."""
        a = np.array([0, 5, 0])

        # Subsample same number of items that are in input (without
        # replacement).
        np.testing.assert_equal(self.module.subsample(a, 5), a)

        # Can only choose from one bin.
        exp = np.array([0, 2, 0])
        np.testing.assert_equal(self.module.subsample(a, 2), exp)
        np.testing.assert_equal(self.module.subsample(a, 2, replace=True), exp)

        # Subsample zero items.
        a = [3, 0, 1]
        exp = np.array([0, 0, 0])
        np.testing.assert_equal(self.module.subsample(a, 0), exp)
        np.testing.assert_equal(self.module.subsample(a, 0, replace=True), exp)

    def test_subsample_without_replacement(self):
        """Should return a random subsample (without replacement)."""
        # Selecting 2 counts from the vector 1000 times yields each of the two
        # possible results at least once each.
        a = np.array([2, 0, 1])
        actual = set()
        for i in range(1000):
            obs = self.module.subsample(a, 2)
            actual.add(tuple(obs))
        self.assertEqual(actual, {(1, 0, 1), (2, 0, 0)})

        obs = self.module.subsample(a, 2)
        self.assertTrue(np.array_equal(obs, np.array([1, 0, 1])) or
                        np.array_equal(obs, np.array([2, 0, 0])))

    def test_subsample_with_replacement(self):
        """Should return a random subsample (with replacement)."""
        # Can choose from all in first bin, all in last bin (since we're
        # sampling with replacement), or split across bins.
        a = np.array([2, 0, 1])
        actual = set()
        for i in range(1000):
            obs = self.module.subsample(a, 2, replace=True)
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
            obs = self.module.subsample(a, 35, replace=True)
            self.assertEqual(obs.sum(), 35)
            actual.add(tuple(obs))
        self.assertTrue(len(actual) > 10)

    def test_subsample_with_replacement_equal_n(self):
        """Returns random subsample (w/ replacement) when n == counts.sum()."""
        a = np.array([0, 0, 3, 4, 2, 1])
        actual = set()
        for i in range(1000):
            obs = self.module.subsample(a, 10, replace=True)
            self.assertEqual(obs.sum(), 10)
            actual.add(tuple(obs))
        self.assertTrue(len(actual) > 1)

    def test_subsample_invalid_input(self):
        """Should raise an error on invalid input."""
        # Negative n.
        with self.assertRaises(ValueError):
            self.module.subsample([1, 2, 3], -1)

        # Floats.
        with self.assertRaises(TypeError):
            self.module.subsample([1, 2.3, 3], 2)

        # Wrong number of dimensions.
        with self.assertRaises(ValueError):
            self.module.subsample([[1, 2, 3], [4, 5, 6]], 2)

        # Input has too few counts.
        with self.assertRaises(ValueError):
            self.module.subsample([0, 5, 0], 6)


class PySubsampleTests(SubsampleTests, unittest.TestCase):
    module = py_subsample


@unittest.skipIf(cy_subsample is None,
                 "Accelerated subsample module unavailable.")
class CySubsampleTests(SubsampleTests, unittest.TestCase):
    module = cy_subsample

if __name__ == '__main__':
    import nose
    nose.runmodule()
