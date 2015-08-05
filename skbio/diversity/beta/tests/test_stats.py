# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
from six import StringIO

from skbio import DistanceMatrix, TreeNode
from skbio.diversity.beta._stats import unweighted_unifrac, weighted_unifrac

class StatsTests(TestCase):

    def setUp(self):
        self.b1 = np.array(
           [[1, 3, 0, 1, 0],
            [0, 2, 0, 4, 4],
            [0, 0, 6, 2, 1],
            [0, 0, 1, 1, 1],
            [5, 3, 5, 0, 0],
            [0, 0, 0, 3, 5]])
        self.sids1 = list('ABCDEF')
        self.oids1 = ['OTU%d' % i for i in range(1,6)]
        self.t1 = TreeNode.read(StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:'
            '1.0):1.0):0.0,(OTU4:0.75,OTU5:0.75):1.25):0.0)root;'))

    def test_unweighted_unifrac_identity(self):
        for i in range(len(self.b1)):
            actual = unweighted_unifrac(
                self.b1[i], self.b1[i], self.oids1, self.t1)
            expected = 0.0
            self.assertEqual(actual, expected)

    def test_unweighted_unifrac_symmetry(self):
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = unweighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                expected = unweighted_unifrac(
                    self.b1[j], self.b1[i], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_unweighted_unifrac_non_overlapping(self):
        # these communities only share the root node
        actual = unweighted_unifrac(
            self.b1[4], self.b1[5], self.oids1, self.t1)
        expected = 1.0
        self.assertEqual(actual, expected)
        actual = unweighted_unifrac(
            [1, 1, 1, 0, 0], [0, 0, 0, 1, 1], self.oids1, self.t1)
        expected = 1.0
        self.assertEqual(actual, expected)

    def test_unweighted_unifrac(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation skbio's initial
        # unweighted unifrac implementation
        # sample A versus all
        actual = unweighted_unifrac(
            self.b1[0], self.b1[1], self.oids1, self.t1)
        expected = 0.238095238095
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.b1[0], self.b1[2], self.oids1, self.t1)
        expected = 0.52
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.b1[0], self.b1[3], self.oids1, self.t1)
        expected = 0.52
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.b1[0], self.b1[4], self.oids1, self.t1)
        expected = 0.545454545455
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.b1[0], self.b1[5], self.oids1, self.t1)
        expected = 0.619047619048
        self.assertAlmostEqual(actual, expected)
        # sample B versus remaining
        actual = unweighted_unifrac(
            self.b1[1], self.b1[2], self.oids1, self.t1)
        expected = 0.347826086957
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.b1[1], self.b1[3], self.oids1, self.t1)
        expected = 0.347826086957
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.b1[1], self.b1[4], self.oids1, self.t1)
        expected = 0.68
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.b1[1], self.b1[5], self.oids1, self.t1)
        expected = 0.421052631579
        self.assertAlmostEqual(actual, expected)
        # sample C versus remaining
        actual = unweighted_unifrac(
            self.b1[2], self.b1[3], self.oids1, self.t1)
        expected = 0.0
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.b1[2], self.b1[4], self.oids1, self.t1)
        expected = 0.68
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.b1[2], self.b1[5], self.oids1, self.t1)
        expected = 0.421052631579
        self.assertAlmostEqual(actual, expected)
        # sample D versus remaining
        actual = unweighted_unifrac(
            self.b1[3], self.b1[4], self.oids1, self.t1)
        expected = 0.68
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.b1[3], self.b1[5], self.oids1, self.t1)
        expected = 0.421052631579
        self.assertAlmostEqual(actual, expected)
        # sample E versus remaining
        actual = unweighted_unifrac(
            self.b1[4], self.b1[5], self.oids1, self.t1)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_identity(self):
        for i in range(len(self.b1)):
            actual = weighted_unifrac(
                self.b1[i], self.b1[i], self.oids1, self.t1)
            expected = 0.0
            self.assertEqual(actual, expected)

    def test_weighted_unifrac_symmetry(self):
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = weighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                expected = weighted_unifrac(
                    self.b1[j], self.b1[i], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_non_overlapping(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation skbio's initial
        # weighted unifrac implementation
        # these communities only share the root node
        actual = weighted_unifrac(
            self.b1[4], self.b1[5], self.oids1, self.t1)
        expected = 4.0
        self.assertEqual(actual, expected)

    def test_weighted_unifrac(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation skbio's initial
        # weighted unifrac implementation
        actual = weighted_unifrac(
            self.b1[0], self.b1[1], self.oids1, self.t1)
        expected = 2.4
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[0], self.b1[2], self.oids1, self.t1)
        expected = 1.86666666667
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[0], self.b1[3], self.oids1, self.t1)
        expected = 2.53333333333
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[0], self.b1[4], self.oids1, self.t1)
        expected = 1.35384615385
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[0], self.b1[5], self.oids1, self.t1)
        expected = 3.2
        self.assertAlmostEqual(actual, expected)
        # sample B versus remaining
        actual = weighted_unifrac(
            self.b1[1], self.b1[2], self.oids1, self.t1)
        expected = 2.26666666667
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[1], self.b1[3], self.oids1, self.t1)
        expected = 0.933333333333
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[1], self.b1[4], self.oids1, self.t1)
        expected = 3.2
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[1], self.b1[5], self.oids1, self.t1)
        expected = 0.8375
        self.assertAlmostEqual(actual, expected)
        # sample C versus remaining
        actual = weighted_unifrac(
            self.b1[2], self.b1[3], self.oids1, self.t1)
        expected = 1.33333333333
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[2], self.b1[4], self.oids1, self.t1)
        expected = 1.89743589744
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[2], self.b1[5], self.oids1, self.t1)
        expected = 2.66666666667
        self.assertAlmostEqual(actual, expected)
        # sample D versus remaining
        actual = weighted_unifrac(
            self.b1[3], self.b1[4], self.oids1, self.t1)
        expected = 2.66666666667
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[3], self.b1[5], self.oids1, self.t1)
        expected = 1.33333333333
        self.assertAlmostEqual(actual, expected)
        # sample E versus remaining
        actual = weighted_unifrac(
            self.b1[4], self.b1[5], self.oids1, self.t1)
        expected = 4.0
        self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_identity_normalized(self):
        for i in range(len(self.b1)):
            actual = weighted_unifrac(
                self.b1[i], self.b1[i], self.oids1, self.t1, normalized=True)
            expected = 0.0
            self.assertEqual(actual, expected)

    def test_weighted_unifrac_symmetry_normalized(self):
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = weighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1,
                    normalized=True)
                expected = weighted_unifrac(
                    self.b1[j], self.b1[i], self.oids1, self.t1,
                    normalized=True)
                self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_non_overlapping_normalized(self):
        # these communities only share the root node
        actual = weighted_unifrac(
        self.b1[4], self.b1[5], self.oids1, self.t1, normalized=True)
        expected = 1.0
        self.assertEqual(actual, expected)
        actual = weighted_unifrac(
            [1, 1, 1, 0, 0], [0, 0, 0, 1, 1], self.oids1, self.t1,
            normalized=True)
        expected = 1.0
        self.assertEqual(actual, expected)

    def test_weighted_unifrac_normalized(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation skbio's initial
        # weighted unifrac implementation
        actual = weighted_unifrac(
            self.b1[0], self.b1[1], self.oids1, self.t1, normalized=True)
        expected = 0.6
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[0], self.b1[2], self.oids1, self.t1, normalized=True)
        expected = 0.466666666667
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[0], self.b1[3], self.oids1, self.t1, normalized=True)
        expected = 0.633333333333
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[0], self.b1[4], self.oids1, self.t1, normalized=True)
        expected = 0.338461538462
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[0], self.b1[5], self.oids1, self.t1, normalized=True)
        expected = 0.8
        self.assertAlmostEqual(actual, expected)
        # sample B versus remaining
        actual = weighted_unifrac(
            self.b1[1], self.b1[2], self.oids1, self.t1, normalized=True)
        expected = 0.566666666667
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[1], self.b1[3], self.oids1, self.t1, normalized=True)
        expected = 0.233333333333
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[1], self.b1[4], self.oids1, self.t1, normalized=True)
        expected = 0.8
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[1], self.b1[5], self.oids1, self.t1, normalized=True)
        expected = 0.209375
        self.assertAlmostEqual(actual, expected)
        # sample C versus remaining
        actual = weighted_unifrac(
            self.b1[2], self.b1[3], self.oids1, self.t1, normalized=True)
        expected = 0.333333333333
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[2], self.b1[4], self.oids1, self.t1, normalized=True)
        expected = 0.474358974359
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[2], self.b1[5], self.oids1, self.t1, normalized=True)
        expected = 0.666666666667
        self.assertAlmostEqual(actual, expected)
        # sample D versus remaining
        actual = weighted_unifrac(
            self.b1[3], self.b1[4], self.oids1, self.t1, normalized=True)
        expected = 0.666666666667
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.b1[3], self.b1[5], self.oids1, self.t1, normalized=True)
        expected = 0.333333333333
        self.assertAlmostEqual(actual, expected)
        # sample E versus remaining
        actual = weighted_unifrac(
            self.b1[4], self.b1[5], self.oids1, self.t1, normalized=True)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)
