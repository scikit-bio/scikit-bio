# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from io import StringIO
from unittest import main, TestCase

import numpy as np

from skbio import TreeNode
from skbio.tree import DuplicateNodeError, MissingNodeError
from skbio.diversity.beta import unweighted_unifrac, weighted_unifrac
from skbio.diversity.beta._unifrac import (_unweighted_unifrac,
                                           _weighted_unifrac,
                                           _weighted_unifrac_branch_correction)


class UnifracTests(TestCase):

    def setUp(self):
        self.b1 = np.array(
            [[1, 3, 0, 1, 0],
             [0, 2, 0, 4, 4],
             [0, 0, 6, 2, 1],
             [0, 0, 1, 1, 1],
             [5, 3, 5, 0, 0],
             [0, 0, 0, 3, 5]])
        self.sids1 = list('ABCDEF')
        self.oids1 = ['OTU%d' % i for i in range(1, 6)]
        self.t1 = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        self.t1_w_extra_tips = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,(OTU5:0.25,(OTU6:0.5,OTU7:0.5):0.5):0.5):1.25):0.0'
                     ')root;'))

        self.t2 = TreeNode.read(
            StringIO('((OTU1:0.1, OTU2:0.2):0.3, (OTU3:0.5, OTU4:0.7):1.1)'
                     'root;'))
        self.oids2 = ['OTU%d' % i for i in range(1, 5)]

    def test_unweighted_otus_out_of_order(self):
        # UniFrac API does not assert the observations are in tip order of the
        # input tree
        shuffled_ids = self.oids1[:]
        shuffled_b1 = self.b1.copy()

        shuffled_ids[0], shuffled_ids[-1] = shuffled_ids[-1], shuffled_ids[0]
        shuffled_b1[:, [0, -1]] = shuffled_b1[:, [-1, 0]]

        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = unweighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                expected = unweighted_unifrac(
                    shuffled_b1[i], shuffled_b1[j], shuffled_ids, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_weighted_otus_out_of_order(self):
        # UniFrac API does not assert the observations are in tip order of the
        # input tree
        shuffled_ids = self.oids1[:]
        shuffled_b1 = self.b1.copy()

        shuffled_ids[0], shuffled_ids[-1] = shuffled_ids[-1], shuffled_ids[0]
        shuffled_b1[:, [0, -1]] = shuffled_b1[:, [-1, 0]]

        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = weighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                expected = weighted_unifrac(
                    shuffled_b1[i], shuffled_b1[j], shuffled_ids, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_unweighted_extra_tips(self):
        # UniFrac values are the same despite unobserved tips in the tree
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = unweighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1_w_extra_tips)
                expected = unweighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_weighted_extra_tips(self):
        # UniFrac values are the same despite unobserved tips in the tree
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = weighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1_w_extra_tips)
                expected = weighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_unweighted_minimal_trees(self):
        # two tips
        tree = TreeNode.read(StringIO('(OTU1:0.25, OTU2:0.25)root;'))
        actual = unweighted_unifrac([1, 0], [0, 0], ['OTU1', 'OTU2'],
                                    tree)
        expected = 1.0
        self.assertEqual(actual, expected)

    def test_weighted_minimal_trees(self):
        # two tips
        tree = TreeNode.read(StringIO('(OTU1:0.25, OTU2:0.25)root;'))
        actual = weighted_unifrac([1, 0], [0, 0], ['OTU1', 'OTU2'], tree)
        expected = 0.25
        self.assertEqual(actual, expected)

    def test_unweighted_root_not_observed(self):
        # expected values computed with QIIME 1.9.1 and by hand
        # root node not observed, but branch between (OTU1, OTU2) and root
        # is considered shared
        actual = unweighted_unifrac([1, 1, 0, 0], [1, 0, 0, 0],
                                    self.oids2, self.t2)
        # for clarity of what I'm testing, compute expected as it would
        # based on the branch lengths. the values that compose shared was
        # a point of confusion for me here, so leaving these in for
        # future reference
        expected = 0.2 / (0.1 + 0.2 + 0.3)  # 0.3333333333
        self.assertAlmostEqual(actual, expected)

        # root node not observed, but branch between (OTU3, OTU4) and root
        # is considered shared
        actual = unweighted_unifrac([0, 0, 1, 1], [0, 0, 1, 0],
                                    self.oids2, self.t2)
        # for clarity of what I'm testing, compute expected as it would
        # based on the branch lengths. the values that compose shared was
        # a point of confusion for me here, so leaving these in for
        # future reference
        expected = 0.7 / (1.1 + 0.5 + 0.7)  # 0.3043478261
        self.assertAlmostEqual(actual, expected)

    def test_weighted_root_not_observed(self):
        # expected values computed by hand, these disagree with QIIME 1.9.1
        # root node not observed, but branch between (OTU1, OTU2) and root
        # is considered shared
        actual = weighted_unifrac([1, 0, 0, 0], [1, 1, 0, 0],
                                  self.oids2, self.t2)
        expected = 0.15
        self.assertAlmostEqual(actual, expected)

        # root node not observed, but branch between (OTU3, OTU4) and root
        # is considered shared
        actual = weighted_unifrac([0, 0, 1, 1], [0, 0, 1, 0],
                                  self.oids2, self.t2)
        expected = 0.6
        self.assertAlmostEqual(actual, expected)

    def test_weighted_normalized_root_not_observed(self):
        # expected values computed by hand, these disagree with QIIME 1.9.1
        # root node not observed, but branch between (OTU1, OTU2) and root
        # is considered shared
        actual = weighted_unifrac([1, 0, 0, 0], [1, 1, 0, 0],
                                  self.oids2, self.t2, normalized=True)
        expected = 0.1764705882
        self.assertAlmostEqual(actual, expected)

        # root node not observed, but branch between (OTU3, OTU4) and root
        # is considered shared
        actual = weighted_unifrac([0, 0, 1, 1], [0, 0, 1, 0],
                                  self.oids2, self.t2, normalized=True)
        expected = 0.1818181818
        self.assertAlmostEqual(actual, expected)

    def test_unweighted_unifrac_identity(self):
        for i in range(len(self.b1)):
            actual = unweighted_unifrac(
                self.b1[i], self.b1[i], self.oids1, self.t1)
            expected = 0.0
            self.assertAlmostEqual(actual, expected)

    def test_unweighted_unifrac_symmetry(self):
        for i in range(len(self.b1)):
            for j in range(len(self.b1)):
                actual = unweighted_unifrac(
                    self.b1[i], self.b1[j], self.oids1, self.t1)
                expected = unweighted_unifrac(
                    self.b1[j], self.b1[i], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_invalid_input(self):
        # Many of these tests are duplicated from
        # skbio.diversity.tests.test_base, but I think it's important to
        # confirm that they are being run when *unifrac is called.

        # tree has duplicated tip ids
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU2:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(DuplicateNodeError, unweighted_unifrac,
                          u_counts, v_counts, otu_ids, t)
        self.assertRaises(DuplicateNodeError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)

        # unrooted tree as input
        t = TreeNode.read(StringIO('((OTU1:0.1, OTU2:0.2):0.3, OTU3:0.5,'
                                   'OTU4:0.7);'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)

        # otu_ids has duplicated ids
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU2']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)

        # len of vectors not equal
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        u_counts = [1, 2, 3]
        v_counts = [1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)

        # negative counts
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, -3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        u_counts = [1, 2, 3]
        v_counts = [1, 1, -1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)

        # tree with no branch lengths
        t = TreeNode.read(
            StringIO('((((OTU1,OTU2),OTU3)),(OTU4,OTU5));'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)

        # tree missing some branch lengths
        t = TreeNode.read(
            StringIO('(((((OTU1,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)

        # otu_ids not present in tree
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU42']
        self.assertRaises(MissingNodeError, unweighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        self.assertRaises(MissingNodeError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)

    def test_unweighted_unifrac_non_overlapping(self):
        # these communities only share the root node
        actual = unweighted_unifrac(
            self.b1[4], self.b1[5], self.oids1, self.t1)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            [1, 1, 1, 0, 0], [0, 0, 0, 1, 1], self.oids1, self.t1)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)

    def test_unweighted_unifrac_zero_counts(self):
        actual = unweighted_unifrac(
            [1, 1, 1, 0, 0], [0, 0, 0, 0, 0], self.oids1, self.t1)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], self.oids1, self.t1)
        expected = 0.0
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            [], [], [], self.t1)
        expected = 0.0
        self.assertAlmostEqual(actual, expected)

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
            self.assertAlmostEqual(actual, expected)

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
        self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_zero_counts(self):
        actual = weighted_unifrac(
            [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], self.oids1, self.t1)
        expected = 0.0
        self.assertAlmostEqual(actual, expected)
        # calculated the following by hand, as QIIME 1.9.1 tells the user
        # that values involving empty vectors will be uninformative, and
        # returns 1.0
        actual = weighted_unifrac(
            [1, 1, 1, 0, 0], [0, 0, 0, 0, 0], self.oids1, self.t1)
        expected = 2.0
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            [], [], [], self.t1)
        expected = 0.0
        self.assertAlmostEqual(actual, expected)

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
            self.assertAlmostEqual(actual, expected)

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
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            [1, 1, 1, 0, 0], [0, 0, 0, 1, 1], self.oids1, self.t1,
            normalized=True)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_zero_counts_normalized(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation skbio's initial
        # weighted unifrac implementation
        actual = weighted_unifrac(
            [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], self.oids1, self.t1,
            normalized=True)
        expected = 0.0
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            [1, 1, 1, 0, 0], [0, 0, 0, 0, 0], self.oids1, self.t1,
            normalized=True)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            [], [], [], self.t1, normalized=True)
        expected = 0.0
        self.assertAlmostEqual(actual, expected)

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

    def test_weighted_unifrac_branch_correction(self):
        # for ((a:1, b:2)c:3,(d:4,e:5)f:6)root;"
        tip_ds = np.array([4, 5, 10, 11, 0, 0, 0])[:, np.newaxis]
        u_counts = np.array([1, 1, 0, 0, 2, 0, 2])
        v_counts = np.array([0, 2, 1, 0, 2, 1, 3])
        u_sum = 2  # counts at the tips
        v_sum = 3
        exp = np.array([2.0,
                        5.0 * (.5 + (2.0/3.0)),
                        10.0 * (1.0 / 3.0),
                        0.0]).sum()
        obs = _weighted_unifrac_branch_correction(
            tip_ds, u_counts/u_sum, v_counts/v_sum)
        self.assertEqual(obs, exp)

    def test_unweighted_unifrac_pycogent_adapted(self):
        # adapted from PyCogent unit tests
        m = np.array([[1, 0, 1], [1, 1, 0], [0, 1, 0], [0, 0, 1], [0, 1, 0],
                      [0, 1, 1], [1, 1, 1], [0, 1, 1], [1, 1, 1]])
        # lengths from ((a:1,b:2):4,(c:3,(d:1,e:1):2):3)
        bl = np.array([1, 2, 1, 1, 3, 2, 4, 3, 0], dtype=float)
        self.assertEqual(_unweighted_unifrac(m[:, 0], m[:, 1], bl), 10/16.0)
        self.assertEqual(_unweighted_unifrac(m[:, 0], m[:, 2], bl), 8/13.0)
        self.assertEqual(_unweighted_unifrac(m[:, 1], m[:, 2], bl), 8/17.0)

    def test_weighted_unifrac_pycogent_adapted(self):
        # lengths from ((a:1,b:2):4,(c:3,(d:1,e:1):2):3)
        bl = np.array([1, 2, 1, 1, 3, 2, 4, 3, 0], dtype=float)

        # adapted from PyCogent unit tests
        m = np.array([[1, 0, 1],  # a
                      [1, 1, 0],  # b
                      [0, 1, 0],  # d
                      [0, 0, 1],  # e
                      [0, 1, 0],  # c
                      [0, 1, 1],  # parent of (d, e)
                      [2, 1, 1],  # parent of a, b
                      [0, 2, 1],  # parent of c (d, e)
                      [2, 3, 2]])  # root

        # sum just the counts at the tips
        m0s = m[:5, 0].sum()
        m1s = m[:5, 1].sum()
        m2s = m[:5, 2].sum()

        # scores computed by educational implementation
        self.assertAlmostEqual(
            _weighted_unifrac(m[:, 0], m[:, 1], m0s, m1s, bl)[0], 7.5)
        self.assertAlmostEqual(
            _weighted_unifrac(m[:, 0], m[:, 2], m0s, m2s, bl)[0], 6.0)
        self.assertAlmostEqual(
            _weighted_unifrac(m[:, 1], m[:, 2], m1s, m2s, bl)[0], 4.5)


if __name__ == '__main__':
    main()
