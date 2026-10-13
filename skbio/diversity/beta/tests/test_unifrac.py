# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import warnings
from io import StringIO
from unittest import main, TestCase, skipIf
from unittest.mock import patch

import numpy as np

from skbio import TreeNode, DistanceMatrix
from skbio.tree import DuplicateNodeError, MissingNodeError
from skbio.diversity import beta_diversity
from skbio.diversity.beta import (unweighted_unifrac, weighted_unifrac,
                                  generalized_unifrac)
from skbio.diversity.beta._unifrac import (_unweighted_unifrac,
                                           _weighted_unifrac,
                                           _weighted_unifrac_branch_correction,
                                           _unweighted_unifrac_pdist_numba,
                                           _weighted_unifrac_pdist_numba,
                                           _generalized_unifrac,
                                           _generalized_unifrac_pdist_numba,
                                           _setup_pairwise_unifrac,
                                           NUMBA_AVAILABLE)
from skbio.diversity._driver import _UNIFRAC_FAST_ENGINE
from skbio.diversity.beta.tests._fixtures import (QiimeTinyTestMixin,
                                                  patch_gpu_backend)
from skbio.util import numba_code

# Measured max abs deviation, CPU-numba vs ssu-ascii-fixture, across all 8
# qiime-191-tt generalized/variance_adjust fixtures and every sample pair:
# 1.332e-7 (generalized_unifrac, alpha=0.5, f2 vs p2). GPU-vs-CPU-numba
# deviation is far tighter (3.331e-16, measured on real NVIDIA/AMD hardware).
# See docs/superpowers/specs/2026-10-02-ssu-unifrac-numba-phase1-design.md.
# 1.5e-6 gives an ~11x safety margin over the larger (ssu-fixture) figure.
SSU_FIXTURE_TOLERANCE = 1.5e-6


class UnifracTests(QiimeTinyTestMixin, TestCase):

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

    def test_unweighted_taxa_out_of_order(self):
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

    def test_weighted_taxa_out_of_order(self):
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
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(DuplicateNodeError, unweighted_unifrac,
                          u_counts, v_counts, taxa, t)
        self.assertRaises(DuplicateNodeError, weighted_unifrac, u_counts,
                          v_counts, taxa, t)

        # unrooted tree as input
        t = TreeNode.read(StringIO('((OTU1:0.1, OTU2:0.2):0.3, OTU3:0.5,'
                                   'OTU4:0.7);'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, taxa, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, taxa, t)

        # taxa has duplicated ids
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        taxa = ['OTU1', 'OTU2', 'OTU2']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, taxa, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, taxa, t)

        # len of vectors not equal
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2]
        v_counts = [1, 1, 1]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, taxa, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, taxa, t)
        u_counts = [1, 2, 3]
        v_counts = [1, 1]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, taxa, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, taxa, t)
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        taxa = ['OTU1', 'OTU2']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, taxa, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, taxa, t)

        # negative counts
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, -3]
        v_counts = [1, 1, 1]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, taxa, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, taxa, t)
        u_counts = [1, 2, 3]
        v_counts = [1, 1, -1]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, taxa, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, taxa, t)

        # tree with no branch lengths
        t = TreeNode.read(
            StringIO('((((OTU1,OTU2),OTU3)),(OTU4,OTU5));'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, taxa, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, taxa, t)

        # tree missing some branch lengths
        t = TreeNode.read(
            StringIO('(((((OTU1,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, taxa, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, taxa, t)

        # taxa not present in tree
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        taxa = ['OTU1', 'OTU2', 'OTU42']
        self.assertRaises(MissingNodeError, unweighted_unifrac, u_counts,
                          v_counts, taxa, t)
        self.assertRaises(MissingNodeError, weighted_unifrac, u_counts,
                          v_counts, taxa, t)

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
            self.b1[0], self.b1[1], self.oids1, self.t1, validate=False)
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

    # Fixture-matched coverage for normalized=False and variance_adjust=True
    # on unweighted_unifrac now lives in test_unifrac_gpu.py: those options
    # are GPU-only in this release, with CPU/numba support planned for a
    # future PR (see the NotImplementedError tests below).

    def test_unweighted_unifrac_normalized_false_requires_gpu(self):
        with self.assertRaises(NotImplementedError):
            unweighted_unifrac(
                [1, 0, 1], [0, 1, 1], ['a', 'b', 'c'], self.t1,
                normalized=False,
            )

    def test_unweighted_unifrac_variance_adjust_requires_gpu(self):
        with self.assertRaises(NotImplementedError):
            unweighted_unifrac(
                [1, 0, 1], [0, 1, 1], ['a', 'b', 'c'], self.t1,
                variance_adjust=True,
            )

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

    # Fixture-matched coverage for weighted_unifrac's variance_adjust now
    # lives in test_unifrac_gpu.py: that option is GPU-only in this release,
    # with CPU/numba support planned for a future PR (see the
    # NotImplementedError test below).

    def test_weighted_unifrac_variance_adjust_requires_gpu(self):
        with self.assertRaises(NotImplementedError):
            weighted_unifrac(
                self.b1[0], self.b1[1], self.oids1, self.t1,
                variance_adjust=True,
            )

    # generalized_unifrac ships GPU-only in this release; CPU/numba support
    # is planned for a future PR. The CPU kernel (_generalized_unifrac /
    # _generalized_unifrac_pdist_numba) stays in skbio.diversity.beta._unifrac
    # as the reference that future PR can build from, so the fixture- and
    # kernel-correctness checks below exercise it directly rather than
    # through the public generalized_unifrac() function, which now always
    # requires engine='gpu'.

    def _assert_pairs_match_ssu_fixture(self, condensed, sample_ids, fname):
        """Spot-check a few sample pairs against a qiime-191-tt fixture matrix.

        A handful of representative pairs rather than the whole matrix: these
        checks pin the CPU kernel against the ssu reference, and the GPU and
        array-API paths are compared against this kernel entry-by-entry in
        test_unifrac_gpu.py/test_unifrac_xp.py.
        """
        obs = DistanceMatrix(condensed, sample_ids)
        expected = self._load_dm_fixture(fname)
        for i, j in [(0, 1), (2, 5), (3, 7)]:
            self.assertAlmostEqual(
                obs[sample_ids[i], sample_ids[j]],
                expected[sample_ids[i], sample_ids[j]],
                delta=SSU_FIXTURE_TOLERANCE
            )

    @numba_code
    def test_generalized_unifrac_matches_ssu_fixture(self):
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        for alpha, fname in [(0.5, 'generalized_unifrac_alpha0.5_dm.txt'),
                             (1.0, 'generalized_unifrac_alpha1.0_dm.txt')]:
            condensed = _generalized_unifrac_pdist_numba(
                table, taxa, tree, alpha=alpha, validate=True,
            )
            self._assert_pairs_match_ssu_fixture(condensed, sample_ids, fname)

    @numba_code
    def test_generalized_unifrac_variance_adjust_matches_ssu_fixture(self):
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        condensed = _generalized_unifrac_pdist_numba(
            table, taxa, tree, alpha=1.0, variance_adjust=True, validate=True,
        )
        self._assert_pairs_match_ssu_fixture(
            condensed, sample_ids, 'generalized_unifrac_alpha1.0_vaw_dm.txt')

    @numba_code
    def test_generalized_unifrac_both_empty_is_zero(self):
        u_node_counts, v_node_counts, u_total_count, v_total_count, tree_index = (
            _setup_pairwise_unifrac(
                [0, 0, 0], [0, 0, 0], self.oids1[:3], self.t1,
                True, normalized=True, unweighted=False,
            )
        )
        obs = _generalized_unifrac(
            u_node_counts, v_node_counts, u_total_count, v_total_count,
            tree_index["length"], 1.0, False,
        )
        self.assertEqual(obs, 0.0)

    def test_generalized_unifrac_both_empty_is_zero_via_real_dispatch(self):
        # Companion to test_generalized_unifrac_both_empty_is_zero above:
        # that test only pins the dead/unreachable CPU reference kernel
        # (_generalized_unifrac) directly. This exercises the same
        # both-samples-empty case through the live dispatch path
        # (_generalized_unifrac_gpu_or_xp) that real callers of
        # generalized_unifrac(engine='gpu')/beta_diversity actually hit --
        # the fused kernel when a GPU backend is usable, otherwise the
        # array-API fallback, either way without needing real GPU hardware
        # to run this test.
        from skbio.diversity.beta._unifrac_gpu import _generalized_unifrac_gpu_or_xp
        table = np.array([[0, 0, 0], [0, 0, 0]])
        obs = _generalized_unifrac_gpu_or_xp(
            table, self.oids1[:3], self.t1, alpha=1.0, variance_adjust=False,
        )
        np.testing.assert_allclose(obs, [0.0])

    def test_generalized_unifrac_cpu_raises_not_implemented(self):
        # generalized_unifrac has no CPU implementation in this release; it
        # ships GPU-only regardless of numba availability (CPU/numba support
        # is planned for a future PR).
        with self.assertRaises(NotImplementedError):
            generalized_unifrac([1, 0, 1], [0, 1, 1], ['a', 'b', 'c'], self.t1)

    def test_generalized_unifrac_cpu_raises_before_alpha_check(self):
        # The "requires GPU" check runs before alpha-range validation on the
        # default (non-'gpu') engine path, mirroring the numba-before-alpha
        # ordering this replaced (see test_generalized_unifrac_alpha_out_of_
        # range_raises below for the GPU-available alpha-validation path).
        with self.assertRaises(NotImplementedError):
            generalized_unifrac(
                [1, 0, 1], [0, 1, 1], ['a', 'b', 'c'], self.t1, alpha=1.5,
            )

    def test_generalized_unifrac_alpha_out_of_range_raises(self):
        with patch_gpu_backend("cuda"):
            with self.assertRaises(ValueError):
                generalized_unifrac(
                    [1, 0, 1], [0, 1, 1], ['a', 'b', 'c'], self.t1,
                    alpha=1.5, engine='gpu',
                )

    def test_beta_diversity_generalized_unifrac_cpu_raises_before_alpha_check(self):
        # beta_diversity counterpart of the regression test above.
        with self.assertRaises(NotImplementedError):
            beta_diversity(
                "generalized_unifrac", self.b1, ids=self.sids1,
                taxa=self.oids1, tree=self.t1, alpha=1.5,
            )

    @numba_code
    def test_generalized_unifrac_pdist_numba_matches_single_pair(self):
        # generalized_unifrac has no cython/pairwise_func counterpart to
        # compare the numba pdist kernel against, so cross-check
        # _generalized_unifrac_pdist_numba directly against the single-pair
        # _generalized_unifrac path (already validated against the SSU
        # fixtures above) for every pair in the tiny-test table, across both
        # alpha values and with/without variance_adjust.
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        n = table.shape[0]
        for alpha in (0.5, 1.0):
            for variance_adjust in (False, True):
                condensed = _generalized_unifrac_pdist_numba(
                    table, taxa, tree, alpha=alpha,
                    variance_adjust=variance_adjust, validate=True,
                )
                obs = DistanceMatrix(condensed, sample_ids)
                for i in range(n):
                    for j in range(i + 1, n):
                        (
                            u_node_counts, v_node_counts,
                            u_total_count, v_total_count, tree_index,
                        ) = _setup_pairwise_unifrac(
                            table[i], table[j], taxa, tree,
                            True, normalized=True, unweighted=False,
                        )
                        expected = _generalized_unifrac(
                            u_node_counts, v_node_counts,
                            u_total_count, v_total_count,
                            tree_index["length"], alpha, variance_adjust,
                        )
                        self.assertAlmostEqual(
                            obs[sample_ids[i], sample_ids[j]], expected, places=10,
                        )

    @numba_code
    def test_unweighted_unifrac_engine_numba_matches_cython(self):
        dm_cy = beta_diversity(
            "unweighted_unifrac", self.b1, ids=self.sids1, taxa=self.oids1,
            tree=self.t1, engine="cython")
        dm_nb = beta_diversity(
            "unweighted_unifrac", self.b1, ids=self.sids1, taxa=self.oids1,
            tree=self.t1, engine="numba")

        self.assertIsInstance(dm_nb, DistanceMatrix)
        self.assertEqual(list(dm_nb.ids), list(dm_cy.ids))
        np.testing.assert_allclose(
            dm_nb.data, dm_cy.data, rtol=1e-12, atol=1e-12)

    @numba_code
    def test_unifrac_fast_resolves_numba(self):
        # Comparing results cannot show which engine "fast" picked: the two
        # unifrac engines agree well within any tolerance a test could use.
        # The choice is a named constant, so assert that instead.
        self.assertEqual(_UNIFRAC_FAST_ENGINE, "numba")

    @skipIf(NUMBA_AVAILABLE, "covers the branch taken when numba is absent")
    def test_unifrac_fast_resolves_cython(self):
        # The counterpart to the test above, for the lane with no numba.
        self.assertEqual(_UNIFRAC_FAST_ENGINE, "cython")

    @numba_code
    def test_unifrac_engine_fast_is_accepted(self):
        # Checks that "fast" is plumbed through and gives the same answer,
        # not which engine ran. The numba and cython unifrac kernels agree to
        # well inside this tolerance, which the two
        # engine_numba_matches_cython tests in this file already show, so no
        # comparison of their results can tell the two engines apart. Which
        # engine "fast" resolves to is covered in skbio/tests/test_config.py.
        for metric in ("unweighted_unifrac", "weighted_unifrac"):
            dm_fast = beta_diversity(
                metric, self.b1, ids=self.sids1, taxa=self.oids1,
                tree=self.t1, engine="fast")
            dm_nb = beta_diversity(
                metric, self.b1, ids=self.sids1, taxa=self.oids1,
                tree=self.t1, engine="numba")
            np.testing.assert_allclose(dm_fast.data, dm_nb.data,
                                       rtol=1e-12, atol=1e-12)

    @skipIf(NUMBA_AVAILABLE, "covers the branch taken when numba is absent")
    def test_unifrac_engine_fast_is_cython_without_numba(self):
        # The counterpart to test_unifrac_engine_fast_is_accepted above.
        # Without numba installed, "fast" resolves to "cython" and runs the
        # identical cython call, so unlike that one this is an exact
        # comparison, and it does pin which engine ran.
        for metric in ("unweighted_unifrac", "weighted_unifrac"):
            dm_fast = beta_diversity(
                metric, self.b1, ids=self.sids1, taxa=self.oids1,
                tree=self.t1, engine="fast")
            dm_cy = beta_diversity(
                metric, self.b1, ids=self.sids1, taxa=self.oids1,
                tree=self.t1, engine="cython")
            np.testing.assert_array_equal(dm_fast.data, dm_cy.data)

    @numba_code
    def test_unweighted_unifrac_engine_numba_larger_random(self):
        rng = np.random.default_rng(0)
        counts = rng.integers(0, 10, size=(12, 5))
        # cover the "both samples empty" (observed == 0) branch: rows 0 and 1
        # are both all-zero, guaranteed rather than left to chance
        counts[0] = 0
        counts[1] = 0
        ids = [f"S{i}" for i in range(counts.shape[0])]

        dm_cy = beta_diversity(
            "unweighted_unifrac", counts, ids=ids, taxa=self.oids1,
            tree=self.t1, engine="cython")
        dm_nb = beta_diversity(
            "unweighted_unifrac", counts, ids=ids, taxa=self.oids1,
            tree=self.t1, engine="numba")

        np.testing.assert_allclose(
            dm_nb.data, dm_cy.data, rtol=1e-12, atol=1e-12)
        # rows 0 and 1 are both all-zero, so their pairwise distance must be 0
        self.assertEqual(dm_nb.data[0, 1], 0.0)

    @numba_code
    def test_unweighted_unifrac_engine_numba_single_sample(self):
        dm_nb = beta_diversity(
            "unweighted_unifrac", self.b1[:1], ids=self.sids1[:1],
            taxa=self.oids1, tree=self.t1, engine="numba")

        self.assertEqual(dm_nb.shape, (1, 1))
        self.assertEqual(dm_nb.data[0, 0], 0.0)

    @numba_code
    def test_weighted_unifrac_engine_numba_matches_cython(self):
        dm_cy = beta_diversity(
            "weighted_unifrac", self.b1, ids=self.sids1, taxa=self.oids1,
            tree=self.t1, engine="cython")
        dm_nb = beta_diversity(
            "weighted_unifrac", self.b1, ids=self.sids1, taxa=self.oids1,
            tree=self.t1, engine="numba")

        self.assertIsInstance(dm_nb, DistanceMatrix)
        self.assertEqual(list(dm_nb.ids), list(dm_cy.ids))
        # weighted UniFrac involves a division and a sum whose ordering
        # differs between NumPy's pairwise summation (Cython/pdist path)
        # and the kernel's sequential accumulation, so use a looser
        # tolerance than the unweighted comparison.
        np.testing.assert_allclose(
            dm_nb.data, dm_cy.data, rtol=1e-10, atol=1e-12)

    @numba_code
    def test_weighted_unifrac_engine_numba_normalized_matches_cython(self):
        dm_cy = beta_diversity(
            "weighted_unifrac", self.b1, ids=self.sids1, taxa=self.oids1,
            tree=self.t1, engine="cython", normalized=True)
        dm_nb = beta_diversity(
            "weighted_unifrac", self.b1, ids=self.sids1, taxa=self.oids1,
            tree=self.t1, engine="numba", normalized=True)

        self.assertIsInstance(dm_nb, DistanceMatrix)
        self.assertEqual(list(dm_nb.ids), list(dm_cy.ids))
        np.testing.assert_allclose(
            dm_nb.data, dm_cy.data, rtol=1e-10, atol=1e-12)
        self.assertTrue((dm_nb.data <= 1.0 + 1e-9).all())

    @numba_code
    def test_weighted_unifrac_engine_numba_larger_random(self):
        rng = np.random.default_rng(0)
        counts = rng.integers(0, 10, size=(12, 5))
        # cover both the unnormalized zero-total path and the normalized
        # both-empty 0/0 guard
        counts[0] = 0
        counts[1] = 0
        ids = [f"S{i}" for i in range(counts.shape[0])]

        for normalized in (False, True):
            dm_cy = beta_diversity(
                "weighted_unifrac", counts, ids=ids, taxa=self.oids1,
                tree=self.t1, engine="cython", normalized=normalized)
            dm_nb = beta_diversity(
                "weighted_unifrac", counts, ids=ids, taxa=self.oids1,
                tree=self.t1, engine="numba", normalized=normalized)

            np.testing.assert_allclose(
                dm_nb.data, dm_cy.data, rtol=1e-10, atol=1e-12)

        # both rows 0 and 1 are all-zero -> normalized 0/0 guard -> 0.0
        self.assertEqual(dm_nb.data[0, 1], 0.0)

    @numba_code
    def test_weighted_unifrac_engine_numba_single_sample(self):
        dm_nb = beta_diversity(
            "weighted_unifrac", self.b1[:1], ids=self.sids1[:1],
            taxa=self.oids1, tree=self.t1, engine="numba")

        self.assertEqual(dm_nb.shape, (1, 1))
        self.assertEqual(dm_nb.data[0, 0], 0.0)

    def test_beta_diversity_generalized_unifrac_cpu_raises_not_implemented(self):
        # generalized_unifrac ships GPU-only in this release; see
        # test_beta_diversity_generalized_unifrac_gpu_engine_raises_without_gpu
        # below for the (hardware-gated) engine='gpu' coverage.
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        with self.assertRaises(NotImplementedError):
            beta_diversity(
                "generalized_unifrac", table, ids=sample_ids, taxa=taxa,
                tree=tree, alpha=0.5,
            )

    def test_beta_diversity_unweighted_unifrac_variance_adjust_cpu_raises_not_implemented(
        self,
    ):
        # variance_adjust is GPU-only for unweighted_unifrac in this release.
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        with self.assertRaises(NotImplementedError):
            beta_diversity(
                "unweighted_unifrac", table, ids=sample_ids, taxa=taxa,
                tree=tree, variance_adjust=True,
            )

    def test_beta_diversity_weighted_unifrac_variance_adjust_cpu_raises_not_implemented(
        self,
    ):
        # variance_adjust is GPU-only for weighted_unifrac in this release.
        table, taxa, tree, sample_ids = self._load_qiime_191_tt()
        for normalized in (False, True):
            with self.assertRaises(NotImplementedError):
                beta_diversity(
                    "weighted_unifrac", table, ids=sample_ids, taxa=taxa,
                    tree=tree, normalized=normalized, variance_adjust=True,
                )

    def test_beta_diversity_engine_invalid(self):
        with self.assertRaisesRegex(
                ValueError, "engine='julia' is not supported"):
            beta_diversity(
                "unweighted_unifrac", self.b1, ids=self.sids1,
                taxa=self.oids1, tree=self.t1, engine="julia")

    def test_beta_diversity_generalized_unifrac_engine_invalid(self):
        # generalized_unifrac has no cython/numba choice to resolve (it has
        # no cython path at all), so unlike unweighted_unifrac/
        # weighted_unifrac above, this branch never called anything that
        # validated `engine`; a typo such as 'bogus' or 'cuda' used to fall
        # through and silently run the CPU numba path.
        for bad in ('bogus', 'cuda', 'GPU', 'numba', 'cython'):
            with self.assertRaisesRegex(ValueError, 'engine'):
                beta_diversity(
                    "generalized_unifrac", self.b1, ids=self.sids1,
                    taxa=self.oids1, tree=self.t1, engine=bad)

    @numba_code
    def test_beta_diversity_engine_numba_with_pairwise_func_is_used(self):
        calls = []

        def recording_pdist(counts, metric, **kwargs):
            calls.append(metric)
            from scipy.spatial.distance import pdist
            return pdist(counts, metric=metric, **kwargs)

        beta_diversity(
            "unweighted_unifrac", self.b1, ids=self.sids1, taxa=self.oids1,
            tree=self.t1, engine="numba", pairwise_func=recording_pdist)

        self.assertEqual(len(calls), 1)

    @numba_code
    def test_beta_diversity_engine_fast_with_pairwise_func_is_quiet(self):
        # engine="numba" warns when a pairwise_func stops the numba kernels
        # being used, because the caller asked for numba and did not get it.
        # "fast" asks scikit-bio to choose, so nothing the caller asked for is
        # overridden and there is nothing to warn about; it just uses the
        # pairwise_func, the same as engine="cython" and the default do.
        #
        # There is no no-numba counterpart to this one, unlike the two tests
        # above. Without numba "fast" resolves to cython, which cannot reach
        # the warning at all, so the test would pass whatever this code did.
        calls = []

        def recording_pdist(counts, metric, **kwargs):
            calls.append(metric)
            from scipy.spatial.distance import pdist
            return pdist(counts, metric=metric, **kwargs)

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            beta_diversity(
                "unweighted_unifrac", self.b1, ids=self.sids1,
                taxa=self.oids1, tree=self.t1, engine="fast",
                pairwise_func=recording_pdist)

        self.assertEqual(len(calls), 1)
        # Look for the engine warning specifically rather than asserting that
        # nothing at all was raised, so an unrelated warning from a dependency
        # cannot fail this on some other environment.
        self.assertEqual(
            [str(w.message) for w in caught if "engine=" in str(w.message)], [])

    @numba_code
    def test_beta_diversity_engine_numba_bogus_kwarg_raises(self):
        for engine in ("cython", "numba"):
            with self.assertRaises(TypeError):
                beta_diversity(
                    "unweighted_unifrac", self.b1, ids=self.sids1,
                    taxa=self.oids1, tree=self.t1, engine=engine,
                    bogus_kwarg=True)

    # -- keyword-only signatures ---------------------------------------

    def test_unifrac_options_are_keyword_only(self):
        # Before `normalized`/`variance_adjust`/`alpha`/`engine` were added,
        # `validate` was the last positional-or-keyword parameter, so
        # `unweighted_unifrac(u, v, taxa, tree, False)` meant validate=False.
        # Those options are keyword-only so such a call raises instead of
        # silently binding False to a different parameter.
        for func in (unweighted_unifrac, weighted_unifrac, generalized_unifrac):
            with self.assertRaises(TypeError):
                func(self.b1[0], self.b1[1], self.oids1, self.t1, False)

    @numba_code
    def test_unifrac_pdist_numba_options_are_keyword_only(self):
        # The three pdist kernels took their options in three different
        # orders, so positional calls were easy to get wrong; everything
        # after `tree` is keyword-only now.
        for func in (_unweighted_unifrac_pdist_numba,
                     _weighted_unifrac_pdist_numba,
                     _generalized_unifrac_pdist_numba):
            with self.assertRaises(TypeError):
                func(self.b1, self.oids1, self.t1, True)

    def test_unifrac_invalid_engine_raises(self):
        for func in (unweighted_unifrac, weighted_unifrac, generalized_unifrac):
            for bad in ('cuda', 'GPU', 'numba', 'cython'):
                with self.assertRaisesRegex(ValueError, 'engine'):
                    func(self.b1[0], self.b1[1], self.oids1, self.t1,
                         engine=bad)

    # -- beta_diversity kwarg handling ---------------------------------

    def test_beta_diversity_unweighted_unifrac_unnormalized_cpu_raises_not_implemented(
        self,
    ):
        # normalized=False for unweighted_unifrac is GPU-only in this
        # release, for every CPU-side engine choice.
        for engine in (None, 'cython', 'numba', 'fast'):
            with self.assertRaises(NotImplementedError):
                beta_diversity(
                    "unweighted_unifrac", self.b1, ids=self.sids1,
                    taxa=self.oids1, tree=self.t1, normalized=False,
                    engine=engine)

    def test_beta_diversity_variance_adjust_with_pairwise_func_raises(self):
        # variance_adjust has no CPU implementation at all in this release
        # (GPU-only), so it must raise regardless of pairwise_func.
        def recording_pdist(counts, metric, **kwargs):
            from scipy.spatial.distance import pdist
            return pdist(counts, metric=metric, **kwargs)

        for metric in ("unweighted_unifrac", "weighted_unifrac"):
            with self.assertRaisesRegex(NotImplementedError, "variance_adjust"):
                beta_diversity(
                    metric, self.b1, ids=self.sids1, taxa=self.oids1,
                    tree=self.t1, variance_adjust=True,
                    pairwise_func=recording_pdist)

    def test_beta_diversity_gpu_engine_bogus_kwarg_raises(self):
        # Checked before GPU backend detection, so this runs with or without
        # a GPU present.
        for metric in ("unweighted_unifrac", "weighted_unifrac",
                       "generalized_unifrac"):
            with self.assertRaises(TypeError):
                beta_diversity(
                    metric, self.b1, ids=self.sids1, taxa=self.oids1,
                    tree=self.t1, engine="gpu", bogus_kwarg=True)

    def test_beta_diversity_gpu_engine_pairwise_func_raises(self):
        def not_a_real_pdist(counts, metric, **kwargs):
            return np.zeros(counts.shape[0] * (counts.shape[0] - 1) // 2)

        for metric in ("unweighted_unifrac", "weighted_unifrac",
                       "generalized_unifrac"):
            with self.assertRaisesRegex(ValueError, "pairwise_func"):
                beta_diversity(
                    metric, self.b1, ids=self.sids1, taxa=self.oids1,
                    tree=self.t1, engine="gpu",
                    pairwise_func=not_a_real_pdist)

    # -- engine='gpu' wiring -------------------------------------------
    #
    # engine='gpu' now always produces a result: a fused kernel when a GPU
    # backend is usable, or the array-API fallback (_unifrac_xp.py)
    # otherwise. These tests force the no-backend case by monkeypatching
    # _detect_gpu_backend(), deterministically exercising the fallback path
    # regardless of what hardware this machine actually has, and check the
    # result against the CPU reference. GPU-available fused-kernel
    # correctness coverage for the underlying *_gpu drivers lives in
    # test_unifrac_gpu.py.

    def test_unweighted_unifrac_gpu_engine_falls_back_without_gpu(self):
        with patch_gpu_backend(None):
            obs = unweighted_unifrac(
                self.b1[0], self.b1[1], self.oids1, self.t1, engine='gpu')
        expected = unweighted_unifrac(
            self.b1[0], self.b1[1], self.oids1, self.t1)
        self.assertAlmostEqual(obs, expected)

    def test_weighted_unifrac_gpu_engine_falls_back_without_gpu(self):
        with patch_gpu_backend(None):
            obs = weighted_unifrac(
                self.b1[0], self.b1[1], self.oids1, self.t1, engine='gpu')
        expected = weighted_unifrac(
            self.b1[0], self.b1[1], self.oids1, self.t1)
        self.assertAlmostEqual(obs, expected)

    def test_generalized_unifrac_gpu_engine_falls_back_without_gpu(self):
        with patch_gpu_backend(None):
            obs = generalized_unifrac(
                self.b1[0], self.b1[1], self.oids1, self.t1, engine='gpu')
        u_node_counts, v_node_counts, u_total_count, v_total_count, tree_index = (
            _setup_pairwise_unifrac(
                self.b1[0], self.b1[1], self.oids1, self.t1,
                True, normalized=True, unweighted=False,
            )
        )
        expected = _generalized_unifrac(
            u_node_counts, v_node_counts, u_total_count, v_total_count,
            tree_index["length"], 1.0, False,
        )
        self.assertAlmostEqual(obs, expected)

    def test_beta_diversity_unweighted_unifrac_gpu_engine_falls_back_without_gpu(
        self,
    ):
        with patch_gpu_backend(None):
            obs = beta_diversity(
                "unweighted_unifrac", self.b1, ids=self.sids1,
                taxa=self.oids1, tree=self.t1, engine="gpu")
        expected = beta_diversity(
            "unweighted_unifrac", self.b1, ids=self.sids1,
            taxa=self.oids1, tree=self.t1)
        np.testing.assert_allclose(obs.data, expected.data, atol=1e-10)

    def test_beta_diversity_weighted_unifrac_gpu_engine_falls_back_without_gpu(
        self,
    ):
        with patch_gpu_backend(None):
            obs = beta_diversity(
                "weighted_unifrac", self.b1, ids=self.sids1,
                taxa=self.oids1, tree=self.t1, engine="gpu")
        expected = beta_diversity(
            "weighted_unifrac", self.b1, ids=self.sids1,
            taxa=self.oids1, tree=self.t1)
        np.testing.assert_allclose(obs.data, expected.data, atol=1e-10)

    def test_beta_diversity_generalized_unifrac_gpu_engine_falls_back_without_gpu(
        self,
    ):
        with patch_gpu_backend(None):
            obs = beta_diversity(
                "generalized_unifrac", self.b1, ids=self.sids1,
                taxa=self.oids1, tree=self.t1, engine="gpu", alpha=0.5)
        expected_condensed = _generalized_unifrac_pdist_numba(
            self.b1, taxa=self.oids1, tree=self.t1, alpha=0.5)
        expected = DistanceMatrix(expected_condensed, self.sids1)
        np.testing.assert_allclose(obs.data, expected.data, atol=1e-10)

    def test_weighted_unifrac_gpu_engine_falls_back_on_kernel_failure(self):
        # Even when a GPU backend IS detected, a fused kernel that fails to
        # build/run must fall back to the array-API path rather than
        # propagating the exception (mirrors the PERMANOVA/Mantel
        # `_mark_gpu_unavailable` precedent).
        from skbio.diversity.beta import _unifrac_gpu

        self.addCleanup(_unifrac_gpu._unavailable_backends.discard, "cuda")
        with patch_gpu_backend("cuda"), patch(
            "skbio.diversity.beta._unifrac_gpu._weighted_unifrac_gpu",
            side_effect=RuntimeError("kernel build failed"),
        ):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                obs = weighted_unifrac(
                    self.b1[0], self.b1[1], self.oids1, self.t1, engine='gpu')
        expected = weighted_unifrac(
            self.b1[0], self.b1[1], self.oids1, self.t1)
        self.assertAlmostEqual(obs, expected)
        self.assertIn("cuda", _unifrac_gpu._unavailable_backends)


if __name__ == '__main__':
    main()
