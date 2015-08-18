# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main
from io import StringIO
import os

import numpy as np
import pandas as pd

from skbio import TreeNode, DistanceMatrix
from skbio.util import get_data_path
from skbio.tree import DuplicateNodeError, MissingNodeError
from skbio.diversity.beta import unweighted_unifrac, weighted_unifrac


class TestUniFrac(TestCase):

    def setUp(self):
        self.table1 = np.array(
           [[1, 3, 0, 1, 0],
            [0, 2, 0, 4, 4],
            [0, 0, 6, 2, 1],
            [0, 0, 1, 1, 1],
            [5, 3, 5, 0, 0],
            [0, 0, 0, 3, 5]])
        self.sids1 = list('ABCDEF')
        self.oids1 = ['OTU%d' % i for i in range(1, 6)]
        self.t1 = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        self.t1_w_extra_tips = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,(OTU5:0.25,(OTU6:0.5,OTU7:0.5):0.5):0.5):1.25):0.0'
                     u')root;'))

        self.t2 = TreeNode.read(
            StringIO(u'((OTU1:0.1, OTU2:0.2):0.3, (OTU3:0.5, OTU4:0.7):1.1)'
                     u'root;'))
        self.oids2 = ['OTU%d' % i for i in range(1, 5)]

        # the following table and tree are derived from the QIIME 1.9.1
        # "tiny-test" data
        tt_table_fp = get_data_path(
            os.path.join('qiime-191-tt', 'otu-table.tsv'), 'data')
        tt_tree_fp = get_data_path(
            os.path.join('qiime-191-tt', 'tree.nwk'), 'data')

        self.q_table = pd.read_csv(tt_table_fp, sep='\t', skiprows=1,
                                   index_col=0)
        self.q_tree = TreeNode.read(tt_tree_fp)

    def test_unweighted_unifrac_qiime_tiny_test(self):
        dm_fp = get_data_path(
            os.path.join('qiime-191-tt', 'unweighted_unifrac_dm.txt'), 'data')
        expected = DistanceMatrix.read(dm_fp)
        for sid1 in self.q_table.columns:
            for sid2 in self.q_table.columns:
                actual = unweighted_unifrac(
                    self.q_table[sid1], self.q_table[sid2],
                    otu_ids=self.q_table.index, tree=self.q_tree)
                self.assertAlmostEqual(actual, expected[sid1, sid2])

    def test_weighted_unifrac_qiime_tiny_test(self):
        dm_fp = get_data_path(
            os.path.join('qiime-191-tt', 'weighted_unifrac_dm.txt'), 'data')
        expected = DistanceMatrix.read(dm_fp)
        for sid1 in self.q_table.columns:
            for sid2 in self.q_table.columns:
                actual = weighted_unifrac(
                    self.q_table[sid1], self.q_table[sid2],
                    otu_ids=self.q_table.index, tree=self.q_tree)
                self.assertAlmostEqual(actual, expected[sid1, sid2],
                                       msg="%s, %s" % (sid1, sid2))

    def test_weighted_normalized_unifrac_qiime_tiny_test(self):
        dm_fp = get_data_path(
            os.path.join('qiime-191-tt', 'weighted_normalized_unifrac_dm.txt'),
            'data')
        expected = DistanceMatrix.read(dm_fp)
        for sid1 in self.q_table.columns:
            for sid2 in self.q_table.columns:
                actual = weighted_unifrac(
                    self.q_table[sid1], self.q_table[sid2],
                    otu_ids=self.q_table.index, tree=self.q_tree,
                    normalized=True)
                self.assertAlmostEqual(actual, expected[sid1, sid2])

    def test_unweighted_extra_tips(self):
        # UniFrac values are the same despite unobserved tips in the tree
        for i in range(len(self.table1)):
            for j in range(len(self.table1)):
                actual = unweighted_unifrac(
                    self.table1[i], self.table1[j], self.oids1,
                    self.t1_w_extra_tips)
                expected = unweighted_unifrac(
                    self.table1[i], self.table1[j], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_weighted_extra_tips(self):
        # UniFrac values are the same despite unobserved tips in the tree
        for i in range(len(self.table1)):
            for j in range(len(self.table1)):
                actual = weighted_unifrac(
                    self.table1[i], self.table1[j], self.oids1,
                    self.t1_w_extra_tips)
                expected = weighted_unifrac(
                    self.table1[i], self.table1[j], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_unweighted_minimal_trees(self):
        # expected values computed by hand
        # zero tips
        tree = TreeNode.read(StringIO(u'root;'))
        actual = unweighted_unifrac([], [], [], tree)
        expected = 0.0
        self.assertEqual(actual, expected)

        # two tips
        tree = TreeNode.read(StringIO(u'(OTU1:0.25, OTU2:0.25)root;'))
        actual = unweighted_unifrac([1, 0], [0, 0], ['OTU1', 'OTU2'], tree)
        expected = 1.0
        self.assertEqual(actual, expected)

    def test_weighted_minimal_trees(self):
        # expected values computed by hand
        # zero tips
        tree = TreeNode.read(StringIO(u'root;'))
        actual = weighted_unifrac([], [], [], tree)
        expected = 0.0
        self.assertEqual(actual, expected)

        # two tips
        tree = TreeNode.read(StringIO(u'(OTU1:0.25, OTU2:0.25)root;'))
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

    def test_unweighted_unifrac_kwargs(self):
        # confirm that **kwargs can be passed
        actual = unweighted_unifrac(self.table1[0], self.table1[0], self.oids1,
                                    self.t1, not_a_known_parameter=42)
        self.assertAlmostEqual(actual, 0.0)

    def test_weighted_unifrac_kwargs(self):
        # confirm that **kwargs can be passed
        actual = weighted_unifrac(self.table1[0], self.table1[0], self.oids1,
                                  self.t1, not_a_known_parameter=42)
        self.assertAlmostEqual(actual, 0.0)

    def test_unweighted_unifrac_identity(self):
        for i in range(len(self.table1)):
            actual = unweighted_unifrac(
                self.table1[i], self.table1[i], self.oids1, self.t1)
            expected = 0.0
            self.assertAlmostEqual(actual, expected)

    def test_unweighted_unifrac_symmetry(self):
        for i in range(len(self.table1)):
            for j in range(len(self.table1)):
                actual = unweighted_unifrac(
                    self.table1[i], self.table1[j], self.oids1, self.t1)
                expected = unweighted_unifrac(
                    self.table1[j], self.table1[i], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_invalid_input(self):
        # Many of these tests are duplicated from
        # skbio.diversity.tests.test_base, but I think it's important to
        # confirm that they are being run when *unifrac is called.

        # tree has duplicated tip ids
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU2:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(DuplicateNodeError, unweighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        self.assertRaises(DuplicateNodeError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)

        # unrooted tree as input
        t = TreeNode.read(StringIO(u'((OTU1:0.1, OTU2:0.2):0.3, OTU3:0.5,'
                                   u'OTU4:0.7);'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts,
                          v_counts, otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts,
                          v_counts, otu_ids, t)

        # otu_ids has duplicated ids
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU2']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts, v_counts,
                          otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts, v_counts,
                          otu_ids, t)

        # len of vectors not equal
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts, v_counts,
                          otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts, v_counts,
                          otu_ids, t)
        u_counts = [1, 2, 3]
        v_counts = [1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts, v_counts,
                          otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts, v_counts,
                          otu_ids, t)
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts, v_counts,
                          otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts, v_counts,
                          otu_ids, t)

        # negative counts
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, -3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts, v_counts,
                          otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts, v_counts,
                          otu_ids, t)
        u_counts = [1, 2, 3]
        v_counts = [1, 1, -1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts, v_counts,
                          otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts, v_counts,
                          otu_ids, t)

        # tree with no branch lengths
        t = TreeNode.read(
            StringIO(u'((((OTU1,OTU2),OTU3)),(OTU4,OTU5));'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts, v_counts,
                          otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts, v_counts,
                          otu_ids, t)

        # tree missing some branch lengths
        t = TreeNode.read(
            StringIO(u'(((((OTU1,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        u_counts = [1, 2, 3]
        v_counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, unweighted_unifrac, u_counts, v_counts,
                          otu_ids, t)
        self.assertRaises(ValueError, weighted_unifrac, u_counts, v_counts,
                          otu_ids, t)

        # otu_ids not present in tree
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
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
            self.table1[4], self.table1[5], self.oids1, self.t1)
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
            self.table1[0], self.table1[1], self.oids1, self.t1)
        expected = 0.238095238095
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.table1[0], self.table1[2], self.oids1, self.t1)
        expected = 0.52
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.table1[0], self.table1[3], self.oids1, self.t1)
        expected = 0.52
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.table1[0], self.table1[4], self.oids1, self.t1)
        expected = 0.545454545455
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.table1[0], self.table1[5], self.oids1, self.t1)
        expected = 0.619047619048
        self.assertAlmostEqual(actual, expected)
        # sample B versus remaining
        actual = unweighted_unifrac(
            self.table1[1], self.table1[2], self.oids1, self.t1)
        expected = 0.347826086957
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.table1[1], self.table1[3], self.oids1, self.t1)
        expected = 0.347826086957
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.table1[1], self.table1[4], self.oids1, self.t1)
        expected = 0.68
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.table1[1], self.table1[5], self.oids1, self.t1)
        expected = 0.421052631579
        self.assertAlmostEqual(actual, expected)
        # sample C versus remaining
        actual = unweighted_unifrac(
            self.table1[2], self.table1[3], self.oids1, self.t1)
        expected = 0.0
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.table1[2], self.table1[4], self.oids1, self.t1)
        expected = 0.68
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.table1[2], self.table1[5], self.oids1, self.t1)
        expected = 0.421052631579
        self.assertAlmostEqual(actual, expected)
        # sample D versus remaining
        actual = unweighted_unifrac(
            self.table1[3], self.table1[4], self.oids1, self.t1)
        expected = 0.68
        self.assertAlmostEqual(actual, expected)
        actual = unweighted_unifrac(
            self.table1[3], self.table1[5], self.oids1, self.t1)
        expected = 0.421052631579
        self.assertAlmostEqual(actual, expected)
        # sample E versus remaining
        actual = unweighted_unifrac(
            self.table1[4], self.table1[5], self.oids1, self.t1)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_identity(self):
        for i in range(len(self.table1)):
            actual = weighted_unifrac(
                self.table1[i], self.table1[i], self.oids1, self.t1)
            expected = 0.0
            self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_symmetry(self):
        for i in range(len(self.table1)):
            for j in range(len(self.table1)):
                actual = weighted_unifrac(
                    self.table1[i], self.table1[j], self.oids1, self.t1)
                expected = weighted_unifrac(
                    self.table1[j], self.table1[i], self.oids1, self.t1)
                self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_non_overlapping(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation skbio's initial
        # weighted unifrac implementation
        # these communities only share the root node
        actual = weighted_unifrac(
            self.table1[4], self.table1[5], self.oids1, self.t1)
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
            self.table1[0], self.table1[1], self.oids1, self.t1)
        expected = 2.4
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[0], self.table1[2], self.oids1, self.t1)
        expected = 1.86666666667
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[0], self.table1[3], self.oids1, self.t1)
        expected = 2.53333333333
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[0], self.table1[4], self.oids1, self.t1)
        expected = 1.35384615385
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[0], self.table1[5], self.oids1, self.t1)
        expected = 3.2
        self.assertAlmostEqual(actual, expected)
        # sample B versus remaining
        actual = weighted_unifrac(
            self.table1[1], self.table1[2], self.oids1, self.t1)
        expected = 2.26666666667
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[1], self.table1[3], self.oids1, self.t1)
        expected = 0.933333333333
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[1], self.table1[4], self.oids1, self.t1)
        expected = 3.2
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[1], self.table1[5], self.oids1, self.t1)
        expected = 0.8375
        self.assertAlmostEqual(actual, expected)
        # sample C versus remaining
        actual = weighted_unifrac(
            self.table1[2], self.table1[3], self.oids1, self.t1)
        expected = 1.33333333333
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[2], self.table1[4], self.oids1, self.t1)
        expected = 1.89743589744
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[2], self.table1[5], self.oids1, self.t1)
        expected = 2.66666666667
        self.assertAlmostEqual(actual, expected)
        # sample D versus remaining
        actual = weighted_unifrac(
            self.table1[3], self.table1[4], self.oids1, self.t1)
        expected = 2.66666666667
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[3], self.table1[5], self.oids1, self.t1)
        expected = 1.33333333333
        self.assertAlmostEqual(actual, expected)
        # sample E versus remaining
        actual = weighted_unifrac(
            self.table1[4], self.table1[5], self.oids1, self.t1)
        expected = 4.0
        self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_identity_normalized(self):
        for i in range(len(self.table1)):
            actual = weighted_unifrac(
                self.table1[i], self.table1[i], self.oids1, self.t1,
                normalized=True)
            expected = 0.0
            self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_symmetry_normalized(self):
        for i in range(len(self.table1)):
            for j in range(len(self.table1)):
                actual = weighted_unifrac(
                    self.table1[i], self.table1[j], self.oids1, self.t1,
                    normalized=True)
                expected = weighted_unifrac(
                    self.table1[j], self.table1[i], self.oids1, self.t1,
                    normalized=True)
                self.assertAlmostEqual(actual, expected)

    def test_weighted_unifrac_non_overlapping_normalized(self):
        # these communities only share the root node
        actual = weighted_unifrac(
            self.table1[4], self.table1[5], self.oids1, self.t1,
            normalized=True)
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
            self.table1[0], self.table1[1], self.oids1, self.t1,
            normalized=True)
        expected = 0.6
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[0], self.table1[2], self.oids1, self.t1,
            normalized=True)
        expected = 0.466666666667
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[0], self.table1[3], self.oids1, self.t1,
            normalized=True)
        expected = 0.633333333333
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[0], self.table1[4], self.oids1, self.t1,
            normalized=True)
        expected = 0.338461538462
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[0], self.table1[5], self.oids1, self.t1,
            normalized=True)
        expected = 0.8
        self.assertAlmostEqual(actual, expected)
        # sample B versus remaining
        actual = weighted_unifrac(
            self.table1[1], self.table1[2], self.oids1, self.t1,
            normalized=True)
        expected = 0.566666666667
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[1], self.table1[3], self.oids1, self.t1,
            normalized=True)
        expected = 0.233333333333
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[1], self.table1[4], self.oids1, self.t1,
            normalized=True)
        expected = 0.8
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[1], self.table1[5], self.oids1, self.t1,
            normalized=True)
        expected = 0.209375
        self.assertAlmostEqual(actual, expected)
        # sample C versus remaining
        actual = weighted_unifrac(
            self.table1[2], self.table1[3], self.oids1, self.t1,
            normalized=True)
        expected = 0.333333333333
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[2], self.table1[4], self.oids1, self.t1,
            normalized=True)
        expected = 0.474358974359
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[2], self.table1[5], self.oids1, self.t1,
            normalized=True)
        expected = 0.666666666667
        self.assertAlmostEqual(actual, expected)
        # sample D versus remaining
        actual = weighted_unifrac(
            self.table1[3], self.table1[4], self.oids1, self.t1,
            normalized=True)
        expected = 0.666666666667
        self.assertAlmostEqual(actual, expected)
        actual = weighted_unifrac(
            self.table1[3], self.table1[5], self.oids1, self.t1,
            normalized=True)
        expected = 0.333333333333
        self.assertAlmostEqual(actual, expected)
        # sample E versus remaining
        actual = weighted_unifrac(
            self.table1[4], self.table1[5], self.oids1, self.t1,
            normalized=True)
        expected = 1.0
        self.assertAlmostEqual(actual, expected)

if __name__ == "__main__":
    main()
