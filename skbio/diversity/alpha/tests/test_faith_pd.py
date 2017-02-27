# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from io import StringIO
import os

import numpy as np
import pandas as pd

from skbio import TreeNode
from skbio.util import get_data_path
from skbio.tree import DuplicateNodeError, MissingNodeError
from skbio.diversity.alpha import faith_pd


class FaithPDTests(TestCase):

    def setUp(self):
        self.counts = np.array([0, 1, 1, 4, 2, 5, 2, 4, 1, 2])
        self.b1 = np.array([[1, 3, 0, 1, 0],
                            [0, 2, 0, 4, 4],
                            [0, 0, 6, 2, 1],
                            [0, 0, 1, 1, 1]])
        self.sids1 = list('ABCD')
        self.oids1 = ['OTU%d' % i for i in range(1, 6)]
        self.t1 = TreeNode.read(StringIO(
            '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):'
            '0.0,(OTU4:0.75,OTU5:0.75):1.25):0.0)root;'))
        self.t1_w_extra_tips = TreeNode.read(
           StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                    '0.75,(OTU5:0.25,(OTU6:0.5,OTU7:0.5):0.5):0.5):1.25):0.0'
                    ')root;'))

    def test_faith_pd_none_observed(self):
        actual = faith_pd(np.array([], dtype=int),
                          np.array([], dtype=int),
                          self.t1)
        expected = 0.0
        self.assertAlmostEqual(actual, expected)
        actual = faith_pd([0, 0, 0, 0, 0], self.oids1, self.t1)
        expected = 0.0
        self.assertAlmostEqual(actual, expected)

    def test_faith_pd_all_observed(self):
        actual = faith_pd([1, 1, 1, 1, 1], self.oids1, self.t1)
        expected = sum(n.length for n in self.t1.traverse()
                       if n.length is not None)
        self.assertAlmostEqual(actual, expected)

        actual = faith_pd([1, 2, 3, 4, 5], self.oids1, self.t1)
        expected = sum(n.length for n in self.t1.traverse()
                       if n.length is not None)
        self.assertAlmostEqual(actual, expected)

    def test_faith_pd(self):
        # expected results derived from QIIME 1.9.1, which
        # is a completely different implementation skbio's initial
        # phylogenetic diversity implementation
        actual = faith_pd(self.b1[0], self.oids1, self.t1)
        expected = 4.5
        self.assertAlmostEqual(actual, expected)
        actual = faith_pd(self.b1[1], self.oids1, self.t1)
        expected = 4.75
        self.assertAlmostEqual(actual, expected)
        actual = faith_pd(self.b1[2], self.oids1, self.t1)
        expected = 4.75
        self.assertAlmostEqual(actual, expected)
        actual = faith_pd(self.b1[3], self.oids1, self.t1)
        expected = 4.75
        self.assertAlmostEqual(actual, expected)

    def test_faith_pd_extra_tips(self):
        # results are the same despite presences of unobserved tips in tree
        actual = faith_pd(self.b1[0], self.oids1, self.t1_w_extra_tips)
        expected = faith_pd(self.b1[0], self.oids1, self.t1)
        self.assertAlmostEqual(actual, expected)
        actual = faith_pd(self.b1[1], self.oids1, self.t1_w_extra_tips)
        expected = faith_pd(self.b1[1], self.oids1, self.t1)
        self.assertAlmostEqual(actual, expected)
        actual = faith_pd(self.b1[2], self.oids1, self.t1_w_extra_tips)
        expected = faith_pd(self.b1[2], self.oids1, self.t1)
        self.assertAlmostEqual(actual, expected)
        actual = faith_pd(self.b1[3], self.oids1, self.t1_w_extra_tips)
        expected = faith_pd(self.b1[3], self.oids1, self.t1)
        self.assertAlmostEqual(actual, expected)

    def test_faith_pd_minimal(self):
        # two tips
        tree = TreeNode.read(StringIO('(OTU1:0.25, OTU2:0.25)root;'))
        actual = faith_pd([1, 0], ['OTU1', 'OTU2'], tree)
        expected = 0.25
        self.assertEqual(actual, expected)

    def test_faith_pd_qiime_tiny_test(self):
        # the following table and tree are derived from the QIIME 1.9.1
        # "tiny-test" data
        tt_table_fp = get_data_path(
            os.path.join('qiime-191-tt', 'otu-table.tsv'), 'data')
        tt_tree_fp = get_data_path(
            os.path.join('qiime-191-tt', 'tree.nwk'), 'data')

        self.q_table = pd.read_csv(tt_table_fp, sep='\t', skiprows=1,
                                   index_col=0)
        self.q_tree = TreeNode.read(tt_tree_fp)

        expected_fp = get_data_path(
            os.path.join('qiime-191-tt', 'faith-pd.txt'), 'data')
        expected = pd.read_csv(expected_fp, sep='\t', index_col=0)
        for sid in self.q_table.columns:
            actual = faith_pd(self.q_table[sid],
                              otu_ids=self.q_table.index,
                              tree=self.q_tree)
            self.assertAlmostEqual(actual, expected['PD_whole_tree'][sid])

    def test_faith_pd_root_not_observed(self):
        # expected values computed by hand
        tree = TreeNode.read(
            StringIO('((OTU1:0.1, OTU2:0.2):0.3, (OTU3:0.5, OTU4:0.7):1.1)'
                     'root;'))
        otu_ids = ['OTU%d' % i for i in range(1, 5)]
        # root node not observed, but branch between (OTU1, OTU2) and root
        # is considered observed
        actual = faith_pd([1, 1, 0, 0], otu_ids, tree)
        expected = 0.6
        self.assertAlmostEqual(actual, expected)

        # root node not observed, but branch between (OTU3, OTU4) and root
        # is considered observed
        actual = faith_pd([0, 0, 1, 1], otu_ids, tree)
        expected = 2.3
        self.assertAlmostEqual(actual, expected)

    def test_faith_pd_invalid_input(self):
        # tree has duplicated tip ids
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU2:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(DuplicateNodeError, faith_pd, counts, otu_ids,
                          t)

        # unrooted tree as input
        t = TreeNode.read(StringIO('((OTU1:0.1, OTU2:0.2):0.3, OTU3:0.5,'
                                   'OTU4:0.7);'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, faith_pd, counts, otu_ids, t)

        # otu_ids has duplicated ids
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU2']
        self.assertRaises(ValueError, faith_pd, counts, otu_ids, t)

        # len of vectors not equal
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, faith_pd, counts, otu_ids, t)
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2']
        self.assertRaises(ValueError, faith_pd, counts, otu_ids, t)

        # negative counts
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, -3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, faith_pd, counts, otu_ids, t)

        # tree with no branch lengths
        t = TreeNode.read(
            StringIO('((((OTU1,OTU2),OTU3)),(OTU4,OTU5));'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, faith_pd, counts, otu_ids, t)

        # tree missing some branch lengths
        t = TreeNode.read(
            StringIO('(((((OTU1,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, faith_pd, counts, otu_ids, t)

        # otu_ids not present in tree
        t = TreeNode.read(
            StringIO('(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU42']
        self.assertRaises(MissingNodeError, faith_pd, counts, otu_ids, t)


if __name__ == "__main__":
    main()
