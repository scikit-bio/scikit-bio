from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

from skbio.core.distance import DistanceMatrix
from skbio.core.tree import TreeNode
from skbio.core.tree.nj import (nj, _compute_q, _compute_collapsed_dm,
    _lowest_index, _otu_to_new_node, _pair_members_to_new_node)

class NjTests(TestCase):

    def setUp(self):
        data1 = [[0,  5,  9,  9,  8],
                 [5,  0, 10, 10,  9],
                 [9, 10,  0,  8,  7],
                 [9, 10,  8,  0,  3],
                 [8,  9,  7,  3,  0]]
        ids1 = list('abcde')
        self.dm1 = DistanceMatrix(data1, ids1)
        # this newick string was confirmed against http://www.trex.uqam.ca/
        # which generated the following (isomorphic) newick string:
        # (d:2.0000,e:1.0000,(c:4.0000,(a:2.0000,b:3.0000):3.0000):2.0000);
        self.expected1_str = "(d:2, (c:4, (b:3, a:2):3):2, e:1);"
        self.expected1_TreeNode = TreeNode.from_newick(self.expected1_str)

    def test_nj(self):
        self.assertEqual(nj(self.dm1, result_constructor=str),
                         self.expected1_str)
        # what is the correct way to compare TreeNode objects for equality?
        actual_TreeNode = nj(self.dm1)
        self.assertEqual(\
            actual_TreeNode.compare_tip_distances(self.expected1_TreeNode), 0.0)

    def test_nj_trivial(self):
        data = [[0, 3, 2],
                [3, 0, 3],
                [2, 3, 0]]
        dm = DistanceMatrix(data, list('abc'))
        expected_str = "(b:2, a:1, c:1);"
        self.assertEqual(nj(dm, result_constructor=str), expected_str)

    def test_compute_q(self):
        expected_data = [[  0, -50, -38, -34, -34],
                         [-50,   0, -38, -34, -34],
                         [-38, -38,   0, -40, -40],
                         [-34, -34, -40,   0, -48],
                         [-34, -34, -40, -48,   0]]
        expected_ids = list('abcde')
        expected = DistanceMatrix(expected_data, expected_ids)
        self.assertEqual(_compute_q(self.dm1), expected)

    def test_compute_collapsed_dm(self):
        expected_data = [[0,  7,  7,  6],
                         [7,  0,  8,  7],
                         [7,  8,  0,  3],
                         [6,  7,  3,  0]]
        expected_ids = ['x', 'c', 'd', 'e']
        expected = DistanceMatrix(expected_data, expected_ids)
        self.assertEqual(_compute_collapsed_dm(self.dm1, 'a', 'b', True, 'x'),
                         expected)

    def test_lowest_index(self):
        self.assertEqual(_lowest_index(self.dm1), (4, 3))
        self.assertEqual(_lowest_index(_compute_q(self.dm1)), (1, 0))

    def test_otu_to_new_node(self):
        self.assertEqual(_otu_to_new_node(self.dm1, 'a', 'b', 'c', True), 7)

    def test_pair_members_to_new_node(self):
        self.assertEqual(_pair_members_to_new_node(self.dm1, 'a', 'b', True),
                         (2, 3))
