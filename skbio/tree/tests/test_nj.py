# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

from skbio import DistanceMatrix, TreeNode, nj
from skbio.tree._nj import (
    _compute_q, _compute_collapsed_dm, _lowest_index,
    _pair_members_to_new_node)


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
        self.expected1_str = ("(d:2.000000, (c:4.000000, (b:3.000000,"
                              " a:2.000000):3.000000):2.000000, e:1.000000);")
        self.expected1_TreeNode = TreeNode.read(
                io.StringIO(self.expected1_str))

        # this example was pulled from the Phylip manual
        # http://evolution.genetics.washington.edu/phylip/doc/neighbor.html
        data2 = [[0.0000, 1.6866, 1.7198, 1.6606, 1.5243, 1.6043, 1.5905],
                 [1.6866, 0.0000, 1.5232, 1.4841, 1.4465, 1.4389, 1.4629],
                 [1.7198, 1.5232, 0.0000, 0.7115, 0.5958, 0.6179, 0.5583],
                 [1.6606, 1.4841, 0.7115, 0.0000, 0.4631, 0.5061, 0.4710],
                 [1.5243, 1.4465, 0.5958, 0.4631, 0.0000, 0.3484, 0.3083],
                 [1.6043, 1.4389, 0.6179, 0.5061, 0.3484, 0.0000, 0.2692],
                 [1.5905, 1.4629, 0.5583, 0.4710, 0.3083, 0.2692, 0.0000]]
        ids2 = ["Bovine", "Mouse", "Gibbon", "Orang", "Gorilla", "Chimp",
                "Human"]
        self.dm2 = DistanceMatrix(data2, ids2)
        self.expected2_str = ("(Mouse:0.76891, (Gibbon:0.35793, (Orang:0.28469"
                              ", (Gorilla:0.15393, (Chimp:0.15167, Human:0.117"
                              "53):0.03982):0.02696):0.04648):0.42027, Bovine:"
                              "0.91769);")
        self.expected2_TreeNode = TreeNode.read(
                io.StringIO(self.expected2_str))

        data3 = [[0, 5, 4, 7, 6, 8],
                 [5, 0, 7, 10, 9, 11],
                 [4, 7, 0, 7, 6, 8],
                 [7, 10, 7, 0, 5, 8],
                 [6, 9, 6, 5, 0, 8],
                 [8, 11, 8, 8, 8, 0]]
        ids3 = map(str, range(6))
        self.dm3 = DistanceMatrix(data3, ids3)
        self.expected3_str = ("((((0:1.000000,1:4.000000):1.000000,2:2.000000"
                              "):1.250000,5:4.750000):0.750000,3:2.750000,4:2."
                              "250000);")
        self.expected3_TreeNode = TreeNode.read(
                io.StringIO(self.expected3_str))

        # this dm can yield negative branch lengths
        data4 = [[0,  5,  9,  9,  800],
                 [5,  0, 10, 10,  9],
                 [9, 10,  0,  8,  7],
                 [9, 10,  8,  0,  3],
                 [800,  9,  7,  3,  0]]
        ids4 = list('abcde')
        self.dm4 = DistanceMatrix(data4, ids4)

    def test_nj_dm1(self):
        self.assertEqual(nj(self.dm1, result_constructor=str),
                         self.expected1_str)
        # what is the correct way to compare TreeNode objects for equality?
        actual_TreeNode = nj(self.dm1)
        self.assertEqual(actual_TreeNode.compare_tip_distances(
            self.expected1_TreeNode), 0.0)

    def test_nj_dm2(self):
        actual_TreeNode = nj(self.dm2)
        self.assertAlmostEqual(actual_TreeNode.compare_tip_distances(
            self.expected2_TreeNode), 0.0)

    def test_nj_dm3(self):
        actual_TreeNode = nj(self.dm3)
        self.assertAlmostEqual(actual_TreeNode.compare_tip_distances(
            self.expected3_TreeNode), 0.0)

    def test_nj_zero_branch_length(self):
        # no nodes have negative branch length when we disallow negative
        # branch length. self is excluded as branch length is None
        tree = nj(self.dm4)
        for n in tree.postorder(include_self=False):
            self.assertTrue(n.length >= 0)
        # only tips associated with the large distance in the input
        # have positive branch lengths when we allow negative branch
        # length
        tree = nj(self.dm4, False)
        self.assertTrue(tree.find('a').length > 0)
        self.assertTrue(tree.find('b').length < 0)
        self.assertTrue(tree.find('c').length < 0)
        self.assertTrue(tree.find('d').length < 0)
        self.assertTrue(tree.find('e').length > 0)

    def test_nj_trivial(self):
        data = [[0, 3, 2],
                [3, 0, 3],
                [2, 3, 0]]
        dm = DistanceMatrix(data, list('abc'))
        expected_str = "(b:2.000000, a:1.000000, c:1.000000);"
        self.assertEqual(nj(dm, result_constructor=str), expected_str)

    def test_nj_error(self):
        data = [[0, 3],
                [3, 0]]
        dm = DistanceMatrix(data, list('ab'))
        self.assertRaises(ValueError, nj, dm)

    def test_compute_q(self):
        expected_data = [[0, -50, -38, -34, -34],
                         [-50,   0, -38, -34, -34],
                         [-38, -38,   0, -40, -40],
                         [-34, -34, -40,   0, -48],
                         [-34, -34, -40, -48,   0]]
        expected_ids = list('abcde')
        expected = DistanceMatrix(expected_data, expected_ids)
        self.assertEqual(_compute_q(self.dm1), expected)

        data = [[0, 3, 2],
                [3, 0, 3],
                [2, 3, 0]]
        dm = DistanceMatrix(data, list('abc'))
        # computed this manually
        expected_data = [[0, -8, -8],
                         [-8,  0, -8],
                         [-8, -8,  0]]
        expected = DistanceMatrix(expected_data, list('abc'))
        self.assertEqual(_compute_q(dm), expected)

    def test_compute_collapsed_dm(self):
        expected_data = [[0,  7,  7,  6],
                         [7,  0,  8,  7],
                         [7,  8,  0,  3],
                         [6,  7,  3,  0]]
        expected_ids = ['x', 'c', 'd', 'e']
        expected1 = DistanceMatrix(expected_data, expected_ids)
        self.assertEqual(_compute_collapsed_dm(self.dm1, 'a', 'b', True, 'x'),
                         expected1)

        # computed manually
        expected_data = [[0, 4, 3],
                         [4, 0, 3],
                         [3, 3, 0]]
        expected_ids = ['yy', 'd', 'e']
        expected2 = DistanceMatrix(expected_data, expected_ids)
        self.assertEqual(
            _compute_collapsed_dm(expected1, 'x', 'c', True, 'yy'), expected2)

    def test_lowest_index(self):
        self.assertEqual(_lowest_index(self.dm1), (4, 3))
        self.assertEqual(_lowest_index(_compute_q(self.dm1)), (1, 0))

    def test_pair_members_to_new_node(self):
        self.assertEqual(_pair_members_to_new_node(self.dm1, 'a', 'b', True),
                         (2, 3))
        self.assertEqual(_pair_members_to_new_node(self.dm1, 'a', 'c', True),
                         (4, 5))
        self.assertEqual(_pair_members_to_new_node(self.dm1, 'd', 'e', True),
                         (2, 1))

    def test_pair_members_to_new_node_zero_branch_length(self):
        # the values in this example don't really make sense
        # (I'm not sure how you end up with these distances between
        # three sequences), but that doesn't really matter for the sake
        # of this test
        data = [[0, 4, 2],
                [4, 0, 38],
                [2, 38, 0]]
        ids = ['a', 'b', 'c']
        dm = DistanceMatrix(data, ids)
        self.assertEqual(_pair_members_to_new_node(dm, 'a', 'b', True), (0, 4))
        # this makes it clear why negative branch lengths don't make sense...
        self.assertEqual(
            _pair_members_to_new_node(dm, 'a', 'b', False), (-16, 20))


if __name__ == "__main__":
    main()
