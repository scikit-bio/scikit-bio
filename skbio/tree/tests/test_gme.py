# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
import numpy as np
from unittest import TestCase, main

from skbio import DistanceMatrix, TreeNode
from skbio.tree._gme import (
    _average_distance_k,
    _lower_subtree_list, _upper_subtree_list,
    _edge_attachment_length,
    _average_distance, _tip_or_root, _tip_or_root_upper,
    _average_distance_upper, _subtree_count,
    _average_subtree_distance, _average_distance_matrix, 
    _ols_edge, gme)


class GmeTests(TestCase):

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
        self.expected1_str = ("((b:3.0,(c:4.0,(d:2.0,e:1.0):2.0):3.0):2.0)a;")
        self.expected1_TreeNode = TreeNode.read(
                io.StringIO(self.expected1_str))
        self.gme_starting1_tree = np.array([[1, 2, 0, 0], [1, 2, 0, 0]])
        self.gme_starting1_ordered = np.array([[1, 2, 0], [1, 2, 0]])
        self.gme_starting1_tips = np.array([[1, 2], [1, 2]])
        self.gme_attach_node = 3

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
        self.expected2_str = ("((Mouse:0.7689100000000002,(Gibbon:0.37287500000000007,"
                              "(Orang:0.2705916666666668,(Gorilla:0.1525416666666667,"
                              "(Chimp:0.14516249999999986,Human:0.12403750000000013)"
                              ":0.04120833333333351):0.0393958333333333)"
                              ":0.032630555555555496):0.42026874999999997)"
                              ":0.9176899999999998)Bovine;")
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
        self.expected3_str = ("((1:4.0,(2:2.0,((3:3.0,4:2.0):0.75,5:4.75):1.25):1.0):1.0)0;")
        self.expected3_TreeNode = TreeNode.read(
                io.StringIO(self.expected3_str))

        # this dm can yield negative branch lengths for OLS-based estimation
        data4 = [[0,  5,  9,  9,  800],
                 [5,  0, 10, 10,  9],
                 [9, 10,  0,  8,  7],
                 [9, 10,  8,  0,  3],
                 [800,  9,  7,  3,  0]]
        ids4 = list('abcde')
        self.dm4 = DistanceMatrix(data4, ids4)

    def test_gme_dm1(self):
        actual_TreeNode = gme(self.dm1)
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            self.expected1_TreeNode), 0.0, places=10)

    def test_gme_dm2(self):
        actual_TreeNode = gme(self.dm2)
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            self.expected2_TreeNode), 0.0, places=10)

    def test_gme_dm3(self):
        actual_TreeNode = gme(self.dm3)
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            self.expected3_TreeNode), 0.0, places=10)

    def test_gme_zero_branch_length(self):
        # OLS-based edge estimation can produce negative branch 
        # lengths when some dm values are much larger than
        # others, analogous to negative branch lengths produced by nj
        tree = gme(self.dm4)
        self.assertTrue(tree.find('b').length < 0)
        self.assertTrue(tree.find('c').length < 0)
        self.assertTrue(tree.find('d').length < 0)
        self.assertTrue(tree.find('e').length > 0)

    def test_gme_error(self):
        data = [[0, 3],
                [3, 0]]
        dm = DistanceMatrix(data, list('ab'))
        self.assertRaises(ValueError, gme, dm)

    def test_average_distance_k(self):
        nodesubtree = 0
        self.assertAlmostEqual(_average_distance_k(
            self.gme_attach_node, nodesubtree, 
            self.gme_starting1_ordered,
            self.gme_starting1_tips, self.dm1, 0
            ), 9.0, places=10)

    def test_average_distance_k_upper(self):
        nodesubtree = 0
        self.assertAlmostEqual(_average_distance_k(
            self.gme_attach_node, nodesubtree, 
            self.gme_starting1_ordered,
            self.gme_starting1_tips, self.dm1, 0,
            upper=True
            ), 9.0, places=10)

    def test_lower_subtree_list(self):
        expected_list = [10.0, 8.0, 9.0]
        self.assertEqual(_lower_subtree_list(
            self.gme_attach_node, self.gme_starting1_ordered,
            self.gme_starting1_tips, self.dm1, 0), expected_list)

    def test_upper_subtree_list(self):
        expected_list = [8.5, 9.5, 9.0]
        self.assertEqual(_upper_subtree_list(
            self.gme_attach_node, self.gme_starting1_ordered,
            self.gme_starting1_tips, self.dm1, 0), expected_list)

    def test_edge_attachment_length(self):
        child = 2
        lowerlist = _lower_subtree_list(
            self.gme_attach_node, self.gme_starting1_ordered,
            self.gme_starting1_tips, self.dm1, 0)
        upperlist = _upper_subtree_list(
            self.gme_attach_node, self.gme_starting1_ordered,
            self.gme_starting1_tips, self.dm1, 0)
        adm = _average_distance_matrix(
            self.dm1, self.gme_starting1_ordered, self.gme_starting1_tips, 0)
        self.assertEqual(_edge_attachment_length(
            child, lowerlist, upperlist,
            self.gme_starting1_ordered, self.gme_starting1_tips, adm), -1.5)

    def test_average_distance(self):
        ordered = np.array([[2, 0, 4, 1, 8, 3, 7], [1, 0, 2, 0, 3, 0, 4]])
        tips = np.array([[2, 4, 8, 7], [1, 2, 3, 4]])
        node_1 = 2
        node_2 = 3
        self.assertAlmostEqual(_average_distance(
            node_1, node_2, self.dm1, ordered, tips),
                               9.5, places=10)

    def test_average_distance_upper(self):
        ordered = np.array([[2, 0, 4, 1, 8, 3, 7], [1, 0, 2, 0, 3, 0, 4]])
        tips = np.array([[2, 4, 8, 7], [1, 2, 3, 4]])
        node_1 = 2
        node_2 = 3
        self.assertAlmostEqual(_average_distance_upper(
            node_1, node_2, self.dm1, ordered, tips, 0),
                               5.0, places=10)

    def test_tip_or_root(self):
        ordered = np.array([[2, 0, 4, 1, 8, 3, 7], [1, 0, 2, 0, 3, 0, 4]])
        tips = np.array([[2, 4, 8, 7], [1, 2, 3, 4]])
        node_internal = 3
        node_leaf = 2
        self.assertEqual(len(_tip_or_root(node_internal, ordered, tips)), 2)
        self.assertEqual(len(_tip_or_root(node_leaf, ordered, tips)), 1)
       
    def test_tip_or_root_upper(self):
        ordered = np.array([[2, 0, 4, 1, 8, 3, 7], [1, 0, 2, 0, 3, 0, 4]])
        tips = np.array([[2, 4, 8, 7], [1, 2, 3, 4]])
        node_internal = 3
        node_leaf = 2
        self.assertEqual(len(_tip_or_root_upper(node_internal, ordered, tips, 0)), 3)
        self.assertEqual(len(_tip_or_root_upper(node_leaf, ordered, tips, 0)), 1)

    def test_subtree_count(self):
        internal_node = 1
        leaf = 8
        root = 0
        tips = np.array([[2, 4, 8, 7], [1, 2, 3, 4]])
        self.assertEqual(_subtree_count(internal_node, tips), 3)
        self.assertEqual(_subtree_count(leaf, tips), 1)
        self.assertEqual(_subtree_count(root, tips), 4)

    def test_average_subtree_distance(self):
        # computed manually
        tips = np.array([[3, 6, 4, 5], [1, 2, 3, 4]])
        ordered = np.array([[3, 6, 2, 0, 4, 1, 5], [1, 2, 0, 0, 3, 0, 4]])
        a = 2
        b = 3
        a1 = 5
        a2 = 6
        self.assertAlmostEqual(_average_subtree_distance(
            a, b, a1, a2, self.dm1, ordered, tips),
                               9.5, places=10)

    def test_average_distance_matrix_trivial(self):
        # In this case, the average distance matrix is equivalent to
        # the original distance matrix
        data = [[0, 3, 2],
                [3, 0, 3],
                [2, 3, 0]]
        ids = list('abc')
        dm = DistanceMatrix(data, ids)
        tips = np.array([[2, 1], [1, 2]])
        ordered = np.array([[2, 0, 1], [1, 0, 2]])
        index = [0, 1, 2]
        actual_adm = _average_distance_matrix(dm, ordered, tips, 0)
        for i in index:
            for j in index:
                if j < i:
                    self.assertEqual(dm[i][j], actual_adm[i][j])
                    self.assertEqual(dm[j][i], actual_adm[j][i])

    def test_average_distance_matrix(self):
        # computed manually
        tips = np.array([[3, 6, 4, 5], [1, 2, 3, 4]])
        ordered = np.array([[3, 6, 2, 0, 4, 1, 5], [1, 2, 0, 0, 3, 0, 4]])
        expected_adm = [[0.0, 10.0, 9.5, 9.0, 10.0, 8.0, 5.0],
                        [10.0, 0.0, 9.0, 3.0, 8.0, 9.0, 9.0],
                        [9.5, 9.0, 0.0, 6.0, 5.5, 7.5, 8.5],
                        [9.0, 3.0, 6.0, 0.0, 7.0, 6.666666666666667, 8.0],
                        [10.0, 8.0, 5.5, 7.0, 0.0, 6.666666666666667, 9.0],
                        [8.0, 9.0, 7.5, 6.666666666666667, 6.666666666666667, 0.0, 7.0],
                        [5.0, 9.0, 8.5, 8.0, 9.0, 7.0, 0.0]]
        actual_adm = _average_distance_matrix(self.dm1, ordered, tips, 0)
        index = [0, 1, 2, 3, 4, 5, 6]
        for i in index:
            for j in index:
                if j < 1:
                    self.assertAlmostEqual(expected_adm[i][j], actual_adm[i][j])
                    self.assertAlmostEqual(expected_adm[j][i], actual_adm[j][i])

    def test_edge_estimation(self):
        data = [[0, 3, 2],
                [3, 0, 3],
                [2, 3, 0]]
        ids = list('abc')
        dm = DistanceMatrix(data, ids)
        tips = np.array([[1, 2], [1, 2]])
        ordered = np.array([[1, 2, 0], [1, 2, 0]])
        pre_estimation_str = "((b,c))a;"
        expected_str = "((b:2.0,c:1.0):1.0)a;"
        actual_TreeNode = TreeNode.read(io.StringIO(pre_estimation_str))
        _ols_edge(actual_TreeNode, dm, ordered, tips, 0)
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            expected_TreeNode), 0.0, places=10)

if __name__ == "__main__":
    main()
