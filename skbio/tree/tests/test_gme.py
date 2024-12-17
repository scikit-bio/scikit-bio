# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np

from skbio import DistanceMatrix, TreeNode
from skbio.tree._gme import (
    _average_distance_k,
    _average_distance_k_upper,
    _lower_subtree_list,
    _upper_subtree_list,
    _edge_attachment_length,
    _average_distance,
    _tip_or_root,
    _average_distance_upper,
    _subtree_count,
    _average_subtree_distance,
    _average_distance_matrix, 
    _ols_edge,
    gme,
)


class GmeTests(TestCase):

    def setUp(self):
        # Example 1
        # This newick string was confirmed against http://www.trex.uqam.ca/ which
        # generated the following (isomorphic) newick string:
        # (d:2.0000,e:1.0000,(c:4.0000,(a:2.0000,b:3.0000):3.0000):2.0000);
        data1 = [[0,  5,  9,  9,  8],
                 [5,  0, 10, 10,  9],
                 [9, 10,  0,  8,  7],
                 [9, 10,  8,  0,  3],
                 [8,  9,  7,  3,  0]]
        ids1 = list('abcde')
        self.dm1 = DistanceMatrix(data1, ids1)
        self.tree1 = TreeNode.read([
            "((b:3.0,(c:4.0,(e:1.0,d:2.0):2.0):3.0):2.0)a;"
        ])

        self.init1 = TreeNode.read(["((b,c))a;"])
        self.attach1 = TreeNode.read(["()d;"])

        self.init1.name = 0
        self.init1.children[0].children[0].name = 1
        self.init1.children[0].children[1].name = 2
        self.attach1.name = 3

        # Example 2
        # This example was pulled from the Phylip manual:
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
        self.tree2 = TreeNode.read([
            "((Mouse:0.7689100000000002,(Gibbon:0.37287500000000007,"
            "(Orang:0.2705916666666668,(Gorilla:0.1525416666666667,"
            "(Human:0.12403749999999991,Chimp:0.14516249999999997)"
            ":0.04120833333333351):0.0393958333333333):0.032630555555555496)"
            ":0.42026874999999997):0.9176899999999998)Bovine;"
        ])

        # Example 3
        data3 = [[0, 5, 4, 7, 6, 8],
                 [5, 0, 7, 10, 9, 11],
                 [4, 7, 0, 7, 6, 8],
                 [7, 10, 7, 0, 5, 8],
                 [6, 9, 6, 5, 0, 8],
                 [8, 11, 8, 8, 8, 0]]
        ids3 = map(str, range(6))
        self.dm3 = DistanceMatrix(data3, ids3)
        self.expected3_str = ()
        self.tree3 = TreeNode.read([
            "((1:4.0,(2:2.0,(5:4.75,(4:2.0,3:3.0):0.75):1.25):1.0):1.0)0;"
        ])

        # Example 4
        # This dm can yield negative branch lengths for OLS-based estimation.
        data4 = [[0,  5,  9,  9,  800],
                 [5,  0, 10, 10,  9],
                 [9, 10,  0,  8,  7],
                 [9, 10,  8,  0,  3],
                 [800,  9,  7,  3,  0]]
        ids4 = list('abcde')
        self.dm4 = DistanceMatrix(data4, ids4)

        # Example 5
        self.tree5 = TreeNode.read([
            "((((e:1.0,d:2.0):2.0,c:4.0):3.0,b:3.0):2.0)a;"
        ])
        for node in self.tree5.traverse(include_self=True):
            if (name := node.name) is not None:
                node.name = "abcde".index(name)

        # Example 6
        self.tree6 = TreeNode.read(["(((b,d),(e,c)))a;"])
        for node in self.tree6.traverse(include_self=True):
            if (name := node.name) is not None:
                node.name = "abcde".index(name)

    def test_gme_dm1(self):
        obs = gme(self.dm1)
        self.assertAlmostEqual(obs.compare_cophenet(self.tree1), 0.0, places=10)

    def test_gme_dm2(self):
        obs = gme(self.dm2)
        self.assertAlmostEqual(obs.compare_cophenet(self.tree2), 0.0, places=10)

    def test_gme_dm3(self):
        obs = gme(self.dm3)
        self.assertAlmostEqual(obs.compare_cophenet(self.tree3), 0.0, places=10)

    def test_gme_zero_branch_length(self):
        # OLS-based edge estimation can produce negative branch 
        # lengths when some dm values are much larger than
        # others, analogous to negative branch lengths produced by nj
        tree = gme(self.dm4)
        self.assertLess(tree.find("b").length, 0)
        self.assertLess(tree.find("c").length, 0)
        self.assertLess(tree.find("d").length, 0)
        self.assertGreater(tree.find("e").length, 0)

    def test_gme_error(self):
        data = [[0, 3],
                [3, 0]]
        dm = DistanceMatrix(data, list("ab"))
        self.assertRaises(ValueError, gme, dm)

    def test_average_distance_k(self):
        (stem,) = self.init1.children
        obs = _average_distance_k(self.attach1, stem, self.dm1.data)
        self.assertAlmostEqual(obs, 9.0, places=10)

    def test_average_distance_k_upper(self):
        (stem,) = self.init1.children
        obs = _average_distance_k_upper(self.attach1, stem, self.dm1.data)
        self.assertAlmostEqual(obs , 9.0, places=10)

    def test_lower_subtree_list(self):
        ordered = list(self.init1.postorder(include_self=False))
        obs = _lower_subtree_list(ordered, self.attach1, self.dm1.data)
        exp = [10.0, 8.0, 9.0]
        self.assertListEqual(obs, exp)

    def test_upper_subtree_list(self):
        ordered = list(self.init1.postorder(include_self=False))
        obs = _upper_subtree_list(ordered, self.attach1, 3, self.dm1)
        exp = [8.5, 9.5, 9.0]
        self.assertListEqual(obs, exp)

    def test_edge_attachment_length(self):
        ordered = list(self.init1.postorder(include_self=False))
        dm = self.dm1.data
        child = self.init1.children[0].children[1]  # taxon "c"
        lowerlist = _lower_subtree_list(ordered, self.attach1, dm)
        upperlist = _upper_subtree_list(ordered, self.attach1, 3, dm)
        adm = _average_distance_matrix(self.init1, dm)
        obs = _edge_attachment_length(child, lowerlist, upperlist, ordered, 3, adm)
        self.assertEqual(obs, -1.5)

    def test_average_distance(self):
        node1 = self.tree5.find(1)
        node2 = self.tree5.find(3).parent
        obs = _average_distance(node1, node2, self.dm1.data)
        self.assertAlmostEqual(obs, 9.5, places=10)

    def test_tip_or_root(self):
        node_internal = self.tree5.find(3).parent
        node_leaf = self.tree5.find(1)
        root = self.tree5
        self.assertEqual(len(_tip_or_root(node_internal)), 2)
        self.assertEqual(_tip_or_root(node_leaf)[0], node_leaf.name)
        self.assertEqual(_tip_or_root(root)[0], root.name)

    def test_average_distance_upper(self):
        # computed manually
        dm = np.array([[0, 0.02, 0.18, 0.34, 0.55],
                       [0.02, 0, 0.19, 0.35, 0.55],
                       [0.18, 0.19, 0, 0.34, 0.54],
                       [0.34, 0.35, 0.34, 0, 0.62],
                       [0.55, 0.55, 0.54, 0.62, 0]])
        tree = TreeNode.read(["((3,(0,(2,1))))4;"])
        for node in tree.traverse(include_self=True):
            if node.name is not None:
                node.name = int(node.name)
        node1 = tree.find(2).parent
        node2 = tree.find(0).parent.parent
        obs = _average_distance_upper(node1, node2, dm)
        self.assertAlmostEqual(obs, 0.545, places=10)

    def test_subtree_count(self):
        tree = self.tree5
        internal_node = tree.find(3).parent.parent
        leaf = tree.find(3)
        self.assertEqual(_subtree_count(internal_node), 3)
        self.assertEqual(_subtree_count(leaf), 1)
        self.assertEqual(_subtree_count(tree), 1)

    def test_average_subtree_distance(self):
        # computed manually
        tree = self.tree6
        a = tree.find(4).parent
        b = tree.find(1)
        a1 = tree.find(4)
        a2 = tree.find(2)
        obs = _average_subtree_distance(a, b, a1, a2, self.dm1)
        self.assertAlmostEqual(obs, 9.5, places=10)

    def test_average_distance_matrix_trivial(self):
        # In this case, the average distance matrix is equivalent to
        # the original distance matrix
        dm = np.array([[0, 3, 2],
                       [3, 0, 3],
                       [2, 3, 0]], dtype=float)
        tree = TreeNode.read(["((2,1))0;"])
        for node in tree.traverse(include_self=True):
            if node.name is not None:
                node.name = int(node.name)
        obs = _average_distance_matrix(tree, dm)
        np.testing.assert_array_equal(obs, dm)

    def test_average_distance_matrix(self):
        # computed manually
        obs = _average_distance_matrix(self.tree6, self.dm1)
        exp = np.array([[0.0, 10.0, 8.0, 9.0, 10.0, 9.5, 5.0],
                        [10.0, 0.0, 6.666666666666667, 3.0, 8.0, 5.5, 9.0],
                        [8.0, 6.666666666666667, 0.0, 6.0, 9.0, 7.5, 7.0],
                        [9.0, 3.0, 6.0, 0.0, 7.0, 6.666666666666667, 8.0],
                        [10.0, 8.0, 9.0, 7.0, 0.0, 9.0, 9.0],
                        [9.5, 5.5, 7.5, 6.666666666666667, 9.0, 0.0, 8.5],
                        [5.0, 9.0, 7.0, 8.0, 9.0, 8.5, 0.0]])
        np.testing.assert_array_equal(obs, exp)

    def test_edge_estimation(self):
        obs = TreeNode.read(["((c,b))a;"])
        exp = TreeNode.read(["((c:1.0,b:2.0):1.0)a;"])
        for t in (obs, exp):
            for node in t.traverse(include_self=True):
                if node.name is not None:
                    node.name = "abc".index(node.name)
        dm = np.array([[0, 3, 2],
                       [3, 0, 3],
                       [2, 3, 0]], dtype=float)
        _ols_edge(obs, dm)
        self.assertAlmostEqual(obs.compare_cophenet(
            exp), 0.0, places=10)


if __name__ == "__main__":
    main()
