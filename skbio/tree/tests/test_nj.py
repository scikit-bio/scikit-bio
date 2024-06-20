# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

from skbio import DistanceMatrix, TreeNode, nj
from skbio.tree._nj import (
    _compute_q, _compute_collapsed_dm, _lowest_index,
    _pair_members_to_new_node, nni, _perform_swap,
    _average_distance, _tip_or_root,
    _average_distance_upper, _subtree_count,
    _swap_length, _swap_heap, _average_subtree_distance,
    _average_distance_matrix, _edge_estimation)


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
        # For nni testing an arbitrary tree is given alongside the distance
        # matrix. Tree topologies are equivalent to that of the unrooted tree
        # of the above newick string.
        self.pre1_nni_str = ("(((b,d),(e,c)))a;")
        self.pre1_nni_TreeNode = TreeNode.read(
                io.StringIO(self.pre1_nni_str))
        self.post1_nni_str = ("((((e:1.0,d:2.0):2.0,c:4.0):3.0,b:3.0):2.0)a;")
        self.post1_nni_TreeNode = TreeNode.read(
                io.StringIO(self.post1_nni_str))

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
        self.pre2_nni_str = ("(((Mouse,Gorilla),(Gibbon,(Bovine,(Orang"
                             ",Chimp)))))Human;")
        self.pre2_nni_TreeNode = TreeNode.read(
                io.StringIO(self.pre2_nni_str))
        self.post2_nni_str = ("((((((Bovine:0.9117125,Mouse:0.7748875):0.42773"
                              "33,Gibbon:0.3504666):0.0408666,Orang:0.2809083)"
                              ":0.0345694,Gorilla:0.1475249):0.0414812,Chimp:0"
                              ".1470600):0.1221399)Human;")
        self.post2_nni_TreeNode = TreeNode.read(
                io.StringIO(self.post2_nni_str))

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
        self.pre3_nni_str = ("((1,(((5,2),4),3)))0;")
        self.pre3_nni_TreeNode = TreeNode.read(
                io.StringIO(self.pre3_nni_str))
        self.post3_nni_str = ("((1:4.0,((5:4.75,(4:2.0,3:3.0):0.75):1.25"
                              ",2:2.0):1.0):1.0)0;")
        self.post3_nni_TreeNode = TreeNode.read(
                io.StringIO(self.post3_nni_str))

        # this dm can yield negative branch lengths for both nj and nni
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
        # precision error on ARM: 1.6653345369377348e-16 != 0.0
        self.assertAlmostEqual(actual_TreeNode.compare_tip_distances(
            self.expected1_TreeNode), 0.0, places=10)

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

    def test_nni_dm1(self):
        self.assertEqual(nj(self.dm1, result_constructor=str),
                         self.expected1_str)
        actual_TreeNode = nni(self.pre1_nni_TreeNode, self.dm1, inplace=False)
        self.assertAlmostEqual(actual_TreeNode.compare_tip_distances(
            self.post1_nni_TreeNode), 0.0, places=10)

    def test_nni_dm2(self):
        # Resulting tree topology is equivalent to result from nj, however,
        # resulting edge lengths are almost equal to 2 places.
        actual_TreeNode = nni(self.pre2_nni_TreeNode, self.dm2, inplace=False)
        self.assertAlmostEqual(actual_TreeNode.compare_tip_distances(
            self.post2_nni_TreeNode), 0.0)

    def test_nni_dm3(self):
        actual_TreeNode = nni(self.pre3_nni_TreeNode, self.dm3, inplace=False)
        self.assertAlmostEqual(actual_TreeNode.compare_tip_distances(
            self.post3_nni_TreeNode), 0.0)
        
    def test_nni_trivial(self):
        # No swaps are performed, but edge lengths are assigned.
        data = [[0, 3, 2],
                [3, 0, 3],
                [2, 3, 0]]
        dm = DistanceMatrix(data, list('abc'))
        pre_str = "((c,b))a;"
        pre_TreeNode = TreeNode.read(
                io.StringIO(pre_str))
        expected_str = "((c:1.0,b:2.0):1.0)a;"
        expected_TreeNode = TreeNode.read(
                io.StringIO(expected_str))
        self.assertEqual(str(nni(pre_TreeNode, dm, inplace=False)),
                         str(expected_TreeNode))

    def test_nni_binary_flag(self):
        data = [[0, 3],
                [3, 0]]
        dm = DistanceMatrix(data, list('ab'))
        pre_str = "((b))a;"
        pre_TreeNode = TreeNode.read(io.StringIO(pre_str))
        msg = "Could not perform NNI. Tree needs to be a binary tree."
        with self.assertRaises(TypeError) as cm:
            nni(pre_TreeNode, dm)
        self.assertEqual(str(cm.exception), msg)

    def test_nni_leaf_root_flag(self):
        pre_str = "((b,d),(e,c))a;"
        pre_TreeNode = TreeNode.read(io.StringIO(pre_str))
        msg = "Could not perform NNI. Tree needs to be rooted at a leaf node."
        with self.assertRaises(TypeError) as cm:
            nni(pre_TreeNode, self.dm1)
        self.assertEqual(str(cm.exception), msg)

    def test_perform_swap(self):
        # Swapping the leaf nodes a tree without edge lengths.
        pre_str = "(((b,d),(e,c)))a;"
        actual_TreeNode = TreeNode.read(
            io.StringIO(pre_str))
        node1 = actual_TreeNode.find('b')
        node2 = actual_TreeNode.find('c')
        expected_str = "(((d,c),(e,b)))a;"
        expected_TreeNode = TreeNode.read(
            io.StringIO(expected_str))
        _perform_swap(node1, node2)
        self.assertEqual(str(actual_TreeNode),
                         str(expected_TreeNode))

    def test_average_distance(self):
        expected_str = ("((((e:1.0,d:2.0):2.0,c:4.0):3.0,b:3.0):2.0)a;")
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        node1 = expected_TreeNode.find('b')
        node2 = expected_TreeNode.find('d').parent
        self.assertAlmostEqual(_average_distance(node1, node2, self.dm1),
                               9.5, places=10)

    def test_tip_or_root(self):
        expected_str = ("((((e:1.0,d:2.0):2.0,c:4.0):3.0,b:3.0):2.0)a;")
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        node_internal = expected_TreeNode.find('d').parent
        node_leaf = expected_TreeNode.find('b')
        root = expected_TreeNode.root()
        self.assertEqual(len(_tip_or_root(node_internal)), 2)
        self.assertEqual(str(_tip_or_root(node_leaf)[0]),
                             str(node_leaf.name))
        self.assertEqual(str(_tip_or_root(root)[0]),
                             str(root.name))

    def test_average_distance_upper(self):
        # computed manually
        data = [[0, 0.02, 0.18, 0.34, 0.55],
                [0.02, 0, 0.19, 0.35, 0.55],
                [0.18, 0.19, 0, 0.34, 0.54],
                [0.34, 0.35, 0.34, 0, 0.62],
                [0.55, 0.55, 0.54, 0.62, 0]]
        ids = ['human','monkey','pig','rat','chicken']
        dm = DistanceMatrix(data, ids)
        expected_str = "((rat,(human,(pig,monkey))))chicken;"
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        node1 = expected_TreeNode.find('pig').parent
        node2 = expected_TreeNode.find('human').parent.parent
        self.assertAlmostEqual(_average_distance_upper(node1, node2, dm), 0.545, places=10)

    def test_subtree_count(self):
        expected_str = ("((((e:1.0,d:2.0):2.0,c:4.0):3.0,b:3.0):2.0)a;")
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        internal_node = expected_TreeNode.find('d').parent.parent
        leaf = expected_TreeNode.find('d')
        root = expected_TreeNode.root()
        self.assertEqual(_subtree_count(internal_node), 3)
        self.assertEqual(_subtree_count(leaf), 1)
        self.assertEqual(_subtree_count(root), 1)

    def test_swap_length(self):
        # results in a positive integer
        # computed manually
        expected_str = ("(((b,d),(e,c)))a;")
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        adm = _average_distance_matrix(expected_TreeNode, self.dm1)
        self.assertAlmostEqual(_swap_length(
            2, 1, 1, 1, 6, 3, 0, 1, adm), 2.5, places=10)

    def test_swap_heap(self):
        # swap length is stored into the maxheap as a negative integer
        expected_str = ("(((b,d),(e,c)))a;")
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        adm = _average_distance_matrix(expected_TreeNode, self.dm1)
        self.assertAlmostEqual(_swap_heap(expected_TreeNode, adm)[0][0],
                               -2.0, places=10)

    def test_average_subtree_distance(self):
        # computed manually
        expected_str = ("(((b,d),(e,c)))a;")
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        a = expected_TreeNode.find('e').parent
        b = expected_TreeNode.find('b')
        a1 = expected_TreeNode.find('e')
        a2 = expected_TreeNode.find('c')
        self.assertAlmostEqual(_average_subtree_distance(a, b, a1, a2, self.dm1),
                               9.5, places=10)

    def test_average_distance_matrix_trivial(self):
        # In this case, the average distance matrix is equivalent to
        # the original distance matrix
        data = [[0, 3, 2],
                [3, 0, 3],
                [2, 3, 0]]
        ids = list('abc')
        dm = DistanceMatrix(data, ids)
        expected_str = "((c,b))a;"
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        index = [0, 1, 2]
        actual_adm = _average_distance_matrix(expected_TreeNode, dm)
        for i in index:
            for j in index:
                if j < i:
                    self.assertEqual(dm[i][j], actual_adm[i][j])
                    self.assertEqual(dm[j][i], actual_adm[j][i])

    def test_average_distance_matrix(self):
        # computed manually
        expected_str = ("(((b,d),(e,c)))a;")
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        expected_adm = [[0.0, 10.0, 8.0, 9.0, 10.0, 9.5, 5.0],
                        [10.0, 0.0, 6.666666666666667, 3.0, 8.0, 5.5, 9.0],
                        [8.0, 6.666666666666667, 0.0, 6.0, 9.0, 7.5, 7.0],
                        [9.0, 3.0, 6.0, 0.0, 7.0, 6.666666666666667, 8.0],
                        [10.0, 8.0, 9.0, 7.0, 0.0, 9.0, 9.0],
                        [9.5, 5.5, 7.5, 6.666666666666667, 9.0, 0.0, 8.5],
                        [5.0, 9.0, 7.0, 8.0, 9.0, 8.5, 0.0]]
        actual_adm = _average_distance_matrix(expected_TreeNode, self.dm1)
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
        pre_estimation_str = "((c,b))a;"
        expected_str = "((c:1.0,b:2.0):1.0)a;"
        actual_TreeNode = TreeNode.read(io.StringIO(pre_estimation_str))
        _edge_estimation(actual_TreeNode, dm)
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        self.assertAlmostEqual(actual_TreeNode.compare_tip_distances(
            expected_TreeNode), 1, places=10)


if __name__ == "__main__":
    main()
