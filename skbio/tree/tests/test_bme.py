# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

from skbio import DistanceMatrix, TreeNode
from skbio.tree._bme import (
    _balanced_lower, _balanced_upper,
    _balanced_attach_length, _balanced_average_matrix,
    _bal_ols_edge, bme)


class BmeTests(TestCase):

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
        self.expected1_str = ("((b:3.0,(c:4.0,(e:1.25,d:1.75):2.0):3.0):2.0)a;")
        self.expected1_TreeNode = TreeNode.read(
                io.StringIO(self.expected1_str))
        self.bme_starting1_str = ("((b,c))a;")
        self.bme_starting1_TreeNode = TreeNode.read(
                io.StringIO(self.bme_starting1_str))
        self.bme_attach_node = TreeNode.read(["()d;"])

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
        self.expected2_str = ("((Mouse:0.7580687500000001,(Gibbon:0.38040000000000007,(Orang:0.2467625,"
                              "(Gorilla:0.16609999999999997,(Human:0.12582499999999996,Chimp:0.143375):"
                              "0.027650000000000063):0.06488749999999999):0.03341250000000007):"
                              "0.4108687500000001):0.92853125)Bovine;")
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
        self.expected3_str = ("((1:4.0,(2:2.0,(4:2.875,(5:4.25,3:3.75):-0.375):1.5):1.0):1.0)0;")
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

    def test_bme_dm1(self):
        actual_TreeNode = bme(self.dm1)
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            self.expected1_TreeNode), 0.0)

    def test_bme_dm2(self):
        actual_TreeNode = bme(self.dm2)
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            self.expected2_TreeNode), 0.0)

    def test_bme_dm3(self):
        actual_TreeNode = bme(self.dm3)
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            self.expected3_TreeNode), 0.0)

    def test_bme_zero_branch_length(self):
        # OLS-based edge estimation can produce negative branch 
        # lengths when some dm values are much larger than
        # others, analogous to negative branch lengths produced by nj
        tree = bme(self.dm4)
        self.assertTrue(tree.find('b').length < 0)
        self.assertTrue(tree.find('c').length < 0)
        self.assertTrue(tree.find('d').length > 0)
        self.assertTrue(tree.find('e').length > 0)

    def test_bme_error(self):
        data = [[0, 3],
                [3, 0]]
        dm = DistanceMatrix(data, list('ab'))
        self.assertRaises(ValueError, bme, dm)

    def test_balanced_lower(self):
        ordered = list(self.bme_starting1_TreeNode.postorder(include_self=False))
        expected_list = [10.0, 8.0, 9.0]
        self.assertEqual(_balanced_lower(ordered, self.bme_attach_node, self.dm1), expected_list)

    def test_balanced_upper(self):
        ordered = list(self.bme_starting1_TreeNode.postorder(include_self=False))
        expected_list = [8.5, 9.5, 9.0]
        self.assertEqual(_balanced_upper(ordered, self.bme_attach_node, self.dm1), expected_list)

    def test_balanced_attach_length(self):
        ordered = list(self.bme_starting1_TreeNode.postorder(include_self=False))
        child = self.bme_starting1_TreeNode.find('c')
        lowerlist = _balanced_lower(ordered, self.bme_attach_node, self.dm1)
        upperlist = _balanced_upper(ordered, self.bme_attach_node, self.dm1)
        adm = _balanced_average_matrix(self.bme_starting1_TreeNode, self.dm1)
        self.assertEqual(_balanced_attach_length(child, lowerlist, upperlist, ordered, adm), -1.5)

    def test_balanced_average_matrix_trivial(self):
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
        actual_adm = _balanced_average_matrix(expected_TreeNode, dm)
        for i in index:
            for j in index:
                if j < i:
                    self.assertEqual(dm[i][j], actual_adm[i][j])
                    self.assertEqual(dm[j][i], actual_adm[j][i])

    def test_balanced_average_matrix(self):
        # computed manually
        expected_str = ("(((b,d),(e,c)))a;")
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        expected_adm = [[0.0, 10.0, 4.75, 9.0, 10.0, 9.5, 5.0],
                        [10.0, 0.0, 2.75, 3.0, 8.0, 5.5, 9.0],
                        [4.75, 2.75, 0.0, 6.0, 9.0, 7.5, 7.0],
                        [9.0, 3.0, 6.0, 0.0, 7.0, 3.0, 8.0],
                        [10.0, 8.0, 9.0, 7.0, 0.0, 4.5, 9.0],
                        [9.5, 5.5, 7.5, 3.0, 4.5, 0.0, 8.5],
                        [5.0, 9.0, 7.0, 8.0, 9.0, 8.5, 0.0]]
        actual_adm = _balanced_average_matrix(expected_TreeNode, self.dm1)
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
        _bal_ols_edge(actual_TreeNode, dm)
        expected_TreeNode = TreeNode.read(io.StringIO(expected_str))
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            expected_TreeNode), 0.0, places=10)

if __name__ == "__main__":
    main()
