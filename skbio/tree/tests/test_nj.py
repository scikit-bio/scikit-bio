# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

import numpy.testing as npt

from skbio import DistanceMatrix, TreeNode
from skbio.tree._nj import nj, _tree_from_linkmat


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
        # https://phylipweb.github.io/phylip/doc/neighbor.html
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

        # this dm can yield negative branch lengths for nj
        data4 = [[0,  5,  9,  9,  800],
                 [5,  0, 10, 10,  9],
                 [9, 10,  0,  8,  7],
                 [9, 10,  8,  0,  3],
                 [800,  9,  7,  3,  0]]
        ids4 = list('abcde')
        self.dm4 = DistanceMatrix(data4, ids4)

    def test_nj_dm1(self):
        actual_TreeNode = nj(self.dm1)
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            self.expected1_TreeNode), 0.0)

    def test_nj_dm2(self):
        actual_TreeNode = nj(self.dm2)
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            self.expected2_TreeNode), 0.0)

    def test_nj_dm3(self):
        actual_TreeNode = nj(self.dm3)
        self.assertAlmostEqual(actual_TreeNode.compare_cophenet(
            self.expected3_TreeNode), 0.0)

    def test_nj_inplace(self):
        dm = self.dm3.copy()
        obs = nj(dm)
        exp = self.expected3_TreeNode
        self.assertAlmostEqual(obs.compare_cophenet(exp), 0.0)
        npt.assert_almost_equal(dm.data, self.dm3.data)

        obs = nj(dm, inplace=True)
        self.assertAlmostEqual(obs.compare_cophenet(exp), 0.0)
        with self.assertRaises(AssertionError):
            npt.assert_almost_equal(dm.data, self.dm3.data)

    def test_nj_zero_branch_length(self):
        # no nodes have negative branch length when we disallow negative
        # branch length. self is excluded as branch length is None
        tree = nj(self.dm4)
        for n in tree.postorder(include_self=False):
            self.assertTrue(n.length >= 0)
        # only tips associated with the large distance in the input
        # have positive branch lengths when we allow negative branch
        # length
        tree = nj(self.dm4, neg_as_zero=False)
        self.assertTrue(tree.find('a').length > 0)
        self.assertTrue(tree.find('b').length < 0)
        self.assertTrue(tree.find('c').length < 0)
        self.assertTrue(tree.find('d').length < 0)
        self.assertTrue(tree.find('e').length > 0)

        # deprecated functionality
        t2 = nj(self.dm4, disallow_negative_branch_length=False)
        self.assertAlmostEqual(tree.compare_cophenet(t2), 0.0)

    def test_nj_trivial(self):
        data = [[0, 3, 2],
                [3, 0, 3],
                [2, 3, 0]]
        dm = DistanceMatrix(data, list('abc'))
        exp = TreeNode.read(["(b:2.000000, a:1.000000, c:1.000000);"])
        self.assertAlmostEqual(nj(dm).compare_cophenet(exp), 0.0)

        # deprecated functionality
        self.assertEqual(nj(dm, result_constructor=str), "(a:1.0,b:2.0,c:1.0);\n")

    def test_nj_error(self):
        data = [[0, 3],
                [3, 0]]
        dm = DistanceMatrix(data, list('ab'))
        self.assertRaises(ValueError, nj, dm)

    def test_tree_from_linkmat(self):
        taxa = ['human', 'chimp', 'monkey', 'pig', 'mouse', 'rat', 'chicken']
        lm = [
            [5, 4, 0.065041, 0.060976],   # rat, mouse
            [6, 7, 0.457317, 0.229675],   # chicken, rat-mouse
            [3, 8, 0.161924, 0.050474],   # pig, chicken-rat-mouse
            [2, 9, 0.018293, 0.162602],   # monkey, pig-...-mouse
            [10, 1, 0.006098, 0.001016],  # chimp, monkey-...-mouse
            [0, 11, 0.003049, 0.0],       # human, chimp-...-mouse
        ]
        obs = str(_tree_from_linkmat(lm, taxa, rooted=True))
        exp = ("(human:0.003049,((monkey:0.018293,(pig:0.161924,(chicken:0.457317,"
               "(rat:0.065041,mouse:0.060976):0.229675):0.050474):0.162602):0.006098,"
               "chimp:0.001016):0.0);\n")
        self.assertEqual(str(obs), exp)
        obs = str(_tree_from_linkmat(lm, taxa, rooted=False))
        exp = ("(human:0.003049,(monkey:0.018293,(pig:0.161924,(chicken:0.457317,"
               "(rat:0.065041,mouse:0.060976):0.229675):0.050474):0.162602):0.006098,"
               "chimp:0.001016);\n")
        self.assertEqual(str(obs), exp)


if __name__ == "__main__":
    main()
