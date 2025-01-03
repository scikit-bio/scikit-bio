# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio import DistanceMatrix, TreeNode
from skbio.tree._me import (
    _check_tree,
    _to_treenode,
    _from_treenode,
    _allocate_arrays,
    _init_tree,
    _insert_taxon_treenode,
    _insert_taxon,
    _avgdist_matrix_naive,
    _avgdist_taxon_naive,
    _gme,
    gme,
    _bme,
    bme,
)
from skbio.tree._c_me import (
    _preorder,
    _postorder,
    _avgdist_matrix,
    _bal_avgdist_matrix,
    _avgdist_taxon,
    _bal_avgdist_taxon,
    _avgdist_d2_insert,
    _bal_avgdist_insert,
    _ols_lengths,
    _ols_lengths_d2,
    _bal_lengths,
    _ols_min_branch_d2,
    _bal_min_branch,
)


class MeTests(TestCase):

    def setUp(self):
        # Example 1
        # This newick string was confirmed against http://www.trex.uqam.ca/ which
        # generated the following (isomorphic) newick string:
        # (d:2.0000,e:1.0000,(c:4.0000,(a:2.0000,b:3.0000):3.0000):2.0000);
        #           /-b
        # -a-------|
        #          |          /-c
        #           \--------|
        #                    |          /-e
        #                     \--------|
        #                               \-d
        self.dm1 = np.array([
            [0,  5,  9,  9,  8],
            [5,  0, 10, 10,  9],
            [9, 10,  0,  8,  7],
            [9, 10,  8,  0,  3],
            [8,  9,  7,  3,  0],
        ], dtype=float)
        self.taxa1 = list("abcde")
        self.nwk1 = "(b:3.0,(c:4.0,(e:1.0,d:2.0):2.0):3.0)a:2.0;"
        self.tree1 = np.array([
            [1, 2, 0, 0, 4, 0, 0, 6],  # 0: a (root)
            [0, 1, 0, 2, 1, 1, 1, 0],  # 1: b
            [3, 4, 0, 1, 3, 1, 2, 5],  # 2: (c,(e,d))
            [0, 2, 2, 4, 1, 2, 3, 1],  # 3: c
            [5, 6, 2, 3, 2, 2, 4, 4],  # 4: (e,d)
            [0, 4, 4, 6, 1, 3, 5, 2],  # 5: e
            [0, 3, 4, 5, 1, 3, 6, 3],  # 6: d
        ])
        self.preodr1 = np.array([0, 1, 2, 3, 4, 5, 6])
        self.postodr1 = np.array([1, 3, 5, 6, 4, 2, 0])
        self.lens1 = np.array([2, 3, 3, 4, 2, 1, 2], dtype=float)

        # An unfinished tree for example 1 (missing e)
        #           /-b
        # -a-------|
        #          |          /-c
        #           \--------|
        #                     \-d
        self.nwk1m1 = "(b:3.0,(c:4.0,d:4.0):3.0)a:2.0;"
        self.tree1m1 = np.array([
            [1, 2, 0, 0, 3, 0, 0, 4],  # 0: a (root)
            [0, 1, 0, 2, 1, 1, 1, 0],  # 1: b
            [3, 4, 0, 1, 2, 1, 2, 3],  # 2: (c,d)
            [0, 2, 2, 4, 1, 2, 3, 1],  # 3: c
            [0, 3, 2, 3, 1, 2, 4, 2],  # 4: d
            [0, 0, 0, 0, 0, 0, 0, 0],  # 5: empty
            [0, 0, 0, 0, 0, 0, 0, 0],  # 6: empty
        ])
        self.preodr1m1 = np.array([0, 1, 2, 3, 4, 0, 0])
        self.postodr1m1 = np.array([1, 3, 4, 2, 0, 0, 0])

        # An alternative tree for example 1 (less ladder-like)
        #                     /-b
        #           /--------|
        #          |          \-d
        # -a-------|
        #          |          /-e
        #           \--------|
        #                     \-c
        self.nwk1v2 = "((b,d),(e,c))a;"
        self.tree1v2 = np.array([
            [1, 2, 0, 0, 4, 0, 0, 6],  # 0: a (root)
            [3, 4, 0, 2, 2, 1, 1, 2],  # 1: (b,d)
            [5, 6, 0, 1, 2, 1, 4, 5],  # 2: (e,c)
            [0, 1, 1, 4, 1, 2, 2, 0],  # 3: b
            [0, 3, 1, 3, 1, 2, 3, 1],  # 4: d
            [0, 4, 2, 6, 1, 2, 5, 3],  # 5: e
            [0, 2, 2, 5, 1, 2, 6, 4],  # 6: c
        ])
        self.preodr1v2 = np.array([0, 1, 3, 4, 2, 5, 6])
        self.postodr1v2 = np.array([3, 4, 1, 5, 6, 2, 0])

        # Example 2
        # This example was adopted from the Phylip manual:
        # https://phylipweb.github.io/phylip/doc/neighbor.html
        #           /-Mouse
        # -Bovine--|
        #          |          /-Gibbon
        #           \--------|
        #                    |          /-Orang
        #                     \--------|
        #                              |          /-Gorilla
        #                               \--------|
        #                                        |          /-Human
        #                                         \--------|
        #                                                   \-Chimp
        self.dm2 = np.array([
            [0.0000, 1.6866, 1.7198, 1.6606, 1.5243, 1.6043, 1.5905],
            [1.6866, 0.0000, 1.5232, 1.4841, 1.4465, 1.4389, 1.4629],
            [1.7198, 1.5232, 0.0000, 0.7115, 0.5958, 0.6179, 0.5583],
            [1.6606, 1.4841, 0.7115, 0.0000, 0.4631, 0.5061, 0.4710],
            [1.5243, 1.4465, 0.5958, 0.4631, 0.0000, 0.3484, 0.3083],
            [1.6043, 1.4389, 0.6179, 0.5061, 0.3484, 0.0000, 0.2692],
            [1.5905, 1.4629, 0.5583, 0.4710, 0.3083, 0.2692, 0.0000]
        ])
        self.taxa2 = [
            "Bovine", "Mouse", "Gibbon", "Orang", "Gorilla", "Chimp", "Human"
        ]
        self.nwk2 = (
            "(Mouse:0.76891,(Gibbon:0.372875,(Orang:0.2705917,(Gorilla:0.1525417,"
            "(Human:0.1240375,Chimp:0.1451625):0.0412083):0.03939583):0.0326305)"
            ":0.42026875)Bovine:0.91769;"
        )
        self.tree2 = np.array([
            [ 1,  2,  0,  0,  6,  0,  0, 10],
            [ 0,  1,  0,  2,  1,  1,  1,  0],
            [ 3,  4,  0,  1,  5,  1,  2,  9],
            [ 0,  2,  2,  4,  1,  2,  3,  1],
            [ 5,  6,  2,  3,  4,  2,  4,  8],
            [ 0,  3,  4,  6,  1,  3,  5,  2],
            [ 7,  8,  4,  5,  3,  3,  6,  7],
            [ 0,  4,  6,  8,  1,  4,  7,  3],
            [ 9, 10,  6,  7,  2,  4,  8,  6],
            [ 0,  6,  8, 10,  1,  5,  9,  4],
            [ 0,  5,  8,  9,  1,  5, 10,  5],
        ])
        self.preodr2 = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        self.postodr2 = np.array([1, 3, 5, 7, 9, 10, 8, 6, 4, 2, 0])
        self.lens2 = np.array([
            0.91769, 0.76891, 0.42026875, 0.372875, 0.0326305, 0.2705917, 0.03939583,
            0.1525417, 0.0412083, 0.1240375, 0.1451625,
        ])

        # Example 3
        # This example uses the same distance matrix as example 2, but with simplified
        # taxon names and different (less ladder-like) tree topology.
        #                               /-d
        #                     /--------|
        #           /--------|          \-c
        #          |         |
        # -a-------|          \-e
        #          |
        #          |          /-b
        #           \--------|
        #                     \-f
        self.dm3 = self.dm2.copy()
        self.taxa3 = list("abcdefg")
        self.nwk3 = "(((d,c),e),(b,f))a;"
        self.tree3 = np.array([
            [1, 2, 0, 0, 5, 0, 0, 8],  # 0, root (a)
            [3, 4, 0, 2, 3, 1, 1, 4],  # 1, ((d,c),e)
            [5, 6, 0, 1, 2, 1, 6, 7],  # 2, (b,f)
            [7, 8, 1, 4, 2, 2, 2, 2],  # 3, (d,c)
            [0, 4, 1, 3, 1, 2, 5, 3],  # 4, e
            [0, 1, 2, 6, 1, 2, 7, 5],  # 5, b
            [0, 5, 2, 5, 1, 2, 8, 6],  # 6, f
            [0, 3, 3, 8, 1, 3, 3, 0],  # 7, d
            [0, 2, 3, 7, 1, 3, 4, 1],  # 8, c
            [0, 0, 0, 0, 0, 0, 0, 0],  # reserved for g's link
            [0, 0, 0, 0, 0, 0, 0, 0],  # reserved for g
        ])
        self.preodr3 = np.array([0, 1, 3, 7, 8, 4, 2, 5, 6, 0, 0])
        self.postodr3 = np.array([7, 8, 3, 4, 1, 5, 6, 2, 0, 0, 0])

    def test_check_tree(self):
        """Check the integrity of a tree structure."""
        _check_tree(self.tree1, self.preodr1, self.postodr1)
        _check_tree(self.tree1m1, self.preodr1m1, self.postodr1m1)
        _check_tree(self.tree1v2, self.preodr1v2, self.postodr1v2)
        _check_tree(self.tree2, self.preodr2, self.postodr2)
        _check_tree(self.tree3, self.preodr3, self.postodr3)
        # TODO: add incorrect examples

    def test_to_treenode(self):
        """Convert an array-based tree structure into a TreeNode object."""
        # entire tree
        obs = _to_treenode(self.tree1, self.taxa1)
        exp = TreeNode.read([self.nwk1])
        self.assertEqual(obs.compare_rfd(exp), 0)

        # with branch lengths
        obs = _to_treenode(self.tree1, self.taxa1, self.lens1)
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        # unroot
        obs = _to_treenode(self.tree1, self.taxa1, self.lens1, unroot=True)
        exp = TreeNode.read(["(b:3.0,(c:4.0,(e:1.0,d:2.0):2.0):3.0,a:2.0);"])
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        # incomplete tree
        obs = _to_treenode(self.tree1m1, self.taxa1)
        exp = TreeNode.read([self.nwk1m1])
        self.assertEqual(obs.compare_rfd(exp), 0)

        # another example
        obs = _to_treenode(self.tree2, self.taxa2, self.lens2)
        exp = TreeNode.read([self.nwk2])
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

    def test_from_treenode(self):
        """Convert a TreeNode object into an array-based tree structure."""
        # a complete tree with branch lengths
        obs = _from_treenode(TreeNode.read([self.nwk1]), self.taxa1)
        _check_tree(*obs[:3])
        npt.assert_array_almost_equal(obs[3], self.lens1)

        # an incomplete tree without branch lengths
        obs = _from_treenode(TreeNode.read([self.nwk3]), self.taxa3)
        _check_tree(*obs[:3])
        self.assertTrue((obs[3] == 0).all())

    def test_preorder(self):
        """Perform preorder traversal."""
        n = self.tree1.shape[0]
        stack = np.full(n, 0)

        # entire tree
        _preorder(obs := np.full(n, 0), self.tree1, stack)
        npt.assert_array_equal(obs, self.preodr1)

        # given clade
        _preorder(obs := np.full(n, 0), self.tree1, stack, start=2)
        exp = np.array([2, 3, 4, 5, 6, 0, 0])
        npt.assert_array_equal(obs, exp)

        # alternative tree
        _preorder(obs := np.full(n, 0), self.tree1v2, stack)
        npt.assert_array_equal(obs, self.preodr1v2)

        # another example
        stack = np.full(n := self.tree3.shape[0], 0)
        _preorder(obs := np.full(n, 0), self.tree3, stack)
        npt.assert_array_equal(obs, self.preodr3)

    def test_postorder(self):
        """Perform postorder traversal."""
        n = self.tree1.shape[0]
        stack = np.full(n, 0)

        # entire tree
        _postorder(obs := np.full(n, 0), self.tree1, stack)
        npt.assert_array_equal(obs, self.postodr1)

        # given clade
        _postorder(obs := np.full(n, 0), self.tree1, stack, start=2)
        exp = np.array([3, 5, 6, 4, 2, 0, 0])
        npt.assert_array_equal(obs, exp)

        # alternative tree
        _postorder(obs := np.full(n, 0), self.tree1v2, stack)
        npt.assert_array_equal(obs, self.postodr1v2)

        # another example
        stack = np.full(n := self.tree3.shape[0], 0)
        _postorder(obs := np.full(n, 0), self.tree3, stack)
        npt.assert_array_equal(obs, self.postodr3)

    def test_allocate_arrays(self):
        """Allocate memory space for arrays."""
        n = 5
        tree, preodr, postodr, ads, adk, lens, stack = _allocate_arrays(n)
        self.assertTupleEqual(tree.shape, (5, 8))
        self.assertTupleEqual(preodr.shape, (5,))
        self.assertTupleEqual(postodr.shape, (5,))
        self.assertTupleEqual(ads.shape, (5, 2))
        self.assertTupleEqual(adk.shape, (5, 2))
        self.assertTupleEqual(lens.shape, (5,))
        self.assertTupleEqual(stack.shape, (5,))

        # make sure arrays are C-continuous
        for arr in tree, preodr, postodr, ads, adk, lens, stack:
            self.assertTrue(arr.flags.c_contiguous)

        _, _, _, ads, _, _, _ = _allocate_arrays(n, matrix=True)
        self.assertTupleEqual(ads.shape, (5, 5))
        self.assertTrue(ads.flags.c_contiguous)

    def test_init_tree(self):
        """Initialize triplet tree."""
        dm = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        n = 3
        tree, preodr, postodr, ads, _, _, _ = _allocate_arrays(n)
        _init_tree(dm, tree, preodr, postodr, ads)
        npt.assert_array_equal(tree, np.array([
            [1, 2, 0, 0, 2, 0, 0, 2],
            [0, 1, 0, 2, 1, 1, 1, 0],
            [0, 2, 0, 1, 1, 1, 2, 1],
        ]))
        npt.assert_array_equal(preodr, np.array([0, 1, 2]))
        npt.assert_array_equal(postodr, np.array([1, 2, 0]))
        npt.assert_array_equal(ads, np.array([
            [0, 0],
            [1, 3],
            [2, 3],
        ]))

        tree, preodr, postodr, ads, _, _, _ = _allocate_arrays(n, matrix=True)
        _init_tree(dm, tree, preodr, postodr, ads, matrix=True)
        np.fill_diagonal(ads, 0)
        npt.assert_array_equal(ads, np.array([
            [0, 1, 2],
            [1, 0, 3],
            [2, 3, 0],
        ]))

    def test_insert_taxon_treenode(self):
        """Insert a taxon between a target node and its parent."""
        # insert into a regular branch
        obs = TreeNode.read([self.nwk1m1])
        _insert_taxon_treenode("e", obs.find("d"), obs)
        exp = TreeNode.read([self.nwk1])
        self.assertEqual(obs.compare_rfd(exp), 0)

        # insert into the root branch
        obs = TreeNode.read([self.nwk1m1])
        _insert_taxon_treenode("e", obs, obs)
        exp = TreeNode.read(["((b,(c,d)),e)a;"])
        self.assertEqual(obs.compare_rfd(exp), 0)

    def test_insert_taxon(self):
        """Insert a taxon between a target node and its parent."""
        # insert taxon e (4) above d (4)
        obs = self.tree1m1.copy(), self.preodr1m1.copy(), self.postodr1m1.copy()
        _insert_taxon(4, 4, *obs)
        exp = np.array([
            [1, 2, 0, 0, 4, 0, 0, 6],
            [0, 1, 0, 2, 1, 1, 1, 0],
            [3, 5, 0, 1, 3, 1, 2, 5],
            [0, 2, 2, 5, 1, 2, 3, 1],
            [0, 3, 5, 6, 1, 3, 5, 2],
            [4, 6, 2, 3, 2, 2, 4, 4],
            [0, 4, 5, 4, 1, 3, 6, 3],
        ]), np.array([
            0, 1, 2, 3, 5, 4, 6
        ]), np.array([
            1, 3, 4, 6, 5, 2, 0
        ])
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # no depth
        obs = self.tree1m1.copy(), self.preodr1m1.copy(), self.postodr1m1.copy()
        _insert_taxon(4, 4, *obs)
        exp[0][:, 5] = [0, 1, 1, 2, 3, 2, 3]
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # insert taxon e (4) above b (1)
        obs = self.tree1m1.copy(), self.preodr1m1.copy(), self.postodr1m1.copy()
        _insert_taxon(4, 1, *obs)
        exp = np.array([
            [5, 2, 0, 0, 4, 0, 0, 6],
            [0, 1, 5, 6, 1, 2, 2, 0],
            [3, 4, 0, 5, 2, 1, 4, 5],
            [0, 2, 2, 4, 1, 2, 5, 3],
            [0, 3, 2, 3, 1, 2, 6, 4],
            [1, 6, 0, 2, 2, 1, 1, 2],
            [0, 4, 5, 1, 1, 2, 3, 1],
        ]), np.array([
            0, 5, 1, 6, 2, 3, 4
        ]), np.array([
            1, 6, 5, 3, 4, 2, 0
        ])
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # insert taxon e (4) above (c,d) (2)
        obs = self.tree1m1.copy(), self.preodr1m1.copy(), self.postodr1m1.copy()
        _insert_taxon(4, 2, *obs)
        exp = np.array([
            [1, 5, 0, 0, 4, 0, 0, 6],
            [0, 1, 0, 5, 1, 1, 1, 0],
            [3, 4, 5, 6, 2, 2, 3, 3],
            [0, 2, 2, 4, 1, 3, 4, 1],
            [0, 3, 2, 3, 1, 3, 5, 2],
            [2, 6, 0, 1, 3, 1, 2, 5],
            [0, 4, 5, 2, 1, 2, 6, 4],
        ]), np.array([
            0, 1, 5, 2, 3, 4, 6
        ]), np.array([
            1, 3, 4, 2, 6, 5, 0
        ])
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # insert taxon e (4) into the root branch (0)
        obs = self.tree1m1.copy(), self.preodr1m1.copy(), self.postodr1m1.copy()
        _insert_taxon(4, 0, *obs)
        exp = np.array([
            [5, 6, 0, 0, 4, 0, 0, 6],
            [0, 1, 5, 2, 1, 2, 2, 0],
            [3, 4, 5, 1, 2, 2, 3, 3],
            [0, 2, 2, 4, 1, 3, 4, 1],
            [0, 3, 2, 3, 1, 3, 5, 2],
            [1, 2, 0, 6, 3, 1, 1, 4],
            [0, 4, 0, 5, 1, 1, 6, 5],
        ]), np.array([
            0, 5, 1, 2, 3, 4, 6
        ]), np.array([
            1, 3, 4, 2, 5, 6, 0
        ])
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # another example; all possible insertions
        tree, taxa = self.tree3, self.taxa3
        m = tree[0, 4] + 1
        n = m * 2 - 3
        taxamap = {}
        for node in self.postodr3:
            if tree[node, 0] == 0:
                taxamap[node] = [taxa[tree[node, 1]]]
            else:
                taxamap[node] = taxamap[tree[node, 0]] + taxamap[tree[node, 1]]
        obj = TreeNode.read([self.nwk3])

        for i in range(n):
            obs = tree.copy(), self.preodr3.copy(), self.postodr3.copy()
            _insert_taxon(m, i, *obs)
            obs = _to_treenode(obs[0], taxa)

            exp = obj.copy()
            _insert_taxon_treenode("g", exp.lca(taxamap[i]), exp)

            self.assertEqual(obs.compare_rfd(exp), 0)

    def test_avgdist_matrix_naive(self):
        """Calculate an average distance matrix using a naive method."""
        # Test if the algorithm reproduces the manually calculated result.
        n = self.tree1.shape[0]
        obs = np.zeros((n, n), dtype=float)
        _avgdist_matrix_naive(obs, self.dm1, self.tree1, self.postodr1)
        exp = np.array([
            [ 0.   ,  5.   ,  8.667,  9.   ,  8.5  ,  8.   ,  9.   ],
            [ 5.   ,  0.   ,  9.667, 10.   ,  9.5  ,  9.   , 10.   ],
            [ 8.667,  9.667,  0.   ,  9.5  ,  9.   ,  8.5  ,  9.5  ],
            [ 9.   , 10.   ,  9.5  ,  0.   ,  7.5  ,  7.   ,  8.   ],
            [ 8.5  ,  9.5  ,  9.   ,  7.5  ,  0.   ,  8.   ,  9.   ],
            [ 8.   ,  9.   ,  8.5  ,  7.   ,  8.   ,  0.   ,  3.   ],
            [ 9.   , 10.   ,  9.5  ,  8.   ,  9.   ,  3.   ,  0.   ],
        ])
        npt.assert_array_equal(obs.round(3), exp)

        # incomplete tree
        obs = np.zeros((n, n), dtype=float)
        _avgdist_matrix_naive(obs, self.dm1, self.tree1m1, self.postodr1m1)
        exp = np.array([
            [ 0. ,  5. ,  9. ,  9. ,  9. ,  0. ,  0. ],
            [ 5. ,  0. , 10. , 10. , 10. ,  0. ,  0. ],
            [ 9. , 10. ,  0. ,  9.5,  9.5,  0. ,  0. ],
            [ 9. , 10. ,  9.5,  0. ,  8. ,  0. ,  0. ],
            [ 9. , 10. ,  9.5,  8. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
        ])
        npt.assert_array_almost_equal(obs, exp)

        # alternative tree
        obs = np.zeros((n, n), dtype=float)
        _avgdist_matrix_naive(obs, self.dm1, self.tree1v2, self.postodr1v2)
        exp = np.array([
            [ 0.  ,  7.  ,  8.5 ,  5.  ,  9.  ,  8.  ,  9.  ],
            [ 7.  ,  0.  ,  7.5 ,  8.  ,  6.67,  6.  ,  9.  ],
            [ 8.5 ,  7.5 ,  0.  ,  9.5 ,  5.5 ,  6.67,  9.  ],
            [ 5.  ,  8.  ,  9.5 ,  0.  , 10.  ,  9.  , 10.  ],
            [ 9.  ,  6.67,  5.5 , 10.  ,  0.  ,  3.  ,  8.  ],
            [ 8.  ,  6.  ,  6.67,  9.  ,  3.  ,  0.  ,  7.  ],
            [ 9.  ,  9.  ,  9.  , 10.  ,  8.  ,  7.  ,  0.  ],
        ])
        npt.assert_array_equal(obs.round(2), exp)

    def test__avgdist_matrix(self):
        """Calculate an average distance matrix."""
        # Test if the algorithm produces the same result as the native method does.
        dm, tree, preodr, postodr = self.dm1, self.tree1, self.preodr1, self.postodr1
        n = tree.shape[0]
        obs = np.zeros((n, n), dtype=float)
        _avgdist_matrix(obs, dm, tree, preodr, postodr)
        exp = np.zeros((n, n), dtype=float)
        _avgdist_matrix_naive(exp, dm, tree, postodr)
        npt.assert_array_almost_equal(obs, exp)

    def test_avgdist_matrix(self):
        """Calculate an average distance matrix."""
        # Test if the algorithm produces the same result as the native method does.
        dm, tree, preodr, postodr = self.dm1, self.tree1, self.preodr1, self.postodr1
        n = tree.shape[0]
        obs = np.zeros((n, n), dtype=float)
        _avgdist_matrix(obs, dm, tree, preodr, postodr)
        exp = np.zeros((n, n), dtype=float)
        _avgdist_matrix_naive(exp, dm, tree, postodr)
        npt.assert_array_almost_equal(obs, exp)

        # make sure all cells are populated
        _avgdist_matrix(obs := np.full((n, n), np.nan), dm, tree, preodr, postodr)
        np.fill_diagonal(obs, 0)
        _avgdist_matrix_naive(exp := np.zeros((n, n)), dm, tree, postodr)
        npt.assert_array_almost_equal(obs, exp)

        # incomplete tree
        tree, preodr, postodr = self.tree1m1, self.preodr1m1, self.postodr1m1
        _avgdist_matrix(obs := np.zeros((n, n)), dm, tree, preodr, postodr)
        _avgdist_matrix_naive(exp := np.zeros((n, n)), dm, tree, postodr)
        npt.assert_array_almost_equal(obs, exp)

        # example 1 v2
        tree, preodr, postodr = self.tree1v2, self.preodr1v2, self.postodr1v2
        _avgdist_matrix(obs := np.zeros((n, n)), dm, tree, preodr, postodr)
        _avgdist_matrix_naive(exp := np.zeros((n, n)), dm, tree, postodr)
        npt.assert_array_almost_equal(obs, exp)

        # example 2 (complete)
        dm, tree, preodr, postodr = self.dm2, self.tree2, self.preodr2, self.postodr2
        n = tree.shape[0]
        _avgdist_matrix(obs := np.zeros((n, n)), dm, tree, preodr, postodr)
        _avgdist_matrix_naive(exp := np.zeros((n, n)), dm, tree, postodr)
        npt.assert_array_almost_equal(obs, exp)

        # example 3 (incomplete)
        dm, tree, preodr, postodr = self.dm3, self.tree3, self.preodr3, self.postodr3
        n = tree.shape[0]
        _avgdist_matrix(obs := np.zeros((n, n)), dm, tree, preodr, postodr)
        _avgdist_matrix_naive(exp := np.zeros((n, n)), dm, tree, postodr)
        npt.assert_array_almost_equal(obs, exp)

    def test_bal_avgdist_matrix(self):
        """Calculate a balanced average distance matrix."""
        dm, tree, preodr, postodr = self.dm1, self.tree1, self.preodr1, self.postodr1
        n = tree.shape[0]
        _bal_avgdist_matrix(obs := np.zeros((n, n)), dm, tree, preodr, postodr)
        exp = np.array([
            [ 0.  ,  5.  ,  8.75,  9.  ,  8.5 ,  8.  ,  9.  ],
            [ 5.  ,  0.  ,  9.75, 10.  ,  9.5 ,  9.  , 10.  ],
            [ 8.75,  9.75,  0.  ,  9.5 ,  9.  ,  8.5 ,  9.5 ],
            [ 9.  , 10.  ,  9.5 ,  0.  ,  7.5 ,  7.  ,  8.  ],
            [ 8.5 ,  9.5 ,  9.  ,  7.5 ,  0.  ,  7.75,  8.75],
            [ 8.  ,  9.  ,  8.5 ,  7.  ,  7.75,  0.  ,  3.  ],
            [ 9.  , 10.  ,  9.5 ,  8.  ,  8.75,  3.  ,  0.  ],
        ])
        npt.assert_array_almost_equal(obs, exp)

        _bal_avgdist_matrix(obs := np.full((n, n), np.nan), dm, tree, preodr, postodr)
        np.fill_diagonal(obs, 0)
        npt.assert_array_almost_equal(obs, exp)

        dm, tree, preodr, postodr = self.dm3, self.tree3, self.preodr3, self.postodr3
        n = tree[0, 4] * 2 - 1
        _bal_avgdist_matrix(obs := np.zeros((n, n)), dm, tree, preodr, postodr)
        exp = np.array([
            [0.   , 1.607, 1.645, 1.69 , 1.524, 1.687, 1.604, 1.661, 1.72 ],
            [1.607, 0.   , 0.965, 1.362, 1.211, 1.475, 0.455, 1.328, 1.395],
            [1.645, 0.965, 0.   , 1.033, 0.897, 1.581, 1.03 , 0.995, 1.071],
            [1.69 , 1.362, 1.033, 0.   , 0.529, 1.504, 0.562, 0.895, 0.995],
            [1.524, 1.211, 0.897, 0.529, 0.   , 1.446, 0.348, 0.463, 0.596],
            [1.687, 1.475, 1.581, 1.504, 1.446, 0.   , 1.439, 1.484, 1.523],
            [1.604, 0.455, 1.03 , 0.562, 0.348, 1.439, 0.   , 0.506, 0.618],
            [1.661, 1.328, 0.995, 0.895, 0.463, 1.484, 0.506, 0.   , 0.712],
            [1.72 , 1.395, 1.071, 0.995, 0.596, 1.523, 0.618, 0.712, 0.   ],
        ])
        npt.assert_array_equal(obs.round(3), exp)

    def test_avgdist_taxon_naive(self):
        """Calculate taxon-to-subtree average distances using a naive method."""
        # Test if the algorithm reproduces the manually calculated result.
        n = self.tree1.shape[0]
        taxon = self.dm1.shape[0] - 1
        obs = np.zeros((n, 2), dtype=float)
        _avgdist_taxon_naive(obs, taxon, self.dm1, self.tree1m1, self.postodr1m1)
        exp = np.array([
            [6.333, 8.   ],
            [9.   , 6.   ],
            [5.   , 8.5  ],
            [7.   , 6.667],
            [3.   , 8.   ],
            [0.   , 0.   ],
            [0.   , 0.   ],
        ])
        npt.assert_array_equal(obs.round(3), exp)

    def test_avgdist_taxon(self):
        """Calculate taxon-to-subtree average distances."""
        # Test if the algorithm produces the same result as the naive method does.
        dm = self.dm1
        taxon = dm.shape[0] - 1
        tree, preodr, postodr = self.tree1m1, self.preodr1m1, self.postodr1m1
        n = tree.shape[0]
        _avgdist_taxon(obs := np.zeros((n, 2)), taxon, dm, tree, preodr, postodr)
        _avgdist_taxon_naive(exp := np.zeros((n, 2)), taxon, dm, tree, postodr)
        npt.assert_array_almost_equal(obs, exp)

    def test_bal_avgdist_taxon(self):
        """Calculate taxon-to-subtree balanced average distances."""
        dm = self.dm1
        taxon = dm.shape[0] - 1
        tree, preodr, postodr = self.tree1m1, self.preodr1m1, self.postodr1m1
        n = tree.shape[0]
        _bal_avgdist_taxon(obs := np.zeros((n, 2)), taxon, dm, tree, preodr, postodr)
        exp = np.array([
            [7.  , 8.  ],
            [9.  , 6.5 ],
            [5.  , 8.5 ],
            [7.  , 5.75],
            [3.  , 7.75],
            [0.  , 0.  ],
            [0.  , 0.  ],
        ])
        npt.assert_array_almost_equal(obs, exp)

        dm = self.dm3
        taxon = dm.shape[0] - 1
        tree, preodr, postodr = self.tree3, self.preodr3, self.postodr3
        n = tree.shape[0]
        _bal_avgdist_taxon(obs := np.zeros((n, 2)), taxon, dm, tree, preodr, postodr)
        exp = np.array([
            [0.639, 1.59 ],
            [0.411, 1.228],
            [0.866, 1.001],
            [0.515, 0.768],
            [0.308, 0.871],
            [1.463, 0.635],
            [0.269, 1.232],
            [0.471, 0.663],
            [0.558, 0.62 ],
            [0.   , 0.   ],
            [0.   , 0.   ],
        ])
        npt.assert_array_equal(obs.round(3), exp)

    def test_avgdist_d2_insert(self):
        """Update distant-2 subtree average distances after taxon insertion."""
        n = self.tree1.shape[0]
        taxon = self.dm1.shape[0] - 1
        # The following values were taken from the full-scale distance matrix.
        ad2 = np.array([
            [ 0. ,  0. ],  # 0: a (root)
            [ 5. , 10. ],  # 1: b
            [ 9. , 10. ],  # 2: (c,d)
            [ 9.5,  8. ],  # 3: c
            [ 9.5,  8. ],  # 4: d
            [ 0. ,  0. ],  # 5: empty
            [ 0. ,  0. ],  # 6: empty
        ])
        # Insert e as a sibling of d. This should recover tree1.
        target = 4
        _avgdist_taxon(
            adk := np.zeros((n, 2), dtype=float), taxon, self.dm1, self.tree1m1,
            self.preodr1m1, self.postodr1m1
        )
        _avgdist_d2_insert(
            obs := ad2.copy(), target, adk, self.tree1m1, self.preodr1m1
        )
        # The following values were taken from the full-scale distance matrix.
        exp = np.array([
            [0.   , 0.   ],  # 0: a (root)
            [5.   , 9.667],  # 1: b
            [8.667, 9.667],  # 2: (c,(e,d))
            [9.5  , 7.5  ],  # 3: c
            [9.   , 3.   ],  # 4: d
            [9.   , 7.5  ],  # 5: (e,d)
            [8.   , 3.   ],  # 6: e
        ])
        npt.assert_array_equal(obs.round(3), exp)

        # another example; all possible insertions
        dm, tree, preodr, postodr = self.dm3, self.tree3, self.preodr3, self.postodr3
        n = tree.shape[0]
        m = tree[0, 4] + 1
        ran_ = np.arange(n)

        # get distant-2 values from the full matrix
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, preodr, postodr)
        ad2 = np.ascontiguousarray(
            np.vstack([adm[ran_, tree[ran_, 2]], adm[ran_, tree[ran_, 3]]]).T)
        _avgdist_taxon(adk := np.zeros((n, 2)), m, dm, tree, preodr, postodr)

        for i in range(n - 2):
            # calculate distant-2 values using the algorithm
            _avgdist_d2_insert(obs := ad2.copy(), i, adk, tree, preodr)

            # insert taxon and calculate full matrix
            tree_, pre_, post_ = tree.copy(), preodr.copy(), postodr.copy()
            _insert_taxon(m, i, tree_, pre_, post_)
            _avgdist_matrix(adm := np.zeros((n, n)), dm, tree_, pre_, post_)
            exp = np.ascontiguousarray(
                np.vstack([adm[ran_, tree_[ran_, 2]], adm[ran_, tree_[ran_, 3]]]).T)

            npt.assert_array_almost_equal(obs, exp)

    def test_bal_avgdist_insert(self):
        """Update balanced average distance matrix after taxon insertion."""
        # Test if the algorithm produces the same result as calculated from the full
        # matrix after taxon insertion.
        tree, preodr, postodr = self.tree1m1, self.preodr1m1, self.postodr1m1
        powers = 2.0 ** -np.arange(5)
        stack = np.zeros(7, dtype=int)
        adm = np.array([
            [ 0. ,  5. ,  9. ,  9. ,  9. ,  0. ,  0. ],
            [ 5. ,  0. , 10. , 10. , 10. ,  0. ,  0. ],
            [ 9. , 10. ,  0. ,  9.5,  9.5,  0. ,  0. ],
            [ 9. , 10. ,  9.5,  0. ,  8. ,  0. ,  0. ],
            [ 9. , 10. ,  9.5,  8. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
        ])
        adk = np.array([
            [7.  , 8.  ],
            [9.  , 6.5 ],
            [5.  , 8.5 ],
            [7.  , 5.75],
            [3.  , 7.75],
            [0.  , 0.  ],
            [0.  , 0.  ],
        ])
        # Insert e as a sibling of d. This should recover tree1.
        target = 4
        _bal_avgdist_insert(
            obs := adm.copy(), target, adk, tree, preodr, postodr, powers, stack
        )
        exp = np.array([
            [ 0.  ,  5.  ,  8.75,  9.  ,  9.  ,  8.5 ,  8.  ],
            [ 5.  ,  0.  ,  9.75, 10.  , 10.  ,  9.5 ,  9.  ],
            [ 8.75,  9.75,  0.  ,  9.5 ,  9.5 ,  9.  ,  8.5 ],
            [ 9.  , 10.  ,  9.5 ,  0.  ,  8.  ,  7.5 ,  7.  ],
            [ 9.  , 10.  ,  9.5 ,  8.  ,  0.  ,  8.75,  3.  ],
            [ 8.5 ,  9.5 ,  9.  ,  7.5 ,  8.75,  0.  ,  7.75],
            [ 8.  ,  9.  ,  8.5 ,  7.  ,  3.  ,  7.75,  0.  ],
        ])
        npt.assert_array_almost_equal(obs, exp)

        # another example; all possible insertions
        dm, tree, preodr, postodr = self.dm3, self.tree3, self.preodr3, self.postodr3
        n = tree.shape[0]
        m = tree[0, 4] + 1
        ran_ = np.arange(n)
        _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree, preodr, postodr)
        _bal_avgdist_taxon(adk := np.zeros((n, 2)), m, dm, tree, preodr, postodr)
        powers = 2.0 ** -np.arange(m)
        stack = np.zeros(n, dtype=int)

        for i in range(n - 2):
            # update matrix using the algorithm
            _bal_avgdist_insert(
                obs := adm.copy(), i, adk, tree, preodr, postodr, powers, stack
            )

            # insert taxon and calculate full matrix
            tree_, pre_, post_ = tree.copy(), preodr.copy(), postodr.copy()
            _insert_taxon(m, i, tree_, pre_, post_)
            _bal_avgdist_matrix(exp := np.zeros((n, n)), dm, tree_, pre_, post_)

            npt.assert_array_almost_equal(obs, exp)

    def test_ols_lengths(self):
        """Calculate tree branch lengths using an OLS framework."""
        # Example 1: This should recover the same branch lengths as in the original
        # tree.
        dm, tree, preodr, postodr = self.dm1, self.tree1, self.preodr1, self.postodr1
        n = tree.shape[0]
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, preodr, postodr)
        _ols_lengths(obs := np.zeros(n), adm, tree, preodr)
        npt.assert_array_almost_equal(obs, self.lens1)

        # Example 2: The output is close but not precisely identical to those in the
        # original tree.
        dm, tree, preodr, postodr = self.dm2, self.tree2, self.preodr2, self.postodr2
        n = tree.shape[0]
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, preodr, postodr)
        _ols_lengths(obs := np.zeros(n), adm, tree, preodr)
        exp = np.array([
            0.91769, 0.76891, 0.42026875, 0.35793125, 0.04316597, 0.28054444,
            0.03137847, 0.15226875, 0.04148125, 0.12214, 0.14706])
        npt.assert_array_almost_equal(obs, exp)

    def test_ols_lengths_d2(self):
        """Calculate tree branch lengths using an OLS framework."""
        # Test if result matches that calculated from the full matrix.
        # example 1
        dm, tree, preodr, postodr = self.dm1, self.tree1, self.preodr1, self.postodr1
        n = tree.shape[0]
        ran_ = np.arange(n)
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, preodr, postodr)
        ad2 = np.ascontiguousarray(
            np.vstack([adm[ran_, tree[ran_, 2]], adm[ran_, tree[ran_, 3]]]).T)
        _ols_lengths_d2(obs := np.zeros(n), ad2, tree, preodr)
        _ols_lengths(exp := np.zeros(n), adm, tree, preodr)
        npt.assert_array_almost_equal(obs, exp)

        # example 2
        dm, tree, preodr, postodr = self.dm2, self.tree2, self.preodr2, self.postodr2
        n = tree.shape[0]
        ran_ = np.arange(n)
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, preodr, postodr)
        ad2 = np.ascontiguousarray(
            np.vstack([adm[ran_, tree[ran_, 2]], adm[ran_, tree[ran_, 3]]]).T)
        _ols_lengths_d2(obs := np.zeros(n), ad2, tree, preodr)
        _ols_lengths(exp := np.zeros(n), adm, tree, preodr)
        npt.assert_array_almost_equal(obs, exp)

    def test_bal_lengths(self):
        """Calculate tree branch lengths using a balanced framework."""
        # Example 1: In this simple example, OLS and balanced frameworks produce the
        # same branch lengths.
        dm, tree, preodr, postodr = self.dm1, self.tree1, self.preodr1, self.postodr1
        n = tree.shape[0]
        _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree, preodr, postodr)
        _bal_lengths(obs := np.zeros(n), adm, tree, preodr)
        npt.assert_array_almost_equal(obs, self.lens1)

        # Example 2: Also slightly different from the original branch lengths.
        dm, tree, preodr, postodr = self.dm2, self.tree2, self.preodr2, self.postodr2
        n = tree.shape[0]
        _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree, preodr, postodr)
        _bal_lengths(obs := np.zeros(n), adm, tree, preodr)
        exp = np.array([
            0.92853125, 0.75806875, 0.41086875, 0.36733125, 0.04648125, 0.28469375,
            0.02695625, 0.15393125, 0.03981875, 0.11678125, 0.15241875])
        npt.assert_array_almost_equal(obs, exp)

    def test_ols_min_branch_d2(self):
        """Find the branch with minimum length change using an OLS framework."""
        # Test if result matches ground truth.
        dm = self.dm1
        tree, preodr, postodr = self.tree1m1, self.preodr1m1, self.postodr1m1
        n = tree.shape[0]
        m = tree[0, 4] + 1
        ran_ = np.arange(n)
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, preodr, postodr)
        ad2 = np.ascontiguousarray(np.vstack([
            adm[ran_, tree[ran_, 2]], adm[ran_, tree[ran_, 3]]]).T)
        _avgdist_taxon(adk := np.zeros((n, 2)), m, dm, tree, preodr, postodr)
        res = _ols_min_branch_d2(obs := np.zeros(n), ad2, adk, tree, preodr)
        self.assertEqual(res, 4)
        exp = np.array([0, 0, -2, -2, -3, 0, 0], dtype=float)
        npt.assert_array_almost_equal(obs, exp)

        # Test if result matches that calculated from the entire tree.
        dm, tree, preodr, postodr = self.dm3, self.tree3, self.preodr3, self.postodr3
        n = tree.shape[0]
        m = tree[0, 4] + 1
        ran_ = np.arange(n)
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, preodr, postodr)
        ad2 = np.ascontiguousarray(np.vstack([
            adm[ran_, tree[ran_, 2]], adm[ran_, tree[ran_, 3]]]).T)
        _avgdist_taxon(adk := np.zeros((n, 2)), m, dm, tree, preodr, postodr)
        res = _ols_min_branch_d2(obs := np.zeros(n), ad2, adk, tree, preodr)
        # For each branch, insert taxon, calculate full matrix, then calculate and sum
        # all branch lengths. The difference between each sum and the sum by the root
        # branch is the length change value calculated by the algorithm. 
        exp = np.zeros(n)
        for i in range(n - 2):
            tree_, pre_, post_ = tree.copy(), preodr.copy(), postodr.copy()
            _insert_taxon(m, i, tree_, pre_, post_)
            _avgdist_matrix(adm := np.zeros((n, n)), dm, tree_, pre_, post_)
            _ols_lengths(lens := np.zeros(n), adm, tree_, pre_)
            exp[i] = lens.sum()
        exp[:n - 2] -= exp[0]

        npt.assert_array_almost_equal(obs, exp)
        self.assertEqual(res, exp[:n - 2].argmin())

    def test_bal_min_branch(self):
        """Find the branch with minimum length change using a balacned framework."""
        # Test if result matches ground truth.
        dm = self.dm1
        tree, preodr, postodr = self.tree1m1, self.preodr1m1, self.postodr1m1
        n = tree.shape[0]
        m = tree[0, 4] + 1
        _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree, preodr, postodr)
        _bal_avgdist_taxon(adk := np.zeros((n, 2)), m, dm, tree, preodr, postodr)
        res = _bal_min_branch(obs := np.zeros(n), adm, adk, tree, preodr)
        self.assertEqual(res, 4)
        exp = np.array([0, 0, -2, -2, -3, 0, 0], dtype=float)
        npt.assert_array_almost_equal(obs, exp)

        # Test if result matches that calculated from the entire tree.
        dm, tree, preodr, postodr = self.dm3, self.tree3, self.preodr3, self.postodr3
        n = tree.shape[0]
        m = tree[0, 4] + 1
        _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree, preodr, postodr)
        _bal_avgdist_taxon(adk := np.zeros((n, 2)), m, dm, tree, preodr, postodr)
        res = _bal_min_branch(obs := np.zeros(n), adm, adk, tree, preodr)
        exp = np.zeros(n)
        for i in range(n - 2):
            tree_, pre_, post_ = tree.copy(), preodr.copy(), postodr.copy()
            _insert_taxon(m, i, tree_, pre_, post_)
            _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree_, pre_, post_)
            _bal_lengths(lens := np.zeros(n), adm, tree_, pre_)
            exp[i] = lens.sum()
        exp[:n - 2] -= exp[0]

        npt.assert_array_almost_equal(obs, exp)
        self.assertEqual(res, exp[:n - 2].argmin())

    def test_gme(self):
        """The entire tree building workflow."""
        obs = gme(DistanceMatrix(self.dm1, self.taxa1))
        exp = TreeNode.read([self.nwk1])
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        obs = gme(DistanceMatrix(self.dm2, self.taxa2))
        exp = TreeNode.read([self.nwk2])
        self.assertEqual(obs.compare_rfd(exp), 0)

    def test_bme(self):
        """The entire tree building workflow."""
        obs = bme(DistanceMatrix(self.dm1, self.taxa1))
        exp = TreeNode.read([self.nwk1])
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        obs = bme(DistanceMatrix(self.dm2, self.taxa2))
        exp = TreeNode.read([self.nwk2])
        self.assertEqual(obs.compare_rfd(exp), 0)


if __name__ == "__main__":
    main()
