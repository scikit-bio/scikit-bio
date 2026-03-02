# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from heapq import heapify

import numpy as np
import numpy.testing as npt

from skbio import DistanceMatrix, TreeNode
from skbio.util import get_data_path
from skbio.tree._me import (
    gme,
    bme,
    nni,
    _to_treenode,
    _from_treenode,
    _root_from_treenode,
    _insert_taxon_treenode,
    _avgdist_matrix_naive,
    _avgdist_taxon_naive,
    _swap_branches_treenode,
    _swap_branches,
)
from skbio.tree._c_me import (
    _preorder,
    _postorder,
    _calc_tacts,
    _calc_sizes,
    _calc_depths,
    _calc_pairs,
    _insert_taxon,
    _avgdist_matrix,
    _bal_avgdist_matrix,
    _avgdist_taxon,
    _avgdist_taxon_c,
    _bal_avgdist_taxon,
    _bal_avgdist_taxon_c,
    _avgdist_d2_insert,
    _bal_avgdist_insert,
    # _bal_avgdist_insert_p,
    _ols_lengths,
    _ols_lengths_d2,
    _bal_lengths,
    _ols_min_branch_d2,
    _bal_min_branch,
    _avgdist_swap,
    _bal_avgdist_swap,
    _ols_all_swaps,
    _ols_corner_swaps,
    _bal_all_swaps,
    _bal_insert_plan,
    _bal_avgdist_chunk,
    _bal_update_spine,
    _bal_avgdist_nest,
    _bal_avgdist_flat,
)


# def _zero_to_nan(dm):
#     valid = ~np.eye(dm.shape[0], dtype=bool)
#     dm[(dm == 0) & valid] = np.nan


# def _fill_pairs(dm):
#     nans = np.isnan(dm)
#     dm[nans] = dm.T[nans]


def _check_tree(tree, n=None):
    """Check the integrity of an array-based tree structure.

    A tree is represented by a 2D integer array of four columns:

        0. Left child index, or 0 if a tip.
        1. Right child index, or taxon index if a tip.
        2. Parent index, or 0 if root.
        3. Sibling index, or 0 if root.

    """
    if n is None:
        n = tree.shape[0]
    for i in range(n):
        left, right = tree[i, 0], tree[i, 1]
        if left == 0:
            assert i == 0 or right != 0
        else:
            assert tree[left, 2] == tree[right, 2] == i
            assert tree[left, 3] == right
            assert tree[right, 3] == left


def _half_adm(adm, order, n=None):
    """Replace the invalid half of an average distance matrix with zero.

    A cell [i, j] is filled only when i < j in preorder traversal.

    """
    if n is None:
        n = adm.shape[0]
    for i in range(n):
        for j in range(i + 1, n):
            adm[order[j], order[i]] = 0


def _full_adm(adm, order, n=None):
    """Fill the invalid half of an average distance matrix with values.

    A cell [j, i] that suffices i < j in preorder traversal is filled with [i, j].

    """
    if n is None:
        n = adm.shape[0]
    for i in range(n):
        a = order[i]
        for j in range(i + 1, n):
            b = order[j]
            adm[b, a] = adm[a, b]


def _halve_adm(adm, tree, n=None):
    """Replace the invalid half of an average distance matrix with zero.

    This helper function is useful in the tests of `_bal_avgdist_insert_p`, which
    inserts a new taxon into an existing tree under the balanced framework. This
    algorithm has been improved to fill only half of the matrix.

    Specifically, only ancestor-descendant pairs and right-left pairs are filled. In
    other words, only late-early pairs in postorder are filled. This function replaces
    early-late pairs with zero, just to mimic the output of `_bal_avgdist_insert_p`.

    """
    shape = adm.shape if n is None else (n, n)
    for i, j in np.ndindex(*shape):
        if tree[i, 7] < tree[j, 7]:
            adm[i, j] = 0


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
            [1, 2, 0, 0],  # 0: a (root)
            [0, 1, 0, 2],  # 1: b
            [3, 4, 0, 1],  # 2: (c,(e,d))
            [0, 2, 2, 4],  # 3: c
            [5, 6, 2, 3],  # 4: (e,d)
            [0, 4, 4, 6],  # 5: e
            [0, 3, 4, 5],  # 6: d
        ])
        self.lens1 = np.array([2, 3, 3, 4, 2, 1, 2], dtype=float)

        # preorder and postorder and indices for reverse mapping
        self.preodr1 = np.array([0, 1, 2, 3, 4, 5, 6])
        self.preidx1 = np.array([0, 1, 2, 3, 4, 5, 6])
        self.posodr1 = np.array([1, 3, 5, 6, 4, 2, 0])
        self.posidx1 = np.array([6, 0, 5, 1, 4, 2, 3])

        # per-node information
        # number of descending tips (i.e., taxa) (incl. self)
        self.tacts1 = np.array([4, 1, 3, 1, 2, 1, 1])

        # number of descending nodes (incl. self)
        self.sizes1 = np.array([7, 1, 5, 1, 3, 1, 1])

        # number of branches to root (i.e., depth)
        self.depths1 = np.array([0, 1, 1, 2, 2, 3, 3])

        # number of ancestor-descendant pairs
        self.pairs1 = np.array([12, 0, 6, 0, 2, 0, 0])

        # An unfinished tree for example 1 (missing e)
        #           /-b
        # -a-------|
        #          |          /-c
        #           \--------|
        #                     \-d
        self.nwk1m1 = "(b:3.0,(c:4.0,d:4.0):3.0)a:2.0;"
        self.tree1m1 = np.array([
            [1, 2, 0, 0],  # 0: a (root)
            [0, 1, 0, 2],  # 1: b
            [3, 4, 0, 1],  # 2: (c,d)
            [0, 2, 2, 4],  # 3: c
            [0, 3, 2, 3],  # 4: d
            [0, 0, 0, 0],  # 5: empty
            [0, 0, 0, 0],  # 6: empty
        ])
        self.preodr1m1 = np.array([0, 1, 2, 3, 4, 0, 0])
        self.preidx1m1 = np.array([0, 1, 2, 3, 4, 0, 0])
        self.posodr1m1 = np.array([1, 3, 4, 2, 0, 0, 0])
        self.posidx1m1 = np.array([4, 0, 3, 1, 2, 0, 0])
        self.tacts1m1 = np.array([3, 1, 2, 1, 1, 0, 0])
        self.sizes1m1 = np.array([5, 1, 3, 1, 1, 0, 0])
        self.depths1m1 = np.array([0, 1, 1, 2, 2, 0, 0])
        self.pairs1m1 = np.array([6, 0, 2, 0, 0, 0, 0])

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
            [1, 2, 0, 0],  # 0: a (root)
            [3, 4, 0, 2],  # 1: (b,d)
            [5, 6, 0, 1],  # 2: (e,c)
            [0, 1, 1, 4],  # 3: b
            [0, 3, 1, 3],  # 4: d
            [0, 4, 2, 6],  # 5: e
            [0, 2, 2, 5],  # 6: c
        ])
        self.preodr1v2 = np.array([0, 1, 3, 4, 2, 5, 6])
        self.preidx1v2 = np.array([0, 1, 4, 2, 3, 5, 6])
        self.posodr1v2 = np.array([3, 4, 1, 5, 6, 2, 0])
        self.posidx1v2 = np.array([6, 2, 5, 0, 1, 3, 4])
        self.tacts1v2 = np.array([4, 2, 2, 1, 1, 1, 1])
        self.sizes1v2 = np.array([7, 3, 3, 1, 1, 1, 1])
        self.depths1v2 = np.array([0, 1, 1, 2, 2, 2, 2])
        self.pairs1v2 = np.array([10, 2, 2, 0, 0, 0, 0])

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
            [ 1,  2,  0,  0],
            [ 0,  1,  0,  2],
            [ 3,  4,  0,  1],
            [ 0,  2,  2,  4],
            [ 5,  6,  2,  3],
            [ 0,  3,  4,  6],
            [ 7,  8,  4,  5],
            [ 0,  4,  6,  8],
            [ 9, 10,  6,  7],
            [ 0,  6,  8, 10],
            [ 0,  5,  8,  9],
        ])
        self.lens2 = np.array([
            0.91769, 0.76891, 0.42026875, 0.372875, 0.0326305, 0.2705917, 0.03939583,
            0.1525417, 0.0412083, 0.1240375, 0.1451625,
        ])
        self.preodr2 = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        self.preidx2 = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        self.posodr2 = np.array([1, 3, 5, 7, 9, 10, 8, 6, 4, 2, 0])
        self.posidx2 = np.array([10, 0, 9, 1, 8, 2, 7, 3, 6, 4, 5])
        self.tacts2 = np.array([6, 1, 5, 1, 4, 1, 3, 1, 2, 1, 1])
        self.sizes2 = np.array([11, 1, 9, 1, 7, 1, 5, 1, 3, 1, 1])
        self.depths2 = np.array([0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5])
        self.pairs2 = np.array([30, 0, 20, 0, 12, 0, 6, 0, 2, 0, 0])

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
            [1, 2, 0, 0],  # 0, root (a)
            [3, 4, 0, 2],  # 1, ((d,c),e)
            [5, 6, 0, 1],  # 2, (b,f)
            [7, 8, 1, 4],  # 3, (d,c)
            [0, 4, 1, 3],  # 4, e
            [0, 1, 2, 6],  # 5, b
            [0, 5, 2, 5],  # 6, f
            [0, 3, 3, 8],  # 7, d
            [0, 2, 3, 7],  # 8, c
            [0, 0, 0, 0],  # reserved for g's link
            [0, 0, 0, 0],  # reserved for g
        ])
        self.preodr3 = np.array([0, 1, 3, 7, 8, 4, 2, 5, 6, 0, 0])
        self.preidx3 = np.array([0, 1, 6, 2, 5, 7, 8, 3, 4, 0, 0])
        self.posodr3 = np.array([7, 8, 3, 4, 1, 5, 6, 2, 0, 0, 0])
        self.posidx3 = np.array([8, 4, 7, 2, 3, 5, 6, 0, 1, 0, 0])
        self.tacts3 = np.array([5, 3, 2, 2, 1, 1, 1, 1, 1, 0, 0])
        self.sizes3 = np.array([9, 5, 3, 3, 1, 1, 1, 1, 1, 0, 0])
        self.depths3 = np.array([0, 1, 1, 2, 2, 2, 2, 3, 3, 0, 0])
        self.pairs3 = np.array([16, 6, 2, 2, 0, 0, 0, 0, 0, 0, 0])

        # Example 4
        # This example is based on example 3, with the last taxon g inserted as the
        # sibling of e.
        #                               /-d
        #                     /--------|
        #                    |          \-c
        #           /--------|
        #          |         |          /-e
        #          |          \--------|
        # -a-------|                    \-g
        #          |
        #          |          /-b
        #           \--------|
        #                     \-f
        self.dm4 = self.dm3.copy()
        self.taxa4 = self.taxa3.copy()
        self.nwk4 = "(((d,c),(e,g)),(b,f))a;"
        self.tree4 = np.array([
            [ 1,  2,  0,  0],  # 0, root (a)
            [ 3,  9,  0,  2],  # 1, ((d,c),(e,g))
            [ 5,  6,  0,  1],  # 2, (b,f)
            [ 7,  8,  1,  9],  # 3, (d,c)
            [ 0,  4,  9, 10],  # 4, e
            [ 0,  1,  2,  6],  # 5, b
            [ 0,  5,  2,  5],  # 6, f
            [ 0,  3,  3,  8],  # 7, d
            [ 0,  2,  3,  7],  # 8, c
            [ 4, 10,  1,  3],  # 9, (e,g)
            [ 0,  6,  9,  4],  # 10, g
        ])
        self.preodr4 = np.array([0, 1, 3, 7, 8, 9, 4, 10, 2, 5, 6])
        self.preidx4 = np.array([0, 1, 8, 2, 6, 9, 10, 3, 4, 5, 7])
        self.posodr4 = np.array([7, 8, 3, 4, 10, 9, 1, 5, 6, 2, 0])
        self.posidx4 = np.array([10, 6, 9, 2, 3, 7, 8, 0, 1, 5, 4])
        self.tacts4 = np.array([6, 4, 2, 2, 1, 1, 1, 1, 1, 2, 1])
        self.sizes4 = np.array([11, 7, 3, 3, 1, 1, 1, 1, 1, 3, 1])
        self.depths4 = np.array([0, 1, 1, 2, 3, 2, 2, 3, 3, 2, 3])
        self.pairs4 = np.array([22, 10, 2, 2, 0, 0, 0, 0, 0, 2, 0])

        # Example 5
        # This dm can yield negative branch lengths.
        self.dm5 = np.array([
            [0,  5,  9,  9,  800],
            [5,  0, 10, 10,  9],
            [9, 10,  0,  8,  7],
            [9, 10,  8,  0,  3],
            [800,  9,  7,  3,  0],
        ], dtype=float)
        self.taxa5 = list("abcde")
        self.nwk5 = "(b:-129,(c:-95,(d:-130,e:133):101):102,a:134);"

        # Example 6
        # This is a relatively large example, used to demonstrate the parallelizaton
        # planning process in BME.
        self.nwk6 = "(((3,(5,6)4)2,(((10,11)9,(13,14)12)8,(16,17)15)7)1,18)0;"
        #                               /-3
        #                     /2-------|
        #                    |         |          /-5
        #                    |          \4-------|
        #                    |                    \-6
        #                    |
        #           /1-------|                              /-10
        #          |         |                    /9-------|
        #          |         |                   |          \-11
        #          |         |          /8-------|
        #          |         |         |         |          /-13
        #          |         |         |          \12------|
        # -0-------|          \7-------|                    \-14
        #          |                   |
        #          |                   |          /-16
        #          |                    \15------|
        #          |                              \-17
        #          |
        #           \-18

    def test_gme(self):
        """The entire tree building workflow."""
        obs = gme(DistanceMatrix(self.dm1, self.taxa1))
        exp = TreeNode.read([self.nwk1])
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        obs = gme(DistanceMatrix(self.dm2, self.taxa2))
        exp = TreeNode.read([self.nwk2])
        self.assertEqual(obs.compare_rfd(exp), 0)

        # float32
        obs = gme(DistanceMatrix(self.dm2.astype("float32"), self.taxa2))
        exp = TreeNode.read([self.nwk2])
        self.assertEqual(obs.compare_rfd(exp), 0)

        # non-contiguous data
        obs = gme(DistanceMatrix(np.asfortranarray(self.dm2), self.taxa2))
        exp = TreeNode.read([self.nwk2])
        self.assertEqual(obs.compare_rfd(exp), 0)

        # clip to zero
        dm = DistanceMatrix(self.dm5, self.taxa5)
        obs = gme(dm, neg_as_zero=False)
        exp = TreeNode.read([self.nwk5])
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)
        for taxon in ("b", "c", "d"):
            self.assertTrue(obs.find(taxon).length < 0)
        obs = gme(dm, neg_as_zero=True)
        for taxon in ("b", "c", "d"):
            self.assertEqual(obs.find(taxon).length, 0)

    def test_gme_condensed(self):
        """The entire tree building workflow."""
        obs = gme(DistanceMatrix(self.dm1, self.taxa1, condensed=True))
        exp = TreeNode.read([self.nwk1])
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        obs = gme(DistanceMatrix(self.dm2, self.taxa2, condensed=True))
        exp = TreeNode.read([self.nwk2])
        self.assertEqual(obs.compare_rfd(exp), 0)

        # clip to zero
        dm = DistanceMatrix(self.dm5, self.taxa5, condensed=True)
        obs = gme(dm, neg_as_zero=False)
        exp = TreeNode.read([self.nwk5])
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)
        for taxon in ("b", "c", "d"):
            self.assertTrue(obs.find(taxon).length < 0)
        obs = gme(dm, neg_as_zero=True)
        for taxon in ("b", "c", "d"):
            self.assertEqual(obs.find(taxon).length, 0)

    def test_gme_real(self):
        """Test on a real dataset with 100 taxa."""
        dm = DistanceMatrix.read(get_data_path("mp100.phy"))
        obs = gme(dm, neg_as_zero=False)
        exp = TreeNode.read(get_data_path("mp100.gme.nwk"))
        self.assertEqual(obs.compare_rfd(exp), 0.0)
        self.assertAlmostEqual(obs.compare_cophenet(exp), 0.0)

    def test_bme(self):
        """The entire tree building workflow."""
        # simple integer example
        obs = bme(DistanceMatrix(self.dm1, self.taxa1))
        exp = TreeNode.read([self.nwk1])
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        # simple float example
        obs = bme(DistanceMatrix(self.dm2, self.taxa2))
        exp = TreeNode.read([self.nwk2])
        self.assertEqual(obs.compare_rfd(exp), 0)

        # float32
        obs = bme(DistanceMatrix(self.dm2.astype("float32"), self.taxa2))
        exp = TreeNode.read([self.nwk2])
        self.assertEqual(obs.compare_rfd(exp), 0)

        # non-contiguous data
        obs = bme(DistanceMatrix(np.asfortranarray(self.dm2), self.taxa2))
        exp = TreeNode.read([self.nwk2])
        self.assertEqual(obs.compare_rfd(exp), 0)

        # parallelization
        obs = bme(DistanceMatrix(self.dm2, self.taxa2), parallel=True)
        exp = TreeNode.read([self.nwk2])
        self.assertEqual(obs.compare_rfd(exp), 0)

        # clip to zero
        dm = DistanceMatrix(self.dm5, self.taxa5)
        obs = bme(dm, neg_as_zero=False)
        for taxon in ("b", "c", "d"):
            self.assertTrue(obs.find(taxon).length < 0)
        obs = bme(dm, neg_as_zero=True)
        for taxon in ("b", "c", "d"):
            self.assertEqual(obs.find(taxon).length, 0)

    def test_bme_condensed(self):
        """The entire tree building workflow."""
        obs = bme(DistanceMatrix(self.dm1, self.taxa1, condensed=True))
        exp = TreeNode.read([self.nwk1])
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        obs = bme(DistanceMatrix(self.dm2, self.taxa2, condensed=True))
        exp = TreeNode.read([self.nwk2])
        self.assertEqual(obs.compare_rfd(exp), 0)

        # clip to zero
        dm = DistanceMatrix(self.dm5, self.taxa5, condensed=True)
        obs = bme(dm, neg_as_zero=False)
        for taxon in ("b", "c", "d"):
            self.assertTrue(obs.find(taxon).length < 0)
        obs = bme(dm, neg_as_zero=True)
        for taxon in ("b", "c", "d"):
            self.assertEqual(obs.find(taxon).length, 0)

    def test_bme_real(self):
        """Test on a real dataset with 100 taxa."""
        dm = DistanceMatrix.read(get_data_path("mp100.phy"))
        obs = bme(dm, neg_as_zero=False)
        exp = TreeNode.read(get_data_path("mp100.bme.nwk"))
        self.assertEqual(obs.compare_rfd(exp), 0.0)
        self.assertAlmostEqual(obs.compare_cophenet(exp), 0.0)

        # disable parallelization (no effect, since tree is small)
        obs = bme(dm, neg_as_zero=False, parallel=False)
        self.assertAlmostEqual(obs.compare_cophenet(exp), 0.0)

        # always enable parallelization
        obs = bme(dm, neg_as_zero=False, parallel=True)
        self.assertAlmostEqual(obs.compare_cophenet(exp), 0.0)

        # specify a threshold for parallelization
        # (50 will not enter phase 3)
        obs = bme(dm, neg_as_zero=False, parallel=50)
        self.assertAlmostEqual(obs.compare_cophenet(exp), 0.0)
        # (15 will enter phase 3)
        obs = bme(dm, neg_as_zero=False, parallel=15)
        self.assertAlmostEqual(obs.compare_cophenet(exp), 0.0)

    def test_nni(self):
        """The entire tree rearrangement workflow."""
        # Test if NNI can convert an incorrect tree into the groud truth tree.

        # Example 1: In this simple example, FastNNI and BNNI should produce the same
        # topology and branch lengths.
        tree = TreeNode.read([self.nwk1v2])
        tree.append(TreeNode(tree.name, tree.length))
        tree.name, tree.length = None, None

        dm = DistanceMatrix(self.dm1, self.taxa1)

        exp = TreeNode.read([self.nwk1])
        exp.append(TreeNode(exp.name, exp.length))
        exp.name, exp.length = None, None

        obs = nni(tree, dm, balanced=False)
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        obs = nni(tree, dm, balanced=True)
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        dm = DistanceMatrix(self.dm1.astype("float32"), self.taxa1)
        for balanced in True, False:
            obs = nni(tree, dm, balanced=balanced)
            self.assertEqual(obs.compare_rfd(exp), 0)
            self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        # Example 4: The same, except for branch lengths.
        tree = TreeNode.read([self.nwk4])
        tree.append(TreeNode(tree.name, tree.length))
        tree.name, tree.length = None, None

        dm = DistanceMatrix(self.dm4, self.taxa4)

        exp = TreeNode.read([self.nwk2])
        exp.append(TreeNode(exp.name, exp.length))
        exp.name, exp.length = None, None
        taxonmap = dict(zip(self.taxa2, self.taxa4))
        for tip in exp.tips():
            tip.name = taxonmap[tip.name]

        obs = nni(tree, dm, balanced=False)
        self.assertEqual(obs.compare_rfd(exp), 0)

        obs = nni(tree, dm, balanced=True)
        self.assertEqual(obs.compare_rfd(exp), 0)

        # clip to zero
        dm = DistanceMatrix(self.dm5, self.taxa5)
        tree = TreeNode.read(["(((a,b),d),(c,e));"])
        obs = nni(tree, dm, balanced=False, neg_as_zero=False)
        exp = TreeNode.read([self.nwk5])
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)
        for taxon in ("b", "c", "d"):
            self.assertTrue(obs.find(taxon).length < 0)
        obs = nni(tree, dm, balanced=False, neg_as_zero=True)
        for taxon in ("b", "c", "d"):
            self.assertEqual(obs.find(taxon).length, 0)

    def test_nni_condensed(self):
        """The entire tree rearrangement workflow."""
        # Test if NNI can convert an incorrect tree into the groud truth tree.

        # Example 1: In this simple example, FastNNI and BNNI should produce the same
        # topology and branch lengths.
        tree = TreeNode.read([self.nwk1v2])
        tree.append(TreeNode(tree.name, tree.length))
        tree.name, tree.length = None, None

        dm = DistanceMatrix(self.dm1, self.taxa1, condensed=True)

        exp = TreeNode.read([self.nwk1])
        exp.append(TreeNode(exp.name, exp.length))
        exp.name, exp.length = None, None

        obs = nni(tree, dm, balanced=False)
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        obs = nni(tree, dm, balanced=True)
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)

        # Example 4: The same, except for branch lengths.
        tree = TreeNode.read([self.nwk4])
        tree.append(TreeNode(tree.name, tree.length))
        tree.name, tree.length = None, None

        dm = DistanceMatrix(self.dm4, self.taxa4, condensed=True)

        exp = TreeNode.read([self.nwk2])
        exp.append(TreeNode(exp.name, exp.length))
        exp.name, exp.length = None, None
        taxonmap = dict(zip(self.taxa2, self.taxa4))
        for tip in exp.tips():
            tip.name = taxonmap[tip.name]

        obs = nni(tree, dm, balanced=False)
        self.assertEqual(obs.compare_rfd(exp), 0)

        obs = nni(tree, dm, balanced=True)
        self.assertEqual(obs.compare_rfd(exp), 0)

        # clip to zero
        dm = DistanceMatrix(self.dm5, self.taxa5, condensed=True)
        tree = TreeNode.read(["(((a,b),d),(c,e));"])
        obs = nni(tree, dm, balanced=False, neg_as_zero=False)
        exp = TreeNode.read([self.nwk5])
        self.assertEqual(obs.compare_rfd(exp), 0)
        self.assertAlmostEqual(obs.compare_cophenet(exp, ignore_self=True), 0)
        for taxon in ("b", "c", "d"):
            self.assertTrue(obs.find(taxon).length < 0)
        obs = nni(tree, dm, balanced=False, neg_as_zero=True)
        for taxon in ("b", "c", "d"):
            self.assertEqual(obs.find(taxon).length, 0)

    def test_nni_real(self):
        """Test on a real dataset with 100 taxa."""
        # FastNNI
        tree = TreeNode.read(get_data_path("mp100.gme.nwk"))
        dm = DistanceMatrix.read(get_data_path("mp100.phy"))
        obs = nni(tree, dm, balanced=False, neg_as_zero=False)
        exp = TreeNode.read(get_data_path("mp100.gme.nni.nwk"))
        self.assertNotEqual(obs.compare_rfd(tree), 0.0)
        self.assertEqual(obs.compare_rfd(exp), 0.0)
        self.assertAlmostEqual(obs.compare_cophenet(exp), 0.0)

        # BNNI
        tree = TreeNode.read(get_data_path("mp100.bme.nwk"))
        dm = DistanceMatrix.read(get_data_path("mp100.phy"))
        obs = nni(tree, dm, balanced=True, neg_as_zero=False)
        exp = TreeNode.read(get_data_path("mp100.bme.nni.nwk"))
        self.assertNotEqual(obs.compare_rfd(tree), 0.0)
        self.assertEqual(obs.compare_rfd(exp), 0.0)
        self.assertAlmostEqual(obs.compare_cophenet(exp), 0.0)

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
        # a complete tree
        obj = TreeNode.read([self.nwk1])
        taxmap = {x: i for i, x in enumerate(self.taxa1)}
        n = len(taxmap) * 2 - 3
        obs = np.empty((n, 4), dtype=int)
        num = _from_treenode(obj, taxmap, obs)
        self.assertEqual(num, n)

        # check if tree topology is correct
        _check_tree(obs)

        # check if nodes are in preorder
        stack = np.empty(n, dtype=int)
        _preorder(order := np.empty(n, dtype=int), obs, stack)
        npt.assert_array_equal(order, np.arange(n))

        # check if tree topology matches expectation
        npt.assert_array_equal(obs, self.tree1)

        # an incomplete tree
        obj = TreeNode.read([self.nwk1m1])
        num = _from_treenode(obj, taxmap, obs)
        self.assertEqual(num, n - 2)
        _check_tree(obs, n - 2)
        obs[-2:] = 0  # fill unused cells before comparison
        npt.assert_array_equal(obs, self.tree1m1)

        # start from index 2
        _from_treenode(obj, taxmap, obs, pos=2)
        exp = np.array([[3, 4, 0, 0],
                        [0, 1, 2, 4],
                        [5, 6, 2, 3],
                        [0, 2, 4, 6],
                        [0, 3, 4, 5]])
        npt.assert_array_equal(obs[2:], exp)

        # non-binary tree
        obj = TreeNode.read(["((b,c,d),e)a;"])
        with self.assertRaises(ValueError) as cm:
            _from_treenode(obj, taxmap, obs)
        self.assertEqual(str(cm.exception), "Tree is not strictly bifurcating.")

        # another example
        obj = TreeNode.read([self.nwk3])
        taxmap = {x: i for i, x in enumerate(self.taxa3)}
        n = len(taxmap) * 2 - 3
        obs = np.empty((n, 4), dtype=int)
        _from_treenode(obj, taxmap, obs)
        _check_tree(obs, n - 2)
        # note: self.tree3 is not in preorder so we can't directly compare

    def test_root_from_treenode(self):
        """Convert TreeNode into tree array rooted at 1st taxon."""
        # Example 1: first move the root taxon ("a") to the child place, making the
        # tree trifurcating. This is the typical output of tree-building methods.
        obj = TreeNode.read([self.nwk1])
        obj.append(TreeNode(obj.name))
        obj.name = None

        # Now convert to array to see if it matches the reference.
        obs = _root_from_treenode(obj, self.taxa1)
        _check_tree(obs)
        npt.assert_array_equal(obs, self.tree1)

        # Test a different order of taxa. This time, taxon "c" will become the root.
        taxa = list("cadeb")
        obs = _root_from_treenode(obj, taxa)
        exp = np.array([
            [1, 4, 0, 0],  # c
            [2, 3, 0, 4],  # (e,d)
            [0, 3, 1, 3],  # e
            [0, 2, 1, 2],  # d
            [5, 6, 0, 1],  # (b,a)
            [0, 4, 4, 6],  # b
            [0, 1, 4, 5],  # a
        ])
        _check_tree(obs)
        npt.assert_array_equal(obs, exp)

        # Root the tree at an internal branch and test to recover the same result.
        obj = obj.find("e").parent
        obj = obj.root_at(above=True)
        self.assertEqual(len(obj.children), 2)
        obs = _root_from_treenode(obj, self.taxa1)
        _check_tree(obs)
        npt.assert_array_equal(obs, self.tree1)

        # Example 4: Test all possible roots (taxa).
        obj = TreeNode.read([self.nwk4])
        obj.append(TreeNode(obj.name))
        obj.name = None
        taxa = self.taxa4
        m = len(taxa)
        n = m * 2 - 3

        for i in range(m):
            # move target taxon to the first place
            taxon_ = taxa[i]
            taxa_ = [taxon_] + taxa[:i] + taxa[i + 1:]

            # use the algorithm to generate the array
            obs = _root_from_treenode(obj, taxa_)
            _check_tree(obs)

            # root the TreeNode object instead (slower) and generate the array
            obj_ = obj.root_at(taxon_)
            obj_.prune()

            taxmap = {x: i for i, x in enumerate(taxa_)}
            _from_treenode(obj_, taxmap, exp := np.empty((n, 4), dtype=int))
            npt.assert_array_equal(obs, exp)

        # non-binary trees
        taxa = list("abcde")
        nwks = [
            "((a,b,c),(d,e));",  # internal node has 3 children
            "((a),b,(c,(d,e)));",  # internal node has 1 child
            # "(((b,c),((d,e),a)));",  # root has 1 child  # TODO: This doesn't work.
            "((a,b),c,d,e);",  # root has 4 children
        ]
        for nwk in nwks:
            obj = TreeNode.read([nwk])
            self.assertRaises(ValueError, _root_from_treenode, obj, taxa)

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
        npt.assert_array_equal(obs, self.posodr1)

        # given clade
        _postorder(obs := np.full(n, 0), self.tree1, stack, start=2)
        exp = np.array([3, 5, 6, 4, 2, 0, 0])
        npt.assert_array_equal(obs, exp)

        # alternative tree
        _postorder(obs := np.full(n, 0), self.tree1v2, stack)
        npt.assert_array_equal(obs, self.posodr1v2)

        # another example
        stack = np.full(n := self.tree3.shape[0], 0)
        _postorder(obs := np.full(n, 0), self.tree3, stack)
        npt.assert_array_equal(obs, self.posodr3)

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
        # insert taxon e (4) above d (order=4, size=1)
        obs = self.tree1m1.copy(), self.preodr1m1.copy()
        _insert_taxon(4, 4, 1, *obs)
        exp = np.array([
            [1, 2, 0, 0],
            [0, 1, 0, 2],
            [3, 5, 0, 1],
            [0, 2, 2, 5],
            [0, 3, 5, 6],
            [4, 6, 2, 3],
            [0, 4, 5, 4],
        ]), np.array([
            0, 1, 2, 3, 5, 4, 6
        ])
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # insert taxon e (4) above b (order=1, size=1)
        obs = self.tree1m1.copy(), self.preodr1m1.copy()
        _insert_taxon(4, 1, 1, *obs)
        exp = np.array([
            [5, 2, 0, 0],
            [0, 1, 5, 6],
            [3, 4, 0, 5],
            [0, 2, 2, 4],
            [0, 3, 2, 3],
            [1, 6, 0, 2],
            [0, 4, 5, 1],
        ]), np.array([
            0, 5, 1, 6, 2, 3, 4
        ])
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # insert taxon e (4) above (c,d) (order=2, size=3)
        obs = self.tree1m1.copy(), self.preodr1m1.copy()
        _insert_taxon(4, 2, 3, *obs)
        exp = np.array([
            [1, 5, 0, 0],
            [0, 1, 0, 5],
            [3, 4, 5, 6],
            [0, 2, 2, 4],
            [0, 3, 2, 3],
            [2, 6, 0, 1],
            [0, 4, 5, 2],
        ]), np.array([
            0, 1, 5, 2, 3, 4, 6
        ])
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # insert taxon e (4) into the root branch (idx=0, size=5)
        obs = self.tree1m1.copy(), self.preodr1m1.copy()
        _insert_taxon(4, 0, 5, *obs)
        exp = np.array([
            [5, 6, 0, 0],
            [0, 1, 5, 2],
            [3, 4, 5, 1],
            [0, 2, 2, 4],
            [0, 3, 2, 3],
            [1, 2, 0, 6],
            [0, 4, 0, 5],
        ]), np.array([
            0, 5, 1, 2, 3, 4, 6
        ])
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # another example: all possible insertions
        taxa, tree, order, sizes = self.taxa3, self.tree3, self.preodr3, self.sizes3
        n = tree.shape[0] - 2
        k = len(taxa) - 1
        taxamap = {}

        for node in self.preodr3[:-2][::-1]:
            if tree[node, 0] == 0:
                taxamap[node] = [taxa[tree[node, 1]]]
            else:
                taxamap[node] = taxamap[tree[node, 0]] + taxamap[tree[node, 1]]
        obj = TreeNode.read([self.nwk3])

        for i in range(n):
            tag = order[i]
            _insert_taxon(k, i, sizes[tag], tree_ := tree.copy(), order.copy())
            obs = _to_treenode(tree_, taxa)
            exp = obj.copy()
            _insert_taxon_treenode("g", exp.lca(taxamap[tag]), exp)
            self.assertEqual(obs.compare_rfd(exp), 0)

    def test_avgdist_matrix_naive(self):
        """Calculate an average distance matrix using a naive method."""
        # Test if the algorithm reproduces the manually calculated result.
        n = self.tree1.shape[0]
        obs = np.zeros((n, n), dtype=float)
        _avgdist_matrix_naive(obs, self.dm1, self.tree1, self.preodr1, self.tacts1)
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
        _avgdist_matrix_naive(
            obs, self.dm1, self.tree1m1, self.preodr1m1, self.tacts1m1
        )
        exp = np.array([
            [ 0. ,  5. ,  9. ,  9. ,  9. ,  0. ,  0. ],
            [ 5. ,  0. , 10. , 10. , 10. ,  0. ,  0. ],
            [ 9. , 10. ,  0. ,  9.5,  9.5,  0. ,  0. ],
            [ 9. , 10. ,  9.5,  0. ,  8. ,  0. ,  0. ],
            [ 9. , 10. ,  9.5,  8. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
        ])
        npt.assert_allclose(obs, exp)

        # alternative tree
        obs = np.zeros((n, n), dtype=float)
        _avgdist_matrix_naive(
            obs, self.dm1, self.tree1v2, self.preodr1v2, self.tacts1v2
        )
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

    def test_avgdist_matrix(self):
        """Calculate an average distance matrix."""
        # Test if the algorithm produces the same result as the native method does.
        dm, tree, order, tacts = self.dm1, self.tree1, self.preodr1, self.tacts1
        n = tree.shape[0]
        obs = np.zeros((n, n), dtype=float)
        _avgdist_matrix(obs, dm, tree, order, tacts)
        exp = np.zeros((n, n), dtype=float)
        _avgdist_matrix_naive(exp, dm, tree, order, tacts)
        npt.assert_allclose(obs, exp)

        # make sure all cells are populated
        _avgdist_matrix(obs := np.full((n, n), np.nan), dm, tree, order, tacts)
        np.fill_diagonal(obs, 0)
        _avgdist_matrix_naive(exp := np.zeros((n, n)), dm, tree, order, tacts)
        npt.assert_allclose(obs, exp)

        # incomplete tree
        tree, order, tacts = self.tree1m1, self.preodr1m1, self.tacts1m1
        _avgdist_matrix(obs := np.zeros((n, n)), dm, tree, order, tacts)
        _avgdist_matrix_naive(exp := np.zeros((n, n)), dm, tree, order, tacts)
        npt.assert_allclose(obs, exp)

        # example 1 v2
        tree, order, tacts = self.tree1v2, self.preodr1v2, self.tacts1v2
        _avgdist_matrix(obs := np.zeros((n, n)), dm, tree, order, tacts)
        _avgdist_matrix_naive(exp := np.zeros((n, n)), dm, tree, order, tacts)
        npt.assert_allclose(obs, exp)

        # example 2 (complete)
        dm, tree, order, tacts = self.dm2, self.tree2, self.preodr2, self.tacts2
        n = tree.shape[0]
        _avgdist_matrix(obs := np.zeros((n, n)), dm, tree, order, tacts)
        _avgdist_matrix_naive(exp := np.zeros((n, n)), dm, tree, order, tacts)
        npt.assert_allclose(obs, exp)

        # example 3 (incomplete)
        dm, tree, order, tacts = self.dm3, self.tree3, self.preodr3, self.tacts3
        n = tree.shape[0]
        _avgdist_matrix(obs := np.zeros((n, n)), dm, tree, order, tacts)
        _avgdist_matrix_naive(exp := np.zeros((n, n)), dm, tree, order, tacts)
        npt.assert_allclose(obs, exp)

    def test_bal_avgdist_matrix(self):
        """Calculate a balanced average distance matrix."""
        dm, tree, order, sizes = self.dm1, self.tree1, self.preodr1, self.sizes1
        n = sizes[0]
        _bal_avgdist_matrix(obs := np.zeros((n, n)), dm, tree, order, sizes)
        exp = np.array([
            [ 0.  ,  5.  ,  8.75,  9.  ,  8.5 ,  8.  ,  9.  ],
            [ 5.  ,  0.  ,  9.75, 10.  ,  9.5 ,  9.  , 10.  ],
            [ 8.75,  9.75,  0.  ,  9.5 ,  9.  ,  8.5 ,  9.5 ],
            [ 9.  , 10.  ,  9.5 ,  0.  ,  7.5 ,  7.  ,  8.  ],
            [ 8.5 ,  9.5 ,  9.  ,  7.5 ,  0.  ,  7.75,  8.75],
            [ 8.  ,  9.  ,  8.5 ,  7.  ,  7.75,  0.  ,  3.  ],
            [ 9.  , 10.  ,  9.5 ,  8.  ,  8.75,  3.  ,  0.  ],
        ])
        npt.assert_allclose(obs, exp)

        _bal_avgdist_matrix(obs := np.full((n, n), np.nan), dm, tree, order, sizes)
        np.fill_diagonal(obs, 0)
        npt.assert_allclose(obs, exp)

        dm, tree, order, sizes = self.dm3, self.tree3, self.preodr3, self.sizes3
        n = sizes[0]
        _bal_avgdist_matrix(obs := np.zeros((n, n)), dm, tree, order, sizes)
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
        obs = np.zeros((2, n), dtype=float)
        _avgdist_taxon_naive(
            obs, taxon, self.dm1, self.tree1m1, self.preodr1m1, self.tacts1m1
        )
        exp = np. array([[6.333, 9.   , 5.   , 7.   , 3.   , 0.   , 0.   ],
                         [8.   , 6.   , 8.5  , 6.667, 8.   , 0.   , 0.   ]])
        npt.assert_array_equal(obs.round(3), exp)

    def test_avgdist_taxon(self):
        """Calculate taxon-to-subtree average distances."""
        # Test if the algorithm produces the same result as the naive method does.
        dm = self.dm1
        taxon = dm.shape[0] - 1
        tree, order, tacts = self.tree1m1, self.preodr1m1, self.tacts1m1
        n = tree.shape[0]
        _avgdist_taxon(obs := np.zeros((2, n)), taxon, dm, tree, order, tacts)
        _avgdist_taxon_naive(exp := np.zeros((2, n)), taxon, dm, tree, order, tacts)
        npt.assert_allclose(obs, exp)

        # non-contiguous data
        dmf = np.asfortranarray(dm)
        _avgdist_taxon(obs := np.zeros((2, n)), taxon, dmf, tree, order, tacts)
        npt.assert_allclose(obs, exp)

        # contiguous function
        _avgdist_taxon_c(obs := np.zeros((2, n)), taxon, dm, tree, order, tacts)
        npt.assert_allclose(obs, exp)

    def test_bal_avgdist_taxon(self):
        """Calculate taxon-to-subtree balanced average distances."""
        dm = self.dm1
        tree, order = self.tree1m1, self.preodr1m1
        k = dm.shape[0] - 1
        n = tree.shape[0]
        obs = np.zeros((2, n))
        _bal_avgdist_taxon(n - 2, k, dm, obs, tree, order)
        exp = np.array([[7.  , 9.  , 5.  , 7.  , 3.  , 0.  , 0.  ],
                        [8.  , 6.5 , 8.5 , 5.75, 7.75, 0.  , 0.  ]])
        npt.assert_allclose(obs, exp)

        # non-contiguous data
        _bal_avgdist_taxon(n - 2, k, np.asfortranarray(dm), obs, tree, order)
        npt.assert_allclose(obs, exp)

        # contiguous function
        _bal_avgdist_taxon_c(n - 2, k, dm, obs, tree, order)
        npt.assert_allclose(obs, exp)

        dm = self.dm3
        tree, order = self.tree3, self.preodr3
        k = dm.shape[0] - 1
        n = tree.shape[0]
        obs = np.zeros((2, n))
        _bal_avgdist_taxon(n - 2, k, dm, obs, tree, order)
        exp = np.array([[0.639, 0.411, 0.866, 0.515, 0.308, 1.463, 0.269, 0.471, 0.558,
                         0.   , 0.   ],
                        [1.59 , 1.228, 1.001, 0.768, 0.871, 0.635, 1.232, 0.663, 0.62 ,
                         0.   , 0.   ]])
        npt.assert_array_equal(obs.round(3), exp)

    def test_avgdist_d2_insert(self):
        """Update distant-2 subtree average distances after taxon insertion."""
        n = self.tree1.shape[0]
        taxon = self.dm1.shape[0] - 1
        # The following values were taken from the full-scale distance matrix.
        ad2 = np.array([
            [ 0. ,  0. ],  # 0: a (root)
            [10. ,  5. ],  # 1: b
            [10. ,  9. ],  # 2: (c,d)
            [ 8. ,  9.5],  # 3: c
            [ 8. ,  9.5],  # 4: d
            [ 0. ,  0. ],  # 5: empty
            [ 0. ,  0. ],  # 6: empty
        ]).T
        # Insert e as a sibling of d. This should recover tree1.
        target = 4
        _avgdist_taxon(
            adk := np.zeros((2, n), dtype=float), taxon, self.dm1, self.tree1m1,
            self.preodr1m1, self.tacts1m1
        )
        _avgdist_d2_insert(
            obs := ad2.copy(), target, adk, self.tree1m1, self.preodr1m1, self.tacts1m1
        )
        # The following values were taken from the full-scale distance matrix.
        exp = np.array([
            [0.   , 0.   ],  # 0: a (root)
            [9.667, 5.   ],  # 1: b
            [9.667, 8.667],  # 2: (c,(e,d))
            [7.5  , 9.5  ],  # 3: c
            [3.   , 9.   ],  # 4: d
            [7.5  , 9.  ],  # 5: (e,d)
            [3.   , 8.   ],  # 6: e
        ]).T
        npt.assert_array_equal(obs.round(3), exp)

        # another example: all possible insertions
        dm, tree, order, tacts = self.dm3, self.tree3, self.preodr3, self.tacts3
        n = tree.shape[0]
        k = tacts[0] + 1
        ran_ = np.arange(n)

        # get distant-2 values from the full matrix
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, tacts)
        ad2 = np.ascontiguousarray(
            np.vstack([adm[ran_, tree[ran_, 3]], adm[ran_, tree[ran_, 2]]])
        )
        _avgdist_taxon(adk := np.zeros((2, n)), k, dm, tree, order, tacts)

        for i in range(n - 2):
            tree_, order_, tacts_ = tree.copy(), order.copy(), tacts.copy()

            # calculate distant-2 values using the algorithm
            _avgdist_d2_insert(obs := ad2.copy(), i, adk, tree_, order_, tacts_)

            # insert taxon and calculate full matrix
            _insert_taxon(k, i, tacts_[order_[i]] * 2 - 1, tree_, order_)
            _avgdist_matrix(adm := np.zeros((n, n)), dm, tree_, order_, tacts_)
            exp = np.ascontiguousarray(
                np.vstack([adm[ran_, tree_[ran_, 3]], adm[ran_, tree_[ran_, 2]]])
            )

            npt.assert_allclose(obs, exp)

    def test_bal_avgdist_insert_1(self):
        """Update balanced average distance matrix after taxon insertion.

        Phase 1 (serial).

        """
        # Test if the algorithm produces the same result as calculated from the full
        # matrix after taxon insertion.
        tree, order, sizes, depths = (
            self.tree1m1, self.preodr1m1, self.sizes1m1, self.depths1m1
        )
        n = tree.shape[0]
        npots = 2.0 ** -np.arange((n + 1) // 2)
        diffs = np.zeros(n)
        ancs = np.zeros(n, dtype=int)
        ancx = np.zeros(n, dtype=int)
        adm = np.array([
            [ 0. ,  5. ,  9. ,  9. ,  9. ,  0. ,  0. ],
            [ 0. ,  0. , 10. , 10. , 10. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  9.5,  9.5,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  8. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
        ])
        adk = np.array([
            [7.  , 9.  , 5.  , 7.  , 3.  , 0.  , 0.  ],
            [8.  , 6.5 , 8.5 , 5.75, 7.75, 0.  , 0.  ],
        ])
        # Insert e as a sibling of d. This should recover tree1.
        itag = 4
        _bal_avgdist_insert(
            n - 2, itag, obs := adm, adk, tree, order, sizes, depths, diffs, npots,
            ancs, ancx
        )
        exp = np.array([
            [ 0.  ,  5.  ,  8.75,  9.  ,  9.  ,  8.5 ,  8.  ],
            [ 0.  ,  0.  ,  9.75, 10.  , 10.  ,  9.5 ,  9.  ],
            [ 0.  ,  0   ,  0.  ,  9.5 ,  9.5 ,  9.  ,  8.5 ],
            [ 0.  ,  0.  ,  0.  ,  0.  ,  8.  ,  7.5 ,  7.  ],
            [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  3.  ],
            [ 0.  ,  0.  ,  0.  ,  0.  ,  8.75,  0.  ,  7.75],
            [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
        ])
        npt.assert_allclose(obs, exp)

        # another example: all possible insertions
        dm, tree, order, sizes, depths = (
            self.dm3, self.tree3, self.preodr3, self.sizes3, self.depths3
        )
        n = tree.shape[0]
        m = dm.shape[0] - 1
        _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, sizes)
        _half_adm(adm, order, n - 2)
        _bal_avgdist_taxon(n - 2, m, dm, adk := np.zeros((2, n)), tree, order)
        npots = 2.0 ** -np.arange(m)
        diffs = np.zeros(n)
        ancs = np.zeros(n, dtype=int)
        ancx = np.zeros(n, dtype=int)

        for i in range(n - 2):
            # update matrix using the algorithm
            tree_, order_, sizes_ = tree.copy(), order.copy(), sizes.copy()
            _bal_avgdist_insert(
                n - 2, i, obs := adm.copy(), adk.copy(), tree_, order_, sizes_,
                depths.copy(), diffs, npots, ancs, ancx
            )
            # insert taxon and calculate full matrix
            _insert_taxon(m, i, sizes_[order_[i]], tree_, order_)
            _bal_avgdist_matrix(exp := np.zeros((n, n)), dm, tree_, order_, sizes_)
            _half_adm(exp, order_, n)

            npt.assert_allclose(obs, exp)

    def test_bal_avgdist_insert_2(self):
        """Update balanced average distance matrix after taxon insertion.

        Phase 2 (flat serial, nested parallel).

        """
        # Test if the algorithm produces the same result as calculated from the full
        # matrix after taxon insertion.
        tree, order, sizes, depths, pairs = (
            self.tree1m1, self.preodr1m1, self.sizes1m1, self.depths1m1,
            self.pairs1m1
        )
        n = tree.shape[0]
        npots = 2.0 ** -np.arange((n + 1) // 2)
        diffs = np.zeros(n)
        ancs = np.zeros(n, dtype=int)
        ancx = np.zeros(n, dtype=int)
        segs = np.empty(n + 1, dtype=int)
        lvls = np.empty(n, dtype=int)
        oops = np.empty(n, dtype=int)
        chunks = np.empty(n + 1, dtype=int)
        chunks[0] = 0
        chusegs = np.empty(n, dtype=int)
        chusegs[0] = 0
        adm = np.array([
            [ 0. ,  5. ,  9. ,  9. ,  9. ,  0. ,  0. ],
            [ 0. ,  0. , 10. , 10. , 10. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  9.5,  9.5,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  8. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
        ])
        adk = np.array([
            [7.  , 9.  , 5.  , 7.  , 3.  , 0.  , 0.  ],
            [8.  , 6.5 , 8.5 , 5.75, 7.75, 0.  , 0.  ],
        ])
        # Insert e as a sibling of d. This should recover tree1.
        itag = 4
        tag = order[itag]
        depth = depths[tag]
        _bal_insert_plan(
            n - 2, itag, adm, adk, tree, order, sizes, depths, pairs, diffs, ancs,
            ancx, segs, lvls, oops, True)

        # The target node has a depth of 2. Therefore it has two ancestors (parent and
        # root).
        obs = (segs == n - 2).argmax()
        self.assertEqual(obs, 3)
        npt.assert_array_equal(segs[:obs + 1], [0, 2, 4, 5])
        npt.assert_array_equal(lvls[:obs], [1, 0, -1])
        npt.assert_array_equal(oops[:obs], [2, 0, 0])

        # This is an extremely simple case. Only the root node has 2 operations, while
        # all other nodes have 0. Therefore, regardless how many chunks are desired,
        # there will always be one chunk generated, which includes just the root.
        enc = 3  # aim for 3 chunks (get 1)
        nc = _bal_avgdist_chunk(
            n - 2, itag, order, sizes, pairs, segs, lvls, oops, enc, chunks, chusegs)
        self.assertEqual(nc, 1)
        _bal_update_spine(tag, depth, sizes, pairs, ancs)
        _bal_avgdist_nest(
            itag, depth, obs := adm, order, sizes, depths, diffs, npots, ancs, ancx,
            segs, lvls, nc, chunks, chusegs)
        exp = np.array([
            [ 0.  ,  5.  ,  8.75,  9.  ,  9.  ,  8.5 ,  8.  ],
            [ 0.  ,  0.  ,  9.75, 10.  , 10.  ,  9.5 ,  9.  ],
            [ 0.  ,  0   ,  0.  ,  9.5 ,  9.5 ,  9.  ,  8.5 ],
            [ 0.  ,  0.  ,  0.  ,  0.  ,  8.  ,  7.5 ,  7.  ],
            [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  3.  ],
            [ 0.  ,  0.  ,  0.  ,  0.  ,  8.75,  0.  ,  7.75],
            [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
        ])
        npt.assert_allclose(obs, exp)

    def test_bal_avgdist_insert_3(self):
        """Update balanced average distance matrix after taxon insertion.

        Phase 3 (flat and nested both parallel).

        """
        # Test if the algorithm produces the same result as calculated from the full
        # matrix after taxon insertion.
        tree, order, sizes, depths, pairs = (
            self.tree1m1, self.preodr1m1, self.sizes1m1, self.depths1m1,
            self.pairs1m1
        )
        n = tree.shape[0]
        npots = 2.0 ** -np.arange((n + 1) // 2)
        diffs = np.zeros(n)
        ancs = np.zeros(n, dtype=int)
        ancx = np.zeros(n, dtype=int)
        segs = np.empty(n + 1, dtype=int)
        lvls = np.empty(n, dtype=int)
        oops = np.empty(n, dtype=int)
        chunks = np.empty(n + 1, dtype=int)
        chunks[0] = 0
        chusegs = np.empty(n, dtype=int)
        chusegs[0] = 0
        adm = np.array([
            [ 0. ,  5. ,  9. ,  9. ,  9. ,  0. ,  0. ],
            [ 0. ,  0. , 10. , 10. , 10. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  9.5,  9.5,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  8. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
            [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
        ])
        adk = np.array([
            [7.  , 9.  , 5.  , 7.  , 3.  , 0.  , 0.  ],
            [8.  , 6.5 , 8.5 , 5.75, 7.75, 0.  , 0.  ],
        ])
        # Insert e as a sibling of d. This should recover tree1.
        itag = 4
        tag = order[itag]
        size = sizes[tag]
        depth = depths[tag]
        _bal_insert_plan(
            n - 2, itag, adm, adk, tree, order, sizes, depths, pairs, diffs, ancs,
            ancx, segs, lvls, oops, False)
        enc = 3
        nc = _bal_avgdist_chunk(
            n - 2, itag, order, sizes, pairs, segs, lvls, oops, enc, chunks, chusegs)
        _bal_update_spine(tag, depth, sizes, pairs, ancs)
        _bal_avgdist_flat(n - 2, itag, size, order, adm, adk, diffs, depths)
        _bal_avgdist_nest(
            itag, depth, obs := adm, order, sizes, depths, diffs, npots, ancs, ancx,
            segs, lvls, nc, chunks, chusegs)
        exp = np.array([
            [ 0.  ,  5.  ,  8.75,  9.  ,  9.  ,  8.5 ,  8.  ],
            [ 0.  ,  0.  ,  9.75, 10.  , 10.  ,  9.5 ,  9.  ],
            [ 0.  ,  0   ,  0.  ,  9.5 ,  9.5 ,  9.  ,  8.5 ],
            [ 0.  ,  0.  ,  0.  ,  0.  ,  8.  ,  7.5 ,  7.  ],
            [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  3.  ],
            [ 0.  ,  0.  ,  0.  ,  0.  ,  8.75,  0.  ,  7.75],
            [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
        ])
        npt.assert_allclose(obs, exp)

    def test_bal_insert_plan(self):
        """Navigate a tree and plan on subsequent parallelization."""
        # This is a relatively large tree with 19 nodes and 10 tips.
        obj = TreeNode.read([self.nwk6])
        n = obj.count()
        m = obj.count(tips=True)

        # Generate tree topology array.
        tree = np.empty((n + 2, 4), dtype=int)
        taxmap = {x.name: int(x.name) for x in obj.tips()}
        _from_treenode(obj, taxmap, tree)

        # Generate tree properties based on the topology.
        order = np.empty(n + 2, dtype=int)
        order[:n] = np.arange(n)
        sizes = np.empty(n + 2, dtype=int)
        _calc_sizes(n, tree, order, sizes)
        depths = np.empty(n + 2, dtype=int)
        _calc_depths(n, tree, order, depths)
        pairs = np.empty(n + 2, dtype=int)
        _calc_pairs(n, tree, order, sizes, pairs)

        # Fill the matrices with zeros, as this test primarily concerns parallelization
        # planning, whereas numeric accuracy was assessed in `test_bal_avgdist_insert_
        # 2`.
        adm = np.zeros((n + 2, n + 2), dtype=float)
        adk = np.zeros((2, n + 2), dtype=float)

        # Allocate arrays that will be filled.
        diffs = np.empty(n + 2, dtype=float)
        ancs = np.empty(n + 2, dtype=int)
        ancx = np.empty(n + 2, dtype=int)
        segs = np.empty(n + 3, dtype=int)
        lvls = np.empty(n + 2, dtype=int)
        oops = np.empty(n + 2, dtype=int)

        # We will insert a taxon into the branch above node 8. This node has a depth of
        # 3, meaning it has three ancestors: 7 (parent), 1 (grandparent) and 0 (root).
        # Meanwhile, 8 itself is the left child of its parent (7), therefore it has a
        # right sibling: 15. Moving up, it has a left cousin 2 and a right cousin 18.
        itag = 8
        tag = order[itag]
        size = sizes[tag]
        depth = depths[tag]

        # do the job
        _bal_insert_plan(n, itag, adm, adk, tree, order, sizes, depths, pairs, diffs,
                         ancs, ancx, segs, lvls, oops, True)

        # ancestors
        npt.assert_array_equal(ancs[:depth], [7, 1, 0])
        npt.assert_array_equal(ancx[:depth], [147, 21, 0])  # == ancs * adm.shape[0]

        # segments (ancestors high-to-low, self, right cousins low-to-high)
        n_segs = (segs == n).argmax()
        self.assertEqual(n_segs, 6)
        npt.assert_array_equal(segs[:n_segs + 1], [0, 1, 7, 8, 15, 18, 19])

        # levels up from parent
        npt.assert_array_equal(lvls[:n_segs], [2, 1, 0, -1, 0, 2])

        # numbers of operations
        npt.assert_array_equal(oops[:n_segs], [22, 18, 6, 4, 2, 2])

        # update sizes and pairs
        _bal_update_spine(tag, depth, sizes, pairs, ancs)

        # Update tree topoloy, then calculate sizes, depths and pairs from scratch, and
        # compare them with the incrementally updated results.
        _insert_taxon(m + 1, itag, size, tree, order)
        exp = np.empty(n + 2, dtype=int)
        _calc_sizes(n + 2, tree, order, exp)
        npt.assert_array_equal(sizes, exp)
        exp = np.empty(n + 2, dtype=int)
        _calc_depths(n + 2, tree, order, exp)
        npt.assert_array_equal(depths, exp)
        exp = np.empty(n + 2, dtype=int)
        _calc_pairs(n + 2, tree, order, sizes, exp)
        npt.assert_array_equal(pairs, exp)

    def test_bal_avgdist_chunk(self):
        """Partition nodes into chunks with roughly even workloads."""
        # Use the same example as above.
        # preparation (which is long...)
        obj = TreeNode.read([self.nwk6])
        n = obj.count()
        tree = np.empty((n + 2, 4), dtype=int)
        taxmap = {x.name: int(x.name) for x in obj.tips()}
        _from_treenode(obj, taxmap, tree)
        order = np.empty(n + 2, dtype=int)
        order[:n] = np.arange(n)
        sizes = np.empty(n + 2, dtype=int)
        _calc_sizes(n, tree, order, sizes)
        depths = np.empty(n + 2, dtype=int)
        _calc_depths(n, tree, order, depths)
        pairs = np.empty(n + 2, dtype=int)
        _calc_pairs(n, tree, order, sizes, pairs)
        adm = np.zeros((n + 2, n + 2), dtype=float)
        adk = np.zeros((2, n + 2), dtype=float)
        diffs = np.empty(n + 2, dtype=float)
        ancs = np.empty(n + 2, dtype=int)
        ancx = np.empty(n + 2, dtype=int)
        segs = np.empty(n + 3, dtype=int)
        lvls = np.empty(n + 2, dtype=int)
        oops = np.empty(n + 2, dtype=int)
        itag = 8
        _bal_insert_plan(n, itag, adm, adk, tree, order, sizes, depths, pairs, diffs,
                         ancs, ancx, segs, lvls, oops, True)

        # chunk boundaries and segment status per chunk
        chunks = np.zeros(n + 1, dtype=int)
        chunks[0] = 0
        chusegs = np.zeros(n, dtype=int)
        chusegs[0] = 0

        # We aim at dividing the tree into 4 chunks, and eventually get 5.
        enc = 4
        obs = _bal_avgdist_chunk(
            n, itag, order, sizes, pairs, segs, lvls, oops, enc, chunks, chusegs)
        self.assertEqual(obs, 5)

        # chunk 0: nodes 0..2: two segs: 0, 1
        # chunk 1: nodes 2..3: one seg: 2 (=1)
        # chunk 2: nodes 3..6: one seg: 3 (=1)
        # chunk 3: nodes 6..15: three segs: 6 (=1), 7, 8
        # chunk 4: nodes 15..19: two segs: 15, 18
        npt.assert_array_equal(chunks[:obs + 1], [0, 2, 3, 6, 15, 19])
        npt.assert_array_equal(chusegs[:obs], [0, 2, 2, 2, 4])

    def test_ols_lengths(self):
        """Calculate tree branch lengths using an OLS framework."""
        # Example 1: This should recover the same branch lengths as in the original
        # tree.
        dm, tree, order, tacts = self.dm1, self.tree1, self.preodr1, self.tacts1
        n = tree.shape[0]
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, tacts)
        _ols_lengths(obs := np.zeros(n), adm, tree, tacts)
        npt.assert_allclose(obs, self.lens1)

        # Example 2: The output is close but not precisely identical to those in the
        # original tree.
        dm, tree, order, tacts = self.dm2, self.tree2, self.preodr2, self.tacts2
        n = tree.shape[0]
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, tacts)
        _ols_lengths(obs := np.zeros(n), adm, tree, tacts)
        exp = np.array([
            0.91769, 0.76891, 0.42026875, 0.35793125, 0.04316597, 0.28054444,
            0.03137847, 0.15226875, 0.04148125, 0.12214, 0.14706])
        npt.assert_allclose(obs, exp)

    def test_ols_lengths_d2(self):
        """Calculate tree branch lengths using an OLS framework."""
        # Test if result matches that calculated from the full matrix.
        # example 1
        dm, tree, order, tacts = self.dm1, self.tree1, self.preodr1, self.tacts1
        n = tree.shape[0]
        ran_ = np.arange(n)
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, tacts)
        ad2 = np.ascontiguousarray(
            np.vstack([adm[ran_, tree[ran_, 3]], adm[ran_, tree[ran_, 2]]]))
        _ols_lengths_d2(obs := np.zeros(n), ad2, tree, tacts)
        _ols_lengths(exp := np.zeros(n), adm, tree, tacts)
        npt.assert_allclose(obs, exp)

        # example 2
        dm, tree, order, tacts = self.dm2, self.tree2, self.preodr2, self.tacts2
        n = tree.shape[0]
        ran_ = np.arange(n)
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, tacts)
        ad2 = np.ascontiguousarray(
            np.vstack([adm[ran_, tree[ran_, 3]], adm[ran_, tree[ran_, 2]]]))
        _ols_lengths_d2(obs := np.zeros(n), ad2, tree, tacts)
        _ols_lengths(exp := np.zeros(n), adm, tree, tacts)
        npt.assert_allclose(obs, exp)

    def test_bal_lengths(self):
        """Calculate tree branch lengths using a balanced framework."""
        # Example 1: In this simple example, OLS and balanced frameworks produce the
        # same branch lengths.
        dm, tree, order, sizes = self.dm1, self.tree1, self.preodr1, self.sizes1
        n = sizes[0]
        _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, sizes)
        _bal_lengths(obs := np.zeros(n), adm, tree)
        npt.assert_allclose(obs, self.lens1)

        # Example 2: Also slightly different from the original branch lengths.
        dm, tree, order, sizes = self.dm2, self.tree2, self.preodr2, self.sizes2
        n = sizes[0]
        _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, sizes)
        _bal_lengths(obs := np.zeros(n), adm, tree)
        exp = np.array([
            0.92853125, 0.75806875, 0.41086875, 0.36733125, 0.04648125, 0.28469375,
            0.02695625, 0.15393125, 0.03981875, 0.11678125, 0.15241875])
        npt.assert_allclose(obs, exp)

    def test_ols_min_branch_d2(self):
        """Find the branch with minimum length change using an OLS framework."""
        # Test if result matches ground truth.
        dm, tree, order, tacts = self.dm1, self.tree1m1, self.preodr1m1, self.tacts1m1
        n = tree.shape[0]
        k = dm.shape[0] - 1
        ran_ = np.arange(n)
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, tacts)
        ad2 = np.ascontiguousarray(np.vstack([
            adm[ran_, tree[ran_, 3]], adm[ran_, tree[ran_, 2]]]))
        _avgdist_taxon(adk := np.zeros((2, n)), k, dm, tree, order, tacts)
        res = _ols_min_branch_d2(obs := np.zeros(n), ad2, adk, tree, order, tacts)
        self.assertEqual(res, 4)
        exp = np.array([0, 0, -2, -2, -3, 0, 0], dtype=float)
        # The algorithm omits factor 0.5, therefore we need to x2 here.
        npt.assert_allclose(obs, 2 * exp)

        # Test if result matches that calculated from the entire tree.
        dm, tree, order, tacts = self.dm3, self.tree3, self.preodr3, self.tacts3
        n = tree.shape[0]
        k = dm.shape[0] - 1
        ran_ = np.arange(n)
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, tacts)
        ad2 = np.ascontiguousarray(
            np.vstack([adm[ran_, tree[ran_, 3]], adm[ran_, tree[ran_, 2]]])
        )
        _avgdist_taxon(adk := np.zeros((2, n)), k, dm, tree, order, tacts)
        res = _ols_min_branch_d2(obs := np.zeros(n), ad2, adk, tree, order, tacts)

        # For each branch, insert taxon, calculate full matrix, then calculate and sum
        # all branch lengths. The difference between each sum and the sum by the root
        # branch is the length change value calculated by the algorithm. 
        exp = np.zeros(n)
        for i in range(n - 2):
            size = tacts[order[i]] * 2 - 1
            _insert_taxon(k, i, size, tree_ := tree.copy(), order_ := order.copy())
            _calc_tacts(n, tree_, order_, tacts_ := np.empty(n, dtype=int))
            _avgdist_matrix(adm := np.zeros((n, n)), dm, tree_, order_, tacts_)
            _ols_lengths(lens := np.zeros(n), adm, tree_, tacts_)
            exp[order[i]] = lens.sum()
        exp[:n - 2] -= exp[0]

        npt.assert_allclose(obs, exp * 2)
        self.assertEqual(order[res], exp[:n - 2].argmin())

    def test_bal_min_branch(self):
        """Find the branch with minimum length change using a balacned framework."""
        # Test if result matches ground truth.
        dm = self.dm1
        tree, order, sizes = self.tree1m1, self.preodr1m1, self.sizes1m1
        n = tree.shape[0]
        k = dm.shape[0] - 1
        _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, sizes)
        _bal_avgdist_taxon(n - 2, k, dm, adk := np.zeros((2, n)), tree, order)
        res = _bal_min_branch(n - 2, obs := np.zeros(n), adm, adk, tree, order)
        self.assertEqual(res, 4)
        exp = np.array([0, 0, -2, -2, -3, 0, 0], dtype=float)
        # The algorithm omits factor 0.25, therefore we need to x4 here.
        npt.assert_allclose(obs, 4 * exp)

        # Test if result matches that calculated from the entire tree.
        dm, tree, order, sizes = self.dm3, self.tree3, self.preodr3, self.sizes3
        n = tree.shape[0]
        k = dm.shape[0] - 1
        _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, sizes)
        _bal_avgdist_taxon(n - 2, k, dm, adk := np.zeros((2, n)), tree, order)
        res = _bal_min_branch(n - 2, obs := np.zeros(n), adm, adk, tree, order)
        exp = np.zeros(n)
        for i in range(n - 2):
            size = sizes[order[i]]
            _insert_taxon(k, i, size, tree_ := tree.copy(), order_ := order.copy())
            _calc_sizes(n, tree_, order_, sizes_ := np.empty(n, dtype=int))
            _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree_, order_, sizes_)
            _bal_lengths(lens := np.zeros(n), adm, tree_)
            exp[order[i]] = lens.sum()
        exp[:n - 2] -= exp[0]

        npt.assert_allclose(obs, 4 * exp)
        self.assertEqual(order[res], exp[:n - 2].argmin())

    def test_swap_branches(self):
        # example 1: swap (e,d) with b
        # swap with tacts (for OLS framework)
        tree, order, index, tacts = (
            self.tree1.copy(), self.preodr1.copy(), self.preidx1.copy(),
            self.tacts1.copy())
        n = tree.shape[0]
        stack = np.full(n, 0)
        _swap_branches(2, 1, tree, order, index, stack, tacts=tacts)

        # make sure tree topology, node order and attributes are valid
        _check_tree(tree)
        _preorder(order_ := np.empty(n, dtype=int), tree, stack)
        npt.assert_array_equal(order, order_)
        index_ = np.empty(n, dtype=int)
        index_[order] = np.arange(n)
        npt.assert_array_equal(index, index_)
        _calc_tacts(n, tree, order, tacts_ := np.empty(n, dtype=int))
        npt.assert_array_equal(tacts, tacts_)

        # check new parent mapping
        self.assertEqual(tree[4, 2], 0)
        self.assertEqual(tree[1, 2], 2)

        # swap with sizes and depths (for balanced framework)
        tree, order, index, sizes, depths = (
            self.tree1.copy(), self.preodr1.copy(), self.preidx1.copy(),
            self.sizes1.copy(), self.depths1.copy())
        _swap_branches(2, 1, tree, order, index, stack, sizes=sizes, depths=depths)

        _check_tree(tree)
        _preorder(order_ := np.empty(n, dtype=int), tree, stack)
        npt.assert_array_equal(order, order_)
        index_ = np.empty(n, dtype=int)
        index_[order] = np.arange(n)
        npt.assert_array_equal(index, index_)
        _calc_sizes(n, tree, order, sizes_ := np.empty(n, dtype=int))
        npt.assert_array_equal(sizes, sizes_)
        _calc_depths(n, tree, order, depths_ := np.empty(n, dtype=int))
        npt.assert_array_equal(depths, depths_)

        self.assertEqual(tree[4, 2], 0)
        self.assertEqual(tree[1, 2], 2)

        # neither tacts nor sizes provided
        with self.assertRaises(ValueError) as cm:
            _swap_branches(2, 1, tree, order, index, stack)
        msg = "Must provide either tacts or sizes."
        self.assertEqual(str(cm.exception), msg)

        # both tacts and sizes provided
        with self.assertRaises(ValueError) as cm:
            _swap_branches(2, 1, tree, order, index, stack,
                           tacts=tacts, sizes=sizes)
        msg = "Cannot provide tacts and sizes simultaneously."
        self.assertEqual(str(cm.exception), msg)

        # example 4: all possible swaps
        taxa, tree, order, index, tacts = (
            self.taxa4, self.tree4, self.preodr4, self.preidx4, self.tacts4)
        n = tree.shape[0]
        stack = np.full(n, 0)

        taxamap = {}
        for node in order[::-1]:
            if tree[node, 0] == 0:
                taxamap[node] = [taxa[tree[node, 1]]]
            else:
                taxamap[node] = taxamap[tree[node, 0]] + taxamap[tree[node, 1]]

        for node in range(1, n):
            if tree[node, 0] == 0:
                continue
            node1 = taxamap[tree[node, 3]]
            for side in range(2):
                tree_, order_, index_, tacts_ = (
                    tree.copy(), order.copy(), index.copy(), tacts.copy())
                node2 = taxamap[tree[node, side]]
                exp = _to_treenode(tree_, taxa)
                _swap_branches_treenode(exp.lca(node1), exp.lca(node2))

                _swap_branches(node, side, tree_, order_, index_, stack, tacts=tacts_)
                _check_tree(tree_)
                obs = _to_treenode(tree_, taxa)

                self.assertEqual(obs.compare_rfd(exp), 0)

    def test_avgdist_swap(self):
        """Update average distance matrix after branch swapping."""
        # Test if the algorithm produces the same result as calculated from the full
        # matrix after branch swapping.

        # example 1: swap (e,d) with b
        dm, tree, order, index, tacts = (
            self.dm1, self.tree1, self.preodr1, self.preidx1, self.tacts1)
        n = tree.shape[0]
        stack = np.full(n, 0)
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, tacts)
        _swap_branches(2, 1, tree, order, index, stack, tacts=tacts)
        _avgdist_swap(index[2], 1, obs := adm.copy(), tree, order, tacts)
        _avgdist_matrix(exp := np.zeros((n, n)), dm, tree, order, tacts)
        npt.assert_allclose(obs, exp)

        # example 4: all possible swaps
        dm, tree, order, index, tacts = (
            self.dm4, self.tree4, self.preodr4, self.preidx4, self.tacts4)
        n = tree.shape[0]
        stack = np.full(n, 0)
        _avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, tacts)
        for node in range(1, n):
            if tree[node, 0] == 0:  # omit tips
                continue
            for side in range(2):  # left and right
                tree_, order_, index_, tacts_ = (
                    tree.copy(), order.copy(), index.copy(), tacts.copy())
                _swap_branches(
                    node, side, tree_, order_, index_, stack, tacts=tacts_)
                _avgdist_swap(
                    index_[node], side, obs := adm.copy(), tree_, order_, tacts_)
                _avgdist_matrix(
                    exp := np.zeros((n, n)), dm, tree_, order_, tacts_)
                npt.assert_allclose(obs, exp)

    def test_bal_avgdist_swap(self):
        """Update balanced average distance matrix after branch swapping."""
        # Test if the algorithm produces the same result as calculated from the full
        # matrix after branch swapping.

        # example 1: swap (e,d) with b
        dm, tree, order, index, sizes, depths = (
            self.dm1, self.tree1, self.preodr1, self.preidx1, self.sizes1, self.depths1
        )
        n = tree.shape[0]
        stack = np.full(n, 0)
        npots = np.ldexp(1.0, -np.arange((n + 1) // 2))
        _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, sizes)
        _swap_branches(2, 1, tree, order, index, stack, sizes=sizes, depths=depths)
        _bal_avgdist_swap(
            obs := adm.copy(), 2, 1, tree, order, index, sizes, depths, npots, stack)
        _bal_avgdist_matrix(exp := np.zeros((n, n)), dm, tree, order, sizes)
        npt.assert_allclose(obs, exp)

        # example 4: all possible swaps
        dm, tree, order, index, sizes, depths = (
            self.dm4, self.tree4, self.preodr4, self.preidx4, self.sizes4, self.depths4
        )
        n = tree.shape[0]
        stack = np.full(n, 0)
        npots = np.ldexp(1.0, -np.arange((n + 1) // 2))
        _bal_avgdist_matrix(adm := np.zeros((n, n)), dm, tree, order, sizes)
        for node in range(1, n):
            if tree[node, 0] == 0:
                continue
            for side in range(2):
                tree_, order_, index_, sizes_, depths_ = (
                    tree.copy(), order.copy(), index.copy(), sizes.copy(),
                    depths.copy())
                _swap_branches(
                    node, side, tree_, order_, index_, stack, sizes=sizes_,
                    depths=depths_)
                _bal_avgdist_swap(
                    obs := adm.copy(), node, side, tree_, order_, index_, sizes_,
                    depths_, npots, stack)
                _bal_avgdist_matrix(exp := np.zeros((n, n)), dm, tree_, order_, sizes_)
                npt.assert_allclose(obs, exp)

    def test_ols_all_swaps(self):
        """Evaluate possible swaps at all branches of a tree."""
        dm, tree, order, index, tacts = (
            self.dm4, self.tree4, self.preodr4, self.preidx4, self.tacts4)
        n = tree.shape[0]
        _avgdist_matrix(adm := np.empty((n, n)), dm, tree, order, tacts)
        _ols_lengths(lens := np.empty(n), adm, tree, tacts)
        lensum = lens.sum()

        # Nodes 2 and 3 will have length reduction.
        sides = np.empty(n, dtype=int)
        _ols_all_swaps(lens, adm, tree, tacts, sides)
        exp = np.array([0, 0, -0.47035, -0.03426528, 0, 0, 0, 0, 0, 0, 0])
        npt.assert_allclose(lens, exp)

        # For each of those two branches, perform branch swapping, update average
        # distance matrix, calculate total tree length, and check if the actual
        # tree length gain equals to the algorithm-calculated value.
        stack = np.full(n, 0)
        for node in np.flatnonzero(lens):
            side = sides[node]
            tree_, order_, index_, tacts_ = (
                tree.copy(), order.copy(), index.copy(), tacts.copy())
            _swap_branches(node, side, tree_, order_, index_, stack, tacts=tacts_)
            _avgdist_swap(
                index_[node], side, adm_ := adm.copy(), tree_, order_, tacts_)
            _ols_lengths(lens_ := np.empty(n), adm_, tree_, tacts_)

            # The algorithm omits factor 0.5, therefore we need to x2 here.
            self.assertAlmostEqual(lens[node], 2 * (lens_.sum() - lensum))

    def test_ols_corner_swaps(self):
        """Update swaps of the four corner branches of a swapped branch."""
        dm, tree, order, index, tacts = (
            self.dm4, self.tree4, self.preodr4, self.preidx4, self.tacts4)
        n = tree.shape[0]
        _avgdist_matrix(adm := np.empty((n, n)), dm, tree, order, tacts)
        sides = np.empty(n, dtype=int)
        _ols_all_swaps(lens := np.empty(n), adm, tree, tacts, sides)
        heap = [(lens[i], i, sides[i]) for i in np.flatnonzero(lens)]
        heapify(heap)

        # Perform swaps at the two branches that are known to produce gains (see
        # `test_ols_all_swaps`), and check if the updated result matches the result
        # re-calculated from the entire tree.
        stack = np.full(n, 0)
        for node in np.flatnonzero(lens):
            side = sides[node]
            tree_, order_, index_, tacts_ = (
                tree.copy(), order.copy(), index.copy(), tacts.copy())
            _swap_branches(node, side, tree_, order_, index_, stack, tacts=tacts_)
            _avgdist_swap(
                index_[node], side, adm_ := adm.copy(), tree_, order_, tacts_)
            _ols_all_swaps(exp := np.empty(n), adm_, tree_, tacts_, sides)
            _ols_corner_swaps(
                node, heap, obs := lens.copy(), adm_, tree_, tacts_, sides)
            obs[node] = 0
            npt.assert_allclose(obs, exp)

    def test_bal_all_swaps(self):
        """Evaluate possible swaps at all branches of a tree."""
        # Using the same strategy as `test_ols_all_swaps`.
        dm, tree, order, index, sizes, depths = (
            self.dm4, self.tree4, self.preodr4, self.preidx4, self.sizes4, self.depths4
        )
        n = tree.shape[0]
        npots = np.ldexp(1.0, -np.arange((n + 1) // 2))
        _bal_avgdist_matrix(adm := np.empty((n, n)), dm, tree, order, sizes)
        _bal_lengths(lens := np.empty(n), adm, tree)
        lensum = lens.sum()

        nodes = np.flatnonzero(tree[1:, 0]) + 1
        nb = nodes.size
        gains = np.zeros(nb, dtype=float)
        sides = np.empty(nb, dtype=int)

        _bal_all_swaps(gains, sides, nodes, adm, tree)
        exp = np.array([0, 0.9407, 0.06885, 0])
        npt.assert_allclose(gains, exp)

        stack = np.full(n, 0)
        for branch in range(nb):
            if gains[branch] == 0:
                continue
            node, side = nodes[branch], sides[branch]
            tree_, order_, index_, sizes_, depths_ = (
                tree.copy(), order.copy(), index.copy(), sizes.copy(), depths.copy())
            _swap_branches(
                node, side, tree_, order_, index_, stack, sizes=sizes_, depths=depths_)
            _bal_avgdist_swap(
                adm_ := adm.copy(), node, side, tree_, order_, index_, sizes_, depths_,
                npots, stack)
            _bal_lengths(lens_ := np.empty(n), adm_, tree_)

            # The algorithm omits factor 0.25, therefore we need to x4 here.
            self.assertAlmostEqual(gains[branch], 4 * (lensum - lens_.sum()))


if __name__ == "__main__":
    main()
