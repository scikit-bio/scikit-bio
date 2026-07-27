# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# line length is useful here, so disabling check
# flake8: noqa: E501

import io
import os
import tempfile
from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio.tree import BPTree
from skbio.tree.bp import parse_newick
import skbio.tree.tests.bp.test_bp_cy as tbc


class BPCythonTests(TestCase):
    def test_cython_level(self):
        # exercise the cdef-level tests in test_bp_cy under the test runner
        # rather than at import time
        for name in dir(tbc):
            if name.startswith('test_'):
                with self.subTest(name=name):
                    getattr(tbc, name)()


class BPTests(TestCase):
    def setUp(self):
        #                       0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
        self.fig1_B = np.array([1, 1, 1, 0, 1, 0, 1, 1 ,0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0], dtype=np.uint8)
        self.bptree = BPTree(self.fig1_B)

    def test_init_unbalanced_raises(self):
        with self.assertRaises(ValueError):
            BPTree(np.array([1, 1, 0], dtype=np.uint8))

    def test_rmq(self):
        #       (  (  (  )  (  )  (  (  )  )   )   (   )   (   (   (   )   (   )   )   )   )
        #excess 1  2  3  2  3  2  3  4  3  2   1   2   1   2   3   4   3   4   3   2   1   0
        #i      0  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  20  21

        exp = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 21],
                  [1, 1, 1, 1, 1, 1, 1, 1, 1, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 21],
                     [2, 3, 3, 3, 3, 3, 3, 3, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 21],
                        [3, 3, 3, 3, 3, 3, 3, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 21],
                           [4, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 21],
                              [5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 21],
                                 [6, 6, 6, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 21],
                                    [7, 8, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 21],
                                       [8, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 21],
                                          [9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 21],
                                             [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 21],
                                                 [11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 21],
                                                     [12, 12, 12, 12, 12, 12, 12, 12, 12, 21],
                                                         [13, 13, 13, 13, 13, 13, 13, 20, 21],
                                                             [14, 14, 14, 14, 14, 19, 20, 21],
                                                                 [15, 16, 16, 16, 19, 20, 21],
                                                                     [16, 16, 16, 19, 20, 21],
                                                                         [17, 18, 19, 20, 21],
                                                                             [18, 19, 20, 21],
                                                                                 [19, 20, 21],
                                                                                     [20, 21],
                                                                                         [21]]
        for i in range(len(self.fig1_B)):
            for j in range(i+1, len(self.fig1_B)):
                self.assertEqual(self.bptree.rmq(i, j), exp[i][j - i])

    def test_rMq(self):
        #       (  (  (  )  (  )  (  (  )  )   )   (   )   (   (   (   )   (   )   )   )   )
        #excess 1  2  3  2  3  2  3  4  3  2   1   2   1   2   3   4   3   4   3   2   1   0
        #i      0  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  20  21

        exp = [[0, 1, 2, 2, 2, 2, 2, 7, 7, 7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7, 7],
                  [1, 2, 2, 2, 2, 2, 7, 7, 7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7, 7],
                     [2, 2, 2, 2, 2, 7, 7, 7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7, 7],
                        [3, 4, 4, 4, 7, 7, 7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7, 7],
                           [4, 4, 4, 7, 7, 7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7, 7],
                              [5, 6, 7, 7, 7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7, 7],
                                 [6, 7, 7, 7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7, 7],
                                    [7, 7, 7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7,  7, 7],
                                       [8, 8,  8,  8,  8,  8,  8, 15, 15, 15, 15, 15, 15, 15],
                                          [9,  9,  9,  9,  9, 14, 15, 15, 15, 15, 15, 15, 15],
                                             [10, 11, 11, 11, 14, 15, 15, 15, 15, 15, 15, 15],
                                                 [11, 11, 11, 14, 15, 15, 15, 15, 15, 15, 15],
                                                     [12, 13, 14, 15, 15, 15, 15, 15, 15, 15],
                                                         [13, 14, 15, 15, 15, 15, 15, 15, 15],
                                                             [14, 15, 15, 15, 15, 15, 15, 15],
                                                                 [15, 15, 15, 15, 15, 15, 15],
                                                                     [16, 17, 17, 17, 17, 17],
                                                                         [17, 17, 17, 17, 17],
                                                                             [18, 18, 18, 18],
                                                                                 [19, 19, 19],
                                                                                     [20, 20],
                                                                                         [21]]
        for i in range(len(self.fig1_B)):
            for j in range(i+1, len(self.fig1_B)):
                self.assertEqual(self.bptree.rMq(i, j), exp[i][j - i])

    def test_mincount(self):
        #       (  (  (  )  (  )  (  (  )  )   )   (   )   (   (   (   )   (   )   )   )   )
        #i      0  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  20  21
        #excess 1  2  3  2  3  2  3  4  3  2   1   2   1   2   3   4   3   4   3   2   1   0

        exp = [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  2,  2,  3,  3,  3,  3,  3,  3,  3,  3,  4, 1],
                  [1, 1, 2, 2, 3, 3, 3, 3, 4,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  3, 1],
                     [1, 1, 1, 2, 2, 2, 2, 3,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  3, 1],
                        [1, 1, 2, 2, 2, 2, 3,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  3, 1],
                           [1, 1, 1, 1, 1, 2,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  3, 1],
                              [1, 1, 1, 1, 2,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  3, 1],
                                 [1, 1, 2, 1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  3, 1],
                                    [1, 1, 1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  3, 1],
                                       [1, 1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  3, 1],
                                          [1,  1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  3, 1],
                                              [1,  1,  2,  2,  2,  2,  2,  2,  2,  2,  3, 1],
                                                  [1,  1,  1,  1,  1,  1,  1,  1,  1,  2, 1],
                                                      [1,  1,  1,  1,  1,  1,  1,  1,  2, 1],
                                                          [1,  1,  1,  1,  1,  1,  2,  1, 1],
                                                              [1,  1,  2,  2,  3,  1,  1, 1],
                                                                  [1,  1,  1,  2,  1,  1, 1],
                                                                      [1,  1,  2,  1,  1, 1],
                                                                          [1,  1,  1,  1, 1],
                                                                              [1,  1,  1, 1],
                                                                                  [1,  1, 1],
                                                                                      [1, 1],
                                                                                         [1]]

        for i in range(len(self.fig1_B)):
            for j in range(i+1, len(self.fig1_B)):
                self.assertEqual(self.bptree.mincount(i, j), exp[i][j - i])

    def test_minselect(self):
        """position of the qth minimum in excess(i), excess(i + 1), . . . , excess(j)."""
        exp = {(0, 20, 1): 0,
               (0, 21, 1): 21,
               (0, 20, 2): 10,
               (0, 21, 2): None,
               (0, 20, 3): 12,
               (0, 20, 4): 20,
               (8, 15, 1): 10,
               (8, 15, 2): 12,
               (6, 9, 1): 9}

        for (i, j, q), e in exp.items():
            self.assertEqual(self.bptree.minselect(i, j, q), e)

    def test_preorder_rank(self):
        exp = [1, 2, 3, 3, 4, 4, 5, 6, 6, 5, 2, 7, 7, 8, 9, 10, 10, 11, 11, 9, 8, 1]
        for i, e in enumerate(exp):
            self.assertEqual(self.bptree.preorder_rank(i), e)

    def test_preorder_select(self):
        exp = [0, 1, 2, 4, 6, 7, 11, 13, 14, 15, 17]
        for k, e in enumerate(exp):
            self.assertEqual(self.bptree.preorder_select(k), e)

    def test_postorder_rank(self):
        exp = [11, 5, 1, 1, 2, 2, 4, 3, 3, 4, 5, 6, 6, 10, 9, 7, 7, 8, 8, 9, 10, 11]
        for i, e in enumerate(exp):
            self.assertEqual(self.bptree.postorder_rank(i), e)

    def test_postorder_select(self):
        exp = [2, 4, 7, 6, 1, 11, 15, 17, 14, 13, 0]
        for k, e in enumerate(exp):
            self.assertEqual(self.bptree.postorder_select(k + 1), e)

    def test_is_ancestor(self):
        exp = {(0, 0): False,  # identity test
               (2, 1): False,  # tip test
               (1, 2): True,   # open test
               (1, 3): True,   # close test
               (0, 7): True,   # nested test
               (1, 7): True}   # nested test

        for (i, j), e in exp.items():
            self.assertEqual(self.bptree.is_ancestor(i, j), e)

    def test_count(self):
        # node counts per subtree (count(i) == former subtree_size(i))
        exp = [11, 5, 1, 1, 1, 1, 2, 1, 1, 2, 5, 1, 1, 4, 3, 1, 1, 1, 1, 3, 4, 11]
        for i, e in enumerate(exp):
            self.assertEqual(self.bptree.count(i), e)

        # the index defaults to the root, i.e. the whole tree
        self.assertEqual(self.bptree.count(), exp[0])

        # tip counts per subtree (count(i, tips=True))
        exp_tips = [6, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 2, 2, 1, 1, 1, 1, 2,
                    2, 6]
        for i, e in enumerate(exp_tips):
            self.assertEqual(self.bptree.count(i, tips=True), e)

        # whole-tree tip count (former ntips)
        self.assertEqual(self.bptree.count(tips=True), 6)

    def test_level_ancestor(self):
        exp = {(2, 1): 1,  # first tip to its parent
               (2, 2): 0,  # first tip to root
               (4, 1): 1,  # second tip to its parent
               (5, 1): 1,  # second tip, closing, to its parent
               (7, 1): 6,  # deep tip to its parent
               (7, 2): 1,  # deep tip to its grandparent
               (7, 3): 0,  # deep tip to its great grand parent
               (7, 9999): 0,  # max out at the root
               (10, 0): -1}  # can't be an ancestor of yourself

        for (i, d), e in exp.items():
            self.assertEqual(self.bptree.level_ancestor(i, d), e)

    def _testinator(self, exp, f, verbose=False):
        self.assertEqual(len(exp), len(self.fig1_B))
        for i, e in enumerate(exp):
            if verbose:
                print(i, e)
            self.assertEqual(f(i), e)

    def test_level_next(self):
        #       (   (  (  )  (  )   (   (   )   )   )  (    )   (   (   (   )   (   )   )   )   )
        exp = [-1, 11, 4, 4, 6, 6, 14, 15, 15, 14, 11, 13, 13, -1, -1, 17, 17, -1, -1, -1, -1, -1]
        self.assertEqual(len(exp), len(self.fig1_B))

        for i, e in enumerate(exp):
            self.assertEqual(self.bptree.level_next(i), e)

    def test_close(self):
        exp = [21, 10, 3, 5, 9, 8, 12, 20, 19, 16, 18]
        for i, e in zip(np.argwhere(self.bptree.data == 1).squeeze(), exp):
            npt.assert_equal(self.bptree.close(i), e)

    def test_lca(self):
        # lca(i, j) = parent(rmq(i, j) + 1)
        # unless is_ancestor(i, j)
        # (so lca(i, j) = i) or is_ancestor(j, i) (so lca(i, j) = j),
        nodes = [self.bptree.preorder_select(k) for k in range(self.fig1_B.sum())]
        exp = {(nodes[2], nodes[3]): nodes[1],
               (nodes[2], nodes[5]): nodes[1],
               (nodes[2], nodes[9]): nodes[0],
               (nodes[9], nodes[10]): nodes[8],
               (nodes[1], nodes[8]): nodes[0]}
        for (i, j), e in exp.items():
            self.assertEqual(self.bptree.lca(i, j), e)

    def test_deepest_node(self):
        # deepest_node(i) = rMq(i, close(i)),
        exp = [7, 7, 2, 2, 4, 4, 7, 7, 7, 7, 7, 11, 11, 15, 15, 15, 15, 17, 17, 15, 15, 7]
        self._testinator(exp, self.bptree.deepest_node)

    def test_height(self):
        # height(i) = excess(deepest_node(i)) - excess(i).
        exp = [3, 2, 0, 0, 0, 0, 1, 0, 0, 1, 2, 0, 0, 2, 1, 0, 0, 0, 0, 1, 2, 3]
        self._testinator(exp, self.bptree.height)

    def test_shear(self):
        #       r  2  3     4     5  6             7       8   9  10      11
        #       (  (  (  )  (  )  (  (  )  )   )   (   )   (   (   (   )   (   )   )   )   )
        #i      0  1  2  3  4  5  6  7  8  9  10  11  12  13  14  15  16  17  18  19  20  21
        names = np.array(['r', '2', '3', None, '4', None, '5', '6', None, None, None, '7', None, '8', '9', '10', None,
                          '11', None, None, None, None])
        lengths = np.array([0, 1, 2, 0, 3, 0, 4, 5, 0, 0, 0, 6, 0, 7, 8, 9, 0, 10, 0, 0, 0, 0], dtype=np.double)
        self.bptree.set_names(names)
        self.bptree.set_lengths(lengths)

        in_ = {'4', '6', '7', '10', '11'}
        exp = np.array([1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 0,
                        0, 0], dtype=np.uint32)
        exp_n = np.array(['r', '2', '4', None, '5', '6', None, None, None, '7', None, '8', '9', '10', None, '11', None,
                          None, None, None])
        exp_l = np.array([0, 1, 3, 0, 4, 5, 0, 0, 0, 6, 0, 7, 8, 9, 0, 10, 0, 0, 0, 0], dtype=np.double)
        obs = self.bptree.shear(in_)
        npt.assert_equal(exp, obs.data)

        for i in range(len(obs.data)):
            self.assertEqual(obs.name(i), exp_n[i])
            self.assertEqual(obs.length(i), exp_l[i])

        in_ = {'10', '11'}
        exp = np.array([1, 1, 1, 1, 0, 1, 0, 0, 0, 0], dtype=np.uint32)
        obs = self.bptree.shear(in_).data
        npt.assert_equal(obs, exp)

    def test_shear_raise_tree_is_empty(self):
        names = np.array(['r', '2', '3', None, '4', None, '5', '6', None, None, None, '7', None, '8', '9', '10', None,
                          '11', None, None, None, None])
        lengths = np.array([0, 1, 2, 0, 3, 0, 4, 5, 0, 0, 0, 6, 0, 7, 8, 9, 0, 10, 0, 0, 0, 0], dtype=np.double)
        self.bptree.set_names(names)
        with self.assertRaises(ValueError):
            self.bptree.shear({'not', 'in', 'tree'})

    def test_collapse(self):
        names = np.array(['r', '2', '3', None, '4', None, '5', '6', None, None, None, '7', None, '8', '9', '10', None,
                          '11', None, None, None, None])
        lengths = np.array([0, 1, 2, 0, 3, 0, 4, 5, 0, 0, 0, 6, 0, 7, 8, 9, 0, 10, 0, 0, 0, 0], dtype=np.double)
        self.bptree.set_names(names)
        self.bptree.set_lengths(lengths)

        exp = np.array([1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0],
                       dtype=np.uint8)
        exp_n = ['r', '2', '3', None, '4', None, '6', None, None, '7', None, '9', '10', None, '11', None, None, None]
        exp_l = [0, 1, 2, 0, 3, 0, 9, 0, 0, 6, 0, 15, 9, 0, 10, 0, 0, 0]

        obs = self.bptree.collapse()

        npt.assert_equal(obs.data, exp)
        for i in range(len(obs.data)):
            self.assertEqual(obs.name(i), exp_n[i])
            self.assertEqual(obs.length(i), exp_l[i])

        bp = BPTree(np.array([1, 1, 1, 0, 0, 1, 0, 0], dtype=np.uint8))
        exp = np.array([1, 1, 0, 1, 0, 0])
        obs = bp.collapse().data

        npt.assert_equal(obs, exp)

    def test_name_unset(self):
        for i in range(self.bptree.data.size):
            self.assertEqual(self.bptree.name(i), None)

    def test_length_unset(self):
        for i in range(self.bptree.data.size):
            self.assertEqual(self.bptree.length(i), 0.0)

    def test_name_length_set(self):
        names = np.full(self.bptree.data.size, None, dtype=object)
        lengths = np.zeros(self.bptree.data.size, dtype=np.double)

        names[0] = 'root'
        names[self.bptree.preorder_select(7)] = 'other'

        lengths[1] = 1.23
        lengths[self.bptree.preorder_select(5)] = 5.43

        self.bptree.set_names(names)
        self.bptree.set_lengths(lengths)

        self.assertEqual(self.bptree.name(0), 'root')
        self.assertEqual(self.bptree.name(1), None)
        self.assertEqual(self.bptree.name(13), 'other')
        self.assertEqual(self.bptree.length(1), 1.23)
        self.assertEqual(self.bptree.length(5), 0.0)
        self.assertEqual(self.bptree.length(7), 5.43)

    def test_write_read_roundtrip(self):
        # give the tree names and lengths so the round-trip has data to preserve
        names = np.array(['r', '2', '3', None, '4', None, '5', '6', None, None,
                          None, '7', None, '8', '9', '10', None, '11', None,
                          None, None, None])
        lengths = np.array([0, 1, 2, 0, 3, 0, 4, 5, 0, 0, 0, 6, 0, 7, 8, 9, 0,
                            10, 0, 0, 0, 0], dtype=np.double)
        self.bptree.set_names(names)
        self.bptree.set_lengths(lengths)

        # round-trip through a file-like object
        buf = io.BytesIO()
        self.bptree.write(buf)
        buf.seek(0)
        obs = BPTree.read(buf)

        # topology, names, and lengths are preserved (edges are not stored)
        npt.assert_equal(obs.data, self.bptree.data)
        for i in range(obs.data.size):
            self.assertEqual(obs.name(i), self.bptree.name(i))
            self.assertEqual(obs.length(i), self.bptree.length(i))

        # and the reconstructed tree answers queries identically
        self.assertEqual(obs.count(tips=True), self.bptree.count(tips=True))

    def test_write_read_roundtrip_file(self):
        names = np.array(['r', '2', '3', None, '4', None, '5', '6', None, None,
                          None, '7', None, '8', '9', '10', None, '11', None,
                          None, None, None])
        lengths = np.array([0, 1, 2, 0, 3, 0, 4, 5, 0, 0, 0, 6, 0, 7, 8, 9, 0,
                            10, 0, 0, 0, 0], dtype=np.double)
        self.bptree.set_names(names)
        self.bptree.set_lengths(lengths)

        # np.savez_compressed appends ".npz" to a bare filename, so the path
        # must already carry that suffix to round-trip via read()
        fd, path = tempfile.mkstemp(suffix='.npz')
        os.close(fd)
        try:
            self.bptree.write(path)
            obs = BPTree.read(path)
        finally:
            os.remove(path)

        npt.assert_equal(obs.data, self.bptree.data)
        for i in range(obs.data.size):
            self.assertEqual(obs.name(i), self.bptree.name(i))
            self.assertEqual(obs.length(i), self.bptree.length(i))


if __name__ == '__main__':
    main()
