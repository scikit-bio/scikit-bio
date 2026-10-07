# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

# line length is useful here, so disabling check
# flake8: noqa: E501

import contextlib
import copy
import io
import os
import pickle
import tempfile
import types
import warnings
from unittest import TestCase, main, mock

import array_api_compat as aac
import numpy as np
import numpy.testing as npt

from skbio._base import SkbioObject
from skbio._config import _resolve_engine, get_config, set_config
from skbio.tree import BPTree, TreeNode
from skbio.tree._exception import DuplicateNodeError, MissingNodeError
from skbio.tree.bp import parse_newick
from skbio.tree.bp import _bp_cy, _bp_gpu, _gpu
from skbio.tree.bp._bp import _build_index, _KERNEL_METHODS
from skbio.tree.bp._bp_numba import NUMBA_AVAILABLE
from skbio.stats.distance import _gpu as _stats_gpu
from skbio.util._array import _get_backend_name, _to_numpy
from skbio.util._testing import ArrayAPITestMixin, array_backends, numba_code
import skbio.tree.tests.bp.test_bp_cy as tbc


def _random_topology(n_nodes, rng):
    """Parentheses of a random tree with ``n_nodes`` nodes."""
    # attach each new node under a uniformly chosen earlier node, then emit the
    # parentheses by an iterative depth-first walk
    children = [[] for _ in range(n_nodes)]
    for v in range(1, n_nodes):
        children[rng.integers(v)].append(v)
    out = []
    stack = [(0, False)]
    while stack:
        v, done = stack.pop()
        if done:
            out.append(0)
            continue
        out.append(1)
        stack.append((v, True))
        stack.extend((c, False) for c in reversed(children[v]))
    return np.array(out, dtype=np.uint8)


def _caterpillar(n_tips):
    out = [1, 1, 0] * (n_tips - 1) + [1, 0] + [0] * (n_tips - 1)
    return np.array(out, dtype=np.uint8)


def _naive_lca(B):
    """Lowest common ancestors by parent pointers, independent of the engines.

    The reference for ``lca`` and ``lca_batch``. Checking those against each
    other proves nothing, as both run the same kernel. Returns a function of
    two positions -- either parenthesis of each node, in either order -- giving
    the opening position of their lowest common ancestor.
    """
    B = np.asarray(B)
    node = np.empty(B.size, dtype=np.intp)  # the opening position of each node
    parent = {}
    stack = []
    for p in range(B.size):
        if B[p]:
            parent[p] = stack[-1] if stack else -1
            stack.append(p)
            node[p] = p
        else:
            node[p] = stack.pop()

    def lca(i, j):
        a, b = node[i], node[j]
        ancestors = set()
        while a != -1:
            ancestors.add(a)
            a = parent[a]
        while b not in ancestors:
            b = parent[b]
        return b

    return lca


def _all_topologies(max_nodes):
    """Every tree of 1 to ``max_nodes`` nodes, as a parenthesis array."""
    out = []

    def grow(seq, n_open, n_close, nodes):
        if n_open == nodes and n_close == nodes:
            out.append(np.array(seq, dtype=np.uint8))
            return
        if n_open < nodes and (n_open == n_close == 0 or n_open > n_close):
            grow(seq + [1], n_open + 1, n_close, nodes)
        # a single root: the excess may reach 0 only at the very end
        if n_close < n_open and (n_close + 1 < n_open or n_open == nodes):
            grow(seq + [0], n_open, n_close + 1, nodes)

    for nodes in range(1, max_nodes + 1):
        grow([], 0, 0, nodes)
    return out


def _index_test_topologies():
    rng = np.random.default_rng(42)
    trees = [np.array([1, 0], dtype=np.uint8),  # a single node
             np.array([1, 1, 0, 0], dtype=np.uint8),
             np.array([1] + [1, 0] * 50 + [0], dtype=np.uint8)]  # a star
    trees += [_caterpillar(n) for n in (2, 3, 17, 500)]
    # sizes straddle rmM block and level boundaries
    trees += [_random_topology(n, rng) for n in
              (3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 100, 257, 1000, 4097,
               20000)]
    return trees


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

    def test_minselect_rank_below_one(self):
        # there is no 0th or earlier minimum (not the last ones, as a negative
        # index would give)
        for q in (0, -1, -5):
            self.assertIsNone(self.bptree.minselect(0, 20, q))

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

    def test_lca_same_node(self):
        # a node is its own lowest common ancestor, at either parenthesis (the
        # root's closing one is the last position)
        opening, stack = [], []
        for i, bit in enumerate(self.fig1_B):
            if bit:
                stack.append(i)
                opening.append(i)
            else:
                opening.append(stack.pop())
        for i, exp in enumerate(opening):
            self.assertEqual(self.bptree.lca(i, i), exp)

    def test_two_node_tree(self):
        # (()): a root with a single tip. Its rmM blocks used to hold one
        # parenthesis each, too few for the backward search, which then missed
        # position 0: the root's last child came back as none, and the tip's
        # closing parenthesis resolved to the root.
        bp = BPTree(np.array([1, 1, 0, 0], dtype=np.uint8))
        self.assertEqual(bp.last_child(0), 1)
        self.assertEqual(bp.last_child(3), 1)
        self.assertEqual(bp.first_child(2), 0)  # a tip has no children
        self.assertEqual(bp.count(2), 1)
        self.assertEqual(bp.height(2), 0)
        self.assertEqual(bp.preorder_rank(2), 2)
        self.assertFalse(bp.is_ancestor(2, 0))
        self.assertEqual(bp.lca(2, 2), 1)

    def test_lca_any_parenthesis_and_order(self):
        # Either parenthesis of a node names it, the order of the two does not
        # matter, and the answer is always an opening parenthesis. Checked
        # against parent pointers: every pair of a small tree, else a sample.
        # Before, lca(j, i) with j > i could return the wrong node, and a
        # closing parenthesis could come back as the answer.
        rng = np.random.default_rng(13)
        for B in _index_test_topologies():
            if B.size > 5000:
                continue
            bp = BPTree(B)
            ref = _naive_lca(B)
            if B.size <= 70:
                i, j = np.meshgrid(np.arange(B.size), np.arange(B.size))
                i, j = i.ravel(), j.ravel()
            else:
                i, j = rng.integers(0, B.size, (2, 3000))
            for a, b in zip(i.tolist(), j.tolist()):
                self.assertEqual(bp.lca(a, b), ref(a, b), (B.size, a, b))

    def test_small_trees_exhaustively(self):
        # Every tree of up to 7 nodes, where the rmM blocks are at their
        # smallest and every boundary case is reachable, against a stack scan.
        # Any parenthesis may name a node. Two bugs this caught:
        # - in (()) the block size was 1 and the backward search missed
        #   position 0, so e.g. the root's last_child came back as none;
        # - is_ancestor(close(x), x) was True: a node, named by its closing
        #   parenthesis, counted as its own ancestor.
        topologies = _all_topologies(7)
        self.assertEqual(len(topologies), 1 + 1 + 2 + 5 + 14 + 42 + 132)
        for B in topologies:
            bp = BPTree(B)
            ref = _naive_lca(B)
            stack, node, parent, children = [], {}, {}, {}
            for p, bit in enumerate(B):
                if bit:
                    parent[p] = stack[-1] if stack else -1
                    children.setdefault(parent[p], []).append(p)
                    stack.append(p)
                    node[p] = p
                else:
                    node[p] = stack.pop()
            close = {node[p]: p for p in range(B.size) if not B[p]}

            def ancestors(x):
                out = set()
                while parent[x] != -1:
                    x = parent[x]
                    out.add(x)
                return out

            tag = ''.join('(' if bit else ')' for bit in B)
            for p in range(B.size):
                x = node[p]
                kids = children.get(x, [])
                if B[p]:
                    self.assertEqual(bp.close(p), close[x], (tag, p))
                    self.assertEqual(bp.parent(p), parent[x], (tag, p))
                self.assertEqual(bp.first_child(p), kids[0] if kids else 0,
                                 (tag, p))
                self.assertEqual(bp.last_child(p), kids[-1] if kids else 0,
                                 (tag, p))
                for q in range(B.size):
                    self.assertEqual(bp.lca(p, q), ref(p, q), (tag, p, q))
                    self.assertEqual(bool(bp.is_ancestor(p, q)),
                                     x in ancestors(node[q]), (tag, p, q))

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
        # BPTree.read/write delegate to the skbio.io registry; the default
        # format is newick. Give the tree names and lengths so the round-trip
        # has data to preserve.
        names = np.array(['r', '2', '3', None, '4', None, '5', '6', None, None,
                          None, '7', None, '8', '9', '10', None, '11', None,
                          None, None, None])
        lengths = np.array([0, 1, 2, 0, 3, 0, 4, 5, 0, 0, 0, 6, 0, 7, 8, 9, 0,
                            10, 0, 0, 0, 0], dtype=np.double)
        self.bptree.set_names(names)
        self.bptree.set_lengths(lengths)

        # round-trip through a file-like object
        buf = io.StringIO()
        self.bptree.write(buf)
        buf.seek(0)
        obs = BPTree.read(buf)

        # topology, names, and lengths are preserved
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

        fd, path = tempfile.mkstemp(suffix='.nwk')
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

    def test_to_from_npz_roundtrip(self):
        # to_npz/from_npz are the binary serialization path, separate from the
        # registry-backed read/write. They preserve topology, names, and
        # lengths (edge numbers are not stored).
        names = np.array(['r', '2', '3', None, '4', None, '5', '6', None, None,
                          None, '7', None, '8', '9', '10', None, '11', None,
                          None, None, None])
        lengths = np.array([0, 1, 2, 0, 3, 0, 4, 5, 0, 0, 0, 6, 0, 7, 8, 9, 0,
                            10, 0, 0, 0, 0], dtype=np.double)
        self.bptree.set_names(names)
        self.bptree.set_lengths(lengths)

        buf = io.BytesIO()
        self.bptree.to_npz(buf)
        buf.seek(0)
        obs = BPTree.from_npz(buf)

        npt.assert_equal(obs.data, self.bptree.data)
        for i in range(obs.data.size):
            self.assertEqual(obs.name(i), self.bptree.name(i))
            self.assertEqual(obs.length(i), self.bptree.length(i))
        self.assertEqual(obs.count(tips=True), self.bptree.count(tips=True))


class BPIndexTests(TestCase):
    """The navigation index matches the scan-based construction."""

    def assert_index_equal(self, obs, exp):
        for key in ('e_index', 'k_index_0', 'k_index_1', 'm', 'M', 'r'):
            npt.assert_array_equal(obs[key], exp[key], err_msg=key)
        self.assertEqual(obs['b'], exp['b'])
        self.assertEqual(obs['height'], exp['height'])

    def test_build_index_matches_reference(self):
        for B in _index_test_topologies():
            with self.subTest(size=B.size):
                self.assert_index_equal(_build_index(B), tbc.reference_index(B))

    def test_tree_index_matches_reference(self):
        # the arrays held by the tree (and viewed by the kernel) are the index
        for B in _index_test_topologies()[:8]:
            bp = BPTree(B)
            exp = tbc.reference_index(B)
            for key in ('e_index', 'k_index_0', 'k_index_1', 'm', 'M', 'r'):
                obs = getattr(bp, '_' + key)
                npt.assert_array_equal(obs, exp[key], err_msg=key)
                self.assertEqual(obs.dtype, np.intp)
                self.assertFalse(obs.flags.writeable)



class BPTreeClassTests(TestCase):
    """``BPTree`` is a Python class delegating navigation to its kernel."""

    def setUp(self):
        self.B = np.array([1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0,
                           1, 0, 0, 0, 0], dtype=np.uint8)
        self.names = np.array(['r', '2', '3', None, '4', None, '5', '6', None,
                               None, None, '7', None, '8', '9', '10', None,
                               '11', None, None, None, None])
        self.lengths = np.arange(self.B.size, dtype=np.double)
        self.bp = BPTree(self.B, lengths=self.lengths, names=self.names)

    def test_skbio_object_subclass(self):
        self.assertTrue(issubclass(BPTree, SkbioObject))
        self.assertIn(SkbioObject, BPTree.__mro__)
        self.assertEqual(str(self.bp), repr(self.bp))

    def test_navigation_bound_to_kernel(self):
        # instance lookups resolve to the compiled kernel's bound methods
        for name in _KERNEL_METHODS:
            meth = getattr(self.bp, name)
            self.assertIs(meth.__self__, self.bp._kernel, name)

    def test_class_methods_match_bound(self):
        # the documented class-level methods give the same answers
        for i in range(self.B.size - 1):
            for name in ('close', 'parent', 'depth', 'is_tip', 'first_child',
                         'last_child', 'next_sibling', 'previous_sibling',
                         'preorder_rank', 'postorder_rank', 'deepest_node',
                         'height', 'level_next', 'name', 'length', 'edge',
                         'count'):
                self.assertEqual(getattr(BPTree, name)(self.bp, i),
                                 getattr(self.bp, name)(i), (name, i))
        self.assertEqual(BPTree.lca(self.bp, 2, 7), self.bp.lca(2, 7))
        self.assertEqual(BPTree.count(self.bp, 0, True),
                         self.bp.count(0, tips=True))

    def test_return_types(self):
        self.assertIs(type(self.bp.length(1)), float)
        self.assertIs(type(self.bp.edge(1)), int)
        self.assertIs(type(self.bp.depth(1)), int)
        self.assertIs(type(self.bp.name(1)), str)
        self.assertIs(type(len(self.bp)), int)

    def test_data_read_only_attribute(self):
        npt.assert_array_equal(self.bp.data, self.B)
        with self.assertRaises(AttributeError):
            self.bp.data = self.B

    def test_single_node(self):
        bp = BPTree(np.array([1, 0], dtype=np.uint8))
        self.assertEqual(len(bp), 1)
        self.assertEqual(bp.close(0), 1)
        self.assertEqual(bp.parent(0), -1)
        self.assertEqual(bp.count(tips=True), 1)
        self.assertTrue(bp.is_tip(0))

    def test_invalid_input(self):
        with self.assertRaises(TypeError):
            BPTree([1, 0])
        with self.assertRaises(ValueError):
            BPTree(self.B.astype(np.int64))
        with self.assertRaises(ValueError):
            BPTree(np.array([], dtype=np.uint8))
        with self.assertRaises(ValueError):
            BPTree(self.B, names=self.names[:-1])
        with self.assertRaises(ValueError):
            BPTree(self.B, lengths=self.lengths.astype(np.float32))
        with self.assertRaises(ValueError):
            self.bp.set_edges(np.full(self.B.size, self.B.size, dtype=np.int32))
        with self.assertRaises(TypeError):
            self.bp.shear(['4'])

    def test_bool_topology(self):
        # accepted as uint8, as the compiled class always did (e.g. trees
        # saved from one built from a bool array)
        B = self.B.astype(bool)
        bp = BPTree(B, lengths=self.lengths, names=self.names)
        self.assertEqual(bp.data.dtype, np.uint8)
        self.assertTrue(np.shares_memory(bp.data, B))
        npt.assert_array_equal(bp.data, self.B)
        for i in range(self.B.size):
            self.assertEqual(bp.close(i), self.bp.close(i))
            self.assertEqual(bp.parent(i), self.bp.parent(i))
        # non-contiguous too
        wide = np.zeros(2 * self.B.size, dtype=bool)
        wide[::2] = B
        npt.assert_array_equal(BPTree(wide[::2]).data, self.B)

        # a saved tree whose topology is bool: .npz and pickle
        with io.BytesIO() as fh:
            np.savez_compressed(fh, names=self.names, lengths=self.lengths, B=B)
            fh.seek(0)
            obs = BPTree.from_npz(fh)
        npt.assert_array_equal(obs.data, self.B)
        self.assertEqual(obs.name(2), self.bp.name(2))

        class Saved:  # pickles as a tree pickled with a bool topology does
            def __reduce__(self_):
                return BPTree, (B, self.lengths, self.names)

        obs = pickle.loads(pickle.dumps(Saved()))
        self.assertIsInstance(obs, BPTree)
        npt.assert_array_equal(obs.data, self.B)

    def test_non_contiguous_topology(self):
        wide = np.zeros(2 * self.B.size, dtype=np.uint8)
        wide[::2] = self.B
        bp = BPTree(wide[::2])
        npt.assert_array_equal(bp.data, self.B)
        self.assertEqual(bp.close(1), 10)

    def test_set_attributes_update_kernel(self):
        names = np.full(self.B.size, None, dtype=object)
        names[1] = 'x'
        self.bp.set_names(names)
        self.assertEqual(self.bp.name(1), 'x')

        lengths = np.zeros(self.B.size, dtype=np.double)
        lengths[1] = 2.5
        self.bp.set_lengths(lengths)
        self.assertEqual(self.bp.length(1), 2.5)

        edges = np.arange(self.B.size, dtype=np.int32)
        self.bp.set_edges(edges)
        self.assertEqual(self.bp.edge(4), 4)
        self.assertEqual(self.bp.edge_from_number(4), 4)

    def test_pickle_and_copy(self):
        for obs in (pickle.loads(pickle.dumps(self.bp)), copy.copy(self.bp),
                    copy.deepcopy(self.bp)):
            self.assertIsNot(obs._kernel, self.bp._kernel)
            npt.assert_array_equal(obs.data, self.bp.data)
            for i in range(self.B.size):
                self.assertEqual(obs.name(i), self.bp.name(i))
                self.assertEqual(obs.length(i), self.bp.length(i))
                self.assertEqual(obs.close(i), self.bp.close(i))

    def test_io_descriptors(self):
        self.assertIn('newick', BPTree.read.__doc__)
        self.assertIn('newick', self.bp.write.__doc__)


def _named_tree(n_nodes, rng):
    """Random tree with named tips and random branch lengths."""
    B = _random_topology(n_nodes, rng)
    tip = np.zeros(B.size, dtype=bool)
    tip[:-1] = (B[:-1] == 1) & (B[1:] == 0)
    names = np.full(B.size, None, dtype=object)
    names[tip] = ['t%d' % k for k in range(tip.sum())]
    lengths = rng.random(B.size)
    return BPTree(B, lengths=lengths, names=names)


class _BatchTests:
    """Batch operations of one compute engine (``engine``)."""

    engine = None

    def setUp(self):
        self.rng = np.random.default_rng(11)
        self.trees = [BPTree(np.array([1, 0], dtype=np.uint8)),
                      BPTree(_caterpillar(30))]
        self.trees += [_named_tree(n, self.rng) for n in (2, 7, 64, 1500)]
        # a named caterpillar: the deepest tree for its size
        B = _caterpillar(200)
        names = np.full(B.size, None, dtype=object)
        names[np.flatnonzero(B[:-1] > B[1:])] = ['c%d' % k for k in range(200)]
        self.trees.append(
            BPTree(B, lengths=self.rng.random(B.size), names=names))

    def test_close_parent_batch(self):
        for bp in self.trees:
            pos = np.arange(bp.data.size - 1)
            npt.assert_array_equal(
                bp.close_batch(pos, engine=self.engine),
                [bp.close(int(i)) for i in pos])
            npt.assert_array_equal(
                bp.parent_batch(pos, engine=self.engine),
                [bp.parent(int(i)) for i in pos])

    def test_lca_batch(self):
        for bp in self.trees:
            n = bp.data.size
            # any position (closing parentheses too), pairs in either order;
            # the reference is independent of the engines (see _naive_lca)
            i = self.rng.integers(0, n, 500)
            j = self.rng.integers(0, n, 500)
            ref = _naive_lca(bp.data)
            exp = [ref(a, b) for a, b in zip(i.tolist(), j.tolist())]
            npt.assert_array_equal(bp.lca_batch(i, j, engine=self.engine), exp)
            npt.assert_array_equal(bp.lca_batch(j, i, engine=self.engine), exp)
            # and the single-node method gives the same answers
            npt.assert_array_equal(
                [bp.lca(a, b) for a, b in zip(i.tolist(), j.tolist())], exp)

    def test_level_ancestor_batch(self):
        for bp in self.trees:
            pos = np.arange(bp.data.size - 1)
            d = self.rng.integers(-1, 5, pos.size)
            npt.assert_array_equal(
                bp.level_ancestor_batch(pos, d, engine=self.engine),
                [bp.level_ancestor(int(a), int(b)) for a, b in zip(pos, d)])
            # a scalar level broadcasts
            npt.assert_array_equal(
                bp.level_ancestor_batch(pos, 1, engine=self.engine),
                [bp.level_ancestor(int(a), 1) for a in pos])

    def test_shapes(self):
        bp = self.trees[-1]
        pos = np.flatnonzero(bp.data)[:12].reshape(3, 4)
        obs = bp.parent_batch(pos, engine=self.engine)
        self.assertEqual(obs.shape, (3, 4))
        self.assertEqual(obs.dtype, np.intp)
        npt.assert_array_equal(obs.ravel(),
                               bp.parent_batch(pos.ravel(), engine=self.engine))
        # broadcasting a column against a row
        obs = bp.lca_batch(pos[:, :1], pos[0], engine=self.engine)
        self.assertEqual(obs.shape, (3, 4))
        # a list and an empty input
        self.assertEqual(
            bp.close_batch([0], engine=self.engine).tolist(), [bp.close(0)])
        self.assertEqual(bp.close_batch([], engine=self.engine).shape, (0,))

    def test_invalid_positions(self):
        bp = self.trees[-1]
        with self.assertRaises(IndexError):
            bp.close_batch([bp.data.size], engine=self.engine)
        with self.assertRaises(IndexError):
            bp.lca_batch([0], [-1], engine=self.engine)
        with self.assertRaises(TypeError):
            bp.parent_batch([1.5], engine=self.engine)
        with self.assertRaises(TypeError):
            bp.level_ancestor_batch([1], [1.0], engine=self.engine)

    def test_cophenet_matches_treenode(self):
        for bp in self.trees[2:]:
            self.assertGreater(bp.count(tips=True), 0)
            tn = TreeNode.from_bptree(bp)
            for use_length in (True, False):
                obs = bp.cophenet(use_length=use_length, engine=self.engine)
                exp = tn.cophenet(use_length=use_length)
                self.assertEqual(obs.ids, exp.ids)
                npt.assert_allclose(obs.data, exp.data, rtol=1e-12, atol=1e-12)

    def test_cophenet_endpoints(self):
        bp = self.trees[-2]
        tn = TreeNode.from_bptree(bp)
        endpoints = ['t5', 't0', 't30', 't2']
        obs = bp.cophenet(endpoints, engine=self.engine)
        exp = tn.cophenet(endpoints)
        self.assertEqual(obs.ids, tuple(endpoints))
        npt.assert_allclose(obs.data, exp.data, rtol=1e-12, atol=1e-12)

        # a large shuffled subset, on the deep tree as well
        for bp in (self.trees[-2], self.trees[-1]):
            tn = TreeNode.from_bptree(bp)
            names = [tip.name for tip in tn.tips()]
            endpoints = list(self.rng.permutation(names)[: len(names) // 2])
            obs = bp.cophenet(endpoints, engine=self.engine)
            exp = tn.cophenet(endpoints)
            self.assertEqual(obs.ids, exp.ids)
            npt.assert_allclose(obs.data, exp.data, rtol=1e-12, atol=1e-12)

    def test_cophenet_errors(self):
        bp = BPTree.read(["((a:1,b:2)c:3,(d:4,e:5)f:6)root;"])
        with self.assertRaises(MissingNodeError):
            bp.cophenet(['a', 'x'], engine=self.engine)
        with self.assertRaises(DuplicateNodeError):
            bp.cophenet(['a', 'b', 'a'], engine=self.engine)
        with self.assertRaises(ValueError):
            bp.cophenet(['a', 'c'], engine=self.engine)
        dup = BPTree.read(["((a:1,a:2)c:3,d:4)root;"])
        with self.assertRaises(DuplicateNodeError):
            dup.cophenet(engine=self.engine)

    def test_cophenet_few_tips(self):
        B = np.array([1, 1, 0, 1, 0, 0], dtype=np.uint8)
        obs = BPTree(B).cophenet(engine=self.engine)
        self.assertEqual(obs.shape, (0, 0))
        names = np.array([None, 'a', None, None, None, None], dtype=object)
        obs = BPTree(B, names=names).cophenet(engine=self.engine)
        self.assertEqual(obs.shape, (1, 1))
        self.assertEqual(obs.ids, ('a',))


class BPBatchCythonTests(_BatchTests, TestCase):
    engine = 'cython'

    def test_engine_resolution(self):
        bp = self.trees[-1]
        with self.assertRaises(ValueError):
            bp.close_batch([0], engine='julia')
        # "fast" and the global default resolve to an available engine
        npt.assert_array_equal(bp.close_batch([0], engine='fast'),
                               [bp.close(0)])
        npt.assert_array_equal(bp.close_batch([0]), [bp.close(0)])

    def test_cophenet_engine_checked_for_few_tips(self):
        # with fewer than two tips there are no pairs to compute, but the
        # engine is still checked, as it is for larger trees
        B = np.array([1, 1, 0, 1, 0, 0], dtype=np.uint8)
        named = np.array([None, 'a', None, None, None, None], dtype=object)
        for names in (None, named):
            with self.assertRaises(ValueError):
                BPTree(B, names=names).cophenet(engine='julia')


@numba_code
class BPBatchNumbaTests(_BatchTests, TestCase):
    engine = 'numba'

    def test_engines_identical(self):
        # the engines run the same algorithm: identical results, bit for bit
        for bp in self.trees:
            n = bp.data.size
            i = self.rng.integers(0, n, 300)
            j = self.rng.integers(0, n, 300)
            npt.assert_array_equal(bp.lca_batch(i, j, engine='numba'),
                                   bp.lca_batch(i, j, engine='cython'))
            for use_length in (True, False):
                npt.assert_array_equal(
                    bp.cophenet(use_length=use_length, engine='numba').data,
                    bp.cophenet(use_length=use_length, engine='cython').data)

    def test_global_engine(self):
        bp = self.trees[-1]
        default = get_config('compute_engine')
        try:
            set_config('compute_engine', 'numba')
            npt.assert_array_equal(bp.close_batch([0]), [bp.close(0)])
        finally:
            set_config('compute_engine', default)


def _engines():
    try:
        import numba  # noqa: F401
    except ImportError:
        return ('cython', 'fast')
    return ('cython', 'numba', 'fast')


class BPDeviceTreeTests(TestCase, ArrayAPITestMixin):
    """A tree of a non-NumPy array (e.g. on a GPU) matches the NumPy tree.

    With ``SKBIO_DEVICE=cuda`` and CuPy or a CUDA PyTorch, the batch operations
    and ``cophenet`` with ``engine="numba"`` run on the GPU.
    """

    def setUp(self):
        rng = np.random.default_rng(7)
        self.host = [_named_tree(n, rng) for n in (2, 9, 300)]
        B = _caterpillar(60)
        names = np.full(B.size, None, dtype=object)
        names[np.flatnonzero(B[:-1] > B[1:])] = ['c%d' % k for k in range(60)]
        self.host.append(BPTree(B, lengths=rng.random(B.size), names=names))
        self.rng = rng

    def on(self, xp, device, bp, dtype=None):
        B = self.make_array(xp, device, bp.data, dtype=dtype or xp.uint8)
        return BPTree(B, lengths=bp._lengths, names=bp._names)

    def assert_ran_on_gpu(self, bp, xp, device):
        """``engine="numba"`` ran the GPU kernels, rather than falling back.

        A fallback gives the same results, so without this a GPU run whose
        kernels cannot compile (e.g. numba-cuda on NumPy 2.5) would
        pass while testing the CPU. It applies to a CUDA device, with a backend
        Numba can read (CuPy, PyTorch) and Numba installed.
        """
        name = _get_backend_name(xp)
        if device in (None, "cpu") or name not in ("cupy", "torch"):
            return
        if not NUMBA_AVAILABLE:
            return
        self.assertIsNotNone(
            bp._gpu_module("numba"),
            "No Numba GPU module for a %s tree on %s: is numba-cuda-mlir (or "
            "numba-cuda) installed, "
            "and can it use this GPU?" % (name, device))
        self.assertNotIn(
            name, _gpu._unavailable,
            "The BPTree GPU kernels failed on %s and fell back to the CPU." % name)
        self.assertIn(name, _bp_gpu._kernels,
                      "The BPTree GPU kernels were never built for %s." % name)
        self.assertEqual(bp._fast_engine(), "numba",
                         "engine='fast' does not choose the GPU for a %s tree." % name)

    @array_backends("numpy", "jax", "torch", "cupy")
    def test_construction(self, xp, device):
        for host in self.host:
            bp = self.on(xp, device, host)
            self.assert_type_preserved(bp.data, xp, device)
            npt.assert_array_equal(_to_numpy(bp.data), host.data)
            exp = tbc.reference_index(host.data)
            for key in ('e_index', 'k_index_0', 'k_index_1', 'm', 'M', 'r'):
                npt.assert_array_equal(getattr(bp, '_' + key), exp[key])
            # every single-node method on the host copy
            for i in range(host.data.size - 1):
                for name in ('close', 'parent', 'depth', 'is_tip', 'name',
                             'length', 'first_child', 'next_sibling',
                             'preorder_rank', 'height', 'count'):
                    self.assertEqual(getattr(bp, name)(i),
                                     getattr(host, name)(i), (name, i))
        # a bool array is held as uint8
        bp = self.on(xp, device, self.host[1], dtype=xp.bool)
        self.assertEqual(bp.data.dtype, xp.uint8)
        self.assert_type_preserved(bp.data, xp, device)
        npt.assert_array_equal(_to_numpy(bp.data), self.host[1].data)

    @array_backends("numpy", "jax", "torch", "cupy")
    def test_invalid_input(self, xp, device):
        with self.assertRaises(ValueError):
            BPTree(self.make_array(xp, device, [1, 1, 0], dtype=xp.uint8))
        with self.assertRaises(ValueError):
            BPTree(self.make_array(xp, device, [1, 0], dtype=xp.int32))

    @array_backends("numpy", "jax", "torch", "cupy")
    def test_batch(self, xp, device):
        for host in self.host:
            bp = self.on(xp, device, host)
            n = host.data.size
            pos = np.arange(n - 1)
            i = self.rng.integers(0, n, 400)
            j = self.rng.integers(0, n, 400)
            d = self.rng.integers(-1, 5, pos.size)
            dev = lambda a: self.make_array(xp, device, a, dtype=xp.int64)
            for engine in _engines():
                cases = [
                    ('close_batch', (pos,), host.close_batch(pos)),
                    ('parent_batch', (pos,), host.parent_batch(pos)),
                    ('lca_batch', (i, j), host.lca_batch(i, j)),
                    ('level_ancestor_batch', (pos, d),
                     host.level_ancestor_batch(pos, d)),
                ]
                for name, args, exp in cases:
                    with self.subTest(engine=engine, op=name, n=n):
                        # queries on the device: a result on the device
                        obs = getattr(bp, name)(*map(dev, args), engine=engine)
                        self.assert_type_preserved(obs, xp, device)
                        npt.assert_array_equal(_to_numpy(obs), exp)
                        # host queries: a host result
                        obs = getattr(bp, name)(*args, engine=engine)
                        self.assertIsInstance(obs, np.ndarray)
                        npt.assert_array_equal(obs, exp)
            self.assert_ran_on_gpu(bp, xp, device)

    @array_backends("numpy", "jax", "torch", "cupy")
    def test_batch_shapes_and_errors(self, xp, device):
        host = self.host[2]
        bp = self.on(xp, device, host)
        pos = np.flatnonzero(host.data)[:12].reshape(3, 4)
        for engine in _engines():
            obs = bp.lca_batch(self.make_array(xp, device, pos[:, :1],
                                               dtype=xp.int64),
                               pos[0], engine=engine)
            self.assertEqual(tuple(obs.shape), (3, 4))
            self.assert_type_preserved(obs, xp, device)
            npt.assert_array_equal(_to_numpy(obs),
                                   host.lca_batch(pos[:, :1], pos[0]))
            empty = self.make_array(xp, device, [], dtype=xp.int64)
            self.assertEqual(tuple(bp.close_batch(empty, engine=engine).shape),
                             (0,))
            with self.assertRaises(IndexError):
                bp.close_batch(self.make_array(xp, device, [host.data.size],
                                               dtype=xp.int64), engine=engine)
            with self.assertRaises(TypeError):
                bp.parent_batch(self.make_array(xp, device, [1.5]),
                                engine=engine)
            with self.assertRaises(TypeError):
                bp.level_ancestor_batch(
                    self.make_array(xp, device, [1], dtype=xp.int64),
                    self.make_array(xp, device, [1.0]), engine=engine)
        self.assert_ran_on_gpu(bp, xp, device)

    @array_backends("numpy", "jax", "torch", "cupy")
    def test_cophenet(self, xp, device):
        for host in self.host[1:]:
            bp = self.on(xp, device, host)
            names = [host.name(int(p)) for p in np.flatnonzero(host.data)
                     if host.is_tip(int(p))]
            subset = list(self.rng.permutation(names)[: len(names) // 2 + 1])
            for engine in _engines():
                for endpoints in (None, subset):
                    for use_length in (True, False):
                        exp = host.cophenet(endpoints, use_length=use_length)
                        obs = bp.cophenet(endpoints, use_length=use_length,
                                          engine=engine)
                        self.assertEqual(obs.ids, exp.ids)
                        self.assert_type_preserved(obs.data, xp, device)
                        # the same algorithm on every device: bit-identical
                        npt.assert_array_equal(_to_numpy(obs.data), exp.data)
            self.assert_ran_on_gpu(bp, xp, device)
        # no pairs
        obs = self.on(xp, device, BPTree(np.array([1, 1, 0, 0], np.uint8)))
        self.assertEqual(obs.cophenet().shape, (0, 0))

    @array_backends("numpy", "jax", "torch", "cupy")
    def test_to_device(self, xp, device):
        # the device as the backend names it (e.g. a cupy.cuda.Device)
        target = aac.device(self.make_array(xp, device, [0], dtype=xp.uint8))
        for host in self.host[1:]:
            bp = host.to_device(xp, target)
            self.assert_type_preserved(bp.data, xp, device)
            npt.assert_array_equal(_to_numpy(bp.data), host.data)
            # the host index and attributes are shared, not rebuilt
            for key in ('_data', '_e_index', '_k_index_0', '_k_index_1', '_m',
                        '_M', '_r', '_names', '_lengths', '_edges'):
                self.assertIs(getattr(bp, key), getattr(host, key), key)
            self.assertIsNot(bp._kernel, host._kernel)
            for i in range(host.data.size - 1):
                self.assertEqual(bp.parent(i), host.parent(i))
                self.assertEqual(bp.name(i), host.name(i))
            n = host.data.size
            i = self.rng.integers(0, n, 300)
            j = self.rng.integers(0, n, 300)
            dev = lambda a: self.make_array(xp, device, a, dtype=xp.int64)
            for engine in _engines():
                obs = bp.lca_batch(dev(i), dev(j), engine=engine)
                self.assert_type_preserved(obs, xp, device)
                npt.assert_array_equal(_to_numpy(obs), host.lca_batch(i, j))
                obs = bp.cophenet(engine=engine)
                self.assert_type_preserved(obs.data, xp, device)
                npt.assert_array_equal(_to_numpy(obs.data), host.cophenet().data)
            self.assert_ran_on_gpu(bp, xp, device)
            # and back to the host, from where the tree lives now
            back = bp.to_device(np)
            self.assertIsInstance(back.data, np.ndarray)
            self.assertIs(back.data, host.data)
            npt.assert_array_equal(back.lca_batch(i, j), host.lca_batch(i, j))

        # setting an attribute of one tree leaves the other as it was (on a
        # tree of its own: the fixtures are shared by the backends' subtests)
        h = self.host[2]
        host = BPTree(h.data, lengths=h._lengths.copy(), names=h._names.copy())
        bp = host.to_device(xp, target)
        names = np.full(host.data.size, 'x', dtype=object)
        bp.set_names(names)
        self.assertEqual(bp.name(0), 'x')
        self.assertNotEqual(host.name(0), 'x')
        lengths = np.full(host.data.size, 9.0)
        host.set_lengths(lengths)
        self.assertEqual(host.length(1), 9.0)
        self.assertNotEqual(bp.length(1), 9.0)
        # the default device of the backend: the current GPU for CuPy
        if _get_backend_name(xp) == 'cupy':
            self.assert_type_preserved(host.to_device(xp).data, xp, device)

    def test_to_device_invalid_backend(self):
        bp = self.host[1]
        for xp in (object(), np.ndarray, 'numpy'):
            with self.assertRaisesRegex(TypeError, 'not a supported array backend'):
                bp.to_device(xp)

    @array_backends("numpy", "jax", "torch", "cupy")
    def test_non_default_stream(self, xp, device):
        # queries made, cast and read on a stream of their own: PyTorch does
        # not tell Numba which stream produced a tensor, so the kernels must
        # wait for it (see _gpu._sync_stream)
        name = _get_backend_name(xp)
        if device in (None, "cpu") or name not in ("cupy", "torch"):
            self.skipTest("needs a CuPy or PyTorch array on a GPU")
        delay = contextlib.nullcontext()
        if name == "torch":
            import torch
            import array_api_compat.torch as compat
            stream = torch.cuda.Stream()
            context = torch.cuda.stream(stream)
            # widen the window of the race: keep the stream busy before each
            # int32 cast, so that a kernel not waiting for it reads it unfinished
            cast = compat.astype

            def slow_cast(x, dtype, *args, **kwargs):
                if x.dtype == torch.int32:
                    torch.cuda._sleep(100_000_000)
                return cast(x, dtype, *args, **kwargs)

            delay = mock.patch.object(compat, "astype", slow_cast)
        else:
            stream = xp.cuda.Stream(non_blocking=True)
            context = stream
        host = self.host[2]
        bp = self.on(xp, device, host)
        n = host.data.size
        pos = self.rng.integers(0, n, 5000)
        d = self.rng.integers(0, 6, pos.size)
        # warm: the index is on the device, so nothing synchronizes before
        # the launch but the kernels' own wait
        bp.level_ancestor_batch(pos[:2], d[:2], engine="numba")
        for _ in range(3):
            with delay, context:
                obs = bp.level_ancestor_batch(
                    self.make_array(xp, device, pos, dtype=xp.int64),
                    self.make_array(xp, device, d, dtype=xp.int32),
                    engine="numba")
                lca = bp.lca_batch(
                    self.make_array(xp, device, pos, dtype=xp.int32),
                    self.make_array(xp, device, pos[::-1].copy(), dtype=xp.int32),
                    engine="numba")
                dm = bp.cophenet(engine="numba")
            npt.assert_array_equal(_to_numpy(obs), host.level_ancestor_batch(pos, d))
            npt.assert_array_equal(_to_numpy(lca),
                                   host.lca_batch(pos, pos[::-1].copy()))
            npt.assert_array_equal(_to_numpy(dm.data), host.cophenet().data)
        self.assert_ran_on_gpu(bp, xp, device)

    @array_backends("numpy", "jax", "torch", "cupy")
    def test_write(self, xp, device):
        # the writers read the host copy, not ``data``, which may be on a GPU
        for host in self.host[1:]:
            edges = np.arange(host.data.size, dtype=np.int32)
            host = BPTree(host.data, lengths=host._lengths, names=host._names,
                          edges=edges)
            bp = BPTree(self.make_array(xp, device, host.data, dtype=xp.uint8),
                        lengths=host._lengths, names=host._names, edges=edges)
            for fmt in ('newick', 'jplace'):
                exp, obs = io.StringIO(), io.StringIO()
                host.write(exp, format=fmt)
                bp.write(obs, format=fmt)
                self.assertEqual(obs.getvalue(), exp.getvalue(), fmt)

    @array_backends("numpy", "jax", "torch", "cupy")
    def test_host_operations(self, xp, device):
        host = self.host[2]
        bp = self.on(xp, device, host)
        # new trees and exports are made from the host copy
        for tree in (pickle.loads(pickle.dumps(bp)), copy.deepcopy(bp),
                     bp.shear({'t0', 't1', 't2'}), bp.collapse()):
            self.assertIsInstance(tree.data, np.ndarray)
        npt.assert_array_equal(pickle.loads(pickle.dumps(bp)).data, host.data)
        self.assertEqual(bp.to_array()['name'].tolist(),
                         host.to_array()['name'].tolist())


@numba_code
class BPGPUFallbackTests(TestCase):
    """A GPU kernel that cannot run falls back to the CPU engines."""

    def setUp(self):
        self.bp = _named_tree(40, np.random.default_rng(2))
        self.key = 'numpy'
        _gpu._unavailable.discard(self.key)
        # a Numba GPU module lacking everything the kernels need
        broken = types.SimpleNamespace(__name__='broken')
        self.bp._gpu_module = lambda engine: broken if engine == 'numba' else None

    def tearDown(self):
        _gpu._unavailable.discard(self.key)

    def test_fallback_warns_once(self):
        bp = self.bp
        n = bp.data.size
        i = np.arange(n)
        j = i[::-1].copy()
        exp = bp.lca_batch(i, j, engine='cython')
        with self.assertWarnsRegex(UserWarning, 'using the CPU engines'):
            obs = bp.lca_batch(i, j, engine='numba')
        npt.assert_array_equal(obs, exp)
        self.assertIn(self.key, _gpu._unavailable)
        # the permutation tests' kernels on the backend are unaffected
        self.assertNotIn('numpy', _stats_gpu._unavailable)
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            npt.assert_array_equal(
                bp.parent_batch(i[:-1], engine='numba'),
                bp.parent_batch(i[:-1], engine='cython'))
            npt.assert_array_equal(
                bp.cophenet(engine='numba').data,
                bp.cophenet(engine='cython').data)

    def test_cophenet_fallback(self):
        with self.assertWarns(UserWarning):
            obs = self.bp.cophenet(engine='numba')
        npt.assert_array_equal(obs.data, self.bp.cophenet(engine='cython').data)

    def test_other_engines_stay_on_host(self):
        # the GPU is only reached for engine="numba"
        with warnings.catch_warnings():
            warnings.simplefilter('error')
            self.bp.lca_batch([1], [2], engine='cython')
            self.bp.cophenet(engine='cython')
        self.assertNotIn(self.key, _gpu._unavailable)


class BPGPUModuleTests(TestCase):
    """Which Numba GPU module the tree's kernels compile through."""

    def setUp(self):
        try:
            import torch
        except ImportError:
            self.skipTest("PyTorch is not installed.")
        if torch.version.hip is not None:
            self.skipTest("a ROCm build of PyTorch maps to numba.hip.")
        # a CPU tensor reports the same namespace as a CUDA one
        self.arr = torch.zeros(2)
        self.mlir = types.SimpleNamespace(__name__='numba_cuda_mlir.cuda')
        self.cuda = types.SimpleNamespace(__name__='numba.cuda')

    def probe(self, installed):
        def available(package, name):
            return installed.get(package + '.' + name)

        with mock.patch.object(_gpu, '_available_module', side_effect=available):
            return _gpu._numba_gpu_module_for(self.arr)

    def test_prefers_numba_cuda_mlir(self):
        installed = {'numba_cuda_mlir.cuda': self.mlir, 'numba.cuda': self.cuda}
        self.assertIs(self.probe(installed), self.mlir)

    def test_falls_back_to_numba_cuda(self):
        self.assertIs(self.probe({'numba.cuda': self.cuda}), self.cuda)

    def test_none_without_an_extension(self):
        self.assertIsNone(self.probe({}))

    def test_none_after_a_failure(self):
        installed = {'numba_cuda_mlir.cuda': self.mlir}
        _gpu._unavailable.add('torch')
        try:
            self.assertIsNone(self.probe(installed))
        finally:
            _gpu._unavailable.discard('torch')


class BPGPUStreamTests(TestCase):
    """A GPU backend's array is read only after its stream's queued work."""

    def test_sync_before_reading_in_place(self):
        calls = []
        arr = types.SimpleNamespace(__cuda_array_interface__={})
        gpu = types.SimpleNamespace(
            as_cuda_array=lambda a: calls.append('read') or a)
        with mock.patch.object(_bp_gpu, '_sync_stream',
                               side_effect=lambda a: calls.append('sync')):
            self.assertIs(_bp_gpu._on_device(gpu, arr), arr)
        self.assertEqual(calls, ['sync', 'read'])

    def test_host_array_is_uploaded_without_sync(self):
        gpu = types.SimpleNamespace(to_device=lambda a: ('uploaded', a))
        with mock.patch.object(_bp_gpu, '_sync_stream') as sync:
            res = _bp_gpu._on_device(gpu, np.arange(3))
        sync.assert_not_called()
        self.assertEqual(res[0], 'uploaded')

    def test_other_arrays_are_ignored(self):
        # e.g. a Numba device array, which no array API namespace takes
        _gpu._sync_stream(types.SimpleNamespace(__cuda_array_interface__={}))
        _gpu._sync_stream(np.arange(3))


class BPFastEngineTests(TestCase):
    """What engine="fast" resolves to: see the Notes of BPTree."""

    def setUp(self):
        self.bp = _named_tree(30, np.random.default_rng(3))

    def resolve(self, openmp, numba, gpu=None):
        """The fast engine of self.bp with the build and install patched."""
        with mock.patch.object(_bp_cy, 'OPENMP', openmp), \
                mock.patch('skbio.tree.bp._bp.NUMBA_AVAILABLE', numba), \
                mock.patch.object(self.bp, '_gpu_module', return_value=gpu):
            return self.bp._fast_engine()

    def test_openmp_flag(self):
        self.assertIsInstance(_bp_cy.OPENMP, bool)

    def test_rules(self):
        gpu = types.SimpleNamespace(__name__='numba.cuda')
        # 1. a tree on a GPU the Numba kernels can use
        self.assertEqual(self.resolve(True, True, gpu), 'numba')
        self.assertEqual(self.resolve(False, True, gpu), 'numba')
        # 2. multithreaded Cython
        self.assertEqual(self.resolve(True, True), 'cython')
        self.assertEqual(self.resolve(True, False), 'cython')
        # 3. serial Cython, and Numba installed
        self.assertEqual(self.resolve(False, True), 'numba')
        # 4. serial Cython, nothing else
        self.assertEqual(self.resolve(False, False), 'cython')

    def test_host_tree_reaches_no_gpu(self):
        # a NumPy tree is never on a GPU: the rule is decided by the build
        self.assertIsNone(self.bp._gpu_module('numba'))
        exp = 'cython' if _bp_cy.OPENMP or not NUMBA_AVAILABLE else 'numba'
        self.assertEqual(self.bp._fast_engine(), exp)

    def test_fast_runs_the_chosen_engine(self):
        bp = self.bp
        exp = bp._fast_engine()
        resolved = []

        def record(*args, **kwargs):
            engine = _resolve_engine(*args, **kwargs)
            resolved.append(engine)
            return engine

        with mock.patch('skbio.tree.bp._bp._resolve_engine', side_effect=record):
            bp.lca_batch([1, 2], [3, 4], engine='fast')
            bp.cophenet(engine='fast')
            default = get_config('compute_engine')
            try:
                set_config('compute_engine', 'fast')
                bp.parent_batch([1])
            finally:
                set_config('compute_engine', default)
        self.assertTrue(resolved)
        self.assertEqual(set(resolved), {exp})


if __name__ == '__main__':
    main()
