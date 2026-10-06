# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import subprocess
import sys
import textwrap
from unittest import TestCase, main

import numpy as np

from skbio.tree import BPTree
from skbio.tree.bp import _bp_numba
from skbio.tree.bp._bp import _KERNEL_METHODS
from skbio.tree.tests.bp import test_bp_cy as tbc
from skbio.tree.tests.bp.test_bp import _index_test_topologies, _random_topology
from skbio.util._testing import numba_code


# Per-position operations with the same signature on BPTree and the primitives.
_UNARY = ("close", "depth", "parent", "is_tip", "first_child", "last_child",
          "next_sibling", "previous_sibling", "preorder_rank", "postorder_rank",
          "level_next", "deepest_node", "height")

# Operations on pairs of positions ``i <= j``.
_PAIRWISE = ("rmq", "rMq", "lca", "is_ancestor", "mincount")

# Index operations reached through the Cython test module's oracle.
_INDEX = ("open", "enclose", "excess")


class PrimitiveCoverageTests(TestCase):
    def test_every_navigation_operation(self):
        # every navigation method of BPTree has a primitive; node attributes are
        # array reads (names are Python objects) and need none
        attributes = {"name", "length", "edge", "edge_from_number"}
        methods = set(_KERNEL_METHODS) - attributes | {"mincount", "minselect"}
        self.assertLessEqual(methods, set(_bp_numba.Primitives._fields))

    def test_bp_arrays(self):
        bp = BPTree(_random_topology(50, np.random.default_rng(0)))
        T = _bp_numba.bp_arrays(bp)
        self.assertIs(T.B, bp.data)
        self.assertEqual(T.n_internal, (1 << T.height) - 1)
        self.assertEqual(T.size, bp.data.size)
        # a transfer function is applied to every array, not to the geometry
        T = _bp_numba.bp_arrays(bp, asarray=np.array)
        self.assertIsNot(T.B, bp.data)
        np.testing.assert_array_equal(T.e_index, bp._e_index)
        self.assertEqual(T.b, bp._b)


@numba_code
class NumbaPrimitiveTests(TestCase):
    """Each primitive returns what the Cython engine returns, at every input."""

    @classmethod
    def setUpClass(cls):
        cls.P = _bp_numba.CPU
        # every index-test topology except the largest (the primitives are
        # called one at a time from Python here)
        cls.trees = [BPTree(B) for B in _index_test_topologies() if B.size < 10000]
        cls.rng = np.random.default_rng(3)

    def pairs(self, n):
        """All pairs ``i <= j`` of a small tree, else a sample."""
        if n <= 80:
            i, j = np.triu_indices(n)
        else:
            i, j = np.sort(self.rng.integers(0, n, (2, 3000)), axis=0)
        return zip(i.tolist(), j.tolist())

    def check(self, obs, exp, *args):
        self.assertEqual(int(obs), int(exp), args)

    def test_unary(self):
        P = self.P
        for bp in self.trees:
            T = _bp_numba.bp_arrays(bp)
            for name in _UNARY:
                prim, meth = getattr(P, name), getattr(bp, name)
                for i in range(bp.data.size):
                    self.check(prim(T, i), meth(i), name, i)
            self.assertEqual(P.root(T), bp.root())

    def test_count(self):
        P = self.P
        for bp in self.trees:
            T = _bp_numba.bp_arrays(bp)
            for i in range(bp.data.size):
                for tips in (False, True):
                    self.check(P.count(T, i, tips), bp.count(i, tips=tips), i, tips)

    def test_level_ancestor(self):
        P = self.P
        for bp in self.trees:
            T = _bp_numba.bp_arrays(bp)
            for i in range(bp.data.size):
                for d in (-1, 0, 1, 2, int(self.rng.integers(3, 8))):
                    self.check(P.level_ancestor(T, i, d),
                               bp.level_ancestor(i, d), i, d)

    def test_select(self):
        P = self.P
        for bp in self.trees:
            T = _bp_numba.bp_arrays(bp)
            n = len(bp)
            for k in range(n):
                self.check(P.preorder_select(T, k), bp.preorder_select(k), k)
                self.check(P.select(T, 1, k), tbc.kernel_index_op(bp, "select", 1, k))
            # the closing-parenthesis select index counts from 1
            for k in range(n + 1):
                self.check(P.postorder_select(T, k), bp.postorder_select(k), k)
                self.check(P.select(T, 0, k), tbc.kernel_index_op(bp, "select", 0, k))

    def test_index_operations(self):
        P = self.P
        for bp in self.trees:
            T = _bp_numba.bp_arrays(bp)
            for i in range(bp.data.size):
                for name in _INDEX:
                    self.check(getattr(P, name)(T, i),
                               tbc.kernel_index_op(bp, name, i), name, i)
                for t in (0, 1):
                    self.check(P.rank(T, t, i),
                               tbc.kernel_index_op(bp, "rank", t, i), t, i)
                for d in (-2, -1, 0, 1, 2):
                    self.check(P.fwdsearch(T, i, d),
                               tbc.kernel_index_op(bp, "fwdsearch", i, d), i, d)
                    self.check(P.bwdsearch(T, i, d),
                               tbc.kernel_index_op(bp, "bwdsearch", i, d), i, d)

    def test_pairwise(self):
        P = self.P
        for bp in self.trees:
            T = _bp_numba.bp_arrays(bp)
            for i, j in self.pairs(bp.data.size):
                for name in _PAIRWISE:
                    self.check(getattr(P, name)(T, i, j),
                               getattr(bp, name)(i, j), name, i, j)
                self.check(P.is_ancestor(T, j, i), bp.is_ancestor(j, i), j, i)
                # lca takes its pair in either order
                self.check(P.lca(T, j, i), bp.lca(j, i), "lca", j, i)

    def test_minselect(self):
        P = self.P
        for bp in self.trees:
            T = _bp_numba.bp_arrays(bp)
            for i, j in self.pairs(bp.data.size):
                # ranks count from 1: nothing precedes the first
                for q in (-1, 0, 1, 2, 3):
                    exp = bp.minselect(i, j, q)
                    self.check(P.minselect(T, i, j, q), -1 if exp is None else exp,
                               i, j, q)

    def test_compiled_caller(self):
        # the primitives compose inside compiled code, their intended use
        from numba import njit

        lca, depth = self.P.lca, self.P.depth

        @njit
        def lca_depths(T, i, j):
            out = np.empty(i.shape[0], dtype=np.intp)
            for t in range(i.shape[0]):
                out[t] = depth(T, lca(T, i[t], j[t]))
            return out

        bp = self.trees[-1]
        i, j = self.rng.integers(0, bp.data.size, (2, 500))
        exp = [bp.depth(bp.lca(a, b)) for a, b in zip(i.tolist(), j.tolist())]
        np.testing.assert_array_equal(lca_depths(_bp_numba.bp_arrays(bp), i, j), exp)


# Runs in a fresh interpreter: the CUDA simulator must be enabled before Numba
# is first imported. It evaluates the device-function build of the primitives in
# a ``cuda.jit`` kernel, one thread per position, against the CPU build.
_SIMULATOR_SCRIPT = textwrap.dedent("""
    import sys
    import numpy as np
    try:
        from numba import cuda
    except ImportError:
        print("SKIP"); sys.exit(0)

    from skbio.tree import BPTree
    from skbio.tree.bp import _bp_numba
    from skbio.tree.tests.bp.test_bp import _random_topology

    G = _bp_numba.gpu_primitives(cuda)
    assert _bp_numba.gpu_primitives(cuda) is G  # built once per module
    COLUMNS = 31

    @cuda.jit
    def evaluate(T, n, out):
        i = cuda.grid(1)
        if i >= T.size:
            return
        j = T.size - 1 - i
        lo, hi = min(i, j), max(i, j)
        row = (
            G.close(T, i), G.open(T, i), G.enclose(T, i), G.excess(T, i),
            G.rank(T, 0, i), G.rank(T, 1, i), G.select(T, 1, i % n),
            G.select(T, 0, i % (n + 1)), G.fwdsearch(T, i, -1),
            G.bwdsearch(T, i, 0), G.rmq(T, lo, hi), G.rMq(T, lo, hi),
            G.mincount(T, lo, hi), G.minselect(T, lo, hi, 1), G.root(T),
            G.depth(T, i), G.parent(T, i), G.is_tip(T, i), G.first_child(T, i),
            G.last_child(T, i), G.next_sibling(T, i), G.previous_sibling(T, i),
            G.preorder_rank(T, i), G.preorder_select(T, i % n),
            G.postorder_rank(T, i), G.postorder_select(T, i % (n + 1)),
            G.is_ancestor(T, lo, hi), G.count(T, i, True),
            G.level_ancestor(T, i, 1), G.lca(T, i, j), G.height(T, i),
        )
        for c in range(len(row)):
            out[i, c] = row[c]

    bp = BPTree(_random_topology(40, np.random.default_rng(5)))
    n = len(bp)
    T = _bp_numba.bp_arrays(bp, asarray=cuda.to_device)
    out = np.zeros((bp.data.size, COLUMNS), dtype=np.intp)
    evaluate[1, 128](T, n, out)

    C, H = _bp_numba.CPU, _bp_numba.bp_arrays(bp)
    for i in range(bp.data.size):
        j = bp.data.size - 1 - i
        lo, hi = min(i, j), max(i, j)
        exp = (
            C.close(H, i), C.open(H, i), C.enclose(H, i), C.excess(H, i),
            C.rank(H, 0, i), C.rank(H, 1, i), C.select(H, 1, i % n),
            C.select(H, 0, i % (n + 1)), C.fwdsearch(H, i, -1),
            C.bwdsearch(H, i, 0), C.rmq(H, lo, hi), C.rMq(H, lo, hi),
            C.mincount(H, lo, hi), C.minselect(H, lo, hi, 1), C.root(H),
            C.depth(H, i), C.parent(H, i), C.is_tip(H, i), C.first_child(H, i),
            C.last_child(H, i), C.next_sibling(H, i), C.previous_sibling(H, i),
            C.preorder_rank(H, i), C.preorder_select(H, i % n),
            C.postorder_rank(H, i), C.postorder_select(H, i % (n + 1)),
            C.is_ancestor(H, lo, hi), C.count(H, i, True),
            C.level_ancestor(H, i, 1), C.lca(H, i, j), C.height(H, i),
        )
        assert len(exp) == COLUMNS
        assert out[i].tolist() == [int(v) for v in exp], (i, out[i], exp)
    print("OK")
""")


@numba_code
class GPUPrimitiveTests(TestCase):
    def test_simulator(self):
        # the same source compiles as device functions: checked on Numba's CUDA
        # simulator, which needs no GPU (a real device is not tested here)
        env = dict(os.environ, NUMBA_ENABLE_CUDASIM="1")
        res = subprocess.run([sys.executable, "-c", _SIMULATOR_SCRIPT], env=env,
                             capture_output=True, text=True)
        if res.stdout.strip() == "SKIP":
            self.skipTest("numba.cuda is not available.")
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertEqual(res.stdout.strip(), "OK")


if __name__ == "__main__":
    main()
