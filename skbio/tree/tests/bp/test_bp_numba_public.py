# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""The public Numba API of BPTree: skbio.tree.bp.numba, and the tree's side."""

import os
import subprocess
import sys
import textwrap
from unittest import TestCase, main, mock

import array_api_compat as aac
import numpy as np
import numpy.testing as npt

from skbio.tree import BPTree, TreeNode
from skbio.tree._exception import MissingNodeError
from skbio.tree.bp import _bp_numba, _gpu
from skbio.tree.bp import numba as bpn
from skbio.tree.tests.bp.test_bp import _caterpillar, _named_tree
from skbio.util._array import _get_backend_name
from skbio.util._testing import ArrayAPITestMixin, array_backends, numba_code


# The kernels of the module documentation, written as a user would: the
# navigation functions are a module global.
if _bp_numba.NUMBA_AVAILABLE:
    from numba import njit, prange

    nav = bpn.cpu()

    @njit(parallel=True)
    def sample_lcas(T, tips, offsets):
        out = np.empty(offsets.size - 1, np.intp)
        for s in prange(out.size):
            a = tips[offsets[s]]
            for k in range(offsets[s] + 1, offsets[s + 1]):
                a = nav.lca(T, a, tips[k])
            out[s] = a
        return out

    @njit(parallel=True)
    def root_distances(T, lengths, nodes):
        out = np.zeros(nodes.size)
        root = nav.root(T)
        for k in prange(nodes.size):
            v = nodes[k]
            while v != root:
                out[k] += lengths[v]
                v = nav.parent(T, v)
        return out


class KernelInputTests(TestCase):
    """What a tree provides to kernels: arrays, attributes, positions."""

    def setUp(self):
        # "x" names both a tip and an internal node; "c" only an internal node
        self.newick = "((a:1,x:2)c:3,(d:4,e:5)x:6)root:0.5;"
        self.bp = BPTree.read([self.newick])

    def test_lengths_and_edges(self):
        bp = self.bp
        npt.assert_array_equal(bp.lengths, bp._lengths)
        npt.assert_array_equal(bp.edges, bp._edges)
        for arr in (bp.lengths, bp.edges):
            self.assertFalse(arr.flags.writeable)
            with self.assertRaises(ValueError):
                arr[0] = 1
        # views of the tree's arrays, which stay writable for the tree
        self.assertTrue(np.shares_memory(bp.lengths, bp._lengths))
        self.assertTrue(bp._lengths.flags.writeable)
        # replaced by set_lengths
        new = np.arange(bp.data.size, dtype=np.float64)
        bp.set_lengths(new)
        npt.assert_array_equal(bp.lengths, new)
        self.assertEqual(bp.length(1), 1.0)

    def test_tip_positions(self):
        bp = self.bp
        obs = bp.tip_positions()
        self.assertEqual(obs.dtype, np.intp)
        exp = [i for i in range(bp.data.size - 1) if bp.is_tip(i)]
        npt.assert_array_equal(obs, exp)
        self.assertEqual([bp.name(int(i)) for i in obs], ['a', 'x', 'd', 'e'])
        # a single node is a tip
        npt.assert_array_equal(
            BPTree(np.array([1, 0], dtype=np.uint8)).tip_positions(), [0])

    def test_positions(self):
        bp = self.bp
        # read directly: TreeNode.from_bptree drops the root's length
        tn = TreeNode.read([self.newick])
        names = ['e', 'c', 'root', 'x', 'a']
        obs = bp.positions(names)
        self.assertEqual(obs.dtype, np.intp)
        # resolved as TreeNode.find: the tip "x" over the internal node "x"
        self.assertEqual([bp.name(int(i)) for i in obs], names)
        self.assertTrue(bp.is_tip(int(obs[3])))
        for name, pos in zip(names, obs):
            self.assertEqual(bp.length(int(pos)), tn.find(name).length)
        # one name gives an int
        self.assertEqual(bp.positions('c'), 1)
        self.assertIsInstance(bp.positions('c'), int)
        self.assertEqual(bp.positions([]).shape, (0,))
        self.assertEqual(bp.positions([]).dtype, np.intp)
        with self.assertRaises(MissingNodeError):
            bp.positions(['a', 'nope'])
        with self.assertRaises(MissingNodeError):
            bp.positions([None])

    def test_numba_arrays(self):
        bp = self.bp
        T = bp.numba_arrays()
        self.assertEqual(T.size, bp.data.size)
        # read-only views of the tree's own arrays: no copy
        for field, own in (('B', bp._data), ('e_index', bp._e_index),
                           ('k_index_0', bp._k_index_0), ('m', bp._m)):
            arr = getattr(T, field)
            self.assertTrue(np.shares_memory(arr, own), field)
            self.assertFalse(arr.flags.writeable, field)
        # the tree's topology array itself is left as it was
        B = np.array([1, 1, 0, 0], dtype=np.uint8)
        BPTree(B).numba_arrays()
        self.assertTrue(B.flags.writeable)


@numba_code
class NavigationTests(TestCase):
    def setUp(self):
        rng = np.random.default_rng(4)
        self.trees = [BPTree(np.array([1, 0], dtype=np.uint8)),
                      BPTree(_caterpillar(40))]
        self.trees += [_named_tree(n, rng) for n in (2, 9, 200)]
        self.rng = rng

    def test_fields_mirror_methods(self):
        nav = bpn.cpu()
        self.assertIs(bpn.cpu(), nav)
        self.assertIsInstance(nav, bpn.Navigation)
        for name in bpn.Navigation._fields:
            # a BPTree method of the same name, and the compiled primitive
            self.assertTrue(callable(getattr(BPTree, name)), name)
            self.assertIs(getattr(nav, name), getattr(_bp_numba.CPU, name))
        # the index operations stay private
        for name in ('rank', 'select', 'excess', 'open', 'enclose', 'fwdsearch',
                     'bwdsearch'):
            self.assertNotIn(name, bpn.Navigation._fields)

    def test_matches_methods(self):
        # the same results, and sentinels, as the methods; see the conventions
        # in the module documentation
        nav = bpn.cpu()
        for bp in self.trees:
            T = bp.numba_arrays()
            n = bp.data.size
            for i in range(n - 1):
                for name in ('close', 'depth', 'parent', 'first_child',
                             'last_child', 'next_sibling', 'previous_sibling',
                             'preorder_rank', 'postorder_rank', 'deepest_node',
                             'height', 'level_next'):
                    self.assertEqual(getattr(nav, name)(T, i),
                                     getattr(bp, name)(i), (name, i))
                self.assertIs(bool(nav.is_tip(T, i)), bool(bp.is_tip(i)))
                self.assertEqual(nav.count(T, i, True), bp.count(i, tips=True))
                self.assertEqual(nav.level_ancestor(T, i, 2),
                                 bp.level_ancestor(i, 2))
            self.assertEqual(nav.root(T), bp.root())
            for _ in range(100):
                i, j = sorted(self.rng.integers(0, n, 2).tolist())
                self.assertEqual(nav.lca(T, i, j), bp.lca(i, j))
                self.assertIs(bool(nav.is_ancestor(T, i, j)),
                              bool(bp.is_ancestor(i, j)))
                self.assertEqual(nav.rmq(T, i, j), bp.rmq(i, j))
                self.assertEqual(nav.rMq(T, i, j), bp.rMq(i, j))
                self.assertEqual(nav.mincount(T, i, j), bp.mincount(i, j))
                for q in (-1, 0, 1, 2, 3):
                    exp = bp.minselect(i, j, q)
                    self.assertEqual(nav.minselect(T, i, j, q),
                                     -1 if exp is None else exp)

    def test_documented_kernels(self):
        # the lca of each sample's tips, against TreeNode.lca of the set
        bp = self.trees[-1]
        tn = TreeNode.from_bptree(bp)
        tips = bp.tip_positions()
        sizes = self.rng.integers(1, 6, 30)
        offsets = np.concatenate([[0], np.cumsum(sizes)])
        members = self.rng.choice(tips, offsets[-1])
        obs = sample_lcas(bp.numba_arrays(), members, offsets)
        for s in range(sizes.size):
            names = [bp.name(int(p)) for p in members[offsets[s]:offsets[s + 1]]]
            self.assertEqual(bp.name(int(obs[s])), tn.lca(names).name, s)

        # the distance of every node from the root, against TreeNode
        nodes = np.flatnonzero(bp.data)
        obs = root_distances(bp.numba_arrays(), bp.lengths, nodes)
        exp = [n.distance(tn) for n in tn.preorder(include_self=True)]
        npt.assert_allclose(obs, exp, rtol=1e-12, atol=1e-12)

    def test_without_numba(self):
        with mock.patch.object(_bp_numba, 'NUMBA_AVAILABLE', False):
            with self.assertRaisesRegex(ImportError, 'require Numba'):
                bpn.cpu()
            with self.assertRaisesRegex(ImportError, 'require Numba'):
                bpn.gpu(object())


# The GPU example of the module documentation, in a fresh interpreter: on the
# CUDA simulator (NUMBA_ENABLE_CUDASIM=1), or on a real device.
_GPU_SCRIPT = textwrap.dedent("""
    import os
    import sys
    import numpy as np
    if os.environ.get("NUMBA_ENABLE_CUDASIM") == "1":
        try:
            from numba import cuda
        except ImportError:
            print("SKIP"); sys.exit(0)
    else:
        # the Numba GPU module that BPTree's dispatch chose (see test_device)
        import importlib
        cuda = importlib.import_module(os.environ["SKBIO_TEST_GPU_MODULE"])

    from skbio.tree.bp import numba as bpn
    from skbio.tree.tests.bp.test_bp import _named_tree

    gnav = bpn.gpu(cuda)
    assert bpn.gpu(cuda) is gnav

    @cuda.jit
    def sample_lcas_gpu(T, tips, offsets, out):
        s = cuda.grid(1)
        if s < out.shape[0]:
            a = tips[offsets[s]]
            for k in range(offsets[s] + 1, offsets[s + 1]):
                a = gnav.lca(T, a, tips[k])
            out[s] = a

    rng = np.random.default_rng(8)
    bp = _named_tree(60, rng)
    T = bp.numba_arrays(gpu=cuda)
    assert bp.numba_arrays(gpu=cuda) is T  # uploaded once, cached on the tree
    sizes = rng.integers(1, 5, 20)
    offsets = np.concatenate([[0], np.cumsum(sizes)])
    members = rng.choice(bp.tip_positions(), offsets[-1])
    out = cuda.device_array(sizes.size, dtype=np.intp)
    sample_lcas_gpu[1, 32](T, cuda.to_device(members), cuda.to_device(offsets), out)

    exp = []
    for s in range(sizes.size):
        a = int(members[offsets[s]])
        for b in members[offsets[s] + 1:offsets[s + 1]]:
            a = bp.lca(a, int(b))
        exp.append(a)
    assert out.copy_to_host().tolist() == exp, (out.copy_to_host(), exp)

    # every function, compiled as a device function, against the CPU's
    nav = bpn.cpu()
    one = [("close",), ("depth",), ("parent",), ("is_tip",), ("first_child",),
           ("last_child",), ("next_sibling",), ("previous_sibling",),
           ("preorder_rank",), ("postorder_rank",), ("deepest_node",),
           ("height",), ("level_next",), ("count", True), ("count", False),
           ("level_ancestor", 0), ("level_ancestor", 1),
           ("level_ancestor", 3)]
    two = [("lca",), ("is_ancestor",), ("rmq",), ("rMq",), ("mincount",),
           ("minselect", -1), ("minselect", 1), ("minselect", 2)]
    covered = {c[0] for c in one + two}
    covered |= {"root", "preorder_select", "postorder_select"}
    assert covered == set(bpn.Navigation._fields), covered

    @cuda.jit
    def nodes(T, idx, ks, out):
        t = cuda.grid(1)
        if t < idx.shape[0]:
            i = idx[t]
            out[t, 0] = gnav.close(T, i)
            out[t, 1] = gnav.depth(T, i)
            out[t, 2] = gnav.parent(T, i)
            out[t, 3] = gnav.is_tip(T, i)
            out[t, 4] = gnav.first_child(T, i)
            out[t, 5] = gnav.last_child(T, i)
            out[t, 6] = gnav.next_sibling(T, i)
            out[t, 7] = gnav.previous_sibling(T, i)
            out[t, 8] = gnav.preorder_rank(T, i)
            out[t, 9] = gnav.postorder_rank(T, i)
            out[t, 10] = gnav.deepest_node(T, i)
            out[t, 11] = gnav.height(T, i)
            out[t, 12] = gnav.level_next(T, i)
            out[t, 13] = gnav.count(T, i, True)
            out[t, 14] = gnav.count(T, i, False)
            out[t, 15] = gnav.level_ancestor(T, i, 0)
            out[t, 16] = gnav.level_ancestor(T, i, 1)
            out[t, 17] = gnav.level_ancestor(T, i, 3)
            out[t, 18] = gnav.root(T)
            out[t, 19] = gnav.preorder_select(T, ks[t])
            out[t, 20] = gnav.postorder_select(T, ks[t])

    @cuda.jit
    def pairs(T, ii, jj, out):
        t = cuda.grid(1)
        if t < ii.shape[0]:
            i = ii[t]
            j = jj[t]
            out[t, 0] = gnav.lca(T, i, j)
            out[t, 1] = gnav.is_ancestor(T, i, j)
            out[t, 2] = gnav.rmq(T, i, j)
            out[t, 3] = gnav.rMq(T, i, j)
            out[t, 4] = gnav.mincount(T, i, j)
            out[t, 5] = gnav.minselect(T, i, j, -1)
            out[t, 6] = gnav.minselect(T, i, j, 1)
            out[t, 7] = gnav.minselect(T, i, j, 2)

    def call(T, entry, *args):
        return int(getattr(nav, entry[0])(T, *args, *entry[1:]))

    Th = bp.numba_arrays()
    n = Th.size
    idx = np.arange(n - 1)
    ks = idx % (n // 2 + 2)  # every rank, and 0 and n + 1, which are not
    out = cuda.device_array((idx.size, 21), dtype=np.intp)
    nodes[(idx.size + 127) // 128, 128](T, cuda.to_device(idx),
                                        cuda.to_device(ks), out)
    obs = out.copy_to_host()
    for t, i in enumerate(idx.tolist()):
        exp = [call(Th, e, i) for e in one]
        exp += [int(nav.root(Th)), int(nav.preorder_select(Th, int(ks[t]))),
                int(nav.postorder_select(Th, int(ks[t])))]
        assert obs[t].tolist() == exp, (i, obs[t].tolist(), exp)

    ij = np.sort(rng.integers(0, n, (300, 2)), axis=1)
    out = cuda.device_array((ij.shape[0], len(two)), dtype=np.intp)
    pairs[(ij.shape[0] + 127) // 128, 128](
        T, cuda.to_device(ij[:, 0].copy()), cuda.to_device(ij[:, 1].copy()), out)
    obs = out.copy_to_host()
    for t, (i, j) in enumerate(ij.tolist()):
        exp = [call(Th, e, i, j) for e in two]
        assert obs[t].tolist() == exp, ((i, j), obs[t].tolist(), exp)
    print("OK")
""")


def _run_gpu_script(test, env):
    res = subprocess.run([sys.executable, "-c", _GPU_SCRIPT], env=env,
                         capture_output=True, text=True)
    if res.stdout.strip() == "SKIP":
        test.skipTest("The numba.cuda simulator is not available.")
    test.assertEqual(res.returncode, 0, res.stderr)
    test.assertEqual(res.stdout.strip(), "OK")


@numba_code
class GPUNavigationTests(TestCase, ArrayAPITestMixin):
    def test_simulator(self):
        _run_gpu_script(self, dict(os.environ, NUMBA_ENABLE_CUDASIM="1"))

    @array_backends("jax", "torch", "cupy")
    def test_device(self, xp, device):
        # on a real GPU (the GPU workflow), through the Numba GPU module that
        # BPTree's dispatch chooses for the device. JAX has no Numba GPU path
        # and is skipped; it is listed because the harness errors on a GPU lane
        # that runs no backend.
        if device == "cpu" or _get_backend_name(xp) == "jax":
            self.skipTest("needs a device-resident CuPy or PyTorch array")
        gpu = _gpu._numba_gpu_module_for(self.make_array(xp, device, np.zeros(2)))
        self.assertIsNotNone(gpu, "No Numba GPU module for %s on %s: is "
                             "numba-cuda-mlir (or numba-cuda) installed?"
                             % (_get_backend_name(xp), device))
        env = {k: v for k, v in os.environ.items() if k != "NUMBA_ENABLE_CUDASIM"}
        env["SKBIO_TEST_GPU_MODULE"] = gpu.__name__
        _run_gpu_script(self, env)

    @array_backends("jax", "torch", "cupy")
    def test_device_tree(self, xp, device):
        # a tree whose data is on a GPU, created from a device array or moved
        # there with to_device: what a kernel takes from it
        if device == "cpu" or _get_backend_name(xp) == "jax":
            self.skipTest("needs a device-resident CuPy or PyTorch array")
        host = _named_tree(200, np.random.default_rng(9))
        B = self.make_array(xp, device, host.data, dtype=xp.uint8)
        gpu = _gpu._numba_gpu_module_for(B)
        self.assertIsNotNone(gpu, "No Numba GPU module for %s on %s"
                             % (_get_backend_name(xp), device))
        names = [host.name(int(i)) for i in host.tip_positions()]
        tips = host.tip_positions()
        exp_lca = host.lca_batch(tips, tips[::-1])
        q1, q2 = (self.make_array(xp, device, q.copy(), dtype=xp.int64)
                  for q in (tips, tips[::-1]))
        for bp in (BPTree(B, lengths=host._lengths, names=host._names),
                   host.to_device(xp, aac.device(B))):
            # positions and node attributes are NumPy arrays on the host
            for obs, exp in ((bp.tip_positions(), host.tip_positions()),
                             (bp.positions(names), host.positions(names)),
                             (bp.lengths, host.lengths), (bp.edges, host.edges)):
                self.assertIsInstance(obs, np.ndarray)
                npt.assert_array_equal(obs, exp)
            self.assertFalse(bp.lengths.flags.writeable)
            # T for the CPU: views of the tree's host copy
            T = bp.numba_arrays()
            self.assertTrue(np.shares_memory(T.B, bp._data))
            # T for the GPU: uploaded once, and the batch operations reuse it
            Tg = bp.numba_arrays(gpu=gpu)
            self.assertIs(bp.numba_arrays(gpu=gpu), Tg)
            npt.assert_array_equal(Tg.B.copy_to_host(), host.data)
            npt.assert_array_equal(Tg.m.copy_to_host(), host._m)
            obs = bp.lca_batch(q1, q2, engine="numba")
            self.assertIs(bp._gpu_arrays[gpu.__name__], Tg)
            self.assertNotIn(_get_backend_name(xp), _gpu._unavailable)
            npt.assert_array_equal(np.asarray(obs.tolist()), exp_lca)


if __name__ == "__main__":
    main()
