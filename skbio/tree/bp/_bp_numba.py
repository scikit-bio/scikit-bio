# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Derived from improved-octo-waddle (https://github.com/biocore/improved-octo-waddle)
# originally authored by Daniel McDonald, distributed under the Modified BSD License.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Numba compute engine for :class:`skbio.tree.BPTree`.

The navigation primitives of the Cython engine (``_bp_cy._BPKernel``), ported
to Numba functions over the tree's flat arrays, and the batch kernels built on
them.

Primitives
----------
Every navigation operation of ``BPTree`` has a primitive here, together with
the index operations they are built from (``rank``, ``select``, ``excess``,
``fwdsearch``, ``bwdsearch``, ``open``, ``enclose``). Each takes the tree's
arrays as a :class:`BPArrays` ``T`` and returns what the ``BPTree`` method of
the same name returns, with ``-1`` in place of ``None`` (``minselect``) and a
boolean for ``is_tip`` and ``is_ancestor``. Node attributes are not
primitives: a kernel reads ``lengths[i]``, ``edges[i]`` or
``edge_lookup[n]`` directly, and names (Python objects) stay on the host.

The primitives follow the Cython implementation, with two differences in
form but not in result: the range queries over the range min-max (rmM) tree
are iterative (bottom-up) rather than recursive, visiting a different but
equivalent set of nodes that covers the same blocks; and the methods that
recurse once to handle a closing parenthesis instead move to the matching
opening parenthesis first.

They are meant to be called from compiled code: a batch kernel, a per-sample
kernel, or a user's own ``@njit`` function. A single call from Python pays
Numba's argument dispatch (about a microsecond), far more than the operation
itself, which is why ``BPTree``'s per-node methods use the Cython engine.

Single source for CPU and GPU
-----------------------------
The primitives are written once, in :func:`define_primitives`, and compiled
by the decorator it is given.

Numba is an optional dependency: :func:`define_primitives` and
:func:`bp_arrays` are always defined, while :data:`CPU` and the batch kernels
are defined only if it is installed (``NUMBA_AVAILABLE``).
"""

from collections import namedtuple

import numpy as np

try:
    from numba import njit, prange

    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False


# The arrays and rmM-tree geometry of a tree, as passed to the kernels.
BPArrays = namedtuple(
    "BPArrays",
    [
        "B",  # parentheses (uint8)
        "e_index",  # excess at each position
        "k_index_0",  # position of the k-th closing parenthesis
        "k_index_1",  # position of the k-th opening parenthesis
        "m",  # rmM tree: minimum excess per node, heap order
        "M",  # rmM tree: maximum excess per node, heap order
        "r",  # rmM tree: rank (opening parentheses) before each node
        "b",  # block size
        "height",  # rmM tree height
        "n_internal",  # rmM tree internal node count
        "size",  # number of parentheses
    ],
)


def bp_arrays(tree, asarray=None):
    """Collect the arrays of a :class:`~skbio.tree.BPTree` for the kernels.

    Parameters
    ----------
    tree : skbio.tree.BPTree
        The tree.
    asarray : callable, optional
        Applied to each array, e.g. ``numba.cuda.to_device`` to place them on
        a GPU. By default the tree's own (host) arrays are used as they are.

    Returns
    -------
    BPArrays
        The arrays and rmM-tree geometry of the tree.
    """
    arrays = (
        tree._data,
        tree._e_index,
        tree._k_index_0,
        tree._k_index_1,
        tree._m,
        tree._M,
        tree._r,
    )
    if asarray is not None:
        arrays = tuple(asarray(a) for a in arrays)
    return BPArrays(
        *arrays,
        tree._b,
        tree._height,
        (1 << tree._height) - 1,
        tree._size,
    )


# The primitives, by name; see the module docstring.
Primitives = namedtuple(
    "Primitives",
    [
        # index operations
        "rank",
        "select",
        "excess",
        "fwdsearch",
        "bwdsearch",
        "open",
        "close",
        "enclose",
        "rmq",
        "rMq",
        "mincount",
        "minselect",
        # navigation
        "root",
        "depth",
        "parent",
        "is_tip",
        "first_child",
        "last_child",
        "next_sibling",
        "previous_sibling",
        "preorder_rank",
        "preorder_select",
        "postorder_rank",
        "postorder_select",
        "is_ancestor",
        "count",
        "level_ancestor",
        "level_next",
        "lca",
        "deepest_node",
        "height",
    ],
)

# Sentinels of the rmM-tree range queries (an empty range).
_SIZE_MAX = int(np.iinfo(np.intp).max)


def define_primitives(jit, jit_inline=None):
    """Compile the navigation primitives with a Numba decorator.

    Parameters
    ----------
    jit : callable
        The decorator compiling each primitive: ``numba.njit`` for the CPU, or
        ``gpu.jit(device=True)`` for device functions of a Numba GPU module.
    jit_inline : callable, optional
        The decorator for the small leaf primitives that the others call on
        every step (``rank``, ``select``, ``excess``, ``open``, ``close``,
        ``enclose``, ``depth``, ``parent``, ``is_tip``, ``is_ancestor``) and for
        the helpers on the rmM tree's node indices, e.g.
        ``numba.njit(inline="always")``. Defaults to ``jit``.

    Notes
    -----
    Inlining the leaf primitives removes a function call per step from the
    operations built on them. The searches (``fwdsearch``, ``bwdsearch``) and
    the range queries stay out of line: inlining them too copies them into
    every caller, which multiplies compile time (several-fold for ``lca``).

    Returns
    -------
    Primitives
        The compiled primitives. Each calls the others through this closure,
        so a set is compiled as a whole for one target.
    """
    if jit_inline is None:
        jit_inline = jit

    # -- the rmM tree: a complete binary tree in breadth-first (heap) order ----

    @jit_inline
    def bt_is_root(v):
        return v == 0

    @jit_inline
    def bt_is_left_child(v):
        return 0 if v == 0 else v % 2

    @jit_inline
    def bt_is_right_child(v):
        return 0 if v == 0 else 1 - (v % 2)

    @jit_inline
    def bt_parent(v):
        return 0 if v == 0 else (v - 1) // 2

    @jit_inline
    def bt_is_leaf(v, n_internal):
        return v >= n_internal

    @jit
    def tree_min(T, lo, hi):
        """Minimum excess over the rmM leaf blocks ``[lo, hi]`` (bottom-up).

        Only fully covered nodes are read. Every block in ``[lo, hi]`` is a
        real block, as the callers pass interior blocks only (see
        ``_BPKernel._rmq_tree_min``).
        """
        res = _SIZE_MAX
        lo += T.n_internal
        hi += T.n_internal
        while lo <= hi:
            if lo % 2 == 0:  # a right child: take it, continue to its right
                res = min(res, T.m[lo])
                lo += 1
            if hi % 2 == 1:  # a left child: take it, continue to its left
                res = min(res, T.m[hi])
                hi -= 1
            lo = (lo - 1) // 2
            hi = (hi - 1) // 2
        return res

    @jit
    def tree_max(T, lo, hi):
        """Maximum excess over the rmM leaf blocks ``[lo, hi]`` (bottom-up)."""
        res = -_SIZE_MAX
        lo += T.n_internal
        hi += T.n_internal
        while lo <= hi:
            if lo % 2 == 0:
                res = max(res, T.M[lo])
                lo += 1
            if hi % 2 == 1:
                res = max(res, T.M[hi])
                hi -= 1
            lo = (lo - 1) // 2
            hi = (hi - 1) // 2
        return res

    # -- index operations --------------------------------------------------

    @jit_inline
    def rank(T, t, i):
        """Number of ``t`` bits in ``B[0..i]``: preorder rank for ``t = 1``."""
        k = i // T.b
        upper = min(min((k + 1) * T.b, T.size), i + 1)
        r = 0
        for j in range(k * T.b, upper):
            r += T.B[j]
        # opening parentheses before the block
        r += T.r[T.n_internal + k]
        if t:
            return r
        return (i - r) + 1

    @jit_inline
    def select(T, t, k):
        """Position of the ``k``-th ``t`` bit (``k`` from 1)."""
        if t:
            return T.k_index_1[k]
        return T.k_index_0[k]

    @jit_inline
    def excess(T, i):
        """Opening minus closing parentheses in ``B[0..i]``."""
        return T.e_index[i]

    @jit
    def scan_block_forward(T, i, k, d):
        lower = max(max(k, 0) * T.b, i + 1)
        upper = min((k + 1) * T.b, T.size)
        for j in range(lower, upper):
            if T.e_index[j] == d:
                return j
        return -1

    @jit
    def scan_block_backward(T, i, k, d):
        lower = max(k, 0) * T.b - 1
        if lower >= 0:
            lower -= 1
        upper = min((k + 1) * T.b, T.size) - 1
        upper = min(i - 1, upper)
        if upper <= 0:
            return -1
        for j in range(upper, lower, -1):
            if T.e_index[j] == d:
                return j
        return -1

    @jit
    def fwdsearch(T, i, d):
        """Position after ``i`` with excess ``excess(i) + d``, or -1."""
        k = i // T.b
        d += T.e_index[i]
        node = T.n_internal + k
        result = -1
        if T.m[node] <= d <= T.M[node]:
            result = scan_block_forward(T, i, k, d)
        if result == -1:
            while not bt_is_root(node):
                if bt_is_left_child(node):
                    node += 1
                    if T.m[node] <= d <= T.M[node]:
                        break
                node = bt_parent(node)
            if bt_is_root(node):
                return -1
            while not bt_is_leaf(node, T.n_internal):
                node = 2 * node + 1
                if not (T.m[node] <= d <= T.M[node]):
                    node += 1
            k = node - T.n_internal
            result = scan_block_forward(T, i, k, d)
        return result

    @jit
    def bwdsearch(T, i, d):
        """Position before ``i`` with excess ``excess(i) + d``, or -1."""
        k = i // T.b
        d += T.e_index[i]
        result = scan_block_backward(T, i, k, d)
        node = T.n_internal + k
        if result == -1 and bt_is_right_child(node):
            node -= 1
            k = node - T.n_internal
            result = scan_block_backward(T, i, k, d)
            k = i // T.b
            node += 1
        if result == -1:
            while not bt_is_root(node):
                if bt_is_right_child(node):
                    node -= 1
                    if T.m[node] <= d <= T.M[node]:
                        break
                node = bt_parent(node)
            if bt_is_root(node):
                return -1
            while not bt_is_leaf(node, T.n_internal):
                node = 2 * node + 2
                if not (T.m[node] <= d <= T.M[node]):
                    node -= 1
            k = node - T.n_internal
            result = scan_block_backward(T, i, k, d)
        return result

    @jit_inline
    def open(T, i):
        if T.B[i] or i <= 0:
            return i
        return bwdsearch(T, i, 0) + 1

    @jit_inline
    def close(T, i):
        if not T.B[i]:
            return i
        return fwdsearch(T, i, -1)

    @jit_inline
    def enclose(T, i):
        if T.B[i]:
            return bwdsearch(T, i, -2) + 1
        return bwdsearch(T, i - 1, -2) + 1

    @jit
    def rmq(T, i, j):
        """Leftmost position of the minimum excess in ``[i, j]``."""
        if i >= j:
            return i
        b = T.b
        bi = i // b
        bj = j // b
        e_i = T.e_index[i]
        d_star = e_i
        if bi == bj:
            for p in range(i + 1, j + 1):
                if T.e_index[p] < d_star:
                    d_star = T.e_index[p]
        else:
            for p in range(i + 1, min((bi + 1) * b, T.size)):
                if T.e_index[p] < d_star:
                    d_star = T.e_index[p]
            if bi + 1 <= bj - 1:
                tree_v = tree_min(T, bi + 1, bj - 1)
                if tree_v < d_star:
                    d_star = tree_v
            for p in range(bj * b, j + 1):
                if T.e_index[p] < d_star:
                    d_star = T.e_index[p]
        if e_i == d_star:
            return i
        return fwdsearch(T, i, d_star - e_i)

    @jit
    def rMq(T, i, j):
        """Leftmost position of the maximum excess in ``[i, j]``."""
        if i >= j:
            return i
        b = T.b
        bi = i // b
        bj = j // b
        e_i = T.e_index[i]
        m_star = e_i
        if bi == bj:
            for p in range(i + 1, j + 1):
                if T.e_index[p] > m_star:
                    m_star = T.e_index[p]
        else:
            for p in range(i + 1, min((bi + 1) * b, T.size)):
                if T.e_index[p] > m_star:
                    m_star = T.e_index[p]
            if bi + 1 <= bj - 1:
                tree_v = tree_max(T, bi + 1, bj - 1)
                if tree_v > m_star:
                    m_star = tree_v
            for p in range(bj * b, j + 1):
                if T.e_index[p] > m_star:
                    m_star = T.e_index[p]
        if e_i == m_star:
            return i
        return fwdsearch(T, i, m_star - e_i)

    @jit
    def mincount(T, i, j):
        """Number of occurrences of the minimum excess in ``[i, j]``."""
        lo = T.e_index[i]
        c = 0
        for p in range(i, j + 1):
            e = T.e_index[p]
            if e < lo:
                lo = e
                c = 1
            elif e == lo:
                c += 1
        return c

    @jit
    def minselect(T, i, j, q):
        """Position of the ``q``-th (from 1) minimum excess in ``[i, j]``, or -1."""
        lo = T.e_index[i]
        for p in range(i + 1, j + 1):
            if T.e_index[p] < lo:
                lo = T.e_index[p]
        if q < 1:
            return -1
        for p in range(i, j + 1):
            if T.e_index[p] == lo:
                q -= 1
                if q == 0:
                    return p
        return -1

    # -- navigation ----------------------------------------------------------

    @jit
    def root(T):
        return 0

    @jit_inline
    def depth(T, i):
        return T.e_index[i]

    @jit_inline
    def parent(T, i):
        """Parent of node ``i``, or -1 for the root."""
        if i == 0 or i == T.size - 1:
            return -1
        return enclose(T, i)

    @jit_inline
    def is_tip(T, i):
        # a closing parenthesis short-circuits, so ``i + 1`` stays in range
        return T.B[i] == 1 and T.B[i + 1] == 0

    @jit
    def first_child(T, i):
        """First child of node ``i``, or 0 for a tip."""
        if not T.B[i]:
            i = open(T, i)
        if is_tip(T, i):
            return 0
        return i + 1

    @jit
    def last_child(T, i):
        """Last child of node ``i``, or 0 for a tip."""
        if not T.B[i]:
            i = open(T, i)
        if is_tip(T, i):
            return 0
        return open(T, close(T, i) - 1)

    @jit
    def next_sibling(T, i):
        """Next sibling of node ``i``, or 0 if none."""
        if not T.B[i]:
            i = open(T, i)
        pos = close(T, i) + 1
        if pos >= T.size or not T.B[pos]:
            return 0
        return pos

    @jit
    def previous_sibling(T, i):
        """Previous sibling of node ``i``, or 0 if none."""
        if not T.B[i]:
            i = open(T, i)
        if T.B[max(0, i - 1)]:
            return 0
        pos = open(T, i - 1)
        if pos < 0 or not T.B[pos]:
            return 0
        return pos

    @jit
    def preorder_rank(T, i):
        if not T.B[i]:
            i = open(T, i)
        return rank(T, 1, i)

    @jit
    def preorder_select(T, k):
        return select(T, 1, k)

    @jit
    def postorder_rank(T, i):
        if T.B[i]:
            return rank(T, 0, close(T, i))
        return rank(T, 0, i)

    @jit
    def postorder_select(T, k):
        return open(T, select(T, 0, k))

    @jit_inline
    def is_ancestor(T, i, j):
        """Whether node ``i`` is an ancestor of node ``j``."""
        if i == j:
            return False
        if not T.B[i]:
            i = open(T, i)
        return i <= j < close(T, i)

    @jit
    def count(T, i, tips):
        """Number of nodes (or of tips) in the subtree of node ``i``."""
        if not T.B[i]:
            i = open(T, i)
        last = close(T, i)
        if not tips:
            return (last - i + 1) // 2
        c = 0
        j = i
        while j < last:
            if T.B[j] and not T.B[j + 1]:
                c += 1
                j += 1
            j += 1
        return c

    @jit
    def level_ancestor(T, i, d):
        """Ancestor ``d`` levels above node ``i``, or -1 for ``d <= 0``."""
        if d <= 0:
            return -1
        if not T.B[i]:
            i = open(T, i)
        return bwdsearch(T, i, -d - 1) + 1

    @jit
    def level_next(T, i):
        return fwdsearch(T, close(T, i), 1)

    @jit
    def lca(T, i, j):
        """Lowest common ancestor of nodes ``i <= j``."""
        if i == j:
            return open(T, i)
        if is_ancestor(T, i, j):
            return i
        elif is_ancestor(T, j, i):
            return j
        return parent(T, rmq(T, i, j) + 1)

    @jit
    def deepest_node(T, i):
        return rMq(T, open(T, i), close(T, i))

    @jit
    def height(T, i):
        """Height of node ``i``, in edges."""
        return T.e_index[deepest_node(T, i)] - T.e_index[open(T, i)]

    scope = locals()
    return Primitives(*(scope[name] for name in Primitives._fields))


# Device-function primitives, by Numba GPU module (built on first use).
_gpu_primitives = {}


def gpu_primitives(gpu):
    """The primitives compiled as device functions of a Numba GPU module.

    Parameters
    ----------
    gpu : module
        ``numba.cuda`` or ``numba.hip`` (e.g. from
        ``skbio.stats.distance._gpu._numba_gpu_module_for``).

    Returns
    -------
    Primitives
        Device functions, callable from the ``gpu.jit`` kernels of that module.
        Each is compiled when a kernel calling it is.
    """
    name = gpu.__name__
    if name not in _gpu_primitives:
        _gpu_primitives[name] = define_primitives(gpu.jit(device=True))
    return _gpu_primitives[name]


if NUMBA_AVAILABLE:
    #: The primitives compiled for the CPU.
    CPU = define_primitives(njit, njit(inline="always"))

    # used by the batch kernels below
    _close = CPU.close
    _parent = CPU.parent
    _lca = CPU.lca
    _level_ancestor = CPU.level_ancestor

    # -- batch kernels (engine="numba") ---------------------------------------

    @njit(parallel=True)
    def close_batch(T, idx):
        """Kernel of :meth:`skbio.tree.BPTree.close_batch`."""
        out = np.empty(idx.shape[0], dtype=np.intp)
        for t in prange(idx.shape[0]):
            out[t] = _close(T, idx[t])
        return out

    @njit(parallel=True)
    def parent_batch(T, idx):
        """Kernel of :meth:`skbio.tree.BPTree.parent_batch`."""
        out = np.empty(idx.shape[0], dtype=np.intp)
        for t in prange(idx.shape[0]):
            out[t] = _parent(T, idx[t])
        return out

    @njit(parallel=True)
    def lca_batch(T, i, j):
        """Kernel of :meth:`skbio.tree.BPTree.lca_batch`."""
        out = np.empty(i.shape[0], dtype=np.intp)
        for t in prange(i.shape[0]):
            out[t] = _lca(T, min(i[t], j[t]), max(i[t], j[t]))
        return out

    @njit(parallel=True)
    def level_ancestor_batch(T, idx, d):
        """Kernel of :meth:`skbio.tree.BPTree.level_ancestor_batch`."""
        out = np.empty(idx.shape[0], dtype=np.intp)
        for t in prange(idx.shape[0]):
            out[t] = _level_ancestor(T, idx[t], d[t])
        return out

    @njit
    def root_distances(T, lengths):
        """Sum of branch lengths from the root to each node, per position."""
        out = np.zeros(T.size, dtype=np.float64)
        stack = np.zeros(T.size // 2 + 1, dtype=np.float64)
        top = 0
        for i in range(1, T.size):
            if T.B[i]:
                top += 1
                stack[top] = stack[top - 1] + lengths[i]
                out[i] = stack[top]
            else:
                top -= 1
        return out

    @njit
    def _tip_distance_row(tips, slot, parent, end, dist, out, row, n):
        r = slot[row]
        out[r, r] = 0.0
        start = row + 1
        da = dist[tips[row]]
        v = parent[tips[row]]
        while start < n:
            stop = end[v]
            dv = dist[v]
            for col in range(start, stop):
                d = (da - dv) + (dist[tips[col]] - dv)
                c = slot[col]
                out[r, c] = d
                out[c, r] = d
            start = stop
            v = parent[v]

    @njit(parallel=True)
    def tip_distances(tips, slot, parent, end, dist):
        """Pairwise path distances between tips, as a square matrix.

        See ``_bp_cy.tip_distances``; rows ``h`` and ``n - 2 - h`` share an
        iteration to balance the triangular workload.
        """
        n = tips.shape[0]
        n_rows = n - 1 if n > 0 else 0
        out = np.empty((n, n), dtype=np.float64)
        for h in prange((n_rows + 1) // 2):
            _tip_distance_row(tips, slot, parent, end, dist, out, h, n)
            if n_rows - 1 - h != h:
                _tip_distance_row(tips, slot, parent, end, dist, out, n_rows - 1 - h, n)
        if n:
            out[slot[n - 1], slot[n - 1]] = 0.0
        return out
