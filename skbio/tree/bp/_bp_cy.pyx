# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
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

"""Cython compute engine for :class:`skbio.tree.BPTree`.

``BPTree`` is a Python class that owns the tree's arrays. This module holds the
compiled kernels that operate on them:

- :class:`_BPKernel` is a handle holding typed, read-only views of those
  arrays. It implements the per-node navigation primitives. ``BPTree`` binds
  the handle's methods onto each instance, so a navigation call goes straight
  to compiled code without passing through a Python frame.
- The module-level functions are whole-tree kernels (conversion, shearing,
  collapsing), each a single compiled pass over a handle.

"""

### NOTE: some doctext strings are copied and pasted from manuscript
### http://www.dcc.uchile.cl/~gnavarro/ps/tcs16.2.pdf

from libc.math cimport pow

import numpy as np
cimport numpy as cnp
cimport cython
from cython.parallel cimport prange

from ._bp_binary_tree cimport *

cnp.import_array()

cdef extern from "Python.h":
    # largest value a Py_ssize_t can hold; used as a nogil "+infinity" sentinel
    cdef Py_ssize_t PY_SSIZE_T_MAX

# The compiler defines _OPENMP when it compiles with OpenMP, which setup.py
# enables where the compiler supports it (gcc, MSVC, Intel; not Apple clang).
cdef extern from *:
    """
    #ifdef _OPENMP
    #define SKBIO_BP_OPENMP 1
    #else
    #define SKBIO_BP_OPENMP 0
    #endif
    """
    int SKBIO_BP_OPENMP

#: Whether this module was compiled with OpenMP, so that its batch kernels
#: (``prange``) run in parallel rather than serially.
OPENMP = bool(SKBIO_BP_OPENMP)


DOUBLE = np.float64
SIZE = np.intp
BOOL = np.uint8
INT32 = np.int32


cdef inline Py_ssize_t min(Py_ssize_t a, Py_ssize_t b) noexcept nogil:
    if a > b:
        return b
    else:
        return a


cdef inline Py_ssize_t max(Py_ssize_t a, Py_ssize_t b) noexcept nogil:
    if a > b:
        return a
    else:
        return b


@cython.final
cdef class _BPKernel:
    """Compiled navigation over the arrays of a :class:`~skbio.tree.BPTree`.

    Not part of the public API. See ``BPTree`` for the documentation of each
    navigation method.
    """

    def __cinit__(self, const BOOL_t[::1] B,
                  const Py_ssize_t[::1] e_index,
                  const Py_ssize_t[::1] k_index_0,
                  const Py_ssize_t[::1] k_index_1,
                  const Py_ssize_t[::1] m,
                  const Py_ssize_t[::1] M,
                  const Py_ssize_t[::1] r,
                  Py_ssize_t b, Py_ssize_t height,
                  cnp.ndarray names, cnp.ndarray lengths, cnp.ndarray edges,
                  object edge_lookup):
        self._B = B
        self._b_ptr = &B[0]
        self.size = B.shape[0]
        self._e_index = e_index
        self._k_index_0 = k_index_0
        self._k_index_1 = k_index_1
        self._m = m
        self._M = M
        self._r = r
        self._b = b
        self._height = height
        self._n_internal = (<Py_ssize_t>1 << height) - 1
        self._names = names
        self._lengths = lengths
        self._edges = edges
        self._edge_lookup = edge_lookup

    cpdef inline unicode name(self, Py_ssize_t i):
        """Name of node ``i``."""
        return self._names[i]

    cpdef inline DOUBLE_t length(self, Py_ssize_t i):
        """Branch length of node ``i``."""
        return self._lengths[i]

    cpdef inline INT32_t edge(self, Py_ssize_t i):
        """Edge number of node ``i``."""
        return self._edges[i]

    cpdef Py_ssize_t edge_from_number(self, INT32_t n):
        """Node index of edge number ``n``."""
        return self._edge_lookup[n]

    cdef inline Py_ssize_t rank(self, Py_ssize_t t, Py_ssize_t i) noexcept nogil:
        """Determine the rank order of the ith bit t

        Rank is the order of the ith bit observed, from left to right. For
        t=1, this is a preorder traversal of the tree.

        Parameters
        ----------
        t : Py_ssize_t
            The bit value, either 0 or 1 where 0 is a closing parenthesis and
            1 is an opening.
        i : Py_ssize_t
            The position to evaluate

        Returns
        -------
        Py_ssize_t
            The rank order of the position.
        """
        cdef Py_ssize_t k
        cdef Py_ssize_t r = 0
        cdef Py_ssize_t lower_bound
        cdef Py_ssize_t upper_bound
        cdef Py_ssize_t j
        cdef Py_ssize_t node

        k = i // self._b

        lower_bound = k * self._b

        # upper_bound is block boundary or end of tree
        upper_bound = min((k + 1) * self._b, self.size)
        upper_bound = min(upper_bound, i + 1)

        # collect rank from within the block
        for j in range(lower_bound, upper_bound):
            r += self._b_ptr[j]

        # collect the rank at the left end of the block
        node = bt_node_from_left(k, self._height)
        r += self._r[node]

        if t:
            return r
        else:
            return (i - r) + 1

    cdef inline Py_ssize_t select(self, Py_ssize_t t, Py_ssize_t k) noexcept nogil:
        """The position in B of the kth occurrence of the bit t."""
        if t:
            return self._k_index_1[k]
        else:
            return self._k_index_0[k]

    cdef Py_ssize_t excess(self, Py_ssize_t i) noexcept nogil:
        """the number of opening minus closing parentheses in B[1, i]"""
        # same as: self.rank(1, i) - self.rank(0, i)
        return self._e_index[i]

    cpdef inline Py_ssize_t close(self, Py_ssize_t i) noexcept nogil:
        """The position of the closing parenthesis that matches B[i]"""
        if not self._b_ptr[i]:
            # identity: the close of a closed parenthesis is itself
            return i

        return self.fwdsearch(i, -1)

    cdef inline Py_ssize_t open(self, Py_ssize_t i) noexcept nogil:
        """The position of the opening parenthesis that matches B[i]"""
        if self._b_ptr[i] or i <= 0:
            # identity: the open of an open parenthesis is itself
            # the open of 0 is open. A negative index cannot be open, so just return
            return i

        return self.bwdsearch(i, 0) + 1

    cdef inline Py_ssize_t enclose(self, Py_ssize_t i) noexcept nogil:
        """The opening parenthesis of the smallest matching pair that contains position i"""
        if self._b_ptr[i]:
            return self.bwdsearch(i, -2) + 1
        else:
            return self.bwdsearch(i - 1, -2) + 1

    cpdef Py_ssize_t rmq(self, Py_ssize_t i, Py_ssize_t j) noexcept nogil:
        """The leftmost minimum excess in i -> j.

        The minimum excess over [i, j] is found in O(log n): the two partial
        end-blocks are scanned directly and the full interior blocks are
        covered by a range query over the rmM tree. The leftmost position that
        attains it is then recovered with a single ``fwdsearch``.
        """
        cdef:
            Py_ssize_t d_star, p, blk_end, tree_v
            Py_ssize_t bi, bj, e_i
            Py_ssize_t b

        if i >= j:
            return i

        b = self._b
        bi = i // b
        bj = j // b

        # smallest absolute excess over [i, j]
        e_i = self._e_index[i]
        d_star = e_i
        if bi == bj:
            for p in range(i + 1, j + 1):
                if self._e_index[p] < d_star:
                    d_star = self._e_index[p]
        else:
            # remainder of the block containing i (inclusive of i)
            blk_end = min((bi + 1) * b, self.size)
            for p in range(i + 1, blk_end):
                if self._e_index[p] < d_star:
                    d_star = self._e_index[p]
            # full interior blocks (bi, bj) via the rmM tree
            if bi + 1 <= bj - 1:
                tree_v = self._rmq_tree_min(0, 0, self._n_internal,
                                            bi + 1, bj - 1)
                if tree_v < d_star:
                    d_star = tree_v
            # leading part of the block containing j (inclusive of j)
            for p in range(bj * b, j + 1):
                if self._e_index[p] < d_star:
                    d_star = self._e_index[p]

        # leftmost position in [i, j] whose excess equals the minimum
        if e_i == d_star:
            return i
        return self.fwdsearch(i, d_star - e_i)

    cpdef Py_ssize_t rMq(self, Py_ssize_t i, Py_ssize_t j) noexcept nogil:
        """The leftmost maximmum excess in i -> j.

        Symmetric to :meth:`rmq`: an O(log n) range-maximum over the rmM tree
        followed by a single ``fwdsearch`` for the leftmost attaining position.
        """
        cdef:
            Py_ssize_t m_star, p, blk_end, tree_v
            Py_ssize_t bi, bj, e_i
            Py_ssize_t b

        if i >= j:
            return i

        b = self._b
        bi = i // b
        bj = j // b

        # largest absolute excess over [i, j]
        e_i = self._e_index[i]
        m_star = e_i
        if bi == bj:
            for p in range(i + 1, j + 1):
                if self._e_index[p] > m_star:
                    m_star = self._e_index[p]
        else:
            blk_end = min((bi + 1) * b, self.size)
            for p in range(i + 1, blk_end):
                if self._e_index[p] > m_star:
                    m_star = self._e_index[p]
            if bi + 1 <= bj - 1:
                tree_v = self._rmq_tree_max(0, 0, self._n_internal,
                                            bi + 1, bj - 1)
                if tree_v > m_star:
                    m_star = tree_v
            for p in range(bj * b, j + 1):
                if self._e_index[p] > m_star:
                    m_star = self._e_index[p]

        # leftmost position in [i, j] whose excess equals the maximum
        if e_i == m_star:
            return i
        return self.fwdsearch(i, m_star - e_i)

    cdef Py_ssize_t _rmq_tree_min(self, Py_ssize_t node, Py_ssize_t node_lo, Py_ssize_t node_hi,
                              Py_ssize_t lo, Py_ssize_t hi) noexcept nogil:
        """Minimum excess over leaf-blocks [lo, hi].

        Segment-tree range query over the rmM binary tree; ``node`` spans the
        contiguous block range [node_lo, node_hi]. Only fully-covered nodes are
        read, and a covered node satisfies ``node_hi <= hi < n_tip``, so it
        (and its subtree) is always a real, populated block -- the partially
        filled last level is only ever touched through the partial end-blocks,
        which the callers scan directly.
        """
        cdef Py_ssize_t mid
        cdef Py_ssize_t left_v, right_v

        if hi < node_lo or node_hi < lo:  # no overlap
            return PY_SSIZE_T_MAX
        if lo <= node_lo and node_hi <= hi:  # fully covered
            return self._m[node]
        mid = (node_lo + node_hi) / 2
        left_v = self._rmq_tree_min(bt_left_child(node), node_lo, mid, lo, hi)
        right_v = self._rmq_tree_min(bt_right_child(node), mid + 1, node_hi,
                                     lo, hi)
        return left_v if left_v < right_v else right_v

    cdef Py_ssize_t _rmq_tree_max(self, Py_ssize_t node, Py_ssize_t node_lo, Py_ssize_t node_hi,
                              Py_ssize_t lo, Py_ssize_t hi) noexcept nogil:
        """Maximum excess over leaf-blocks [lo, hi].

        The range-maximum counterpart of :meth:`_rmq_tree_min`; see its note on
        why every fully-covered node is a real block.
        """
        cdef Py_ssize_t mid
        cdef Py_ssize_t left_v, right_v

        if hi < node_lo or node_hi < lo:  # no overlap
            return -PY_SSIZE_T_MAX
        if lo <= node_lo and node_hi <= hi:  # fully covered
            return self._M[node]
        mid = (node_lo + node_hi) / 2
        left_v = self._rmq_tree_max(bt_left_child(node), node_lo, mid, lo, hi)
        right_v = self._rmq_tree_max(bt_right_child(node), mid + 1, node_hi,
                                     lo, hi)
        return left_v if left_v > right_v else right_v

    cpdef Py_ssize_t depth(self, Py_ssize_t i) noexcept nogil:
        """The depth of node ``i``."""
        return self._e_index[i]

    cpdef Py_ssize_t root(self) noexcept nogil:
        """The index of the root node of the tree."""
        return 0

    cpdef Py_ssize_t parent(self, Py_ssize_t i) noexcept nogil:
        """The parent of node ``i``, or -1 for the root."""
        if i == self.root() or i == (self.size - 1):
            return -1
        else:
            return self.enclose(i)

    cpdef BOOL_t is_tip(self, Py_ssize_t i) noexcept nogil:
        """Whether node ``i`` is a tip."""
        return self._b_ptr[i] and (not self._b_ptr[i + 1])

    cpdef Py_ssize_t first_child(self, Py_ssize_t i) noexcept nogil:
        """Index of the first child of node ``i``, or 0 for a tip."""
        if self._b_ptr[i]:
            if self.is_tip(i):
                return 0
            else:
                return i + 1
        else:
            return self.first_child(self.open(i))

    cpdef Py_ssize_t last_child(self, Py_ssize_t i) noexcept nogil:
        """Index of the last child of node ``i``, or 0 for a tip."""
        if self._b_ptr[i]:
            if self.is_tip(i):
                return 0
            else:
                return self.open(self.close(i) - 1)
        else:
            return self.last_child(self.open(i))

    cpdef Py_ssize_t next_sibling(self, Py_ssize_t i) noexcept nogil:
        """Index of the next sibling of node ``i``, or 0 if none."""
        cdef Py_ssize_t pos

        if self._b_ptr[i]:
            pos = self.close(i) + 1
        else:
            pos = self.next_sibling(self.open(i))

        if pos >= self.size:
            return 0
        elif self._b_ptr[pos]:
            return pos
        else:
            return 0

    cpdef Py_ssize_t previous_sibling(self, Py_ssize_t i) noexcept nogil:
        """Index of the previous sibling of node ``i``, or 0 if none."""
        cdef Py_ssize_t pos

        if self._b_ptr[i]:
            if self._b_ptr[max(0, i - 1)]:
                return 0

            pos = self.open(i - 1)
        else:
            pos = self.previous_sibling(self.open(i))

        if pos < 0:
            return 0
        elif self._b_ptr[pos]:
            return pos
        else:
            return 0

    cpdef Py_ssize_t preorder_rank(self, Py_ssize_t i) noexcept nogil:
        """Preorder rank of node ``i``."""
        if self._b_ptr[i]:
            return self.rank(1, i)
        else:
            return self.preorder_rank(self.open(i))

    cpdef Py_ssize_t preorder_select(self, Py_ssize_t k) noexcept nogil:
        """Index of the node with preorder rank ``k``."""
        return self.select(1, k)

    cpdef Py_ssize_t postorder_rank(self, Py_ssize_t i) noexcept nogil:
        """Postorder rank of node ``i``."""
        if self._b_ptr[i]:
            return self.rank(0, self.close(i))
        else:
            return self.rank(0, i)

    cpdef Py_ssize_t postorder_select(self, Py_ssize_t k) noexcept nogil:
        """Index of the node with postorder rank ``k``."""
        return self.open(self.select(0, k))

    cpdef BOOL_t is_ancestor(self, Py_ssize_t i, Py_ssize_t j) noexcept nogil:
        """Whether node ``i`` is an ancestor of node ``j``.

        Either parenthesis names a node. Strictly inside ``i``'s span means a
        descendant; ``i``'s own parentheses, its two ends, are not, so a node
        is not its own ancestor whichever parenthesis names it.
        """
        if not self._b_ptr[i]:
            i = self.open(i)

        return i < j < self.close(i)

    cpdef Py_ssize_t count(self, Py_ssize_t i=0, bint tips=False) noexcept nogil:
        """Count of nodes (or tips) in the subtree rooted at node ``i``."""
        cdef:
            Py_ssize_t last, j, c

        if not self._b_ptr[i]:
            i = self.open(i)

        if not tips:
            return (self.close(i) - i + 1) / 2

        # a tip is an open parenthesis immediately followed by a close; count
        # them within the span [i, close(i)] of the subtree
        last = self.close(i)
        j = i
        c = 0
        while j < last:
            if self._b_ptr[j] and not self._b_ptr[j + 1]:
                c += 1
                j += 1
            j += 1

        return c

    cpdef Py_ssize_t level_ancestor(self, Py_ssize_t i, Py_ssize_t d) noexcept nogil:
        """Index of the ancestor ``d`` levels above node ``i``."""
        if d <= 0:
            return -1

        if not self._b_ptr[i]:
            i = self.open(i)

        return self.bwdsearch(i, -d - 1) + 1

    cpdef Py_ssize_t level_next(self, Py_ssize_t i) noexcept nogil:
        """Index of the next node at the same depth as node ``i``."""
        return self.fwdsearch(self.close(i), 1)

    cpdef Py_ssize_t lca(self, Py_ssize_t i, Py_ssize_t j) noexcept nogil:
        """The lowest common ancestor of nodes ``i`` and ``j``.

        Either parenthesis of a node names it, and the order of ``i`` and ``j``
        does not matter. The search below holds only for opening parentheses
        with ``i < j``, so both are brought to that form first; for the usual
        opening-parenthesis arguments that is two byte tests and a compare.
        """
        cdef Py_ssize_t t
        if not self._b_ptr[i]:
            i = self.open(i)
        if not self._b_ptr[j]:
            j = self.open(j)
        if i > j:
            t = i
            i = j
            j = t
        if i == j:
            # a node is its own lowest common ancestor; the search below would
            # return its parent
            return i
        if j < self.close(i):
            # i encloses j. The converse cannot hold: an ancestor opens before
            # its descendants, and j opens after i.
            return i
        return self.parent(self.rmq(i, j) + 1)

    cpdef Py_ssize_t deepest_node(self, Py_ssize_t i) noexcept nogil:
        """Index of the deepest node descending from node ``i``."""
        return self.rMq(self.open(i), self.close(i))

    cpdef Py_ssize_t height(self, Py_ssize_t i) noexcept nogil:
        """The height of node ``i``, in edges."""
        return self.excess(self.deepest_node(i)) - self.excess(self.open(i))

    cdef Py_ssize_t scan_block_forward(self, Py_ssize_t i, Py_ssize_t k, Py_ssize_t b, Py_ssize_t d) noexcept nogil:
        """Scan a block forward from i.

        Parameters
        ----------
        i : int
            The index position to start from in the tree
        k : int
            The block to explore
        b : int
            The block size
        d : int
            The depth to search for

        Returns
        -------
        int
            The index position of the result. -1 is returned if a result is not
            found.
        """
        cdef Py_ssize_t lower_bound
        cdef Py_ssize_t upper_bound
        cdef Py_ssize_t j

        # lower_bound is block boundary or right of i
        lower_bound = max(k, 0) * b
        lower_bound = max(i + 1, lower_bound)

        # upper_bound is block boundary or end of tree
        upper_bound = min((k + 1) * b, self.size)

        for j in range(lower_bound, upper_bound):
            if self._e_index[j] == d:
                return j

        return -1

    cdef Py_ssize_t scan_block_backward(self, Py_ssize_t i, Py_ssize_t k, Py_ssize_t b, Py_ssize_t d) noexcept nogil:
        """Scan a block backward from i.

        Parameters
        ----------
        i : int
            The index position to start from in the tree
        k : int
            The block to explore
        b : int
            The block size
        d : int
            The depth to search for

        Returns
        -------
        int
            The index position of the result. -1 is returned if a result is not
            found.
        """
        cdef Py_ssize_t lower_bound
        cdef Py_ssize_t upper_bound
        cdef Py_ssize_t j

        # range stop is exclusive, so need to set "stop" at -1 of boundary
        lower_bound = max(k, 0) * b - 1

        # include the right most position of the k-1 block so we can identify
        # closures spanning blocks.
        if lower_bound >= 0:
            lower_bound -= 1

        # upper bound is block boundary or left of i, whichever is less
        # note that this is an inclusive boundary since this is a backward search
        upper_bound = min((k + 1) * b, self.size) - 1
        upper_bound = min(i - 1, upper_bound)

        if upper_bound <= 0:
            return -1

        for j in range(upper_bound, lower_bound, -1):
            if self.excess(j) == d:
                return j

        return -1

    cdef Py_ssize_t fwdsearch(self, Py_ssize_t i, Py_ssize_t d) noexcept nogil:
        """Search forward from i for desired excess.

        Parameters
        ----------
        i : int
            The index to search forward from
        d : int
            The excess difference to search for (relative to E[i])

        Returns
        -------
        int
            The index of the result, or -1 if no result was found
        """
        cdef Py_ssize_t k  # the block being interrogated
        cdef Py_ssize_t result = -1 # the result of a scan within a block
        cdef Py_ssize_t node  # the node within the binary tree being examined

        # get the block of parentheses to check
        k = i // self._b

        # desired excess
        d += self._e_index[i]

        # determine which node our block corresponds too
        node = bt_node_from_left(k, self._height)

        # see if our result is in our current block
        if self._m[node] <= d <= self._M[node]:
            result = self.scan_block_forward(i, k, self._b, d)

        # if we do not have a result, we need to begin traversal of the tree
        if result == -1:
            # walk up the tree
            while not bt_is_root(node):
                if bt_is_left_child(node):
                    node = bt_right_sibling(node)
                    if self._m[node] <= d  <= self._M[node]:
                        break
                node = bt_parent(node)

            if bt_is_root(node):
                return -1

            # descend until we hit a leaf node
            while not bt_is_leaf(node, self._height):
                node = bt_left_child(node)

                # evaluate right, if not found, pick left
                if not (self._m[node] <= d <= self._M[node]):
                    node = bt_right_sibling(node)

            # we have found a block with contains our solution. convert from the
            # node index back into the block index
            k = node - <Py_ssize_t>(pow(2, self._height) - 1)

            # scan for a result using the original d
            result = self.scan_block_forward(i, k, self._b, d)

        return result

    cdef Py_ssize_t bwdsearch(self, Py_ssize_t i, Py_ssize_t d) noexcept nogil:
        """Search backward from i for desired excess

        Parameters
        ----------
        i : int
            The index to search forward from
        d : int
            The excess difference to search for (relative to E[i])

        Returns
        -------
        int
            The index of the result, or -1 if no result was found
        """
        cdef Py_ssize_t k  # the block being interrogated
        cdef Py_ssize_t result = -1 # the result of a scan within a block
        cdef Py_ssize_t node  # the node within the binary tree being examined

        # get the block of parentheses to check
        k = i // self._b

        # desired excess
        d += self.excess(i)

        # see if our result is in our current block
        result = self.scan_block_backward(i, k, self._b, d)

        # determine which node our block corresponds too
        node = bt_node_from_left(k, self._height)

        # special case: check sibling
        if result == -1 and bt_is_right_child(node):
            node = bt_left_sibling(node)
            k = node - <Py_ssize_t>(pow(2, self._height) - 1)
            result = self.scan_block_backward(i, k, self._b, d)

           # reset node and k in the event that result == -1
            k = i // self._b
            node = bt_right_sibling(node)

        # if we do not have a result, we need to begin traversal of the tree
        if result == -1:
            while not bt_is_root(node):
                # right nodes cannot contain the solution as we are searching left
                # As such, if we are the right node already, evaluate its sibling.
                if bt_is_right_child(node):
                    node = bt_left_sibling(node)
                    if self._m[node] <= d <= self._M[node]:
                        break

                # if we did not find a valid node, adjust for the relative
                # excess of the current node, and ascend to the parent
                node = bt_parent(node)

            if bt_is_root(node):
                return -1

            # descend until we hit a leaf node
            while not bt_is_leaf(node, self._height):
                node = bt_right_child(node)

                # evaluate right, if not found, pick left
                if not (self._m[node] <= d <= self._M[node]):
                    node = bt_left_sibling(node)

            # we have found a block with contains our solution. convert from the
            # node index back into the block index
            k = node - <Py_ssize_t>(pow(2, self._height) - 1)

            # scan for a result
            result = self.scan_block_backward(i, k, self._b, d)

        return result


# ---------------------------------------------------------------------------
# Whole-tree kernels
# ---------------------------------------------------------------------------


def build_index(const BOOL_t[::1] B, Py_ssize_t b, Py_ssize_t height):
    """Build the navigation index of a balanced-parentheses array.

    A single O(n) pass for the excess array, the select indexes and the
    leaves of the range min-max (rmM) tree of Navarro and Sadakane
    (http://www.dcc.uchile.cl/~gnavarro/ps/talg12.pdf) over blocks of the
    excess array, then one pass over the rmM tree's internal nodes.

    Parameters
    ----------
    B : numpy.ndarray of uint8
        The parentheses, 1 for an opening and 0 for a closing parenthesis.
    b : int
        rmM block size.
    height : int
        rmM tree height.

    Returns
    -------
    dict
        ``e_index`` (excess at each position), ``k_index_1`` and ``k_index_0``
        (select indexes: position of the k-th opening / closing parenthesis,
        with position 0 at k = 0 for the bit that is absent there), ``m``,
        ``M`` and ``r`` (minimum excess, maximum excess and rank of the rmM tree
        nodes in heap order), ``b`` (block size) and ``height`` (rmM tree
        height). Arrays are of intp.

    """
    cdef:
        Py_ssize_t n = B.shape[0]
        Py_ssize_t n_tip = (n + b - 1) // b
        Py_ssize_t n_internal = (<Py_ssize_t>1 << height) - 1
        Py_ssize_t n_total = n_tip + n_internal
        Py_ssize_t n_open, i, j, k, upper, lvl, pos, node, lchild, rchild
        Py_ssize_t excess = 0, rank = 0, min_, max_, ptr_0 = 0, ptr_1 = 0
        Py_ssize_t[::1] e_index, k_index_0, k_index_1, m, M, r
        BOOL_t first_open = B[0] != 0

    n_open = np.count_nonzero(B)
    e_index_arr = np.empty(n, dtype=SIZE)
    # position 0 is in both select indexes: the bit found there, plus a
    # leading 0 for the other bit
    k_index_1_arr = np.empty(n_open + (0 if first_open else 1), dtype=SIZE)
    k_index_0_arr = np.empty(n - n_open + (1 if first_open else 0), dtype=SIZE)
    m_arr = np.zeros(n_total, dtype=SIZE)
    M_arr = np.zeros(n_total, dtype=SIZE)
    r_arr = np.zeros(n_total, dtype=SIZE)
    e_index, k_index_0, k_index_1 = e_index_arr, k_index_0_arr, k_index_1_arr
    m, M, r = m_arr, M_arr, r_arr

    with nogil:
        if first_open:
            k_index_0[0] = 0
            ptr_0 = 1
        else:
            k_index_1[0] = 0
            ptr_1 = 1

        # leaves: one block of b parentheses each
        for k in range(n_tip):
            upper = min((k + 1) * b, n)
            min_ = PY_SSIZE_T_MAX
            max_ = 0
            r[n_internal + k] = rank
            for j in range(k * b, upper):
                # +1 for an opening and -1 for a closing parenthesis
                excess += -1 + (2 * B[j])
                rank += B[j]
                e_index[j] = excess
                if B[j]:
                    k_index_1[ptr_1] = j
                    ptr_1 += 1
                else:
                    k_index_0[ptr_0] = j
                    ptr_0 += 1
                if excess < min_:
                    min_ = excess
                if excess > max_:
                    max_ = excess
            m[n_internal + k] = min_
            M[n_internal + k] = max_

        # internal nodes, in reverse level order starting above the leaves
        for lvl in range(height - 1, -1, -1):
            for pos in range(<Py_ssize_t>1 << lvl):
                node = bt_node_from_left(pos, lvl)
                lchild = bt_left_child(node)
                rchild = bt_right_child(node)

                if lchild >= n_total:
                    continue
                elif rchild >= n_total:
                    m[node] = m[lchild]
                    M[node] = M[lchild]
                else:
                    m[node] = min(m[lchild], m[rchild])
                    M[node] = max(M[lchild], M[rchild])

                r[node] = r[lchild]

    return {
        "e_index": e_index_arr,
        "k_index_0": k_index_0_arr,
        "k_index_1": k_index_1_arr,
        "m": m_arr,
        "M": M_arr,
        "r": r_arr,
        "b": b,
        "height": height,
    }


def from_treenode_arrays(tree):
    """Balanced-parentheses arrays of a :class:`~skbio.tree.TreeNode`.

    Returns
    -------
    tuple of numpy.ndarray
        ``(topology, names, lengths, edges)``, each of length ``2 * n_nodes``.
    """
    cdef:
        Py_ssize_t n_nodes, ptr
        cnp.ndarray[BOOL_t, ndim=1] topo
        cnp.ndarray[object, ndim=1] names
        cnp.ndarray[DOUBLE_t, ndim=1] lengths
        cnp.ndarray[INT32_t, ndim=1] edges
        set seen

    n_nodes = len(list(tree.traverse(include_self=True)))

    topo = np.zeros(n_nodes * 2, dtype=np.uint8)
    names = np.full(n_nodes * 2, None, dtype=object)
    lengths = np.zeros(n_nodes * 2, dtype=np.double)
    edges = np.zeros(n_nodes * 2, dtype=np.int32)

    ptr = 0
    seen = set()
    for n in tree.pre_and_postorder(include_self=True):
        if n not in seen:
            topo[ptr] = 1
            names[ptr] = n.name
            lengths[ptr] = n.length or 0.0
            edges[ptr] = getattr(n, 'edge_num', None) or 0

            if n.is_tip():
                ptr += 1

            seen.add(n)

        ptr += 1
    return topo, names, lengths, edges


def to_array(_BPKernel k):
    """Kernel of :meth:`skbio.tree.BPTree.to_array`."""
    cdef:
        Py_ssize_t i, n
        Py_ssize_t chi_ptr, cur_index
        Py_ssize_t node_idx, first_child, last_child, sib_idx
        cnp.ndarray[DOUBLE_t, ndim=1] length
        cnp.ndarray[UINT32_t, ndim=1] node_ids
        cnp.ndarray[object, ndim=1] name
        dict id_index

    class mock_node:
        def __init__(self, id, is_tip):
            self.is_tip_ = is_tip
            self.id = id

        def is_tip(self):
            return self.is_tip_

    n = k.size // 2

    child_index = np.zeros((n - k.count(tips=True), 3), dtype=np.int64)
    length = np.zeros(n, dtype=np.double)
    node_ids = np.zeros(k.size, dtype=np.uint32)
    name = np.full(n, None, dtype=object)

    # TreeNode.assign_ids, decompose target
    chi_ptr = 0
    cur_index = 0  # the index into node_ids, equivalent to TreeNode.assign_ids
    id_index = dict.fromkeys(set(range(n)))  # map a node's "id" to an object which indicates if it is a leaf or not
    for i in range(n):
        node_idx = k.postorder_select(i + 1)  # the index within the BP of the node

        if not k.is_tip(node_idx):
            first_child = k.first_child(node_idx)
            last_child = k.last_child(node_idx)

            sib_idx = first_child  # the sibling index wtihin the BP of the node
            while sib_idx != 0 and sib_idx <= last_child:
                node_ids[sib_idx] = cur_index
                id_index[cur_index] = mock_node(cur_index, k.is_tip(sib_idx))
                length[cur_index] = k.length(sib_idx)
                name[cur_index] = k.name(sib_idx)

                cur_index += 1
                sib_idx = k.next_sibling(sib_idx)

            child_index[chi_ptr] = [node_idx, node_ids[first_child], node_ids[last_child]]
            chi_ptr += 1

    # make sure to capture root
    id_index[n - 1] = mock_node(cur_index, False)

    node_ids[0] = cur_index
    child_index[:, 0] = node_ids[child_index[:, 0]]
    child_index = child_index[np.argsort(child_index[:, 0])]

    return {'child_index': child_index, 'length': length, 'id_index': id_index,
            'name': name}


def to_node_arrays(_BPKernel k):
    """Kernel of :meth:`skbio.tree.BPTree._to_node_arrays`."""
    cdef:
        Py_ssize_t i, n
        Py_ssize_t node_idx, root
        cnp.ndarray[object, ndim=1] name
        cnp.ndarray[DOUBLE_t, ndim=1] length
        cnp.ndarray[INT32_t, ndim=1] edge
        cnp.ndarray[Py_ssize_t, ndim=1] parent

    n = k.size // 2
    name = np.empty(n, dtype=object)
    length = np.empty(n, dtype=DOUBLE)
    edge = np.empty(n, dtype=INT32)
    parent = np.empty(n, dtype=SIZE)

    root = k.root()
    for i in range(n):
        node_idx = k.preorder_select(i)
        name[i] = k.name(node_idx)
        length[i] = k.length(node_idx)
        edge[i] = k.edge(node_idx)
        if node_idx == root:
            parent[i] = -1
        else:
            # preorder_rank is 1-based, so -1 converts it to a 0-based index
            parent[i] = k.preorder_rank(k.parent(node_idx)) - 1

    return name, length, edge, parent


def edge_lookup(const BOOL_t[::1] B, const INT32_t[::1] edges):
    """Map each edge number to the position of the node that carries it."""
    cdef:
        Py_ssize_t i, n = B.shape[0]
        cnp.ndarray[Py_ssize_t, ndim=1] lookup

    lookup = np.full(n, 0, dtype=SIZE)
    for i in range(n):
        if B[i] == 1:
            lookup[edges[i]] = i
    return lookup


def shear_mask(_BPKernel k, set tips):
    """Mask of the positions kept by :meth:`skbio.tree.BPTree.shear`.

    Returns
    -------
    tuple
        ``(mask, count)``: a uint8 mask over the parentheses, set for every
        retained tip and its ancestors (both parentheses of each), and the
        number of requested tips found.
    """
    cdef:
        Py_ssize_t i, p, count = 0
        cnp.ndarray[BOOL_t, ndim=1] mask_arr
        BOOL_t* mask

    mask_arr = np.zeros(k.size, dtype=BOOL)
    mask = <BOOL_t*>mask_arr.data
    mask[k.root()] = 1
    mask[k.close(k.root())] = 1

    for i in range(k.size):
        # is_tip is only defined on the open parenthesis
        if k.is_tip(i):
            if k.name(i) in tips:  # gil is required for set operation
                with nogil:
                    count += 1
                    mask[i] = 1
                    mask[i + 1] = 1

                    p = k.parent(i)
                    while p != 0 and mask[p] == 0:
                        mask[p] = 1
                        mask[k.close(p)] = 1

                        p = k.parent(p)

    return mask_arr, count


def collapse_mask(_BPKernel k):
    """Mask and merged lengths for :meth:`skbio.tree.BPTree.collapse`.

    Returns
    -------
    tuple of numpy.ndarray
        ``(mask, lengths)``: a uint8 mask over the parentheses, set for every
        retained node (both parentheses), and a copy of the branch lengths in
        which each removed single-child node's length is added to its child.
    """
    cdef:
        Py_ssize_t i, n = k.size // 2
        Py_ssize_t current, first, last
        cnp.ndarray[BOOL_t, ndim=1] mask_arr
        cnp.ndarray[DOUBLE_t, ndim=1] new_lengths
        BOOL_t* mask
        DOUBLE_t* new_lengths_ptr

    mask_arr = np.zeros(k.size, dtype=BOOL)
    mask = <BOOL_t*>mask_arr.data
    mask[k.root()] = 1
    mask[k.close(k.root())] = 1

    new_lengths = k._lengths.copy()
    new_lengths_ptr = <DOUBLE_t*>new_lengths.data

    with nogil:
        for i in range(n):
            current = k.preorder_select(i)

            if k.is_tip(current):
                mask[current] = 1
                mask[k.close(current)] = 1
            else:
                first = k.first_child(current)
                last = k.last_child(current)

                if first == last:
                    new_lengths_ptr[first] = new_lengths_ptr[first] + \
                            new_lengths_ptr[current]
                else:
                    mask[current] = 1
                    mask[k.close(current)] = 1

    return mask_arr, new_lengths


# ---------------------------------------------------------------------------
# Batch kernels (engine="cython")
#
# Each runs the compiled navigation primitives over an array of queries in an
# OpenMP ``prange`` without the GIL; the queries are independent. Inputs are
# validated by the ``BPTree`` methods that call these.
# ---------------------------------------------------------------------------


def close_batch(_BPKernel k, const Py_ssize_t[::1] idx):
    """Kernel of :meth:`skbio.tree.BPTree.close_batch`."""
    cdef:
        Py_ssize_t t, n = idx.shape[0]
        Py_ssize_t[::1] out
    out_arr = np.empty(n, dtype=SIZE)
    out = out_arr
    for t in prange(n, nogil=True):
        out[t] = k.close(idx[t])
    return out_arr


def parent_batch(_BPKernel k, const Py_ssize_t[::1] idx):
    """Kernel of :meth:`skbio.tree.BPTree.parent_batch`."""
    cdef:
        Py_ssize_t t, n = idx.shape[0]
        Py_ssize_t[::1] out
    out_arr = np.empty(n, dtype=SIZE)
    out = out_arr
    for t in prange(n, nogil=True):
        out[t] = k.parent(idx[t])
    return out_arr


def lca_batch(_BPKernel k, const Py_ssize_t[::1] i, const Py_ssize_t[::1] j):
    """Kernel of :meth:`skbio.tree.BPTree.lca_batch`."""
    cdef:
        Py_ssize_t t, n = i.shape[0]
        Py_ssize_t[::1] out
    out_arr = np.empty(n, dtype=SIZE)
    out = out_arr
    for t in prange(n, nogil=True):
        out[t] = k.lca(i[t], j[t])
    return out_arr


def level_ancestor_batch(_BPKernel k, const Py_ssize_t[::1] idx,
                        const Py_ssize_t[::1] d):
    """Kernel of :meth:`skbio.tree.BPTree.level_ancestor_batch`."""
    cdef:
        Py_ssize_t t, n = idx.shape[0]
        Py_ssize_t[::1] out
    out_arr = np.empty(n, dtype=SIZE)
    out = out_arr
    for t in prange(n, nogil=True):
        out[t] = k.level_ancestor(idx[t], d[t])
    return out_arr


def root_distances(_BPKernel k, const DOUBLE_t[::1] lengths):
    """Sum of branch lengths from the root to each node.

    Indexed by position; set at opening parentheses (0 elsewhere). The root's
    own length is not included. Each node's value is its parent's plus its own
    length, so every value is an exact root-to-node path sum.
    """
    cdef:
        Py_ssize_t i, top = 0, n = k.size
        DOUBLE_t[::1] out, stack
    out_arr = np.zeros(n, dtype=DOUBLE)
    out = out_arr
    # distance of each open ancestor; the root is at depth 0
    stack = np.zeros(n // 2 + 1, dtype=DOUBLE)
    with nogil:
        for i in range(1, n):
            if k._b_ptr[i]:
                top += 1
                stack[top] = stack[top - 1] + lengths[i]
                out[i] = stack[top]
            else:
                top -= 1
    return out_arr


def tip_distances(const Py_ssize_t[::1] tips, const Py_ssize_t[::1] slot,
                  const Py_ssize_t[::1] parent, const Py_ssize_t[::1] end,
                  const DOUBLE_t[::1] dist):
    """Pairwise path distances between tips, as a square matrix.

    Parameters
    ----------
    tips : array of intp
        Positions (opening parentheses) of the tips, in tree order.
    slot : array of intp
        Row (and column) of the output for each tip of ``tips``.
    parent : array of intp
        Per-position parent position (at opening parentheses; -1 for the root).
    end : array of intp
        Per-position (at opening parentheses) exclusive end of the node's tips
        in ``tips``: a node's tips are ``tips[start:end[v]]`` for some start.
    dist : array of float64
        Per-position distance from the root (branch lengths or edge counts).

    Returns
    -------
    numpy.ndarray of float64, shape (n, n)
        The symmetric distance matrix, zero on the diagonal.

    Notes
    -----
    The tips of any node are contiguous in tree order. For tip ``a``, walking
    up its ancestors, the tips after ``a`` within ancestor ``v`` but outside
    the child of ``v`` containing ``a`` are exactly those whose lowest common
    ancestor with ``a`` is ``v``: their distances are filled as one run, with
    no search per pair. Each row is independent; rows ``h`` and ``n - 2 - h``
    share a ``prange`` iteration to balance the triangular workload.
    """
    cdef:
        Py_ssize_t n = tips.shape[0]
        Py_ssize_t n_rows = n - 1 if n > 0 else 0
        Py_ssize_t h
        DOUBLE_t[:, ::1] out
    out_arr = np.empty((n, n), dtype=DOUBLE)
    out = out_arr
    for h in prange((n_rows + 1) // 2, nogil=True, schedule='dynamic'):
        _tip_distance_row(tips, slot, parent, end, dist, out, h, n)
        if n_rows - 1 - h != h:
            _tip_distance_row(tips, slot, parent, end, dist, out,
                              n_rows - 1 - h, n)
    if n:
        out[slot[n - 1], slot[n - 1]] = 0
    return out_arr


cdef inline void _tip_distance_row(const Py_ssize_t[::1] tips,
                                   const Py_ssize_t[::1] slot,
                                   const Py_ssize_t[::1] parent,
                                   const Py_ssize_t[::1] end,
                                   const DOUBLE_t[::1] dist,
                                   DOUBLE_t[:, ::1] out,
                                   Py_ssize_t row, Py_ssize_t n) noexcept nogil:
    """Fill the pairs of tip ``row`` with the tips after it."""
    cdef:
        Py_ssize_t col, v, stop, r, c
        Py_ssize_t start = row + 1
        DOUBLE_t da, dv, d
    r = slot[row]
    out[r, r] = 0
    da = dist[tips[row]]
    v = parent[tips[row]]
    while start < n:
        # tips start .. end[v] - 1 have lowest common ancestor v with this tip
        stop = end[v]
        dv = dist[v]
        for col in range(start, stop):
            # two differences (no multiply), so compilers cannot contract it
            # into a fused multiply-add and every engine rounds identically
            d = (da - dv) + (dist[tips[col]] - dv)
            c = slot[col]
            out[r, c] = d
            out[c, r] = d
        start = stop
        v = parent[v]
