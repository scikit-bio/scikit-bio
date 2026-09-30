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

cimport numpy as cnp
cimport cython

ctypedef cnp.uint32_t UINT32_t
ctypedef cnp.int32_t INT32_t
ctypedef cnp.float64_t DOUBLE_t
ctypedef cnp.uint8_t BOOL_t


@cython.final
cdef class _BPKernel:
    cdef:
        # read-only views onto the arrays owned by the Python ``BPTree``
        const BOOL_t[::1] _B
        const BOOL_t* _b_ptr
        const Py_ssize_t[::1] _e_index
        const Py_ssize_t[::1] _k_index_0
        const Py_ssize_t[::1] _k_index_1
        # range min-max (rmM) tree, in heap (breadth-first) order
        const Py_ssize_t[::1] _m
        const Py_ssize_t[::1] _M
        const Py_ssize_t[::1] _r
        Py_ssize_t _b  # rmM block size
        Py_ssize_t _height  # rmM tree height
        Py_ssize_t _n_internal  # rmM internal node count
        Py_ssize_t size
        # node attributes, replaced in place by ``BPTree.set_*``
        public cnp.ndarray _names
        public cnp.ndarray _lengths
        public cnp.ndarray _edges
        public object _edge_lookup

    cdef inline Py_ssize_t rank(self, Py_ssize_t t, Py_ssize_t i) noexcept nogil
    cdef inline Py_ssize_t select(self, Py_ssize_t t, Py_ssize_t k) noexcept nogil
    cdef Py_ssize_t excess(self, Py_ssize_t i) noexcept nogil
    cdef Py_ssize_t fwdsearch(self, Py_ssize_t i, Py_ssize_t d) noexcept nogil
    cdef Py_ssize_t bwdsearch(self, Py_ssize_t i, Py_ssize_t d) noexcept nogil
    cpdef inline Py_ssize_t close(self, Py_ssize_t i) noexcept nogil
    cdef inline Py_ssize_t open(self, Py_ssize_t i) noexcept nogil
    cpdef inline BOOL_t is_tip(self, Py_ssize_t i) noexcept nogil
    cdef inline Py_ssize_t enclose(self, Py_ssize_t i) noexcept nogil
    cpdef Py_ssize_t next_sibling(self, Py_ssize_t i) noexcept nogil
    cpdef Py_ssize_t previous_sibling(self, Py_ssize_t i) noexcept nogil
    cpdef Py_ssize_t last_child(self, Py_ssize_t i) noexcept nogil
    cpdef Py_ssize_t first_child(self, Py_ssize_t i) noexcept nogil
    cpdef Py_ssize_t parent(self, Py_ssize_t i) noexcept nogil
    cpdef Py_ssize_t depth(self, Py_ssize_t i) noexcept nogil
    cpdef Py_ssize_t root(self) noexcept nogil
    cdef Py_ssize_t scan_block_forward(self, Py_ssize_t i, Py_ssize_t k, Py_ssize_t b, Py_ssize_t d) noexcept nogil
    cdef Py_ssize_t scan_block_backward(self, Py_ssize_t i, Py_ssize_t k, Py_ssize_t b, Py_ssize_t d) noexcept nogil

    cpdef inline unicode name(self, Py_ssize_t i)
    cpdef inline DOUBLE_t length(self, Py_ssize_t i)
    cpdef inline INT32_t edge(self, Py_ssize_t i)
    cpdef Py_ssize_t edge_from_number(self, INT32_t n)
    cpdef Py_ssize_t rmq(self, Py_ssize_t i, Py_ssize_t j) noexcept nogil
    cpdef Py_ssize_t rMq(self, Py_ssize_t i, Py_ssize_t j) noexcept nogil
    cdef Py_ssize_t _rmq_tree_min(self, Py_ssize_t node, Py_ssize_t node_lo, Py_ssize_t node_hi, Py_ssize_t lo, Py_ssize_t hi) noexcept nogil
    cdef Py_ssize_t _rmq_tree_max(self, Py_ssize_t node, Py_ssize_t node_lo, Py_ssize_t node_hi, Py_ssize_t lo, Py_ssize_t hi) noexcept nogil
    cpdef Py_ssize_t postorder_select(self, Py_ssize_t k) noexcept nogil
    cpdef Py_ssize_t postorder_rank(self, Py_ssize_t i) noexcept nogil
    cpdef Py_ssize_t preorder_select(self, Py_ssize_t k) noexcept nogil
    cpdef Py_ssize_t preorder_rank(self, Py_ssize_t i) noexcept nogil
    cpdef BOOL_t is_ancestor(self, Py_ssize_t i, Py_ssize_t j) noexcept nogil
    cpdef Py_ssize_t level_ancestor(self, Py_ssize_t i, Py_ssize_t d) noexcept nogil
    cpdef Py_ssize_t count(self, Py_ssize_t i=*, bint tips=*) noexcept nogil
    cpdef Py_ssize_t level_next(self, Py_ssize_t i) noexcept nogil
    cpdef Py_ssize_t height(self, Py_ssize_t i) noexcept nogil
    cpdef Py_ssize_t deepest_node(self, Py_ssize_t i) noexcept nogil
    cpdef Py_ssize_t lca(self, Py_ssize_t i, Py_ssize_t j) noexcept nogil
