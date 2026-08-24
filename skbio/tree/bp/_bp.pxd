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

from ._ba cimport BIT_ARRAY

ctypedef cnp.uint32_t UINT32_t
ctypedef cnp.int32_t INT32_t
ctypedef cnp.float64_t DOUBLE_t
ctypedef cnp.uint8_t BOOL_t


cdef class mM:
    cdef Py_ssize_t b  # block size (Py_ssize_t so block-index * block-size products stay 64-bit)
    cdef int n_tip  # number of tips in the binary tree
    cdef int n_internal  # number of internal nodes in the binary tree
    cdef int n_total  # total number of nodes in the binary tree
    cdef int height  # the height of the binary tree
    cdef int m_idx  # m is minimum excess
    cdef int M_idx  # M is maximum excess
    cdef int r_idx  # rank
    cdef Py_ssize_t[:, ::1] mM
    cdef Py_ssize_t[:] r

    cdef void rmm(self, BOOL_t[:] B, Py_ssize_t B_size) nogil


@cython.final
cdef class BPTree:
    cdef:
        public cnp.ndarray data
        BOOL_t* _b_ptr
        Py_ssize_t[:] _e_index
        Py_ssize_t[:] _k_index_0
        Py_ssize_t[:] _k_index_1
        cnp.ndarray _names
        cnp.ndarray _lengths
        cnp.ndarray _edges
        cnp.ndarray _edge_lookup
        mM _rmm
        Py_ssize_t size

    cdef inline Py_ssize_t rank(self, Py_ssize_t t, Py_ssize_t i) nogil
    cdef inline Py_ssize_t select(self, Py_ssize_t t, Py_ssize_t k) nogil
    cdef Py_ssize_t _excess(self, Py_ssize_t i) nogil
    cdef Py_ssize_t excess(self, Py_ssize_t i) nogil
    cdef Py_ssize_t fwdsearch(self, Py_ssize_t i, Py_ssize_t d) nogil
    cdef Py_ssize_t bwdsearch(self, Py_ssize_t i, Py_ssize_t d) nogil
    cpdef inline Py_ssize_t close(self, Py_ssize_t i) nogil
    cdef inline Py_ssize_t open(self, Py_ssize_t i) nogil
    cpdef inline BOOL_t is_tip(self, Py_ssize_t i) nogil
    cdef inline Py_ssize_t enclose(self, Py_ssize_t i) nogil
    cdef BPTree _mask_from_self(self, BIT_ARRAY* mask, cnp.ndarray[DOUBLE_t, ndim=1] lengths)
    cpdef Py_ssize_t next_sibling(self, Py_ssize_t i) nogil
    cpdef Py_ssize_t previous_sibling(self, Py_ssize_t i) nogil
    cpdef Py_ssize_t last_child(self, Py_ssize_t i) nogil
    cpdef Py_ssize_t first_child(self, Py_ssize_t i) nogil
    cpdef Py_ssize_t parent(self, Py_ssize_t i) nogil
    cpdef Py_ssize_t depth(self, Py_ssize_t i) nogil
    cpdef Py_ssize_t root(self) nogil
    cdef Py_ssize_t scan_block_forward(self, Py_ssize_t i, int k, Py_ssize_t b, Py_ssize_t d) nogil
    cdef Py_ssize_t scan_block_backward(self, Py_ssize_t i, int k, Py_ssize_t b, Py_ssize_t d) nogil
    cdef void _set_edges(self, cnp.ndarray[INT32_t, ndim=1] edges)

    cpdef inline unicode name(self, Py_ssize_t i)
    cpdef inline DOUBLE_t length(self, Py_ssize_t i)
    cpdef inline INT32_t edge(self, Py_ssize_t i)
    cpdef Py_ssize_t edge_from_number(self, INT32_t n)
    cpdef Py_ssize_t rmq(self, Py_ssize_t i, Py_ssize_t j) nogil
    cpdef Py_ssize_t rMq(self, Py_ssize_t i, Py_ssize_t j) nogil
    cdef Py_ssize_t _rmq_tree_min(self, int node, int node_lo, int node_hi, int lo, int hi) nogil
    cdef Py_ssize_t _rmq_tree_max(self, int node, int node_lo, int node_hi, int lo, int hi) nogil
    cpdef Py_ssize_t postorder_select(self, Py_ssize_t k) nogil
    cpdef Py_ssize_t postorder_rank(self, Py_ssize_t i) nogil
    cpdef Py_ssize_t preorder_select(self, Py_ssize_t k) nogil
    cpdef Py_ssize_t preorder_rank(self, Py_ssize_t i) nogil
    cpdef BOOL_t is_ancestor(self, Py_ssize_t i, Py_ssize_t j) nogil
    cpdef Py_ssize_t level_ancestor(self, Py_ssize_t i, Py_ssize_t d) nogil
    cpdef Py_ssize_t count(self, Py_ssize_t i=*, bint tips=*) nogil
    cpdef BPTree shear(self, set tips)
    cpdef BPTree collapse(self)
    cpdef Py_ssize_t level_next(self, Py_ssize_t i) nogil
    cpdef Py_ssize_t height(self, Py_ssize_t i) nogil
    cpdef Py_ssize_t deepest_node(self, Py_ssize_t i) nogil
    cpdef Py_ssize_t lca(self, Py_ssize_t i, Py_ssize_t j) nogil
