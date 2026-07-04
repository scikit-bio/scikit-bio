cimport numpy as cnp
cimport cython

from ._ba cimport BIT_ARRAY

ctypedef cnp.npy_intp SIZE_t
ctypedef cnp.uint32_t UINT32_t
ctypedef cnp.int32_t INT32_t
ctypedef cnp.float64_t DOUBLE_t
ctypedef cnp.uint8_t BOOL_t


cdef class mM:
    cdef int b  # block size
    cdef int n_tip  # number of tips in the binary tree
    cdef int n_internal  # number of internal nodes in the binary tree
    cdef int n_total  # total number of nodes in the binary tree
    cdef int height  # the height of the binary tree
    cdef int m_idx  # m is minimum excess
    cdef int M_idx  # M is maximum excess
    cdef int r_idx  # rank
    cdef SIZE_t[:, ::1] mM
    cdef SIZE_t[:] r

    cdef void rmm(self, BOOL_t[:] B, int B_size) nogil


@cython.final
cdef class BPTree:
    cdef:
        public cnp.ndarray data
        BOOL_t* _b_ptr
        SIZE_t[:] _e_index
        SIZE_t[:] _k_index_0
        SIZE_t[:] _k_index_1
        cnp.ndarray _names
        cnp.ndarray _lengths
        cnp.ndarray _edges
        cnp.ndarray _edge_lookup
        mM _rmm
        SIZE_t size

    cdef inline SIZE_t rank(self, SIZE_t t, SIZE_t i) nogil
    cdef inline SIZE_t select(self, SIZE_t t, SIZE_t k) nogil
    cdef SIZE_t _excess(self, SIZE_t i) nogil
    cdef SIZE_t excess(self, SIZE_t i) nogil
    cdef SIZE_t fwdsearch(self, SIZE_t i, int d) nogil
    cdef SIZE_t bwdsearch(self, SIZE_t i, int d) nogil
    cpdef inline SIZE_t close(self, SIZE_t i) nogil
    cdef inline SIZE_t open(self, SIZE_t i) nogil
    cpdef inline BOOL_t is_tip(self, SIZE_t i) nogil
    cdef inline SIZE_t enclose(self, SIZE_t i) nogil
    cdef BPTree _mask_from_self(self, BIT_ARRAY* mask, cnp.ndarray[DOUBLE_t, ndim=1] lengths)
    cpdef SIZE_t next_sibling(self, SIZE_t i) nogil
    cpdef SIZE_t previous_sibling(self, SIZE_t i) nogil
    cpdef SIZE_t last_child(self, SIZE_t i) nogil
    cpdef SIZE_t first_child(self, SIZE_t i) nogil
    cpdef SIZE_t parent(self, SIZE_t i) nogil
    cpdef SIZE_t depth(self, SIZE_t i) nogil
    cpdef SIZE_t root(self) nogil
    cdef int scan_block_forward(self, int i, int k, int b, int d) nogil
    cdef int scan_block_backward(self, int i, int k, int b, int d) nogil
    cdef void _set_edges(self, cnp.ndarray[INT32_t, ndim=1] edges)

    cpdef inline unicode name(self, SIZE_t i)
    cpdef inline DOUBLE_t length(self, SIZE_t i)
    cpdef inline INT32_t edge(self, SIZE_t i)
    cpdef SIZE_t edge_from_number(self, INT32_t n)
    cpdef SIZE_t rmq(self, SIZE_t i, SIZE_t j) nogil
    cpdef SIZE_t rMq(self, SIZE_t i, SIZE_t j) nogil
    cpdef SIZE_t postorder_select(self, SIZE_t k) nogil
    cpdef SIZE_t postorder_rank(self, SIZE_t i) nogil
    cpdef SIZE_t preorder_select(self, SIZE_t k) nogil
    cpdef SIZE_t preorder_rank(self, SIZE_t i) nogil
    cpdef BOOL_t is_ancestor(self, SIZE_t i, SIZE_t j) nogil
    cpdef SIZE_t level_ancestor(self, SIZE_t i, SIZE_t d) nogil
    cpdef SIZE_t subtree_size(self, SIZE_t i) nogil
    cpdef BPTree shear(self, set tips)
    cpdef BPTree collapse(self)
    cpdef SIZE_t ntips(self) nogil
    cpdef SIZE_t level_next(self, SIZE_t i) nogil
    cpdef SIZE_t height(self, SIZE_t i) nogil
    cpdef SIZE_t deepest_node(self, SIZE_t i) nogil
    cpdef SIZE_t lca(self, SIZE_t i, SIZE_t j) nogil
