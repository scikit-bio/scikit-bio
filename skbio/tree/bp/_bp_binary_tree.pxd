# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

# An implementation of a complete binary tree in breadth first order adapted
# from https://github.com/jfuentess/sea2015/blob/master/binary_trees.h

from libc.math cimport pow, log2, floor
from ._bp cimport SIZE_t


cdef inline SIZE_t bt_is_root(SIZE_t v) nogil:
    """Is v the root"""
    return v == 0


cdef inline SIZE_t bt_is_left_child(SIZE_t v) nogil:
    """Is v a left child of some node"""
    return 0 if bt_is_root(v) else v % 2


cdef inline SIZE_t bt_is_right_child(SIZE_t v) nogil:
    """Is v a right child of some node"""
    return 0 if bt_is_root(v) else 1 - (v % 2)


cdef inline SIZE_t bt_parent(SIZE_t v) nogil:
    """Get the index of the parent of v"""
    return 0 if bt_is_root(v) else (v - 1) // 2


cdef inline SIZE_t bt_left_child(SIZE_t v) nogil:
    """Get the index of the left child of v"""
    return 2 * v + 1


cdef inline SIZE_t bt_right_child(SIZE_t v) nogil:
    """Get the index of the right child of v"""
    return 2 * v + 2


cdef inline SIZE_t bt_left_sibling(SIZE_t v) nogil:
    """Get the index of the left sibling of v"""
    return v - 1


cdef inline SIZE_t bt_right_sibling(SIZE_t v) nogil:
    """Get the index of the right sibling of v"""
    return v + 1


cdef inline SIZE_t bt_is_leaf(SIZE_t v, SIZE_t height) nogil:
    """Determine if v is a leaf"""
    return <SIZE_t>(v >= pow(2, height) - 1)


cdef inline SIZE_t bt_node_from_left(SIZE_t pos, SIZE_t height) nogil:
    """Get the index from the left of a node at a given height"""
    return <SIZE_t>pow(2, height) - 1 + pos


cdef inline SIZE_t bt_offset_from_left(SIZE_t v) nogil:
    """Get the position from left of a node at its level

    This is the inverse of bt_node_from_left
    """
    cdef double leftmost_check

    if bt_is_root(v):
        return 0

    leftmost_check = log2(v + 1)
    if leftmost_check == floor(leftmost_check):
        return 0

    return v - <SIZE_t>pow(2, floor(log2(v))) + 1


cdef inline SIZE_t bt_offset_from_right(SIZE_t v) nogil:
    """Get the position from right of a node at its level"""
    cdef SIZE_t lvl = <SIZE_t>floor(log2(v + 1))
    cdef SIZE_t n_nodes_at_lvl = <SIZE_t>pow(2, lvl)

    return n_nodes_at_lvl - bt_offset_from_left(v) - 1


cdef inline SIZE_t bt_left_leaf(SIZE_t v, SIZE_t height) nogil:
    """Determine the index of a nodes left most leaf"""
    cdef SIZE_t left_tip = <SIZE_t>pow(2, height) - 1
    cdef SIZE_t block_size

    if bt_is_root(v):
        return left_tip

    block_size = <SIZE_t>pow(2, height - floor(log2(v + 1)))

    return left_tip + (block_size * bt_offset_from_left(v))


cdef inline SIZE_t bt_right_leaf(SIZE_t v, SIZE_t height) nogil:
    """Determine the index of a nodes right most leaf"""
    cdef SIZE_t right_tip = <SIZE_t>pow(2, height + 1) - 2
    cdef SIZE_t block_size

    if bt_is_root(v):
        return right_tip

    block_size = <SIZE_t>pow(2, height - floor(log2(v + 1)))

    return right_tip - (block_size * bt_offset_from_right(v))
