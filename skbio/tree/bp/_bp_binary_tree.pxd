# cython: cdivision=True, boundscheck=False, wraparound=False
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

# An implementation of a complete binary tree in breadth first order adapted
# from https://github.com/jfuentess/sea2015/blob/master/binary_trees.h

from libc.math cimport pow, log2, floor


cdef inline Py_ssize_t bt_is_root(Py_ssize_t v) nogil:
    """Is v the root"""
    return v == 0


cdef inline Py_ssize_t bt_is_left_child(Py_ssize_t v) nogil:
    """Is v a left child of some node"""
    return 0 if bt_is_root(v) else v % 2


cdef inline Py_ssize_t bt_is_right_child(Py_ssize_t v) nogil:
    """Is v a right child of some node"""
    return 0 if bt_is_root(v) else 1 - (v % 2)


cdef inline Py_ssize_t bt_parent(Py_ssize_t v) nogil:
    """Get the index of the parent of v"""
    return 0 if bt_is_root(v) else (v - 1) // 2


cdef inline Py_ssize_t bt_left_child(Py_ssize_t v) nogil:
    """Get the index of the left child of v"""
    return 2 * v + 1


cdef inline Py_ssize_t bt_right_child(Py_ssize_t v) nogil:
    """Get the index of the right child of v"""
    return 2 * v + 2


cdef inline Py_ssize_t bt_left_sibling(Py_ssize_t v) nogil:
    """Get the index of the left sibling of v"""
    return v - 1


cdef inline Py_ssize_t bt_right_sibling(Py_ssize_t v) nogil:
    """Get the index of the right sibling of v"""
    return v + 1


cdef inline Py_ssize_t bt_is_leaf(Py_ssize_t v, Py_ssize_t height) nogil:
    """Determine if v is a leaf"""
    return <Py_ssize_t>(v >= pow(2, height) - 1)


cdef inline Py_ssize_t bt_node_from_left(Py_ssize_t pos, Py_ssize_t height) nogil:
    """Get the index from the left of a node at a given height"""
    return <Py_ssize_t>pow(2, height) - 1 + pos


cdef inline Py_ssize_t bt_offset_from_left(Py_ssize_t v) nogil:
    """Get the position from left of a node at its level

    This is the inverse of bt_node_from_left
    """
    cdef double leftmost_check

    if bt_is_root(v):
        return 0

    leftmost_check = log2(v + 1)
    if leftmost_check == floor(leftmost_check):
        return 0

    return v - <Py_ssize_t>pow(2, floor(log2(v))) + 1


cdef inline Py_ssize_t bt_offset_from_right(Py_ssize_t v) nogil:
    """Get the position from right of a node at its level"""
    cdef Py_ssize_t lvl = <Py_ssize_t>floor(log2(v + 1))
    cdef Py_ssize_t n_nodes_at_lvl = <Py_ssize_t>pow(2, lvl)

    return n_nodes_at_lvl - bt_offset_from_left(v) - 1


cdef inline Py_ssize_t bt_left_leaf(Py_ssize_t v, Py_ssize_t height) nogil:
    """Determine the index of a nodes left most leaf"""
    cdef Py_ssize_t left_tip = <Py_ssize_t>pow(2, height) - 1
    cdef Py_ssize_t block_size

    if bt_is_root(v):
        return left_tip

    block_size = <Py_ssize_t>pow(2, height - floor(log2(v + 1)))

    return left_tip + (block_size * bt_offset_from_left(v))


cdef inline Py_ssize_t bt_right_leaf(Py_ssize_t v, Py_ssize_t height) nogil:
    """Determine the index of a nodes right most leaf"""
    cdef Py_ssize_t right_tip = <Py_ssize_t>pow(2, height + 1) - 2
    cdef Py_ssize_t block_size

    if bt_is_root(v):
        return right_tip

    block_size = <Py_ssize_t>pow(2, height - floor(log2(v + 1)))

    return right_tip - (block_size * bt_offset_from_right(v))
