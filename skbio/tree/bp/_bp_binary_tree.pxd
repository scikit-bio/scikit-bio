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


cdef inline Py_ssize_t bt_is_root(Py_ssize_t v) noexcept nogil:
    """Is v the root"""
    return v == 0


cdef inline Py_ssize_t bt_is_left_child(Py_ssize_t v) noexcept nogil:
    """Is v a left child of some node"""
    return 0 if bt_is_root(v) else v % 2


cdef inline Py_ssize_t bt_is_right_child(Py_ssize_t v) noexcept nogil:
    """Is v a right child of some node"""
    return 0 if bt_is_root(v) else 1 - (v % 2)


cdef inline Py_ssize_t bt_parent(Py_ssize_t v) noexcept nogil:
    """Get the index of the parent of v"""
    return 0 if bt_is_root(v) else (v - 1) // 2


cdef inline Py_ssize_t bt_left_child(Py_ssize_t v) noexcept nogil:
    """Get the index of the left child of v"""
    return 2 * v + 1


cdef inline Py_ssize_t bt_right_child(Py_ssize_t v) noexcept nogil:
    """Get the index of the right child of v"""
    return 2 * v + 2


cdef inline Py_ssize_t bt_left_sibling(Py_ssize_t v) noexcept nogil:
    """Get the index of the left sibling of v"""
    return v - 1


cdef inline Py_ssize_t bt_right_sibling(Py_ssize_t v) noexcept nogil:
    """Get the index of the right sibling of v"""
    return v + 1


cdef inline Py_ssize_t bt_is_leaf(Py_ssize_t v, Py_ssize_t height) noexcept nogil:
    """Determine if v is a leaf"""
    return <Py_ssize_t>(v >= pow(2, height) - 1)


cdef inline Py_ssize_t bt_node_from_left(Py_ssize_t pos, Py_ssize_t height) noexcept nogil:
    """Get the index from the left of a node at a given height"""
    return <Py_ssize_t>pow(2, height) - 1 + pos

