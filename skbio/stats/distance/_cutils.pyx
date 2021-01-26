# -----------------------------------------------------------------------------
#  Copyright (c) 2021-2021, scikit-bio development team.
#
#  Distributed under the terms of the Modified BSD License.
#
#  The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

# distutils: extra_compile_args=-fopenmp
# distutils: extra_link_args=-fopenmp

import numpy as np
cimport cython
from cython.parallel import prange

ctypedef fused TReal:
    float
    double

@cython.boundscheck(False)
@cython.wraparound(False)
def distmat_reorder_cy(TReal[:, ::1] in_mat, long[::1] reorder_vec, 
                       TReal[:, ::1] out_mat):
    """
    Reorder the rows and columns of a distance matrix
    given a reorder vector.
    Not all of the columns need to be used.

    For example:
     [ [0, 1, 2, 3] ,
       [1, 0, 4, 5] ,
       [2, 4, 0, 6] ,
       [3, 5, 6, 0] ]
     with
     [1,0,3,2]
     will result in
     [ [0, 1, 5, 4] ,
       [1, 0, 3, 2] ,
       [5, 3, 0, 6] ,
       [4, 2, 6, 0] ]

    Note: No error checking is performed.
          The caller must ensure that all values in reorder_vec are valid

    Parameters
    ----------
    in_mat : 2D array_like
        Distance matrix.
    reorder_vec : 1D_array_like
        List of permutation indexes
    out_mat : 2D array_like
        Output, Distance matrix, must be same size as reorder_vec
    """
    cdef Py_ssize_t in_n = in_mat.shape[0]
    cdef Py_ssize_t in2 = in_mat.shape[1]
    cdef Py_ssize_t out_n = reorder_vec.shape[0]
    cdef Py_ssize_t on2 = out_mat.shape[0]
    cdef Py_ssize_t on3 = out_mat.shape[1]

    assert in_n == in2
    assert out_n == on2
    assert out_n == on3

    cdef Py_ssize_t row,col
    cdef Py_ssize_t vrow

    for row in prange(out_n, nogil=True):
        vrow = reorder_vec[row]
        for col in range(out_n):
           out_mat[row,col] = in_mat[vrow,reorder_vec[col]]

@cython.boundscheck(False)
@cython.wraparound(False)
def distmat_reorder_condensed_cy(TReal[:, ::1] in_mat, long[::1] reorder_vec,
                                  TReal[::1] out_mat_condensed):
    """
    Reorder the rows and columns of a distance matrix
    given a reorder vector.
    Not all of the columns need to be used.

    For example:
     [ [0, 1, 2, 3] ,
       [1, 0, 4, 5] ,
       [2, 4, 0, 6] ,
       [3, 5, 6, 0] ]
     with
     [1,0,3,2]
     will result in
     [ 1, 5, 4 , 3, 2, 6 ]

    Note: No error checking is performed.
          The caller must ensure that all values in reorder_vec are valid

    Parameters
    ----------
    in_mat : 2D array_like
        Distance matrix.
    reorder_vec : 1D_array_like
        List of permutation indexes
    out_mat_condensed : 1D array_like
        Output, condensed distance matrix
    """
    cdef Py_ssize_t in_n = in_mat.shape[0]
    cdef Py_ssize_t in2 = in_mat.shape[1]
    cdef Py_ssize_t out_n = reorder_vec.shape[0]
    cdef Py_ssize_t on2 = out_mat_condensed.shape[0]

    assert in_n == in2
    assert on2 == ((out_n-1)*out_n)/2

    cdef Py_ssize_t row,col
    cdef Py_ssize_t vrow
    cdef Py_ssize_t idx

    for row in prange(out_n-1, nogil=True):
        vrow = reorder_vec[row]
        idx = row*out_n - ((row-1)*row)/2
        for col in range(out_n-row-1):
           out_mat_condensed[idx+col] = in_mat[vrow,reorder_vec[col+row+1]]


