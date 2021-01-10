# -----------------------------------------------------------------------------
#  Copyright (c) 2013--, scikit-bio development team.
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
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

ctypedef fused TReal:
    float
    double

#
#    Compute E matrix from a distance matrix and store in temp centered matrix.
#    Squares and divides by -2 the input elementwise. Eq. 9.20 in
#    Legendre & Legendre 1998.
#    Compute sum of the rows at the same time.
#

# assuming mat is symmetric 

@cython.boundscheck(False)
@cython.wraparound(False)
def E_matrix_means(TReal[:, ::1] mat, TReal[:, ::1] centered, TReal[::1] row_means):
    cdef Py_ssize_t n_samples = mat.shape[0]
    cdef Py_ssize_t d2 = mat.shape[1]
    cdef Py_ssize_t d3 = centered.shape[1]
    cdef Py_ssize_t d4 = centered.shape[1]
    cdef Py_ssize_t d5 = row_means.shape[0]

    assert n_samples == d2
    assert n_samples == d3
    assert n_samples == d4
    assert n_samples == d5

    cdef Py_ssize_t row,col
    cdef TReal row_sum
    cdef TReal el0

    cdef TReal global_sum = 0.0
    for row in prange(n_samples, nogil=True):
        row_sum = 0.0

        for col in range(n_samples):
            el0 = mat[row,col]
            el0 =  -0.5*el0*el0
            centered[row,col] = el0
            # Note: do not use +=, so it is not flagged as a global reduction
            row_sum = row_sum + el0

        global_sum += row_sum
        row_means[row] = row_sum/n_samples

    return (global_sum/n_samples)/n_samples

#
#    Compute F matrix from E matrix.
#    Centring step: for each element, the mean of the corresponding
#    row and column are substracted, and the mean of the whole
#    matrix is added. Eq. 9.21 in Legendre & Legendre 1998.
#    Pseudo-code:
#    row_means = E_matrix.mean(axis=1, keepdims=True)
#    col_means = Transpose(row_means)
#    matrix_mean = E_matrix.mean()
#    E_matrix = E_matrix - row_means - col_means + matrix_mean
#


@cython.boundscheck(False)
@cython.wraparound(False)
def F_matrix_inplace(TReal[::1] row_means, TReal global_mean, TReal[:, ::1] centered):
    cdef Py_ssize_t n_samples = centered.shape[0]
    cdef Py_ssize_t d2 = centered.shape[1]
    cdef Py_ssize_t d3 = row_means.shape[0]

    assert n_samples == d2
    assert n_samples == d3

    cdef Py_ssize_t trow,tcol,row,col
    cdef Py_ssize_t trow_max,tcol_max
    cdef TReal gr_mean

    # use a tiled pattern to maximize locality of row_means
    for trow in prange(0, n_samples, 512, nogil=True):
        trow_max = min(trow+512, n_samples)

        for tcol in range(0, n_samples, 512):
            tcol_max = min(tcol+512, n_samples)

            for row in range(trow, trow_max, 1):
                gr_mean = global_mean - row_means[row]

                for col in range(tcol, tcol_max, 1):
                    # Note: do not use +=, so it is not flagged as a global reduction
                    centered[row,col] = centered[row,col] + (gr_mean - row_means[col])


#
# Centers a distance matrix, single precision
#
# mat       <const float *> Distance matrix to center, must be symmetric.
# centered  <      float *> Output centered matrix. Must be pre-allocated. Can point to mat.
#

@cython.boundscheck(False)
@cython.wraparound(False)
def center_distance_matrix(TReal[:, ::1] mat, TReal[:, ::1] centered):
    cdef Py_ssize_t n_samples = mat.shape[0]
    cdef Py_ssize_t d2 = mat.shape[1]
    cdef Py_ssize_t d3 = centered.shape[1]
    cdef Py_ssize_t d4 = centered.shape[1]

    assert n_samples == d2
    assert n_samples == d3
    assert n_samples == d4
  
    cdef TReal global_mean

    if TReal is float:
        dtype = np.float
    else:
        dtype = np.double

    row_means = np.zeros((n_samples,), dtype=dtype)

    global_mean = E_matrix_means(mat, centered, row_means)
    F_matrix_inplace(row_means, global_mean, centered)

