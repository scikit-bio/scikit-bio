# -----------------------------------------------------------------------------
#  Copyright (c) 2013--, scikit-bio development team.
#
#  Distributed under the terms of the Modified BSD License.
#
#  The full license is in the file LICENSE.txt, distributed with this software.
# -----------------------------------------------------------------------------
import numpy as np
cimport cython
from cython.parallel import prange
from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

ctypedef fused TReal:
    float
    double

@cython.boundscheck(False)
@cython.wraparound(False)
def e_matrix_means_cy(TReal[:, ::1] mat, TReal[:, ::1] centered, TReal[::1] row_means):
    """
    Compute E matrix from a distance matrix, and 
    also compute the means in the process.

    Squares and divides by -2 the input elementwise. Eq. 9.20 in
    Legendre & Legendre 1998.


    Parameters
    ----------
    mat : 2D array_like
        Distance matrix.
    centered : 2D array_like
        Output, E matrix. Must be pre-allocated and same shape as mat.
        Can point to mat (i.e. in-place)
    row_means : 1D_array_like
        Output, Mean values of each row in `centered`
    Returns
    -------
    global_mean : real
        Global mean value
    """
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
    cdef long double row_sum
    cdef TReal el0

    cdef long double global_sum = 0.0
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

    cdef TReal global_mean = (global_sum/n_samples)/n_samples

    return global_mean

@cython.boundscheck(False)
@cython.wraparound(False)
def f_matrix_inplace_cy(TReal[::1] row_means, TReal global_mean, TReal[:, ::1] centered):
    """
    Compute F matrix from E matrix inplace.
    Centering step: for each element, the mean of the corresponding
    row and column are subtracted, and the mean of the whole
    matrix is added. Eq. 9.21 in Legendre & Legendre 1998.

    Modified from :func:`skbio.stats.ordination.f_matrix_inplace` function,

    Parameters
    ----------
    row_means : 1D_array_like
        Mean values of each row in `centered`
    global_mean : real
        Global mean value in `centered`
    centered : 2D array_like, must be symmetric
        In,  a matrix representing the "E matrix" as described above.
        Out, the centered matrix
    """
    cdef Py_ssize_t n_samples = centered.shape[0]
    cdef Py_ssize_t d2 = centered.shape[1]
    cdef Py_ssize_t d3 = row_means.shape[0]

    assert n_samples == d2
    assert n_samples == d3

    cdef Py_ssize_t trow,tcol,row,col
    cdef Py_ssize_t trow_max,tcol_max
    cdef TReal gr_mean

    # use a tiled pattern to maximize locality of row_means
    for trow in prange(0, n_samples, 24, nogil=True):
        trow_max = min(trow+24, n_samples)

        for tcol in range(0, n_samples, 24):
            tcol_max = min(tcol+24, n_samples)

            for row in range(trow, trow_max, 1):
                gr_mean = global_mean - row_means[row]

                for col in range(tcol, tcol_max, 1):
                    # Note: do not use +=, so it is not flagged as a global reduction
                    centered[row,col] = centered[row,col] + (gr_mean - row_means[col])

@cython.boundscheck(False)
@cython.wraparound(False)
def center_distance_matrix_cy(TReal[:, ::1] mat, TReal[:, ::1] centered):
    """
    Centers a distance matrix.

    Note: If the used distance was euclidean, pairwise distances
    needn't be computed from the data table Y because F_matrix =
    Y.dot(Y.T) (if Y has been centered).
    But since we're expecting distance_matrix to be non-euclidian,
    we do the following computation as per
    Numerical Ecology (Legendre & Legendre 1998).

    Parameters
    ----------
    mat : 2D array_like
        Distance matrix.
    centered : 2D array_like
        Output centered matrix. Must be pre-allocated and same shape as mat.
        Can point to mat (i.e. in-place)
    """
    cdef Py_ssize_t n_samples = mat.shape[0]
    cdef Py_ssize_t d2 = mat.shape[1]
    cdef Py_ssize_t d3 = centered.shape[1]
    cdef Py_ssize_t d4 = centered.shape[1]

    assert n_samples == d2
    assert n_samples == d3
    assert n_samples == d4
  
    cdef TReal global_mean

    if TReal is float:
        dtype_real = np.float32
    else:
        dtype_real = np.float64

    row_means_np = np.zeros((n_samples,), dtype=dtype_real)
    cdef TReal[::1] row_means = row_means_np

    global_mean = e_matrix_means_cy(mat, centered, row_means)
    f_matrix_inplace_cy(row_means, global_mean, centered)

