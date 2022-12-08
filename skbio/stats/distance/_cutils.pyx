# -----------------------------------------------------------------------------
#  Copyright (c) 2021-2021, scikit-bio development team.
#
#  Distributed under the terms of the Modified BSD License.
#
#  The full license is in the file LICENSE.txt, distributed with this software.
# -----------------------------------------------------------------------------

import numpy as np
cimport cython
from cython.parallel import prange

ctypedef fused TReal:
    float
    double

@cython.boundscheck(False)
@cython.wraparound(False)
def is_symmetric_and_hollow_cy(TReal[:, ::1] mat):
    """
    Check if mat is symmetric and hollow.
    Equivalent to [not (mat.T != mat).any(), np.trace(mat) == 0]

    Parameters
    ----------
    mat : 2D array_like
        Distance matrix.

    Result:
    -------
    is_symmetric: Boolean
        not (mat.T != mat).any()
    is_hollow: Boolean
        np.trace(mat) == 0
    """
    cdef Py_ssize_t in_n = mat.shape[0]
    cdef Py_ssize_t in2 = mat.shape[1]

    assert in_n == in2

    cdef Py_ssize_t trow,tcol
    cdef Py_ssize_t trow_max,tcol_max
    cdef Py_ssize_t row,col

    cdef TReal testval

    # use int instead of bool for portability
    cdef int is_sym = True
    cdef int is_hollow = True

    # use a tiled approach to maximize memory locality
    for trow in prange(0, in_n, 24, nogil=True):
        trow_max = min(trow+24, in_n)
        for tcol in range(0, in_n, 24):
            tcol_max = min(tcol+24, in_n)
            for row in range(trow, trow_max, 1):
                for col in range(tcol, tcol_max, 1):
                   is_sym &= (mat[row,col]==mat[col,row])
            if (trow==tcol): # diagonal block, only ones that can have col==row
                for col in range(tcol, tcol_max, 1):
                   is_hollow &= (mat[col,col]==0)

    return [(is_sym==True), (is_hollow==True)]


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
           out_mat[row,col] = in_mat[vrow, reorder_vec[col]]

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
     [ 1, 5, 4, 3, 2, 6 ]

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
        idx = row*(out_n-1) - ((row-1)*row)/2
        for col in range(out_n-row-1):
           out_mat_condensed[idx+col] = in_mat[vrow, reorder_vec[col+row+1]]


@cython.boundscheck(False)
@cython.wraparound(False)
def mantel_perm_pearsonr_cy(TReal[:, ::1] x_data, long[:, ::1] perm_order,
                            TReal xmean, TReal normxm,
                            TReal[::1] ym_normalized,
                            TReal[::1] permuted_stats):
    """
    Fused permute, fma, pearsonr for mantel.

    Replaces the following python code:
    def _mantel_perm_pearsonr_one(x_flat, xmean, normxm, ym_normalized):
        # inline pearsonr, condensed from scipy.stats.pearsonr
        # and reusing some of the known values
        xm_normalized = (x_flat - xmean)/normxm
        one_stat = np.dot(xm_normalized, ym_normalized)
        one_stat = max(min(one_stat, 1.0), -1.0)
        return one_stat

    perm_gen = (_mantel_perm_pearsonr_one(distmat_reorder_condensed(x._data, perm_order[p,:]),
                                          xmean, normxm, ym_normalized)
                for p in range(permutations))
    permuted_stats = np.fromiter(perm_gen, np.float, count=permutations)

    Parameters
    ----------
    x_data : 2D array_like
        Distance matrix.
    perm_order : 2D array_like
        List of permutation orders.
    xmean: real
        Mean value of condensed x_data
    normxm: real
        Norm of pre-processed xm
    ym_normalized : 1D_array_like
        Normalized condensed y_data
    permuted_stats : 1D array_like
        Output, Pearson stats
    """
    cdef Py_ssize_t in_n = x_data.shape[0]
    cdef Py_ssize_t in2 = x_data.shape[1]
    cdef Py_ssize_t perms_n = perm_order.shape[0]
    cdef Py_ssize_t out_n = perm_order.shape[1]
    cdef Py_ssize_t y_n = ym_normalized.shape[0]
    cdef Py_ssize_t on2 = permuted_stats.shape[0]

    assert in_n == in2
    assert y_n == ((out_n-1)*out_n)/2
    assert perms_n == on2

    cdef Py_ssize_t p
    cdef Py_ssize_t row,col,icol
    cdef Py_ssize_t vrow
    cdef Py_ssize_t idx

    cdef TReal mul = 1.0/normxm
    cdef TReal add = -xmean/normxm

    cdef TReal my_ps
    cdef TReal yval
    cdef TReal xval

    for p in prange(perms_n, nogil=True):
        my_ps = 0.0
        for row in range(out_n-1):
            vrow = perm_order[p, row]
            idx = row*(out_n-1) - ((row-1)*row)/2
            for icol in range(out_n-row-1):
               col = icol+row+1
               yval = ym_normalized[idx+icol]
               xval = x_data[vrow, perm_order[p, col]]*mul + add
               # do not use += to avoid having prange consider it for reduction
               my_ps = yval*xval + my_ps

        # Presumably, if abs(one_stat) > 1, then it is only some small artifact of
        # floating point arithmetic.
        if my_ps>1.0:
            my_ps = 1.0
        elif my_ps<-1.0:
            my_ps = -1.0
        permuted_stats[p] = my_ps


@cython.boundscheck(False)
@cython.wraparound(False)
def permanova_f_stat_sW_cy(TReal[:, ::1] distance_matrix,
                           Py_ssize_t[::1] group_sizes,
                           Py_ssize_t[::1] grouping):
    """Compute PERMANOVA pseudo-F partial statistic."""
    cdef Py_ssize_t in_n = distance_matrix.shape[0]
    cdef Py_ssize_t in2 = distance_matrix.shape[1]
    cdef Py_ssize_t in3 = grouping.shape[0]

    assert in_n == in2
    assert in_n == in3

    cdef double s_W = 0.0

    cdef Py_ssize_t group_idx
    cdef double local_s_W
    cdef double val

    cdef Py_ssize_t row, col, rowi, coli

    for rowi in prange(in_n/2, nogil=True):
        # since columns get shorter, combine first and last
        row=rowi
        local_s_W = 0.0
        group_idx = grouping[row]
        for coli in range(in_n-row-1):
            col = coli+row+1
            if grouping[col] == group_idx:
                val = distance_matrix[row,col]
                local_s_W = local_s_W + val * val
        s_W += local_s_W/group_sizes[group_idx]

        row = in_n-rowi-2
        if row!=rowi: # don't double count
            local_s_W = 0.0
            group_idx = grouping[row]
            for coli in range(in_n-row-1):
                col = coli+row+1
                if grouping[col] == group_idx:
                    val = distance_matrix[row,col]
                    local_s_W = local_s_W + val * val
            s_W += local_s_W/group_sizes[group_idx]

    return s_W
