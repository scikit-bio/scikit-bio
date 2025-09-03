# -----------------------------------------------------------------------------
#  Copyright (c) 2021-2025, scikit-bio development team.
#
#  Distributed under the terms of the Modified BSD License.
#
#  The full license is in the file LICENSE.txt, distributed with this software.
# -----------------------------------------------------------------------------

# Define NPY_NO_DEPRECATED_API to avoid deprecation warnings
cdef extern from *:
    """
    #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
    """

import numpy as np
cimport numpy as cnp
cnp.import_array()
cimport cython
from cython.parallel import prange
from libc.math cimport sqrt, fabs

ctypedef cnp.npy_intp intp_t
ctypedef cnp.float64_t float64_t
ctypedef cnp.float32_t float32_t

ctypedef fused TReal:
    float
    double

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def permanova_f_stat_sW_opt_cy(TReal[:, ::1] distance_matrix,
                               Py_ssize_t[::1] group_sizes,
                               Py_ssize_t[::1] grouping):
    """Optimized PERMANOVA pseudo-F partial statistic computation.
    
    This is an optimized version that:
    1. Uses better cache locality patterns through blocking
    2. Reduces redundant computations
    3. Minimizes memory access patterns
    4. Processes only upper triangle of symmetric matrix
    """
    cdef Py_ssize_t n = distance_matrix.shape[0]
    cdef Py_ssize_t n_groups = group_sizes.shape[0]
    
    # Pre-allocate group-wise accumulators
    cdef cnp.ndarray[double, ndim=1] group_sums_arr = np.zeros(n_groups, dtype=np.float64)
    cdef double[:] group_sums = group_sums_arr
    
    cdef double s_W = 0.0
    cdef Py_ssize_t i, j, gi, gj, g
    cdef double val, val_sq
    cdef Py_ssize_t block_size = 32  # Cache-friendly block size
    cdef Py_ssize_t bi, bj, block_i_end, block_j_end
    
    # Process in blocks for better cache locality
    for bi in range(0, n, block_size):
        block_i_end = min(bi + block_size, n)
        
        for bj in range(bi, n, block_size):
            block_j_end = min(bj + block_size, n)
            
            # Process block
            for i in range(bi, block_i_end):
                gi = grouping[i]
                
                # Only process upper triangle
                for j in range(max(i + 1, bj), block_j_end):
                    gj = grouping[j]
                    
                    # Check if same group
                    if gi == gj:
                        val = distance_matrix[i, j]
                        val_sq = val * val
                        group_sums[gi] += val_sq
    
    # Final reduction with group size normalization
    for g in range(n_groups):
        if group_sizes[g] > 0:
            s_W += group_sums[g] / <double>group_sizes[g]
    
    return s_W


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef double permanova_f_stat_sW_inline(TReal[:, ::1] distance_matrix,
                                       Py_ssize_t[::1] group_sizes,
                                       Py_ssize_t[::1] grouping) nogil:
    """Inline version for nogil contexts."""
    cdef Py_ssize_t n = distance_matrix.shape[0]
    cdef Py_ssize_t n_groups = group_sizes.shape[0]
    
    cdef double s_W = 0.0
    cdef double group_sum
    cdef Py_ssize_t i, j, gi, gj, g
    cdef double val
    
    # For each group, compute within-group sum of squared distances
    for g in range(n_groups):
        if group_sizes[g] == 0:
            continue
            
        group_sum = 0.0
        
        # Only process upper triangle
        for i in range(n - 1):
            if grouping[i] != g:
                continue
            
            for j in range(i + 1, n):
                if grouping[j] == g:
                    val = distance_matrix[i, j]
                    group_sum += val * val
        
        s_W += group_sum / <double>group_sizes[g]
    
    return s_W


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def permanova_compute_f_stat_opt_cy(Py_ssize_t sample_size, 
                                    Py_ssize_t num_groups,
                                    TReal[:, ::1] distance_matrix,
                                    Py_ssize_t[::1] group_sizes,
                                    double s_T,
                                    Py_ssize_t[::1] grouping):
    """Compute complete PERMANOVA F-statistic in Cython.
    
    This combines the s_W computation with F-stat calculation
    to avoid Python overhead.
    """
    cdef double s_W = permanova_f_stat_sW_opt_cy(distance_matrix, group_sizes, grouping)
    cdef double s_A = s_T - s_W
    
    # Compute F-statistic with division safety
    cdef double numerator = s_A / <double>(num_groups - 1)
    cdef double denominator = s_W / <double>(sample_size - num_groups)
    
    if denominator == 0:
        return 0.0
    
    return numerator / denominator


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def permanova_batch_f_stats_opt_cy(Py_ssize_t sample_size,
                                   Py_ssize_t num_groups, 
                                   TReal[:, ::1] distance_matrix,
                                   Py_ssize_t[::1] group_sizes,
                                   double s_T,
                                   Py_ssize_t[:, ::1] groupings_batch):
    """Compute F-statistics for multiple permutations in batch.
    
    This processes multiple permutations to improve efficiency.
    """
    cdef Py_ssize_t n_perms = groupings_batch.shape[0]
    cdef Py_ssize_t n = distance_matrix.shape[0]
    
    # Output array for F-statistics
    cdef cnp.ndarray[double, ndim=1] f_stats = np.empty(n_perms, dtype=np.float64)
    cdef double[:] f_stats_view = f_stats
    
    cdef Py_ssize_t p
    cdef double s_W, s_A, f_stat
    
    # Process permutations (can be parallelized with OpenMP)
    for p in prange(n_perms, nogil=True):
        # Compute s_W for this permutation
        s_W = permanova_f_stat_sW_inline(
            distance_matrix, 
            group_sizes,
            groupings_batch[p, :]
        )
        
        s_A = s_T - s_W
        
        # Compute F-statistic
        if (sample_size - num_groups) > 0 and s_W > 0:
            f_stat = (s_A / <double>(num_groups - 1)) / (s_W / <double>(sample_size - num_groups))
        else:
            f_stat = 0.0
        
        f_stats_view[p] = f_stat
    
    return f_stats


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def compute_s_T_opt_cy(TReal[:, ::1] distance_matrix):
    """Optimized computation of s_T (total sum of squared distances)."""
    cdef Py_ssize_t n = distance_matrix.shape[0]
    cdef double s_T = 0.0
    cdef Py_ssize_t i, j
    cdef double val
    
    # Only process upper triangle (matrix is symmetric)
    for i in prange(n - 1, nogil=True):
        for j in range(i + 1, n):
            val = distance_matrix[i, j]
            s_T += val * val
    
    # Normalize by sample size
    s_T = s_T / <double>n
    
    return s_T