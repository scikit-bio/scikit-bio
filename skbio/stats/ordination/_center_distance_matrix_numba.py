# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Optional CPU Numba implementation of distance-matrix double-centering."""

import numpy as np

try:
    from numba import njit, prange

    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False


if NUMBA_AVAILABLE:

    @njit(parallel=True)
    def e_matrix_means_nb(mat, centered, row_means):
        """Apply E-matrix transform and collect row/global means in one pass."""
        n = mat.shape[0]
        global_sum = np.float64(0.0)

        for row in prange(n):
            row_sum = np.float64(0.0)

            for col in range(n):
                val = np.float64(mat[row, col])
                el = np.float64(-0.5) * val * val
                centered[row, col] = el
                row_sum += el

            row_means[row] = row_sum / np.float64(n)
            global_sum += row_sum

        return (global_sum / np.float64(n)) / np.float64(n)

    @njit(parallel=True)
    def f_matrix_inplace_nb(row_means, global_mean, centered):
        """Double-center E-matrix in-place."""
        n = centered.shape[0]

        n_blocks = (n + 23) // 24
        for brow in prange(n_blocks):
            trow = brow * 24
            trow_max = min(trow + 24, n)

            for tcol in range(0, n, 24):
                tcol_max = min(tcol + 24, n)

                for row in range(trow, trow_max):
                    gr_mean = global_mean - row_means[row]

                    for col in range(tcol, tcol_max):
                        centered[row, col] += gr_mean - row_means[col]

    def center_distance_matrix_nb(mat, centered):
        """Drop-in style replacement for ``center_distance_matrix_cy``.

        Parameters
        ----------
        mat : ndarray (n, n), float32 or float64, C-contiguous
            Input distance matrix.
        centered : ndarray (n, n), same dtype as ``mat``
            Pre-allocated output array. May alias ``mat`` for in-place use.
        """
        n = mat.shape[0]
        row_means = np.zeros(n, dtype=mat.dtype)
        global_mean = e_matrix_means_nb(mat, centered, row_means)
        f_matrix_inplace_nb(row_means, np.float64(global_mean), centered)
