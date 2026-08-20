# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

"""Optional CPU Numba implementation of distance-matrix double-centering.

Mirrors the Cython implementation in ``_cutils.pyx`` (``e_matrix_means_cy``,
``f_matrix_inplace_cy`` and ``center_distance_matrix_cy``) so the two engines
are numerically interchangeable.
"""

import numpy as np

try:
    from numba import njit, prange

    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False


if NUMBA_AVAILABLE:

    @njit(parallel=True)
    def e_matrix_means_nb(mat, centered, row_means):
        """Apply the E-matrix transform and collect row/global means in one pass.

        Writes ``-0.5 * d**2`` into ``centered`` and returns the global mean of
        that transformed matrix. Values are promoted to float64 before squaring
        to match the Cython ``double`` accumulator (matters for float32 input).
        """
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
        """Double-center the E-matrix in-place.

        Uses the same 24-wide block tiling as ``f_matrix_inplace_cy`` so the
        two engines share a memory-access pattern.
        """
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
        """Drop-in replacement for ``center_distance_matrix_cy``.

        Parameters
        ----------
        mat : ndarray of shape (n, n)
            Input distance matrix, C-contiguous, float32 or float64.
        centered : ndarray of shape (n, n)
            Pre-allocated output array of the same dtype as ``mat``. May alias
            ``mat`` for in-place centering.

        """
        n = mat.shape[0]
        row_means = np.zeros(n, dtype=mat.dtype)
        global_mean = e_matrix_means_nb(mat, centered, row_means)
        f_matrix_inplace_nb(row_means, np.float64(global_mean), centered)
