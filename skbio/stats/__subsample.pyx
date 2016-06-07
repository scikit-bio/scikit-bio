# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
cimport numpy as cnp


def _subsample_counts_without_replacement(
    cnp.ndarray[cnp.int64_t, ndim=1] counts, n, counts_sum):
    cdef:
        cnp.ndarray[cnp.int64_t, ndim=1] result, permuted, unpacked
        cnp.int64_t cnt
        Py_ssize_t unpacked_idx, i, j

    unpacked = np.empty(counts_sum, dtype=int)
    unpacked_idx = 0
    for i in range(counts.shape[0]):
        cnt = counts[i]
        for j in range(cnt):
            unpacked[unpacked_idx] = i
            unpacked_idx += 1

    permuted = np.random.permutation(unpacked)[:n]

    result = np.zeros_like(counts)
    for idx in range(permuted.shape[0]):
        result[permuted[idx]] += 1

    return result
