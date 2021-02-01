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
def is_symmetric_cy(TReal[:, ::1] mat):
    """
    Check is is symmetric
    Equivalent to not (mat.T != mat).any()

    Parameters
    ----------
    mat : 2D array_like
        Distance matrix.

    Result:
    -------
    is_symmetric: Boolean
        not (mat.T != mat).any()
    """
    cdef Py_ssize_t in_n = mat.shape[0]
    cdef Py_ssize_t in2 = mat.shape[1]

    assert in_n == in2

    cdef Py_ssize_t trow,tcol
    cdef Py_ssize_t trow_max,tcol_max
    cdef Py_ssize_t row,col

    # use int instead of bool for portabiltiy
    cdef int is_sym = True

    # use a tiled approach to maximize memory locality
    for trow in prange(0, in_n, 64, nogil=True):
        trow_max = min(trow+64, in_n)
        for tcol in range(0, in_n, 64):
            tcol_max = min(tcol+64, in_n)
            for row in range(trow, trow_max, 1):
                for col in range(tcol, tcol_max, 1):
                   is_sym &= mat[row,col]==mat[col,row]

    return (is_sym==True)
