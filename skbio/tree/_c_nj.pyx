# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

cimport cython
from libc.float cimport DBL_MAX


@cython.boundscheck(False)
@cython.wraparound(False)
def nj_minq_cy(double[:, :] dm, double[:] sums):
    r"""Find the minimum value in a Q-matrix during the NJ algorithm.

    Parameters
    ----------
    dm : (N, N) ndarray
        Distance matrix.
    sums : (N,) ndarray
        Distance sum vector.

    Returns
    -------
    (int, int)
        Row and column indices of the minimum value.

    See Also
    --------
    skbio.tree.nj

    Notes
    -----
    The agglomerative clustering algorithm neighbor joining (NJ) involves converting
    the original distance matrix into a "Q-matrix" and finding the location of the
    minimum _q_ value in the matrix.

    .. math::

        Q(i, j) = (n - 2) d(i, j) - \sum d(i) - \sum d(j)

    This function calculates _q_ values and finds the minimum as the calculation goes,
    and avoids creating the entire Q-matrix, thereby saving compute.

    """
    cdef Py_ssize_t n = dm.shape[0]
    cdef Py_ssize_t i, j
    cdef int n_2 = n - 2

    # current minimum q value and its location
    cdef Py_ssize_t min_i, min_j
    cdef double min_q = DBL_MAX

    # current q value and minimum q value plus \sum d(i)
    cdef double q_plus, min_q_plus

    # loop the upper-right triangle of the distance matrix
    # i < j is guaranteed
    for i in range(n):
        min_q_plus = min_q + sums[i]
        for j in range(i + 1, n):
            q_plus = dm[i, j] * n_2 - sums[j]
            if q_plus < min_q_plus:
                min_q_plus = q_plus
                min_q = q_plus - sums[i]
                min_i, min_j = i, j

    return min_i, min_j
