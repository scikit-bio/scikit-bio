# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.util._decorator import experimental


@experimental(as_of="0.4.0")
def mean_and_std(a, axis=None, weights=None, with_mean=True, with_std=True,
                 ddof=0):
    """Compute the weighted average and standard deviation along the
    specified axis.

    Parameters
    ----------
    a : array_like
        Calculate average and standard deviation of these values.
    axis : int, optional
        Axis along which the statistics are computed. The default is
        to compute them on the flattened array.
    weights : array_like, optional
        An array of weights associated with the values in `a`. Each
        value in `a` contributes to the average according to its
        associated weight. The weights array can either be 1-D (in
        which case its length must be the size of `a` along the given
        axis) or of the same shape as `a`. If `weights=None`, then all
        data in `a` are assumed to have a weight equal to one.
    with_mean : bool, optional, defaults to True
        Compute average if True.
    with_std : bool, optional, defaults to True
        Compute standard deviation if True.
    ddof : int, optional, defaults to 0
        It means delta degrees of freedom. Variance is calculated by
        dividing by `n - ddof` (where `n` is the number of
        elements). By default it computes the maximum likelyhood
        estimator.
    Returns
    -------
    average, std
        Return the average and standard deviation along the specified
        axis. If any of them was not required, returns `None` instead
    """
    if not (with_mean or with_std):
        raise ValueError("Either the mean or standard deviation need to be"
                         " computed.")
    a = np.asarray(a)
    if weights is None:
        avg = a.mean(axis=axis) if with_mean else None
        std = a.std(axis=axis, ddof=ddof) if with_std else None
    else:
        avg = np.average(a, axis=axis, weights=weights)
        if with_std:
            if axis is None:
                variance = np.average((a - avg)**2, weights=weights)
            else:
                # Make sure that the subtraction to compute variance works for
                # multidimensional arrays
                a_rolled = np.rollaxis(a, axis)
                # Numpy doesn't have a weighted std implementation, but this is
                # stable and fast
                variance = np.average((a_rolled - avg)**2, axis=0,
                                      weights=weights)
            if ddof != 0:  # Don't waste time if variance doesn't need scaling
                if axis is None:
                    variance *= a.size / (a.size - ddof)
                else:
                    variance *= a.shape[axis] / (a.shape[axis] - ddof)
            std = np.sqrt(variance)
        else:
            std = None
        avg = avg if with_mean else None
    return avg, std


@experimental(as_of="0.4.0")
def scale(a, weights=None, with_mean=True, with_std=True, ddof=0, copy=True):
    """Scale array by columns to have weighted average 0 and standard
    deviation 1.

    Parameters
    ----------
    a : array_like
        2D array whose columns are standardized according to the
        weights.
    weights : array_like, optional
        Array of weights associated with the columns of `a`. By
        default, the scaling is unweighted.
    with_mean : bool, optional, defaults to True
        Center columns to have 0 weighted mean.
    with_std : bool, optional, defaults to True
        Scale columns to have unit weighted std.
    ddof : int, optional, defaults to 0
        If with_std is True, variance is calculated by dividing by `n
        - ddof` (where `n` is the number of elements). By default it
        computes the maximum likelyhood stimator.
    copy : bool, optional, defaults to True
        Whether to perform the standardization in place, or return a
        new copy of `a`.

    Returns
    -------
    2D ndarray
        Scaled array.

    Notes
    -----
    Wherever std equals 0, it is replaced by 1 in order to avoid
    division by zero.
    """
    if copy:
        a = a.copy()
    a = np.asarray(a, dtype=np.float64)
    avg, std = mean_and_std(a, axis=0, weights=weights, with_mean=with_mean,
                            with_std=with_std, ddof=ddof)
    if with_mean:
        a -= avg
    if with_std:
        std[std == 0] = 1.0
        a /= std
    return a


@experimental(as_of="0.4.0")
def svd_rank(M_shape, S, tol=None):
    """Matrix rank of `M` given its singular values `S`.

    See `np.linalg.matrix_rank` for a rationale on the tolerance
    (we're not using that function because it doesn't let us reuse a
    precomputed SVD)."""
    if tol is None:
        tol = S.max() * max(M_shape) * np.finfo(S.dtype).eps
    return np.sum(S > tol)


@experimental(as_of="0.4.0")
def corr(x, y=None):
    """Computes correlation between columns of `x`, or `x` and `y`.

    Correlation is covariance of (columnwise) standardized matrices,
    so each matrix is first centered and scaled to have variance one,
    and then their covariance is computed.

    Parameters
    ----------
    x : 2D array_like
        Matrix of shape (n, p). Correlation between its columns will
        be computed.
    y : 2D array_like, optional
        Matrix of shape (n, q). If provided, the correlation is
        computed between the columns of `x` and the columns of
        `y`. Else, it's computed between the columns of `x`.

    Returns
    -------
    correlation
        Matrix of computed correlations. Has shape (p, p) if `y` is
        not provided, else has shape (p, q).
    """
    x = np.asarray(x)
    if y is not None:
        y = np.asarray(y)
        if y.shape[0] != x.shape[0]:
            raise ValueError("Both matrices must have the same number of rows")
        x, y = scale(x), scale(y)
    else:
        x = scale(x)
        y = x
    # Notice that scaling was performed with ddof=0 (dividing by n,
    # the default), so now we need to remove it by also using ddof=0
    # (dividing by n)
    return x.T.dot(y) / x.shape[0]


@experimental(as_of="0.4.0")
def e_matrix(distance_matrix):
    """Compute E matrix from a distance matrix.

    Squares and divides by -2 the input elementwise. Eq. 9.20 in
    Legendre & Legendre 1998."""
    return distance_matrix * distance_matrix / -2


def f_matrix(E_matrix):
    """Compute F matrix from E matrix.

    Centring step: for each element, the mean of the corresponding
    row and column are substracted, and the mean of the whole
    matrix is added. Eq. 9.21 in Legendre & Legendre 1998."""
    row_means = E_matrix.mean(axis=1, keepdims=True)
    col_means = E_matrix.mean(axis=0, keepdims=True)
    matrix_mean = E_matrix.mean()
    return E_matrix - row_means - col_means + matrix_mean


def center_distance_matrix(distance_matrix, inplace=False):
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
    distance_matrix : 2D array_like
        Distance matrix.
    inplace : bool, optional
        Whether or not to center the given distance matrix in-place, which
        is more efficient in terms of memory and computation.
    """
    if inplace:
        return _f_matrix_inplace(_e_matrix_inplace(distance_matrix))
    else:
        return f_matrix(e_matrix(distance_matrix))


def _e_matrix_inplace(distance_matrix):
    """
    Compute E matrix from a distance matrix inplace.
    Squares and divides by -2 the input element-wise. Eq. 9.20 in
    Legendre & Legendre 1998.

    Modified from :func:`skbio.stats.ordination.e_matrix` function,
    performing row-wise operations to avoid excessive memory allocations.

    Parameters
    ----------
    distance_matrix : 2D array_like
        Distance matrix.
    """
    distance_matrix = distance_matrix.astype(np.float)

    for i in np.arange(len(distance_matrix)):
        distance_matrix[i] = (distance_matrix[i] * distance_matrix[i]) / -2
    return distance_matrix


def _f_matrix_inplace(e_matrix):
    """
    Compute F matrix from E matrix inplace.
    Centering step: for each element, the mean of the corresponding
    row and column are subtracted, and the mean of the whole
    matrix is added. Eq. 9.21 in Legendre & Legendre 1998.

    Modified from :func:`skbio.stats.ordination.f_matrix` function,
    performing row-wise operations to avoid excessive memory allocations.

    Parameters
    ----------
    e_matrix : 2D array_like
        A matrix representing the "E matrix" as described above.
    """
    e_matrix = e_matrix.astype(np.float)

    row_means = np.zeros(len(e_matrix), dtype=float)
    col_means = np.zeros(len(e_matrix), dtype=float)
    matrix_mean = 0.0

    for i in np.arange(len(e_matrix)):
        row_means[i] = e_matrix[i].mean()
        matrix_mean += e_matrix[i].sum()
        col_means += e_matrix[i]
    matrix_mean /= len(e_matrix) ** 2
    col_means /= len(e_matrix)

    for i in np.arange(len(e_matrix)):
        v = e_matrix[i]
        v -= row_means[i]
        v -= col_means
        v += matrix_mean
        e_matrix[i] = v
    return e_matrix
