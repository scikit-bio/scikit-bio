r"""
Spatial Statistics (:mod:`skbio.stats.spatial`)
===============================================

.. currentmodule:: skbio.stats.spatial

This module provides functions for spatial analysis.

Functions
---------

.. autosummary::
   :toctree: generated/

   procrustes

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np


def procrustes(data1, data2):
    r"""Procrustes analysis, a similarity test for two data sets

    Each input matrix is a set of points or vectors (the rows of the matrix).
    The dimension of the space is the number of columns of each matrix. Given
    two identially sized matrices, procrustes standardizes both such that:

    - trace(AA') = 1  (A' is the transpose, and the product is a standard
      matrix product).
    - Both sets of points are centered around the origin.

    Procrustes ([1]_, [2]_) then applies the optimal transform to the second
    matrix (including scaling/dilation, rotations, and reflections) to minimize
    M^2 = sum(square(mtx1 - mtx2)), or the sum of the squares of the pointwise
    differences between the two input datasets.

    If two data sets have different dimensionality (different number of
    columns), simply add columns of zeros the the smaller of the two.

    This function was not designed to handle datasets with different numbers of
    datapoints (rows).

    Parameters
    ----------
    data1 : array_like
        matrix, n rows represent points in k (columns) space data1 is the
        reference data, after it is standardised, the data from data2 will
        be transformed to fit the pattern in data1 (must have >1 unique
        points).

    data2 : array_like
        n rows of data in k space to be fit to data1.  Must be the  same
        shape (numrows, numcols) as data1 (must have >1 unique points).


    Returns
    -------
    mtx1 : array_like
        a standardized version of data1
    mtx2 : array_like
        the orientation of data2 that best fits data1. Centered, but not
        necessarily trace(mtx2*mtx2') = 1
    disparity : array_like
        M^2 defined above


    Notes
    -----

    - The disparity should not depend on the order of the input matrices, but
      the output matrices will, as only the first output matrix is guaranteed
      to be scaled such that ``trace(AA') = 1``.

    - Duplicate datapoints are generally ok, duplicating a data point will
      increase its effect on the procrustes fit.

    - The disparity scales as the number of points per input matrix.

    References
    ----------

    .. [1] Krzanowski, W. J. (2000). "Principles of Multivariate analysis".
    .. [2] Gower, J. C. (1975). "Generalized procrustes analysis".

    Examples
    --------

    >>> import numpy as np
    >>> from skbio.stats.spatial import procrustes
    >>> a = np.array([[1, 3], [1, 2], [1, 1], [2, 1]], 'd')
    >>> b = np.array([[4, -2], [4, -4], [4, -6], [2, -6]], 'd')
    >>> mtx1, mtx2, disparity = procrustes(a, b)
    >>> print(round(disparity))
    0.0

    """
    num_rows, num_cols = np.shape(data1)
    if (num_rows, num_cols) != np.shape(data2):
        raise ValueError("input matrices must be of same shape")
    if num_rows == 0 or num_cols == 0:
        raise ValueError("input matrices must be >0 rows, >0 cols")

    # standardize each matrix
    mtx1 = _center(data1)
    mtx2 = _center(data2)

    if (not np.any(mtx1)) or (not np.any(mtx2)):
        raise ValueError("input matrices must contain >1 unique points")

    mtx1 = _normalize(mtx1)
    mtx2 = _normalize(mtx2)

    # transform mtx2 to minimize disparity (sum( (mtx1[i,j] - mtx2[i,j])^2) )
    mtx2 = _match_points(mtx1, mtx2)

    disparity = _get_disparity(mtx1, mtx2)

    return mtx1, mtx2, disparity


def _center(mtx):
    """Translate all data (rows of the matrix) to center on the origin

    Parameters
    ----------
    mtx : array_like
        Matrix to translate the data for.

    Returns
    -------
    result : array_like ('d') array
        Shifted version of the input data.  The new matrix is such that the
        center of mass of the row vectors is centered at the origin.

    """
    result = np.array(mtx, 'd')
    result -= np.mean(result, 0)
    # subtract each column's mean from each element in that column
    return result


def _normalize(mtx):
    """change scaling of data (in rows) such that trace(mtx*mtx') = 1

    Parameters
    ----------
    mtx : array_like
        Matrix to scale the data for.

    Notes
    -----
    mtx' denotes the transpose of mtx

    """
    mtx = np.asarray(mtx, dtype=float)
    return mtx / np.linalg.norm(mtx)


def _match_points(mtx1, mtx2):
    """Returns a transformed mtx2 that matches mtx1.

    Returns
    -------

    A new matrix which is a transform of mtx2.  Scales and rotates a copy of
    mtx 2.  See procrustes docs for details.

    """
    u, s, vh = np.linalg.svd(np.dot(np.transpose(mtx1), mtx2))
    q = np.dot(np.transpose(vh), np.transpose(u))
    new_mtx2 = np.dot(mtx2, q)
    new_mtx2 *= np.sum(s)

    return new_mtx2


def _get_disparity(mtx1, mtx2):
    """Measures the dissimilarity between two data sets

    Returns
    -------

    M^2 = sum(square(mtx1 - mtx2)), the pointwise sum of squared differences

    """
    return(np.sum(np.square(mtx1 - mtx2)))
