# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from ._cutils import is_symmetric_and_hollow_cy
from ._cutils import distmat_reorder_cy, distmat_reorder_condensed_cy


def is_symmetric_and_hollow(mat):
    """Check if a Distance Matrix is symmetric and hollow.

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
    # is_symmetric_and_hollow_cy is optimized
    # for the common cas of c_contiguous.
    # For all other cases, make a copy.
    if not mat.flags.c_contiguous:
        mat = np.asarray(mat, order="C")

    return is_symmetric_and_hollow_cy(mat)


def is_symmetric(mat):
    """Check if a Distance Matrix is symmetric.

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
    # the is_hollow check is really cheap,
    # so can reuse is_symmetric_and_hollow
    return is_symmetric_and_hollow(mat)[0]


def is_hollow(mat):
    """Check if a Distance Matrix is hollow.

    Equivalent to np.trace(mat) == 0

    Parameters
    ----------
    mat : 2D array_like
        Distance matrix.

    Result:
    -------
    is_hollow: Boolean
        np.trace(mat) == 0

    """
    # is_symmetric_and_hollow_cy spends most
    # of its time in symetry check, just use numpy
    return np.trace(mat) == 0


def distmat_reorder_buf(in_mat, reorder_vec, out_mat, validate=False):
    """Reorder the rows and columns of a distance matrix given a reorder vector.

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

    Parameters
    ----------
    in_mat : 2D array_like
        Distance matrix
    reorder_vec : 1D_array_like
        List of permutation indexes
    out_mat : 2D array_like
        Output, Distance matrix,
        must be in c_order and same size as reorder_vec
    validate: boolean
        Optional, if True, validate reorder_vec content, detaults to False

    """
    np_reorder = np.asarray(reorder_vec, dtype=np.intp)
    if validate:
        maxsize = in_mat.shape[0]
        bad_cnt = np.where((np_reorder < 0) or (np_reorder >= maxsize))[0].size
        if bad_cnt > 0:
            raise ValueError("Invalid reorder_vec")

    if not in_mat.flags.c_contiguous:
        in_mat = np.asarray(in_mat, order="C")

    distmat_reorder_cy(in_mat, np_reorder, out_mat)


def distmat_reorder(in_mat, reorder_vec, validate=False):
    """Reorder the rows and columns of a distance matrix given a reorder vector.

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

    Parameters
    ----------
    in_mat : 2D array_like
        Distance matrix, must be in c_order
    reorder_vec : 1D_array_like
        List of permutation indexes
    validate: boolean
        Optional, if True, validate reorder_vec content, detaults to False

    Returns
    -------
    out_mat : 2D array_like
        Distance matrix

    """
    np_reorder = np.asarray(reorder_vec, dtype=np.intp)
    if validate:
        maxsize = in_mat.shape[0]
        bad_cnt = np.where((np_reorder < 0) or (np_reorder >= maxsize))[0].size
        if bad_cnt > 0:
            raise ValueError("Invalid reorder_vec")

    if not in_mat.flags.c_contiguous:
        in_mat = np.asarray(in_mat, order="C")

    out_mat = np.empty([np_reorder.size, np_reorder.size], in_mat.dtype)
    distmat_reorder_cy(in_mat, np_reorder, out_mat)
    return out_mat


def distmat_reorder_condensed(in_mat, reorder_vec, validate=False):
    """Reorder the rows and columns of a distance matrix given a reorder vector.

    Not all of the columns need to be used.

    For example:
     [ [0, 1, 2, 3] ,
       [1, 0, 4, 5] ,
       [2, 4, 0, 6] ,
       [3, 5, 6, 0] ]
     with
     [1,0,3,2]
     will result in
     [ 1, 5, 4 , 3, 2, 6 ]

    Parameters
    ----------
    in_mat : 2D array_like
        Distance matrix, must be in c_order
    reorder_vec : 1D_array_like
        List of permutation indexes
    validate: boolean
        Optional, if True, validate reorder_vec content, detaults to False

    Returns
    -------
    out_mat_condensed : 1D array_like
        Condensed distance matrix

    """
    np_reorder = np.asarray(reorder_vec, dtype=np.intp)
    if validate:
        maxsize = in_mat.shape[0]
        bad_cnt = np.where((np_reorder < 0) or (np_reorder >= maxsize))[0].size
        if bad_cnt > 0:
            raise ValueError("Invalid reorder_vec")

    if not in_mat.flags.c_contiguous:
        in_mat = np.asarray(in_mat, order="C")

    csize = ((np_reorder.size - 1) * np_reorder.size) // 2
    out_mat_condensed = np.empty([csize], in_mat.dtype)
    distmat_reorder_condensed_cy(in_mat, np_reorder, out_mat_condensed)
    return out_mat_condensed
