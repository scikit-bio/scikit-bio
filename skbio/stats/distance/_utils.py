# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import array_api_compat as _aac

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

    Notes
    -----
    A non-NumPy array-API buffer (e.g. GPU-resident) is checked via its own array
    namespace; a NumPy array uses the parallel Cython kernel.

    This function uses parallel computation for improved performance.
    See the :install:`parallelization guide <#parallelization>` for information on
    controlling the number of threads used.

    """
    if _aac.is_array_api_obj(mat) and not _aac.is_numpy_array(mat):
        # Non-NumPy (e.g. GPU-resident) buffer: check symmetry and hollowness via
        # the array-API namespace. A NaN makes the inequality True, so this also
        # rejects NaNs, matching the Cython path.
        xp = _aac.array_namespace(mat)
        is_symmetric = not bool(xp.any(xp.matrix_transpose(mat) != mat))
        is_hollow = not bool(xp.any(xp.linalg.diagonal(mat) != 0))
        return is_symmetric, is_hollow

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

    Notes
    -----
    This function uses parallel computation for improved performance.
    See the :install:`parallelization guide <#parallelization>` for information on
    controlling the number of threads used.

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
    # of its time in symmetry check, just use numpy
    return np.trace(mat) == 0


def _validate_order(order, mat):
    if np.any((order < 0) | (order >= mat.shape[0])):
        raise ValueError("Invalid reorder_vec")


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
        Optional, if True, validate reorder_vec content, defaults to False

    Returns
    -------
    out_mat : 2D array_like
        Distance matrix

    Notes
    -----
    A non-NumPy array-API buffer (e.g. GPU-resident) is reordered on its own device
    via the array namespace; a NumPy array uses the parallel Cython kernel.

    This function uses parallel computation for improved performance.
    See the :install:`parallelization guide <#parallelization>` for information on
    controlling the number of threads used.

    """
    if _aac.is_array_api_obj(in_mat) and not _aac.is_numpy_array(in_mat):
        # Non-NumPy (e.g. GPU-resident) buffer: reorder rows and columns with an
        # xp fancy-index gather, keeping the result on the input's device.
        if validate:
            _validate_order(np.asarray(reorder_vec, dtype=np.intp), in_mat)
        xp = _aac.array_namespace(in_mat)
        order = xp.asarray(
            np.asarray(reorder_vec, dtype=int),
            device=getattr(in_mat, "device", None),
        )
        return in_mat[order][:, order]

    np_reorder = np.asarray(reorder_vec, dtype=np.intp)
    if validate:
        _validate_order(np_reorder, in_mat)

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
        Optional, if True, validate reorder_vec content, defaults to False

    Returns
    -------
    out_mat_condensed : 1D array_like
        Condensed distance matrix

    Notes
    -----
    A non-NumPy array-API buffer (e.g. GPU-resident) is reordered on its own device
    via the array namespace; a NumPy array uses the parallel Cython kernel.

    This function uses parallel computation for improved performance.
    See the :install:`parallelization guide <#parallelization>` for information on
    controlling the number of threads used.

    """
    if _aac.is_array_api_obj(in_mat) and not _aac.is_numpy_array(in_mat):
        # Non-NumPy (e.g. GPU-resident) buffer: reorder then take the row-major
        # upper triangle (the condensed_form order), on the input's device.
        # Note: the boolean-mask gather relies on library-specific boolean
        # indexing (supported by numpy/torch/jax/cupy), not the strict array API.
        if validate:
            _validate_order(np.asarray(reorder_vec, dtype=np.intp), in_mat)
        xp = _aac.array_namespace(in_mat)
        dev = getattr(in_mat, "device", None)
        order = xp.asarray(np.asarray(reorder_vec, dtype=int), device=dev)
        reordered = in_mat[order][:, order]
        n = reordered.shape[0]
        idx = xp.arange(n, device=dev)
        return reordered[idx[:, None] < idx[None, :]]

    np_reorder = np.asarray(reorder_vec, dtype=np.intp)
    if validate:
        _validate_order(np_reorder, in_mat)

    if not in_mat.flags.c_contiguous:
        in_mat = np.asarray(in_mat, order="C")

    csize = ((np_reorder.size - 1) * np_reorder.size) // 2
    out_mat_condensed = np.empty([csize], in_mat.dtype)
    distmat_reorder_condensed_cy(in_mat, np_reorder, out_mat_condensed)
    return out_mat_condensed
