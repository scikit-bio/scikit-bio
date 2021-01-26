# ----------------------------------------------------------------------------
# Copyright (c) 2021-2021, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from ._cutils import distmat_reorder_cy


def distmat_reorder_buf(in_mat, reorder_vec, out_mat, validate = false):
    """
    Reorder the rows and columns of a distance matrix
    given a reorder vector.
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
    out_mat : 2D array_like
        Output, Distance matrix, 
        must be in c_order and same size as reorder_vec
    validate: boolean
        Optional, if true, validate reorder_vec content, detaults to false
    """
    np_reorder = np.asarray(reorder_vec, dtype=np.int)
    if validate:
        maxsize = in_mat.shape[0]
        bad_cnt = np.where((np_reorder < 0) or (np_reorder >= maxsize))[0].size
        if bad_cnt > 0:
            raise ValueError("Invalid reorder_vec")

    distmat_reorder_cy(in_mat, np_reorder, out_mat)

def distmat_reorder(in_mat, reorder_vec, validate = false):
    """
    Reorder the rows and columns of a distance matrix
    given a reorder vector.
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
        Optional, if true, validate reorder_vec content, detaults to false

    Returns
    -------
    out_mat : 2D array_like
        Distance matrix
    """
    np_reorder = np.asarray(reorder_vec, dtype=np.int)
    if validate:
        maxsize = in_mat.shape[0]
        bad_cnt = np.where((np_reorder < 0) or (np_reorder >= maxsize))[0].size
        if bad_cnt > 0:
            raise ValueError("Invalid reorder_vec")

    out_mat = np.empty([np_reorder.size, np_reorder.size], in_mat.dtype)
    distmat_reorder_cy(in_mat, np_reorder, out_mat)
    return out_mat

