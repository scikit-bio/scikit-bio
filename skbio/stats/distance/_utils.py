# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from ._cutils import is_symmetric_and_hollow_cy


def is_symmetric_and_hollow(mat):
    """
    Check if a Distance Matrix is symmetric and hollow.
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
        mat = np.asarray(mat, order='C')

    return is_symmetric_and_hollow_cy(mat)


def is_symmetric(mat):
    """
    Check if a Distance Matrix is symmetric.
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
    """
    Check if a Distance Matrix is hollow.
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
    return (np.trace(mat) == 0)
