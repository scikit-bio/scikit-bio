# ----------------------------------------------------------------------------
# Copyright (c) 2021-2021, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from ._cutils import is_symmetric_cy


def is_symmetric(mat):
    """
    Check if a Distance Matrix is symmetric
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
    # is_symmetric_cy is optimized for the common cas of c_contiguous
    # for all other cases, make a copy
    if not mat.flags.c_contiguous:
        mat = np.asarray(mat, order='C')

    return is_symmetric_cy(mat)
