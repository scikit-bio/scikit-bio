#!/usr/bin/env python
r"""
Subsampling (:mod:`skbio.maths.subsample`)
==========================================

.. currentmodule:: skbio.maths.subsample

This module provides functionality for subsampling from vectors of counts.

Functions
---------

.. autosummary::
   :toctree: generated/

   subsample

"""
from __future__ import division

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
cimport numpy as cnp


def subsample(counts, n, replace=False):
    """Randomly subsample from a vector of counts, with or without replacement.

    Parameters
    ----------
    counts : 1-D array_like
        Vector of counts (integers) to randomly subsample from.
    n : int
        Number of items to subsample from `counts`. Must be less than or equal
        to the sum of `counts`.
    replace : bool, optional
        If ``True``, subsample with replacement. If ``False`` (the default),
        subsample without replacement.

    Returns
    -------
    subsampled : ndarray
        Subsampled vector of counts where the sum of the elements equals `n`
        (i.e., ``subsampled.sum() == n``). Will have the same shape as
        `counts`.

    Raises
    ------
    TypeError
        If `counts` cannot be safely converted to an integer datatype.
    ValueError
        If `n` is less than zero or greater than the sum of `counts`.

    Notes
    -----
    If subsampling is performed without replacement (``replace=False``), a copy
    of `counts` is returned if `n` is equal to the number of items in `counts`,
    as all items will be chosen from the original vector.

    If subsampling is performed with replacement (``replace=True``) and `n` is
    equal to the number of items in `counts`, the subsampled vector that is
    returned may not necessarily be the same vector as `counts`.

    Examples
    --------
    Subsample 4 items (without replacement) from a vector of counts:

    >>> import numpy as np
    >>> from skbio.maths.subsample import subsample
    >>> a = np.array([4, 5, 0, 2, 1])
    >>> sub = subsample(a, 4)
    >>> sub.sum()
    4
    >>> sub.shape
    (5,)

    Trying to subsample an equal number of items (without replacement) results
    in the same vector as our input:

    >>> subsample([0, 3, 0, 1], 4)
    array([0, 3, 0, 1])

    Subsample 5 items (with replacement):

    >>> sub = subsample([1, 0, 1, 2, 2, 3, 0, 1], 5, replace=True)
    >>> sub.sum()
    5
    >>> sub.shape
    (8,)

    """
    cdef:
        cnp.ndarray[cnp.int64_t, ndim=1] int_counts, result, permuted, unpacked
        cnp.int64_t cnt
        Py_ssize_t unpacked_idx, i, j

    if n < 0:
        raise ValueError("n cannot be negative.")

    counts = np.asarray(counts)
    int_counts = counts.astype(int, casting='safe')

    if int_counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")

    counts_sum = int_counts.sum()
    if n > counts_sum:
        raise ValueError("Cannot subsample more items than exist in input "
                         "counts vector.")

    if replace:
        probs = int_counts / counts_sum
        result = np.random.multinomial(n, probs)
    else:
        unpacked = np.empty(counts_sum, dtype=int)
        unpacked_idx = 0
        for i in range(int_counts.shape[0]):
            cnt = int_counts[i]
            for j in range(cnt):
                unpacked[unpacked_idx] = i
                unpacked_idx += 1

        permuted = np.random.permutation(unpacked)[:n]

        result = np.zeros_like(int_counts)
        for idx in range(permuted.shape[0]):
            result[permuted[idx]] += 1

    return result
