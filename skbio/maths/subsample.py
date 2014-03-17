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

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import numpy as np


def subsample(counts, n):
    """Subsample from a vector of counts.

    Returns `counts` if `n` is equal to or larger than the number of items in
    `counts`.

    Parameters
    ----------
    counts : 1-D array_like
        Vector of counts (integers) to subsample from.
    n : int
        Number of items to subsample from `counts`.

    Returns
    -------
    subsampled : ndarray
        Subsampled vector of counts where ``subsampled.sum() == n``.

    Examples
    --------
    Subsample 4 items from a vector of counts:

    >>> import numpy as np
    >>> from skbio.maths.subsample import subsample
    >>> a = np.array([4, 5, 0, 2, 1])
    >>> sub = subsample(a, 4)
    >>> sub.sum()
    4

    Trying to subsample an equal or greater number of items results in the same
    vector as our input:

    >>> subsample([0, 3, 0, 1], 8)
    array([0, 3, 0, 1])

    """
    counts = np.asarray(counts)
    counts = counts.astype(int, casting='safe')

    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")

    if counts.sum() <= n:
        return counts

    nz = counts.nonzero()[0]
    unpacked = np.concatenate([np.repeat(np.array(i,), counts[i]) for i in nz])
    permuted = np.random.permutation(unpacked)[:n]

    result = np.zeros(len(counts), dtype=int)
    for p in permuted:
        result[p] += 1

    return result
