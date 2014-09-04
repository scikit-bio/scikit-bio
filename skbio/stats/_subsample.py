# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from warnings import warn

import numpy as np

from skbio.util import EfficiencyWarning
try:
    from .__subsample import _subsample_without_replacement
except ImportError:
    pass


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

    Raises
    ------
    EfficiencyWarning
        If the accelerated code isn't present or hasn't been compiled.

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
    >>> from skbio.stats import subsample
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
    if n < 0:
        raise ValueError("n cannot be negative.")

    counts = np.asarray(counts)
    counts = counts.astype(int, casting='safe')

    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")

    counts_sum = counts.sum()
    if n > counts_sum:
        raise ValueError("Cannot subsample more items than exist in input "
                         "counts vector.")

    if replace:
        probs = counts / counts_sum
        result = np.random.multinomial(n, probs)
    else:
        if counts_sum == n:
            result = counts
        else:
            try:
                result = _subsample_without_replacement(counts, n, counts_sum)
            except NameError:
                warn("Accelerated subsampling without replacement isn't"
                     " available.", EfficiencyWarning)

                nz = counts.nonzero()[0]
                unpacked = np.concatenate([np.repeat(np.array(i,), counts[i])
                                           for i in nz])
                permuted = np.random.permutation(unpacked)[:n]

                result = np.zeros(len(counts), dtype=int)
                for p in permuted:
                    result[p] += 1

    return result
