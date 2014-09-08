# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import sys
from warnings import warn
from heapq import heappush, heappop

import numpy as np
from future.utils import viewitems
from collections import defaultdict
from copy import copy

from skbio.util import EfficiencyWarning
try:
    from .__subsample import _subsample_without_replacement
except ImportError:
    pass


def uneven_subsample(iter_, maximum, minimum=1, random_buf_size=100000,
                     bin_f=None):
    """Get a random subset of items per bin

    Parameters
    ----------
    iter_ : Iterable
        The iterable to walk over.
    maximum : unsigned int
        The maximum number of items per bin.
    minimum : unsigned int, optional
        The minimum number of items per bin. The default is 1.
    random_buf_size : unsigned int, optional
        The size of the random value buffer. THe default is 100000.
    bin_f : function, optional
        Method to determine what bin an item is associated with. If None (the
        default), then all items are considered to be part of the same bin.

    Notes
    -----
    Randomly get ``maximum`` items for each bin. If the bin has less than
    ``maximum``, only those bins that have > ``minimum`` items are
    returned.

    This method will at most hold ``maximum`` * N data, where N is the number
    of bins.

    All items associated to a bin have an equal probability of being retained.

    If ``maximum`` is equal to ``minimum``, then this method should be the
    same as ``subsample``.

    Raises
    ------
    ValueError
        If ``minimum`` is > ``maximum``.
    ValueError
        If ``minimum`` < 1 or if ``maximum`` < 1.

    Returns
    -------
    generator
        (bin, item)

    Examples
    --------
    Randomly keep up to 2 sequences per sample from a set of demultiplexed
    sequences:

    >>> from skbio.stats import uneven_subsample
    >>> import numpy as np
    >>> np.random.seed(123)
    >>> sequences = [('sampleA', 'AATTGG'),
    ...              ('sampleB', 'ATATATAT'),
    ...              ('sampleC', 'ATGGCC'),
    ...              ('sampleB', 'ATGGCT'),
    ...              ('sampleB', 'ATGGCG'),
    ...              ('sampleA', 'ATGGCA')]
    >>> bin_f = lambda item: item[0]
    >>> for bin_, item in sorted(uneven_subsample(sequences, 2, bin_f=bin_f)):
    ...     print(bin_, item[1])
    sampleA AATTGG
    sampleA ATGGCA
    sampleB ATATATAT
    sampleB ATGGCG
    sampleC ATGGCC
    """
    if minimum > maximum:
        raise ValueError("minimum cannot be > maximum!")
    if minimum < 1 or maximum < 1:
        raise ValueError("minimum and maximum must be > 0!")
    if bin_f is None:
        bin_f = lambda x: True

    # buffer some random values
    random_values = np.random.randint(0, sys.maxint, random_buf_size)
    random_idx = 0

    result = defaultdict(list)
    for item in iter_:
        bin_ = bin_f(item)
        heap = result[bin_]

        # pull a random value, and recompute random values if we've consumed
        # our buffer
        random_value = random_values[random_idx]
        random_idx += 1
        if random_idx >= random_buf_size:
            random_values = np.random.randint(0, sys.maxint, random_buf_size)
            random_idx = 0

        # push our item on to the heap and drop the smallest if necessary
        heappush(heap, (random_value, copy(item)))
        if len(heap) > maximum:
            heappop(heap)

    # yield items
    for bin_, heap in viewitems(result):
        if len(heap) < minimum:
            continue

        for _, item in heap:
            yield (bin_, item)


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
