# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import sys
from heapq import heappush, heappop
from collections import defaultdict
from copy import copy

import numpy as np
import scipy.sparse as sparse
from skbio.util import get_rng

from biom import subsample as biom_subsample


def isubsample(items, maximum, minimum=1, buf_size=1000, bin_f=None):
    """Randomly subsample items from bins, without replacement.

    Randomly subsample items without replacement from an unknown number of
    input items, that may fall into an unknown number of bins. This method is
    intended for data that either a) cannot fit into memory or b) subsampling
    collections of arbitrary datatypes.

    Parameters
    ----------
    items : Iterable
        The items to evaluate.
    maximum : unsigned int
        The maximum number of items per bin.
    minimum : unsigned int, optional
        The minimum number of items per bin. The default is 1.
    buf_size : unsigned int, optional
        The size of the random value buffer. This buffer holds the random
        values assigned to each item from items. In practice, it is unlikely
        that this value will need to change. Increasing it will require more
        resident memory, but potentially reduce the number of function calls
        made to the PRNG, whereas decreasing it will result in more function
        calls and lower memory overhead. The default is 1000.
    bin_f : function, optional
        Method to determine what bin an item is associated with. If None (the
        default), then all items are considered to be part of the same bin.
        This function will be provided with each entry in items, and must
        return a hashable value indicating the bin that that entry should be
        placed in.

    Returns
    -------
    generator
        (bin, item)

    Raises
    ------
    ValueError
        If ``minimum`` is > ``maximum``.
    ValueError
        If ``minimum`` < 1 or if ``maximum`` < 1.

    See Also
    --------
    subsample_counts

    Notes
    -----
    Randomly get up to ``maximum`` items for each bin. If the bin has less than
    ``maximum``, only those bins that have >= ``minimum`` items are
    returned.

    This method will at most hold ``maximum`` * N data, where N is the number
    of bins.

    All items associated to a bin have an equal probability of being retained.

    Examples
    --------
    Randomly keep up to 2 sequences per sample from a set of demultiplexed
    sequences:

    >>> from skbio.stats import isubsample
    >>> import numpy as np
    >>> np.random.seed(123)
    >>> seqs = [('sampleA', 'AATTGG'),
    ...         ('sampleB', 'ATATATAT'),
    ...         ('sampleC', 'ATGGCC'),
    ...         ('sampleB', 'ATGGCT'),
    ...         ('sampleB', 'ATGGCG'),
    ...         ('sampleA', 'ATGGCA')]
    >>> bin_f = lambda item: item[0]
    >>> for bin_, item in sorted(isubsample(seqs, 2, bin_f=bin_f)):
    ...     print(bin_, item[1])
    sampleA AATTGG
    sampleA ATGGCA
    sampleB ATATATAT
    sampleB ATGGCG
    sampleC ATGGCC

    Now, let's set the minimum to 2:

    >>> bin_f = lambda item: item[0]
    >>> for bin_, item in sorted(isubsample(seqs, 2, 2, bin_f=bin_f)):
    ...     print(bin_, item[1])
    sampleA AATTGG
    sampleA ATGGCA
    sampleB ATATATAT
    sampleB ATGGCG

    """
    if minimum > maximum:
        raise ValueError("minimum cannot be > maximum.")
    if minimum < 1 or maximum < 1:
        raise ValueError("minimum and maximum must be > 0.")
    if bin_f is None:

        def bin_f(x):
            return True

    # buffer some random values
    random_values = np.random.randint(0, sys.maxsize, buf_size, dtype=np.int64)
    random_idx = 0

    result = defaultdict(list)
    for item in items:
        bin_ = bin_f(item)
        heap = result[bin_]

        # pull a random value, and recompute random values if we've consumed
        # our buffer
        random_value = random_values[random_idx]
        random_idx += 1
        if random_idx >= buf_size:
            random_values = np.random.randint(0, sys.maxsize, buf_size, dtype=np.int64)
            random_idx = 0

        # push our item on to the heap and drop the smallest if necessary
        heappush(heap, (random_value, copy(item)))
        if len(heap) > maximum:
            heappop(heap)

    # yield items
    for bin_, heap in result.items():
        if len(heap) < minimum:
            continue

        for _, item in heap:
            yield (bin_, item)


def subsample_counts(counts, n, replace=False, seed=None):
    """Randomly subsample from a vector of counts, with or without replacement.

    Parameters
    ----------
    counts : 1-D array_like
        Vector of counts (integers or floats) to randomly subsample from.
    n : int
        Number of items to subsample from `counts`. Must be less than or equal
        to the sum of `counts`.
    replace : bool, optional
        If ``True``, subsample with replacement. If ``False`` (the default),
        subsample without replacement.
    seed : int or np.random.Generator, optional
        A user-provided random seed or random generator instance.

    Returns
    -------
    subsampled : ndarray
        Subsampled vector of counts where the sum of the elements equals `n`
        (i.e., ``subsampled.sum() == n``). Will have the same shape as
        `counts`.

    Raises
    ------
    ValueError
        If `n` is less than zero or greater than the sum of `counts`
        when `replace=False`.
    EfficiencyWarning
        If the accelerated code isn't present or hasn't been compiled.

    See Also
    --------
    isubsample
    skbio.diversity.alpha

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
    >>> from skbio.stats import subsample_counts
    >>> a = np.array([4, 5, 0, 2, 1])
    >>> sub = subsample_counts(a, 4)
    >>> sub.sum()
    4
    >>> sub.shape
    (5,)

    Trying to subsample an equal number of items (without replacement) results
    in the same vector as our input:

    >>> subsample_counts([0, 3, 0, 1], 4)
    array([0, 3, 0, 1])

    Subsample 5 items (with replacement):

    >>> sub = subsample_counts([1, 0, 1, 2, 2, 3, 0, 1], 5, replace=True)
    >>> sub.sum()
    5
    >>> sub.shape
    (8,)

    """
    if n < 0:
        raise ValueError("n cannot be negative.")

    # cast to float as that's what the biom subsample method currently requires
    counts = np.asarray(counts, dtype=float)

    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")

    # csr_matrix will report ndim of 2 if vector
    counts = sparse.csr_matrix(counts)

    counts_sum = counts.sum()
    if n > counts_sum and not replace:
        raise ValueError(
            "Cannot subsample more items than exist in input "
            "counts vector when `replace=False`."
        )

    rng = get_rng(seed)
    biom_subsample(counts, n, replace, rng)

    return np.atleast_1d(counts.astype(int).toarray().squeeze())
