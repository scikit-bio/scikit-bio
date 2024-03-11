# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import partial

import numpy as np
from scipy.stats import rankdata

from ._base import _preprocess_input, _run_monte_carlo_stats, _build_results


def anosim(distance_matrix, grouping, column=None, permutations=999):
    """Test for significant differences between groups using ANOSIM.

    Analysis of Similarities (ANOSIM) is a non-parametric method that tests
    whether two or more groups of objects (e.g., samples) are significantly
    different based on a categorical factor. The ranks of the distances in the
    distance matrix are used to calculate an R statistic, which ranges between
    -1 (anti-grouping) to +1 (strong grouping), with an R value of 0 indicating
    random grouping.

    Statistical significance is assessed via a permutation test. The assignment
    of objects to groups (`grouping`) is randomly permuted a number of times
    (controlled via `permutations`). An R statistic is computed for each
    permutation and the p-value is the proportion of permuted R statisics that
    are equal to or greater than the original (unpermuted) R statistic.

    Parameters
    ----------
    distance_matrix : DistanceMatrix
        Distance matrix containing distances between objects (e.g., distances
        between samples of microbial communities).
    grouping : 1-D array_like or pandas.DataFrame
        Vector indicating the assignment of objects to groups. For example,
        these could be strings or integers denoting which group an object
        belongs to. If `grouping` is 1-D ``array_like``, it must be the same
        length and in the same order as the objects in `distance_matrix`. If
        `grouping` is a ``DataFrame``, the column specified by `column` will be
        used as the grouping vector. The ``DataFrame`` must be indexed by the
        IDs in `distance_matrix` (i.e., the row labels must be distance matrix
        IDs), but the order of IDs between `distance_matrix` and the
        ``DataFrame`` need not be the same. All IDs in the distance matrix must
        be present in the ``DataFrame``. Extra IDs in the ``DataFrame`` are
        allowed (they are ignored in the calculations).
    column : str, optional
        Column name to use as the grouping vector if `grouping` is a
        ``DataFrame``. Must be provided if `grouping` is a ``DataFrame``.
        Cannot be provided if `grouping` is 1-D ``array_like``.
    permutations : int, optional
        Number of permutations to use when assessing statistical
        significance. Must be greater than or equal to zero. If zero,
        statistical significance calculations will be skipped and the p-value
        will be ``np.nan``.

    Returns
    -------
    pandas.Series
        Results of the statistical test, including ``test statistic`` and
        ``p-value``.

    See Also
    --------
    permanova

    Notes
    -----
    See [1]_ for the original method reference. The general algorithm and
    interface are similar to ``vegan::anosim``, available in R's vegan package
    [2]_.

    The p-value will be ``np.nan`` if `permutations` is zero.

    References
    ----------
    .. [1] Clarke, KR. "Non-parametric multivariate analyses of changes in
       community structure." Australian journal of ecology 18.1 (1993):
       117-143.

    .. [2] http://cran.r-project.org/web/packages/vegan/index.html

    Examples
    --------
    Load a 4x4 distance matrix and grouping vector denoting 2 groups of
    objects:

    >>> from skbio import DistanceMatrix
    >>> dm = DistanceMatrix([[0, 1, 1, 4],
    ...                      [1, 0, 3, 2],
    ...                      [1, 3, 0, 3],
    ...                      [4, 2, 3, 0]],
    ...                     ['s1', 's2', 's3', 's4'])
    >>> grouping = ['Group1', 'Group1', 'Group2', 'Group2']

    Run ANOSIM using 99 permutations to calculate the p-value:

    >>> import numpy as np
    >>> # make output deterministic; not necessary for normal use
    >>> np.random.seed(0)
    >>> from skbio.stats.distance import anosim
    >>> anosim(dm, grouping, permutations=99)
    method name               ANOSIM
    test statistic name            R
    sample size                    4
    number of groups               2
    test statistic              0.25
    p-value                     0.67
    number of permutations        99
    Name: ANOSIM results, dtype: object

    The return value is a ``pandas.Series`` object containing the results of
    the statistical test.

    To suppress calculation of the p-value and only obtain the R statistic,
    specify zero permutations:

    >>> anosim(dm, grouping, permutations=0)
    method name               ANOSIM
    test statistic name            R
    sample size                    4
    number of groups               2
    test statistic              0.25
    p-value                      NaN
    number of permutations         0
    Name: ANOSIM results, dtype: object

    You can also provide a ``pandas.DataFrame`` and a column denoting the
    grouping instead of a grouping vector. The following ``DataFrame``'s
    ``Group`` column specifies the same grouping as the vector we used in the
    previous examples:

    >>> # make output deterministic; not necessary for normal use
    >>> np.random.seed(0)
    >>> import pandas as pd
    >>> df = pd.DataFrame.from_dict(
    ...     {'Group': {'s2': 'Group1', 's3': 'Group2', 's4': 'Group2',
    ...                's5': 'Group3', 's1': 'Group1'}})
    >>> anosim(dm, df, column='Group', permutations=99)
    method name               ANOSIM
    test statistic name            R
    sample size                    4
    number of groups               2
    test statistic              0.25
    p-value                     0.67
    number of permutations        99
    Name: ANOSIM results, dtype: object

    The results match the first example above.

    Note that when providing a ``DataFrame``, the ordering of rows and/or
    columns does not affect the grouping vector that is extracted. The
    ``DataFrame`` must be indexed by the distance matrix IDs (i.e., the row
    labels must be distance matrix IDs).

    If IDs (rows) are present in the ``DataFrame`` but not in the distance
    matrix, they are ignored. The previous example's ``s5`` ID illustrates this
    behavior: note that even though the ``DataFrame`` had 5 objects, only 4
    were used in the test (see the "Sample size" row in the results above to
    confirm this). Thus, the ``DataFrame`` can be a superset of the distance
    matrix IDs. Note that the reverse is not true: IDs in the distance matrix
    *must* be present in the ``DataFrame`` or an error will be raised.

    """
    sample_size, num_groups, grouping, tri_idxs, distances = _preprocess_input(
        distance_matrix, grouping, column
    )

    divisor = sample_size * ((sample_size - 1) / 4)
    ranked_dists = rankdata(distances, method="average")

    test_stat_function = partial(_compute_r_stat, tri_idxs, ranked_dists, divisor)
    stat, p_value = _run_monte_carlo_stats(test_stat_function, grouping, permutations)

    return _build_results(
        "ANOSIM", "R", sample_size, num_groups, stat, p_value, permutations
    )


def _compute_r_stat(tri_idxs, ranked_dists, divisor, grouping):
    """Compute ANOSIM R statistic (between -1 and +1)."""
    # Create a matrix where True means that the two objects are in the same
    # group. This ufunc requires that grouping is a numeric vector (e.g., it
    # won't work with a grouping vector of strings).
    grouping_matrix = np.equal.outer(grouping, grouping)

    # Extract upper triangle from the grouping matrix. It is important to
    # extract the values in the same order that the distances are extracted
    # from the distance matrix (see ranked_dists). Extracting the upper
    # triangle (excluding the diagonal) preserves this order.
    grouping_tri = grouping_matrix[tri_idxs]

    # within
    r_W = np.mean(ranked_dists[grouping_tri])

    # between
    r_B = np.mean(ranked_dists[np.invert(grouping_tri)])

    return (r_B - r_W) / divisor
