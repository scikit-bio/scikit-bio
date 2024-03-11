# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import partial

import numpy as np

from ._base import (
    _preprocess_input_sng,
    _run_monte_carlo_stats,
    _build_results,
    DistanceMatrix,
)
from ._cutils import permanova_f_stat_sW_cy


def permanova(distance_matrix, grouping, column=None, permutations=999):
    """Test for significant differences between groups using PERMANOVA.

    Permutational Multivariate Analysis of Variance (PERMANOVA) is a
    non-parametric method that tests whether two or more groups of objects
    (e.g., samples) are significantly different based on a categorical factor.
    It is conceptually similar to ANOVA except that it operates on a distance
    matrix, which allows for multivariate analysis. PERMANOVA computes a
    pseudo-F statistic.

    Statistical significance is assessed via a permutation test. The assignment
    of objects to groups (`grouping`) is randomly permuted a number of times
    (controlled via `permutations`). A pseudo-F statistic is computed for each
    permutation and the p-value is the proportion of permuted pseudo-F
    statisics that are equal to or greater than the original (unpermuted)
    pseudo-F statistic.

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
    anosim

    Notes
    -----
    See [1]_ for the original method reference, as well as ``vegan::adonis``,
    available in R's vegan package [2]_.

    The p-value will be ``np.nan`` if `permutations` is zero.

    References
    ----------
    .. [1] Anderson, Marti J. "A new method for non-parametric multivariate
       analysis of variance." Austral Ecology 26.1 (2001): 32-46.

    .. [2] http://cran.r-project.org/web/packages/vegan/index.html

    Examples
    --------
    See :mod:`skbio.stats.distance.anosim` for usage examples (both functions
    provide similar interfaces).

    """
    if not isinstance(distance_matrix, DistanceMatrix):
        raise TypeError("Input must be a DistanceMatrix.")
    sample_size = distance_matrix.shape[0]

    num_groups, grouping = _preprocess_input_sng(
        distance_matrix.ids, sample_size, grouping, column
    )

    # Calculate number of objects in each group.
    group_sizes = np.bincount(grouping)
    s_T = (distance_matrix[:] ** 2).sum() / sample_size
    # we are going over the whole matrix, instead of just upper triangle
    # so cut in half
    s_T /= 2.0

    test_stat_function = partial(
        _compute_f_stat, sample_size, num_groups, distance_matrix, group_sizes, s_T
    )
    stat, p_value = _run_monte_carlo_stats(test_stat_function, grouping, permutations)

    return _build_results(
        "PERMANOVA", "pseudo-F", sample_size, num_groups, stat, p_value, permutations
    )


def _compute_f_stat(
    sample_size, num_groups, distance_matrix, group_sizes, s_T, grouping
):
    """Compute PERMANOVA pseudo-F statistic."""
    # Calculate s_W for each group, accounting for different group sizes.
    s_W = permanova_f_stat_sW_cy(distance_matrix.data, group_sizes, grouping)

    s_A = s_T - s_W
    return (s_A / (num_groups - 1)) / (s_W / (sample_size - num_groups))
