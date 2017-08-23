# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import partial

import numpy as np
import pandas as pd
from scipy.stats import f_oneway
from scipy.spatial.distance import euclidean

import hdmedians as hd

from ._base import (_preprocess_input, _run_monte_carlo_stats, _build_results)

from skbio.stats.ordination import pcoa
from skbio.util._decorator import experimental


@experimental(as_of="0.5.1")
def permdisp(distance_matrix, grouping, column=None, test='median',
             permutations=999):

    """Test for Homogeneity of Multivariate Groups Disperisons using Marti
    Anderson's PERMDISP2 procedure. PERMDISP is a multivariate analogue of
    Levene's test for homogeneity of multivariate variances. Non-euclidean
    distances are handled by reducing the original distances to principle
    coordinates. PERMDISP calculates an F-statistic to assess whether the
    dispersions between groups is significant.

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
    test : {'centroid', 'median'}
        determines whether the analysis is done using centroid or spaitial
        median.
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

    Raises
    ------
    TypeError
        If, when using the spatial median test, the pcoa ordination is not of
        type np.float32 or np.float64, the spatial median function will fail
        and the centroid test should be used instead
    ValueError
        If the test is not centroid or median. test is set to median by default
    TypeError
        If the distance matrix is not an instance of a distance matrix

    See Also
    --------
    permanova

    Notes
    -----
    See [1]_ for the original method reference, as well as
    ``vegan::betadisper``, available in R's vegan package [2]_.

    References
    ----------
    .. [1] Anderson, Marti J. "Distance-Based Tests for Homogeneity of
        Multivariate Dispersions." Biometrics 62 (2006):245-253

    .. [2] http://cran.r-project.org/web/packages/vegan/index.html

    Examples
    --------
    Load a 4x4 distance matrix and grouping vector denoting 2 groups of
    objects:

    >>> from skbio import DistanceMatrix
    >>> dm = DistanceMatrix([0, 1, 1, 4],
                            [1, 0, 3, 2],
                            [1, 3, 0, 3],
                            [4, 2, 3, 0]],
                            ['s1', 's2', 's3', 's4'])
    >>> grouping = ['G1', 'G1', 'G2', 'G2']

    Run PERMDISP using 99 permutations to caluculate the p-value:

    >>> from skbio.stats.distance import permdisp
    >>> import numpy as np
    >>> #make output deterministic; not necessary for normal use
    >>> np.random.seed(0)
    >>> permdisp(dm, grouping, permutations=99)
    


    """
    ordination = pcoa(distance_matrix)

    sample_size, num_groups, grouping, tri_idxs, distances = _preprocess_input(
        distance_matrix, grouping, column)

    if test == 'centroid':
        test_stat_function = partial(_cen_oneway, ordination)
    elif test == 'median':
        test_stat_function = partial(_med_oneway, ordination)
    else:
        raise ValueError('Test must be centroid or median')

    stat, p_value = _run_monte_carlo_stats(test_stat_function, grouping,
                                           permutations)

    return _build_results('PERMDISP', 'F-value', sample_size, num_groups,
                          stat, p_value, permutations)


def _eu_dist(x, vector): # not explicitly tested
    """
    return a series of Euclidean distances from the aggregated series,
    sliced to exclude the grouping column to an established centroid or
    spatial median vector
    """
    return pd.Series([euclidean(x[:-1],
                     vector.loc[x.grouping]), x.grouping],
                     index=['distance', 'grouping'])


def _compute_centroid_groups(ordination, grouping):

    groups = []

    ordination.samples['grouping'] = grouping

    centroids = ordination.samples.groupby('grouping').aggregate(_centroid)

    grouped = ordination.samples.apply(_eu_dist, axis=1,
                                       vector=centroids).groupby('grouping')
    for _, group in grouped:
        groups.append(group['distance'].tolist())

    return groups


def _centroid(x): # not explicitly tested
    return x.sum()/len(x)


def _cen_oneway(ordination, grouping): # not explicitly tested
    stat, _ = f_oneway(*(_compute_centroid_groups(ordination, grouping)))
    return stat


def _config_med(x): # not explicitly tested
    """slice the vector up to the last value to exclude grouping column
    and transpose the vector to be compatible with hd.geomedian
    """
    X = x.values[:, :-1]
    return np.array(hd.geomedian(X.T))


def _compute_median_groups(ordination, grouping):

    groups = []

    ordination.samples['grouping'] = grouping

    medians = ordination.samples.groupby('grouping').aggregate(_config_med)

    grouped = ordination.samples.apply(_eu_dist, axis=1,
                                         vector=medians).groupby('grouping')

    for _, group in grouped:
        groups.append(group['distance'].tolist())

    return groups


def _med_oneway(ordination, grouping): # not explicitly tested
    stat, _ = f_oneway(*(_compute_median_groups(ordination, grouping)))
    return stat
