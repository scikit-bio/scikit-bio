# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import partial

import numpy as np
import pandas as pd
from scipy.stats import f_oneway
from scipy.spatial.distance import cdist
from ._cutils import geomedian_axis_one


from ._base import (
    _preprocess_input_sng,
    _run_monte_carlo_stats,
    _build_results,
    DistanceMatrix,
)

from skbio.stats.ordination import pcoa, OrdinationResults


def permdisp(
    distance_matrix,
    grouping,
    column=None,
    test="median",
    permutations=999,
    method="eigh",
    number_of_dimensions=10,
):
    """Test for Homogeneity of Multivariate Groups Disperisons.

    PERMDISP is a multivariate analogue of Levene's test for homogeneity of
    multivariate variances. Distances are handled by reducing the
    original distances to principal coordinates. PERMDISP calculates an
    F-statistic to assess whether the dispersions between groups is significant


    Parameters
    ----------
    distance_matrix : DistanceMatrix or OrdinationResults
        Distance matrix containing distances between objects (e.g., distances
        between samples of microbial communities) or
        result of pcoa on such a matrix.
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
    method : str, optional
        Eigendecomposition method to use in performing PCoA.
        By default, uses SciPy's `eigh`, which computes exact
        eigenvectors and eigenvalues for all dimensions. The alternate
        method, `fsvd`, uses faster heuristic eigendecomposition but loses
        accuracy. The magnitude of accuracy lost is dependent on dataset.
        Note that using `fsvd` is still considered experimental and
        should be used with care.
        Not used if distance_matrix is a OrdinationResults object.
    number_of_dimensions : int, optional
        Dimensions to reduce the distance matrix to if using the `fsvd` method.
        Not used if the `eigh` method is being selected.

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
        If the test is not centroid or median,
        or if method is not eigh or fsvd
    TypeError
        If the distance matrix is not an instance of a
        ``skbio.DistanceMatrix``.
    ValueError
        If there is only one group
    ValueError
        If a list and a column name are both provided
    ValueError
        If a list is provided for `grouping` and it's length does not match
        the number of ids in distance_matrix
    ValueError
        If all of the values in the grouping vector are unique
    KeyError
        If there are ids in grouping that are not in distance_matrix

    See Also
    --------
    permanova
    anosim

    Notes
    -----
    This function uses Marti Anderson's PERMDISP2 procedure.

    The significance of the results from this function will be the same as the
    results found in vegan's betadisper, however due to floating point
    variability the F-statistic results may vary slightly.

    See [1]_ for the original method reference, as well as
    ``vegan::betadisper``, available in R's vegan package [2]_.

    References
    ----------
    .. [1] Anderson, M. J. (2006). Distance-based tests for homogeneity of multivariate
       dispersions. Biometrics, 62(1), 245-253.

    .. [2] http://cran.r-project.org/web/packages/vegan/index.html

    Examples
    --------
    Load a 6x6 distance matrix and grouping vector denoting 2 groups of
    objects:

    >>> from skbio import DistanceMatrix
    >>> dm = DistanceMatrix([[0,    0.5,  0.75, 1, 0.66, 0.33],
    ...                       [0.5,  0,    0.25, 0.33, 0.77, 0.61],
    ...                       [0.75, 0.25, 0,    0.1, 0.44, 0.55],
    ...                       [1,    0.33, 0.1,  0, 0.75, 0.88],
    ...                       [0.66, 0.77, 0.44, 0.75, 0, 0.77],
    ...                       [0.33, 0.61, 0.55, 0.88, 0.77, 0]],
    ...                       ['s1', 's2', 's3', 's4', 's5', 's6'])
    >>> grouping = ['G1', 'G1', 'G1', 'G2', 'G2', 'G2']

    Run PERMDISP using 99 permutations to caluculate the p-value:

    >>> from skbio.stats.distance import permdisp
    >>> import numpy as np
    >>> #make output deterministic, should not be included during normal use
    >>> np.random.seed(0)
    >>> permdisp(dm, grouping, permutations=99)
    method name               PERMDISP
    test statistic name        F-value
    sample size                      6
    number of groups                 2
    test statistic     ... 1.03...
    p-value            ...
    number of permutations          99
    Name: PERMDISP results, dtype: object

    The return value is a ``pandas.Series`` object containing the results of
    the statistical test.

    To suppress calculation of the p-value and only obtain the F statistic,
    specify zero permutations:

    >>> permdisp(dm, grouping, permutations=0)
    method name               PERMDISP
    test statistic name        F-value
    sample size                      6
    number of groups                 2
    test statistic      ... 1.03...
    p-value                        NaN
    number of permutations           0
    Name: PERMDISP results, dtype: object

    PERMDISP computes variances based on two types of tests, using either
    centroids or spatial medians, also commonly referred to as a geometric
    median. The spatial median is thought to yield a more robust test
    statistic, and this test is used by default. Spatial medians are computed
    using an iterative algorithm to find the optimally minimum point from all
    other points in a group while centroids are computed using a deterministic
    formula. As such the two different tests yeild slightly different F
    statistics.

    >>> np.random.seed(0)
    >>> permdisp(dm, grouping, test='centroid', permutations=6)
    method name               PERMDISP
    test statistic name        F-value
    sample size                      6
    number of groups                 2
    test statistic     ... 3.67...
    p-value            ... 0.42...
    number of permutations           6
    Name: PERMDISP results, dtype: object

    You can also provide a ``pandas.DataFrame`` and a column denoting the
    grouping instead of a grouping vector. The following DataFrame's
    Grouping column specifies the same grouping as the vector we used in the
    previous examples.:

    >>> import pandas as pd
    >>> df = pd.DataFrame.from_dict(
    ...      {'Grouping': {'s1': 'G1', 's2': 'G1', 's3': 'G1', 's4': 'G2',
    ...                    's5': 'G2', 's6': 'G2'}})
    >>> permdisp(dm, df, 'Grouping', permutations=6, test='centroid')
    method name               PERMDISP
    test statistic name        F-value
    sample size                      6
    number of groups                 2
    test statistic      ... 3.67...
    p-value             ... 0.42...
    number of permutations           6
    Name: PERMDISP results, dtype: object

    Note that when providing a ``DataFrame``, the ordering of rows and/or
    columns does not affect the grouping vector that is extracted. The
    ``DataFrame`` must be indexed by the distance matrix IDs (i.e., the row
    labels must be distance matrix IDs).

    If IDs (rows) are present in the ``DataFrame`` but not in the distance
    matrix, they are ignored. The previous example's ``s7`` ID illustrates this
    behavior: note that even though the ``DataFrame`` had 7 objects, only 6
    were used in the test (see the "Sample size" row in the results above to
    confirm this). Thus, the ``DataFrame`` can be a superset of the distance
    matrix IDs. Note that the reverse is not true: IDs in the distance matrix
    *must* be present in the ``DataFrame`` or an error will be raised.

    PERMDISP should be used to determine whether the dispersions between the
    groups in your distance matrix are significantly separated.
    A non-significant test result indicates that group dispersions are similar
    to each other. PERMANOVA or ANOSIM should then be used in conjunction to
    determine whether clustering within groups is significant.

    """
    if test not in ["centroid", "median"]:
        raise ValueError("Test must be centroid or median")

    if isinstance(distance_matrix, OrdinationResults):
        ordination = distance_matrix
        ids = ordination.samples.axes[0].to_list()
        sample_size = len(ids)
        distance_matrix = None  # not used anymore, avoid using by mistake
    elif isinstance(distance_matrix, DistanceMatrix):
        if method == "eigh":
            # eigh does not natively support specifying number_of_dimensions
            # and pcoa expects it to be 0
            number_of_dimensions = 0
        elif method != "fsvd":
            raise ValueError("Method must be eigh or fsvd")

        ids = distance_matrix.ids
        sample_size = distance_matrix.shape[0]

        ordination = pcoa(
            distance_matrix, method=method, number_of_dimensions=number_of_dimensions
        )
    else:
        raise TypeError("Input must be a DistanceMatrix or OrdinationResults.")

    samples = ordination.samples

    num_groups, grouping = _preprocess_input_sng(ids, sample_size, grouping, column)

    test_stat_function = partial(_compute_groups, samples, test)

    stat, p_value = _run_monte_carlo_stats(test_stat_function, grouping, permutations)

    return _build_results(
        "PERMDISP", "F-value", sample_size, num_groups, stat, p_value, permutations
    )


def _compute_groups(samples, test_type, grouping):
    groups = []

    samples["grouping"] = grouping
    if test_type == "centroid":
        centroids = samples.groupby("grouping").aggregate("mean")
    elif test_type == "median":
        grouping_cols = samples.columns.to_list()
        centroids = samples.groupby("grouping")[grouping_cols].apply(_config_med)

    for label, df in samples.groupby("grouping"):
        groups.append(
            cdist(
                df.values[:, :-1].astype("float64"),
                [centroids.loc[label].values],
                metric="euclidean",
            )
        )

    stat, _ = f_oneway(*groups)
    stat = stat[0]

    return stat


def _config_med(x):
    """Slice and transpose the vector.

    Slice the vector up to the last value to exclude grouping column
    and transpose the vector to be compatible with hd.geomedian.
    """
    X = x.values[:, :-1]
    return pd.Series(np.array(geomedian_axis_one(X.T)), index=x.columns[:-1])
