# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from ._cutils import geomedian_axis_one
from ._base import (
    _preprocess_input_sng,
    _run_monte_carlo_stats,
    _build_results,
    DistanceMatrix,
)
from skbio.stats.ordination import pcoa, OrdinationResults
from skbio.util._array import ingest_array
from skbio.util._decorator import params_aliased

if TYPE_CHECKING:  # pragma: no cover
    from numpy.typing import ArrayLike
    from skbio.util._typing import SeedLike


@params_aliased(
    [
        ("dimensions", "number_of_dimensions", "0.7.0", False),
        ("distmat", "distance_matrix", "0.7.0", False),
    ]
)
def permdisp(
    distmat: DistanceMatrix | OrdinationResults,
    grouping: pd.DataFrame | ArrayLike,
    column: str | None = None,
    test: str = "median",
    permutations: int = 999,
    method: str = "eigh",
    dimensions: int = 10,
    seed: SeedLike | None = None,
    warn_neg_eigval: bool | float = 0.01,
) -> pd.Series:
    r"""Test for Homogeneity of Multivariate Groups Disperisons.

    PERMDISP is a multivariate analog of Levene's test for homogeneity of multivariate
    variances. Distances are handled by reducing the original distances to principal
    coordinates. PERMDISP calculates an F-statistic to assess whether the dispersions
    between groups is significant.

    Parameters
    ----------
    distmat : DistanceMatrix or OrdinationResults
        Distance matrix containing distances between objects (e.g., distances
        between samples of microbial communities) or result of pcoa on such a matrix.
    grouping : 1-D array_like or pandas.DataFrame
        Vector indicating the assignment of objects to groups. For example,
        these could be strings or integers denoting which group an object
        belongs to. If `grouping` is 1-D ``array_like``, it must be the same
        length and in the same order as the objects in `distmat`. If
        `grouping` is a ``DataFrame``, the column specified by `column` will be
        used as the grouping vector. The ``DataFrame`` must be indexed by the
        IDs in `distmat` (i.e., the row labels must be distance matrix
        IDs), but the order of IDs between `distmat` and the
        ``DataFrame`` need not be the same. All IDs in the distance matrix must
        be present in the ``DataFrame``. Extra IDs in the ``DataFrame`` are
        allowed (they are ignored in the calculations).
    column : str, optional
        Column name to use as the grouping vector if `grouping` is a
        ``DataFrame``. Must be provided if `grouping` is a ``DataFrame``.
        Cannot be provided if `grouping` is 1-D ``array_like``.
    test : {'centroid', 'median'}, optional
        Determines whether the analysis is done using centroid or spatial
        median (default).
    permutations : int, optional
        Number of permutations to use when assessing statistical
        significance. Must be greater than or equal to zero. If zero,
        statistical significance calculations will be skipped and the p-value
        will be ``np.nan``.
    method : {'eigh', 'fsvd'}, optional
        Matrix decomposition method to use. Options are "eigh" (eigendecomposition,
        default) and "fsvd" (fast singular value decomposition). See
        :func:`~skbio.stats.ordination.pcoa <pcoa>` for details. Not used if
        distmat is a OrdinationResults object.
    dimensions : int, optional
        Dimensions to reduce the distance matrix to if using the `fsvd` method.
        Not used if the `eigh` method is being selected.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

        .. versionadded:: 0.6.3

    warn_neg_eigval : bool or float, optional
        Raise a warning if any negative eigenvalue of large magnitude is generated
        during PCoA. See :func:`skbio.stats.ordination.pcoa <pcoa>` for details.

        .. versionadded:: 0.6.3

    Returns
    -------
    pandas.Series
        Results of the statistical test, including ``test statistic`` and ``p-value``.

    Raises
    ------
    TypeError
        If, when using the spatial median test, the pcoa ordination is not of
        type np.float32 or np.float64, the spatial median function will fail
        and the centroid test should be used instead
    ValueError
        If the test is not centroid or median, or if method is not eigh or fsvd.
    TypeError
        If the distance matrix is not an instance of ``DistanceMatrix``.
    ValueError
        If there is only one group.
    ValueError
        If a list and a column name are both provided.
    ValueError
        If a list is provided for `grouping` and it's length does not match.
        the number of ids in distmat.
    ValueError
        If all of the values in the grouping vector are unique.
    KeyError
        If there are ids in grouping that are not in distmat.

    See Also
    --------
    permanova
    anosim

    Notes
    -----
    This function uses parallel computation for improved performance.
    See the :install:`parallelization guide <#parallelization>` for information on
    controlling the number of threads used.

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
    Load a 6x6 distance matrix and grouping vector denoting 2 groups of objects:

    >>> from skbio import DistanceMatrix
    >>> dm = DistanceMatrix([[0,    0.5,  0.75, 1, 0.66, 0.33],
    ...                       [0.5,  0,    0.25, 0.33, 0.77, 0.61],
    ...                       [0.75, 0.25, 0,    0.1, 0.44, 0.55],
    ...                       [1,    0.33, 0.1,  0, 0.75, 0.88],
    ...                       [0.66, 0.77, 0.44, 0.75, 0, 0.77],
    ...                       [0.33, 0.61, 0.55, 0.88, 0.77, 0]],
    ...                       ['s1', 's2', 's3', 's4', 's5', 's6'])
    >>> grouping = ['G1', 'G1', 'G1', 'G2', 'G2', 'G2']

    Run PERMDISP using 99 permutations to calculate the p-value. The seed is to make
    the output deterministic. You may skip it if that's not necessary.

    >>> from skbio.stats.distance import permdisp
    >>> permdisp(dm, grouping, permutations=99, seed=42) # doctest: +ELLIPSIS
    method name               PERMDISP
    test statistic name        F-value
    sample size                      6
    number of groups                 2
    test statistic             1.03296
    p-value                       ...
    number of permutations          99
    Name: PERMDISP results, dtype: object

    The return value is a ``pandas.Series`` object containing the results of the
    statistical test.

    To suppress calculation of the p-value and only obtain the F statistic, specify
    zero permutations:

    >>> permdisp(dm, grouping, permutations=0)
    method name               PERMDISP
    test statistic name        F-value
    sample size                      6
    number of groups                 2
    test statistic             1.03296
    p-value                        NaN
    number of permutations           0
    Name: PERMDISP results, dtype: object

    PERMDISP computes variances based on two types of tests, using either
    centroids or spatial medians, also commonly referred to as a geometric
    median. The spatial median is thought to yield a more robust test
    statistic, and this test is used by default. Spatial medians are computed
    using an iterative algorithm to find the optimally minimum point from all
    other points in a group while centroids are computed using a deterministic
    formula. As such the two different tests yield slightly different F
    statistics.

    >>> permdisp(dm, grouping, test='centroid', permutations=6, seed=42)
    method name               PERMDISP
    test statistic name        F-value
    sample size                      6
    number of groups                 2
    test statistic            3.670816
    p-value                   0.285714
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
    >>> permdisp(dm, df, 'Grouping', permutations=6, test='centroid', seed=42)
    method name               PERMDISP
    test statistic name        F-value
    sample size                      6
    number of groups                 2
    test statistic            3.670816
    p-value                   0.285714
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
    if test not in ("centroid", "median"):
        raise ValueError("Test must be centroid or median.")

    if isinstance(distmat, OrdinationResults):
        ordination = distmat
        ids = _extract_sample_ids(ordination, grouping)
        sample_size = len(ids)
    elif isinstance(distmat, DistanceMatrix):
        if method == "eigh":
            # eigh does not natively support specifying dimensions
            # and pcoa expects it to be 0
            dimensions = 0
        elif method != "fsvd":
            raise ValueError("Method must be eigh or fsvd.")

        ids = distmat.ids
        sample_size = distmat.shape[0]

        ordination = pcoa(
            distmat,
            method=method,
            dimensions=dimensions,
            warn_neg_eigval=warn_neg_eigval,
        )
    else:
        raise TypeError("Input must be a DistanceMatrix or OrdinationResults.")

    sample_data = _extract_sample_data(ordination, sample_size)

    num_groups, grouping = _preprocess_input_sng(ids, sample_size, grouping, column)

    test_stat_function = partial(_compute_groups, sample_data, test)

    stat, p_value = _run_monte_carlo_stats(
        test_stat_function, grouping, permutations, seed
    )

    return _build_results(
        "PERMDISP", "F-value", sample_size, num_groups, stat, p_value, permutations
    )


def _compute_groups(samples, test_type, grouping):
    xp, data = ingest_array(samples)
    groups = []

    grouping_array = np.asarray(grouping)
    if np.issubdtype(grouping_array.dtype, np.number):
        group_codes = grouping_array
    else:
        group_codes = _encode_grouping_labels(grouping)

    for group_id in np.unique(group_codes):
        group_data = data[group_codes == group_id]

        if test_type == "centroid":
            center = xp.mean(group_data, axis=0)
        else:
            if isinstance(group_data, np.ndarray):
                center_np = np.asarray(
                    geomedian_axis_one(group_data.T), dtype=np.float64
                )
                center = xp.asarray(center_np, dtype=group_data.dtype)
            else:
                center = _spatial_median(group_data, xp)

        diff = group_data - center
        distances = _vector_norm(xp, diff, axis=1)
        groups.append(distances)

    return _f_oneway_stat(xp, groups)


def _extract_sample_ids(ordination, grouping):
    sample_ids = ordination.sample_ids
    if sample_ids is not None:
        return list(sample_ids)

    samples = ordination.samples
    if hasattr(samples, "axes"):
        return samples.axes[0].to_list()

    if isinstance(grouping, (pd.DataFrame, pd.Series)):
        raise ValueError(
            "OrdinationResults with array-backed samples require `sample_ids` "
            "when grouping is provided as a DataFrame or Series."
        )

    _, sample_data = ingest_array(samples, to_numpy=True)
    return [str(i) for i in range(sample_data.shape[0])]


def _extract_sample_data(ordination, sample_size):
    samples = ordination.samples
    xp, sample_data = ingest_array(samples)

    if sample_data.ndim != 2:
        raise ValueError("Ordination sample coordinates must be two-dimensional.")
    if sample_data.shape[0] != sample_size:
        raise ValueError(
            "Ordination sample coordinate count must match the number of sample IDs."
        )

    try:
        dtype = np.dtype(sample_data.dtype)
    except (TypeError, AttributeError):
        dtype = None

    if dtype is None or not np.issubdtype(dtype, np.floating):
        sample_data = xp.asarray(sample_data, dtype=xp.float64)

    return sample_data


def _vector_norm(xp, arr, axis=None):
    try:
        return xp.linalg.vector_norm(arr, axis=axis)
    except AttributeError:
        return xp.sqrt(xp.sum(arr * arr, axis=axis))


def _f_oneway_stat(xp, groups):
    n_groups = len(groups)
    total_n = 0
    mean_num = 0.0
    group_means = []
    ss_within = None

    for group in groups:
        n_i = group.shape[0]
        mean_i = xp.mean(group)
        group_means.append((n_i, mean_i))
        mean_num = mean_num + mean_i * n_i

        centered = group - mean_i
        ss_i = xp.sum(centered * centered)
        ss_within = ss_i if ss_within is None else ss_within + ss_i
        total_n += n_i

    overall_mean = mean_num / total_n
    ss_between = None
    for n_i, mean_i in group_means:
        diff = mean_i - overall_mean
        term = n_i * diff * diff
        ss_between = term if ss_between is None else ss_between + term

    df_between = n_groups - 1
    df_within = total_n - n_groups
    ms_between = ss_between / df_between
    ms_within = ss_within / df_within

    ms_between_value = _to_python_scalar(ms_between)
    ms_within_value = _to_python_scalar(ms_within)

    if ms_within_value == 0.0:
        if ms_between_value == 0.0:
            return np.nan
        return np.inf

    return _to_python_scalar(ms_between / ms_within)


def _spatial_median(data, xp, eps=1e-7, maxiters=500):
    n = data.shape[0]
    y = xp.mean(data, axis=0)

    if n == 1:
        return y

    for _ in range(maxiters):
        diffs = data - y
        dists = _vector_norm(xp, diffs, axis=1)
        mask = dists > eps
        nzeros = int(_to_python_scalar(xp.sum(xp.logical_not(mask))))

        if nzeros == n:
            break

        dinv = xp.where(mask, 1.0 / dists, 0.0)
        dinvs = xp.sum(dinv)
        dinvs_value = _to_python_scalar(dinvs)
        if dinvs_value == 0.0:
            break

        weights = dinv / dinvs
        weighted = data * xp.expand_dims(weights, axis=1)
        T = xp.sum(weighted, axis=0)

        if nzeros == 0:
            y1 = T
        else:
            R = (T - y) * dinvs
            r = _to_python_scalar(_vector_norm(xp, R))
            if r > eps:
                rinv = nzeros / r
            else:
                rinv = 0.0
            y1 = max(0.0, 1.0 - rinv) * T + min(1.0, rinv) * y

        if _to_python_scalar(_vector_norm(xp, y - y1)) < eps:
            return y1

        y = y1

    return y


def _to_python_scalar(value):
    _, value_np = ingest_array(value, to_numpy=True)
    return float(np.asarray(value_np, dtype=np.float64).reshape(-1)[0])


def _encode_grouping_labels(grouping):
    labels = np.asarray(grouping, dtype=object)
    codes = np.empty(labels.shape[0], dtype=np.intp)
    mapping = {}
    next_code = 0

    for idx, label in enumerate(labels):
        code = mapping.get(label)
        if code is None:
            code = next_code
            mapping[label] = code
            next_code += 1
        codes[idx] = code

    return codes

