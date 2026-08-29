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
from scipy.stats import f_oneway

from ._cutils import geomedian_axis_one
from ._base import (
    _preprocess_input_sng,
    _run_monte_carlo_stats,
    _build_results,
    DistanceMatrix,
)
from skbio.stats.ordination import pcoa, OrdinationResults
from skbio.util._decorator import params_aliased
from skbio._config import _resolve_engine

try:
    from numba import njit

    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False

if TYPE_CHECKING:  # pragma: no cover
    from numpy.typing import ArrayLike
    import pandas as pd
    from skbio.util._typing import SeedLike

# Must match the defaults of geomedian_axis_one's eps/maxiters in _cutils.pyx, so
# the two engines resolve to the same behavior when the caller doesn't override them.
_GEOMEDIAN_EPS = 1e-7
_GEOMEDIAN_MAXITERS = 500


if NUMBA_AVAILABLE:

    @njit
    def _geomedian_axis_one_nb(X, eps, maxiters):
        """Compute geometric median along axis 1 using Numba.

        This is a Numba port of the Weiszfeld algorithm implemented by the
        Cython ``geomedian_axis_one`` (see ``_cutils.pyx``). It always
        operates in float64. All loops run sequentially (no ``prange``):
        profiling showed that at the problem sizes involved here (tens to
        low hundreds of points per group, called thousands of times per
        ``permdisp`` call for the Monte Carlo permutations), the fixed
        per-entry dispatch overhead of Numba's parallel regions dominates
        the actual compute time, making ``parallel=True`` several times
        *slower* than a plain sequential ``@njit`` function.
        """
        p = X.shape[0]
        n = X.shape[1]

        # initial guess = row means (mean along axis 1)
        y = np.empty(p, dtype=np.float64)
        for j in range(p):
            s = 0.0
            for i in range(n):
                s += X[j, i]
            y[j] = s / n

        if n == 1:
            return y

        D = np.empty(n, dtype=np.float64)
        Dinv = np.empty(n, dtype=np.float64)
        W = np.empty(n, dtype=np.float64)
        T = np.empty(p, dtype=np.float64)
        y1 = np.empty(p, dtype=np.float64)

        for _ in range(maxiters):
            for i in range(n):
                d = 0.0
                for j in range(p):
                    tmp = X[j, i] - y[j]
                    d += tmp * tmp
                Di = np.sqrt(d)
                D[i] = Di
                if abs(Di) > eps:
                    Dinv[i] = 1.0 / Di
                else:
                    Dinv[i] = 0.0

            # accumulate in the same order as Cython's _sum, for bit-identical results
            Dinvs = 0.0
            for i in range(n):
                Dinvs += Dinv[i]

            for i in range(n):
                W[i] = Dinv[i] / Dinvs

            for j in range(p):
                total = 0.0
                for i in range(n):
                    if abs(D[i]) > eps:
                        total += W[i] * X[j, i]
                T[j] = total

            nzeros = n
            for i in range(n):
                if abs(D[i]) > eps:
                    nzeros -= 1

            if nzeros == 0:
                for j in range(p):
                    y1[j] = T[j]
            elif nzeros == n:
                break
            else:
                r = 0.0
                for j in range(p):
                    Rj = (T[j] - y[j]) * Dinvs
                    r += Rj * Rj
                r = np.sqrt(r)
                if r > eps:
                    rinv = nzeros / r
                else:
                    rinv = 0.0
                a = max(0.0, 1.0 - rinv)
                b = min(1.0, rinv)
                for j in range(p):
                    y1[j] = a * T[j] + b * y[j]

            dist = 0.0
            for j in range(p):
                tmp = y[j] - y1[j]
                dist += tmp * tmp
            dist = np.sqrt(dist)
            if dist < eps:
                break

            for j in range(p):
                y[j] = y1[j]

        return y


def _geomedian_axis_one(X, engine):
    """Compute the geometric median along axis 1 using the given engine.

    Parameters
    ----------
    X : array_like
        Array of shape ``(p, n)`` (``p`` dimensions, ``n`` points).
    engine : str
        Already-resolved compute engine, either ``"cython"`` or ``"numba"``.

    Returns
    -------
    numpy.ndarray
        The length-``p`` geometric median.
    """
    if engine == "numba":
        # X is typically a transposed view (group_data.T); asarray (unlike
        # ascontiguousarray) only casts dtype and leaves an already-float64 array
        # as-is, so this avoids forcing a full transpose-copy on every call. Numba
        # handles the resulting non-contiguous 2-D array directly.
        Xc = np.asarray(X, dtype=np.float64)
        return np.asarray(
            _geomedian_axis_one_nb(Xc, _GEOMEDIAN_EPS, _GEOMEDIAN_MAXITERS)
        )
    return np.asarray(geomedian_axis_one(X))


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
    engine: str | None = None,
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

    engine : {"cython", "numba"}, optional
        Compute engine to use for the geometric median (``test="median"``).
        ``"cython"`` (default) uses the Cython implementation. ``"numba"``
        uses the optional Numba implementation and requires Numba to be
        installed. If not provided, the global default is used (see
        :func:`skbio.set_config`). Ignored when ``test="centroid"``.

        .. versionadded:: 0.7.4

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
    ValueError
        If ``engine`` is not one of the supported compute engines.
    ImportError
        If ``engine="numba"`` is specified but Numba is not installed.

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

    This illustrative distance matrix is non-Euclidean. The examples below disable
    PCoA's negative-eigenvalue warning to keep the focus on PERMDISP usage. This does
    not alter the results.

    Run PERMDISP using 99 permutations to calculate the p-value. The seed is to make
    the output deterministic. You may skip it if that's not necessary.

    >>> from skbio.stats.distance import permdisp
    >>> permdisp(dm, grouping, permutations=99, seed=42,
    ...          dimensions=dm.shape[0],
    ...          warn_neg_eigval=False)  # doctest: +ELLIPSIS
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

    >>> permdisp(dm, grouping, permutations=0, dimensions=dm.shape[0],
    ...          warn_neg_eigval=False)
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

    >>> permdisp(dm, grouping, test='centroid', permutations=6, seed=42,
    ...          dimensions=dm.shape[0], warn_neg_eigval=False)
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
    >>> permdisp(dm, df, 'Grouping', permutations=6, test='centroid', seed=42,
    ...          dimensions=dm.shape[0], warn_neg_eigval=False)
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
        ids = ordination.samples.axes[0].to_list()
        sample_size = len(ids)
    elif isinstance(distmat, DistanceMatrix):
        if method not in ("eigh", "fsvd"):
            raise ValueError("Method must be eigh or fsvd.")
        # The `dimensions` argument is passed through to pcoa unchanged.
        # Since #2285 the eigh path in pcoa honors `dimensions` and can
        # yield substantial speedups on large distance matrices; the
        # previous code here forced `dimensions=0`, which both discarded
        # that speedup and tripped pcoa's "EIGH: since no value for
        # dimensions..." RuntimeWarning on every call (#2430). Users who
        # deliberately pass `dimensions=0` will still see that warning,
        # as intended.

        ids = distmat.ids
        sample_size = distmat.shape[0]

        ordination = pcoa(
            distmat,
            method=method,
            dimensions=dimensions,
            seed=seed,
            warn_neg_eigval=warn_neg_eigval,
        )
    else:
        raise TypeError("Input must be a DistanceMatrix or OrdinationResults.")

    sample_data = ordination.samples.to_numpy(copy=False)

    num_groups, grouping = _preprocess_input_sng(ids, sample_size, grouping, column)

    # Only resolve engine for test="median": it's the only test type that
    # reads it (see _compute_groups below), and the docstring documents
    # engine as ignored for test="centroid". Resolving unconditionally would
    # raise ImportError/ValueError for an unavailable/invalid engine even
    # when test="centroid" never uses it.
    if test == "median":
        engine = _resolve_engine(engine, ("cython", "numba"))

    test_stat_function = partial(_compute_groups, sample_data, test, engine)

    stat, p_value = _run_monte_carlo_stats(
        test_stat_function, grouping, permutations, seed
    )

    return _build_results(
        "PERMDISP", "F-value", sample_size, num_groups, stat, p_value, permutations
    )


def _compute_groups(samples, test_type, engine, grouping):
    data = np.asarray(samples)
    groups = []

    grouping_array = np.asarray(grouping)
    for group_id in np.unique(grouping_array):
        group_data = data[grouping_array == group_id]

        if test_type == "centroid":
            center = group_data.mean(axis=0)
        else:  # median
            center = _geomedian_axis_one(group_data.T, engine)

        # Distances from each sample in this group to the group center.
        groups.append(np.linalg.norm(group_data - center, axis=1))

    stat, _ = f_oneway(*groups)
    return float(np.ravel(stat)[0])
