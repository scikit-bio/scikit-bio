# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import annotations

import math
from functools import partial
from warnings import warn
from typing import TYPE_CHECKING

import numpy as np

from ._base import (
    _preprocess_input_sng,
    _run_monte_carlo_stats,
    _build_results,
    DistanceMatrix,
)
from ._cutils import permanova_f_stat_sW_cy, permanova_f_stat_sW_condensed_cy
from skbio.binaries import (
    permanova_available as _skbb_permanova_available,
    permanova as _skbb_permanova,
)
from skbio.util import get_rng
from skbio.util._decorator import params_aliased
from skbio._config import _resolve_engine

try:
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False

if TYPE_CHECKING:  # pragma: no cover
    from numpy.typing import ArrayLike
    import pandas as pd
    from skbio.util._typing import SeedLike


if NUMBA_AVAILABLE:

    @njit(parallel=True)
    def _permanova_f_stat_sW_nb(distance_matrix, group_sizes, grouping):
        """Compute s_W for full (non-condensed) distance matrix using Numba.

        Replaces the following natural Numba code, which triggers a parfor
        cycle bug in Numba 0.65.x when the scalar reduction variable is read
        inside the inner loop::

            @njit(parallel=True)
            def _permanova_f_stat_sW_nb(
                distance_matrix, group_sizes, grouping
            ):
                n = distance_matrix.shape[0]
                s_W = 0.0
                for row_idx in prange(n - 1):
                    group_idx = grouping[row_idx]
                    local_sum = 0.0
                    for col_idx in range(row_idx + 1, n):
                        if grouping[col_idx] == group_idx:
                            val = distance_matrix[row_idx, col_idx]
                            local_sum += val * val
                    s_W += local_sum / group_sizes[group_idx]
                return s_W

        Instead, each parallel iteration writes its row-pair contribution to
        a partials array and the final sum is performed sequentially. The
        loop also pairs row_idx with its mirror (n - row_idx - 2) so each
        prange iteration handles a comparable amount of work, matching the
        Cython implementation's load-balancing pattern.

        Parameters
        ----------
        distance_matrix : np.ndarray, shape (n, n)
            Full symmetric distance matrix.
        group_sizes : np.ndarray, shape (num_groups,)
            Number of objects in each group (precomputed via np.bincount).
        grouping : np.ndarray, shape (n,)
            Integer group label for each object (values 0 to num_groups - 1).

        Returns
        -------
        float
            The within-group sum of squares s_W.

        """
        n = distance_matrix.shape[0]
        n_half = n // 2
        partials = np.zeros(n_half, np.float64)

        for row_idx in prange(n_half):
            group_idx = grouping[row_idx]
            local_sum = 0.0
            for col_idx in range(row_idx + 1, n):
                if grouping[col_idx] == group_idx:
                    # promote to float64 before squaring to match the Cython
                    # ``double`` accumulator (matters for float32 input)
                    val = np.float64(distance_matrix[row_idx, col_idx])
                    local_sum += val * val

            group_sum = local_sum / group_sizes[group_idx]

            mirror_row = n - row_idx - 2
            if mirror_row != row_idx:
                group_idx = grouping[mirror_row]
                local_sum = 0.0
                for col_idx in range(mirror_row + 1, n):
                    if grouping[col_idx] == group_idx:
                        val = np.float64(distance_matrix[mirror_row, col_idx])
                        local_sum += val * val
                group_sum += local_sum / group_sizes[group_idx]

            partials[row_idx] = group_sum

        return partials.sum()

    @njit(parallel=True)
    def _permanova_f_stat_sW_condensed_nb(condensed_matrix, group_sizes, grouping):
        """Compute s_W for condensed distance matrix using Numba.

        Same as ``_permanova_f_stat_sW_nb`` but accepts ``condensed_matrix``
        in condensed form (1-D array of length ``n * (n - 1) // 2``) instead
        of the full 2-D matrix. Uses the standard condensed-index formula to
        map ``(row_idx, col_idx)`` pairs back to the 1-D array. The same
        Numba 0.65.x parfor cycle workaround applies (partials array + final
        sequential sum); see the full version for the equivalent natural
        Numba code.

        Parameters
        ----------
        condensed_matrix : np.ndarray, shape (k,)
            Condensed (upper triangle) distance matrix where
            k = n * (n - 1) // 2.
        group_sizes : np.ndarray, shape (num_groups,)
            Number of objects in each group.
        grouping : np.ndarray, shape (n,)
            Integer group label for each object.

        Returns
        -------
        float
            The within-group sum of squares s_W.

        """
        k = condensed_matrix.shape[0]
        n = int((1.0 + math.sqrt(1.0 + 8.0 * k)) / 2.0)
        n_half = n // 2
        partials = np.zeros(n_half, np.float64)

        for row_idx in prange(n_half):
            group_idx = grouping[row_idx]
            local_sum = 0.0
            for col_idx in range(row_idx + 1, n):
                if grouping[col_idx] == group_idx:
                    condensed_idx = (
                        row_idx * n + col_idx
                        - ((row_idx + 2) * (row_idx + 1)) // 2
                    )
                    val = np.float64(condensed_matrix[condensed_idx])
                    local_sum += val * val

            group_sum = local_sum / group_sizes[group_idx]

            mirror_row = n - row_idx - 2
            if mirror_row != row_idx:
                group_idx = grouping[mirror_row]
                local_sum = 0.0
                for col_idx in range(mirror_row + 1, n):
                    if grouping[col_idx] == group_idx:
                        condensed_idx = (
                            mirror_row * n + col_idx
                            - ((mirror_row + 2) * (mirror_row + 1)) // 2
                        )
                        val = condensed_matrix[condensed_idx]
                        local_sum += val * val
                group_sum += local_sum / group_sizes[group_idx]

            partials[row_idx] = group_sum

        return partials.sum()

    @njit(parallel=True)
    def _permanova_f_stat_sW_rowtile_nb(
        distance_matrix, groupings_T, inv_group_sizes, partials
    ):
        """Compute partial s_W for a chunk of permutations in one matrix pass.

        Parallelises over distance-matrix rows (each paired with its mirror row
        ``n - row - 2`` for load balance, as in ``_permanova_f_stat_sW_nb``)
        rather than over permutations. Each row is owned by one worker, so each
        distance is read from memory once per chunk and reused across all
        permutations of the chunk, instead of streaming the whole matrix once
        per permutation. The kernel is memory-bound, so this reuse -- not the
        arithmetic -- is what determines its speed; it mirrors the cache reuse
        of scikit-bio-binaries' ``pmn_f_stat_sW_cpu`` in
        ``permanova_dyn_impl.hpp``.

        The permutation loop is kept branchless (a uniform read-and-accumulate
        over the grouping windows, with no data-dependent control flow), and
        ``groupings_T`` is sample-major ``int32`` so each permutation's labels
        are read as a contiguous, low-bandwidth window.

        Parameters
        ----------
        distance_matrix : np.ndarray, shape (n, n)
            Full symmetric distance matrix.
        groupings_T : np.ndarray, shape (n, n_perm), int32
            Permuted groupings, sample-major (transposed) and C-contiguous.
        inv_group_sizes : np.ndarray, shape (num_groups,)
            Reciprocal group sizes (``1.0 / group_sizes``), precomputed so the
            hot loop multiplies instead of dividing.
        partials : np.ndarray, shape (n // 2, C) with C >= n_perm
            Reused scratch (allocated once by the caller). Each iteration zeros
            and fills the first ``n_perm`` entries of its own row ``partials[i]``
            with the s_W contribution of the mirror-pair it handles. The caller
            sums the first ``n_perm`` columns over axis 0 to obtain s_W per
            permutation.

        """
        n = distance_matrix.shape[0]
        n_perm = groupings_T.shape[1]
        n_half = n // 2

        # Each worker owns one row (and its mirror for load balance), so each
        # distance is read once per chunk and reused across all permutations.
        for row_idx in prange(n_half):
            # partials is reused across chunks, so zero this row's accumulator.
            local_sW = partials[row_idx]
            for p in range(n_perm):
                local_sW[p] = 0.0
            rsum = np.empty(n_perm, np.float64)

            for p in range(n_perm):
                rsum[p] = 0.0
            for col in range(row_idx + 1, n):
                # promote to float64 before squaring to match the Cython
                # ``double`` accumulator (matters for float32 input)
                val = np.float64(distance_matrix[row_idx, col])
                val2 = val * val
                # branchless: scale val2 by the 0/1 group-match instead of an
                # if, so the loop is a uniform read-accumulate over the grouping
                # window (no data-dependent control flow) -- val2 stays in a
                # register and is reused for every permutation
                for p in range(n_perm):
                    rsum[p] += val2 * np.float64(
                        groupings_T[col, p] == groupings_T[row_idx, p]
                    )
            for p in range(n_perm):
                local_sW[p] += rsum[p] * inv_group_sizes[groupings_T[row_idx, p]]

            mirror_row = n - row_idx - 2
            if mirror_row != row_idx:
                for p in range(n_perm):
                    rsum[p] = 0.0
                for col in range(mirror_row + 1, n):
                    val = np.float64(distance_matrix[mirror_row, col])
                    val2 = val * val
                    for p in range(n_perm):
                        rsum[p] += val2 * np.float64(
                            groupings_T[col, p] == groupings_T[mirror_row, p]
                        )
                for p in range(n_perm):
                    local_sW[p] += (
                        rsum[p] * inv_group_sizes[groupings_T[mirror_row, p]]
                    )

    def _run_permanova_rowtile_nb(
        distmat, grouping, group_sizes, s_T, num_groups, sample_size,
        permutations, seed
    ):
        """Run PERMANOVA with the row-tile (single matrix pass) strategy.

        Generates the ``permutations + 1`` groupings in fixed-size chunks and,
        for each chunk, computes s_W via ``_permanova_f_stat_sW_rowtile_nb``,
        which reads the distance matrix only once per chunk (parallelism is over
        matrix rows, not permutations).  This is the cross-thread cache-reuse
        path of scikit-bio-binaries; it is markedly faster than the
        per-permutation kernel on large matrices, which are memory-bound.

        Groupings are transposed to sample-major ``int32`` once per chunk so the
        kernel's permutation loop reads contiguous, low-bandwidth windows; the
        transpose costs a small fraction of the matrix scan.  Permutations are
        drawn in the same RNG order as the sequential Monte Carlo path, so
        p-values are unchanged.

        Parameters
        ----------
        distmat : DistanceMatrix
            Full (non-condensed) distance matrix.
        grouping : np.ndarray, shape (n,)
            Integer group labels (0-indexed).
        group_sizes : np.ndarray, shape (num_groups,)
            Number of objects per group.
        s_T : float
            Total sum of squares, precomputed by the caller.
        num_groups : int
            Number of distinct groups.
        sample_size : int
            Number of objects (n).
        permutations : int
            Number of Monte Carlo permutations.
        seed : int, Generator, or None
            Random seed passed to ``get_rng``.

        Returns
        -------
        stat : float
            Observed pseudo-F statistic.
        p_value : float
            Permutation p-value.

        """
        if permutations < 0:
            raise ValueError(
                "Number of permutations must be greater than or equal to zero."
            )

        rng = get_rng(seed)
        n_total = permutations + 1
        out_sW = np.empty(n_total, np.float64)
        inv_group_sizes = 1.0 / group_sizes.astype(np.float64)
        n_half = sample_size // 2

        # CHUNK is the number of permutations processed per pass over the
        # matrix. It is not a parallelism knob (the kernel parallelises over
        # matrix rows); it sizes each worker's private working set -- its rsum
        # accumulator and the grouping windows it reads -- so it is a fixed
        # constant chosen to keep that set in L1, not scaled with the thread
        # count. It also bounds peak memory to (CHUNK x n) regardless of the
        # permutation count.
        CHUNK = 256
        width = min(CHUNK, n_total)
        # Permutation-major int32 buffer; rng.permutation fills a row, so the
        # writes are sequential. Allocated once and reused across chunks.
        buf = np.empty((width, sample_size), dtype=np.int32)
        # Kernel scratch, allocated once and reused (the kernel zeros its own
        # rows); max width is CHUNK.
        partials = np.empty((n_half, width), np.float64)
        for start in range(0, n_total, CHUNK):
            end = min(start + CHUNK, n_total)
            chunk_size = end - start
            for i in range(start, end):
                if i == 0:
                    buf[0] = grouping
                else:
                    buf[i - start] = rng.permutation(grouping)

            # Transpose to sample-major so the kernel reads each permutation's
            # labels as a contiguous window; buf is permutation-major and buf.T
            # is a non-contiguous view, so this materialises the C-contiguous
            # copy the kernel needs.
            groupings_T = np.ascontiguousarray(buf[:chunk_size].T)
            _permanova_f_stat_sW_rowtile_nb(
                distmat.data, groupings_T, inv_group_sizes, partials
            )
            # sum each row-pair's partial contribution into s_W per permutation
            out_sW[start:end] = partials[:, :chunk_size].sum(axis=0)

        s_A = s_T - out_sW
        f_stats = (s_A / (num_groups - 1)) / (out_sW / (sample_size - num_groups))
        stat = f_stats[0]
        if permutations == 0:
            return stat, np.nan
        p_value = (1.0 + np.sum(f_stats[1:] >= stat)) / (1.0 + permutations)
        return stat, p_value


@params_aliased([("distmat", "distance_matrix", "0.7.0", False)])
def permanova(
    distmat: DistanceMatrix,
    grouping: pd.DataFrame | ArrayLike,
    column: str | None = None,
    permutations: int = 999,
    seed: SeedLike | None = None,
    engine: str | None = None,
) -> pd.Series:
    r"""Test for significant differences between groups using PERMANOVA.

    Permutational Multivariate Analysis of Variance (PERMANOVA) is a non-parametric
    method that tests whether two or more groups of objects (e.g., samples) are
    significantly different based on a categorical factor. It is conceptually
    similar to ANOVA except that it operates on distances between objects via a
    distance matrix, which allows for multivariate analysis. Unlike classical
    Multivariate Analysis of Variance (MANOVA), PERMANOVA makes no assumptions
    about the distribution of the underlying data. As such, rather than computing
    a true `F` statistic based in known distributions of variables, it computes a
    pseudo-`F` statistic whose significance can be assessed by a permutation test.

    The pseudo-`F` statistic is the ratio of between-group variance to within-group
    variance, defined in [1]_ analogously to the `F` statistic in ANOVA:

    .. math::
        F = \frac{{SS}_{between}/(g - 1)}{{SS}_{within}/(n - g)}

    It is computed from the sums of squares :math:`{SS}_{between}` and
    :math:`{SS}_{within}` divided by their corresponding degrees of freedom, where
    :math:`n` is the number of distinct objects and :math:`g` is the number of groups.

    Statistical significance is assessed via a permutation test. Objects in the distance
    matrix are assigned to groups (`grouping`) based on a categorical factor. This
    assignment of groups is permuted a number of times (controlled via `permutations`),
    and a pseudo-`F` statistic is computed for each permutation. Under the null
    hypothesis that the groupings of objects have no effect on the distribution of
    the underlying data, the pseudo-`F` statistics of these permutations should be
    identically distributed for a given distance matrix. The probability of a given
    pseudo-`F` statistic being at least as extreme as an observed one is then the
    proportion of permuted pseudo-`F` statistics (:math:`F^{\pi}`) that are greater
    than or equal to the observed (unpermuted) one (:math:`F`):

    .. math::
        p = \frac{1 + \text{no. of } F^{\pi} \geq F}{1 + \text{no. of permutations}}

    Parameters
    ----------
    distmat : DistanceMatrix
        Distance matrix containing distances between objects (e.g., distances
        between samples of microbial communities).
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
    permutations : int, optional
        Number of permutations to use when assessing statistical
        significance. Must be greater than or equal to zero. If zero,
        statistical significance calculations will be skipped and the `p`-value
        will be ``np.nan``.
    seed : int, Generator or RandomState, optional
        A user-provided random seed or random generator instance. See
        :func:`details <skbio.util.get_rng>`.

        .. versionadded:: 0.6.3
    engine : {"cython", "numba"}, optional
        Compute engine to use. ``"cython"`` (default) uses the Cython
        implementation. ``"numba"`` uses the optional Numba implementation
        and requires Numba to be installed. If not provided, the global
        default is used (see :func:`skbio.set_config`).

    Returns
    -------
    pandas.Series
        Results of the statistical test, including ``test statistic`` and
        ``p-value``.

    See Also
    --------
    anosim
    statsmodels.multivariate.manova.MANOVA
    scipy.stats.permutation_test

    Notes
    -----
    This function uses parallel computation for improved performance.
    See the :install:`parallelization guide <#parallelization>` for information on
    controlling the number of threads used.

    Low-level acceleration is available for this function. See
    :install:`scikit-bio-binaries <#acceleration>` for more information.

    See [1]_ for the original method reference, as well as ``vegan::adonis``,
    available in R's vegan package [2]_.

    The precision of the `p`-value is dependent on the number of permutations. The
    default precision is :math:`0.001=1/(1+999)` from the default value
    ``permutations=999``. The unpermuted grouping always contributes the first
    permutation to the numerator and denominator of the `p`-value, so 1 is added
    to both. This circumvents the risk of the probability being zero by chance
    even when it is nonzero. It is suggested in [1]_ that at least 1000
    permutations should be performed for a confidence level of 0.05, and
    5000 permutations should be performed for a confidence level of 0.01. The
    `p`-value will be ``np.nan`` if ``permutations`` is zero.

    A related statistic reported by some implementations (such as
    ``vegan::adonis``) is the :math:`R^2` value, which describes the proportion
    of variance in the data explained by the grouping:

    .. math::
        R^2 = \frac{{SS}_{between}}{{SS}_{total}}

    This is not currently computed by this function, but it may be derived from the
    outputs using the following formula:

    .. math::
        R^2 = \frac{1}{1 + \frac{n - g}{(g - 1)F}}

    where :math:`F` is the pseudo-`F` statistic, :math:`n` is the number of
    objects, and :math:`g` is the number of groups.

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
    if not isinstance(distmat, DistanceMatrix):
        raise TypeError("Input must be a DistanceMatrix.")

    sample_size = distmat.shape[0]

    num_groups, grouping = _preprocess_input_sng(
        distmat.ids, sample_size, grouping, column
    )

    engine = _resolve_engine(engine, ("cython", "numba"))

    if _skbb_permanova_available(
        distmat, grouping, permutations, seed
    ):  # pragma: no cover
        # unlikely to throw here, but just in case
        try:
            stat, p_value = _skbb_permanova(distmat, grouping, permutations, seed)
            return _build_results(
                "PERMANOVA",
                "pseudo-F",
                sample_size,
                num_groups,
                stat,
                p_value,
                permutations,
            )
        except Exception as e:
            warn(
                "Attempted to use binaries.permanova but failed, "
                "using regular logic instead.",
                RuntimeWarning,
            )
    # if we got here, we could not use skbb
    # Calculate number of objects in each group.
    group_sizes = np.bincount(grouping)
    s_T = (distmat.data**2).sum() / sample_size
    if not distmat._flags["CONDENSED"]:
        # we are going over the whole matrix, instead of just upper triangle
        # so cut in half
        s_T /= 2.0

    if engine == "numba" and not distmat._flags["CONDENSED"]:
        stat, p_value = _run_permanova_rowtile_nb(
            distmat, grouping, group_sizes, s_T,
            num_groups, sample_size, permutations, seed
        )
    else:
        test_stat_function = partial(
            _compute_f_stat, sample_size, num_groups, distmat, group_sizes, s_T, engine
        )
        stat, p_value = _run_monte_carlo_stats(
            test_stat_function, grouping, permutations, seed
        )

    return _build_results(
        "PERMANOVA", "pseudo-F", sample_size, num_groups, stat, p_value, permutations
    )


def _compute_f_stat(
    sample_size, num_groups, distance_matrix, group_sizes, s_T, engine, grouping
):
    """Compute PERMANOVA pseudo-F statistic."""
    # Calculate s_W for each group, accounting for different group sizes.
    if distance_matrix._flags["CONDENSED"]:
        if engine == "numba":
            s_W = _permanova_f_stat_sW_condensed_nb(
                distance_matrix.data, group_sizes, grouping
            )
        else:
            s_W = permanova_f_stat_sW_condensed_cy(
                distance_matrix.data, group_sizes, grouping
            )
    else:
        if engine == "numba":
            s_W = _permanova_f_stat_sW_nb(
                distance_matrix.data, group_sizes, grouping
            )
        else:
            s_W = permanova_f_stat_sW_cy(
                distance_matrix.data, group_sizes, grouping
            )

    s_A = s_T - s_W
    return (s_A / (num_groups - 1)) / (s_W / (sample_size - num_groups))
