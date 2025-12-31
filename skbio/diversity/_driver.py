# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import annotations

from functools import partial
from itertools import chain
from inspect import signature, getmembers, isfunction
from typing import Any, TYPE_CHECKING

import numpy as np
import pandas as pd

from skbio.diversity import alpha
from skbio.diversity.alpha._pd import _setup_pd, _faith_pd, _phydiv
from skbio.diversity.beta._unifrac import (
    _setup_multiple_unweighted_unifrac,
    _setup_multiple_weighted_unifrac,
    _normalize_weighted_unifrac_by_default,
)
from skbio.stats.distance import DistanceMatrix
from skbio.diversity._util import (
    _validate_counts_matrix,
    _get_phylogenetic_kwargs,
    _qualify_counts,
)
from skbio.util._decorator import deprecated
from skbio.table._tabular import _ingest_table

if TYPE_CHECKING:  # pragma: no cover
    from collections.abc import Iterable, Callable
    from numpy.typing import ArrayLike
    from skbio.util._typing import TableLike


# Qualitative (absence/presence-based) metrics, which require Boolean input.
_qualitative_metrics = {
    "dice",
    "jaccard",
    "matching",
    "rogerstanimoto",
    "russellrao",
    "sokalmichener",
    "sokalsneath",
    "yule",
    "unweighted_unifrac",
}

# Beta diversity metrics implemented in SciPy's pdist.
_pdist_metrics = {
    "euclidean",
    "cityblock",
    "braycurtis",
    "canberra",
    "chebyshev",
    "correlation",
    "cosine",
    "dice",
    "hamming",
    "jaccard",
    "jensenshannon",
    "mahalanobis",
    "manhattan",  # aliases to "cityblock" in beta_diversity
    "matching",
    "minkowski",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalmichener",
    "sokalsneath",
    "sqeuclidean",
    "yule",
}


def get_alpha_diversity_metrics() -> list[str]:
    r"""List scikit-bio's alpha diversity metrics.

    The alpha diversity metrics listed here can be passed as metrics to
    :func:`alpha_diversity`.

    Returns
    -------
    list of str
        Alphabetically sorted list of alpha diversity metrics implemented in
        scikit-bio.

    See Also
    --------
    alpha
    alpha_diversity
    get_beta_diversity_metrics

    """
    return [x[0] for x in getmembers(alpha, isfunction)]


def get_beta_diversity_metrics() -> list[str]:
    """List scikit-bio's beta diversity metrics.

    The beta diversity metrics listed here can be passed as metrics to
    :func:`beta_diversity`.

    Returns
    -------
    list of str
        Alphabetically sorted list of beta diversity metrics implemented in
        scikit-bio.

    See Also
    --------
    beta_diversity
    get_alpha_diversity_metrics
    scipy.spatial.distance.pdist

    Notes
    -----
    SciPy implements many additional beta diversity metrics that are not
    included in this list. See documentation for
    SciPy's :func:`~scipy.spatial.distance.pdist` for more details.

    """
    return sorted(_pdist_metrics.union(["unweighted_unifrac", "weighted_unifrac"]))


def alpha_diversity(
    metric: str | Callable,
    counts: TableLike,
    ids: ArrayLike | None = None,
    validate: bool = True,
    **kwargs: Any,
) -> pd.Series:
    r"""Compute alpha diversity for one or more samples.

    Parameters
    ----------
    metric : str or callable
        The alpha diversity metric to apply to the sample(s). See
        :mod:`~skbio.diversity.alpha` for available metrics. Passing metric as a string
        is preferable as this often results in an optimized version of the metric
        being used.
    counts : table_like of shape (n_samples, n_taxa) or (n_taxa,)
        Vector or matrix containing count/abundance data of one or multiple samples.
        See :ref:`supported formats <table_like>`.
    ids : array_like of shape (n_samples,), optional
        Identifiers for each sample in ``counts``. If not provided, will extract sample
        IDs from ``counts``, if available, or assign integer identifiers in the order
        that samples were provided.
    validate: bool, optional
        If True (default), validate the input data before applying the alpha diversity
        metric. See :mod:`skbio.diversity` for the details of validation.
    kwargs : dict, optional
        Metric-specific parameters. Refer to the documentation of the chosen metric.
        A special parameter is ``taxa``, needed by some phylogenetic metrics. If not
        provided, will extract taxa (feature IDs) from ``counts``, if available, and
        pass to the metric.

    Returns
    -------
    pd.Series
        Values of ``metric`` for all samples provided in ``counts``. The index
        will be ``ids``, if provided.

    Raises
    ------
    ValueError, MissingNodeError, DuplicateNodeError
        If validation fails. Exact error will depend on what was invalid.
    TypeError
        If invalid method-specific parameters are provided.

    See Also
    --------
    alpha
    get_alpha_diversity_metrics
    beta_diversity

    """
    taxa = kwargs.pop("taxa", None)
    counts, ids, taxa = _ingest_table(counts, ids, taxa)
    if validate:
        counts = _validate_counts_matrix(counts)

    # phylogenetic diversity (call on entire matrix)
    if metric in ("faith_pd", "phydiv"):
        taxa, tree, kwargs = _get_phylogenetic_kwargs(kwargs, taxa)
        is_faith_pd = metric == "faith_pd"
        counts, lengths = _setup_pd(
            counts, taxa, tree, validate, rooted=is_faith_pd, single_sample=False
        )
        if is_faith_pd:
            func = _faith_pd
        else:
            func = _phydiv
            kwargs.setdefault("rooted", tree._is_rooted())
            kwargs.setdefault("weight", False)
        func = partial(func, branch_lengths=lengths, **kwargs)

    # other metrics (call on each sample)
    else:
        if callable(metric):
            func = metric
        elif isinstance(metric, str):
            func = getattr(alpha, metric, None)
            if func is None or not isfunction(func):
                raise ValueError(
                    f'"{metric}" is not an available alpha diversity metric name. '
                    "Refer to `get_alpha_diversity_metrics` for a list of available "
                    "metrics."
                )
        else:
            raise ValueError(f"Invalid metric provided: {metric!r}.")

        # add "taxa" back to parameters
        if taxa is not None and "taxa" in signature(func).parameters:
            kwargs["taxa"] = taxa
        func = partial(func, **kwargs)

    # kwargs is provided here so an error is raised on extra kwargs
    results = [func(c, **kwargs) for c in counts]
    return pd.Series(results, index=ids)


def beta_diversity(
    metric: str | Callable,
    counts: TableLike,
    ids: ArrayLike | None = None,
    validate: bool = True,
    pairwise_func: Callable | None = None,
    **kwargs: Any,
) -> DistanceMatrix:
    r"""Compute distances between all pairs of samples.

    Parameters
    ----------
    metric : str or callable
        The beta diversity metric, i.e., a pairwise distance function to apply to the
        sample(s). See :mod:`~skbio.diversity.beta` and SciPy's
        :func:`~scipy.spatial.distance.pdist` for available metrics. Passing metric as
        a string is preferable as this often results in an optimized version of the
        metric being used.
    counts : table_like of shape (n_samples, n_taxa) or (n_taxa,)
        Vector or matrix containing count/abundance data of one or multiple samples.
        See :ref:`supported formats <table_like>`.
    ids : array_like of shape (n_samples,), optional
        Identifiers for each sample in ``counts``. If not provided, will extract sample
        IDs from ``counts``, if available, or assign integer identifiers in the order
        that samples were provided.
    validate: bool, optional
        If True (default), validate the input data before applying the alpha diversity
        metric. See :mod:`skbio.diversity` for the details of validation.
    pairwise_func : callable, optional
        The function to use for computing pairwise distances. Must take ``counts`` and
        ``metric`` and return a square, hollow, 2-D float array of dissimilarities.
        Examples of functions that can be provided are SciPy's
        :func:`~scipy.spatial.distance.pdist` (default) and scikit-learn's
        :func:`~sklearn.metrics.pairwise_distances`.
    kwargs : dict, optional
        Metric-specific parameters. Refer to the documentation of the chosen metric.
        A special parameter is ``taxa``, needed by some phylogenetic metrics. If not
        provided, will extract taxa (feature IDs) from ``counts``, if available, and
        pass to the metric.

    Returns
    -------
    :class:`~skbio.stats.distance.DistanceMatrix`
        Distances between all pairs of samples (i.e., rows). The number of
        rows and columns will be equal to the number of rows in ``counts``.

    Raises
    ------
    ValueError, MissingNodeError, DuplicateNodeError
        If validation fails. Exact error will depend on what was invalid.
    Any Exception
        If invalid method-specific parameters are provided.

    See Also
    --------
    beta
    get_beta_diversity_metrics
    alpha_diversity
    scipy.spatial.distance.pdist
    sklearn.metrics.pairwise_distances

    """
    taxa = kwargs.pop("taxa", None)
    counts, ids, taxa = _ingest_table(counts, ids, taxa)
    if 0 in counts.shape:
        # if the input counts are empty, return an empty DistanceMatrix.
        # this check is not necessary for scipy.spatial.distance.pdist but
        # it is necessary for sklearn.metrics.pairwise_distances where the
        # latter raises an exception over empty data.
        return DistanceMatrix(np.zeros((len(ids), len(ids))), ids)
    if validate:
        counts = _validate_counts_matrix(counts)
    if metric in _qualitative_metrics:
        counts = _qualify_counts(counts)
    if metric in ("unweighted_unifrac", "weighted_unifrac"):
        taxa, tree, kwargs = _get_phylogenetic_kwargs(kwargs, taxa)

    if metric == "unweighted_unifrac":
        metric, counts = _setup_multiple_unweighted_unifrac(
            counts, taxa=taxa, tree=tree, validate=validate
        )
    elif metric == "weighted_unifrac":
        # get the value for normalized. if it was not provided, it will fall
        # back to the default value inside of _weighted_unifrac_pdist_f
        normalized = kwargs.pop("normalized", _normalize_weighted_unifrac_by_default)
        metric, counts = _setup_multiple_weighted_unifrac(
            counts, taxa=taxa, tree=tree, normalized=normalized, validate=validate
        )
    elif metric == "manhattan":
        metric = "cityblock"
    elif metric == "mahalanobis":
        nrow, ncol = counts.shape
        if nrow < ncol:
            raise ValueError(
                "Metric 'mahalanobis' requires more samples than features. "
                f"The input has {nrow} samples and {ncol} features."
            )
    elif callable(metric):
        # add "taxa" back to parameters
        if taxa is not None and "taxa" in signature(metric).parameters:
            kwargs["taxa"] = taxa
        metric = partial(metric, **kwargs)
        # remove all values from kwargs, since they have already been provided
        # through the partial
        kwargs = {}
    elif metric not in _pdist_metrics:
        raise ValueError(
            f'"{metric}" is not an available beta diversity metric name. '
            "Refer to `get_beta_diversity_metrics` for a list of available metrics."
        )

    if pairwise_func is None:
        from scipy.spatial.distance import pdist

        pairwise_func = pdist

    distances = pairwise_func(counts, metric=metric, **kwargs)
    return DistanceMatrix(distances, ids)


@deprecated(
    "0.5.0",
    msg="The return type is unstable. Developer caution is advised. The resulting "
    "DistanceMatrix object will include zeros when distance has not been calculated, "
    "and therefore can be misleading.",
)
def partial_beta_diversity(
    metric: str | Callable,
    counts: TableLike,
    ids: ArrayLike | None,
    id_pairs: Iterable[tuple[str, str]],
    validate: bool = True,
    **kwargs: Any,
) -> DistanceMatrix:
    r"""Compute distances only between specified ID pairs.

    Parameters
    ----------
    metric : str or callable
        The beta diversity metric to apply to the samples. See :func:`beta_diversity`
        for details.
    counts : table_like of shape (n_samples, n_taxa)
        Matrix containing count/abundance data of the samples. See
        :ref:`supported formats <table_like>`.
    ids : iterable of strs
        Identifiers for each sample in ``counts``.
    id_pairs : iterable of tuple of (str, str)
        Pairs of sample IDs to compare (e.g., ``[('a', 'b'), ('a', 'c'), ...])``. If
        specified, the they must be a subset of ``ids``.
    validate : bool, optional
        Validate the input data. See ``beta_diversity`` for details.
    kwargs : dict, optional
        Metric-specific parameters. See ``beta_diversity`` for details.

    Returns
    -------
    :class:`~skbio.stats.distance.DistanceMatrix`
        Distances between pairs of samples indicated by ``id_pairs``. Pairwise
        distances not defined by id_pairs will be 0.0. Use this resulting
        DistanceMatrix with caution as 0.0 is a valid distance.

    Raises
    ------
    ValueError
        If ``ids`` are not specified.
        If ``id_pairs`` are not a subset of ``ids``.
        If ``metric`` is not a callable or is unresolvable string by scikit-bio.
        If duplicates are observed in ``id_pairs``.

    See Also
    --------
    beta_diversity
    get_beta_diversity_metrics

    """
    taxa = kwargs.pop("taxa", None)
    counts, ids, taxa = _ingest_table(counts, ids, taxa)
    if validate:
        counts = _validate_counts_matrix(counts)
    if metric in _qualitative_metrics:
        counts = _qualify_counts(counts)

    id_pairs = list(id_pairs)
    all_ids_in_pairs = set(chain.from_iterable(id_pairs))
    if not all_ids_in_pairs.issubset(ids):
        raise ValueError("`id_pairs` are not a subset of `ids`")

    hashes = {i for i in id_pairs}.union({i[::-1] for i in id_pairs})
    if len(hashes) != len(id_pairs) * 2:
        raise ValueError("A duplicate or a self-self pair was observed.")

    if metric == "unweighted_unifrac":
        taxa, tree, kwargs = _get_phylogenetic_kwargs(kwargs, taxa)
        func, counts = _setup_multiple_unweighted_unifrac(
            counts, taxa=taxa, tree=tree, validate=validate
        )
    elif metric == "weighted_unifrac":
        # get the value for normalized. if it was not provided, it will fall
        # back to the default value inside of _weighted_unifrac_pdist_f
        normalized = kwargs.pop("normalized", _normalize_weighted_unifrac_by_default)
        taxa, tree, kwargs = _get_phylogenetic_kwargs(kwargs, taxa)
        func, counts = _setup_multiple_weighted_unifrac(
            counts, taxa=taxa, tree=tree, normalized=normalized, validate=validate
        )
    elif callable(metric):
        func = partial(metric, **kwargs)
        # remove all values from kwargs, since they have already been provided
        # through the partial
        kwargs = {}
    else:
        raise ValueError(
            "partial_beta_diversity is only compatible with "
            "optimized unifrac methods and callable functions."
        )

    dm = np.zeros((len(ids), len(ids)))
    id_index = {id_: idx for idx, id_ in enumerate(ids)}
    id_pairs_indexed = ((id_index[u], id_index[v]) for u, v in id_pairs)

    for u, v in id_pairs_indexed:
        dm[u, v] = func(counts[u], counts[v], **kwargs)

    return DistanceMatrix(dm + dm.T, ids)
