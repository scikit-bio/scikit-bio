# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import partial
from itertools import chain
from inspect import signature, getmembers, isfunction

import numpy as np
import pandas as pd

import skbio
from skbio.diversity.alpha._pd import _faith_pd, _phydiv, _setup_pd
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


def get_alpha_diversity_metrics():
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
    return [x[0] for x in getmembers(skbio.diversity.alpha, isfunction)]


def get_beta_diversity_metrics():
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


def alpha_diversity(metric, counts, ids=None, validate=True, **kwargs):
    r"""Compute alpha diversity for one or more samples.

    Parameters
    ----------
    metric : str or callable
        The alpha diversity metric to apply to the sample(s). Passing metric as
        a string is preferable as this often results in an optimized version of
        the metric being used.
    counts : 1D or 2D array_like of ints or floats, Table
        Vector or matrix containing count/abundance data. If a matrix, each row
        should contain counts of taxa in a given sample.
    ids : iterable of strs, optional
        Identifiers for each sample in ``counts``. If not provided, will extract sample
        IDs from ``counts``, if available, or assign integer identifiers in the order
        that samples were provided.
    validate: bool, optional
        If ``False``, validation of the input won't be performed. This step can
        be slow, so if validation is run elsewhere it can be disabled here.
        However, invalid input data can lead to invalid results or error
        messages that are hard to interpret, so this step should not be
        bypassed if you're not certain that your input data are valid. See
        :mod:`skbio.diversity` for the description of what validation entails
        so you can determine if you can safely disable validation.
    kwargs : kwargs, optional
        Metric-specific parameters.

    Returns
    -------
    pd.Series
        Values of ``metric`` for all vectors provided in ``counts``. The index
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
        metric = partial(func, branch_lengths=lengths, **kwargs)

    # other metrics (call on each sample)
    else:
        if isinstance(metric, str):
            metric = getattr(skbio.diversity.alpha, metric, None)
            if metric is None or not isfunction(metric):
                raise ValueError(
                    f'"{metric}" is not an available alpha diversity metric name. '
                    "Refer to `get_alpha_diversity_metrics` for a list of available "
                    "metrics."
                )
        elif not callable(metric):
            raise ValueError(f"Invalid metric provided: {metric!r}.")

        # add "taxa" back to parameters
        if taxa is not None and "taxa" in signature(metric).parameters:
            kwargs["taxa"] = taxa
        metric = partial(metric, **kwargs)

    # kwargs is provided here so an error is raised on extra kwargs
    results = [metric(c, **kwargs) for c in counts]
    return pd.Series(results, index=ids)


@deprecated(
    "0.5.0",
    msg="The return type is unstable. Developer caution is advised. The resulting "
    "DistanceMatrix object will include zeros when distance has not been calculated, "
    "and therefore can be misleading.",
)
def partial_beta_diversity(metric, counts, ids, id_pairs, validate=True, **kwargs):
    r"""Compute distances only between specified ID pairs.

    Parameters
    ----------
    metric : str or callable
        The pairwise distance function to apply. If ``metric`` is a string, it
        must be resolvable by scikit-bio (e.g., UniFrac methods), or must be
        callable.
    counts : 2D array_like of ints or floats
        Matrix containing count/abundance data where each row contains counts
        of taxa in a given sample.
    ids : iterable of strs
        Identifiers for each sample in ``counts``.
    id_pairs : iterable of tuple
        An iterable of tuples of IDs to compare (e.g., ``[('a', 'b'), ('a',
        'c'), ...])``. If specified, the set of IDs described must be a subset
        of ``ids``.
    validate : bool, optional
        See :func:`beta_diversity` for details.
    kwargs : kwargs, optional
        Metric-specific parameters.

    Returns
    -------
    skbio.DistanceMatrix
        Distances between pairs of samples indicated by id_pairs. Pairwise
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
        metric, counts = _setup_multiple_unweighted_unifrac(
            counts, taxa=taxa, tree=tree, validate=validate
        )
    elif metric == "weighted_unifrac":
        # get the value for normalized. if it was not provided, it will fall
        # back to the default value inside of _weighted_unifrac_pdist_f
        normalized = kwargs.pop("normalized", _normalize_weighted_unifrac_by_default)
        taxa, tree, kwargs = _get_phylogenetic_kwargs(kwargs, taxa)
        metric, counts = _setup_multiple_weighted_unifrac(
            counts, taxa=taxa, tree=tree, normalized=normalized, validate=validate
        )
    elif callable(metric):
        metric = partial(metric, **kwargs)
        # remove all values from kwargs, since they have already been provided
        # through the partial
        kwargs = {}
    else:
        raise ValueError(
            "partial_beta_diversity is only compatible with "
            "optimized unifrac methods and callable functions."
        )

    dm = np.zeros((len(ids), len(ids)), dtype=float)
    id_index = {id_: idx for idx, id_ in enumerate(ids)}
    id_pairs_indexed = ((id_index[u], id_index[v]) for u, v in id_pairs)

    for u, v in id_pairs_indexed:
        dm[u, v] = metric(counts[u], counts[v], **kwargs)

    return DistanceMatrix(dm + dm.T, ids)


def beta_diversity(
    metric, counts, ids=None, validate=True, pairwise_func=None, **kwargs
):
    r"""Compute distances between all pairs of samples.

    Parameters
    ----------
    metric : str or callable
        The pairwise distance function to apply. See the scipy ``pdist`` docs
        and the scikit-bio functions linked under *See Also* for available
        metrics. Passing metrics as a strings is preferable as this often
        results in an optimized version of the metric being used.
    counts : 2D array_like of ints or floats, 2D pandas DataFrame, Table
        Matrix containing count/abundance data where each row contains counts
        of taxa in a given sample.
    ids : iterable of strs, optional
        Identifiers for each sample in ``counts``. By default, samples will be
        assigned integer identifiers in the order that they were provided
        (where the type of the identifiers will be ``str``).
    validate : bool, optional
        If ``False``, validation of the input won't be performed. This step can
        be slow, so if validation is run elsewhere it can be disabled here.
        However, invalid input data can lead to invalid results or error
        messages that are hard to interpret, so this step should not be
        bypassed if you're not certain that your input data are valid. See
        :mod:`skbio.diversity` for the description of what validation entails
        so you can determine if you can safely disable validation.
    pairwise_func : callable, optional
        The function to use for computing pairwise distances. This function
        must take ``counts`` and ``metric`` and return a square, hollow, 2-D
        ``numpy.ndarray`` of dissimilarities (floats). Examples of functions
        that can be provided are ``scipy.spatial.distance.pdist`` and
        ``sklearn.metrics.pairwise_distances``. By default,
        ``scipy.spatial.distance.pdist`` will be used.
    kwargs : kwargs, optional
        Metric-specific parameters.

    Returns
    -------
    skbio.DistanceMatrix
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
