# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import functools
import itertools
from warnings import warn

import numpy as np
import scipy.spatial.distance
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
    _quantitative_to_qualitative_counts,
    _table_to_numpy,
    _validate_table,
)
from skbio.util._warning import _warn_deprecated


def _get_alpha_diversity_metric_map():
    return {
        "ace": skbio.diversity.alpha.ace,
        "chao1": skbio.diversity.alpha.chao1,
        "chao1_ci": skbio.diversity.alpha.chao1_ci,
        "berger_parker_d": skbio.diversity.alpha.berger_parker_d,
        "brillouin_d": skbio.diversity.alpha.brillouin_d,
        "dominance": skbio.diversity.alpha.dominance,
        "doubles": skbio.diversity.alpha.doubles,
        "enspie": skbio.diversity.alpha.enspie,
        "esty_ci": skbio.diversity.alpha.esty_ci,
        "faith_pd": skbio.diversity.alpha.faith_pd,
        "fisher_alpha": skbio.diversity.alpha.fisher_alpha,
        "gini_index": skbio.diversity.alpha.gini_index,
        "goods_coverage": skbio.diversity.alpha.goods_coverage,
        "inv_simpson": skbio.diversity.alpha.inv_simpson,
        "hill": skbio.diversity.alpha.hill,
        "heip_e": skbio.diversity.alpha.heip_e,
        "kempton_taylor_q": skbio.diversity.alpha.kempton_taylor_q,
        "lladser_ci": skbio.diversity.alpha.lladser_ci,
        "lladser_pe": skbio.diversity.alpha.lladser_pe,
        "margalef": skbio.diversity.alpha.margalef,
        "mcintosh_d": skbio.diversity.alpha.mcintosh_d,
        "mcintosh_e": skbio.diversity.alpha.mcintosh_e,
        "menhinick": skbio.diversity.alpha.menhinick,
        "michaelis_menten_fit": skbio.diversity.alpha.michaelis_menten_fit,
        "observed_features": skbio.diversity.alpha.observed_features,
        "observed_otus": skbio.diversity.alpha.observed_otus,
        "osd": skbio.diversity.alpha.osd,
        "phydiv": skbio.diversity.alpha.phydiv,
        "pielou_e": skbio.diversity.alpha.pielou_e,
        "renyi": skbio.diversity.alpha.renyi,
        "robbins": skbio.diversity.alpha.robbins,
        "shannon": skbio.diversity.alpha.shannon,
        "simpson": skbio.diversity.alpha.simpson,
        "simpson_d": skbio.diversity.alpha.simpson_d,
        "simpson_e": skbio.diversity.alpha.simpson_e,
        "singles": skbio.diversity.alpha.singles,
        "sobs": skbio.diversity.alpha.sobs,
        "strong": skbio.diversity.alpha.strong,
        "tsallis": skbio.diversity.alpha.tsallis,
    }


def get_alpha_diversity_metrics():
    """List scikit-bio's alpha diversity metrics.

    The alpha diversity metrics listed here can be passed as metrics to
    ``skbio.diversity.alpha_diversity``.

    Returns
    -------
    list of str
        Alphabetically sorted list of alpha diversity metrics implemented in
        scikit-bio.

    See Also
    --------
    alpha_diversity
    get_beta_diversity_metrics

    """
    metrics = _get_alpha_diversity_metric_map()
    return sorted(metrics.keys())


def get_beta_diversity_metrics():
    """List scikit-bio's beta diversity metrics.

    The beta diversity metrics listed here can be passed as metrics to
    ``skbio.diversity.beta_diversity``.

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
    ``scipy.spatial.distance.pdist`` for more details.

    """
    return sorted(["unweighted_unifrac", "weighted_unifrac"])


def alpha_diversity(metric, counts, ids=None, validate=True, **kwargs):
    """Compute alpha diversity for one or more samples.

    Parameters
    ----------
    metric : str, callable
        The alpha diversity metric to apply to the sample(s). Passing metric as
        a string is preferable as this often results in an optimized version of
        the metric being used.
    counts : 1D or 2D array_like of ints or floats, Table
        Vector or matrix containing count/abundance data. If a matrix, each row
        should contain counts of taxa in a given sample.
    ids : iterable of strs, optional
        Identifiers for each sample in ``counts``. By default, samples will be
        assigned integer identifiers in the order that they were provided.
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
    skbio.diversity
    skbio.diversity.alpha
    skbio.diversity.get_alpha_diversity_metrics
    skbio.diversity.beta_diversity

    """
    if isinstance(counts, skbio.Table):
        counts, ids = _validate_table(counts, ids, kwargs)

    metric_map = _get_alpha_diversity_metric_map()

    if validate:
        counts = _validate_counts_matrix(counts, ids=ids)

    if metric == "faith_pd":
        taxa, tree, kwargs = _get_phylogenetic_kwargs(counts, **kwargs)
        counts_by_node, branch_lengths = _setup_pd(
            counts, taxa, tree, validate, rooted=True, single_sample=False
        )
        counts = counts_by_node
        metric = functools.partial(_faith_pd, branch_lengths=branch_lengths, **kwargs)

    elif metric == "phydiv":
        taxa, tree, kwargs = _get_phylogenetic_kwargs(counts, **kwargs)
        counts_by_node, branch_lengths = _setup_pd(
            counts, taxa, tree, validate, rooted=False, single_sample=False
        )
        counts = counts_by_node
        if "rooted" not in kwargs:
            kwargs["rooted"] = len(tree.root().children) == 2
        if "weight" not in kwargs:
            kwargs["weight"] = False
        metric = functools.partial(_phydiv, branch_lengths=branch_lengths, **kwargs)

    elif callable(metric):
        metric = functools.partial(metric, **kwargs)
    elif metric in metric_map:
        metric = functools.partial(metric_map[metric], **kwargs)
    else:
        raise ValueError("Unknown metric provided: %r." % metric)

    # kwargs is provided here so an error is raised on extra kwargs
    results = [metric(c, **kwargs) for c in counts]
    return pd.Series(results, index=ids)


def partial_beta_diversity(metric, counts, ids, id_pairs, validate=True, **kwargs):
    """Compute distances only between specified ID pairs.

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
        See ``skbio.diversity.beta_diversity`` for details.
    kwargs : kwargs, optional
        Metric-specific parameters.

    Returns
    -------
    skbio.DistanceMatrix
        Distances between pairs of samples indicated by id_pairs. Pairwise
        distances not defined by id_pairs will be 0.0. Use this resulting
        DistanceMatrix with caution as 0.0 is a valid distance.

    Warnings
    --------
    ``partial_beta_diversity`` is deprecated as of ``0.5.0``. The return type is
    unstable. Developer caution is advised. The resulting DistanceMatrix object will
    include zeros when distance has not been calculated, and therefore can be
    misleading.

    Raises
    ------
    ValueError
        If ``ids`` are not specified.
        If ``id_pairs`` are not a subset of ``ids``.
        If ``metric`` is not a callable or is unresolvable string by
        scikit-bio.
        If duplicates are observed in ``id_pairs``.

    See Also
    --------
    skbio.diversity.beta_diversity
    skbio.diversity.get_beta_diversity_metrics

    """
    # @deprecated
    _warn_deprecated(
        partial_beta_diversity,
        "0.5.0",
        msg="The return type is unstable. Developer caution is advised. The resulting "
        "DistanceMatrix object will include zeros when distance has not been "
        "calculated, and therefore can be misleading.",
    )

    if validate:
        counts = _validate_counts_matrix(counts, ids=ids)

    id_pairs = list(id_pairs)
    all_ids_in_pairs = set(itertools.chain.from_iterable(id_pairs))
    if not all_ids_in_pairs.issubset(ids):
        raise ValueError("`id_pairs` are not a subset of `ids`")

    hashes = {i for i in id_pairs}.union({i[::-1] for i in id_pairs})
    if len(hashes) != len(id_pairs) * 2:
        raise ValueError("A duplicate or a self-self pair was observed.")

    if metric == "unweighted_unifrac":
        counts = _quantitative_to_qualitative_counts(counts)
        taxa, tree, kwargs = _get_phylogenetic_kwargs(counts, **kwargs)
        metric, counts_by_node = _setup_multiple_unweighted_unifrac(
            counts, taxa=taxa, tree=tree, validate=validate
        )
        counts = counts_by_node
    elif metric == "weighted_unifrac":
        # get the value for normalized. if it was not provided, it will fall
        # back to the default value inside of _weighted_unifrac_pdist_f
        normalized = kwargs.pop("normalized", _normalize_weighted_unifrac_by_default)
        taxa, tree, kwargs = _get_phylogenetic_kwargs(counts, **kwargs)
        metric, counts_by_node = _setup_multiple_weighted_unifrac(
            counts, taxa=taxa, tree=tree, normalized=normalized, validate=validate
        )
        counts = counts_by_node
    elif callable(metric):
        metric = functools.partial(metric, **kwargs)
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


# The following two lists are adapted from sklearn.metrics.pairwise. Metrics
# that are not available in SciPy (only in sklearn) have been removed from
# the list of _valid_beta_metrics here (those are: manhatten, wminkowski,
# nan_euclidean, and haversine)
_valid_beta_metrics = [
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
]


_qualitative_beta_metrics = [
    "dice",
    "jaccard",
    "matching",
    "rogerstanimoto",
    "russellrao",
    "sokalmichener",
    "sokalsneath",
    "yule",
]


def beta_diversity(
    metric, counts, ids=None, validate=True, pairwise_func=None, **kwargs
):
    """Compute distances between all pairs of samples.

    Parameters
    ----------
    metric : str, callable
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
    iTypeError
        If invalid method-specific parameters are provided.

    See Also
    --------
    skbio.diversity
    skbio.diversity.beta
    skbio.diversity.get_beta_diversity_metrics
    skbio.diversity.alpha_diversity
    scipy.spatial.distance.pdist
    sklearn.metrics.pairwise_distances

    """
    if isinstance(counts, skbio.Table):
        counts, ids = _validate_table(counts, ids, kwargs)

    if isinstance(counts, pd.DataFrame) and ids is None:
        ids = list(counts.index)

    if validate:
        counts = _validate_counts_matrix(counts, ids=ids)

    if 0 in counts.shape:
        # if the input counts are empty, return an empty DistanceMatrix.
        # this check is not necessary for scipy.spatial.distance.pdist but
        # it is necessary for sklearn.metrics.pairwise_distances where the
        # latter raises an exception over empty data.
        return DistanceMatrix(np.zeros((len(ids), len(ids))), ids)

    if metric == "unweighted_unifrac":
        counts = _quantitative_to_qualitative_counts(counts)
        taxa, tree, kwargs = _get_phylogenetic_kwargs(counts, **kwargs)
        metric, counts_by_node = _setup_multiple_unweighted_unifrac(
            counts, taxa=taxa, tree=tree, validate=validate
        )
        counts = counts_by_node
    elif metric == "weighted_unifrac":
        # get the value for normalized. if it was not provided, it will fall
        # back to the default value inside of _weighted_unifrac_pdist_f
        normalized = kwargs.pop("normalized", _normalize_weighted_unifrac_by_default)
        taxa, tree, kwargs = _get_phylogenetic_kwargs(counts, **kwargs)
        metric, counts_by_node = _setup_multiple_weighted_unifrac(
            counts, taxa=taxa, tree=tree, normalized=normalized, validate=validate
        )
        counts = counts_by_node
    elif metric == "manhattan":
        metric = "cityblock"
    elif callable(metric):
        metric = functools.partial(metric, **kwargs)
        # remove all values from kwargs, since they have already been provided
        # through the partial
        kwargs = {}
    elif metric in _qualitative_beta_metrics:
        counts = _quantitative_to_qualitative_counts(counts)
    elif metric not in _valid_beta_metrics:
        raise ValueError(
            "Metric %s is not available. "
            "Only the following metrics can be passed as strings to "
            "beta_diversity as we know whether each of these should be "
            "treated as a qualitative or quantitative metric. Other metrics "
            "can be provided as functions.\n Available metrics are: %s"
            % (metric, ", ".join(_valid_beta_metrics))
        )
    else:
        # metric is a string that scikit-bio doesn't know about, for
        # example one of the SciPy metrics
        pass

    if pairwise_func is None:
        pairwise_func = scipy.spatial.distance.pdist

    distances = pairwise_func(counts, metric=metric, **kwargs)
    return DistanceMatrix(distances, ids)
