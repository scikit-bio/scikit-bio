# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import functools
import itertools
import operator

import numpy as np
import scipy.spatial.distance
import pandas as pd

import skbio
from skbio.diversity.alpha._faith_pd import _faith_pd, _setup_faith_pd
from skbio.diversity.beta._unifrac import (
    _setup_multiple_unweighted_unifrac, _setup_multiple_weighted_unifrac,
    _normalize_weighted_unifrac_by_default)
from skbio.util._decorator import experimental
from skbio.stats.distance import DistanceMatrix
from skbio.diversity._util import (_validate_counts_matrix,
                                   _get_phylogenetic_kwargs)


def _get_alpha_diversity_metric_map():
    return {
        'ace': skbio.diversity.alpha.ace,
        'chao1': skbio.diversity.alpha.chao1,
        'chao1_ci': skbio.diversity.alpha.chao1_ci,
        'berger_parker_d': skbio.diversity.alpha.berger_parker_d,
        'brillouin_d': skbio.diversity.alpha.brillouin_d,
        'dominance': skbio.diversity.alpha.dominance,
        'doubles': skbio.diversity.alpha.doubles,
        'enspie': skbio.diversity.alpha.enspie,
        'esty_ci': skbio.diversity.alpha.esty_ci,
        'faith_pd': skbio.diversity.alpha.faith_pd,
        'fisher_alpha': skbio.diversity.alpha.fisher_alpha,
        'goods_coverage': skbio.diversity.alpha.goods_coverage,
        'heip_e': skbio.diversity.alpha.heip_e,
        'kempton_taylor_q': skbio.diversity.alpha.kempton_taylor_q,
        'margalef': skbio.diversity.alpha.margalef,
        'mcintosh_d': skbio.diversity.alpha.mcintosh_d,
        'mcintosh_e': skbio.diversity.alpha.mcintosh_e,
        'menhinick': skbio.diversity.alpha.menhinick,
        'michaelis_menten_fit': skbio.diversity.alpha.michaelis_menten_fit,
        'observed_otus': skbio.diversity.alpha.observed_otus,
        'osd': skbio.diversity.alpha.osd,
        'pielou_e': skbio.diversity.alpha.pielou_e,
        'robbins': skbio.diversity.alpha.robbins,
        'shannon': skbio.diversity.alpha.shannon,
        'simpson': skbio.diversity.alpha.simpson,
        'simpson_e': skbio.diversity.alpha.simpson_e,
        'singles': skbio.diversity.alpha.singles,
        'strong': skbio.diversity.alpha.strong,
        'gini_index': skbio.diversity.alpha.gini_index,
        'lladser_pe': skbio.diversity.alpha.lladser_pe,
        'lladser_ci': skbio.diversity.alpha.lladser_ci}


@experimental(as_of="0.4.1")
def get_alpha_diversity_metrics():
    """ List scikit-bio's alpha diversity metrics

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


@experimental(as_of="0.4.1")
def get_beta_diversity_metrics():
    """ List scikit-bio's beta diversity metrics

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
    return sorted(['unweighted_unifrac', 'weighted_unifrac'])


@experimental(as_of="0.4.1")
def alpha_diversity(metric, counts, ids=None, validate=True, **kwargs):
    """ Compute alpha diversity for one or more samples

    Parameters
    ----------
    metric : str, callable
        The alpha diversity metric to apply to the sample(s). Passing metric as
        a string is preferable as this often results in an optimized version of
        the metric being used.
    counts : 1D or 2D array_like of ints or floats
        Vector or matrix containing count/abundance data. If a matrix, each row
        should contain counts of OTUs in a given sample.
    ids : iterable of strs, optional
        Identifiers for each sample in ``counts``. By default, samples will be
        assigned integer identifiers in the order that they were provided.
    validate: bool, optional
        If `False`, validation of the input won't be performed. This step can
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
    metric_map = _get_alpha_diversity_metric_map()

    if validate:
        counts = _validate_counts_matrix(counts, ids=ids)

    if metric == 'faith_pd':
        otu_ids, tree, kwargs = _get_phylogenetic_kwargs(counts, **kwargs)
        counts_by_node, branch_lengths = _setup_faith_pd(
            counts, otu_ids, tree, validate, single_sample=False)
        counts = counts_by_node
        metric = functools.partial(_faith_pd, branch_lengths=branch_lengths)
    elif callable(metric):
        metric = functools.partial(metric, **kwargs)
    elif metric in metric_map:
        metric = functools.partial(metric_map[metric], **kwargs)
    else:
        raise ValueError('Unknown metric provided: %r.' % metric)

    # kwargs is provided here so an error is raised on extra kwargs
    results = [metric(c, **kwargs) for c in counts]
    return pd.Series(results, index=ids)


def _partial_pw(ids, id_pairs, counts, metric, **kwargs):
    """Compute distances only between specified ID pairs

    Parameters
    ----------
    metric : callable
        The pairwise distance function to apply.
    counts : 2D array_like of ints or floats
        Matrix containing count/abundance data where each row contains counts
        of OTUs in a given sample.
    ids : iterable of strs
        Identifiers for each sample in ``counts``.
    id_pairs : iterable of tuple, optional
        An iterable of tuples of IDs to compare (e.g., ``[('a', 'b'), ('a',
        'c'), ...])``. If specified, the set of IDs described must be a subset
        of ``ids``. ``id_pairs`` is a mutually exclusive argument with
        ``pairwise_func``. ``metric`` must be resolvable by scikit-bio (e.g.,
        UniFrac methods), or must be callable.
    kwargs : kwargs, optional
        Metric-specific parameters.

    Returns
    -------
    np.array
        A square, hollow, matrix of distances between ID pairs.

    Raises
    ------
    ValueError
        If ``ids`` are not specified.
        If ``id_pairs`` are not a subset of ``ids``.
        If ``metric`` is not a callable.
        If duplicates are observed in ``id_pairs``.

    See Also
    --------
    skbio.diversity.get_beta_diversity_metrics

    """
    if ids is None:
        raise ValueError("`ids` must be specified if `id_pairs` is specified")

    all_ids_in_pairs = set(itertools.chain.from_iterable(id_pairs))
    if not all_ids_in_pairs.issubset(ids):
        raise ValueError("`id_pairs` are not a subset of `ids`")

    hashes = {i for i in id_pairs}.union({i[::-1] for i in id_pairs})
    if len(hashes) != len(id_pairs) * 2:
        raise ValueError("Duplicate ID pairs observed.")

    if isinstance(metric, str):
        # The mechanism for going from str to callable in scipy's pdist is
        # not exposed
        raise ValueError("`metric` must be callable")

    dm = np.zeros((len(ids), len(ids)), dtype=float)
    id_index = {id_: idx for idx, id_ in enumerate(ids)}
    id_pairs_indexed = ((id_index[u], id_index[v]) for u, v in id_pairs)

    for u, v in id_pairs_indexed:
        dm[u, v] = metric(counts[u], counts[v], **kwargs)

    return dm + dm.T


def _generate_id_blocks(ids, k=4):
    """Generate blocks of IDs that map into a DistanceMatrix

    Parameters
    ----------
    ids : Iterable
        An iterable of IDs of whatever type.
    k : int
        The size of a block to generate IDs for

    Notes
    -----

    This method is intended to facilitate partial beta diversity calculations.
    Blocks of IDs are generated from the upper triangle of the subsequent
    distance matrix. For instance, given the following distance matrix with
    IDs {A, B, C, D, E}:

      A B C D E
    A 0 # # # #
    B # 0 # # #
    C # # 0 # #
    D # # # 0 #
    E # # # # 0

    The goal of this method is to generate tuples of IDs of size k over the
    upper triangle which correspond to blocks of the matrix to compute. Given
    a k=3, the following ID tuples would be generated:

    ((A, B, C), (A, B, C))
    ((A, B, C), (D, E))
    ((D, E), (D, E))

    This method is not responsible for describing which specific pairs of IDs
    are to be computed, only the subset of the matrix of interest.

    Returns
    -------
    tuple of 1D np.array
        Index 0 contans the row IDs, and index 1 contains the column IDs
    """
    n = len(ids)
    ids_idx = np.arange(n)
    for row_start in range(0, n, k):
        for col_start in range(row_start, n, k):
            row_ids = ids_idx[row_start:row_start + k]
            col_ids = ids_idx[col_start:col_start + k]

            yield (row_ids, col_ids)


def _block_party(counts, row_ids=None, col_ids=None, **kwargs):
    """Subset counts to relevant rows and columns

    Parameters
    ----------
    counts : 2D array_like of ints or floats
        Matrix containing count/abundance data where each row contains counts
        of OTUs in a given sample.
    row_ids : 1D np.ndarray of int
        Block row IDs to keep in the counts matrix.
    col_ids : 1D np.ndarray of int
        Block column IDs to keep in the counts matrix. Note, these correspond
        to rows in the counts matrix, but columns in a subsequent distance
        matrix.

    Returns
    -------
    2D np.ndarray
        The subset counts block containing only the rows which exist in row_ids
        and col_ids, and does not contain any zero'd columns.
    dict
        kwargs with any relevant filtering (e.g., filtering a phylogenetic tree
        to only reflect the OTUs present in the counts matrix).
    """
    ids_to_keep = np.unique(np.hstack([row_ids, col_ids]))

    # create a view of the relevant samples
    counts_block = counts[ids_to_keep]

    # remove from the block any empty observations
    # NOTE: this will perform an implicit copy
    nonzero_cols = (counts_block != 0).any(axis=0)
    counts_block = counts_block[:, nonzero_cols]

    block_kwargs = kwargs.copy()
    if 'tree' in kwargs and 'otu_ids' in kwargs:
        block_kwargs['otu_ids'] = np.asarray(kwargs['otu_ids'])[nonzero_cols]
        block_kwargs['tree'] = kwargs['tree'].shear(block_kwargs['otu_ids'])

    block_kwargs['ids'] = ids_to_keep

    return counts_block, block_kwargs


def _pairs_to_compute(rids, cids):
    """Determine the pairs of samples to compute distances between

    Parameters
    ----------
    rids : Iterable
        The row IDs in the partial pairwise computation.

    cids : Iterable
        The column IDs in the partial pairwise computation.

    Notes
    -----
    The diagonal and any pairs in the lower triangle of the distance matrix
    are ignored.

    Returns
    -------
    list of tuple
        The ID pairs to compute distances between.
    """
    if (rids == cids).all():
        return [(i, j) for idx, i in enumerate(rids) for j in rids[idx+1:]]
    else:
        if set(rids).intersection(set(cids)):
            raise ValueError("Attempting to compute a lower triangle")
        return [(i, j) for i in rids for j in cids if i != j]


def _block_and_id_pairs(ids, k):
    """Generate pairs of IDs to compute distances between"""
    for row_ids, col_ids in _generate_id_blocks(ids, k):
        id_pairs = _pairs_to_compute(row_ids, col_ids)
        yield row_ids, col_ids, id_pairs


def _block_compute(metric, counts, **kwargs):
    """Compute a block within the resulting distance matrix

    Parameters
    ----------
    metric : function or string
        The distance metric to use

    counts : np.ndarray
        The counts array to compute over

    Returns
    -------
    PartialDistanceMatrix
        TODO: Unsure how this works. `beta_diversity` will return an
        incomplete distance matrix corresponding to just the partial distances
        with in the block to compute
    """
    block_counts, block_kw = _block_party(counts, **kwargs)

    return beta_diversity(metric, block_counts, **block_kw)


def _block_computable(*args, **kwargs):
    """Construct a function to compute a block in the resulting matrix

    Returns
    -------
    function
        An executable function which returns a `PartialDistanceMatrix`.
    """
    ids = kwargs.get('ids')
    blocksize = kwargs.get('blocksize')
    for row_ids, col_ids, id_pairs in _block_and_id_pairs(ids, blocksize):
        kw = kwargs.copy()
        kw['id_pairs'] = id_pairs
        kw['row_ids'] = row_ids
        kw['col_ids'] = col_ids

        yield functools.partial(_block_compute(*args, **kw))


def _map_block_decomposition(iterable):
    for func in iterable:
        yield func()


def _reduce_block_decomposition(iterable):
    return functools.reduce(operator.add, list(iterable))


def block_beta_diversity(*args, **kwargs):
    """Perform a block-decomposition beta diversity calculation

    Parameters
    ----------
    reduce_f : function, optional
        A method to reduce `PartialDistanceMatrix` objects into a single
        `DistanceMatrix`. The expected signature is:

        `DistanceMatrix <- f(Iterable of PartialDistanceMatrix`

    map_f: function, optional
        A method that accepts a `_computable`. The expected signature is:
        `PartialDistanceMatrix <- f(Iterable of _computable)`
    """
    reduce_f = kwargs.get('reduce_f', _reduce_block_decomposition)
    map_f = kwargs.get('map_f', _map_block_decomposition)

    return reduce_f(map_f(_block_computable(*args, **kwargs)))


@experimental(as_of="0.4.0")
def beta_diversity(metric, counts, ids=None, validate=True, pairwise_func=None,
                   id_pairs=None, **kwargs):
    """Compute distances between all pairs of samples

    Parameters
    ----------
    metric : str, callable
        The pairwise distance function to apply. See the scipy ``pdist`` docs
        and the scikit-bio functions linked under *See Also* for available
        metrics. Passing metrics as a strings is preferable as this often
        results in an optimized version of the metric being used.
    counts : 2D array_like of ints or floats
        Matrix containing count/abundance data where each row contains counts
        of OTUs in a given sample.
    ids : iterable of strs, optional
        Identifiers for each sample in ``counts``. By default, samples will be
        assigned integer identifiers in the order that they were provided
        (where the type of the identifiers will be ``str``).
    validate : bool, optional
        If `False`, validation of the input won't be performed. This step can
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
    id_pairs : iterable of tuple, optional
        An iterable of tuples of IDs to compare (e.g., ``[('a', 'b'), ('a',
        'c'), ...])``. If specified, the set of IDs described must be a subset
        of ``ids``. ``id_pairs`` is a mutually exclusive argument with
        ``pairwise_func``. ``metric`` must be resolvable by scikit-bio (e.g.,
        UniFrac methods), or must be callable.
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
    TypeError
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
    if validate:
        counts = _validate_counts_matrix(counts, ids=ids)

    if metric == 'unweighted_unifrac':
        otu_ids, tree, kwargs = _get_phylogenetic_kwargs(counts, **kwargs)
        metric, counts_by_node = _setup_multiple_unweighted_unifrac(
                counts, otu_ids=otu_ids, tree=tree, validate=validate)
        counts = counts_by_node
    elif metric == 'weighted_unifrac':
        # get the value for normalized. if it was not provided, it will fall
        # back to the default value inside of _weighted_unifrac_pdist_f
        normalized = kwargs.pop('normalized',
                                _normalize_weighted_unifrac_by_default)
        otu_ids, tree, kwargs = _get_phylogenetic_kwargs(counts, **kwargs)
        metric, counts_by_node = _setup_multiple_weighted_unifrac(
                counts, otu_ids=otu_ids, tree=tree, normalized=normalized,
                validate=validate)
        counts = counts_by_node
    elif callable(metric):
        metric = functools.partial(metric, **kwargs)
        # remove all values from kwargs, since they have already been provided
        # through the partial
        kwargs = {}
    else:
        # metric is a string that scikit-bio doesn't know about, for
        # example one of the SciPy metrics
        pass

    if id_pairs is not None:
        if pairwise_func is not None:
            raise ValueError("`pairwise_func` is not compatible with "
                             "`id_pairs`")

        pairwise_func = functools.partial(_partial_pw, ids, id_pairs)

    if pairwise_func is None:
        pairwise_func = scipy.spatial.distance.pdist

    distances = pairwise_func(counts, metric=metric, **kwargs)
    return DistanceMatrix(distances, ids)
