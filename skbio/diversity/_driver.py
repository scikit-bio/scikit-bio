# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from functools import partial

from scipy.spatial.distance import pdist
import pandas as pd
import numpy as np

import skbio
from skbio.util._decorator import experimental
from skbio.stats.distance import DistanceMatrix
from skbio.tree import DuplicateNodeError, MissingNodeError
from skbio.diversity._phylogenetic import _nodes_by_counts


def _get_alpha_diversity_metrics():
    return {'ace': skbio.diversity.alpha.ace,
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
            'lladser_ci': skbio.diversity.alpha.lladser_ci
            }


@experimental(as_of="0.4.0-dev")
def alpha_diversity(metric, counts, ids=None, **kwargs):
    """ Compute alpha diversity for one or more count vectors

    Parameters
    ----------
    metric: str, callable
        The alpha diversity metric to apply to the count vector(s).
        Passing metric as a string is preferable as this often results in an
        optimized version of the metric being used.
    counts : 1D or 2D array_like of ints or floats
        Vector or matrix containing count/abundance data. If a matrix, each row
        should contain counts of observations in a given sample.
    ids : iterable of strs, optional
        Identifiers for each sample in ``counts``.

    Returns
    -------
    pd.Series
        Values of metric for all vectors provided in ``counts``. The index
        will be ``ids``, if provided.

    Raises
    ------
    ValueError
        If ``len(ids) != len(counts)``, or if ``otu_ids`` and ``tree`` are not
        provided when ``metric=faith_pd``.

    See Also
    --------
    skbio.diversity.beta.beta_diversity

    Notes
    -----
    The value that you provide for for ``metric`` can be either a string (e.g.,
    ``"faith_pd"``) or a function
    (e.g., ``skbio.diversity.alpha.faith_pd``). The metric should
    generally be passed as a string, as this often uses an optimized version
    of the metric. For example, passing  ``"faith_pd"`` (a string) will be
    tens of times faster than passing the function
    ``skbio.diversity.alpha.faith_pd``. The latter may be faster if computing
    alpha diversity for only one or a few distances, but in these cases the
    difference in runtime is negligible, so it's safer to just err on the side
    of passing ``metric`` as a string.

    """
    metrics = _get_alpha_diversity_metrics()

    counts = np.atleast_2d(counts)
    num_samples = len(counts)
    if ids is not None and num_samples != len(ids):
        raise ValueError(
            "Number of rows in counts must be equal to number of provided "
            "ids.")

    if metric == 'faith_pd':
        if "otu_ids" not in kwargs:
            raise ValueError("otu_ids is required when metric is faith_pd.")
        if "tree" not in kwargs:
            raise ValueError("tree is required when metric=faith_pd")
        counts_by_node, tree_index, branch_lengths = \
            _vectorize_counts_and_tree(counts, kwargs['otu_ids'],
                                       kwargs['tree'])
        counts = counts_by_node
        metric = partial(skbio.diversity.alpha._faith_pd._faith_pd,
                         branch_lengths=branch_lengths)
    elif callable(metric):
        metric = partial(metric, **kwargs)
    elif metric in metrics:
        metric = metrics[metric]
    else:
        raise ValueError('Unknown metric provided: %r.' % metric)

    results = map(metric, counts)
    return pd.Series(results, index=ids)


@experimental(as_of="0.4.0")
def beta_diversity(metric, counts, ids=None, **kwargs):
    """Compute distances between all pairs of columns in a counts matrix

    Parameters
    ----------
    metric : str, callable
        The pairwise distance function to apply. See the scipy ``pdist`` docs
        and the scikit-bio functions linked under *See Also* for available
        metrics. Passing metrics as a string is preferable as this often
        results in an optimized version of the metric being used.
    counts : 2D array_like of ints or floats
        Matrix containing count/abundance data where each row contains counts
        of observations in a given sample.
    ids : iterable of strs, optional
        Identifiers for each sample in ``counts``.

    Returns
    -------
    skbio.DistanceMatrix
        Distances between all pairs of samples (i.e., rows). The number of
        row and columns will be equal to the number of rows in ``counts``.

    Raises
    ------
    ValueError
        If ``len(ids) != len(counts)``.

    See Also
    --------
    unweighted_unifrac
    weighted_unifrac
    scipy.spatial.distance.pdist
    skbio.diversity.alpha.alpha_diversity

    Notes
    -----
    The value that you provide for for ``metric`` can be either a string (e.g.,
    ``"unweighted_unifrac"``) or a function
    (e.g., ``skbio.diversity.beta.unweighted_unifrac``). The metric should
    generally be passed as a string, as this often uses an optimized version
    of the metric. For example, passing  ``"unweighted_unifrac"`` (a string)
    will be hundreds of times faster than passing the function
    ``skbio.diversity.beta.unweighted_unifrac``. The latter is faster if
    computing only one or a few distances, but in these cases the difference in
    runtime is negligible, so it's safer to just err on the side of passing
    ``metric`` as a string.

    """
    num_samples = len(counts)
    if ids is not None and num_samples != len(ids):
        raise ValueError(
            "Number of rows in counts must be equal to number of provided "
            "ids.")

    if metric == 'unweighted_unifrac':
        f = skbio.diversity.beta._unifrac._unweighted_unifrac_pdist_f
        metric, counts, _ = f(counts, otu_ids=kwargs['otu_ids'],
                              tree=kwargs['tree'])
    elif metric == 'weighted_unifrac':
        normalized = kwargs.get('normalized', False)
        f = skbio.diversity.beta._unifrac._weighted_unifrac_pdist_f
        metric, counts, _ = f(counts, otu_ids=kwargs['otu_ids'],
                              tree=kwargs['tree'], normalized=normalized)
    elif callable(metric):
        metric = partial(metric, **kwargs)
    else:
        pass

    distances = pdist(counts, metric)
    return DistanceMatrix(distances, ids)


def _validate_counts_vector(counts, suppress_cast=False):
    """Validate and convert input to an acceptable counts vector type.

    Note: may not always return a copy of `counts`!

    """
    counts = np.asarray(counts)

    if not suppress_cast:
        counts = counts.astype(int, casting='safe', copy=False)

    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")
    elif (counts < 0).any():
        raise ValueError("Counts vector cannot contain negative values.")

    return counts


def _validate_counts_vectors(*args, **kwargs):
    results = []
    lens = []
    # py2-compatible mechanism for specifying a keyword argument when also
    # passing *args derived from SO answer:
    # http://stackoverflow.com/a/15302038/3424666
    suppress_cast = kwargs.pop('suppress_cast', False)
    for counts in args:
        results.append(_validate_counts_vector(counts, suppress_cast))
        lens.append(len(counts))
    if len(set(lens)) > 1:
        raise ValueError("Input vectors u_counts and v_counts must be of "
                         "equal length.")

    return results


def _validate_otu_ids_and_tree(counts, otu_ids, tree):
    # all otu_ids are unique
    # len(otu_ids) == len(counts)
    len_otu_ids = len(otu_ids)
    set_otu_ids = set(otu_ids)
    if len_otu_ids != len(set_otu_ids):
        raise ValueError("OTU IDs vector cannot contain duplicated ids.")
    if len(counts) != len_otu_ids:
        raise ValueError("OTU IDs vector must be the same length as counts "
                         "vector(s).")

    # the tree is rooted
    if len(tree.root().children) > 2:
        # this is an imperfect check for whether the tree is rooted or not.
        # can this be improved?
        raise ValueError("Tree must be rooted.")

    # all nodes (except the root node) have corresponding branch lengths
    # all tip names in tree are unique
    # all otu_ids correspond to tip names in tree
    branch_lengths = []
    tip_names = []
    for e in tree.traverse():
        if not e.is_root():
            branch_lengths.append(e.length)
        if e.is_tip():
            tip_names.append(e.name)
    set_tip_names = set(tip_names)
    if len(tip_names) != len(set_tip_names):
        raise DuplicateNodeError("All tip names must be unique.")
    if np.array([l is None for l in branch_lengths]).any():
        raise ValueError("All non-root nodes in tree must have a branch "
                         "length.")
    missing_tip_names = set_otu_ids - set_tip_names
    if missing_tip_names != set():
        raise MissingNodeError("All otu_ids must be present as tip names in "
                               "tree. Tree is missing tips with names: %s"
                               % " ".join(missing_tip_names))


def _vectorize_counts_and_tree(counts, otu_ids, tree, tree_index=None):
    """ Index tree and convert counts to np.array in corresponding order
    """
    if tree_index is None:
        tree_index = tree.to_array(nan_length_value=0.0)
    otu_ids = np.asarray(otu_ids)
    counts = np.atleast_2d(counts)
    counts_by_node = _nodes_by_counts(counts, otu_ids, tree_index)
    branch_lengths = tree_index['length']

    # branch_lengths is just a reference to the array inside of tree_index,
    # but it's used so much that it's convenient to just pull it out here.
    return counts_by_node.T, tree_index, branch_lengths
