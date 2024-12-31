# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from itertools import combinations

import numpy as np
import scipy.spatial.distance as spdist
from scipy.spatial.distance import squareform

from skbio.stats.distance import DistanceMatrix
from skbio.util import get_rng


def unitcorr(a, b):
    """Calculate unit correlation distance.

    d = (1 - r) / 2, in which r is the Pearson's correlation coefficient.

    """
    return spdist.correlation(a, b) * 0.5


def _check_input(trees, ids):
    """Validate input trees and IDs."""
    if not (n := len(trees)):
        raise ValueError(f"No tree is provided.")
    if ids is None:
        return
    if (m := len(ids)) != n:
        raise ValueError(f"Numbers of IDs and trees don't match.")
    if m > len(set(ids)):
        raise ValueError(f"IDs contain duplicates.")


def _check_dist_metric(metric):
    """Validate distance metric."""
    if isinstance(metric, str):
        if metric == "unitcorr":
            metric = unitcorr
        else:
            metric = getattr(spdist, metric)
    elif not callable(metric):
        raise ValueError("`metric` must be a string or callable.")
    return metric


def _check_shuffler(shuffler):
    """Validate sample shuffler."""
    if not callable(shuffler):
        shuffler = get_rng(shuffler).shuffle
    return shuffler


def _calc_dists(trees, ids, func, kwargs, shared_by_all):
    """Calculate paiwise distances among trees."""
    if shared_by_all:
        result = func(trees, **kwargs)
        result = squareform(np.array(result, dtype=float), checks=False)
        return DistanceMatrix(result, ids, validate=False)
    else:
        return DistanceMatrix.from_iterable(
            trees,
            metric=lambda x, y: func((x, y), **kwargs)[0],
            keys=ids,
            validate=False,
        )


def _setdiff(a, b, proportion):
    r"""Quantify the difference between two sets.

    For Robinson-Foulds (RF) distance calculation.

    """
    # calculate symmetric difference
    result = len(a ^ b)

    # normalize result to unit range [0, 1]
    # if total is 0, return 1 (dist = 1 means saturation)
    if proportion:
        total = len(a) + len(b)
        result = result / total if total else 1.0

    # cast result to float
    else:
        result = float(result)

    return result


def _mapdiff(a, b, metric):
    r"""Quantify the difference between two mappings.

    For weighted Robinson-Foulds (wRF) distance calculation.

    """
    # unique keys are assigned value 0
    union = frozenset(a).union(b)
    La = [a.get(x, 0.0) for x in union]
    Lb = [b.get(x, 0.0) for x in union]
    return metric(La, Lb)


def _topo_dists(
    trees,
    rooted=False,
    shared_only=True,
    proportion=False,
    include_tips=False,
    weighted=False,
    metric=None,
):
    r"""Calculate the topological distances among multiple trees.

    This function calculates the Robinson-Foulds (RF) distance or its variants. These
    metrics rely on the compatibility of branches between two trees.

    See :meth:`TreeNode.compare_rfd` and :meth:`~TreeNode.compare_wrfd` for details of
    the metrics and parameters.

    Parameters
    ----------
    trees : list of TreeNode
        The trees to compare.
    rooted : bool, optional
        Rooted or unrooted.
    shared_only : bool, optional
        Refine to shared taxa.
    proportion : bool, optional
        Normalize to proportion.
    include_tips : bool, optional
        Include singletons.
    weighted : bool, optional
        Weight by branch length.
    metric : callable, optional
        Distance metric (must provide if weighted).

    Returns
    -------
    list of float
        Pairwise path-length distances or variants in condensed form.

    See Also
    --------
    rf_dists
    wrfd_dists
    TreeNode.compare_rfd
    TreeNode.compare_wrfd
    TreeNode.compare_subsets
    TreeNode.compare_biparts

    """
    kwargs = dict(include_tips=include_tips, map_to_length=weighted)

    # get topological units reflecting shared taxa
    if shared_only:
        taxon_sets = [t.subset() for t in trees]
        shared = frozenset.intersection(*taxon_sets)
        n_shared = len(shared)

        # Decide the "within" parameter (None or shared taxa) for each method call,
        # which is more efficient when within=None.
        withins = [shared if len(s) > n_shared else None for s in taxon_sets]

        # get topological units from all trees
        if rooted:
            sets = [t.subsets(within=w, **kwargs) for t, w in zip(trees, withins)]
        else:
            sets = [
                t.biparts(within=w, full=s, **kwargs)
                for t, w, s in zip(trees, withins, taxon_sets)
            ]

    # get all topological units
    else:
        sets = [(t.subsets if rooted else t.biparts)(**kwargs) for t in trees]

    # unweighted (symmetric difference between sets)
    if not weighted:
        result = [_setdiff(a, b, proportion) for a, b in combinations(sets, 2)]

    # branch length weighted (vector distance)
    else:
        result = [_mapdiff(a, b, metric) for a, b in combinations(sets, 2)]

    return result


def rf_dists(trees, ids=None, shared_by_all=True, proportion=False, rooted=False):
    r"""Calculate Robinson-Foulds (RF) distances among trees.

    .. versionadded:: 0.6.3

    Parameters
    ----------
    trees : list of TreeNode
        Input trees.
    ids : list of str, optional
        Unique identifiers of input trees. If omitted, will use incremental integers
        "0", "1", "2",...
    shared_by_all : bool, optional
        Calculate the distance between each pair of trees based on taxa shared across
        all trees (True, default), or shared between the current pair of trees (False).
    proportion : bool, optional
        Whether to return the RF distance as count (False, default) or proportion
        (True).
    rooted : bool, optional
        Whether to consider the trees as unrooted (False, default) or rooted (True).

    Returns
    -------
    DistanceMatrix
        Matrix of the Robinson-Foulds distances.

    See Also
    --------
    TreeNode.compare_rfd

    Notes
    -----
    The Robinson-Foulds (RF) distance [1]_, a.k.a. symmetric difference, is the number
    of bipartitions differing between two trees.

    This function is equivalent to :meth:`TreeNode.compare_rfd` for two trees. Refer
    to the latter for details about the metric and its parameters. However, the current
    function extends the operation to an arbitrary number of trees and returns a
    distance matrix for them.

    This function is optimized for calculation based on taxa shared across all trees.
    One can instead set ``shared_by_all`` to False to calculate based on taxa shared
    between each pair of trees, which is however less efficient since bipartitions need
    to be re-inferred during each comparison.

    References
    ----------
    .. [1] Robinson, D. F., & Foulds, L. R. (1981). Comparison of phylogenetic trees.
       Mathematical biosciences, 53(1-2), 131-147.

    Examples
    --------
    >>> from skbio import TreeNode
    >>> trees = [TreeNode.read([x]) for x in (
    ...     "(((a,b),c),d,e);",
    ...     "((a,(b,c)),d,e);",
    ...     "((a,b),(c,d),e);",
    ...     "(a,b,(c,(d,e)));",
    ... )]
    >>> dm = rf_dists(trees, ids=list("ABCD"))
    >>> print(dm)
    4x4 distance matrix
    IDs:
    'A', 'B', 'C', 'D'
    Data:
    [[ 0.  2.  2.  0.]
     [ 2.  0.  4.  2.]
     [ 2.  4.  0.  2.]
     [ 0.  2.  2.  0.]]

    """
    _check_input(trees, ids)
    kwargs = dict(
        rooted=rooted,
        shared_only=True,
        proportion=proportion,
        include_tips=False,
        weighted=False,
        metric=None,
    )
    return _calc_dists(trees, ids, _topo_dists, kwargs, shared_by_all)


def wrf_dists(
    trees,
    ids=None,
    shared_by_all=True,
    metric="cityblock",
    rooted=False,
    include_tips=True,
):
    r"""Calculate weighted Robinson-Foulds (wRF) distances or variants among trees.

    .. versionadded:: 0.6.3

    Parameters
    ----------
    trees : list of TreeNode
        Input trees.
    ids : list of str, optional
        Unique identifiers of input trees. If omitted, will use incremental integers
        "0", "1", "2",...
    shared_by_all : bool, optional
        Calculate the distance between each pair of trees based on taxa shared across
        all trees (True, default), or shared between the current pair of trees (False).
    metric : str or callable, optional
        The distance metric to use. Can be a preset, a distance function name under
        :mod:`scipy.spatial.distance`, or a custom function that takes two vectors and
        returns a number. See :meth:`~TreeNode.compare_wrfd` for details.
    rooted : bool, optional
        Whether to consider the trees as unrooted (False, default) or rooted (True).
    include_tips : bool, optional
        Whether to include single-taxon biparitions (terminal branches) in the
        calculation. Default is True.

    Returns
    -------
    DistanceMatrix
        Matrix of weighted Robinson-Foulds distances or variants.

    See Also
    --------
    TreeNode.compare_wrfd

    Notes
    -----
    The weighted Robinson-Foulds (wRF) distance [1]_ is the sum of differences of
    branch lengths of matching bipartitions between a pair of trees.

    This function is equivalent to :meth:`TreeNode.compare_wrfd` for two trees. Refer
    to the latter for details about the metric and its variants, and the parameter
    settings for calculating them. However, the current function extends the operation
    to an arbitrary number of trees and returns a distance matrix for them.

    A restriction of the current function compared to ``compare_wrfd`` is that
    ``metric`` must be symmetric (i.e., :math:`d(x, y) = d(y, x)`), and equals zero
    from a vector to itself (i.e., :math:`d(x, x) = 0`). It does not have to suffice
    non-negativity or triangle inequality though.

    This function is optimized for calculation based on taxa shared across all trees.
    One can instead set ``shared_by_all`` to False to calculate based on taxa shared
    between each pair of trees, which is however less efficient as bipartitions need to
    be re-inferred during each comparison.

    References
    ----------
    .. [1] Robinson, D. F., & Foulds, L. R. (1979) Comparison of weighted labelled
       trees. In Combinatorial Mathematics VI: Proceedings of the Sixth Australian
       Conference on Combinatorial Mathematics, Armidale, Australia (pp. 119-126).

    Examples
    --------
    >>> from skbio import TreeNode
    >>> trees = [TreeNode.read([x]) for x in (
    ...     "((a:1,b:2):1,c:4,((d:4,e:5):2,f:6):1);",
    ...     "((a:3,(b:2,c:2):1):3,d:8,(e:5,f:6):2);",
    ...     "((a:1,c:6):2,(b:3,(d:2,e:3):1):2,f:7);",
    ... )]
    >>> dm = wrf_dists(trees, ids=list("ABC"))
    >>> print(dm)
    3x3 distance matrix
    IDs:
    'A', 'B', 'C'
    Data:
    [[  0.  16.  15.]
     [ 16.   0.  27.]
     [ 15.  27.   0.]]

    """
    _check_input(trees, ids)
    metric = _check_dist_metric(metric)
    kwargs = dict(
        rooted=rooted,
        shared_only=True,
        include_tips=include_tips,
        weighted=True,
        metric=metric,
    )
    return _calc_dists(trees, ids, _topo_dists, kwargs, shared_by_all)


def _path_dists(trees, sample, metric, shuffler, use_length, ignore_self):
    r"""Calculate path-length distances or variants among multiple trees.

    This function calculates several tree distance metrics based on the lengths of
    paths connecting pairs of taxa, a.k.a., cophenetic distances.

    See :meth:`TreeNode.compare_cophenet` for details of the metrics and parameters.

    Parameters
    ----------
    trees : list of TreeNode
        The trees to compare.
    sample : int, optional
        Number of taxa to sample.
    metric : callable, optional
        Distance function.
    shuffler : callable, optional
        Shuffling function.
    use_length : bool, optional
        Use branch lengths.
    ignore_self : bool, optional
        Ignore taxon to itself.

    Returns
    -------
    list of float
        Pairwise path-length distances or variants in condensed form.

    See Also
    --------
    path_dists
    TreeNode.compare_cophenet

    """
    # Get shared taxa among trees.
    shared = frozenset.intersection(*(t.subset() for t in trees))

    # Sort shared taxa. This is for the stability of sampling and distance calculation
    # as Python set is unordered.
    shared = sorted(shared)

    # Randomly sample a subset of taxa.
    if sample is not None:
        if (n_shared := len(shared)) < sample:
            raise ValueError(
                f"Sample size ({sample}) is greater than the number of shared taxa "
                f"({n_shared})."
            )
        shuffler(shared)
        shared = shared[:sample]

    # Generate cophenetic (tip-to-tip) distance matrices.
    paths = [x.cophenet(endpoints=shared, use_length=use_length) for x in trees]

    # Convert square matrices into condensed 1-D arrays, such that diagonal (d(x, x) =
    # 0) and symmetric (d(x, y) = d(y, x)) comparisons are omitted.
    if ignore_self:
        paths = [x.condensed_form() for x in paths]

    # Otherwise, create 1-D views of the entire square matrics. `ravel` copies the data
    # only when necessary, and in this case doesn't.
    else:
        paths = [x.data.ravel() for x in paths]

    # Calculate pairwise distances among matrices.
    return [metric(a, b) for a, b in combinations(paths, 2)]


def path_dists(
    trees,
    ids=None,
    shared_by_all=True,
    metric="euclidean",
    use_length=True,
    sample=None,
    shuffler=None,
):
    r"""Calculate path-length distances or variants among trees.

    .. versionadded:: 0.6.3

    Parameters
    ----------
    trees : list of TreeNode
        Input trees.
    ids : list of str, optional
        Unique identifiers of input trees. If omitted, will use incremental integers
        "0", "1", "2",...
    shared_by_all : bool, optional
        Calculate the distance between each pair of trees based on taxa shared across
        all trees (True, default), or shared between the current pair of trees (False).
    metric : str or callable, optional
        The distance metric to use. Can be a preset, a distance function name under
        :mod:`scipy.spatial.distance`, or a custom function that takes two vectors and
        returns a number. See :meth:`~TreeNode.compare_cophenet` for details.
    use_length : bool, optional
        Whether to calculate the sum of branch lengths (True, default) or the
        number of branches (False) connecting each pair of tips.
    sample : int, optional
        Randomly subsample this number of shared taxa for calculation. This is useful
        when comparing very large trees.
    shuffler : int, np.random.Generator or callable, optional
        The shuffling function to use if ``sample`` is specified. Default is
        :meth:`~numpy.random.Generator.shuffle`. If an integer is provided, a random
        generator will be constructed using this number as the seed.

    Returns
    -------
    DistanceMatrix
        Matrix of the path-length distances or variants.

    See Also
    --------
    TreeNode.compare_cophenet

    Notes
    -----
    The path-length distance [1]_ is the square root of the sum of squared differences
    of path lengths among all pairs of taxa between two trees.

    This function is equivalent to :meth:`TreeNode.compare_cophenet` for two trees.
    Refer to the latter for details about the metric and its variants, and parameter
    settings for calculating them. However, the current function extends the operation
    to an arbitrary number of trees and returns a distance matrix for them. It is named
    so because the term "cophenetic distance" refers to the distance between two taxa
    in a tree instead.

    A restriction of the current function compared to ``compare_cophenet`` is that
    ``metric`` must be symmetric (i.e., :math:`d(x, y) = d(y, x)`), and equals zero
    from a vector to itself (i.e., :math:`d(x, x) = 0`). It does not have to suffice
    non-negativity or triangle inequality though.

    This function is optimized for calculation based on taxa shared across all trees.
    One can instead set ``shared_by_all`` to False to calculate based on taxa shared
    between each pair of trees, which is however less efficient as the path lengths
    need to be re-calculated during each comparison.

    References
    ----------
    .. [1] Lapointe, F. J., & Cucumel, G. (1997). The average consensus procedure:
       combination of weighted trees containing identical or overlapping sets of
       taxa. Systematic Biology, 46(2), 306-312.

    Examples
    --------
    >>> from skbio import TreeNode
    >>> trees = [TreeNode.read([x]) for x in (
    ...     "((a:1,b:2):1,c:4,((d:4,e:5):2,f:6):1);",
    ...     "((a:3,(b:2,c:2):1):3,d:8,(e:5,f:6):2);",
    ...     "((a:1,c:6):2,(b:3,(d:2,e:3):1):2,f:7);",
    ... )]
    >>> dm = path_dists(trees, ids=list("ABC"))
    >>> print(dm)
    3x3 distance matrix
    IDs:
    'A', 'B', 'C'
    Data:
    [[  0.          13.7113092   11.87434209]
     [ 13.7113092    0.          19.5192213 ]
     [ 11.87434209  19.5192213    0.        ]]

    """
    _check_input(trees, ids)
    metric = _check_dist_metric(metric)
    if sample is not None:
        shuffler = _check_shuffler(shuffler)
    kwargs = dict(
        sample=sample,
        metric=metric,
        shuffler=shuffler,
        use_length=use_length,
        ignore_self=True,
    )
    return _calc_dists(trees, ids, _path_dists, kwargs, shared_by_all)
