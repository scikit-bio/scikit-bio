# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from functools import partial
from itertools import combinations

import numpy as np
import scipy.spatial.distance as spdist

from skbio.stats.distance import DistanceMatrix
from skbio.util import get_rng


def unitcorr(a, b):
    """Calculate unit correlation distance."""
    return spdist.correlation(a, b) * 0.5


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


def _check_ids(trees, ids):
    """Check tree IDs."""
    if ids is None:
        return
    if (n := len(ids)) != len(trees):
        raise ValueError(f"ID number does not match tree number.")
    if n > len(set(ids)):
        raise ValueError(f"IDs contain duplicates.")


def _shared(trees):
    """Get shared taxa among multiple trees."""
    taxon_sets = [t.subset() for t in trees]
    shared = frozenset.intersection(*taxon_sets)
    if not shared:
        raise ValueError("No taxon is shared across all trees.")
    return shared, taxon_sets


def _withins(trees):
    """Set the "within" parameter based on shared taxa among multiple trees."""
    shared, taxon_sets = _shared(trees)
    n_shared = len(shared)
    return [shared if len(s) > n_shared else None for s in taxon_sets]


def _sample_taxa(sample, taxa, shuffler):
    """Sample a given number of taxa using a given shuffler."""
    if (n_taxa := len(taxa)) < sample:
        raise ValueError(
            f"{sample} taxa are to be sampled whereas only {n_taxa} taxa are shared "
            "between the trees."
        )
    shuffler(taxa)
    return taxa[:sample]


def _topo_dist(
    tree1,
    tree2,
    method="biparts",
    shared_only=True,
    proportion=False,
    symmetric=True,
    include_single=False,
    weighted=False,
    metric=None,
):
    r"""Calculate the topological distance between two trees.

    This function calculates the Robinson-Foulds (RF) distance or its derivates.

    Parameters
    ----------
    tree1 : TreeNode
        The first tree to compare.
    tree2 : TreeNode
        The second tree to compare.
    method : str, optional
        "subsets" or "biparts"
    shared_only : bool, optional
        Refine to shared taxa.
    proportion : bool, optional
        Normalize to proportion.
    symmetric : bool, optional
        Symmetric difference.
    include_single : bool, optional
        Include singletons.
    weighted : bool, optional
        Weight by branch length.
    metric : callable, optional
        Distance metric (must provide if weighted).

    Returns
    -------
    float
        Difference between the two trees.

    See Also
    --------
    rf_dists
    wrfd_dists
    TreeNode.compare_rfd
    TreeNode.compare_wrfd
    TreeNode.compare_subsets
    TreeNode.compare_biparts

    """
    topo1, topo2 = getattr(tree1, method), getattr(tree2, method)
    kwargs = dict(include_single=include_single, map_to_length=weighted)
    if shared_only:
        set1, set2 = tree1.subset(), tree2.subset()
        n_shared = len(shared := set1 & set2)
        if len(set1) > n_shared:
            sets1 = topo1(within=shared, **kwargs)
        else:
            sets1 = topo1(**kwargs)
        if len(set2) > n_shared:
            sets2 = topo2(within=shared, **kwargs)
        else:
            sets2 = topo2(**kwargs)
    else:
        sets1, sets2 = topo1(**kwargs), topo2(**kwargs)

    # unweighted (set difference)
    if not weighted:
        if symmetric:
            result = sets1 ^ sets2  # symmetric difference
        else:
            result = sets1 - sets2  # difference
        result = len(result)

        # normalize result to unit range [0, 1]
        # if total is 0, return 1 (dist = 1 means saturation)
        if proportion:
            total = len(sets1) + (symmetric and len(sets2))
            result = result / total if total else 1.0

        # cast result to float
        else:
            result = float(result)

    # branch length weighted (vector distance)
    else:
        union = frozenset(sets1).union(sets2)
        L1 = [sets1.get(x, 0.0) for x in union]
        L2 = [sets2.get(x, 0.0) for x in union]
        result = metric(L1, L2)

    return result


def rf_dists(trees, ids=None, pairwise=False, proportion=False, rooted=False):
    r"""Calculate Robinson-Foulds (RF) distances among trees.

    .. versionadded:: 0.6.3

    Parameters
    ----------
    trees : list of TreeNode
        Input trees.
    ids : list of str, optional
        Unique identifiers of input trees. If omitted, will use incremental integers
        "0", "1", "2",...
    pairwise : bool, optional
        Calculate the distance between each pair of trees based on taxa shared across
        all trees (False, default), or shared between the current pair of trees (True).
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
    One can instead set ``pairwise`` to True to calculate based on taxa shared between
    each pair of trees, which is however less efficient since bipartitions need to be
    re-inferred during each comparison.

    References
    ----------
    .. [1] Robinson, D. F., & Foulds, L. R. (1981). Comparison of phylogenetic
        trees. Mathematical biosciences, 53(1-2), 131-147.

    Examples
    --------
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
    _check_ids(trees, ids)
    method = "subsets" if rooted else "biparts"

    if pairwise:
        metric = partial(_topo_dist, method=method, proportion=proportion)
        return DistanceMatrix.from_iterable(
            trees, metric=metric, keys=ids, validate=False
        )

    withins = _withins(trees)
    sets = [getattr(x, method)(within=y) for x, y in zip(trees, withins)]
    lens = [len(s) for s in sets]

    result = np.zeros((n_trees := len(trees), n_trees))
    for i, j in combinations(range(n_trees), 2):
        rf = len(sets[i] ^ sets[j])
        if proportion:
            total = lens[i] + lens[j]
            rf = rf / total if total else 1.0
        result[i, j] = result[j, i] = rf

    return DistanceMatrix(result, ids, validate=False)


def wrf_dists(
    trees,
    ids=None,
    pairwise=False,
    metric="cityblock",
    rooted=False,
    include_single=True,
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
    pairwise : bool, optional
        Calculate the distance between each pair of trees based on taxa shared across
        all trees (False, default), or shared between the current pair of trees (True).
    metric : str or callable, optional
        The distance metric to use. Can be a preset, a distance function name under
        :mod:`scipy.spatial.distance`, or a custom function that takes two vectors and
        returns a number. See :meth:`~TreeNode.compare_wrfd` for details.
    rooted : bool, optional
        Whether to consider the trees as unrooted (False, default) or rooted (True).
    include_single : bool, optional
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
    One can instead set ``pairwise`` to True to calculate based on taxa shared between
    each pair of trees, which is however less efficient as bipartitions need to be
    re-inferred during each comparison.

    References
    ----------
    .. [1] Robinson, D. F., & Foulds, L. R. (1979) Comparison of weighted labelled
        trees. In Combinatorial Mathematics VI: Proceedings of the Sixth Australian
        Conference on Combinatorial Mathematics, Armidale, Australia (pp. 119-126).

    Examples
    --------
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
    _check_ids(trees, ids)
    metric = _check_dist_metric(metric)
    method = "subsets" if rooted else "biparts"

    if pairwise:
        metric = partial(
            _topo_dist,
            method=method,
            weighted=True,
            metric=metric,
            include_single=include_single,
        )
        return DistanceMatrix.from_iterable(
            trees, metric=metric, keys=ids, validate=False
        )

    withins = _withins(trees)
    maps = [
        getattr(x, method)(within=y, include_single=include_single, map_to_length=True)
        for x, y in zip(trees, withins)
    ]
    gets = [x.get for x in maps]
    sets = [frozenset(x) for x in maps]

    result = np.zeros((n_trees := len(trees), n_trees))
    for i, j in combinations(range(n_trees), 2):
        union = sets[i] | sets[j]
        Li = [gets[i](x, 0.0) for x in union]
        Lj = [gets[j](x, 0.0) for x in union]
        result[i, j] = result[j, i] = metric(Li, Lj)

    return DistanceMatrix(result, ids, validate=False)


def _path_dist(
    tree1,
    tree2,
    sample=None,
    metric="unitcorr",
    shuffler=None,
    use_length=True,
    ignore_self=False,
):
    tipmap1 = {n.name: n for n in tree1.tips()}
    tipmap2 = {n.name: n for n in tree2.tips()}
    shared = [x for x in tipmap1 if x in tipmap2]
    if not shared:
        raise ValueError("No tips are in common between the two trees.")

    if sample is not None:
        shuffler = _check_shuffler(shuffler)
        shared = _sample_taxa(sample, shared, shuffler)

    tips1 = [tipmap1[x] for x in shared]
    tips2 = [tipmap2[x] for x in shared]

    dm1 = tree1.cophenet(endpoints=tips1, use_length=use_length)
    dm2 = tree2.cophenet(endpoints=tips2, use_length=use_length)

    if ignore_self:
        dm1 = dm1.condensed_form()
        dm2 = dm2.condensed_form()
    else:
        dm1 = dm1.data.flat
        dm2 = dm2.data.flat

    metric = _check_dist_metric(metric)
    return metric(dm1, dm2)


def path_dists(
    trees,
    ids=None,
    pairwise=False,
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
    pairwise : bool, optional
        Calculate the distance between each pair of trees based on taxa shared across
        all trees (False, default), or shared between the current pair of trees (True).
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
    One can instead set ``pairwise`` to True to calculate based on taxa shared between
    each pair of trees, which is however less efficient as the path lengths need to be
    re-calculated during each comparison.

    References
    ----------
    .. [1] Lapointe, F. J., & Cucumel, G. (1997). The average consensus procedure:
        combination of weighted trees containing identical or overlapping sets of
        taxa. Systematic Biology, 46(2), 306-312.

    Examples
    --------
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
    _check_ids(trees, ids)

    if sample is not None:
        shuffler = _check_shuffler(shuffler)

    if pairwise:
        metric = partial(
            _path_dist,
            sample=sample,
            metric=metric,
            shuffler=shuffler,
            use_length=use_length,
            ignore_self=True,
        )
        return DistanceMatrix.from_iterable(
            trees, metric=metric, keys=ids, validate=False
        )

    metric = _check_dist_metric(metric)

    shared = sorted(_shared(trees)[0])

    if sample is not None:
        shared = _sample_taxa(sample, shared, shuffler)

    paths = [
        x.cophenet(endpoints=shared, use_length=use_length).condensed_form()
        for x in trees
    ]

    result = np.zeros((n_trees := len(trees), n_trees))
    for i, j in combinations(range(n_trees), 2):
        result[i, j] = result[j, i] = metric(paths[i], paths[j])

    return DistanceMatrix(result, ids, validate=False)
