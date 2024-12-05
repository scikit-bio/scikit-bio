# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from itertools import combinations

import numpy as np

from skbio.stats.distance import DistanceMatrix


def rf_dists(trees, ids=None, shared_by_all=True, proportion=False, rooted=False):
    r"""Calculate pairwise Robinson-Foulds (RF) distances among trees.

    .. versionadded:: 0.6.3

    Parameters
    ----------
    trees : list of TreeNode
        Input trees.
    ids : list of str, optional
        Identifiers of input trees. Must be unique. If omitted, will use incremental
        integers "0", "1", "2",...
    shared_by_all : bool, optional
        Distance between any pair of trees is calculated based on taxa that are shared
        across all trees (True, default), or shared between the current pair of trees
        (False).
    proportion : bool, optional
        Whether to return the RF distance as count (False, default) or proportion
        (True).
    rooted : bool, optional
        Whether to consider the trees as unrooted (False, default) or rooted (True).

    Returns
    -------
    DistanceMatrix
        Robinson-Foulds distance matrix.

    Notes
    -----
    The Robinson-Foulds (RF) distance, a.k.a. symmetric difference, is a measure of
    topological dissimilarity between two trees. It was originally described in
    [1]_. It is calculated as the number of bipartitions that differ between two
    unrooted trees.

    .. math::

        \text{RF}(T_1, T_2) = |S_1 \triangle S_2| = |(S_1 \setminus S_2) \cup (S_2
        \setminus S_1)|

    where :math:`S_1` and :math:`S_2` are the sets of bipartitions of trees
    :math:`T_1` and :math:`T_2`, respectively.

    For rooted trees, the RF distance is calculated as the number of unshared
    clades (subsets of taxa) [2]_.

    By specifying ``proportion=True``, a unit distance will be returned, ranging
    from 0 (identical) to 1 (completely different).

    See Also
    --------
    TreeNode.compare_rfd

    References
    ----------
    .. [1] Robinson, D. F., & Foulds, L. R. (1981). Comparison of phylogenetic
        trees. Mathematical biosciences, 53(1-2), 131-147.

    .. [2] Bogdanowicz, D., & Giaro, K. (2013). On a matching distance between
        rooted phylogenetic trees. International Journal of Applied Mathematics
        and Computer Science, 23(3), 669-684.

    """
    n_trees = len(trees)
    method = "subsets" if rooted else "biparts"
    result = np.zeros((n_trees, n_trees))
    taxon_sets = [t.subset() for t in trees]
    if shared_by_all:
        within = frozenset.intersection(*taxon_sets)
    for i, j in combinations(range(n_trees), 2):
        if not shared_by_all:
            within = taxon_sets[i] & taxon_sets[j]
        result[i, j] = result[j, i] = _compare_topology(
            trees[i], trees[j], method, within=within, proportion=proportion
        )
    return DistanceMatrix(result, ids, validate=False)


def _compare_topology(
    tree1,
    tree2,
    method="subsets",
    within=None,
    shared_only=True,
    proportion=False,
    symmetric=True,
    include_single=False,
    weighted=False,
    metric="euclidean",
):
    r"""Calculate the topological difference between two trees.

    This function calculates the Robinson-Foulds (RF) distance or its derivates.

    Parameters
    ----------
    tree1 : TreeNode
        The first tree to compare with.
    tree2 : TreeNode
        The second tree to compare with.
    method : str, optional
        Subsets or bipartitions.
    within : set of str, optional
        Refine to given taxa.
    shared_only : bool, optional
        Refine to shared taxa.
    proportion : bool, optional
        Normalize to fraction.
    symmetric : bool, optional
        Symmetric difference.
    include_single : bool, optional
        Include singletons.
    weighted : bool, optional
        Weight by branch length.
    metric : str or callable, optional
        Pairwise distance metric.

    Returns
    -------
    float
        Difference between self and other.

    See Also
    --------
    compare_rfd
    compare_wrfd
    compare_subsets
    compare_biparts

    """
    topo1, topo2 = getattr(tree1, method), getattr(tree2, method)
    kwargs = dict(include_single=include_single, map_to_length=weighted)
    if within:
        sets1, sets2 = topo1(within=within, **kwargs), topo2(within=within, **kwargs)
    elif shared_only:
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
            result = sets1.symmetric_difference(sets2)
        else:
            result = sets1.difference(sets2)
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
        if isinstance(metric, str):
            result = getattr(spdist, metric)(L1, L2)
        else:
            result = metric(L1, L2)

    return result
