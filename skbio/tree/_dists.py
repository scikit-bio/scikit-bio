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

from skbio.tree import TreeNode
from skbio.stats.distance import DistanceMatrix


def _withins(trees):
    """Return shared taxa among multiple trees and whether each tree has unique taxa."""
    taxon_sets = [t.subset() for t in trees]
    shared = frozenset.intersection(*taxon_sets)
    if not shared:
        raise ValueError("No taxon is shared across all trees.")
    n_shared = len(shared)
    return [shared if len(s) > n_shared else None for s in taxon_sets]


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

    This function is optimized for calculation based on taxa shared across all trees.
    One can instead set ``shared_by_all`` to False to calculate based on taxa shared
    between each pair of trees, which is however less efficient as bipartitions need
    to be re-inferred during each comparison.

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
    if not shared_by_all:
        metric = partial(TreeNode.compare_rfd, proportion=proportion, rooted=rooted)
        return DistanceMatrix.from_iterable(trees, metric=metric, keys=ids)

    n_trees = len(trees)
    method = "subsets" if rooted else "biparts"
    withins = _withins(trees)
    sets = [getattr(x, method)(within=y) for x, y in zip(trees, withins)]
    lens = [len(s) for s in sets]
    result = np.zeros((n_trees, n_trees))
    for i, j in combinations(range(n_trees), 2):
        rf = len(sets[i] ^ sets[j])
        if proportion:
            total = lens[i] + lens[j]
            rf = rf / total if total else 1.0
        result[i, j] = result[j, i] = rf
    return DistanceMatrix(result, ids, validate=False)
