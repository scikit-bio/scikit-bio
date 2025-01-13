# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.tree import TreeNode
from skbio.util._decorator import params_aliased
from skbio.util._warning import _warn_deprecated
from ._c_nj import nj_minq_cy
from ._utils import _check_dm


@params_aliased([("neg_as_zero", "disallow_negative_branch_length", "0.6.3", True)])
def nj(
    dm,
    neg_as_zero=True,
    result_constructor=None,
    inplace=False,
):
    r"""Perform neighbor joining (NJ) for phylogenetic reconstruction.

    Parameters
    ----------
    dm : skbio.DistanceMatrix
        Input distance matrix containing pairwise distances among taxa.
    neg_as_zero : bool, optional
        If True (default), convert negative branch lengths into zeros.
    result_constructor : function, optional
        Function to apply to construct the result object. This must take a
        newick-formatted string as input. Deprecated and to be removed in a
        future release.

        .. deprecated:: 0.6.3

    inplace : bool, optional
        If True, the input distance matrix will be manipulated in-place to reduce
        memory consumption, at the cost of losing the original data. Default is False.

        .. versionadded:: 0.6.3

    Returns
    -------
    TreeNode
        Reconstructed phylogenetic tree.

        .. versionchanged:: 0.6.3
            The NJ algorithm has been optimized. The output may be slightly different
            from the previous version in root placement, node ordering, and the numeric
            precision of branch lengths. However, the overall tree topology and branch
            lengths should remain the same.

    See Also
    --------
    bme
    upgma
    TreeNode.root_at_midpoint

    Notes
    -----
    Neighbor joining (NJ) was initially described by Saitou and Nei (1987) [1]_. It is
    a simple and efficient agglomerative clustering method that builds a phylogenetic
    tree based on a distance matrix.

    This function implements the canonical NJ method. The algorithm has been optimized
    to improve efficiency, but no heuristic was involved to further accelerate it at
    the cost of accuracy. Therefore, one is guaranteed to obtain an optimal tree. The
    algorithm is *O*\(*n*:sup:`3`) in time and *O*\(*n*:sup:`2`) in space.

    NJ creates an **unrooted** tree with varying tip heights. This contrasts UPGMA
    (:func:`upgma`), which always produces ultrametric trees. One may subsequently use
    choice of strategies such as midpoint rooting (:meth:`~TreeNode.root_at_midpoint`)
    or outgroup rooting (:meth:`~TreeNode.root_by_outgroup`) to convert the result into
    a rooted tree.

    NJ is most accurate when distances are **additive** -- the distance between two
    taxa in the matrix equals to the sum of branch lengths connecting them in the tree.
    When this assumption is violated, which is common in real studies, negative branch
    lengths may be produced, challenging interpretation and subsequent analyses. To
    address this issue, this function converts negative branch lengths into zeros by
    default, but this behavior can be disabled by setting ``neg_as_zero`` to False.

    Gascuel and Steel (2006) provide a detailed overview of neighbor joining in terms
    of its biological relevance and limitations [2]_. They proved that NJ is a greedy
    heuristic to the balanced minimum evolution (BME) problem. An alternative method
    to this problem is provided by :func:`bme`.

    The example presented here is derived from the Wikipedia page on neighbor joining
    [3]_.

    References
    ----------
    .. [1] Saitou, N., & Nei, M. (1987). The neighbor-joining method: a new method for
       reconstructing phylogenetic trees. Mol Biol Evol, 4(4), 406-425.

    .. [2] Gascuel, O., & Steel, M. (2006). Neighbor-joining revealed. Mol Biol Evol,
       23(11), 1997-2000.

    .. [3] http://en.wikipedia.org/wiki/Neighbour_joining

    Examples
    --------
    Define a new distance matrix object describing the distances between five
    taxa: a, b, c, d, and e.

    >>> from skbio import DistanceMatrix
    >>> from skbio.tree import nj

    >>> data = [[0,  5,  9,  9,  8],
    ...         [5,  0, 10, 10,  9],
    ...         [9, 10,  0,  8,  7],
    ...         [9, 10,  8,  0,  3],
    ...         [8,  9,  7,  3,  0]]
    >>> ids = list('abcde')
    >>> dm = DistanceMatrix(data, ids)

    Construct the neighbor joining tree representing the relationship between
    those taxa. This is returned as a TreeNode object.

    >>> tree = nj(dm)
    >>> print(tree.ascii_art())
                                  /-a
                        /--------|
              /--------|          \-b
             |         |
             |          \-c
    ---------|
             |--e
             |
              \-d

    """
    if result_constructor is not None:
        msg = (
            "`result_constructor` is deprecated and will be removed in a future "
            "release."
        )
        _warn_deprecated(nj, "0.6.3", msg)

    _check_dm(dm)
    taxa = list(dm.ids)

    dm_ = dm.data
    if not inplace:
        dm_ = dm_.copy()

    lm = _nj(dm_)

    tree = _tree_from_linkmat(lm, taxa, rooted=False, neg_as_zero=neg_as_zero)
    if result_constructor is not None:
        tree = result_constructor(str(tree))
    return tree


def _nj(dm):
    r"""Perform neighbor joining (NJ) for phylogenetic reconstruction.

    Parameters
    ----------
    dm : (N, N) ndarray
        Input distance matrix.

    Returns
    -------
    ndarray of shape (N - 1, 4)
        Linkage matrix representing the tree.

    Notes
    -----
    This function manipulates the distance matrix in-place. Therefore, one should make
    a copy prior to running the function if the original distance matrix needs to be
    preserved.

    """
    # This function re-uses the original array space during iteration, without creating
    # additional intermediate arrays. Therefore, it is memory-efficient, and avoids the
    # time overhead of allocating new memory space.

    # This function only operates on arrays of numbers, therefore it can be further
    # Cythonized. However, Cythonization did not bring significant performance gain in
    # tests, likely because all operations already utilize NumPy APIs. That being said,
    # further optimization and testing should be convenient.

    N = n = dm.shape[0]  # dimension
    sums = dm.sum(axis=0)  # distance sums
    idxs = np.arange(N)  # cluster indices
    lm = np.empty((N - 1, 4))  # linkage matrix

    # Iteratively merge taxa until there are three left.
    while n > 3:
        # Create memory views of currently relevant array areas.
        dm_ = dm[:n, :n]
        sums_ = sums[:n]
        idxs_ = idxs[:n]

        # Find the minimum value of the Q-matrix and return its position (i, j).
        #   Q(i, j) = (n - 2) d(i, j) - \sum d(i) - \sum d(j)
        # The function call avoids constructing the entire Q-matrix, but instead
        # computes values and finds the minimum as the computation goes.
        i, j = nj_minq_cy(dm_, sums_)

        # Get half of the original distance at (i, j).
        d_ij_ = dm[i, j] / 2

        # Taxa i and j will be merged into a cluster {i, j}. The updated distance from
        # cluster to any other taxon k is:
        #   d({i, j}, k) = (d(i, k) + d(j, k) - d(i, j)) / 2
        # We first compute (d(i, k) + d(j, k)) / 2 and save the results in row i.
        dm_[i] += dm_[j]
        dm_[i] /= 2

        # Compute branch lengths of taxa i and j.
        #   \delta = (\sum d(i) - \sum d(j)) / (2(n - 2))
        #   L(i) = d(i, j) / 2 + \delta
        #   L(j) = d(i, j) / 2 - \delta
        delta_ = (sums_[i] - sums_[j]) / (2 * n - 4)
        L_i = d_ij_ + delta_
        L_j = d_ij_ - delta_

        # The previously calculated sums can be updated for re-use. Specifically, for
        # taxon k, there is:
        #   new sum = old sum - d(i, k) - d(j, k) + d({i,j}, k)
        #           = old sum - d(i, k) - d(j, k) + (d(i, k) + d(j, k) - d(i, j)) / 2
        #           = old sum - (d(i, k) + d(j, k)) / 2 - d(i, j) / 2
        # We already have (d(i, k) + d(j, k)) / 2 stored in row i, therefore:
        sums_[:] -= dm_[i]
        sums_[:] -= d_ij_

        # Now complete the calculation of the updated distances d({i, j}, k).
        dm_[i] -= d_ij_

        # Update column i to match row i.
        dm_[:, i] = dm_[i]

        # Because two taxa have been merged into one cluster, we will shrink the
        # distance matrix from (n, n) to (n - 1, n - 1). Specifically, we will move
        # the last row/column (index: n - 1) to row/column j.
        n_1 = n - 1
        dm_[j] = dm_[n_1]
        dm_[:, j] = dm_[:, n_1]

        # Also move the last sum to j.
        sums_[j] = sums_[n_1]

        # Then calculate the updated sum at i (now cluster {i, j}), which is the sum
        # of the updated distances.
        sums_[i] = dm_[i, :n_1].sum()

        # Store the taxa and branch lengths to the linkage matrix.
        lm[N - n] = idxs_[i], idxs_[j], L_i, L_j

        # Update cluster indices. Specifically, position i will have the new cluster
        # index. Meanwhile, position j will be replaced with the last cluster.
        idxs_[i] = 2 * N - n
        idxs_[j] = idxs_[n_1]

        n -= 1

    # Perform final calculation on the three remaining taxa. They will become children
    # of the root node, and the entire tree is unrooted.
    L_0 = (dm[0, 1] + dm[0, 2] - dm[1, 2]) / 2
    lm[N - 3] = idxs[1], idxs[2], dm[0, 1] - L_0, dm[0, 2] - L_0
    lm[N - 2] = idxs[0], 2 * N - 3, L_0, 0

    return lm


def _tree_from_linkmat(lm, taxa, rooted=True, neg_as_zero=True):
    r"""Convert a linkage matrix into a tree structure.

    Parameters
    ----------
    lm : (N, 4) array_like
        Linkage matrix.
    taxa : list of str
        Taxon names.
    rooted : bool, optional
        If True (default), will generate a rooted binary tree, with two children
        attached to the root node. If False, will generate an unrooted tree with
        three children attached to the root node.
    neg_as_zero : bool, optional
        Convert negative branch lengths to zero.

    Returns
    -------
    TreeNode
        Converted phylogenetic tree.

    See Also
    --------
    TreeNode.from_linkage_matrix

    Notes
    -----
    The linkage matrix data structure resembles that used in SciPy's hierarchical
    clustering, but also stores the branch lengths of both children, which could be
    unequal in a phylogenetic tree.

    In SciPy's linkage matrix, the four elements per row are (see:
    https://stackoverflow.com/questions/9838861/):

        cluster1, cluster2, length1+2, # original taxa

    In the current data structure, they are:

        cluster1, cluster2, length1, length2

    """
    # allocate node list
    nodes = [TreeNode(name=x) for x in taxa] + [TreeNode() for _ in range(len(lm))]

    # build tree incrementally
    idx = len(taxa)
    for c1, c2, L1, L2 in lm if rooted else lm[:-2]:
        c1, c2 = nodes[int(c1)], nodes[int(c2)]
        if neg_as_zero:
            c1.length = L1 if L1 >= 0 else 0.0
            c2.length = L2 if L2 >= 0 else 0.0
        else:
            c1.length, c2.length = L1, L2
        nodes[idx].extend([c1, c2], uncache=False)
        idx += 1

    # final treatment of an unroot tree (root is trifurcating)
    if not rooted:
        # this code assumes that the first element of the last row is the node
        # see the end of _nj
        (c0, _, L0, _), (c1, c2, L1, L2) = lm[-1], lm[-2]
        c0, c1, c2 = nodes[int(c0)], nodes[int(c1)], nodes[int(c2)]
        if neg_as_zero:
            c0.length = L0 if L0 >= 0 else 0.0
            c1.length = L1 if L1 >= 0 else 0.0
            c2.length = L2 if L2 >= 0 else 0.0
        else:
            c0.length, c1.length, c2.length = L0, L1, L2
        nodes[-1].extend([c0, c1, c2], uncache=False)
    return nodes[-1]
