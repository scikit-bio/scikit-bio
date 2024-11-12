# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io

import numpy as np

from skbio.stats.distance import DistanceMatrix
from skbio.tree import TreeNode
from skbio.util._warning import _warn_deprecated
from ._cutils import nj_minq_cy


def nj(
    dm,
    clip_to_zero=True,
    result_constructor=None,
    disallow_negative_branch_length=None,
):
    r"""Perform neighbor joining (NJ) for phylogenetic reconstruction.

    Parameters
    ----------
    dm : skbio.DistanceMatrix
        Input distance matrix containing pairwise distances among taxa.
    clip_to_zero : bool, optional
        If True (default), convert negative branch lengths into zeros.

        .. versionadded:: 0.6.3

    result_constructor : function, optional
        Function to apply to construct the result object. This must take a
        newick-formatted string as input. Deprecated and to be removed in a
        future release.

        .. deprecated:: 0.6.3

    disallow_negative_branch_length : bool, optional
        Alias of ``clip_to_zero`` for backward compatibility. Deprecated and to be
        removed in a future release.

        .. deprecated:: 0.6.3

    Returns
    -------
    TreeNode
        Reconstructed phylogenetic tree.

        .. versionchanged:: 0.6.3
            The NJ algorithm has been optimized. The output may be slightly different
            from the previous one in root placement, node ordering, and the numeric
            precision of branch lengths. However, the overall tree topology and branch
            lengths should remain the same.

    See Also
    --------
    bme
    upgma
    TreeNode.root_at_midpoint

    Notes
    -----
    Neighbor joining (NJ) was initially described in Saitou and Nei (1987) [1]_. The
    example presented here is derived from the Wikipedia page on neighbor joining [2]_.
    Gascuel and Steel (2006) provide a detailed overview of neighbor joining in terms
    of its biological relevance and limitations [3]_.

    Neighbor joining, by definition, creates unrooted trees with varying tip heights,
    which contrasts UPGMA (:func:`upgma`). One strategy for rooting the resulting tree
    is midpoint rooting, which is accessible as :meth:`TreeNode.root_at_midpoint`.

    Note that the tree constructed using neighbor joining is not rooted at a tip,
    unlike minimum evolution (:func:`bme`), so re-rooting is required before tree
    re-arrangement operations such as nearest neighbor interchange (NNI) (:func:`nni`)
    can be performed.

    References
    ----------
    .. [1] Saitou N, and Nei M. (1987) "The neighbor-joining method: a new
       method for reconstructing phylogenetic trees." Molecular Biology and
       Evolution. PMID: 3447015.

    .. [2] http://en.wikipedia.org/wiki/Neighbour_joining

    .. [3] Gascuel O, and Steel M. (2006) "Neighbor-Joining Revealed" Molecular
       Biology and Evolution, Volume 23, Issue 11, November 2006,
       Pages 1997-2000, https://doi.org/10.1093/molbev/msl072

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
             |          \-b
             |
    ---------|--c
             |
             |          /-d
              \--------|
                        \-e

    """
    # @deprecated
    if disallow_negative_branch_length is not None:
        _warn_deprecated(
            nj,
            "0.6.3",
            msg=(
                "`disallow_negative_branch_length` has been renamed as `clip_to_zero`."
                "The old name will be removed in a future release."
            ),
        )
        clip_to_zero = disallow_negative_branch_length
    if result_constructor is not None:
        _warn_deprecated(
            nj,
            "0.6.3",
            msg=(
                "`result_constructor` is deprecated and will be removed in a future "
                "release."
            ),
        )

    if dm.shape[0] < 3:
        raise ValueError(
            "Distance matrix must be at least 3x3 to generate a neighbor joining tree."
        )
    taxa = list(dm.ids)
    dm_ = dm.data.astype(float)  # make a copy of the distance matrix data
    lm = _nj(dm_)
    tree = _tree_from_linkmat(lm, taxa, rooted=False, clip_to_zero=clip_to_zero)
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
    # Further optimization and testing should be convenient.

    N = n = dm.shape[0]  # dimension
    sums = dm.sum(axis=0)  # distance sums
    idxs = np.arange(N)  # cluster indices
    lm = np.empty((N - 1, 4))  # linkage matrix

    # Iteratively merge taxa until there are three left.
    while n > 3:
        # Find the minimum value of the Q-matrix and return its position (i, j).
        #   Q(i, j) = (n - 2) d(i, j) - \sum d(i) - \sum d(j)
        # The function call avoids constructing the entire Q-matrix, but instead
        # computes values and finds the minimum as the computation goes.
        i, j = nj_minq_cy(dm[:n, :n], sums[:n])

        # Get half of the original distance at (i, j).
        d_ij_ = dm[i, j] / 2

        # Taxa i and j will be merged into a cluster, therefore the values in rows
        # and columns i and j are no longer needed, and their memory space can be
        # utilized.

        # The updated distance from cluster {i, j} to any other taxon k is:
        #   d({i, j}, k) = (d(i, k) + d(j, k) - d(i, j)) / 2
        # We first compute (d(i, k) + d(j, k)) / 2 and save the results in row i.
        dm[i, :n] += dm[j, :n]
        dm[i, :n] /= 2

        # Compute branch lengths of taxa i and j.
        #   \delta = (\sum d(i) - \sum d(j)) / (2(n - 2))
        #   L(i) = d(i, j) / 2 + \delta
        #   L(j) = d(i, j) / 2 - \delta
        delta_ = (sums[i] - sums[j]) / (2 * n - 4)
        L_i = d_ij_ + delta_
        L_j = d_ij_ - delta_

        # The previously calculated sums can be re-used, but they need to be updated.
        # For taxon k, there is:
        #   new sum = old sum - d(i, k) - d(j, k) + d({i, j}, k)
        #           = old sum - d(i, k) - d(j, k) + (d(i, k) + d(j, k) - d(i, j)) / 2
        #           = old sum - (d(i, k) + d(j, k)) / 2 - d(i, j) / 2
        # We already have (d(i, k) + d(j, k)) / 2 stored in row i, therefore:
        sums[:n] -= dm[i, :n]
        sums[:n] -= d_ij_

        # Now complete the calculation of the updated distances d({i, j}, k).
        dm[i, :n] -= d_ij_

        # Because two taxa have been merged into one cluster, we will shrink the
        # distance matrix from (n, n) to (n - 1, n - 1). Specifically, we will delete
        # row/column j, and update row/column i.

        # Every row/column beyond j will be moved -1 position. This involves moving
        # the right block leftward, moving the bottom block upward, and moving the
        # bottom-right block upleftward.
        dm[j : n - 1, :j] = dm[j + 1 : n, :j]
        dm[:j, j : n - 1] = dm[:j, j + 1 : n]
        dm[j : n - 1, j : n - 1] = dm[j + 1 : n, j + 1 : n]

        # Also move the sums beyond j leftward.
        sums[j : n - 1] = sums[j + 1 : n]

        # Then update row/column i. Because the updated distances were already stored
        # in row i, we only need to update column i.
        dm[: n - 1, i] = dm[i, : n - 1]

        # Then calculate the updated sum at i (now cluster {i, j}), which is the sum
        # of updated distances.
        sums[i] = dm[i, : n - 1].sum()

        # Store the taxa and branch lengths to the linkage matrix.
        lm[N - n] = idxs[i], idxs[j], L_i, L_j

        # Update cluster indices. Specifically, position i will have the new cluster
        # index. Meanwhile, indices beyond j will be moved leftward.
        idxs[i] = 2 * N - n
        idxs[j : n - 1] = idxs[j + 1 : n]

        n -= 1

    # Perform final calculation on the three remaining taxa. They will become children
    # of the root node, and the entire tree is unrooted.
    L_0 = (dm[0, 1] + dm[0, 2] - dm[1, 2]) / 2
    lm[N - 3] = idxs[1], idxs[2], dm[0, 1] - L_0, dm[0, 2] - L_0
    lm[N - 2] = idxs[0], 2 * N - 3, L_0, 0

    return lm


def _tree_from_linkmat(lm, taxa, rooted=True, clip_to_zero=True):
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
    clip_to_zero : bool, optional
        Convert negative branch lengths to zero.

    Returns
    -------
    TreeNode
        Converted phylogenetic tree.

    Notes
    -----
    The linkage matrix data structure resembles that used in SciPy's hierarchical
    clustering, but also stores the branch lengths of both children, which could be
    unequal in a phylogenetic tree.

    In SciPy's linkage matrix, the four elements per row are (see:
    https://stackoverflow.com/questions/9838861/):

        taxon1, taxon2, length1+2, # original taxa

    In the current data structure, they are:

        taxon1, taxon2, length1, length2

    """
    nodes = [TreeNode(name=x) for x in taxa] + [TreeNode() for _ in range(len(lm))]
    idx = len(taxa)
    for c1, c2, l1, l2 in lm:
        c1, c2 = nodes[int(c1)], nodes[int(c2)]
        if clip_to_zero:
            c1.length = l1 if l1 >= 0 else 0.0
            c2.length = l2 if l2 >= 0 else 0.0
        else:
            c1.length, c2.length = l1, l2
        nodes[idx].extend([c1, c2], uncache=False)
        idx += 1
    tree = nodes[-1]
    if not rooted:
        tree.unroot(uncache=False)
    return tree
