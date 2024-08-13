# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import heapq as hq

from skbio.tree._gme import _average_distance_matrix, _ols_edge
from skbio.tree._bme import _balanced_average_matrix, _bal_ols_edge


def nni(tree, dm, allow_edge_estimation=True, inplace=True, balanced=True):
    r"""Perform nearest neighbor interchange (NNI) on a phylogenetic tree.

    Parameters
    ----------
    tree : skbio.TreeNode
        Input phylogenetic tree to be rearranged.
    dm : skbio.DistanceMatrix
        Input distance matrix containing distances between taxa.
    allow_edge_estimation : bool, optional
        Whether to perform an OLS-based estimation of edge values (``True``, default)
        or return a tree without edge values assigned (``False``).
    inplace : bool, optional
        Whether manipulate the tree in place (``True``, default) or return a
        copy of the tree (``False``).
    balanced : bool, optional
        Whether to use the minimum evolution framework or the balanced
        minimum evolution framework. The definition of average distance
        between subtrees either ignores subtree size (``True``, default)
        or calculates based on subtree size (``False``).

    Returns
    -------
    TreeNode
        Rearranged phylogenetic tree (if ``inplace`` is ``True``).

    Notes
    -----
    NNI algorithm for minimum evolution problem on phylogenetic trees. It rearranges
    an initial tree topology by performing subtree exchanges such that the distance
    is minimized. This implementation is based on the FastNNI algorithm and the
    BNNI algorithm (BNNI)[1]_.

    The two versions of NNI are due to the relationship with minimum evolution
    and the use of average distances between subtrees. FastNNI is based on the
    Minimum Evolution (ME) problem while BNNI is based on the Balanced Minimum
    Evolution (BME) problem. In BME the sizes of subtrees are ignored while ME
    considers the size of subtrees in the calculation of average distance.

    For both versions of NNI, the input tree is required to be binary and rooted
    at a leaf node such that there is a unique descendant from the root.

    References
    ----------
    .. [1] Desper R, Gascuel O. Fast and accurate phylogeny reconstruction
       algorithms based on the minimum-evolution principle. J Comput Biol.
       2002;9(5):687-705. doi: 10.1089/106652702761034136. PMID: 12487758.

    Examples
    --------
    Define a new distance matrix object describing the distances between five
    taxa: human, monkey, pig, rat, and chicken.

    >>> from skbio import DistanceMatrix
    >>> from skbio.tree import nni

    >>> dm = DistanceMatrix([[0, 0.02,  0.18,  0.34,  0.55],
    ...                      [0.02,  0, 0.19, 0.35,  0.55],
    ...                      [0.18, 0.19,  0,  0.34,  0.54],
    ...                      [0.34, 0.35,  0.34,  0,  0.62],
    ...                      [0.55,  0.55,  0.54,  0.62,  0]],
    ...                      ['human','monkey','pig','rat','chicken'])

    Also, provide a tree topology to be rearranged. The tree provided is
    required to be a binary tree rooted at a leaf node.

    Note that the tree provided does not require to have assigned edge lengths.

    >>> from skbio.tree import TreeNode

    >>> tree = TreeNode.read(["(((human,chicken),(rat,monkey)))pig;"])
    >>> print(tree.ascii_art())
                                  /-human
                        /--------|
                       |          \-chicken
    -pig----- /--------|
                       |          /-rat
                        \--------|
                                  \-monkey

    Perform nearest neighbor interchange (NNI), here the BNNI version is used.
    By default, the tree is rearrangede in place.

    >>> nni(tree, dm)
    >>> print(tree.ascii_art())
                                  /-rat
                        /--------|
                       |          \-chicken
    -pig----- /--------|
                       |          /-monkey
                        \--------|
                                  \-human

    Besides rearranging the tree, estimated edge lengths are assigned to the
    tree.

    >>> rat = tree.find('rat')
    >>> print(rat.length)
    0.21

    """
    # Initialize and populate the average distance matrix
    if not inplace:
        tree = tree.copy()
    if len(tree.root().children) != 1:
        raise TypeError(
            "Could not perform NNI. " "Tree needs to be rooted at a leaf node."
        )
    for node in tree.non_tips():
        if len(node.children) != 2:
            raise TypeError("Could not perform NNI. Tree needs to be a binary tree.")
    if balanced:
        adm = _balanced_average_matrix(tree, dm)
    else:
        adm = _average_distance_matrix(tree, dm)
    while True:
        # create heap of possible swaps and then swapping subtrees
        # until no more swaps are possible.
        adm = _average_distance_matrix(tree, dm)
        if balanced:
            heap = _swap_heap(tree, adm)
        else:
            heap = _swap_heap(tree, adm, balanced=False)
        if not heap:
            break
        swap = hq.heappop(heap)
        _perform_swap(swap[1][0], swap[1][1])
    # edge values are added using an OLS framework.
    if allow_edge_estimation:
        if balanced:
            _bal_ols_edge(tree, dm)
        else:
            _ols_edge(tree, dm)
    if not inplace:
        return tree


def _perform_swap(node1, node2):
    """Return a tree after swapping two subtrees."""
    parent1, parent2 = node1.parent, node2.parent
    parent1.append(node2)
    parent2.append(node1)


def _swap_length(a, b, c, d, i, j, k, m, adm):
    """Return the change in overall tree length after a given swap.

    The count of leaves contained in each subtree are denoted 'a, b, c, d' while
    each node defining the subtree has the index 'i, j, k, m', respectively.

    """
    lambda1 = (a * d + b * c) / ((a + b) * (c + d))
    lambda2 = (a * d + b * c) / ((a + c) * (b + d))
    return 0.5 * (
        (lambda1 - 1) * (adm[i][k] + adm[j][m])
        - (lambda2 - 1) * (adm[i][j] + adm[k][m])
        - (lambda1 - lambda2) * (adm[i][m] + adm[j][k])
    )


def _balanced_swap_length(i, j, k, m, adm):
    """Return the change in overall tree length after a given swap.

    Uses the definition of average distance from the balanced minimum evolution
    problem, and only requires a node's index while ignoring subtree size.

    """
    return 0.25 * ((adm[i][j] + adm[k][m]) - (adm[i][k] + adm[j][m]))


def _swap_heap(tree, adm, balanced=True):
    """Return a maxheap ordered by the swap length for all possible swaps.

    The 'balanced' option uses the balanced definition (``True``, default)
    or the classical definition (``False``).

    """
    heap = []
    ordered = list(tree.postorder(include_self=False))
    root = tree.root()
    n_taxa = root.count(tips=True) + 1
    # begin by finding nodes which are the child node of an internal edge
    for node in ordered:
        # ignore tips of the tree
        if node.is_tip():
            continue
        # identify the parent and grandparent nodes
        parent = node
        a = parent.parent
        # identify the index of each neighboring node
        for index, node in enumerate(ordered):
            if node == a:
                i1 = index
        for child in parent.children:
            if child.is_tip():
                continue
            childnode = child
            c, d = childnode.children
            for sibling in childnode.siblings():
                b = sibling
            for index, node in enumerate(ordered):
                if node == b:
                    i2 = index
                elif node == c:
                    i3 = index
                elif node == d:
                    i4 = index
            # count the tips of the subtrees defined by the neighboring nodes
            sub_tips = []
            for subtree in [b, c, d]:
                sub_tips.append(1 if subtree.is_tip() else subtree.count(tips=True))
            b_, c_, d_ = sub_tips
            a_ = n_taxa - b_ - c_ - d_
            # calculate the swap length for the two possible swaps given the edge
            if balanced:
                swap_1 = _swap_length(a_, b_, c_, d_, i1, i2, i3, i4, adm)
                swap_2 = _swap_length(a_, b_, d_, c_, i1, i2, i4, i3, adm)
            else:
                swap_1 = _balanced_swap_length(i1, i2, i3, i4, adm)
                swap_2 = _balanced_swap_length(i1, i2, i4, i3, adm)
            # store the best possible swap into a maxheap
            if swap_1 > swap_2 and swap_1 > 0:
                swap = -1 * swap_1
                hq.heappush(heap, (swap, (b, c)))
            elif swap_2 > swap_1 and swap_2 > 0:
                swap = -1 * swap_2
                hq.heappush(heap, (swap, (b, d)))
    return heap
