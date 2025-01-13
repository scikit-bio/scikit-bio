# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from heapq import heapify, heappop

import numpy as np

from skbio.tree import TreeNode
from ._c_me import (
    _preorder,
    _postorder,
    _insert_taxon,
    _avgdist_taxon,
    _bal_avgdist_taxon,
    _avgdist_d2_insert,
    _bal_avgdist_insert,
    _ols_lengths,
    _ols_lengths_d2,
    _bal_lengths,
    _ols_min_branch_d2,
    _bal_min_branch,
    _avgdist_matrix,
    _bal_avgdist_matrix,
    _avgdist_swap,
    _bal_avgdist_swap,
    _ols_all_swaps,
    _ols_corner_swaps,
    _bal_all_swaps,
    _bal_avgdist_insert_p,
)
from ._utils import _check_dm, _check_dm_tree


def gme(dm, neg_as_zero=True):
    r"""Perform greedy minimum evolution (GME) for phylogenetic reconstruction.

    .. versionadded:: 0.6.2

    .. versionchanged:: 0.6.3
        Computational efficiency significantly improved.

    Parameters
    ----------
    dm : skbio.DistanceMatrix
        Input distance matrix containing distances between taxa.
    neg_as_zero : bool, optional
        If True (default), convert negative branch lengths into zeros.

    Returns
    -------
    TreeNode
        Reconstructed phylogenetic tree.

    See Also
    --------
    bme
    nni
    nj

    Notes
    -----
    Greedy Minimum Evolution (GME) [1]_ is a distance-based algorithm for phylogenetic
    reconstruction. It utilizes the minimum evolution (ME) principle [2]_ for selecting
    a tree topology with the lowest sum of branch lengths, as calculated using an
    ordinary least squares (OLS) framework [3]_.

    GME is *O*\(*n*:sup:`2`) in time and *O*\(*n*) in space, making it a more scalable
    method than multiple alternatives (e.g., :func:`nj` and :func:`bme`). Therefore it
    is suitable for reconstructing very large phylogenetic trees.

    GME generates an unrooted tree with variable tip heights and potentially some
    negative branch lengths. The first taxon in the distance matrix is always placed as
    a child of the root node. See the notes of :func:`nj` for how to deal with unrooted
    trees and negative branch lengths.

    A GME-generated tree may be further improved by executing the FastNNI algorithm
    implemented in :func:`nni` (with ``balanced=False``).

    A similar but less scalable algorithm using the balanced instead of OLS framework
    is provided in :func:`bme`.

    The same methods underlying :func:`gme`, :func:`bme` and :func:`nni` were also
    provided by the software package FastME [4]_.

    .. note::
        These scikit-bio functions were implemented following the original paper [1]_.
        It is not guaranteed that they will precisely mirror FastME's output. Although
        in practices they usually generate identical or equally optimal phylogenetic
        trees as FastME does.

    References
    ----------
    .. [1] Desper, R., & Gascuel, O. (2002). Fast and accurate phylogeny reconstruction
       algorithms based on the minimum-evolution principle. J Comput Biol, 9(5),
       687-705.

    .. [2] Rzhetsky, A., & Nei, M. (1993). Theoretical foundation of the
       minimum-evolution method of phylogenetic inference. Mol Biol Evol, 10(5),
       1073-1095.

    .. [3] Cavalli-Sforza, L. L., & Edwards, A. W. (1967). Phylogenetic analysis.
       Models and estimation procedures. American journal of human genetics,
       19(3 Pt 1), 233.

    .. [4] Lefort, V., Desper, R., & Gascuel, O. (2015). FastME 2.0: a comprehensive,
       accurate, and fast distance-based phylogeny inference program. Mol Biol Evol,
       32(10), 2798-2800.

    Examples
    --------
    Define a new distance matrix object describing the distances between five taxa:
    human, monkey, pig, rat, and chicken.

    >>> from skbio import DistanceMatrix
    >>> from skbio.tree import gme

    >>> dm = DistanceMatrix([[0, 0.02,  0.18,  0.34,  0.55],
    ...                      [0.02,  0, 0.19, 0.35,  0.55],
    ...                      [0.18, 0.19,  0,  0.34,  0.54],
    ...                      [0.34, 0.35,  0.34,  0,  0.62],
    ...                      [0.55,  0.55,  0.54,  0.62,  0]],
    ...                     ['human','monkey','pig','rat','chicken'])

    Perform Greedy Minimum Evoltuion (GME) and construct the minimum evolution tree
    representing the relationship between those taxa. This is returned as a TreeNode
    object.

    >>> tree = gme(dm)
    >>> print(tree.ascii_art())
              /-monkey
             |
             |          /-pig
             |---------|
    ---------|         |          /-rat
             |          \--------|
             |                    \-chicken
             |
              \-human

    """
    _check_dm(dm)

    # reconstruct tree topology and branch lengths using GME
    tree, lens = _gme(dm.data)

    if neg_as_zero:
        lens[lens < 0] = 0

    return _to_treenode(tree, dm.ids, lens, unroot=True)


def bme(dm, neg_as_zero=True, **kwargs):
    r"""Perform balanced minimum evolution (BME) for phylogenetic reconstruction.

    .. versionadded:: 0.6.3

    Parameters
    ----------
    dm : skbio.DistanceMatrix
        Input distance matrix containing distances between taxa.
    neg_as_zero : bool, optional
        If True (default), convert negative branch lengths into zeros.

    Returns
    -------
    TreeNode
        Reconstructed phylogenetic tree.

    See Also
    --------
    gme
    nni
    nj

    Notes
    -----
    Balanced Minimum Evolution (BME) [1]_ is a refinement of the distance-based minimum
    evolution problem where the average distances between subtrees are independent of
    the sizes of the subtrees. This is referred to as a balanced (or simply BME)
    framework [2]_, as in contrast to the OLS framework used by GME (:func:`gme`).

    The BME algorithm implemented here uses a similar greedy algorithm as implemented
    in :func:`gme`, but less scalable due to the need to update subtree distances as
    the tree topology changes. The algorithm is sub-*O*\(*n*:sup:`3`) in time and
    *O*\(*n*:sup:`2`) in space.

    Refer to :func:`gme` for the format of the output tree and subsequent treatments.

    A BME-generated tree may be further improved by executing the BNNI algorithm
    implemented in :func:`nni` (with ``balanced=True``).

    The same method was provided by FastME [3]_. See :func:`gme` for notes on this.

    .. note::
        Experimental feature: Add ``parallel=True`` will enable parallelization,
        which may increase the performance of the algorithm. This feature may not
        be stable and may be modified without notice in the future.

    References
    ----------
    .. [1] Desper, R., & Gascuel, O. (2002). Fast and accurate phylogeny reconstruction
       algorithms based on the minimum-evolution principle. J Comput Biol, 9(5),
       687-705.

    .. [2] Pauplin, Y. (2000). Direct calculation of a tree length using a distance
       matrix. J Mol Evol, 51, 41-47.

    .. [3] Lefort, V., Desper, R., & Gascuel, O. (2015). FastME 2.0: a comprehensive,
       accurate, and fast distance-based phylogeny inference program. Mol Biol Evol,
       32(10), 2798-2800.

    Examples
    --------
    Define a new distance matrix object describing the distances between five taxa:
    human, monkey, pig, rat, and chicken.

    >>> from skbio import DistanceMatrix
    >>> from skbio.tree import bme

    >>> dm = DistanceMatrix([[0, 0.02,  0.18,  0.34,  0.55],
    ...                      [0.02,  0, 0.19, 0.35,  0.55],
    ...                      [0.18, 0.19,  0,  0.34,  0.54],
    ...                      [0.34, 0.35,  0.34,  0,  0.62],
    ...                      [0.55,  0.55,  0.54,  0.62,  0]],
    ...                     ['human','monkey','pig','rat','chicken'])

    Perform Balanced Minimum Evoltuion (BME) and construct the minimum evolution tree
    representing the relationship between those taxa. This is returned as a TreeNode
    object.

    >>> tree = bme(dm)
    >>> print(tree.ascii_art())
              /-monkey
             |
             |          /-pig
             |---------|
    ---------|         |          /-rat
             |          \--------|
             |                    \-chicken
             |
              \-human

    """
    _check_dm(dm)

    # reconstruct tree topology and branch lengths using BME
    tree, lens = _bme(dm.data, **kwargs)

    if neg_as_zero:
        lens[lens < 0] = 0

    return _to_treenode(tree, dm.ids, lens, unroot=True)


def nni(tree, dm, balanced=True, neg_as_zero=True):
    r"""Perform nearest neighbor interchange (NNI) to improve a phylogenetic tree.

    .. versionadded:: 0.6.2

    .. versionchanged:: 0.6.3
        Computational efficiency significantly improved.

    Parameters
    ----------
    tree : skbio.TreeNode
        Input phylogenetic tree to be rearranged.
    dm : skbio.DistanceMatrix
        Input distance matrix containing distances between taxa.
    balanced : bool, optional
        Use the OLS framework (False) or the balanced framework (True, default).
    neg_as_zero : bool, optional
        If True (default), convert negative branch lengths into zeros.

    Returns
    -------
    TreeNode
        Rearranged phylogenetic tree.

    Notes
    -----
    Nearest neighbor interchange (NNI) is a method for tree rearrangement aiming at
    optimizing a given tree topology according to certain criteria. It iteratively
    swaps neighboring branches that could result in a better tree, until no such swaps
    remain in the optimized tree.

    This function performs NNI for the minimum evolution (ME) problem on phylogenetic
    trees. The implementation is based on [1]_. Two versions of the method are
    provided: the FastNNI algorithm (``balanced=False``) based on an ordinary least
    squares (OLS) framework (see also :func:`gme`) and the BNNI algorithm
    (``balanced=True``) based on a balanced framework (see also :func:`bme`). The
    former is more scalable than the latter.

    The input tree may be rooted or unrooted, and it must be strictly bifurcating. One
    can apply :meth:`~TreeNode.is_bifurcating` or :meth:`~TreeNode.bifurcate` to check
    or make a bifurcating tree. The set of taxa in the tree (accessible via
    :meth:`~TreeNode.subset`) must match those in the distance matrix.

    The output tree is unrooted, with the first taxon in the distance matrix placed as
    a child of the root node. Branch lengths are assigned according to the framework of
    choice. See also :func:`nj` regarding rooting an unrooted tree and dealing with
    potential negative branch lengths.

    The same methods were provided by FastME [2]_. See :func:`gme` for notes on this.

    References
    ----------
    .. [1] Desper, R., & Gascuel, O. (2002). Fast and accurate phylogeny reconstruction
       algorithms based on the minimum-evolution principle. J Comput Biol, 9(5),
       687-705.

    .. [2] Lefort, V., Desper, R., & Gascuel, O. (2015). FastME 2.0: a comprehensive,
       accurate, and fast distance-based phylogeny inference program. Mol Biol Evol,
       32(10), 2798-2800.

    Examples
    --------
    >>> from skbio import DistanceMatrix, TreeNode
    >>> from skbio.tree import nni

    Define a distance matrix describing the distances between five taxa.

    >>> dm = DistanceMatrix([[0, 0.02,  0.18,  0.34,  0.55],
    ...                      [0.02,  0, 0.19, 0.35,  0.55],
    ...                      [0.18, 0.19,  0,  0.34,  0.54],
    ...                      [0.34, 0.35,  0.34,  0,  0.62],
    ...                      [0.55,  0.55,  0.54,  0.62,  0]],
    ...                     ['human','monkey','pig','rat','chicken'])

    Provide a tree to be rearranged. The tree must be bifurcating. Branches lengths are
    not required.

    >>> tree = TreeNode.read(["((human,chicken),(rat,monkey),pig);"])
    >>> print(tree.ascii_art())
                        /-human
              /--------|
             |          \-chicken
             |
    ---------|          /-rat
             |---------|
             |          \-monkey
             |
              \-pig

    Perform nearest neighbor interchange (NNI). This will rearrange tree topolgy to
    better approach the minimum evolution criterion.

    >>> tree = nni(tree, dm)
    >>> print(tree.ascii_art())
                        /-pig
              /--------|
             |         |          /-chicken
             |          \--------|
    ---------|                    \-rat
             |
             |--monkey
             |
              \-human

    Besides rearranging tree topology, estimated branch lengths are assigned to the
    tree.

    >>> rat = tree.find('rat')
    >>> print(rat.length)
    0.20875

    """
    _check_dm_tree(dm, tree)

    # generate tree array
    taxa = dm.ids
    tree, preodr, postodr = _root_from_treenode(tree, taxa)

    # allocate lengths
    lens = np.empty(len(tree), dtype=float)

    # perform BNNI or FastNNI
    func = _bnni if balanced else _fastnni
    func(dm.data, tree, preodr, postodr, lens)

    if neg_as_zero:
        lens[lens < 0] = 0

    # generate TreeNode object
    return _to_treenode(tree, taxa, lens, unroot=True)


def _gme(dm):
    r"""Perform greedy minimum evolution (GME) for phylogenetic reconstruction.

    Parameters
    ----------
    dm : ndarray of float of shape (m, m)
        Input distance matrix containing distances between taxa.

    Returns
    -------
    tree : ndarray of int of shape (2m - 3, 4)
        Reconstructed tree topology by GME.
    lens : ndarray of float of shape (2m - 3,)
        Estimated branch lengths by OLS.

    Notes
    -----
    The GME algorithm involves growing an unrooted tree by iteratively inserting taxa.
    In the beginning it initializes a triplet tree with the first three taxa:

          0
          |
          x
         / \
        1   2

    Taxon 0 is placed at the root position, with its unique descendant x serving as
    connector for taxa 1 and 2 as its left and right children. Next, taxon 3 will
    be inserted into a branch that suffices the minimum evolution criterion. Let's say
    the branch connecting x and 2. The tree becomes:

          0
          |
          x
         / \
        1   y
           / \
          2   3

    So on so forth until all taxa are included in the tree.

    Finally, branch lengths are estimated using an ordinary least squares (OLS)
    framework.

    ---

    This implementation uses an array-based tree structure to improve efficiency. See
    :func:`_check_tree` for details of this structure. Basically, the root node (taxon
    0) is omitted, making the tree strictly binary. Instead, the unique descendant (x)
    is now treated as the root and marked with taxon 0:

          0
         / \
        1   y
           / \
          2   3

    The root is always the first row in the tree array. Other nodes are appended to the
    tree array in the same order as they are added to the tree.

    ---

    For an input distance matrix with m taxa, the output tree will have 2m - 3 nodes.
    This memory space is pre-allocated to improve efficiency. The algorithm doesn't
    remake arrays during iteration, but repeatedly uses the already allocated space.
    During an interation involving x nodes, the first x positions of each array are
    occupied, while the remaining positions are disregarded. There is no need to reset
    array values after each iteration.

    Since they are created de novo, all arrays are C-continuous. This permits further
    optimization in the Cython code.

    """
    # number of taxa
    m = dm.shape[0]

    # number of nodes in the final tree
    n = 2 * m - 3

    # Pre-allocate memory spaces:
    # tree structure
    tree, preodr, postodr = _allocate_tree(n)

    # average distances between distant-2 subtrees
    ad2 = np.empty((n, 2), dtype=float)

    # average distances from a taxon to each subtree
    adk = np.empty((n, 2), dtype=float)

    # branch lengths or length changes
    lens = np.empty((n,), dtype=float)

    # Initialize 3-taxon tree.
    _init_tree(dm, tree, preodr, postodr, ad2, matrix=False)

    # Iteratively add taxa to the tree.
    for k in range(3, m):
        # Calculate average distances from new taxon to existing subtrees.
        _avgdist_taxon(adk, k, dm, tree, preodr, postodr)

        # Find a branch with minimum length change.
        target = _ols_min_branch_d2(lens, ad2, adk, tree, preodr)

        # Update average distances between distant-2 subtrees.
        _avgdist_d2_insert(ad2, target, adk, tree, preodr)

        # Insert new taxon into tree.
        _insert_taxon(k, target, tree, preodr, postodr, use_depth=False)

    # Calculate branch lengths using an OLS framework.
    _ols_lengths_d2(lens, ad2, tree)

    return tree, lens


def _bme(dm, parallel=False):
    r"""Perform balanced minimum evolution (BME) for phylogenetic reconstruction.

    Parameters
    ----------
    dm : ndarray of float of shape (m, m)
        Input distance matrix containing distances between taxa.

    Returns
    -------
    tree : ndarray of int of shape (2m - 3, 4)
        Reconstructed tree topology by BME.
    lens : ndarray of float of shape (2m - 3,)
        Estimated branch lengths by a balanced framework.

    Notes
    -----
    The BME algorithm resembles GME but uses a balanced instead of OLS framework.
    See :func:`_gme` for details of the latter.

    The main differences are: 1) BME requires a full matrix stored and iteratively
    updated, which is less efficient than GME in time and space. 2) BME doesn't rely
    on subtree sizes, which is simpler and more efficient than GME in individual
    calculations (but this cannot compensate for the former).

    """
    func = _bal_avgdist_insert_p if parallel else _bal_avgdist_insert

    # numbers of taxa and nodes in the tree
    m = dm.shape[0]
    n = 2 * m - 3

    # Pre-allocate memory spaces:
    # tree structure
    tree, preodr, postodr = _allocate_tree(n)

    # average distances between all subtrees
    # (This is the dominant factor in BME, and what makes it more expensive than GME.)
    adm = np.empty((n, n), dtype=float)

    # average distances from a taxon to each subtree
    adk = np.empty((n, 2), dtype=float)

    # branch lengths or length changes
    lens = np.empty((n,), dtype=float)

    # a stack for traversal operations
    stack = np.empty((n,), dtype=int)

    # initialize 3-taxon tree
    _init_tree(dm, tree, preodr, postodr, adm, matrix=True)

    # Pre-calculate negative powers of 2.
    powers = np.ldexp(1.0, -np.arange(m))

    # Iteratively add taxa to the tree.
    for k in range(3, m):
        # Calculate balanced average distances from new taxon to existing subtrees.
        _bal_avgdist_taxon(adk, k, dm, tree, preodr, postodr)

        # Find the branch with minimum length change.
        target = _bal_min_branch(lens, adm, adk, tree, preodr)

        # Update balanced average distance matrix between all subtrees.
        func(adm, target, adk, tree, preodr, postodr, powers, stack)

        # Insert new taxon into tree.
        _insert_taxon(k, target, tree, preodr, postodr, use_depth=True)

    # Calculate branch lengths using a balanced framework.
    _bal_lengths(lens, adm, tree)

    return tree, lens


def _fastnni(dm, tree, preodr, postodr, lens):
    r"""Perform fast nearest neighbor interchange (FastNNI) on a tree.

    To improve the tree under the minimum evolution (ME) criterion using an OLS
    framework on a distance matrix between taxa.

    This algorithm was implemented following Section 2.2 and Appendix 4 of Desper and
    Gascuel (2002). Basically, it creates a full average distance matrix between all
    pairs of subtrees, evaluates all possible swaps, then it iteratively perform
    swaps that lead to the maximum reduction of overall branch lengths. During this
    process, it needs to update the evaluation of the four corner branches after one
    branch is swapped.

    ---

    The swaps that can reduce overall branch length are stored in a heap. The paper
    suggests using a maxheap, but the following code uses Python's heapq which is a
    minheap. Therefore the length changes are negated (in contrast to `_bnni`). Each
    item in the heap is a tuple consisting of:

        0. Length change (negative; smaller is better)
        1. Corresponding node index
        2. Corresponding side (0: left, 1: right) child to be swapped with sibling

    For example, (-0.15, 7, 1) means that swapping the right child and the sibling of
    node 7 will reduce the overall branch length by 0.15.

    Meanwhile, the length change of each node is stored in `lens`. For root and tips,
    this value is 0. For internal nodes, this value is negative or 0. Only negative
    values are stored in the heap. 0 means this node is not useful.

    Additionally, the "side" information is stored in column 7 of the tree array.

    When a swap is updated, the new length change (if negative) is pushed into the
    heap, but the old length change won't be deleted from the heap (which would be
    inefficient as O(n)). Instead, a "lazy deletion" mechanism is adopted: When the
    negative-most swap is popped out of the heap, it is checked against the stored
    length change (in `lens`) and side (in `tree[:, 7]`). If they don't match, then
    it means that it is outdated and should be skipped. Otherwise, it is considered
    relevant and the swap will be performed.

    TODO: heapq is a Python module and may not be efficient. Tests showed that the
    code isn't significantly faster than a naive `np.argmax` on then entire `lens`.
    Consider optimization.

    """
    n = tree.shape[0]
    stack = np.empty(n, dtype=int)

    # Calculate average distances between all pairs of subtrees.
    adm = np.empty((n, n))
    _avgdist_matrix(adm, dm, tree, preodr, postodr)

    # Calculate length changes of all possible swaps.
    _ols_all_swaps(lens, tree, adm)

    # Create a heap that stores negative length changes and their nodes and sides.
    heap = [(lens[i], i, tree[i, 7]) for i in np.nonzero(lens)[0]]
    heapify(heap)

    # Iteratively swap branches until there is no more beneficial swap.
    while heap:
        # Pop the swap that reduces the most length.
        L, target, side = heappop(heap)

        # Skip if this swap is outdated.
        if L != lens[target] or side != tree[target, 7]:
            continue

        # Reset this swap.
        lens[target] = 0

        # Swap the branches in the tree.
        _swap_branches(target, side, tree, preodr, stack, use_depth=False)

        # Update average distances after swapping.
        _avgdist_swap(adm, target, side, tree)

        # Update length reductions of swaps of the four corner branches.
        _ols_corner_swaps(target, heap, lens, tree, adm)

    # Calculate branch lengths using an OLS framework.
    _ols_lengths(lens, adm, tree)


def _bnni(dm, tree, preodr, postodr, lens):
    r"""Perform balanced nearest neighbor interchange (BNNI) on a tree.

    To improve the tree under the balanced minimum evolution (BME) criterion given a
    distance matrix between taxa.

    Implemented according to Section 3.2 and Appendix 5 of Desper and Gascuel (2002).
    Basically, BNNI is similar to FastNNI (see `_fastnni`), but it needs to update
    more than the four corner branches after each swap, as balanced average distances
    are dependent on the topology of the tree.

    ---

    The original paper doesn't explain how to store and update positive swaps. Because
    a large proportion of swaps need to be re-evaluated during every iteration, one
    would expect that a heap structure as suggested for the FastNNI algorithm may not
    be efficient. The current code simply stores all swaps (positive or not) and runs
    `np.argmax` on all of them during each iteration.

    """
    n = tree.shape[0]
    stack = np.empty(n, dtype=int)

    # Calculate balanced average distances between all pairs of subtrees.
    adm = np.empty((n, n))
    _bal_avgdist_matrix(adm, dm, tree, preodr, postodr)

    # Pre-calculate negative powers of 2.
    powers = np.ldexp(1.0, -np.arange(dm.shape[0]))

    # Initialize branch swapping information.
    gains, sides, nodes = _init_swaps(tree)

    # Iteratively swap branches until there is no more beneficial swap.
    while True:
        # Calculate length reductions of all possible swaps.
        _bal_all_swaps(gains, sides, nodes, adm, tree)

        # Find the swap with the maximum length reduction, and stop if non-positive.
        branch = gains.argmax()
        if (gain := gains[branch]) <= 0:
            break
        side = sides[branch]
        target = nodes[branch]

        # Swap the branches in the tree.
        _swap_branches(target, side, tree, preodr, stack, use_depth=True)

        # Update balanced average distances after swapping.
        _bal_avgdist_swap(adm, target, side, tree, preodr, powers, stack)

    # Calculate branch lengths using a balanced framework.
    _bal_lengths(lens, adm, tree)


def _check_tree(tree, preodr, postodr):
    r"""Check the integrity of an array-based tree structure.

    The greedy algorithms implemented in this module use a special array-based tree
    structure to improve efficiency.

    In this 2-D array, rows represent nodes in the order of addition to the tree,
    with the root placed at row 0. There are eight columns:

        0. Left child index, or 0 if a tip.
        1. Right child index, or taxon index if a tip.
        2. Parent index, or 0 if root.
        3. Sibling index, or 0 if root.
        4. Size, i.e., number of taxa descending from the node, or 1 if a tip.
        5. Depth, i.e., number of branches constituting the path to root.
        6. Preorder index.
        7. Postorder index.

    Meanwhile, node indices in preorder and postorder are stored in two separate 1-D
    arrays to facilitate traversals.

    Note: Columns 2-7 and the two separate arrays can be derived from columns 0 and 1.
    They are explicitly stored because of their frequent usage during the calculation.
    Also, they are updated iteratively as the tree grows, which is more efficient than
    de novo calculation on an existing tree.

    ---

    For example, a tree with four taxa is like (see :func:`_gme`):

          0
         / \
        1   y
           / \
          2   3

    This tree can be represented by the following 2-D array:

        [[1, 2, 0, 0, 3, 0, 0, 4],
         [0, 1, 0, 2, 1, 1, 1, 0],
         [3, 4, 0, 1, 2, 1, 2, 3],
         [0, 2, 2, 4, 1, 2, 3, 1],
         [0, 3, 2, 3, 1, 2, 4, 2]]

    The separate pre- and postorder arrays are:

        [0, 1, 2, 3, 4]
        [1, 3, 4, 2, 0]

    ---

    This function ensures that a tree array is valid. It is for test purpose, but not
    used in the actual greedy algorithms.

    """
    assert tree.shape[1] == 8
    n = tree[0, 4] * 2 - 1
    assert (tree[:n, 6].argsort() == preodr[:n]).all()
    assert (tree[:n, 7].argsort() == postodr[:n]).all()

    for i in range(n):
        left, right, parent, sibling, size, depth, preidx, postidx = tree[i]
        if left == 0:
            assert i == 0 or right != 0
        else:
            assert tree[left, 2] == tree[right, 2] == i
            assert tree[left, 3] == right
            assert tree[right, 3] == left

    sizes = np.zeros((n,), dtype=int)
    for i in range(n):
        node = postodr[i]
        if tree[node, 0] == 0:
            sizes[node] = 1
        else:
            sizes[node] = sizes[tree[node, :2]].sum()
    assert (sizes == tree[:n, 4]).all()

    depths = np.zeros((n,), dtype=int)
    for i in range(1, n):
        node = preodr[i]
        depths[node] = depths[tree[node, 2]] + 1
    assert (depths == tree[:n, 5]).all()

    stack = np.zeros((n,), dtype=int)
    _preorder(exp := np.zeros((n,), dtype=int), tree, stack)
    assert (exp == preodr[:n]).all()
    _postorder(exp := np.zeros((n,), dtype=int), tree, stack)
    assert (exp == postodr[:n]).all()


def _allocate_tree(n):
    r"""Pre-allocate memory space for an array-based tree structure.

    Parameters
    ----------
    n : int
        Total number of nodes in the tree.

    Returns
    -------
    (n, 8) ndarray of int
        Tree structure.
    (n,) ndarray of int
        Nodes in preorder.
    (n,) ndarray of int
        Nodes in postorder.

    See Also
    --------
    _check_tree

    """
    tree = np.empty((n, 8), dtype=int)
    preodr = np.empty((n,), dtype=int)
    postodr = np.empty((n,), dtype=int)
    return tree, preodr, postodr


def _to_treenode(tree, taxa, lens=None, unroot=False):
    r"""Convert an array-based tree structure into a TreeNode object.

    Parameters
    ----------
    tree : (n, 2+) array_like of int
        Tree structure.
    taxa : list of str
        Taxon names in order.
    lens : (n,) array_like of float, optional
        Branch lengths.
    unroot : bool, optional
        If True, generate an unrooted tree in which the root node has three children,
        the last of which is the first taxon. Otherwise leave the first taxon at the
        root node.

    Returns
    -------
    TreeNode
        Converted TreeNode object.

    See Also
    --------
    TreeNode

    Notes
    -----
    This function only needs the first two columns (left child, right child / taxon
    name) of a tree array to generate a TreeNode object.

    """
    nodes = [TreeNode(taxa[0])] + [
        TreeNode(None if a else taxa[b]) for a, b in tree[1:, :2]
    ]
    if lens is not None:
        for node, L in zip(nodes, lens):
            node.length = L
    for node, (left, right, *_) in zip(nodes, tree):
        if left:
            node.extend([nodes[left], nodes[right]], uncache=False)
    tree = nodes[0]
    if unroot:
        tree.append(TreeNode(taxa[0], length=tree.length), uncache=False)
        tree.name = None
        tree.length = None
    return tree


def _from_treenode(obj, taxmap, tree, preodr, postodr, pos=0, depth=0):
    r"""Convert a TreeNode object into an array-based tree structure.

    Parameters
    ----------
    tree : TreeNode
        Input tree. It must be strictly bifurcating.
    taxmap : dict of str : int
        Mapping of taxon names to indices.
    tree : (n, 8) ndarray of int
        Tree structure.
    preodr : (n,) ndarray of int
        Nodes in preorder.
    postodr : (n,) ndarray of int
        Nodes in postorder.
    pos : int, optional
        Position (axis 0) in the tree array for new data.
    depth : int, optional
        Depth of the root node.

    Raises
    ------
    ValueError
        If tree is not strictly bifurcating.

    See Also
    --------
    _allocate_tree
    _root_from_treenode

    Notes
    -----
    This function fills in pre-allocated tree arrays. It doesn't return any value.

    Nodes in the resulting tree array are ordered by preorder traversal. This may not
    agree with the order they were added to the original tree, if that ever happened.
    Therefore, one shouldn't directly compare the original and generated tree arrays.
    However, the topologies represented by the two tree arrays should be identical.

    The function doesn't fill branch lengths, which aren't required by the downstream
    algorithms. Adding this functionality should be straightforward.

    Parameters ``pos`` and ``depth`` are typically set by :func:`_root_from_treenode`.
    One doesn't need to specify them if working with this function alone. ``depth``
    also impacts the position in ``postodr``. See ``_root_from_treenode`` for
    rationales.

    """
    # determine position in postorder
    post_pos = pos - depth

    # create root node
    tree[pos, 2] = 0
    tree[pos, 3] = 0
    tree[pos, 5] = depth

    # if root node is a tip, just wrap up and return
    if not obj.children:
        tree[pos, 0] = 0
        tree[pos, 1] = taxmap[obj.name]
        tree[pos, 4] = 1
        preodr[pos] = tree[pos, 6] = pos
        postodr[post_pos] = pos
        tree[pos, 7] = post_pos
        return

    # perform preorder traversal, index nodes, fill parent (2), depth (5), taxon (1)
    # and length
    obj.i = pos
    pos_1 = pos + 1
    for i, node in enumerate(obj.preorder(include_self=False)):
        node.i = (ni := pos_1 + i)
        tree[ni, 2] = (pi := node.parent.i)
        tree[ni, 5] = tree[pi, 5] + 1
        if not node.children:
            tree[ni, 0] = 0
            tree[ni, 1] = taxmap[node.name]

    # number of nodes in the current tree
    n = i + 2

    # fill preorder and indices
    end = pos + n
    preodr[pos:end] = tree[pos:end, 6] = np.arange(pos, end)

    # perform postorder traversal, fill children (0 and 1), sibling (3), and size (4)
    for i, node in enumerate(obj.postorder()):
        postodr[post_pos + i] = (ni := node.i)
        tree[ni, 7] = post_pos + i
        if not node.children:
            tree[ni, 4] = 1
        else:
            try:
                left, right = node.children
            except ValueError:
                raise ValueError("Tree is not strictly bifurcating.")
            li, ri = left.i, right.i
            tree[ni, 0] = tree[ri, 3] = li
            tree[ni, 1] = tree[li, 3] = ri
            tree[ni, 4] = tree[li, 4] + tree[ri, 4]
            del left.i
            del right.i

    del obj.i


def _root_from_treenode(obj, taxa):
    r"""Convert a TreeNode object into a tree array rooted at the first taxon.

    Parameters
    ----------
    tree : TreeNode
        Input tree. It must be strictly bifurcating.
    taxa : list of str
        Taxon names in order.
    tree : (n, 8) ndarray of int
        Tree structure.

    Returns
    -------
    (n, 8) ndarray of int
        Tree structure.
    (n,) ndarray of int
        Nodes in preorder.
    (n,) ndarray of int
        Nodes in postorder.

    Raises
    ------
    ValueError
        If tree is not strictly bifurcating.

    See Also
    --------
    _from_treenode
    TreeNode.unrooted_move

    Notes
    -----
    This function locates the first taxon in the tree, walks from it to the root, and
    rotate nodes along the path such that the first taxon becomes the new root in the
    generated tree array. For example (p: parent, c: clade):

          root           taxon
         /  | \           /  \
        c3 p2 c4         c1  p2
          /  \      =>      /  \
         c2  p1           c2  root
            /  \              /  \
        taxon  c1            c3   c4

    Note that the original parent node always becomes the right child of the current
    node. The first parent (p1) is gone and replaced by the taxon. The output array is
    always preordered.

    """
    errmsg = "Tree is not strictly bifurcating."
    taxmap = {taxon: i for i, taxon in enumerate(taxa)}

    m = len(taxa)  # number of taxa
    n = m * 2 - 3  # number of nodes

    # allocate tree arrays
    tree, preodr, postodr = _allocate_tree(n)

    # position of the current node
    pos = 0

    # depth of the current node
    depth = 0

    # position of the previous node
    prev_pos = 0

    # position of the other child of the previous node in the path;
    # it will become the sibling of the current node
    sib_pos = 0

    # number of tips under the current node
    size = m - 1

    n_1 = n - 1

    # helper to add the current node to the array
    def _add_current():
        tree[pos, 2] = prev_pos
        tree[pos, 4] = size
        tree[pos, 5] = depth
        tree[pos, 6] = pos
        preodr[pos] = pos
        tree[pos, 7] = (post_pos := n_1 - depth)
        postodr[post_pos] = pos
        tree[prev_pos, 1] = pos
        tree[sib_pos, 3] = pos
        tree[pos, 3] = sib_pos

    # iterate over ancestors of the first taxon
    curr = obj.find(taxa[0])
    parent = None
    while True:
        prev = curr
        curr = curr.parent
        parent = curr.parent
        siblings = [x for x in curr.children if x is not prev]

        if not siblings:  # singleton branch
            raise ValueError(errmsg)
        n_sibs = len(siblings)
        if n_sibs == 1:
            other = siblings[0]

            # Current node is a regular internal node of the tree (the typical case):
            # Add the current node, then add the other child of it.
            #      parent            prev
            #        |                 |
            #      curr     =>       curr
            #     /   \             /   \
            # other   prev      other   parent
            if parent is not None:
                _add_current()

                prev_pos = pos
                pos += 1

                # add the other child clade
                _from_treenode(other, taxmap, tree, preodr, postodr, pos, depth + 1)

                tree[pos - 1, 0] = pos  # curr's left child
                tree[pos, 2] = pos - 1  # parent (now missing sibling)

                sib_pos = pos
                size -= (k := tree[pos, 4])
                pos += k * 2 - 1

            # Current node is the root of a rooted tree (i.e., it has two children):
            # Don't add the current node, instead add the other child.
            #      curr             prev
            #     /   \     =>      /
            # other   prev      other
            else:
                _from_treenode(other, taxmap, tree, preodr, postodr, pos, depth)

                tree[pos, 2] = prev_pos
                tree[prev_pos, 1] = pos
                tree[sib_pos, 3] = pos
                tree[pos, 3] = sib_pos

                break

        elif n_sibs == 2:
            # Current node is the root of an unrooted tree (i.e., it has three
            # children): Add the current node, then add the other two children of it.
            #                         prev
            #       curr               |
            #     /  |  \    =>       curr
            # left right prev        /   \
            #                     left   right
            if parent is None:
                _add_current()
                curr_pos = pos
                pos += 1

                left, right = siblings
                _from_treenode(left, taxmap, tree, preodr, postodr, pos, depth + 1)
                tree[pos, 2] = curr_pos
                left_pos = pos
                tree[curr_pos, 0] = pos
                pos += tree[pos, 4] * 2 - 1

                _from_treenode(right, taxmap, tree, preodr, postodr, pos, depth + 1)
                tree[pos, 2] = curr_pos
                tree[curr_pos, 1] = pos
                tree[left_pos, 3] = pos
                tree[pos, 3] = left_pos

                break

            else:
                raise ValueError(errmsg)
        else:
            raise ValueError(errmsg)

        # move up one level
        depth += 1

    return tree, preodr, postodr


def _init_tree(dm, tree, preodr, postodr, ads, matrix=False):
    """Initialize tree (triplet with taxa 0, 1 and 2)."""
    # triplet tree
    tree[:3] = [
        [1, 2, 0, 0, 2, 0, 0, 2],  # 0: root
        [0, 1, 0, 2, 1, 1, 1, 0],  # 1: left child
        [0, 2, 0, 1, 1, 1, 2, 1],  # 2: right child
    ]

    # nodes in pre- and postorder
    preodr[:3] = [0, 1, 2]
    postodr[:3] = [1, 2, 0]

    # average distance matrix between all subtrees
    if matrix:
        ads[0, 1] = ads[1, 0] = dm[0, 1]
        ads[0, 2] = ads[2, 0] = dm[0, 2]
        ads[1, 2] = ads[2, 1] = dm[1, 2]

    # average distances between distant-2 subtrees
    else:
        ads[1, 0] = dm[0, 1]
        ads[2, 0] = dm[0, 2]
        ads[1, 1] = ads[2, 1] = dm[1, 2]
        ads[0, 0] = ads[0, 1] = 0


def _insert_taxon_treenode(taxon, target, tree):
    r"""Insert a taxon between a target node and its parent.

    This function resembles :func:`_insert_taxon` but operates on a TreeNode object. It
    is for test purpose. It is not used in the actual algorithms.

    """
    link = TreeNode()
    tip = TreeNode(taxon)
    parent = target.parent
    if parent is None:  # root branch
        link.extend(tree.children)
        tree.append(link)
        tree.append(tip)
    else:  # other branch
        is_left = target is parent.children[0]
        link.append(target)
        link.append(tip)
        parent.append(link)
        if is_left:
            parent.children = parent.children[::-1]


def _avgdist_matrix_naive(adm, dm, tree, postodr):
    r"""Calculate a matrix of average distances between all pairs of subtrees.

    This function produces the same result as :func:`_avgdist_matrix`. However, it
    calculates all subtree-to-subtree distances based on the original taxon-to-taxon
    distances, as discussed in Eq. 1 of Desper and Gascuel (2002). instead of adopting
    a recursive strategy (Eq. 2). Therefore, it is significantly slower.

        d(A, B) = \frac{1}{|A||B|} \sum_{i \in A, j \in B} d(i, j)

    This algorithm is implemented as a reference and for comparison purpose. It is not
    used in the actual GME algorithm.

    """
    n = tree[0, 4] * 2 - 1
    taxas = {}
    for node in postodr[:n]:
        left, right = tree[node, :2]
        if not left:
            taxas[node] = frozenset([right])
        else:
            taxas[node] = taxas[left] | taxas[right]

    full = taxas[0] | frozenset([0])

    for a in range(n):
        a_taxa = taxas[a]
        a_taxa_r = full - a_taxa
        # adm[a, a] = dm[list(a_taxa)][:, list(a_taxa_r)].mean()  # self distance
        for b in range(a + 1, n):
            # If a is an ancestor (proper superset) of b, flip a's taxon set as the
            # the upper subtree of a will be considered.
            a_taxa_ = a_taxa_r if a_taxa > (b_taxa := taxas[b]) else a_taxa
            adm[a, b] = adm[b, a] = dm[list(a_taxa_)][:, list(b_taxa)].mean()
            # Because this is postorder traversal, a cannot be a descendant (proper
            # subset) of b.


def _avgdist_taxon_naive(adk, taxon, dm, tree, postodr):
    """Calculate average distances between a new taxon and existing subtrees.

    This function produces the same result as :func:`_avgdist_taxon`, but it
    calculates all taxon-to-subtree distances based on the original taxon-to-taxon
    distances, as discussed in Eq. 1 of Desper and Gascuel (2002), instead of
    adopting a recursive strategy (Eq. 2). Therefore, it is significantly slower.

    This algorithm is implemented as a reference and for comparison purpose, but not
    used in the actual GME algorithm.

    """
    n = tree[0, 4] * 2 - 1
    taxas = {}
    for node in postodr[:n]:
        left, right = tree[node, :2]
        if not left:
            taxas[node] = frozenset([right])
        else:
            taxas[node] = taxas[left] | taxas[right]

    full = taxas[0] | frozenset([0])

    dk = dm[taxon]
    for node in range(n):
        taxa_lower = taxas[node]
        adk[node, 0] = dk[list(taxa_lower)].mean()
        taxa_upper = full - taxa_lower
        adk[node, 1] = dk[list(taxa_upper)].mean()


def _init_swaps(tree):
    """Initialize branch swapping information.

    It will create three 1-D arrays to store information of all internal branches of
    the tree:

    - `gains`: Overall tree length reduction (larger is better).
    - `sides`: Which side (left: 0, right: 1) of the target node should be swapped.
    - `nodes`: Corresponding node index of the branch.

    Meanwhile, column 7 of the tree array will be filled with the branch index of the
    each node, if it corresponds to one. This column was previous used to store
    postorder index, but that information is not needed during swapping.

    For a tree with n taxa, there should be n - 3 internal branches. The association
    between nodes and internal branches are fixed, no matter how tree is swapped.

    """
    # number of internal branches
    n = tree[0, 4] - 2
    gains = np.zeros((n,), dtype=float)
    sides = np.empty((n,), dtype=int)
    nodes = np.empty((n,), dtype=int)

    # identify all internal branches
    branch = 0
    for node in range(1, tree.shape[0]):
        if tree[node, 0]:
            tree[node, 7] = branch
            nodes[branch] = node
            branch += 1

    return gains, sides, nodes


def _swap_branches_treenode(node1, node2):
    """Swap two nodes in a TreeNode object.

    This function is for testing purpose.

    """
    parent1, parent2 = node1.parent, node2.parent
    parent1.append(node2, False)
    parent2.append(node1, False)


def _swap_branches(target, side, tree, preodr, stack, use_depth=True):
    r"""Swap one child of the target node with its sibling.

                 |                         |
               parent                    parent
               /   \                     /   \
           target  sibling    =>     target  child
            /  \                      /  \
        other  child              other  sibling

    Target must be an internal node. Root and tips are not allowed.

    This function currently doesn't handle postorder, as subsequent algorithms don't
    need this information.

    This function is specifically designed for nearest neighbor interchange (NNI). It
    may be generalized to swapping non-neighbor branches, but that is not currently
    implemented.

    """
    # This function can potentially optimized using Cython. See _insert_taxon.

    child = tree[target, side]
    other = tree[target, 1 - side]  # other child
    parent = tree[target, 2]
    sibling = tree[target, 3]

    c_size = tree[child, 4]
    o_size = tree[other, 4]
    s_size = tree[sibling, 4]

    # update connections
    tree[target, side] = sibling
    tree[target, 3] = child
    tree[target, 4] += s_size - c_size

    p_side = int(tree[parent, 0] == target)  # sibling side of parent
    tree[parent, p_side] = child

    tree[sibling, 2] = target
    tree[sibling, 3] = other

    tree[child, 2] = parent
    tree[child, 3] = target

    tree[other, 3] = sibling

    # locate the clades under the relevant nodes
    c_start = tree[child, 6]
    c_width = c_size * 2 - 1
    c_end = c_start + c_width
    c_clade = preodr[c_start:c_end]

    s_start = tree[sibling, 6]
    s_width = s_size * 2 - 1
    s_end = s_start + s_width
    s_clade = preodr[s_start:s_end]

    o_start = tree[other, 6]
    o_width = o_size * 2 - 1
    o_end = o_start + o_width
    o_clade = preodr[o_start:o_end]

    # update depth (if needed)
    if use_depth:
        tree[c_clade, 5] -= 1
        tree[s_clade, 5] += 1

    # update preorder (naive)
    # _preorder(preodr, tree, stack)
    # tree[:n, 6] = np.argsort(preodr[:n])

    # update postorder (naive)
    # _postorder(postodr, tree, stack)
    # tree[:n, 7] = np.argsort(postodr[:n])

    # update preorder
    # sibling is the left child of parent
    if p_side == 0:
        tree[target, 6] += c_width - s_width
        stack[:c_width] = c_clade
        stack[c_width] = target
        c1_width = c_width + 1

        # sibling_clade, target, child_clade, other_clade =>
        # child_clade, target, sibling_clade, other_clade
        if side == 0:
            tree[s_clade, 6] += c_width + 1
            tree[c_clade, 6] -= s_width + 1
            preodr[(sc1_start := s_start + c1_width) : c_end] = s_clade
            preodr[s_start:sc1_start] = stack[:c1_width]

        # sibling_clade, target, other_clade, child_clade =>
        # child_clade, target, other_clade, sibling_clade
        else:
            tree[s_clade, 6] += c_width + o_width + 1
            tree[c_clade, 6] -= s_width + o_width + 1
            tree[o_clade, 6] += c_width - s_width
            stack[c1_width : (c1o_width := c1_width + o_width)] = o_clade
            preodr[(sc1o_start := s_start + c1o_width) : c_end] = s_clade
            preodr[s_start:sc1o_start] = stack[:c1o_width]

    # sibling is the right child of parent
    else:
        stack[:s_width] = s_clade

        # target, child_clade, other_clade, sibling_clade =>
        # target, sibling_clade, other_clade, child_clade
        if side == 0:
            tree[c_clade, 6] += (so_width := s_width + o_width)
            tree[s_clade, 6] -= c_width + o_width
            tree[o_clade, 6] += s_width - c_width
            stack[s_width:so_width] = o_clade
            preodr[(cso_start := c_start + so_width) : s_end] = c_clade
            preodr[c_start:cso_start] = stack[:so_width]

        # target, other_clade, child_clade, sibling_clade =>
        # target, other_clade, sibling_clade, child_clade
        else:
            tree[c_clade, 6] += s_width
            tree[s_clade, 6] -= c_width
            preodr[(cs_start := c_start + s_width) : s_end] = c_clade
            preodr[c_start:cs_start] = stack[:s_width]
