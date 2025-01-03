# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from skbio.tree import TreeNode
from ._c_me import (
    _preorder,
    _postorder,
    _avgdist_taxon,
    _bal_avgdist_taxon,
    _avgdist_d2_insert,
    _bal_avgdist_insert,
    _ols_lengths_d2,
    _bal_lengths,
    _ols_min_branch_d2,
    _bal_min_branch,
)
from ._utils import _check_dm


def gme(dm):
    r"""Perform greedy minimum evolution (GME) for phylogenetic reconstruction.

    .. versionadded:: 0.6.2

    Parameters
    ----------
    dm : skbio.DistanceMatrix
        Input distance matrix containing distances between taxa.

    Returns
    -------
    TreeNode
        Reconstructed phylogenetic tree.

    See Also
    --------
    bme
    nj

    Notes
    -----
    Greedy Minimum Evolution (GME) is a distance-based algorithm for phylogenetic
    reconstruction utilizing the minimum evolution principle for selecting a tree
    topology with the lowest sum of branch lengths according to a given method of
    estimating branch lengths. Ordinary Least Squares (OLS) is a natural framework for
    edge estimation as it is statistically consistent with minimum evolution and is
    used for GME.

    References
    ----------
    .. [1] Desper, R., & Gascuel, O. (2002). Fast and accurate phylogeny reconstruction
       algorithms based on the minimum-evolution principle. J Comput Biol, 9(5),
       687-705.

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

    return _to_treenode(tree, dm.ids, lens, unroot=True)


def bme(dm):
    r"""Perform balanced minimum evolution (BME) for phylogenetic reconstruction.

    .. versionadded:: 0.6.3

    Parameters
    ----------
    dm : skbio.DistanceMatrix
        Input distance matrix containing distances between taxa.

    Returns
    -------
    TreeNode
        Reconstructed phylogenetic tree.

    See Also
    --------
    gme
    nj

    Notes
    -----
    Balanced Minimum Evolution (BME) is a refinement of the distance-based minimum
    evolution problem where average distances between subtrees ignores the size of the
    subtrees. The BME algorithm implemented here uses the same OLS based edge
    estimation used with Greedy Minimum Evolution (GME).

    References
    ----------
    .. [1] Desper, R., & Gascuel, O. (2002). Fast and accurate phylogeny reconstruction
       algorithms based on the minimum-evolution principle. J Comput Biol, 9(5),
       687-705.

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
    tree, lens = _bme(dm.data)

    return _to_treenode(tree, dm.ids, lens, unroot=True)


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

    For an input distance matrix with m taxa, the output tree will have 2m - 3 nodes.
    This memory space is pre-allocated to improve efficiency. All operations consider
    only the required block of indices within each array during each iteration.

    """
    # number of taxa
    m = dm.shape[0]

    # number of nodes in the final tree
    n = 2 * m - 3

    # pre-allocate memory space
    tree, preodr, postodr, ad2, adk, lens, stack = _allocate_arrays(n, False)

    # initialize 3-taxon tree
    _init_tree(dm, tree, preodr, postodr, ad2, False)

    ### iteratively add taxa to the tree
    for k in range(3, m):
        # calculate average distances from new taxon to existing subtrees
        _avgdist_taxon(adk, k, dm, tree, preodr, postodr)

        # find a branch with minimum length change
        target = _ols_min_branch_d2(lens, ad2, adk, tree, preodr)

        # update average distances between distant-2 subtrees
        _avgdist_d2_insert(ad2, target, adk, tree, preodr)

        # insert new taxon into tree
        _insert_taxon(k, target, tree, preodr, postodr, False)

    # calculate branch lengths using an OLS framework
    _ols_lengths_d2(lens, ad2, tree, preodr)

    return tree, lens


def _bme(dm):
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
    # numbers of taxa and nodes in the tree
    m = dm.shape[0]
    n = 2 * m - 3

    # pre-allocate memory space and initialize 3-taxon tree (same as GME but creates a
    # full matrix)
    tree, preodr, postodr, adm, adk, lens, stack = _allocate_arrays(n, True)
    _init_tree(dm, tree, preodr, postodr, adm, True)

    # pre-calculate negative powers of 2
    powers = np.ldexp(1.0, -np.arange(m))

    for k in range(3, m):
        # calculate balanced average distances from new taxon to existing subtrees
        _bal_avgdist_taxon(adk, k, dm, tree, preodr, postodr)

        # find a branch with minimum length change
        target = _bal_min_branch(lens, adm, adk, tree, preodr)

        # update balanced average distance matrix between all subtrees
        _bal_avgdist_insert(adm, target, adk, tree, preodr, postodr, powers, stack)

        # insert new taxon into tree
        _insert_taxon(k, target, tree, preodr, postodr, True)

    # calculate branch lengths using a balanced framework
    _bal_lengths(lens, adm, tree, preodr)

    return tree, lens


# ------------------------------------------------
# Tree structure
# ------------------------------------------------


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
    assert not tree[n:].any()
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


# ------------------------------------------------
# Format conversion
# ------------------------------------------------


def _to_treenode(tree, taxa, lens=None, unroot=False):
    r"""Convert an array-based tree structure into a TreeNode object.

    Parameters
    ----------
    tree : (N, 2+) array_like of int
        Tree structure.
    taxa : list of str
        Taxon names in order.
    lens : (N,) array_like of float, optional
        Branch lengths.
    unroot : bool, optional
        If True, generate an unrooted tree with three children from the root node.

    Returns
    -------
    TreeNode
        Converted TreeNode object.

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


def _from_treenode(tree, taxa):
    r"""Convert a TreeNode object into an array-based tree structure.

    Parameters
    ----------
    tree : TreeNode
        Tree structure.
    taxa : list of str
        Taxon names in order.

    Returns
    -------
    (N, 8) ndarray of int
        Converted tree array.
    (N,) ndarray of int
        Nodes in preorder.
    (N,) ndarray of int
        Nodes in postorder.
    (N,) ndarray of float
        Branch lengths.

    Notes
    -----
    Nodes in the resulting tree array are ordered by a preorder traversal. This may not
    agree with the order they were added to the original tree, if that ever happened.
    Therefore, one shouldn't directly compare the original and generated tree arrays.
    However, the topologies represented by the two tree arrays should be identical.

    """
    # number of nodes in the full tree
    N = len(taxa) * 2 - 3

    # allocate arrays
    result = np.zeros((N, 8), dtype=int)
    preodr = np.zeros((N,), dtype=int)
    postodr = np.zeros((N,), dtype=int)
    lens = np.zeros((N,), dtype=float)

    # perform preorder traversal, mark nodes, fill parent, depth, taxon and length
    ordered = list(tree.preorder())
    for i, node in enumerate(ordered):
        node.i = i
        row = result[i]
        if node.parent is not None:
            row[2] = (pi := node.parent.i)
            row[5] = result[pi, 5] + 1
        if not node.children:
            row[1] = taxa.index(node.name)
        lens[i] = node.length or 0

    # number of nodes in the current tree
    n = len(ordered)

    # fill preorder and indices
    preodr[:n] = result[:n, 6] = np.arange(n)

    # perform postorder traversal, fill children, sibling, and size
    for i, node in enumerate(tree.postorder()):
        postodr[i] = node.i
        row = result[node.i]
        row[7] = i
        if not node.children:
            row[4] = 1
        else:
            left, right = node.children
            row[0] = (rrow := result[right.i])[3] = left.i
            row[1] = (lrow := result[left.i])[3] = right.i
            row[4] = lrow[4] + rrow[4]

    # clean up temporary attribute
    for node in ordered:
        delattr(node, "i")

    return result, preodr, postodr, lens


# ------------------------------------------------
# Tree traversal
# ------------------------------------------------


def _preorder_py(order, tree, stack, start=0):
    r"""Perform preorder traversal.

    This function and :func:`_postorder` use stacks to avoid recursion. The stack
    array is pre-allocated. The output (ordered nodes) is also written into a
    pre-allocated array.

    This function and :func:`_postorder` are not actually used in the greedy algorithms
    in this module, which incrementally grow the tree as well as the orders. The two
    functions are implemented for reference and test purpose.

    """
    stack[0] = start
    order_i = 0  # next index of order
    stack_i = 1  # next index of stack
    while stack_i:
        # pop a node from stack into order
        stack_i -= 1
        curr = stack[stack_i]
        order[order_i] = curr
        # index[curr] = order_i  # preorder index
        order_i += 1
        # append children to stack, right first such that left is processed first
        left, right = tree[curr, :2]
        if left:
            stack[stack_i] = right
            stack[stack_i + 1] = left
            stack_i += 2


def _postorder_py(order, tree, stack, start=0):
    """Perform postorder traversal.

    See also :func:`_preorder`.

    """
    stack[0] = start
    curr = tree[start, 0]
    order_i = 0
    stack_i = 1
    prev = 0  # last visited node
    while stack_i:
        if curr:
            stack[stack_i] = curr
            stack_i += 1
            curr = tree[curr, 0]  # go to left child
        else:
            last = stack[stack_i - 1]
            left, right = tree[last, :2]
            if left and prev != right:
                curr = right  # go to right child
            else:
                order[order_i] = last
                # index[last] = order_i  # postorder index
                order_i += 1
                prev = last
                stack_i -= 1


# ------------------------------------------------
# Initiation
# ------------------------------------------------


def _allocate_arrays(n, matrix=False):
    r"""Pre-allocate memory space for arrays.

    This function creates multiple arrays. Each array has a length of n on axis 1.
    Here, n is the number of nodes in the final tree. Thus, one doesn't need to remake
    arrays during iteration, but repeatedly uses the already allocated space. During an
    interation involving x nodes, the first x positions of each array are occupied,
    while the remaining positions are left as zero. There is no need to reset array
    values to zero after each iteration.

    Since they are created de novo, all arrays are C-continuous. This permits further
    optimization in the Cython code.

    `ads` stores the average distances between subtrees within the current tree. If
    `matrix` is True, an n by n square matrix will be allocated to cover all pairs of
    subtrees. This is required by BME and NNI. Such a matrix represents the primary
    memory bound in the entire algorithm.

    Otherwise, an n by 2 array will be allocated, only to include pairs of "distant-2"
    subtrees. This is adopted by GME to reduce computational complexity.

    """
    # tree structure
    tree = np.empty((n, 8), dtype=int)

    # nodes in pre- and postorder
    preodr = np.empty((n,), dtype=int)
    postodr = np.empty((n,), dtype=int)

    # average distances between subtrees
    x = n if matrix else 2
    ads = np.empty((n, x), dtype=float)

    # average distances from a taxon to each subtree
    adk = np.empty((n, 2), dtype=float)

    # branch lengths or length changes
    lens = np.empty((n,), dtype=float)

    # a stack for traversal operations
    stack = np.empty((n,), dtype=int)

    return tree, preodr, postodr, ads, adk, lens, stack


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


# ------------------------------------------------
# Tree manipulation
# ------------------------------------------------


def _insert_taxon(taxon, target, tree, preodr, postodr, use_depth=True):
    r"""Insert a taxon between a target node and its parent.

    For example, with the following local structure of the original tree:

          A
         / \
        B   C

    With target=B, this function inserts a taxon into the branch A-B. The structure
    becomes:

            A
           / \
        link  C
         / \
        B  taxon

    A special case is that the taxon is inserted into the root branch (node=0). The
    tree becomes:

            A
           / \
        link taxon
         / \
        B   C

    The inserted taxon always becomes the right child.

    """
    # determine tree dimensions
    # typically n = 2 * taxon - 3, but this function doesn't enforce this
    m = tree[0, 4]
    n = m * 2 - 1
    link = n
    tip = n + 1

    # Special case (root branch): taxon k becomes the sibling of all existing taxa
    # except for the root (taxon 0).
    if target == 0:
        # children
        left, right = tree[0, :2]
        tree[left, 2] = tree[right, 2] = link

        # root
        tree[0, :2] = [link, tip]
        tree[0, 4] = m + 1
        tree[0, 7] = n + 1

        # link
        tree[link] = [left, right, 0, tip, m, 1, 1, n - 1]

        # tip
        tree[tip] = [0, taxon, 0, link, 1, 1, n + 1, n]

        # entire tree depth + 1
        tree[1:n, 5] += 1

        # preorder
        tree[1:n, 6] += 1
        preodr[2 : n + 1] = preodr[1:n]
        preodr[1] = link
        preodr[n + 1] = tip

        # postorder
        postodr[n - 1 : n + 2] = [link, tip, 0]

    # Regular case (any other branch): The link becomes the parent of the target node,
    # and child of its original parent. Taxon k becomes the sibling
    else:
        left, right, parent, sibling, size, depth, pre_i, post_i = tree[target]
        side = (tree[parent, 0] != target).astype(int)
        tree[parent, side] = link
        tree[sibling, 3] = link
        tree[target, 2] = link
        tree[target, 3] = tip

        # identify nodes desceding from target
        clade_n = size * 2 - 1  # number of nodes in clade (including target)
        pre_i_after = pre_i + clade_n  # preorder index of node after clade
        clade = preodr[pre_i:pre_i_after]  # range of clade in preorder

        # link
        tree[link] = [
            target,
            tip,
            parent,
            sibling,
            size + 1,
            depth,
            pre_i,
            post_i + 2,
        ]

        # tip
        tree[tip] = [
            0,
            taxon,
            link,
            target,
            1,
            depth + 1,
            pre_i_after + 1,
            post_i + 1,
        ]

        # clade depth +1
        if use_depth:
            tree[clade, 5] += 1

        # preorder shift: nodes after clade +2, tip inserted after clade, nodes within
        # clade +1, link inserted before clade
        pre_after = preodr[pre_i_after:n]
        tree[pre_after, 6] += 2
        preodr[pre_i_after + 2 : n + 2] = pre_after
        preodr[pre_i_after + 1] = tip
        tree[clade, 6] += 1
        preodr[pre_i + 1 : pre_i_after + 1] = clade
        preodr[pre_i] = link

        # postorder shift: all nodes after clade +2, tip and link inserted after clade
        post_after = postodr[post_i + 1 : n]
        tree[post_after, 7] += 2
        postodr[post_i + 3 : n + 2] = post_after
        postodr[post_i + 2] = link
        postodr[post_i + 1] = tip

        # size +1 from link to root
        curr = link
        while curr:
            parent = tree[curr, 2]
            tree[parent, 4] += 1
            curr = parent


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


# ------------------------------------------------
# Calculation
# ------------------------------------------------


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
        for b in range(a + 1, n):
            # If a is an ancestor (proper superset) of b, flip a's taxon set as the
            # the upper subtree of a will be considered.
            if (a_taxa := taxas[a]) > (b_taxa := taxas[b]):
                a_taxa = full - a_taxa
            # Because this is postorder traversal, a cannot be a descendant (proper
            # subset) of b.
            adm[a, b] = adm[b, a] = dm[list(a_taxa)][:, list(b_taxa)].mean()


def _avgdist_matrix_py(adm, dm, tree, preodr, postodr):
    r"""Calculate a matrix of average distances between all pairs of subtrees.

    This function will update adm, a float array of (n, n) representing pairwise
    distances between all nodes (tips, internal nodes and the root) in the tree.

    Implemented according to Appendix 4 of Desper and Gascuel (2002). Basically, this
    algorithm traverses the tree and calculates the average distance between each pair
    of subtrees. Here, a "subtree" is identified by a node, but it can be one of the
    two scenarios:

    0. Lower subtree (L): the subtree descending from a node (including the node).
    1. Upper subtree (U): the subtree branching from the node upward. The root of this
       subtree is the node's parent. Its immediate children are the parent's parent
       and the node's sibling.

    Then it iteratively applies Eq. 2 to calculate the average distance between two
    subtrees based on the average distances from the children of one of the subtrees
    to the other.

        d(A, B) = (|B_1| * d(A, B_1) + |B_2| * d(A, B_2)) / |B|

    """
    # total numbers of taxa and nodes in the tree
    m = tree[0, 4] + 1
    n = 2 * m - 3

    # Calculate the average distance between each pair of subtrees defined by nodes
    # a and b (A4.1 of the paper).

    # Step 1: Calculate non-nested subtree-to-subtree distances (i.e., one is not an
    # ancestor of another. Therefore each subtree is the lower (descending) tree of
    # the node. (A4.1 (a))
    # The paper suggests outer and inner postorder traverals for all nodes. Since
    # distances are symmetric, one can save half of the calculations, as done by the
    # following code.

    # Loop over nodes in postorder.
    # Skip the root, which is always the last node in the postorder traversal.
    for a in postodr[: n - 1]:
        # check if a is a tip
        if (a_size := tree[a, 4]) == 1:
            a_taxon = tree[a, 1]
        else:
            a_taxon = 0  # a can never be root (taxon 0), therefore 0 means none
            a1, a2 = tree[a, :2]
            a1_size, a2_size = tree[a1, 4], tree[a2, 4]

        # Iterate over the ancestors of a, and find the other subtree on the right.
        # We can skip the left subtree because it must have been calculated already.
        # If the right subtree contains a's ancestry, just skip.
        curr = a
        while curr:
            parent = tree[curr, 2]
            if tree[parent, 0] != curr:
                curr = parent
                continue
            sibling = tree[curr, 3]

            # Loop over nodes within the right subtree in postorder.
            # This postorder doesn't need to be re-calculated. Because all nodes within
            # a clade are continuous in postorder, one can take a slice of the full
            # postorder that represent the descending nodes of the current node. The
            # size of the slice is taxon count * 2 - 2. *
            for b in postodr[
                (i := tree[sibling, 7]) - tree[sibling, 4] * 2 + 2 : i + 1
            ]:
                # If both a and b are tips, take the original taxon-to-taxon distance
                # (A4.1 (a) i).
                if (b_size := tree[b, 4]) == 1 and a_taxon:
                    dist = dm[a_taxon, tree[b, 1]]

                # If a is an internal node, and b is either (a tip or an internal node),
                # calculate the average distance based on the two child subtrees of a
                # (A4.1 (a) ii).
                elif not a_taxon:
                    dist = (a1_size * adm[a1, b] + a2_size * adm[a2, b]) / a_size

                # If a is a tip, and b is an internal node, calculate the average
                # distance based on the two child subtrees of b (A4.1 (a) iii).
                else:
                    b1, b2 = tree[b, :2]
                    dist = (
                        tree[b1, 4] * adm[a, b1] + tree[b2, 4] * adm[a, b2]
                    ) / b_size

                adm[a, b] = adm[b, a] = dist

            curr = parent

    # Step 2: Calculate subtree to root (taxon 0) distances (A4.1 (b)).

    # This is done through a postorder traversal.
    for a in postodr[: n - 1]:
        if (a_size := tree[a, 4]) == 1:
            a_taxon = tree[a, 1]
            dist = dm[0, a_taxon]
        else:
            a1, a2 = tree[a, :2]
            dist = (tree[a1, 4] * adm[a1, 0] + tree[a2, 4] * adm[a2, 0]) / a_size

        adm[a, 0] = adm[0, a] = dist

    # Step 3: Calculate nested subtree to subtree distances, in which the first node
    # (a) is a descendant of the second node (b), therefore the first subtree (A) is
    # the lower (descending) tree of node a, whereas the second subtree (B) is the
    # upper (ancestral) tree of node b (A4.1 (c)).

    # This is done through a preorder traversal.
    for a in preodr[1:n]:
        parent, sibling = tree[a, 2:4]

        # The size of (upper) subtree b is the complement of its descendants. Same for
        # the parent subtree.
        a_size = m - tree[a, 4]
        p_size = m - tree[parent, 4]
        s_size = tree[sibling, 4]

        # Iterate over all subtrees below b.
        # The paper says this traversal can be done in any manner. Here, we use the
        # postorder. See * above.
        for b in postodr[(i := tree[a, 7]) - tree[a, 4] * 2 + 2 : i]:
            dist = (s_size * adm[b, sibling] + p_size * adm[b, parent]) / a_size
            adm[a, b] = adm[b, a] = dist


def _bal_avgdist_matrix_py(adm, dm, tree, preodr, postodr):
    r"""Calculate a matrix of balanced average distances between all pairs of subtrees.

    This function resembles :func:`_avgdist_matrix`, but it weighs subtrees equally
    regardless of their sizes. Specifically, it replaces Eq. 2 with Eq. 6. of Desper
    and Gascuel (2002):

        d(A, B) = (d(A, B_1) + d(A, B_2)) / 2

    Same for all functions starting with `bal_`.

    """
    n = 2 * tree[0, 4] - 1

    # Step 1: Calculate non-nested subtree to subtree distances.
    for a in postodr[: n - 1]:
        if tree[a, 0] == 0:
            a_taxon = tree[a, 1]
        else:
            a_taxon = 0
            a1, a2 = tree[a, :2]
        curr = a
        while curr:
            parent = tree[curr, 2]
            if tree[parent, 0] != curr:
                curr = parent
                continue
            sibling = tree[curr, 3]
            for b in postodr[
                (i := tree[sibling, 7]) - tree[sibling, 4] * 2 + 2 : i + 1
            ]:
                if tree[b, 0] == 0 and a_taxon:
                    dist = dm[a_taxon, tree[b, 1]]
                elif not a_taxon:
                    dist = 0.5 * (adm[a1, b] + adm[a2, b])
                else:
                    b1, b2 = tree[b, :2]
                    dist = 0.5 * (adm[a, b1] + adm[a, b2])
                adm[a, b] = adm[b, a] = dist
            curr = parent

    # Step 2: Calculate subtree to root distances.
    for a in postodr[: n - 1]:
        if tree[a, 0] == 0:
            a_taxon = tree[a, 1]
            dist = dm[0, a_taxon]
        else:
            a1, a2 = tree[a, :2]
            dist = 0.5 * (adm[a1, 0] + adm[a2, 0])

        adm[a, 0] = adm[0, a] = dist

    # Step 3: Calculate nested subtree to subtree distances.
    for a in preodr[1:n]:
        parent, sibling = tree[a, 2:4]
        for b in postodr[(i := tree[a, 7]) - tree[a, 4] * 2 + 2 : i]:
            adm[a, b] = adm[b, a] = 0.5 * (adm[b, sibling] + adm[b, parent])


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


def _avgdist_taxon_py(adk, taxon, dm, tree, preodr, postodr):
    """Calculate average distances between a new taxon and existing subtrees.

    This function will update adk, a float array of (n, 2) in which columns 0 and 1
    represent the average distances from the taxon to the lower and upper subtrees of
    each existing node, respectively.

    Implemented according to Appendix 3 of Desper and Gascuel (2002). Basically, this
    algorithm calculates all lower subtree distances via a postorder traversal, then
    calculates all upper subtree distances via a preorder traversal.

    """
    # distance from taxon to all other taxa
    dk = dm[taxon]

    # total numbers of taxa and nodes in the tree
    m = tree[0, 4] + 1
    n = 2 * m - 3

    # Calculate the distance between taxon and the lower subtree of each node.
    for node in postodr[:n]:
        left, right = tree[node, :2]
        if left == 0:
            adk[node, 0] = dk[right]
        else:
            adk[node, 0] = (
                tree[left, 4] * adk[left, 0] + tree[right, 4] * adk[right, 0]
            ) / tree[node, 4]

    # Assign upper distance of root.
    adk[0, 1] = dk[0]

    # Calculate the distance between taxon k and the upper subtree of each node.
    for node in preodr[1:n]:
        adk[node, 1] = (
            (m - tree[parent := tree[node, 2], 4]) * adk[parent, 1]
            + tree[sibling := tree[node, 3], 4] * adk[sibling, 0]
        ) / (m - tree[node, 4])


def _bal_avgdist_taxon_py(adk, taxon, dm, tree, preodr, postodr):
    r"""Calculate balanced average distances between a new taxon and existing subtrees.

    This function resembles :func:`_avgdist_taxon` but uses the balanced framework.

    """
    n = tree[0, 4] * 2 - 1
    dk = dm[taxon]
    for node in postodr[:n]:
        left, right = tree[node, :2]
        if left == 0:
            adk[node, 0] = dk[right]
        else:
            adk[node, 0] = 0.5 * (adk[left, 0] + adk[right, 0])
    adk[0, 1] = dk[0]
    for node in preodr[1:n]:
        adk[node, 1] = 0.5 * (adk[tree[node, 2], 1] + adk[tree[node, 3], 0])


def _avgdist_d2_insert_py(ad2, target, adk, tree, preodr):
    r"""Update average distances between distant-2 subtrees after taxon insertion.

    This function will update ad2, a float array of (n, 2) representing pairwise
    distances between all distant-2 subtrees in the tree. Here, `distant-2 subtrees`
    refer to subtrees that are two branches away from each other. Specifically, there
    are two scenarios:

    - Column 0: Distance between the lower subtree of the current node and the upper
      subtree of its parent.
    - Column 1: Distance between the lower subtree of the current node and the lower
      subtree of its sibling.

    This function assumes that the taxon will be inserted into the branch connecting
    the target node and its parent. After insertion, the taxon will become the sibling
    of the target.

               parent
                /  \
             link  sibling
             /  \
        target  taxon

    This function should be executed *before* calling :func:`_insert_taxon`, which will
    mutate the tree.

    Implemented according to Eq. 8 of Desper and Gascuel (2002).

    """
    m = tree[0, 4] + 1
    n = 2 * m - 3
    link = n
    tip = n + 1

    ###### Special case: insert into the root branch. ######

    if target == 0:
        # k (lower) to root (parent, upper): pre-calculated
        ad2[tip, 0] = adk[0, 1]

        # k (lower) to link (sibling, lower): equals to k to root (lower).
        ad2[link, 1] = ad2[tip, 1] = adk[0, 0]

        # Link (lower) to root (parent, upper): de novo calculation according to the
        # equation in A4.1(b). It is basically the distance between the upper and lower
        # subtrees of the root itself.
        ad2[link, 0] = (
            tree[children := tree[0, :2], 4] * ad2[children, 0]
        ).sum() / tree[0, 4]

        # Calculate all node (lower) to parent (upper, containing k) distances. These
        # parents include the new link.
        ad2[1:n, 0] = (
            adk[1:n, 0] + ad2[1:n, 0] * (p_sizes := m - tree[tree[1:n, 2], 4])
        ) / (p_sizes + 1)

        return

    ###### Regular case: insert into any other branch. ######

    parent, sibling, size = tree[target, 2:5]

    # Temporarily copy the distances of target to link (will edit later).
    ad2[link, :2] = ad2[target, :2]

    # Distance between k (lower) and link (parent, upper) equals to that between k and
    # the upper subtree of target.
    ad2[tip, 0] = adk[target, 1]

    # Distance between target (lower) and link (parent, upper) needs to be calculated
    # using the equation in A4.1(c). Basically, it is the distance between the lower
    # and upper subtrees of the same target.
    ad2[target, 0] = (
        tree[sibling, 4] * ad2[target, 1] + (m - tree[parent, 4]) * ad2[target, 0]
    ) / (m - size)

    # Transfer the pre-calculated distance between target (lower) and k (sibling,
    # lower).
    ad2[target, 1] = ad2[tip, 1] = adk[target, 0]

    # Within the clade below target, calculate the distance between each node (lower)
    # and its parent (upper, containing k).
    clade = preodr[(i := tree[target, 6] + 1) : i + size * 2 - 2]
    ad2[clade, 0] = (
        adk[clade, 0] + ad2[clade, 0] * (p_sizes := m - tree[tree[clade, 2], 4])
    ) / (p_sizes + 1)

    # Iterate over the ancestors of target, starting from link and ending at root.
    curr = link
    while curr:
        # Calculate the distance between each pair of lower (containing k) and upper
        # ancestors.
        ad2[curr, 0] = (adk[parent, 1] + ad2[curr, 0] * size) / (size + 1)

        # Calculate the distance between each ancestor (lower, containing k) and its
        # sibling (lower).
        ad2[curr, 1] = ad2[sibling, 1] = (adk[sibling, 0] + ad2[curr, 1] * size) / (
            size + 1
        )

        # Within the clade below each sibling, calculate the distance between each node
        # (lower) and its parent (upper, containing k).
        clade = preodr[(i := tree[sibling, 6] + 1) : i + tree[sibling, 4] * 2 - 2]
        ad2[clade, 0] = (
            adk[clade, 0] + ad2[clade, 0] * (p_sizes := m - tree[tree[clade, 2], 4])
        ) / (p_sizes + 1)

        curr = parent
        parent, sibling, size = tree[curr, 2:5]


def _bal_avgdist_insert_py(adm, target, adk, tree, preodr, postodr, powers, stack):
    r"""Update balanced average distance matrix after taxon insertion.

    This function resembles :func:`_avgdist_d2_insert` but it 1) uses the balanced
    framework and 2) updates the entire matrix.

    Two additional parameters are provided: `powers` is a pre-calculated array of
    2^(-l) powers (l is the depth difference between two nodes). `stack` is an
    integer array to store ancestral nodes of target.

    """
    m = tree[0, 4] + 1
    n = 2 * m - 3
    link = n
    tip = n + 1

    ###### Special case: insert into the root branch. ######

    if target == 0:
        # Transfer distance between k and root (upper).
        adm[0, tip] = adm[tip, 0] = adk[0, 1]

        # Transfer distances from k to any non-root nodes (lower).
        adm[1:n, tip] = adm[tip, 1:n] = adk[1:n, 0]

        # k to link: equals to k to root (lower).
        adm[link, tip] = adm[tip, link] = adk[0, 0]

        # Root to link: de novo calculation according to the equation in A4.1(b). It is
        # basically the distance between the upper and lower subtrees of the root.
        a1, a2 = tree[0, :2]
        adm[0, link] = adm[link, 0] = 0.5 * (adm[a1, 0] + adm[a2, 0])

        # Non-root to link: use Eq. 8 to calculate the distance between the upper tree
        # with two taxa (0 and k) and any non-root node.
        adm[1:n, link] = adm[link, 1:n] = 0.5 * (adk[1:n, 0] + adm[1:n, 0])

        # Calculate all ancestor (a, upper, containing k) to descendant (b, lower)
        # distances.
        for a in range(1, n):
            if tree[a, 0] > 0:
                descs = postodr[(i := tree[a, 7]) - tree[a, 4] * 2 + 2 : i]
                adm[a, descs] = adm[descs, a] = adm[a, descs] + powers[
                    tree[a, 5] + 1
                ] * (adk[descs, 0] - adm[0, descs])

        return

    ###### Regular case: insert into any other branch. ######

    parent, sibling, size, depth = tree[target, 2:6]

    ### Step 1: Distances that don't need iteration. ###

    # Distance between k (lower) and link (upper) equals to that between k and the
    # upper subtree of target.
    adm[tip, link] = adm[link, tip] = adk[target, 1]

    # Distance between target (lower) and link (upper) needs to be calculated using the
    # equation in A4.1(c). Basically, it is the distance between the lower and upper
    # subtrees of the same target.
    adm[target, link] = adm[link, target] = 0.5 * (
        adm[target, sibling] + adm[target, parent]
    )

    ### Step 2: Distances within the clade below target. ###

    # Locate the clade below target (including target)
    start = (i := tree[target, 7]) - size * 2 + 2

    # Transfer pre-calculated distance between k (lower) and any node within the clade
    # (including target, lower).
    clade = postodr[start : i + 1]
    adm[clade, n + 1] = adm[n + 1, clade] = adk[clade, 0]

    # Distance from any descendant (lower) to link (upper) equals to that to target.
    descs = postodr[start:i]
    adm[descs, n] = adm[n, descs] = adm[descs, target]

    # Within the clade, find all ancestor (a) - descendant (b) pairs, and calculate the
    # distance between the upper subtree of a (containing k) and the lower subtree of b.
    for a in clade:
        if tree[a, 0]:
            descs = postodr[(i := tree[a, 7]) - tree[a, 4] * 2 + 2 : i]
            adm[a, descs] = adm[descs, a] = adm[a, descs] + powers[
                tree[a, 5] - depth + 1
            ] * (adk[descs, 0] - adm[target, descs])

    ### Step 3: Distances among nodes outside the clade. ###

    # Iterate over ancestors of target in ascending order.
    i = 0
    curr = target
    while curr:
        stack[i] = anc = tree[curr, 2]

        # Transfer the pre-calculated distance between k and the ancestor (upper).
        adm[anc, n + 1] = adm[n + 1, anc] = adk[anc, 1]

        # Calculate the distance between link (lower, containing k) and the ancestor
        # (upper).
        adm[anc, n] = adm[n, anc] = 0.5 * (adk[anc, 1] + adm[anc, target])

        # Calculate the distance between each previous ancestor (lower, containing k)
        # and the current ancestor (upper).
        prevs = stack[:i]
        adm[anc, prevs] = adm[prevs, anc] = adm[anc, prevs] + powers[
            depth - tree[prevs, 5] + 1
        ] * (adk[anc, 1] - adm[target, anc])

        # Identify the sibling clade descending from the ancestor.
        cousin = tree[curr, 3]
        clade = postodr[(j := tree[cousin, 7]) - tree[cousin, 4] * 2 + 2 : j + 1]

        # Transfer the pre-calculated distances between k and each descendant
        # (lower).
        adm[clade, tip] = adm[tip, clade] = adk[clade, 0]

        # Calculate the distance between link (lower, containing k) and each
        # descendant (lower).
        adm[clade, link] = adm[link, clade] = 0.5 * (adk[clade, 0] + adm[clade, target])

        # Calculate the distance between each previous ancestor (lower, containing k)
        # and each descendant (lower).
        for a in prevs:
            adm[a, clade] = adm[clade, a] = adm[a, clade] + powers[
                depth - tree[a, 5] + 1
            ] * (adk[clade, 0] - adm[clade, target])

        # Iterate over descendants of each member of the clade, and calculate the
        # distance between the former (upper, containing k) and the latter (lower).
        for a in clade:
            if tree[a, 0]:
                descs = postodr[(j := tree[a, 7]) - tree[a, 4] * 2 + 2 : j]
                adm[a, descs] = adm[descs, a] = adm[a, descs] + powers[
                    depth + tree[a, 5] - 2 * tree[anc, 5]
                ] * (adk[descs, 0] - adm[descs, target])

        curr = anc
        i += 1


# ------------------------------------------------
# Branch lengths
# ------------------------------------------------


def _ols_lengths_py(lens, adm, tree, preodr):
    r"""Calculate branch lengths of a tree based on the OLS framework.

    Using an average distance matrix between all pairs of subtrees.

    Implemented according to Eqs. 3 & 4 of Desper and Gascuel (2002).

    """
    # total numbers of taxa and nodes in the tree
    m = tree[0, 4] + 1
    n = 2 * m - 3

    for node in preodr[1:n]:
        left, right, parent, sibling = tree[node, :4]

        # External (terminal) branch: based on the triplet (iAB) of self (i), parent
        # (A), and sibling (B) (Eq. 4).
        if not left:
            lens[node] = 0.5 * (
                adm[parent, node] + adm[sibling, node] - adm[parent, sibling]
            )

        # Internal branch: based on the quartet (AB|CD) of parent (A, upper), sibling
        # (B), and children (C and D) (Eq. 3).
        else:
            l_size = tree[left, 4]
            r_size = tree[right, 4]
            p_size = m - tree[parent, 4]
            s_size = tree[sibling, 4]
            lambda_ = (p_size * r_size + s_size * l_size) / (
                (p_size + s_size) * (l_size + r_size)
            )
            lens[node] = 0.5 * (
                lambda_ * (adm[parent, left] + adm[sibling, right])
                + (1 - lambda_) * (adm[parent, right] + adm[sibling, left])
                - (adm[parent, sibling] + adm[left, right])
            )

    # root branch
    left, right = tree[0, :2]
    lens[0] = 0.5 * (adm[left, 0] + adm[right, 0] - adm[left, right])


def _ols_lengths_d2_py(lens, ad2, tree, preodr):
    r"""Calculate branch lengths of a tree based on an OLS framework.

    Using only average distances between pairs of distant-2 subtrees.

    This function produces the same result as `_ols_lengths`. The latter
    relies on Eq. 3 of Desper and Gascuel (2002), which involves distances between
    distant-3 subtrees. Therefore, we need to modify the equation such that it only
    takes distances between distant-2 subtrees as input.

    Specifically, with the following local structure:

                |
              parent
              /   \
           node  sibling
           /  \
        left  right

    We will need to calculate the distances between left / right (lower) and parent
    (upper) / sibling (lower). This can be achieved by (according to Fig. 2a):

        d(left(L), parent(U)) = d(left(L), node(U)) + d(node(L), parent(U))
            - d(node(L), node(U))

    Do the same for d(right(L), parent(U)), d(left(L), sibling(L)), and d(right(L),
    sibling(L)). Plug results into Eq. 3, we will get:

        l(node) = 0.5 * (
            d(left(L), parent(U)) + d(right(L), parent(U)) + d(node(L), parent(U))
            + d(node(L), sibling(U)) - d(sibling(L), parent(U)) - d(left(L), right(L))
        ) - d(node(L), node(U))

    Note that this equation is free of the lambda factor.

    Here, d(node(L), node(U)) is the distance between the lower and upper subtrees of
    the same node. This can be calculated using the equation in A4.1(c):

        d(node(L), node(U)) = (|sibling(L)| * d(node(L), sibling(L)) + |parent(U)|
            * d(node(L), parent(U))) / |node(U)|

    Therefore, we will get l(node).

    """
    # total numbers of taxa and nodes in the tree
    m = tree[0, 4] + 1
    n = 2 * m - 3

    for node in preodr[1:n]:
        left, right, parent, sibling = tree[node, :4]

        # External (terminal) branch: based on Eq. 4 of the paper.
        if not left:
            lens[node] = 0.5 * (ad2[node, 0] + ad2[node, 1] - ad2[sibling, 0])

        # Internal branch: based on the equation discussed above.
        else:
            lens[node] = 0.5 * (
                ad2[left, 0]
                + ad2[right, 0]
                + ad2[node, 0]
                + ad2[node, 1]
                - ad2[sibling, 0]
                - ad2[left, 1]
            ) - (
                tree[sibling, 4] * ad2[node, 1] + (m - tree[parent, 4]) * ad2[node, 0]
            ) / (m - tree[node, 4])

    # root branch
    left, right = tree[0, :2]
    lens[0] = 0.5 * (ad2[left, 0] + ad2[right, 0] - ad2[left, 1])


def _bal_lengths_py(lens, adm, tree, preodr):
    r"""Calculate branch lengths of a tree based on the balanced framework.

    Using a balanced average distance matrix between all pairs of subtrees.

    This function resembles :func:`_ols_lengths` but it uses the balanced
    framework.

    Implemented according to Eqs. 3 & 4 and the description on top of pg. 691 of Desper
    and Gascuel (2002).

    """
    m = tree[0, 4] + 1
    n = 2 * m - 3
    for node in preodr[1:n]:
        left, right, parent, sibling = tree[node, :4]
        if not left:
            lens[node] = 0.5 * (
                adm[parent, node] + adm[sibling, node] - adm[parent, sibling]
            )
        else:
            lens[node] = 0.25 * (
                adm[parent, left]
                + adm[sibling, right]
                + adm[parent, right]
                + adm[sibling, left]
            ) - 0.5 * (adm[parent, sibling] + adm[left, right])
    left, right = tree[0, :2]
    lens[0] = 0.5 * (adm[left, 0] + adm[right, 0] - adm[left, right])


# ------------------------------------------------
# Branch search
# ------------------------------------------------


def _ols_min_branch_d2_py(lens, ad2, adk, tree, preodr):
    """Find the branch with the minimum length change after inserting a new taxon.

    It returns the node at the lower end of the branch.

    Implemented according to Eq. 7 of Desper and Gascuel (2002).

    IMPORTANT NOTE: Should there are ties (which are common), this function returns
    the first minimum branch seen during preorder traversal. This behavior resembles
    FastME. Alternatively, one can return the first minimum branch by the order of
    addition using `return lens[:n].argmin()`, which could produce different results
    at the presence of ties. It must be noted that all candidates in a tie are equally
    optimal in the current iteration of the greedy algorithm.

    """
    min_node = 0
    lens[min_node] = min_len = 0

    m = tree[0, 4] + 1
    n = 2 * m - 3
    for node in preodr[1:n]:
        parent, sibling, size = tree[node, 2:5]
        p_size = m - tree[parent, 4]
        s_size = tree[sibling, 4]

        numerator = s_size + size * p_size
        lambda_0 = numerator / (s_size + size) / (p_size + 1)
        lambda_1 = numerator / (s_size + p_size) / (size + 1)

        lens[node] = L = lens[parent] + 0.5 * (
            (lambda_0 - lambda_1) * (adk[sibling, 0] + ad2[node, 0])
            + (lambda_1 - 1) * (ad2[sibling, 1] + adk[parent, 1])
            + (1 - lambda_0) * (ad2[sibling, 0] + adk[node, 0])
        )

        if L < min_len:
            min_len, min_node = L, node

    return min_node


def _bal_min_branch_py(lens, adm, adk, tree, preodr):
    """Find the branch with the minimum length change after inserting a new taxon.

    This function resembles :func:`_ols_min_branch_d2` but it 1) uses the
    balanced framework and 2) calculates based on the entire matrix. See also
    the important note of the latter.

    Implemented according to Eq. 10 of Desper and Gascuel (2002).

    """
    min_node = 0
    lens[min_node] = min_len = 0
    for node in preodr[1 : 2 * tree[0, 4] - 1]:
        parent, sibling, size = tree[node, 2:5]
        lens[node] = L = lens[parent] + 0.25 * (
            adm[sibling, parent] + adk[node, 0] - adm[sibling, node] - adk[parent, 1]
        )
        if L < min_len:
            min_len, min_node = L, node
    return min_node
