# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io

import numpy as np
import heapq as hq

from skbio.stats.distance import DistanceMatrix
from skbio.tree import TreeNode


def nj(dm, disallow_negative_branch_length=True, result_constructor=None):
    r"""Apply neighbor joining for phylogenetic reconstruction.

    Parameters
    ----------
    dm : skbio.DistanceMatrix
        Input distance matrix containing distances between taxa.
    disallow_negative_branch_length : bool, optional
        Neighbor joining can result in negative branch lengths, which don't
        make sense in an evolutionary context. If `True`, negative branch
        lengths will be returned as zero, a common strategy for handling this
        issue that was proposed by the original developers of the algorithm.
    result_constructor : function, optional
        Function to apply to construct the result object. This must take a
        newick-formatted string as input. The result of applying this function
        to a newick-formatted string will be returned from this function. This
        defaults to ``lambda x: TreeNode.read(StringIO(x), format='newick')``.

    Returns
    -------
    TreeNode
        By default, the result object is a `TreeNode`, though this can be
        overridden by passing `result_constructor`.

    See Also
    --------
    TreeNode.root_at_midpoint

    Notes
    -----
    Neighbor joining was initially described in Saitou and Nei (1987) [1]_. The
    example presented here is derived from the Wikipedia page on neighbor
    joining [2]_. Gascuel and Steel (2006) provide a detailed overview of
    Neighbor joining in terms of its biological relevance and limitations [3]_.

    Neighbor joining, by definition, creates unrooted trees. One strategy for
    rooting the resulting trees is midpoint rooting, which is accessible as
    ``TreeNode.root_at_midpoint``.

    References
    ----------
    .. [1] Saitou N, and Nei M. (1987) "The neighbor-joining method: a new
       method for reconstructing phylogenetic trees." Molecular Biology and
       Evolution. PMID: 3447015.
    .. [2] http://en.wikipedia.org/wiki/Neighbour_joining
    .. [3] Gascuel O, and Steel M. (2006) "Neighbor-Joining Revealed" Molecular
       Biology and Evolution, Volume 23, Issue 11, November 2006,
       Pages 1997–2000, https://doi.org/10.1093/molbev/msl072

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
              /-d
             |
             |          /-c
             |---------|
    ---------|         |          /-b
             |          \--------|
             |                    \-a
             |
              \-e

    Again, construct the neighbor joining tree, but instead return the newick
    string representing the tree, rather than the TreeNode object. (Note that
    in this example the string output is truncated when printed to facilitate
    rendering.)

    >>> newick_str = nj(dm, result_constructor=str)
    >>> print(newick_str[:55], "...")
    (d:2.000000, (c:4.000000, (b:3.000000, a:2.000000):3.00 ...

    """
    if dm.shape[0] < 3:
        raise ValueError(
            "Distance matrix must be at least 3x3 to "
            "generate a neighbor joining tree."
        )

    if result_constructor is None:

        def result_constructor(x):
            return TreeNode.read(io.StringIO(x), format="newick")

    # initialize variables
    node_definition = None

    # while there are still more than three distances in the distance matrix,
    # join neighboring nodes.
    while dm.shape[0] > 3:
        # compute the Q matrix
        q = _compute_q(dm)

        # identify the pair of nodes that have the lowest Q value. if multiple
        # pairs have equally low Q values, the first pair identified (closest
        # to the top-left of the matrix) will be chosen. these will be joined
        # in the current node.
        idx1, idx2 = _lowest_index(q)
        pair_member_1 = dm.ids[idx1]
        pair_member_2 = dm.ids[idx2]
        # determine the distance of each node to the new node connecting them.
        pair_member_1_len, pair_member_2_len = _pair_members_to_new_node(
            dm, idx1, idx2, disallow_negative_branch_length
        )
        # define the new node in newick style
        node_definition = "(%s:%f, %s:%f)" % (
            pair_member_1,
            pair_member_1_len,
            pair_member_2,
            pair_member_2_len,
        )
        # compute the new distance matrix, which will contain distances of all
        # other nodes to this new node
        dm = _compute_collapsed_dm(
            dm,
            pair_member_1,
            pair_member_2,
            disallow_negative_branch_length=disallow_negative_branch_length,
            new_node_id=node_definition,
        )

    # When there are three distances left in the distance matrix, we have a
    # fully defined tree. The last node is internal, and its distances are
    # defined by these last three values.
    # First determine the distance between the last two nodes to be joined in
    # a pair...
    pair_member_1 = dm.ids[1]
    pair_member_2 = dm.ids[2]
    pair_member_1_len, pair_member_2_len = _pair_members_to_new_node(
        dm, pair_member_1, pair_member_2, disallow_negative_branch_length
    )
    # ...then determine their distance to the other remaining node, but first
    # handle the trivial case where the input dm was only 3 x 3
    node_definition = node_definition or dm.ids[0]
    internal_len = 0.5 * (
        dm[pair_member_1, node_definition]
        + dm[pair_member_2, node_definition]
        - dm[pair_member_1, pair_member_2]
    )
    if disallow_negative_branch_length and internal_len < 0:
        internal_len = 0

    # ...and finally create the newick string describing the whole tree.
    newick = "(%s:%f, %s:%f, %s:%f);" % (
        pair_member_1,
        pair_member_1_len,
        node_definition,
        internal_len,
        pair_member_2,
        pair_member_2_len,
    )

    # package the result as requested by the user and return it.
    return result_constructor(newick)


def _compute_q(dm):
    """Compute Q matrix, used to identify the next pair of nodes to join."""
    q = np.zeros(dm.shape)
    n = dm.shape[0]
    big_sum = np.array([dm.data.sum(1)] * dm.shape[0])
    big_sum_diffs = big_sum + big_sum.T
    q = (n - 2) * dm.data - big_sum_diffs
    np.fill_diagonal(q, 0)
    return DistanceMatrix(q, dm.ids)


def _compute_collapsed_dm(dm, i, j, disallow_negative_branch_length, new_node_id):
    """Return the distance matrix resulting from joining ids i and j in a node.

    If the input distance matrix has shape ``(n, n)``, the result will have
    shape ``(n-1, n-1)`` as the ids `i` and `j` are collapsed to a single new
    ids.

    """
    in_n = dm.shape[0]
    out_n = in_n - 1
    out_ids = [new_node_id]
    out_ids.extend([e for e in dm.ids if e not in (i, j)])
    result = np.zeros((out_n, out_n))
    # pre-populate the result array with known distances
    ij_indexes = [dm.index(i), dm.index(j)]
    result[1:, 1:] = np.delete(
        np.delete(dm.data, ij_indexes, axis=0), ij_indexes, axis=1
    )
    # calculate the new distances from the current DistanceMatrix
    k_to_u = 0.5 * (dm[i] + dm[j] - dm[i, j])
    # set negative branches to 0 if specified
    if disallow_negative_branch_length:
        k_to_u[k_to_u < 0] = 0
    # drop nodes being joined
    k_to_u = np.delete(k_to_u, ij_indexes)
    # assign the distances to the result array
    result[0] = result[:, 0] = np.concatenate([[0], k_to_u])
    return DistanceMatrix(result, out_ids)


def _lowest_index(dm):
    """Return the index of the lowest value in the input distance matrix.

    If there are ties for the lowest value, the index of top-left most
    occurrence of that value will be returned.

    This should be ultimately be replaced with a new DistanceMatrix object
    method (#228).

    """
    # get the positions of the lowest value
    results = np.vstack(np.where(dm.data == np.amin(dm.condensed_form()))).T
    # select results in the bottom-left of the array
    results = results[results[:, 0] > results[:, 1]]
    # calculate the distances of the results to [0, 0]
    res_distances = np.sqrt(results[:, 0] ** 2 + results[:, 1] ** 2)
    # detect distance ties & return the point which would have
    # been produced by the original function
    if np.count_nonzero(res_distances == np.amin(res_distances)) > 1:
        eqdistres = results[res_distances == np.amin(res_distances)]
        res_coords = eqdistres[np.argmin([r[0] for r in eqdistres])]
    else:
        res_coords = results[np.argmin(res_distances)]

    return tuple([res_coords[0], res_coords[1]])


def _pair_members_to_new_node(dm, i, j, disallow_negative_branch_length):
    """Return the distance between a new node and descendants of that new node.

    Parameters
    ----------
    dm : skbio.DistanceMatrix
        The input distance matrix.
    i, j : str
        Identifiers of entries in the distance matrix to be collapsed (i.e.,
        the descendants of the new node, which is internally represented as
        `u`).
    disallow_negative_branch_length : bool
        Neighbor joining can result in negative branch lengths, which don't
        make sense in an evolutionary context. If `True`, negative branch
        lengths will be returned as zero, a common strategy for handling this
        issue that was proposed by the original developers of the algorithm.

    """
    n = dm.shape[0]
    i_to_j = dm[i, j]
    i_to_u = (0.5 * i_to_j) + ((dm[i].sum() - dm[j].sum()) / (2 * (n - 2)))

    if disallow_negative_branch_length and i_to_u < 0:
        i_to_u = 0

    j_to_u = i_to_j - i_to_u

    if disallow_negative_branch_length and j_to_u < 0:
        j_to_u = 0

    return i_to_u, j_to_u


def nni(tree, dm, inplace=True):
    r"""Perform nearest neighbor interchanges on a phylogenetic tree.

    Parameters
    ----------
    tree : skbio.TreeNode
        Input phylogenetic tree to be rearranged.
    dm : skbio.DistanceMatrix
        Input distance matrix containing distances between taxa.
    inplace : bool, optional
        Whether manipulate the tree in place (``True``, default) or return a
        copy of the tree (``False``).

    Returns
    -------
    TreeNode
        Rearranged phylogenetic tree (if ``inplace`` is ``True``).

    Notes
    -----
    NNI algorithm for minimum evolution problem on phylogenetic trees. It rearranges
    an initial tree topology by performing subtree exchanges such that the distance
    is minimized. This implementation is based on the FastNNI algorithm [1]_.

    References
    ----------
    .. [1] Desper R, Gascuel O. Fast and accurate phylogeny reconstruction
       algorithms based on the minimum-evolution principle. J Comput Biol.
       2002;9(5):687-705. doi: 10.1089/106652702761034136. PMID: 12487758.

    """
    # Needs errors for nonbinary trees and trees rooted on internal node
    # Initialize and populate the average distance matrix
    if not inplace:
        tree = tree.copy()
    adm = _average_distance_matrix(tree, dm)
    while True:
        # create heap of possible swaps and then swapping subtrees
        # until no more swaps are possible.
        adm = _average_distance_matrix(tree, dm)
        heap = _swap_heap(tree, adm)
        if not heap:
            break
        swap = hq.heappop(heap)
        _perform_swap(swap[1][0], swap[1][1])
    # edge values are added using an OLS framework.
    _edge_estimation(tree, dm)
    if not inplace:
        return tree


def _perform_swap(node1, node2):
    """Returns a tree after swapping two subtrees."""
    parent1, parent2 = node1.parent, node2.parent
    parent1.append(node2)
    parent2.append(node1)


def _average_distance(node1, node2, dm):
    """Returns the average distance between the leaves of two subtrees.

    Distances between nodes are calculated using a distance matrix.
    """
    nodelist1 = _tip_or_root(node1)
    nodelist2 = _tip_or_root(node2)
    df = dm.between(nodelist1, nodelist2)
    return df["value"].mean()


def _tip_or_root(node):
    """Get name(s) of a node if it's a tip or root, otherwise its descending tips."""
    if node.is_tip() or node.is_root():
        return [node.name]
    else:
        return [x.name for x in node.tips()]


def _average_distance_upper(node1, node2, dm):
    """Returns the average distance between the leaves of two subtrees.

    Used for subtrees that contain ancestors corresponding to the root of the
    tree as descendants within the subtree.
    """
    nodelist1 = _tip_or_root(node1)
    if node2.is_root():
        nodelist2 = []
    else:
        root2 = node2.root()
        nodelist2 = [root2.name]
        nodelist2.extend(root2.subset() - node2.subset())
    df = dm.between(nodelist1, nodelist2)
    return df["value"].mean()


def _subtree_count(subtree):
    """Returns the number of leaves in a subtree.

    Assumes the root as a leaf node.
    """
    if subtree.is_tip() or subtree.is_root():
        return 1
    else:
        return subtree.count(tips=True)


def _swap_length(a, b, c, d, i, j, k, m, adm):
    """Returns the change in overall tree length after a given swap."""
    lambda1 = (a * d + b * c) / ((a + b) * (c + d))
    lambda2 = (a * d + b * c) / ((a + c) * (b + d))
    return 0.5 * (
        (lambda1 - 1) * (adm[i][k] + adm[j][m])
        - (lambda2 - 1) * (adm[i][j] + adm[k][m])
        - (lambda1 - lambda2) * (adm[i][m] + adm[j][k])
    )


def _swap_heap(tree, adm):
    """ """
    heap = []
    ordered = list(tree.postorder(include_self=False))
    root = tree.root()
    n_taxa = root.count(tips=True) + 1
    for node in ordered:
        if node.is_tip():
            continue
        parent = node
        a = parent.parent
        for index, node in enumerate(ordered):
            if node == a:
                i1 = index
        for child in parent.children:
            if child.is_tip():
                continue
            else:
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
                if c.is_tip():
                    c_ = 1
                else:
                    c_ = c.count(tips=True)
                if d.is_tip():
                    d_ = 1
                else:
                    d_ = d.count(tips=True)
                if b.is_tip():
                    b_ = 1
                else:
                    b_ = b.count(tips=True)
                a_ = n_taxa - b_ - c_ - d_
                swap_1 = _swap_length(a_, b_, c_, d_, i1, i2, i3, i4, adm)
                swap_2 = _swap_length(a_, b_, d_, c_, i1, i2, i4, i3, adm)
                if swap_1 > swap_2 and swap_1 > 0:
                    swap = -1 * swap_1
                    hq.heappush(heap, (swap, (b, c)))
                elif swap_2 > swap_1 and swap_2 > 0:
                    swap = -1 * swap_2
                    hq.heappush(heap, (swap, (b, d)))
    return heap


def _average_subtree_distance(a, b, a1, a2, dm):
    """ """
    return (
        _subtree_count(a1) * _average_distance(a1, b, dm)
        + _subtree_count(a2) * _average_distance(a2, b, dm)
    ) / _subtree_count(a)


def _average_distance_matrix(tree, dm):
    """ """
    ordered = list(tree.postorder(include_self=False))
    n = len(ordered)
    r = tree.root()
    taxa_size = r.count(tips=True) + 1
    adm = np.empty((n, n))
    for i, a in enumerate(ordered):  # consider refactoring with tip-to-tip distance
        if a in tree.children:  # skipping over unique descendant
            continue
        if a.is_tip():  # part (b)
            adm[n - 1, i] = adm[i, n - 1] = dm[a.name, r.name]
        else:
            a1, a2 = a.children
            adm[n - 1, i] = _average_subtree_distance(a, r, a1, a2, dm)
            adm[i, n - 1] = adm[n - 1, i]
        for j in range(i + 1, n - 1):  # part (a)
            b = ordered[j]
            if b in a.ancestors():  # skipping over ancestors
                continue
            if a.is_tip() and b.is_tip():
                adm[i, j] = adm[j, i] = dm[a.name, b.name]
            elif b.is_tip():
                a1, a2 = a.children
                adm[i, j] = adm[j, i] = _average_subtree_distance(a, b, a1, a2, dm)
            else:
                b1, b2 = b.children
                adm[i, j] = adm[j, i] = _average_subtree_distance(b, a, b1, b2, dm)
    for j, b in enumerate(ordered):  # calculating for ancestors
        if b in tree.children:  # skipping over unique descendant
            continue
        s_ = b.siblings()
        for sibling in s_:
            s = sibling
        p = b.parent
        for i, a in enumerate(ordered):
            if b in a.ancestors():
                adm[i, j] = (
                    _subtree_count(s) * _average_distance(a, s, dm)
                    + (taxa_size - _subtree_count(p))
                    * _average_distance_upper(a, p, dm)
                ) / (taxa_size - _subtree_count(b))
                adm[j, i] = adm[i, j]
    return adm


def _edge_estimation(tree, dm):
    """ """
    adm = _average_distance_matrix(tree, dm)
    ordered = list(tree.postorder(include_self=False))
    root = tree.root()
    taxa_size = root.count(tips=True) + 1
    for treenode in ordered:
        if treenode.is_root():
            continue
        elif treenode.parent.is_root():
            for index, node in enumerate(ordered):
                if node == treenode:
                    i1 = index
            a, b = treenode.children
            for index, node in enumerate(ordered):
                if node == a:
                    i2 = index
                elif node == b:
                    i3 = index
            length = 0.5 * (adm[i2][i1] + adm[i3][i1] - adm[i2][i3])
            treenode.length = length
        elif treenode.is_tip():
            parentnode = treenode.parent
            a = parentnode.parent
            if a.is_root():
                for child in a.children:
                    a = child
            for siblingnode in treenode.siblings():
                b = siblingnode
                for index, node in enumerate(ordered):
                    if node == treenode:
                        i1 = index
                    if node == a:
                        i2 = index
                    if node == b:
                        i3 = index
            length = 0.5 * (adm[i2][i1] + adm[i3][i1] - adm[i2][i3])
            treenode.length = length
        else:
            parentnode = treenode.parent
            a = parentnode.parent
            if a.is_root():
                for child in a.children:
                    a = child
            for index, node in enumerate(ordered):
                if node == a:
                    i1 = index
            c, d = treenode.children
            for sibling in treenode.siblings():
                b = sibling
            for index, node in enumerate(ordered):
                if node == b:
                    i2 = index
                elif node == c:
                    i3 = index
                elif node == d:
                    i4 = index
            if c.is_tip():
                c_ = 1
            else:
                c_ = c.count(tips=True)
            if d.is_tip():
                d_ = 1
            else:
                d_ = d.count(tips=True)
            if b.is_tip():
                b_ = 1
            else:
                b_ = b.count(tips=True)
            a_ = taxa_size - b_ - c_ - d_
            lambda1 = (a_ * d_ + b_ * c_) / ((a_ + b_) * (c_ + d_))
            length = 0.5 * (
                (lambda1 * (adm[i1][i3] + adm[i2][i4]))
                + ((1 - lambda1) * (adm[i1][i4] + adm[i2][i3]))
                - (adm[i1][i2] + adm[i3][i4])
            )
            treenode.length = length
