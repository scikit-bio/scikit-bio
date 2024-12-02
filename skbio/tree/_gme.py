# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from operator import itemgetter
from skbio.tree._util import (
    _ordered,
    _parent,
    _sibling,
    _ancestors,
    _subtree,
    _move_subtree,
    _array_to_tree,
)
from ._cutils import num_dist_cy


def _gme(dm):
    """Perform greedy minimum evolution (GME).

    Parameters
    ----------
    dm : (N, N) ndarray
        Input distance matrix.

    Returns
    -------
    (2, N - 1) ndarray
        Output tree as an array.

    """
    # count the number of taxa
    n = dm.shape[0]

    # create initial nodes
    tree = np.array([[1, 2] + [0] * (n - 3), [1, 2] + [0] * (n - 3)])

    # loop over rest of nodes for k = 3 to n
    for k in range(3, n):
        leaves = tree[:, tree[0] != 0]
        taxa = len(leaves[0]) + 1
        ordered = _ordered(tree)

        # Calculate average distance matrix between subtrees of T_(k-1).
        # Note that this should be replaced with a step to update an initial
        # computation of the average distance matrix rather than to recompute
        # the matrix after every iteration as implemented here.
        adm = _average_distance_matrix(dm, ordered, leaves, 0)

        # create average distance lists for subtrees of T_(k-1)
        lowerlist = _lower_subtree_list(k, ordered, leaves, dm, 0)
        upperlist = _upper_subtree_list(k, ordered, leaves, dm, 0)

        # Initialize the list of edge parent nodes for computation.
        # Here, we start with setting the initial length value of the edge
        # defined by the root and its unique descendant to 0.
        edge_list = [[0, 0, 0]]
        for a in edge_list:
            if a[1] != 0:
                continue
            a1 = 2 * a[0] + 1
            a2 = a1 + 1
            value1 = _edge_attachment_length(
                a1, lowerlist, upperlist, ordered, leaves, adm
            )
            value2 = _edge_attachment_length(
                a2, lowerlist, upperlist, ordered, leaves, adm
            )
            a1_id = ordered[1, np.where(ordered[0] == a1)[0][0]]
            a2_id = ordered[1, np.where(ordered[0] == a2)[0][0]]
            e1 = [a1, a1_id, value1 + a[2]]
            e2 = [a2, a2_id, value2 + a[2]]
            edge_list.append(e1)
            edge_list.append(e2)

        # find the edge with minimum length after edge attachment
        minimum_child = sorted(edge_list, key=itemgetter(2))[0][:2]

        # attach new taxa to the edge
        subtree_nodes = _subtree(minimum_child[0], leaves[0])

        # check intermediates
        # print("data is: ")
        # print(tree_array)
        # print(subtree_nodes)
        # print(minimum_child[0])

        _move_subtree(tree, subtree_nodes, minimum_child[0], 2 * minimum_child[0] + 1)
        tree[0, k - 1] = 2 * minimum_child[0] + 2
        tree[1, k - 1] = k

    return tree


def gme(dm, allow_edge_estimation=True):
    r"""Perform greedy minimum evolution (GME) for phylogenetic reconstruction.

    Parameters
    ----------
    dm : skbio.DistanceMatrix
        Input distance matrix containing distances between taxa.
    allow_edge_estimation : bool, optional
        Whether to perform an OLS-based estimation of edge values (``True``, default)
        or return a tree without edge values assigned (``False``).

    Returns
    -------
    TreeNode
        Reconstructed phylogenetic Tree with estimated edge
        values (if ``allow_edge_estimation`` is ``True``).

    Notes
    -----
    Greedy Minimum Evolution (GME) is a distance-based algorithm for phylogenetic
    reconstruction utilizing the minimum evolution principle for selecting a tree
    topology with the lowest sum of branch lengths according to a given method
    of estimating branch lengths. Ordinary Least Squares (OLS) is a natural
    framework for edge estimation as it is statistically consistent with
    minimum evolution and is used for GME.

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
    >>> from skbio.tree import gme

    >>> dm = DistanceMatrix([[0, 0.02,  0.18,  0.34,  0.55],
    ...                      [0.02,  0, 0.19, 0.35,  0.55],
    ...                      [0.18, 0.19,  0,  0.34,  0.54],
    ...                      [0.34, 0.35,  0.34,  0,  0.62],
    ...                      [0.55,  0.55,  0.54,  0.62,  0]],
    ...                      ['human','monkey','pig','rat','chicken'])

    Perform Greedy Minimum Evoltuion (GME) and construct the minimum evolution tree
    representing the relationship between those taxa. This is returned as a TreeNode
    object.

    >>> tree = gme(dm)
    >>> print(tree.ascii_art())
                        /-monkey
    -human--- /--------|
                       |          /-pig
                        \--------|
                                 |          /-rat
                                  \--------|
                                            \-chicken

    Notice that, unlike neighbor joining, the tree is rooted at a taxa/leaf node.
    This will allow it to have nearest neighbor interchange performed on it without
    needing to re-root the tree.

    """
    if dm.shape[0] < 3:
        raise ValueError(
            "Distance matrix must be at least 3x3 to "
            "generate a minimum evolution tree."
        )

    # Perform GME and get an array representing the tree.
    tree = _gme(dm.data)

    # Create TreeNode object from the tree array.
    leaves = tree[:, tree[0] != 0]
    ordered = _ordered(tree)
    tree = _array_to_tree(dm.ids, tree)

    # OLS-based edge estimation
    if allow_edge_estimation:
        _ols_edge(tree, dm, ordered, leaves, 0)

    return tree


def _tip_or_root(node, ordered, tips):
    """Get the taxon index(es) of a node if it's a tip or root, otherwise its
    descending tips.

    Parameters
    ----------
    node : int
        Index of node in the tree.
    ordered : (2, 2 * N - 1) ndarray
        Indices of all nodes in the tree.
        Indices of all nodes in the distance matrix, or 0 for internal nodes.
    tips : (2, N) ndarray
        Indices of tips in the tree.
        Indices of tips in the distance matrix.

    Returns
    -------
    list of int
        Indices of tips.

    """
    idx = ordered[1, np.where(ordered[0] == node)[0][0]]

    # tip
    if idx > 0:
        return [idx]

    # internal node
    idxlst = []
    for x in tips.T:
        if num_dist_cy(node, x[0]) > 0:
            idxlst.append(x[1])
    return idxlst


def _tip_or_root_upper(node, ordered, tips, root):
    """Get name(s) of a node if it's a tip or root, otherwise
    its ascending tips including root.
    """
    idx = ordered[1, np.where(ordered[0] == node)[0][0]]
    if idx > 0:
        return [idx]
    idxlst = []
    for x in tips.T:
        if num_dist_cy(node, x[0]) < 0:
            idxlst.append(x[1])
    return idxlst + [root]


def _average_distance(node1, node2, dm, ordered, tips):
    """Return the average distance between the leaves of two subtrees.

    Distances between nodes are calculated using a distance matrix.

    """
    nodelist1 = _tip_or_root(node1, ordered, tips)
    nodelist2 = _tip_or_root(node2, ordered, tips)
    return dm[nodelist1][:, nodelist2].mean()


def _average_distance_upper(node1, node2, dm, ordered, tips, root):
    """Return the average distance between the tips of two subtrees.

    Used for subtrees which have a set of tips that are the complement
    of the set of tips that are descendants from the node defining
    the subtree.

    Given an internal edge of a binary tree, exactly one adjacent edge
    will connect to a node defining a subtree of this form.

    """
    nodelist1 = _tip_or_root(node1, ordered, tips)
    nodelist2 = _tip_or_root_upper(node2, ordered, tips, root)
    return dm[nodelist1][:, nodelist2].mean()


def _subtree_count(subtree, tips):
    """Return the number of tips in a subtree."""
    count = 0
    for x in tips[0]:
        if num_dist_cy(subtree, x) > 0:
            count += 1
    return count or 1


def _average_subtree_distance(a, b, a1, a2, dm, ordered, tips):
    """Return the average distance between two subtrees."""
    return (
        _subtree_count(a1, tips) * _average_distance(a1, b, dm, ordered, tips)
        + _subtree_count(a2, tips) * _average_distance(a2, b, dm, ordered, tips)
    ) / _subtree_count(a, tips)


def _average_distance_k(node_k, nodesubtree, ordered, leaves, dm, root, upper=False):
    """Return the average distance between a subtree, defined
    by a node and its descendants, and a taxa to attach.

    The subtree here is referred to as a lower subtree.

    """
    if upper:
        nodelistsubtree = _tip_or_root_upper(nodesubtree, ordered, leaves, root)
    else:
        nodelistsubtree = _tip_or_root(nodesubtree, ordered, leaves)
    return dm[node_k, nodelistsubtree].mean()


def _lower_subtree_list(k, ordered, leaves, dm, root):
    """Return the list of values representing the change in tree length
    after attaching a taxa at an edge not the root edge. The subtree
    involved is a lower subtree.

    This list is indexed by the tree's postorder.

    """
    lower_list = []
    for node in ordered.T:
        if node[1] != 0:
            lower_list.append(dm[node[1], k])
        else:
            node1 = 2 * node[0] + 1
            node2 = node1 + 1
            average = (
                _subtree_count(node1, leaves)
                * _average_distance_k(k, node1, ordered, leaves, dm, root)
                + _subtree_count(node2, leaves)
                * _average_distance_k(k, node2, ordered, leaves, dm, root)
            ) / _subtree_count(node[0], leaves)
            lower_list.append(average)
    return lower_list


def _upper_subtree_list(k, ordered, leaves, dm, root):
    """Return the list of values representing the change in tree length
    after attaching a taxa at an edge not the root edge. The subtree
    involved is an upper subtree.

    This list is indexed by the tree's postorder.

    """
    taxa = len(leaves[0]) + 1
    upper_list = []
    for node in ordered.T:
        if node[0] == 0:
            upper_list.append(dm[root, k])
        else:
            p = _parent(node[0])
            s = _sibling(node[0])
            average = (
                (taxa - _subtree_count(p, leaves))
                * _average_distance_k(k, p, ordered, leaves, dm, root, upper=True)
                + _subtree_count(s, leaves)
                * _average_distance_k(k, s, ordered, leaves, dm, root)
            ) / (taxa - _subtree_count(node[0], leaves))
            upper_list.append(average)
    return upper_list


def _average_distance_matrix(dm, ordered, leaves, root):
    """Return the matrix of distances between pairs of subtrees.

    Parameters
    ----------
    dm : (N, N) ndarray
        Distance matrix.
    ordered : (2, 2 * N - 1) ndarray
        Indices of all nodes in the tree.
    leaves : (2, N) ndarray
        Indices of leaf nodes in the tree.
        Indices of leaf nodes in the distance matrix.
    root : int
        Index of root node in the tree.

    Returns
    -------
    (2 * N - 1, 2 * N - 1) ndarray
        Average distance matrix.
    """
    n = len(ordered[0])
    taxa_size = len(leaves[0]) + 1
    adm = np.empty((n, n))
    for i, a in enumerate(ordered.T):
        ancestors = _ancestors(a[0], ordered[0])

        # skip over unique descendant
        # TO-DO: there is only one instance which is always at the end, therefore this
        # check can be skipped and the for loop should be ordered.T[:-1].
        if a[0] == 0:
            continue

        # find the average distance between given node and root
        if a[1] != 0:
            adm[n - 1, i] = adm[i, n - 1] = dm[a[1], 0]
        else:
            a1 = 2 * a[0] + 1
            a2 = a1 + 1
            adm[n - 1, i] = adm[i, n - 1] = (
                _subtree_count(a1, leaves)
                * _average_distance_upper(a1, 0, dm, ordered, leaves, 0)
                + _subtree_count(a2, leaves)
                * _average_distance_upper(a2, 0, dm, ordered, leaves, 0)
            ) / _subtree_count(a[0], leaves)

        # find the average distance between first node and a second node
        # which is above the first node in the postorder as well as an ancestor
        for j in range(i + 1, n - 1):  # part (a)
            b = ordered[:, j]

            # b is an ancestor of a
            if b[0] in ancestors:
                # TO-DO: stop before the last element (see above)
                if b[0] == 0:
                    continue
                s = _sibling(b[0])
                p = _parent(b[0])
                adm[i, j] = adm[j, i] = (
                    _subtree_count(s, leaves)
                    * _average_distance(a[0], s, dm, ordered, leaves)
                    + (taxa_size - _subtree_count(p, leaves))
                    * _average_distance_upper(a[0], p, dm, ordered, leaves, 0)
                ) / (taxa_size - _subtree_count(b[0], leaves))

            # all other situations
            else:
                # both nodes are tips
                if a[1] != 0 and b[1] != 0:
                    adm[i, j] = adm[j, i] = dm[a[1], b[1]]
                # second node is a tip, but not the first node
                elif b[1] != 0:
                    a1 = 2 * a[0] + 1
                    a2 = a1 + 1
                    adm[i, j] = adm[j, i] = _average_subtree_distance(
                        a[0], b[0], a1, a2, dm, ordered, leaves
                    )
                # neither node is a tip
                else:
                    b1 = 2 * b[0] + 1
                    b2 = b1 + 1
                    adm[i, j] = adm[j, i] = _average_subtree_distance(
                        b[0], a[0], b1, b2, dm, ordered, leaves
                    )

    # zero the diagonal
    np.fill_diagonal(adm, 0)

    return adm


def _edge_attachment_length(child, lowerlist, upperlist, ordered, leaves, adm):
    """Return the change in the length of a tree after attaching a new taxa
    at an edge in comparison to attaching at the root edge.

    The value of attaching at the root edge is considered to be 0.

    """
    taxa = len(leaves[0]) + 1
    parent = _parent(child)
    sibling = _sibling(child)
    for node in ordered.T:
        if node[0] == parent:
            c_i = np.where(ordered[0] == node[0])[0][0]
            c_size = taxa - _subtree_count(parent, leaves)
        elif node[0] == child:
            b_i = np.where(ordered[0] == node[0])[0][0]
            b_size = _subtree_count(child, leaves)
        elif node[0] == sibling:
            a_i = np.where(ordered[0] == node[0])[0][0]
            a_size = _subtree_count(sibling, leaves)
    lambda_1 = (a_size + b_size * c_size) / ((a_size + b_size) * (c_size + 1))
    lambda_2 = (a_size + b_size * c_size) / ((a_size + c_size) * (b_size + 1))
    length = 0.5 * (
        (lambda_1 - lambda_2) * (lowerlist[a_i] + adm[b_i, c_i])
        + (lambda_2 - 1) * (adm[a_i, b_i] + upperlist[c_i])
        + (1 - lambda_1) * (adm[a_i, c_i] + lowerlist[b_i])
    )
    return length


def _ols_edge(tree, dm, ordered, leaves, root):
    """Assign estimated edge values to a tree based on a given distance matrix.

    Estimation of edge values is based on an ordinary least squares (OLS) framework.

    """
    adm = _average_distance_matrix(dm, ordered, leaves, root)
    ordered = list(tree.postorder(include_self=False))
    root = tree.root()
    taxa_size = root.count(tips=True) + 1

    ### identify edges by first finding the child node of an edge
    for edge_node in ordered:
        parent = edge_node.parent
        ### skip over root node
        if edge_node.is_root():
            continue
        ### calculate edge length for the edge adjacent to the root
        elif parent.is_root():
            for index, node in enumerate(ordered):
                if node == edge_node:
                    i1 = index
            a, b = edge_node.children
            for index, node in enumerate(ordered):
                if node == a:
                    i2 = index
                elif node == b:
                    i3 = index
            edge_node.length = 0.5 * (adm[i2, i1] + adm[i3, i1] - adm[i2, i3])
        ### calculate edge lengths for external edges
        elif edge_node.is_tip():
            a = parent.parent
            if a.is_root():
                for child in a.children:
                    a = child
            for siblingnode in edge_node.siblings():
                b = siblingnode
                for index, node in enumerate(ordered):
                    if node is edge_node:
                        i1 = index
                    elif node is a:
                        i2 = index
                    elif node is b:
                        i3 = index
            edge_node.length = 0.5 * (adm[i2][i1] + adm[i3][i1] - adm[i2][i3])
        ### calculate edge lengths for internal edges
        else:
            a = parent.parent
            if a.is_root():
                for child in a.children:
                    a = child
            for index, node in enumerate(ordered):
                if node == a:
                    i1 = index
            c, d = edge_node.children
            for sibling in edge_node.siblings():
                b = sibling
            for index, node in enumerate(ordered):
                if node is b:
                    i2 = index
                elif node is c:
                    i3 = index
                elif node is d:
                    i4 = index
            ### count the tips of subtrees which are adjacent to the internal edge
            sub_tips = []
            for subtree in b, c, d:
                sub_tips.append(1 if subtree.is_tip() else subtree.count(tips=True))
            b_, c_, d_ = sub_tips
            a_ = taxa_size - b_ - c_ - d_
            ## calculate the edge length
            lambda1 = (a_ * d_ + b_ * c_) / ((a_ + b_) * (c_ + d_))
            edge_node.length = 0.5 * (
                (lambda1 * (adm[i1, i3] + adm[i2, i4]))
                + ((1 - lambda1) * (adm[i1, i4] + adm[i2, i3]))
                - (adm[i1, i2] + adm[i3, i4])
            )
