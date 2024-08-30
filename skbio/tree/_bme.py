# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from operator import itemgetter
from skbio.tree import TreeNode
from skbio.tree._gme import _average_distance_k, _average_distance_k_upper


def bme(dm, allow_edge_estimation=True):
    r"""Perform balanced minimum evolution (BME) for phylogenetic reconstruction.

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
    Balanced Minimum Evolution (BME) is a refinement of the distance-based minimum
    evolution problem where average distances between subtrees ignores the size of the
    subtrees. The BME algorithm implemented here uses the same OLS based edge estimation
    used with Greedy Minimum Evolution (GME).

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
    >>> from skbio.tree import bme

    >>> dm = DistanceMatrix([[0, 0.02,  0.18,  0.34,  0.55],
    ...                      [0.02,  0, 0.19, 0.35,  0.55],
    ...                      [0.18, 0.19,  0,  0.34,  0.54],
    ...                      [0.34, 0.35,  0.34,  0,  0.62],
    ...                      [0.55,  0.55,  0.54,  0.62,  0]],
    ...                      ['human','monkey','pig','rat','chicken'])

    Perform Balanced Minimum Evoltuion (BME) and construct the minimum evolution tree
    representing the relationship between those taxa. This is returned as a TreeNode
    object.

    >>> tree = bme(dm)
    >>> print(tree.ascii_art())
                        /-monkey
    -human--- /--------|
                       |          /-pig
                        \--------|
                                 |          /-chicken
                                  \--------|
                                            \-rat

    Notice that, unlike neighbor joining, the tree is rooted at a taxa/leaf node.
    This will allow it to have nearest neighbor interchange performed on it without
    needing to re-root the tree.

    """
    if dm.shape[0] < 3:
        raise ValueError(
            "Distance matrix must be at least 3x3 to "
            "generate a minimum evolution tree."
        )

    # count the number of taxa
    n = len(dm.ids)

    # create initial nodes
    root_name, node_a_name, node_b_name = dm.ids[0], dm.ids[1], dm.ids[2]
    root = TreeNode(root_name)
    root_child = TreeNode("")
    node_a = TreeNode(node_a_name)
    node_b = TreeNode(node_b_name)
    root.append(root_child)
    root_child.append(node_a)
    root_child.append(node_b)

    # loop over rest of nodes for k = 4 to n
    for k in range(3, n):
        taxa = root.count(tips=True) + 1
        # calculate average distance matrix between subtrees of T_(k-1)
        # Note that this could be replaced with an updating step to an
        # initially computed average distance matrix, as noted in GME,
        # however the updating step requires all subtrees to be updated
        # in contrast to GME, so instead the matrix is recomputed after
        # every iteration.
        adm = _balanced_average_matrix(root, dm)
        # create node for taxa k
        node_k_name = dm.ids[k]
        connecting_node = TreeNode("")
        node_k = TreeNode(node_k_name)
        connecting_node.append(node_k)
        # connecting_node = TreeNode.read(["(%s);" % (node_k_name)])
        # node_k = connecting_node.find(node_k_name)
        # create average distance lists for subtrees of T_(k-1)
        ordered = list(root.postorder(include_self=False))
        lowerlist = _balanced_lower(ordered, node_k, dm)
        upperlist = _balanced_upper(ordered, node_k, taxa, dm)
        # Initialize the list of edge parent nodes for computation.
        # The tuple is a pair of a TreeNode object and the value of the tree
        # length for attaching at the edge consisting of the node and its parent.
        # Here, we start with setting the initial length value of the edge
        # defined by the root and its unique descendant to 0.
        edge_list = [(root_child, 0)]
        for a, i in edge_list:
            if a.is_tip():
                continue
            a1, a2 = a.children
            value1 = _balanced_attach_length(
                a1, lowerlist, upperlist, ordered, taxa, adm
            )
            value2 = _balanced_attach_length(
                a2, lowerlist, upperlist, ordered, taxa, adm
            )
            e1 = (a1, value1 + i)
            e2 = (a2, value2 + i)
            edge_list.append(e1)
            edge_list.append(e2)
        # find the edge with minimum length after edge attachment
        minimum_child = sorted(edge_list, key=itemgetter(1))[0][0]
        # attach new taxa to the edge
        minimum_parent = minimum_child.parent
        minimum_parent.append(connecting_node)
        connecting_node.append(minimum_child)
    if allow_edge_estimation:
        _bal_ols_edge(root, dm)
    return root


def _balanced_lower(ordered, node, dm):
    """Return the list of values representing the change in tree length
    after attaching a taxa at an edge not the root edge. The subtree
    involved is a lower subtree.

    This list is indexed by the tree's postorder.

    """
    lower_list = []
    for i, a in enumerate(ordered):
        if a.is_tip():
            lower_list.append(dm[a.name, node.name])
        else:
            a1, a2 = a.children
            average = (
                _average_distance_k(node, a1, dm) + _average_distance_k(node, a2, dm)
            ) * 0.5
            lower_list.append(average)
    return lower_list


def _balanced_upper(ordered, node, taxa, dm):
    """Return the list of values representing the change in tree length
    after attaching a taxa at an edge not the root edge. The subtree
    involved is an upper subtree.

    This list is indexed by the tree's postorder.

    """
    upper_list = []
    for a in ordered:
        if a.parent.is_root():
            upper_list.append(dm[a.parent.name, node.name])
        else:
            p = a.parent
            s = a.siblings()[0]
            average = (
                _average_distance_k_upper(node, p, dm)
                + _average_distance_k(node, s, dm)
            ) * 0.5
            upper_list.append(average)
    return upper_list


def _balanced_attach_length(child, lowerlist, upperlist, ordered, taxa, adm):
    """Return the change in the length of a tree after attaching a new taxa
    at an edge in comparison to attaching at the root edge.

    The value of attaching at the root edge is considered to be 0. This uses
    the balanced framework for average distance.

    """
    parent = child.parent
    sibling = child.siblings()[0]
    for i, a in enumerate(ordered):
        if a is parent:
            c_i = i
        elif a is child:
            b_i = i
        elif a is sibling:
            a_i = i
    length = 0.25 * (
        (adm[a_i, c_i] + lowerlist[b_i]) - (adm[a_i, b_i] + upperlist[c_i])
    )
    return length


def _balanced_average_matrix(tree, dm):
    """Return the matrix of distances between pairs of subtrees.

    Here, the definition of average distance coincides with the balanced minimum
    evolution problem.

    """
    # start from unique descendant
    name = tree.name
    (tree,) = tree.children
    ordered = list(tree.postorder(include_self=False))
    n = len(ordered)
    adm = np.zeros((n + 1, n + 1))
    for i, a in enumerate(ordered):
        ancestors = a.ancestors()
        # find the average distance between given node and root
        if a.is_tip():
            adm[n, i] = adm[i, n] = dm[a.name, name]
        else:
            a1, a2 = a.children
            for k, x in enumerate(ordered):
                if x is a1:
                    i1 = k
                if x is a2:
                    i2 = k
            adm[n, i] = adm[i, n] = (adm[n, i1] + adm[n, i2]) * 0.5
        # find the average distance between first node and a second node
        # which is above the first node in the postorder as well as an ancestor
        for j in range(i + 1, n):  # part (a)
            b = ordered[j]
            # skipping over ancestors
            if b in ancestors:
                continue
            # both nodes are tips
            if a.is_tip() and b.is_tip():
                adm[i, j] = adm[j, i] = dm[a.name, b.name]
            # second node is a tip, but not the first node
            elif b.is_tip():
                a1, a2 = a.children
                for k, x in enumerate(ordered):
                    if x is a1:
                        i1 = k
                    if x is a2:
                        i2 = k
                adm[i, j] = adm[j, i] = (adm[j, i1] + adm[j, i2]) * 0.5
            # neither node is a tip
            else:
                b1, b2 = b.children
                for k, x in enumerate(ordered):
                    if x is b1:
                        j1 = k
                    if x is b2:
                        j2 = k
                adm[i, j] = adm[j, i] = (adm[j1, i] + adm[j2, i]) * 0.5
    # calculating for second nodes which are ancestors
    for j, b in enumerate(ordered):
        (s,) = b.siblings()
        p = b.parent
        for k, x in enumerate(ordered):
            if x is s:
                j_s = k
            if x is p:
                j_p = k
        for i, a in enumerate(ordered):
            if b in a.ancestors():
                adm[i, j] = adm[j, i] = (adm[i, j_s] + adm[i, j_p]) * 0.5
    return adm


def _bal_ols_edge(tree, dm):
    """Assign estimated edge values to a tree based on a given distance matrix.

    Estimation of edge values is based on an ordinary least squares (OLS) framework.
    Average distances defined here coincide with the balanced minimum evolution problem.

    """
    adm = _balanced_average_matrix(tree, dm)
    ordered = list(tree.postorder(include_self=False))
    root = tree.root()
    # identify edges by first finding the child node of an edge
    for edge_node in ordered:
        parent = edge_node.parent
        # skip over root node
        if edge_node.is_root():
            continue
        # calculate edge length for the edge adjacent to the root
        elif parent.is_root():
            for index, node in enumerate(ordered):
                if node is edge_node:
                    i1 = index
            a, b = edge_node.children
            for index, node in enumerate(ordered):
                if node is a:
                    i2 = index
                elif node is b:
                    i3 = index
            edge_node.length = 0.5 * (adm[i2, i1] + adm[i3, i1] - adm[i2, i3])
        # calculate edge lengths for external edges
        elif edge_node.is_tip():
            a = parent.parent
            if a.is_root():
                for child in a.children:
                    a = child
            b = edge_node.siblings()[0]
            for index, node in enumerate(ordered):
                if node is edge_node:
                    i1 = index
                elif node is a:
                    i2 = index
                elif node is b:
                    i3 = index
            edge_node.length = 0.5 * (adm[i2][i1] + adm[i3][i1] - adm[i2][i3])
        # calculate edge lengths for internal edges
        else:
            a = parent.parent
            if a.is_root():
                a = a.children[0]
            for index, node in enumerate(ordered):
                if node is a:
                    i1 = index
            c, d = edge_node.children
            b = edge_node.siblings()[0]
            for index, node in enumerate(ordered):
                if node is b:
                    i2 = index
                elif node is c:
                    i3 = index
                elif node is d:
                    i4 = index
            # calculate the edge length
            edge_node.length = 0.5 * (
                (0.5 * (adm[i1, i3] + adm[i2, i4]))
                + (0.5 * (adm[i1, i4] + adm[i2, i3]))
                - (adm[i1, i2] + adm[i3, i4])
            )
