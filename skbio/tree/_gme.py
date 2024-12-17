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
                                 |          /-chicken
                                  \--------|
                                            \-rat

    Notice that, unlike neighbor joining, the tree is rooted at a taxa/leaf node.
    This will allow it to have nearest neighbor interchange performed on it without
    needing to re-root the tree.

    """
    if (n := dm.shape[0]) < 3:
        raise ValueError("Distance matrix must be at least 3x3 to generate a tree.")

    # extract taxa and matrix data
    taxa = dm.ids
    dm = dm.data

    # create initial triplet tree: ((left,right)stem)root;
    # - root (r, taxon 0)
    # - stem (d): a unique descendant of root
    # - left (taxon 1) and right (taxon 2): children of stem
    root, left, right = [TreeNode(i) for i in range(3)]
    stem = TreeNode(children=[left, right])
    root.append(stem, uncache=False)

    sortkey = itemgetter(1)

    # Iterate over taxa 3 to n - 1. In each literation, add taxon k to the tree.
    for k in range(3, n):
        n_taxa = k + 1

        # Calculate average distance matrix between subtrees of T_(k-1).
        # Note that this should be replaced with a step to update an initial
        # computation of the average distance matrix rather than to recompute
        # the matrix after every iteration as implemented here.
        adm = _average_distance_matrix(root, dm)

        # create node for taxon k
        node_k = TreeNode(k)
        connector = TreeNode()
        connector.append(node_k, uncache=False)

        # create average distance lists for subtrees of T_(k-1)
        # TODO: avoid repeated postorder traversal
        ordered = list(root.postorder(include_self=False))

        lowerlist = _lower_subtree_list(ordered, node_k, dm)
        upperlist = _upper_subtree_list(ordered, node_k, n_taxa, dm)

        # Initialize the list of edge parent nodes for computation.
        # The tuple is a pair of a TreeNode object and the value of the tree
        # length for attaching at the edge consisting of the node and its parent.
        # Here, we start with setting the initial length value of the edge
        # defined by the root and its unique descendant to 0.
        edge_list = [(stem, 0)]
        for a, i in edge_list:
            if a.children:
                a1, a2 = a.children
                value1 = _edge_attachment_length(
                    a1, lowerlist, upperlist, ordered, n_taxa, adm
                )
                value2 = _edge_attachment_length(
                    a2, lowerlist, upperlist, ordered, n_taxa, adm
                )
                edge_list.extend([(a1, value1 + i), (a2, value2 + i)])

        # find the edge with minimum length after edge attachment
        min_child = sorted(edge_list, key=sortkey)[0][0]

        # attach new taxon to the edge
        min_parent = min_child.parent
        min_parent.append(connector, uncache=False)
        connector.append(min_child, uncache=False)

    # estimate branch lengths
    if allow_edge_estimation:
        _ols_edge(root, dm)

    # replace taxon indices with taxon names, and clean up intermediate attributes
    root.name = taxa[0]
    for node in root.children[0].traverse(include_self=True):
        if not node.children:
            node.name = taxa[node.name]
        delattr(node, "taxa")
        delattr(node, "size")
        delattr(node, "i")

    return root


def _average_distance_k(nodek, nodesubtree, dm):
    """Return the average distance between a subtree, defined
    by a node and its descendants, and a taxa to attach.

    The subtree here is referred to as a lower subtree.

    """
    nodelistk = [nodek.name]
    nodelistsubtree = _tip_or_root(nodesubtree)
    return dm[nodelistk][:, nodelistsubtree].mean()


def _average_distance_k_upper(nodek, nodesubtree, dm):
    """Return the average distance between a subtree, defined
    by all nodes that are not descendants, and a taxa to attach.

    The subtree here is referred to as an upper subtree.

    """
    nodelistk = [nodek.name]
    if nodesubtree.is_root():
        nodelistsubtree = []
    else:
        root = nodesubtree.root()
        nodelistsubtree = [root.name]
        nodelistsubtree.extend(sorted(root.subset() - nodesubtree.subset()))
    return dm[nodelistk][:, nodelistsubtree].mean()
    # df = dm.between(nodelistk, nodelistsubtree)
    # return df["value"].mean()


def _lower_subtree_list(ordered, node, dm):
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
                _subtree_count(a1) * _average_distance_k(node, a1, dm)
                + _subtree_count(a2) * _average_distance_k(node, a2, dm)
            ) / _subtree_count(a)
            lower_list.append(average)
    return lower_list


def _upper_subtree_list(ordered, node, taxa, dm):
    """Return the list of values representing the change in tree length
    after attaching a taxa at an edge not the root edge. The subtree
    involved is an upper subtree.

    This list is indexed by the tree's postorder.

    """
    upper_list = []
    for a in enumerate(ordered):
        a = a[1]
        if a.parent.is_root():
            upper_list.append(dm[a.parent.name, node.name])
        else:
            p = a.parent
            s = a.siblings()[0]
            average = (
                (taxa - _subtree_count(p)) * _average_distance_k_upper(node, p, dm)
                + _subtree_count(s) * _average_distance_k(node, s, dm)
            ) / (taxa - _subtree_count(a))
            upper_list.append(average)
    return upper_list


def _edge_attachment_length(child, lowerlist, upperlist, ordered, taxa, adm):
    """Return the change in the length of a tree after attaching a new taxa
    at an edge in comparison to attaching at the root edge.

    The value of attaching at the root edge is considered to be 0.

    """
    parent = child.parent
    sibling = child.siblings()[0]
    for i, a in enumerate(ordered):
        if a is parent:
            c_i = i
            c_size = taxa - _subtree_count(parent)
        elif a is child:
            b_i = i
            b_size = _subtree_count(child)
        elif a is sibling:
            a_i = i
            a_size = _subtree_count(sibling)
    lambda_1 = (a_size + b_size * c_size) / ((a_size + b_size) * (c_size + 1))
    lambda_2 = (a_size + b_size * c_size) / ((a_size + c_size) * (b_size + 1))
    length = 0.5 * (
        (lambda_1 - lambda_2) * (lowerlist[a_i] + adm[b_i, c_i])
        + (lambda_2 - 1) * (adm[a_i, b_i] + upperlist[c_i])
        + (1 - lambda_1) * (adm[a_i, c_i] + lowerlist[b_i])
    )
    return length


def _average_distance(node1, node2, dm):
    """Return the average distance between the leaves of two subtrees.

    Distances between nodes are calculated using a distance matrix.

    """
    nodelist1 = _tip_or_root(node1)
    nodelist2 = _tip_or_root(node2)
    return dm[nodelist1][:, nodelist2].mean()
    # df = dm.between(nodelist1, nodelist2)
    # return df["value"].mean()


def _tip_or_root(node):
    """Get name(s) of a node if it's a tip or root, otherwise its descending tips."""
    if node.is_tip() or node.is_root():
        return [node.name]
    else:
        return [x.name for x in node.tips()]


def _average_distance_upper(node1, node2, dm):
    """Return the average distance between the leaves of two subtrees.

    Used for subtrees which have a set of tips that are the complement
    of the set of tips that are descendants from the node defining
    the subtree.

    Given an internal edge of a binary tree, exactly one adjacent edge
    will connect to a node defining a subtree of this form.

    """
    nodelist1 = _tip_or_root(node1)
    if node2.is_root():
        nodelist2 = []
    # Second subtree serves as the tree with a set of tips
    # complementary to the set of tips that descend from the
    # corresponding second node.
    else:
        root2 = node2.root()
        nodelist2 = [root2.name]
        nodelist2.extend(sorted(root2.subset() - node2.subset()))
    return dm[nodelist1][:, nodelist2].mean()
    # df = dm.between(nodelist1, nodelist2)
    # return df["value"].mean()


def _subtree_count(subtree):
    """Return the number of leaves in a subtree.

    Assumes the root as a leaf node.

    """
    if subtree.is_tip() or subtree.is_root():
        return 1
    else:
        return subtree.count(tips=True)


def _average_subtree_distance(a, b, a1, a2, dm):
    """Return the average distance between two subtrees."""
    return (
        _subtree_count(a1) * _average_distance(a1, b, dm)
        + _subtree_count(a2) * _average_distance(a2, b, dm)
    ) / _subtree_count(a)


def _average_distance_matrix_eq1(tree, dm):
    """Calculate a matrix of average distances between pairs of subtrees.

    This algorithm produces the same result as `_average_distance_matrix`. However, it
    calculates all subtree-to-subtree distances based on the original taxon-to-taxon
    distances, as discussed in Eq. 1 of Desper and Gascuel (2002). instead of adopting
    a recursive strategy (Eq. 2). Therefore, it is slower. It is implemented as a
    reference, but not used.

    """
    ordered = []
    full = [0]
    for node in tree.postorder(include_self=False):
        if not node.children:
            full.append(node.name)
            node.taxa = frozenset([node.name])
        else:
            node.taxa = frozenset().union(*[x.taxa for x in node.children])
        ordered.append(node.taxa)

    ordered[-1] = frozenset([0])
    full = frozenset(full)

    n = len(ordered)
    adm = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = ordered[i], ordered[j]
            # a is a descendant (proper subset) of b
            if a < b:
                b = full - b
            # a is an ancestor (proper superset) of b
            elif a > b:
                a = full - a
            adm[i, j] = adm[j, i] = dm[list(a)][:, list(b)].mean()

    return adm


def _average_distance_matrix(tree, dm):
    """Calculate a matrix of average distances between pairs of subtrees.

    The tree has a root (taxon 0), its unique descendant (stem), and a strictly binary
    tree below the stem.

    For a tree with k taxa (including root), the resulting distance matrix should have
    2 * k - 3 elements (nodes except for root).

    This algorithm is implemented according to Appendix 4 of Desper and Gascuel (2002).

    """
    # Full set of taxa in the current tree (not all taxa in the input distance matrix).
    # They correspond to indices in the original distance matrix (dm).
    full = [0]

    # Ordered list of nodes. Their indices correspond to the indices in the resulting
    # average distance matrix (adm).
    ordered = []

    # index of the current node in adm
    i = 0

    # Perform postorder traversal of the tree.
    # Attributes are assigned to each node:
    # - i: Index of the node in adm.
    # - size: Number of taxa (tips) descending from the node.
    # - taxa: Indices of descending taxa in dm.
    for node in tree.postorder(include_self=False):
        if not node.children:  # tip (taxon)
            full.append(node.name)
            node.taxa = frozenset([node.name])
            node.size = 1
        else:  # internal node
            node.taxa = frozenset().union(*[x.taxa for x in node.children])
            node.size = sum(x.size for x in node.children)
        ordered.append(node)
        node.i = i
        i += 1

    # The last node must be the stem, therefore it has a single taxon 0.
    tree.children[0].taxa = frozenset([0])

    n = len(ordered)  # number of nodes
    m = len(full)  # number of taxa

    # Initiate the resulting average distance matrix
    adm = np.zeros((n, n), dtype=float)

    # Calculate the average distance between each pair of subtrees defined by nodes
    # a and b (A4.1 of the paper)

    # Step 1: Calculate non-nested subtree to subtree distances (i.e., one is not an
    # ancestor of another. Therefore each subtree is the lower (descending) tree of
    # the node. (A4.1 (a))

    # Loop over nodes in post order.
    # Skip the stem, which is always the last node in the postorder traversal.
    for i in range(n - 1):
        a = ordered[i]

        # check if a is a tip
        if a.size == 1:
            (a_taxon,) = a.taxa
        else:
            a_taxon = 0  # a can never be the root (taxon 0), therefore 0 means none

        # Loop over remaining nodes in post order.
        # The original paper says looping over all nodes (vague), but since distances
        # are symmetric, one can loop over remaining nodes
        for j in range(i + 1, n - 1):
            b = ordered[j]

            # TODO: The original paper suggests a postorder traversal that skips
            # subtree a. However, this is not implemented here. Therefore the current
            # code performs traversal over the entire tree, and select non-nested
            # nodes.

            # exclude nested nodes (i.e., b is ancestral to or descending from a)
            if not b.taxa.isdisjoint(a.taxa):
                continue

            # If both a and b are tips, take the original taxon-to-taxon distance
            # (A4.1 (a) i).
            if b.size == 1 and a_taxon:
                (b_taxon,) = b.taxa
                dist = dm[a_taxon, b_taxon]

            # If a is an internal node, and b is either (a tip or an internal node),
            # calculate the average distance based on the two child subtrees of a
            # (A4.1 (a) ii).
            elif not a_taxon:
                a1, a2 = a.children
                dist = (a1.size * adm[a1.i, b.i] + a2.size * adm[a2.i, b.i]) / a.size

            # If a is a tip, and b is an internal node, calculate the average distance
            # based on the two child subtrees of b (A4.1 (a) iii).
            else:
                b1, b2 = b.children
                dist = (b1.size * adm[a.i, b1.i] + b2.size * adm[a.i, b2.i]) / b.size

            adm[i, j] = adm[j, i] = dist

    # Step 2: Calculate subtree to root (taxon 0) distances (A4.1 (b)).

    # This is done through a postorder traversal
    for i, a in enumerate(ordered[:-1]):
        if a.size == 1:
            (a_taxon,) = a.taxa
            dist = dm[0, a_taxon]
        else:
            a1, a2 = a.children
            dist = (a1.size * adm[a1.i, n - 1] + a2.size * adm[a2.i, n - 1]) / a.size

        adm[a.i, n - 1] = adm[n - 1, a.i] = dist

    # Step 3: Calculate nested subtree to subtree distances, in which the first node
    # (a) is a descendant of the second node (b), therefore the first subtree (A) is
    # the lower (descending) tree of node a, whereas the second subtree (B) is the
    # upper (ancestral) tree of node b (A4.1 (c)).

    # This is done through a preorder traversal.
    for b in tree.children[0].preorder(include_self=False):
        # The upper subtree of b consists of two child subtrees: its parent and its
        # sibling.
        p = b.parent
        s = p.children[0] if b is p.children[1] else p.children[1]

        # The size of (upper) subtree b the complement of its descendants. Same for
        # the parent subtree.
        b_size = m - b.size
        p_size = m - p.size
        s_size = s.size

        # Iterate over all subtrees below b.
        # Likewise, the size of subtree b is the complement.
        # The paper says this traversal can be done in any manner. Here, we use the
        # postorder.
        for a in b.postorder(include_self=False):
            dist = (s_size * adm[a.i, s.i] + p_size * adm[a.i, p.i]) / b_size
            adm[a.i, b.i] = adm[b.i, a.i] = dist

    return adm


def _average_distance_matrix_old(tree, dm):
    """Return the matrix of distances between pairs of subtrees.

    For a tree with k taxa (including root), the resulting distance matrix should have
    2 * k - 3 elements (nodes except for root).

    """
    ordered = list(tree.postorder(include_self=False))

    n = len(ordered)  # number of nodes
    m = (n + 3) // 2  # number of taxa

    # adm = np.empty((n, n))
    adm = np.full((n, n), np.nan)

    # Calculate the average distance between each pair of subtrees defined by nodes
    # a and b

    # skip the stem (the unique descendant of root, which is always the last node in
    # the postorder traversal)
    for i, a in enumerate(ordered[:-1]):
        ancestors = a.ancestors()

        # Find the average distance between a given subtree and the root (i.e., the
        # first taxon).
        # When a is a tip (taxon): the average distance is simply the distance from a
        # to root.
        if not a.children:
            adm[n - 1, i] = adm[i, n - 1] = dm[a.name, 0]

        # when a is an internal node,
        else:
            a1, a2 = a.children
            adm[n - 1, i] = adm[i, n - 1] = _average_subtree_distance(
                a, tree, a1, a2, dm
            )

        # find the average distance between first node and a second node
        # which is above the first node in the postorder as well as an ancestor
        for j in range(i + 1, n - 1):  # part (a)
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
                adm[i, j] = adm[j, i] = _average_subtree_distance(a, b, a1, a2, dm)
            # neither node is a tip
            else:
                b1, b2 = b.children
                adm[i, j] = adm[j, i] = _average_subtree_distance(b, a, b1, b2, dm)

    # calculating for second nodes which are ancestors
    for j, b in enumerate(ordered):
        # skipping over unique descendant
        if b in tree.children:
            continue
        s_ = b.siblings()
        for sibling in s_:
            s = sibling
        p = b.parent
        for i, a in enumerate(ordered):
            if b in a.ancestors():
                adm[i, j] = adm[j, i] = (
                    _subtree_count(s) * _average_distance(a, s, dm)
                    + (m - _subtree_count(p)) * _average_distance_upper(a, p, dm)
                ) / (m - _subtree_count(b))

    np.fill_diagonal(adm, 0)
    return adm


def _ols_edge_old(tree, dm):
    """Assign estimated edge values to a tree based on a given distance matrix.

    Estimation of edge values is based on an ordinary least squares (OLS) framework.

    """
    adm = _average_distance_matrix(tree, dm)
    ordered = list(tree.postorder(include_self=False))
    root = tree.root()
    taxa_size = root.count(tips=True) + 1
    # identify edges by first finding the child node of an edge
    for edge_node in ordered:
        parent = edge_node.parent
        # skip over root node
        if edge_node.is_root():
            continue
        # calculate edge length for the edge adjacent to the root
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
        # calculate edge lengths for external edges
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
        # calculate edge lengths for internal edges
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
            # count the tips of subtrees which are adjacent to the internal edge
            sub_tips = []
            for subtree in b, c, d:
                sub_tips.append(1 if subtree.is_tip() else subtree.count(tips=True))
            b_, c_, d_ = sub_tips
            a_ = taxa_size - b_ - c_ - d_
            # calculate the edge length
            lambda1 = (a_ * d_ + b_ * c_) / ((a_ + b_) * (c_ + d_))
            edge_node.length = 0.5 * (
                (lambda1 * (adm[i1, i3] + adm[i2, i4]))
                + ((1 - lambda1) * (adm[i1, i4] + adm[i2, i3]))
                - (adm[i1, i2] + adm[i3, i4])
            )


def _ols_edge(tree, dm, adm=None):
    """Calculate branch lengths of a tree based on a distance matrix of taxa.

    Estimation of branch lengths is based on an ordinary least squares (OLS) framework.

    This algorithm is implemented according to Eqs. 3 & 4 of Desper and Gascuel (2002).

    """
    # calculate average distance matrix if not provided
    if adm is None:
        adm = _average_distance_matrix(tree, dm)

    m = dm.shape[0]
    for node in tree.children[0].traverse(include_self=False):
        p = node.parent
        a, b = p.children
        s = a if node is b else b

        # External (terminal) branch: based on the triplet (iAB) of self (i), parent
        # (A), and sibling (B) (Eq. 4).
        if not node.children:
            node.length = 0.5 * (adm[p.i, node.i] + adm[s.i, node.i] - adm[p.i, s.i])

        # Internal branch: based on the quartet (AB|CD) of parent (A, upper), sibling
        # (B), and children (C and D) (Eq. 3).
        else:
            a, b = node.children
            lambda_ = ((m - p.size) * b.size + s.size * a.size) / (
                (m - p.size + s.size) * (a.size + b.size)
            )
            node.length = 0.5 * (
                lambda_ * (adm[p.i, a.i] + adm[s.i, b.i])
                + (1 - lambda_) * (adm[p.i, b.i] + adm[s.i, a.i])
                - (adm[p.i, s.i] + adm[a.i, b.i])
            )

    # branch connecting stem and root
    node = tree.children[0]
    a, b = node.children
    node.length = 0.5 * (adm[a.i, -1] + adm[b.i, -1] - adm[a.i, b.i])
