from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from collections import defaultdict
from future.builtins import zip

import numpy as np

from skbio.tree import TreeNode


def _walk_clades(trees, weights):
    """Walk all the clades of all the trees

    Parameters
    ----------
    trees : list of TreeNode
        The trees to walk
    weights : np.array
        Tree weights

    Returns
    -------
    list of tuple
        The clades and support values sorted by support value such that the
        most supported clade is index 0. The tuples are of the form:
        (frozenset, float).
    defaultdict(float)
        The edge lengths, keyed by frozenset of the clade, and valued by the
        weighted average length of the clade by the trees the clade was
        observed in.

    """
    clade_counts = defaultdict(float)
    edge_lengths = defaultdict(float)
    total = weights.sum()

    # get clade counts
    tipnames_f = lambda n: [n.name] if n.is_tip() else []
    for tree, weight in zip(trees, weights):
        tree.cache_attr(tipnames_f, 'tip_names', frozenset)

        for node in tree.postorder():
            tip_names = node.tip_names

            # if node.length is not None, fetch it and weight it
            length = node.length * weight if node.length is not None else None

            clade_counts[tip_names] += weight

            if length is None:
                edge_lengths[tip_names] = None
            else:
                edge_lengths[tip_names] += length / total

    # sort clades by number times observed
    clade_counts = sorted(clade_counts.items(), key=lambda x: len(x[0]),
                          reverse=True)

    return clade_counts, edge_lengths


def _filter_clades(clade_counts, cutoff_threshold):
    """Filter clades that not well supported or are contradicted

    Parameters
    ----------
    clade_counts : list of tuple
        Where the first element in each tuple is the frozenset of the clade,
        and the second element is the support value. It is expected that this
        list is sorted by descending order by support.
    cutoff_threshold : float
        The minimum weighted observation count that a clade must have to be
        considered supported.

    Returns
    -------
    dict
        A dict of the accepted clades, keyed by the frozenset of the clade and
        valued by the support value.
    """
    accepted_clades = {}

    for clade, count in clade_counts:
        conflict = False

        if count <= cutoff_threshold:
            continue

        if len(clade) > 1:
            # check the current clade against all the accepted clades to see if
            # it conflicts. A conflict is defined as:
            # 1. the clades are not disjoint
            # 2. neither clade is a subset of the other
            for accepted_clade in accepted_clades:
                intersect = clade.intersection(accepted_clade)
                subset = clade.issubset(accepted_clade)
                superset = clade.issuperset(accepted_clade)

                if intersect and not (subset or superset):
                    conflict = True

        if conflict is False:
            accepted_clades[clade] = count

    return accepted_clades


def _build_trees(clade_counts, edge_lengths, support_attr):
    """Construct the trees with support

    Parameters
    ----------
    clade_counts : dict
        Keyed by the frozenset of the clade and valued by the support
    edge_lengths : dict
        Keyed by the frozenset of the clade and valued by the weighted length
    support_attr : str
        The name of the attribute to hold the support value

    Returns
    -------
    list of TreeNode
        A list of the constructed trees
    """
    nodes = {}
    queue = [(len(clade), clade) for clade in clade_counts]
    while queue:
        # The values within the queue are updated on each iteration, so it
        # doesn't look like an insertion sort will make sense unfortunately
        queue.sort()
        (clade_size, clade) = queue.pop(0)
        new_queue = []

        # search for ancestors of clade
        for (_, ancestor) in queue:
            if clade.issubset(ancestor):
                # update ancestor such that, in the following example:
                # ancestor == {1, 2, 3, 4}
                # clade == {2, 3}
                # new_ancestor == {1, {2, 3}, 4}
                new_ancestor = (ancestor - clade) | frozenset([clade])

                # update references for counts and lengths
                clade_counts[new_ancestor] = clade_counts.pop(ancestor)
                edge_lengths[new_ancestor] = edge_lengths.pop(ancestor)

                ancestor = new_ancestor

            new_queue.append((len(ancestor), ancestor))

        # if the clade is a tip, then we have a name
        if clade_size == 1:
            name = list(clade)[0]
        else:
            name = None

        # the clade will not be in nodes if it is a tip
        children = [nodes.pop(c) for c in clade if c in nodes]
        length = edge_lengths[clade]

        node = TreeNode(children=children, length=length, name=name)
        setattr(node, support_attr, clade_counts[clade])
        nodes[clade] = node

        queue = new_queue

    return list(nodes.values())


def majority_rule(trees, weights=None, cutoff=0.5, support_attr='support'):
    r"""Determines consensus trees from a list of rooted trees

    Parameters
    ----------
    trees : list of TreeNode
        The trees to operate on
    weights : list or np.array of {int, float}, optional
        If provided, the list must be in index order with `trees`. Each tree
        will receive the corresponding weight. If omitted, all trees will be
        equally weighted.
    cutoff : float, 0.0 <= cutoff <= 1.0
        Any clade that has <= cutoff support will be dropped. If cutoff is
        < 0.5, then it is possible that ties will result. If so, ties are
        broken arbitrarily depending on list sort order.
    support_attr : str
        The attribute to be decorated onto the resulting trees that contain the
        consensus support.

    Returns
    -------
    list of TreeNode
        Multiple trees can be returned in the case of two or more disjoint sets
        of tips represented on input.

    Notes
    -----
    This code was adapted from PyCogent's majority consensus code originally
    written by Matthew Wakefield. The method is based off the original
    description of consensus trees in [1]_. An additional description can be
    found in the Phylip manual [2]_. This method does not support majority rule
    extended.

    Support is computed as a weighted average of the tree weights in which the
    clade was observed in. For instance, if {A, B, C} was observed in 5 trees
    all with a weight of 1, its support would then be 5.

    References
    ----------
    .. [1] Margush T, McMorris FR. (1981) "Consensus n-trees." Bulletin for
           Mathematical Biology 43(2) 239-44.
    .. [2] http://evolution.genetics.washington.edu/phylip/doc/consense.html

    Examples
    --------
    Computing the majority consensus, using the example from the Phylip manual
    with the exception that we are computing majority rule and not majority
    rule extended.

    >>> from skbio.tree import TreeNode
    >>> trees = [
    ... TreeNode.from_newick("(A,(B,(H,(D,(J,(((G,E),(F,I)),C))))));"),
    ... TreeNode.from_newick("(A,(B,(D,((J,H),(((G,E),(F,I)),C)))));"),
    ... TreeNode.from_newick("(A,(B,(D,(H,(J,(((G,E),(F,I)),C))))));"),
    ... TreeNode.from_newick("(A,(B,(E,(G,((F,I),((J,(H,D)),C))))));"),
    ... TreeNode.from_newick("(A,(B,(E,(G,((F,I),(((J,H),D),C))))));"),
    ... TreeNode.from_newick("(A,(B,(E,((F,I),(G,((J,(H,D)),C))))));"),
    ... TreeNode.from_newick("(A,(B,(E,((F,I),(G,(((J,H),D),C))))));"),
    ... TreeNode.from_newick("(A,(B,(E,((G,(F,I)),((J,(H,D)),C)))));"),
    ... TreeNode.from_newick("(A,(B,(E,((G,(F,I)),(((J,H),D),C)))));")]
    >>> consensus = majority_rule(trees, cutoff=0.5)[0]
    >>> print(consensus.ascii_art())
                                  /-E
                                 |
                                 |          /-G
                        /--------|         |
                       |         |         |          /-F
                       |         |         |---------|
                       |          \--------|          \-I
                       |                   |
                       |                   |          /-C
              /--------|                   |         |
             |         |                    \--------|          /-D
             |         |                             |         |
             |         |                              \--------|--J
    ---------|         |                                       |
             |         |                                        \-H
             |         |
             |          \-B
             |
              \-A
    >>> for node in consensus.non_tips():
    ...     support_value = node.support
    ...     names = ' '.join([n.name for n in node.tips()])
    ...     print("Tips: %s, support: %s" % (names, support_value))
    Tips: F I, support: 9.0
    Tips: D J H, support: 6.0
    Tips: C D J H, support: 6.0
    Tips: G F I C D J H, support: 6.0
    Tips: E G F I C D J H, support: 9.0
    Tips: E G F I C D J H B, support: 9.0

    In the next example, multiple trees will be returned which can happen if
    clades are not well supported across the trees. In addition, this can arise
    if not all tips are present across all trees.

    >>> trees = [
    ...     TreeNode.from_newick("((a,b),(c,d),(e,f))"),
    ...     TreeNode.from_newick("(a,(c,d),b,(e,f))"),
    ...     TreeNode.from_newick("((c,d),(e,f),b)"),
    ...     TreeNode.from_newick("(a,(c,d),(e,f))")]
    >>> consensus_trees = majority_rule(trees)
    >>> print(len(consensus_trees))
    4
    >>> for tree in consensus_trees:
    ...     print(tree.ascii_art())
    --b
    --a
              /-f
    ---------|
              \-e
              /-d
    ---------|
              \-c

    """
    if weights is None:
        weights = np.ones(len(trees), dtype=float)
    else:
        weights = np.asarray(weights)
        if len(weights) != len(trees):
            raise ValueError("Number of weights and trees differ!")

    cutoff_threshold = cutoff * weights.sum()

    clade_counts, edge_lengths = _walk_clades(trees, weights)
    clade_counts = _filter_clades(clade_counts, cutoff_threshold)
    trees = _build_trees(clade_counts, edge_lengths, support_attr)

    return trees
