# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
from skbio.tree import TreeNode
from skbio.tree._cutils import num_dist_cy


def _ordered(tree_array):
    """Constructs a postorder of nodes as a numerical array.

    First row contains node locations, and second row contains the id indexes of
    leaf nodes and zero for internal nodes.

    TODO: Pre-allocate ordered list since you know its size.

    """
    # take row of node locations of taxa
    tree_array_top = tree_array[0]
    # remove unassigned taxa
    filtered_tree_array = tree_array[:, tree_array[0] != 0]
    # initialize node locations list for postorder
    ordered = []
    # iterating over leaf nodes
    for x in filtered_tree_array[0]:
        if x % 2 != 0:
            ordered.append(int(x))
        else:
            ordered.append(int(x))
            y = int((x - 2) * 0.5)
            ordered.append(int(y))
            while y % 2 == 0 and y != 0:
                y = int((y - 2) * 0.5)
                ordered.append(int(y))
    ids_index = np.zeros_like(ordered)
    value_map = dict(zip(filtered_tree_array[0], filtered_tree_array[1]))
    for i, b_val in enumerate(ordered):
        if b_val in value_map:
            ids_index[i] = value_map[b_val]
    return np.array([ordered, ids_index])


def _pair_lca(left_node, right_node):
    """Finds the LCA of a pair of nodes."""
    # initialize variables the left node's path back to root
    left_path = []
    # iterating over parents of left node
    # until back to root which is represented as zero
    while left_node != 0:
        if left_node % 2 == 0:
            left_node = int(((left_node - 2) * 0.5))
        else:
            left_node = int(((left_node - 1) * 0.5))
        # left_path contains nodes back to parent in order
        left_path.append(left_node)
    # iterating over parents of right night
    # until reaching a node in the left path which
    # represents the LCA
    while right_node not in left_path:
        if right_node % 2 == 0:
            right_node = int(((right_node - 2) * 0.5))
        else:
            right_node = int(((right_node - 1) * 0.5))
    return right_node


def _subtree_root(taxa):
    """Find the root node of a subtree given the leaves of the subtree."""
    # find LCA of first pair
    lca = _pair_lca(taxa[0], taxa[1])
    # iterate over remaining nodes in list of taxa.
    for taxon in taxa[2 : len(taxa)]:
        # skip taxa that are ancestral to
        # the current lca
        if num_dist_cy(lca, taxon) >= 0:
            continue
        lca = _pair_lca(lca, taxon)
    return lca


def _parent(node):
    """Find parent node of a given node."""
    # case for node with even value location
    if node % 2 == 0:
        parent = int((node - 2) * 0.5)
    # case for node with odd value location
    else:
        parent = int((node - 1) * 0.5)
    return parent


def _sibling(node):
    """Find sibling node of a given node."""
    # case for node with even value location
    if node % 2 == 0:
        sibling = node - 1
    # case for node with odd value location
    else:
        sibling = node + 1
    return sibling


def _ancestors(x, array):
    """Returns list of ancestral nodes for
    a given node.

    """
    # _num_dist checks for ancestral lineage
    result = [v for v in array if num_dist_cy(v, x) > 0]
    return np.array(result)


def _subtree(x, array):
    """Return nodes in array of node locations that
    are in the subtree rooted at a given node.

    """
    result = [v for v in array if num_dist_cy(x, v) >= 0]
    return np.array(result)


def _move_subtree(tree_array, subtree_nodes, old_lca, new_lca):
    """Adjusts the node locations of nodes in a subtree where the root
    is relocated to a new node location value.

    """
    # moves each node in subtree individually
    for leaf in subtree_nodes:
        tree_array[0, np.where(tree_array[0] == leaf)[0][0]] = _move_node(
            leaf, old_lca, new_lca
        )


def _move_node(x, old_lca, new_lca):
    """Adjusts the node location of a given node when also given an initial and final
    location not necessarily the location of the given node. The initial node
    location can be viewed as the location of the root of a subtree that contains
    the given node.

    For use when not considering a subtree, the initial node location can be the
    location of the given node.
    """
    diff = new_lca - old_lca
    # order of arguments in distance function ensures a positive value is returned
    return x + (diff * (2 ** num_dist_cy(old_lca, x)))


def _array_to_tree(taxa, tree_array):
    # initialize variables
    internal = []
    iterating = True
    nodes = {}

    # create leaf_array from reordering tree_array from index values
    sorted_tree_indices = np.argsort(tree_array[1])
    sorted_tree_array = tree_array[:, sorted_tree_indices]
    leaf_array = sorted_tree_array[0][sorted_tree_array[0] != 0]

    # Iterate over nodes in order of greatest value, where values of
    # nodes are according to their index in a complete binary tree.
    # Idea is to have the nodes index in a complete binary tree as a
    # value to avoid needing unnecessary spaces in the array.
    # Also, only leaf node values are needed which further reduces the
    # size of the array to n-1 for a tree with n leaves.
    while iterating:
        node = np.max(np.concatenate((leaf_array, internal)))
        # Create nodes as they are chosen from the array.
        # Can be replaced with Newick strings instead of TreeNode.
        # Set value to 0 to show the node has been chosen and created.
        if node in leaf_array:
            node_index = np.where(leaf_array == node)[0][0]
            node_a = TreeNode.read([taxa[node_index + 1] + ";"])
            leaf_array[np.where(leaf_array == node)[0][0]] = 0
        # Set node to already created internal node, and
        # then set value to 0 to show the node has been chosen.
        else:
            node_a = nodes["node_" + str(node)]
            internal[internal.index(node)] = 0
        # Check if even or odd, then find sibling and parent nodes.
        if node % 2 == 0:
            sibling = node - 1
            parent = int((node - 2) * 0.5)
        else:
            sibling = node + 1
            parent = int((node - 1) * 0.5)
        # Create node for sibling if needed, and updated arrays similar to above
        if sibling in leaf_array:
            sibling_index = np.where(leaf_array == sibling)[0][0]
            node_b = TreeNode.read([taxa[sibling_index + 1] + ";"])
            leaf_array[sibling_index] = 0
        else:
            node_b = nodes["node_" + str(sibling)]
            internal[np.where(internal == sibling)[0][0]] = 0
        # Attach nodes to parent.
        # Can be replaced with code using Newick strings.
        nodes[f"node_{parent}"] = TreeNode.read([";"])
        nodes[f"node_{parent}"].append(node_b)
        nodes[f"node_{parent}"].append(node_a)
        internal.append(parent)
        # Stop iterating when all nodes have been chosen
        if np.max(np.concatenate((leaf_array, internal))) == 0:
            iterating = False
    # Create root node and append the rest of the tree.
    tree = TreeNode.read([taxa[0] + ";"])
    tree.append(nodes["node_0"])
    return tree
