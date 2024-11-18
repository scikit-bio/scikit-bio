# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
from skbio.tree import TreeNode


def _ordered(tree_array):
    tree_array_top = tree_array[0]
    filtered_tree_array = tree_array[:, tree_array[0] != 0]
    ordered = []
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


def pair_lca(left_node, right_node):
    left_path = []
    while left_node != 0:
        if left_node % 2 == 0:
            left_node = int(((left_node - 2) * 0.5))
        else:
            left_node = int(((left_node - 1) * 0.5))
        left_path.append(left_node)
    while right_node not in left_path:
        if right_node % 2 == 0:
            right_node = int(((right_node - 2) * 0.5))
        else:
            right_node = int(((right_node - 1) * 0.5))
    return right_node


def num_dist(anc, desc):
    dist = 0
    while desc != anc:
        if desc < anc:
            dist = -1
            break
        if desc % 2 == 0:
            desc = int(((desc - 2) * 0.5))
        else:
            desc = int(((desc - 1) * 0.5))
        dist += 1
    # if dist == -1:
    # raise ValueError("no ancestral path")
    return dist


def subtree_root(taxa):
    lca = pair_lca(taxa[0], taxa[1])
    for taxon in taxa[2 : len(taxa)]:
        lca = pair_lca(lca, taxon)
    return lca


def _parent(node):
    if node % 2 == 0:
        parent = int((node - 2) * 0.5)
    else:
        parent = int((node - 1) * 0.5)
    return parent


def _sibling(node):
    if node % 2 == 0:
        sibling = node - 1
    else:
        sibling = node + 1
    return sibling


def _ancestors(x, array):
    result = [v for v in array if num_dist(v, x) > 0]
    return np.array(result)


def _subtree(x, array):
    result = [v for v in array if num_dist(x, v) >= 0]
    return np.array(result)


def move_subtree(tree_array, subtree_nodes, old_lca, new_lca):
    for leaf in subtree_nodes:
        tree_array[0, np.where(tree_array[0] == leaf)[0][0]] = move_node(
            leaf, old_lca, new_lca
        )
    return


def move_node(x, old_lca, new_lca):
    diff = new_lca - old_lca
    return x + (diff * (2 ** num_dist(old_lca, x)))


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
    if iterating:
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
