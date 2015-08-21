# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np


def traverse_reduce(bound_indices, f):
    """Applies a[i] = f(a[j:k]) over list of [(a[i], a[j:k])].

    If list is in traversal order, has same effect as consolidating the
    function over the tree, only much faster.

    Note that f(a[j:k]) must return an object that can be broadcast to the
    same shape as a[i], e.g. summing a 2D array to get a vector.
    """
    for i, s in bound_indices:
        i[:] = f(s, 0)


def _skbio_counts_to_envs(otu_ids, *counts):
    """This is a place holder method to help get the API hooked up

    This method should be unnecessary once the implementation is optimized to
    the interface.
    """
    envs = {}
    n_counts = len(counts)

    for packed in zip(otu_ids, *counts):
        # NOTE: this is ducttape to fit the API
        otu_id = packed[0]

        counts = {}
        for env in range(1, n_counts + 1):
            if packed[env]:
                counts[env] = packed[env]

        if counts:
            envs[otu_id] = counts

    return envs


def bool_descendants(bound_indices):
    """For each internal node, sets col to True if any descendant is True."""
    traverse_reduce(bound_indices, np.logical_or.reduce)


def bind_to_array(tree_index, a):
    """Binds tree_index to array a, returning result in list.

    Takes as input list of (node, first_child, last_child)
    returns list of (node_row, child_rows) such that node_row points to the
    row of a that corresponds to the current node, and child_rows points to the
    row or rows of a that correspond to the direct children of the current
    node.

    Order is assumed to be traversal order, i.e. for the typical case of
    postorder traversal iterating over the items in the result and
    consolidating each time should give the same result as postorder
    traversal of the original tree. Should also be able to modify for
    preorder traversal.
    """
    # note: range ends with end+1, not end, b/c end is included
    return [(a[node], a[start:end+1]) for node, start, end in tree_index]


def _fast_unifrac_setup(t, envs, make_subtree=True):
    """Setup shared by fast_unifrac and by significance tests."""
    if make_subtree:
        t2 = t.copy()
        wanted = set(envs.keys())

        def delete_test(node):
            if node.is_tip() and node.name not in wanted:
                return True
            return False
        t2.remove_deleted(delete_test)
        t2.prune()
        t = t2

    # index tree
    node_index, nodes = index_tree(t)
    # get good nodes, defined as those that are in the env file.
    good_nodes = dict([(i.name, envs[i.name])
                       for i in t.tips() if i.name in envs])
    envs = good_nodes
    (count_array, unique_envs,
     env_to_index, node_to_index) = index_envs(envs, node_index)
    env_names = sorted(unique_envs)
    # Note: envs get sorted at the step above
    branch_lengths = get_branch_lengths(node_index)
    if not envs:
        raise ValueError(''.join(["No valid samples/environments found.",
                                  "Check whether tree tips match ",
                                  "otus/taxa present ",
                                  "in samples/environments"]))
    return (envs, count_array, unique_envs, env_to_index,
            node_to_index, env_names, branch_lengths, nodes, t)


def index_tree(t):
    """

    Indexes nodes in-place as n._leaf_index.

    Algorithm is as follows:
    for each node in post-order traversal over tree:
        if the node has children:
            set an index on each child
            for each child with children:
                add the child and its start and end tips to the result

    Parameters
    ----------
    t : skbio.TreeNode
        Phylogenetic tree

    Returns
    -------
    tuple
        Contains {node_id:node}, [node_id,first_child,last_child]
    """
    id_index = {}    # needs to be dict, not list, b/c adding out of order
    child_index = []
    curr_index = 0
    for n in t.traverse(self_before=False, self_after=True):
        for c in n.children:
            c._leaf_index = curr_index
            id_index[curr_index] = c
            curr_index += 1
            if c:    # c has children itself, so need to add to result
                child_index.append((c._leaf_index, c.children[0]._leaf_index,
                                    c.children[-1]._leaf_index))
    # handle root, which should be t itself
    t._leaf_index = curr_index
    id_index[curr_index] = t
    # only want to add to the child_index if t has children...
    if t.children:
        child_index.append((t._leaf_index, t.children[0]._leaf_index,
                            t.children[-1]._leaf_index))
    return id_index, child_index


def index_envs(env_counts, tree_index, array_constructor=int):
    """
    Calculates the taxon counts in each env (what is an env?)

    Parameters
    ----------
    env_counts : np.array
       Should be the output of count_envs(lines).
    tree_index : np.array
       Should be the id_index of index_tree(t).
    array_constructor: function
        int by default (may need to change to float later
        to handle microarray data).

    Returns
    -------
    np.array
        Array of taxon x env with counts of the taxon in each env.
    """
    num_nodes = len(tree_index)
    unique_envs, num_envs = get_unique_envs(env_counts)
    env_to_index = dict([(e, i) for i, e in enumerate(unique_envs)])
    result = np.zeros((num_nodes, num_envs), array_constructor)
    # figure out taxon label to index map
    node_to_index = {}
    for i, node in tree_index.items():
        if node.name is not None:
            node_to_index[node.name] = i
    # walk over env_counts, adding correct slots in array
    for name in env_counts:
        curr_row_index = node_to_index[name]
        for env, count in env_counts[name].items():
            result[curr_row_index, env_to_index[env]] = count
    # return all the data structures we created; will be useful for other tasks
    return result, unique_envs, env_to_index, node_to_index


def get_branch_lengths(tree_index):
    """
    Parameters
    ----------
    tree_index: dict


    Returns
    -------
    np.array
        Array of branch lengths, in tree index order.
    """
    result = np.zeros(len(tree_index), float)
    for i, node in tree_index.items():
        try:
            if node.length is not None:
                result[i] = node.length
        except AttributeError:
            pass
    return result


def get_unique_envs(envs):
    """extract all unique envs from envs dict"""
    result = set()
    for v in envs.values():
        result.update(v.keys())
    # sort envs for convenience in testing and display
    return sorted(result), len(result)
