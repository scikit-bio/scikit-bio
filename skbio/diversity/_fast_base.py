# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from warnings import warn

import numpy as np

from skbio.util import EfficiencyWarning
from skbio.diversity._fast_base_cy import nodes_by_counts


def index_tree(tree):
    """Index a tree to allow for bulk numpy aggregations

    Paramters
    ---------
    tree : TreeNode
        A tree to index

    Returns
    -------
        dict of array
            {id_index: {id: TreeNode},
             child_index: [(node_id, left_child_id, right_child_id)],
             attr_1: array(...),
             ...
             attr_N: array(...)}

    Notes
    -----
    This wraps `TreeNode.to_array`, but replaces any `nan` in the length array
    with 0.0 which can arise from an edge not having a length, notably the
    root node parent edge.
    """
    indexed = tree.to_array()
    length = indexed['length']
    indexed['length'] = np.where(np.isnan(length), 0.0, length)

    return indexed


### preserving for the purposes of comments, this method must be deleted
### prior to merge
def nodes_by_counts_py(counts, tip_ids, indexed_tree):
    """Construct the count array, and the counts up the tree

    Parameters
    ----------
    counts : np.array of int
        A 1D or 2D vector in which each row corresponds to the observed counts
        in an environment. The rows are expected to be in order with respect to
        `tip_ids`.
    tip_ids : np.array of str
        A vector of tip names that correspond to the columns in the `counts`
        matrix.
    indexed_tree : dict
        The result of `index_tree`.

    Returns
    -------
    np.array of int
        The observed counts of every node and the counts if its descendents.

    Notes
    -----

    """
    nodes = indexed_tree['name']

    # allow counts to be a vector
    counts = np.atleast_2d(counts)

    # determine observed IDs
    observed_indices = np.where(np.logical_or.reduce(counts))[0]
    observed_ids = tip_ids[observed_indices]
    observed_ids_set = set(observed_ids)

    # construct mappings of the observed to their positions in the node array
    node_lookup = {n: i for i, n in enumerate(nodes) if n in observed_ids_set}
    otus_in_nodes = np.array([node_lookup[n] for n in observed_ids])

    # count_array has a column per node (not tip) and a row per env. This is
    # transposed on return as the subsequent operations will collapse over
    # the nodes.
    n_count_vectors = counts.shape[0]
    count_array = np.zeros((n_count_vectors, len(nodes)), dtype=int)
    count_array[:, otus_in_nodes] = counts[:, observed_indices]
    count_array = count_array.T

    ### method call no longer possible as it is now a cdef
    # traverse_reduce(indexed_tree['child_index'], count_array)

    return count_array


def _counts_and_length(counts, otu_ids, tree, indexed):
    """Get the counts array and the tree branch lengths"""
    if indexed is None:
        warn("It is recommended to drive the fast unifrac method with an "
             "indexed tree. Please see skbio.diversity.index_tree.",
             EfficiencyWarning)

        indexed = index_tree(tree)

    counts = np.atleast_2d(counts)
    count_array = nodes_by_counts(counts, otu_ids, indexed)

    return count_array, indexed['length']
