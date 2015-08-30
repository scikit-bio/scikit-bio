# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np

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


def _counts_and_length(counts, otu_ids, tree, indexed):
    """Get the counts array and the tree branch lengths"""
    if indexed is None:
        indexed = index_tree(tree)

    counts = np.atleast_2d(counts)
    count_array = nodes_by_counts(counts, otu_ids, indexed)

    return count_array, indexed
