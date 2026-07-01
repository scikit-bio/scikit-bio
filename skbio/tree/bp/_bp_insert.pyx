# encoding: utf-8
# cython: profile=False, boundscheck=False, wraparound=False
# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Derived from improved-octo-waddle (https://github.com/biocore/improved-octo-waddle)
# originally authored by Daniel McDonald, distributed under the Modified BSD License.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from ._bp cimport BPTree
from ._bp_conv import to_skbio_treenode
import pandas as pd
import json
import skbio
cimport cython


def _preallocate_fragment(r):
    return skbio.TreeNode(name=r['fragment'], length=r['pendant_length'])


def _preallocate_empty(r):
    return skbio.TreeNode()


def _insert_setup(placements, bptree, insert_type):
    # insertion setup addresses:
    # * treenode caching
    # * placement ordering
    # * preallocation of objects where "easy"

    sktree = to_skbio_treenode(bptree)
    node_lookup = {n.edge_num: n for n in sktree.traverse(include_self=True)}

    # we are only setup to handle a single placement per fragment, so pull
    # deduplicated following guidance from Prof. Siavash Mirarab. We sort so
    # "better" has a smaller index value
    # fragment -> group the rows by the fragment, fragment order doesn't matter
    # like_weight_ratio -> our first selection criteria, higher is better
    # pendant_length -> our second selection criteria, lower is better
    placements = placements.sort_values(['fragment', 'like_weight_ratio',
                                         'pendant_length'],
                                        ascending=[True, False, True])

    # take the first non-duplicated row per fragment. because of the sort, this
    # is assured to be the highest weight ratio, and the smallest pendant
    # length. Ties are handled arbitrarily.
    placements = placements[~placements['fragment'].duplicated()]

    if insert_type == 'multifurcating':
        placements = placements.sort_values(['edge_num', 'pendant_length'])
    elif insert_type == 'fully_resolved':
        placements = placements.sort_values(['edge_num', 'distal_length'],
                                            ascending=[True, False])
    else:
        raise ValueError()

    # Use .assign() rather than ``placements['node'] = ...`` to attach columns.
    # pandas 3.0's chained-assignment detector flags ``df[key] = value`` when
    # ``sys.getrefcount(df)`` is low and ``df`` is not found as a local in the
    # *caller's* Python frame. That frame introspection cannot see Cython
    # locals, so direct column assignment from this compiled module always
    # raises a spurious ChainedAssignmentError. .assign() performs the
    # assignment inside pandas' own (pure-Python) frame, avoiding the misfire.
    placements = placements.assign(
        node=placements.apply(_preallocate_fragment, axis=1))

    if insert_type == 'fully_resolved':
        placements = placements.assign(
            parent=placements.apply(_preallocate_empty, axis=1))

    return placements, sktree, node_lookup


# pd.DataFrame is not a resolved type so we cannot use it here for cython
def insert_fully_resolved(object placements, BPTree bptree):
    """Update the backbone, fully resolving edges with multiple queries

    Parameters
    ----------
    placements : pd.DataFrame
        jplace data represented as a DataFrame
    bptree : bp.BPTree
        An instance of a BPTree, this is expected to contain edge numbers
        and correspond to the backbone for the jplace data

    Returns
    -------
    skbio.TreeNode
        A tree with the fragments placed
    """
    placements, sktree, node_lookup = \
        _insert_setup(placements, bptree, 'fully_resolved')

    for edge, edge_grp in placements.groupby('edge_num'):
        existing_node = node_lookup[edge]
        current_parent = existing_node.parent

        # break the edge
        # uncache=False: no caches are present during bulk construction
        current_parent.remove(existing_node, uncache=False)
        existing_node.parent = None
        existing_length = existing_node.length

        for _, fragment in edge_grp.iterrows():
            distal_length = fragment['distal_length']
            fragment_node = fragment['node']
            fragment_parent = fragment['parent']

            # update branch lengths
            fragment_parent.length = existing_length - distal_length
            existing_length = distal_length

            # attach the nodes
            fragment_parent.append(fragment_node, uncache=False)
            current_parent.append(fragment_parent, uncache=False)

            # update
            current_parent = fragment_parent

        existing_node.length = existing_length
        current_parent.append(existing_node, uncache=False)
        existing_node.length = distal_length

    return sktree
