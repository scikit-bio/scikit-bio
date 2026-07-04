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

import skbio
import numpy as np
cimport numpy as cnp

from ._bp cimport BPTree


def to_skbio_treenode(BPTree bp):
    """Convert BPTree to TreeNode

    Parameters
    ----------
    bp : BPTree
        A BPTree

    Returns
    -------
    skbio.TreeNode
        The tree represented as an skbio.TreeNode
    """
    cdef int i

    nodes = [skbio.TreeNode() for i in range(bp.data.sum())]

    root = nodes[0]

    for i in range(bp.data.sum()):
        node_idx = bp.preorder_select(i)
        nodes[i].name = bp.name(node_idx)
        nodes[i].length = bp.length(node_idx)
        nodes[i].edge_num = bp.edge(node_idx)

        if node_idx != bp.root():
            # preorder starts at 1 annoyingly
            parent = bp.preorder_rank(bp.parent(node_idx)) - 1
            # uncache=False avoids repeated cache clearing during bulk construction
            nodes[parent].append(nodes[i], uncache=False)

    root.length = None

    return root


def from_skbio_treenode(tree):
    """Convert a skbio TreeNode into BPTree

    Parameters
    ----------
    tree : skbio.TreeNode
        The tree to convert

    Returns
    -------
    tuple
        (BPTree, np.array of str, np.array of double)
    """
    n_nodes = len(list(tree.traverse(include_self=True)))

    topo = np.zeros(n_nodes * 2, dtype=np.uint8)
    names = np.full(n_nodes * 2, None, dtype=object)
    lengths = np.zeros(n_nodes * 2, dtype=np.double)
    edges = np.zeros(n_nodes * 2, dtype=np.int32)

    ptr = 0
    seen = set()
    for n in tree.pre_and_postorder(include_self=True):
        if n not in seen:
            topo[ptr] = 1
            names[ptr] = n.name
            lengths[ptr] = n.length or 0.0
            edges[ptr] = getattr(n, 'edge_num', None) or 0

            if n.is_tip():
                ptr += 1

            seen.add(n)

        ptr += 1
    return BPTree(topo, names=names, lengths=lengths, edges=edges)


def to_skbio_treearray(BPTree bp):
    """Convert to a tree array comparable to TreeNode.to_array

    Parameters
    ----------
    bp : BPTree
        A BPTree

    Returns
    -------
    dict
        Dict with keys 'child_index', 'length', 'id_index', 'name'
    """
    cdef int i

    class mock_node:
        def __init__(self, id, is_tip):
            self.is_tip_ = is_tip
            self.id = id

        def is_tip(self):
            return self.is_tip_

    child_index = np.zeros((int(bp.data.sum()) - bp.ntips(), 3), dtype=np.int64)
    length = np.zeros(bp.data.sum(), dtype=np.double)
    node_ids = np.zeros(bp.data.size, dtype=np.uint32)
    name = np.full(bp.data.sum(), None, dtype=object)

    # TreeNode.assign_ids, decompose target
    chi_ptr = 0
    cur_index = 0  # the index into node_ids, equivalent to TreeNode.assign_ids
    id_index = dict.fromkeys(set(range(bp.data.sum())))  # map a node's "id" to an object which indicates if it is a leaf or not
    for i in range(bp.data.sum()):
        node_idx = bp.postorder_select(i + 1)  # the index within the BP of the node

        if not bp.is_tip(node_idx):
            first_child = bp.first_child(node_idx)
            last_child = bp.last_child(node_idx)

            sib_idx = first_child  # the sibling index wtihin the BP of the node
            while sib_idx != 0 and sib_idx <= last_child:
                node_ids[sib_idx] = cur_index
                id_index[cur_index] = mock_node(cur_index, bp.is_tip(sib_idx))
                length[cur_index] = bp.length(sib_idx)
                name[cur_index] = bp.name(sib_idx)

                cur_index += 1
                sib_idx = bp.next_sibling(sib_idx)

            child_index[chi_ptr] = [node_idx, node_ids[first_child], node_ids[last_child]]
            chi_ptr += 1

    # make sure to capture root
    id_index[bp.data.sum() - 1] = mock_node(cur_index, False)

    node_ids[0] = cur_index
    child_index[:, 0] = node_ids[child_index[:, 0]]
    child_index = child_index[np.argsort(child_index[:, 0])]

    return {'child_index': child_index, 'length': length, 'id_index': id_index,
            'name': name}
