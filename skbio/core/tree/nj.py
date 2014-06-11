from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

r"""

>>> from skbio.core.distance import DistanceMatrix
>>> from skbio.core.tree.nj import nj

>>> data = [[0,  5,  9,  9,  8],
            [5,  0, 10, 10,  9],
            [9, 10,  0,  8,  7],
            [9, 10,  8,  0,  3],
            [8,  9,  7,  3,  0]]
>>> ids = list('abcde')
>>> dm = DistanceMatrix(data, ids)

>>> tree = nj(dm)
>>> print tree.ascii_art()

          /-d
         |
         |          /-c
         |---------|
---------|         |          /-b
         |          \--------|
         |                    \-a
         |
          \-e

>>> newick_str = nj(dm, result_constructor=str)
>>> print newick_str
(d:2, (c:4, (b:3, a:2):3):2, e:1);

"""

import numpy as np

from skbio.core.distance import DistanceMatrix
from skbio.core.tree import TreeNode

def nj(dm, dissallow_negative_branch_length=True,
       result_constructor=TreeNode.from_newick):
    while(dm.shape[0] > 3):
        q = _compute_q(dm)
        idx1, idx2 = _lowest_index(q)
        pair_member_1 = dm.ids[idx1]
        pair_member_2 = dm.ids[idx2]
        pair_member_1_len, pair_member_2_len = _pair_members_to_new_node(dm,
                                                                        idx1,
                                                                        idx2,
                                                                        dissallow_negative_branch_length)
        node_definition = "(%s:%d, %s:%d)" % (pair_member_1,
                                              pair_member_1_len,
                                              pair_member_2,
                                              pair_member_2_len)
        dm = _compute_collapsed_dm(dm, pair_member_1, pair_member_2,
                                  dissallow_negative_branch_length=dissallow_negative_branch_length,
                                  new_node_id=node_definition)

    # Define the last node - this is an internal node
    pair_member_1 = dm.ids[1]
    pair_member_2 = dm.ids[2]
    pair_member_1_len, pair_member_2_len = _pair_members_to_new_node(dm,
                                                                    pair_member_1,
                                                                    pair_member_2,
                                                                    dissallow_negative_branch_length)
    internal_len = _otu_to_new_node(dm, pair_member_1, pair_member_2, node_definition,
                                   dissallow_negative_branch_length=dissallow_negative_branch_length,
                                   )
    newick = "(%s:%d, %s:%d, %s:%d);" % (pair_member_1, pair_member_1_len,
                                         node_definition, internal_len,
                                         pair_member_2, pair_member_2_len)

    return result_constructor(newick)

def _compute_q(dm):
    q = np.zeros(dm.shape)
    n = dm.shape[0]
    for i in range(n):
        for j in range(i):
            q[i, j] = q[j, i] = ((n - 2) * dm[i, j]) - dm[i].sum() - dm[j].sum()
    return DistanceMatrix(q, dm.ids)

def _compute_collapsed_dm(dm, i, j, dissallow_negative_branch_length, new_node_id=None):
    in_n = dm.shape[0]
    out_n = in_n - 1
    new_node_id = new_node_id or "(%s, %s)" % (i, j)
    out_ids = [new_node_id]
    out_ids.extend([e for e in dm.ids if e not in (i, j)])
    result = np.zeros((out_n, out_n))
    for idx1, out_id1 in enumerate(out_ids[1:]):
        result[0, idx1 + 1] = result[idx1 + 1, 0] = \
         _otu_to_new_node(dm, i, j, out_id1, dissallow_negative_branch_length)
        for idx2, out_id2 in enumerate(out_ids[1:idx1+1]):
            result[idx1+1, idx2+1] = result[idx2+1, idx1+1] = dm[out_id1, out_id2]
    return DistanceMatrix(result, out_ids)

def _lowest_index(dm):
    lowest_value = np.inf
    for i in range(dm.shape[0]):
        for j in range(i):
            curr_index = i, j
            curr_value = dm[curr_index]
            if curr_value < lowest_value:
                lowest_value = curr_value
                result = curr_index
    return result

def _otu_to_new_node(dm, i, j, k, dissallow_negative_branch_length):
    k_to_u = 0.5 * (dm[i, k] + dm[j, k] - dm[i, j])

    if dissallow_negative_branch_length and k_to_u < 0:
        k_to_u = 0

    return k_to_u

def _pair_members_to_new_node(dm, i, j, dissallow_negative_branch_length):
    n = dm.shape[0]
    i_to_j = dm[i, j]
    i_to_u = (0.5 * i_to_j) + (1 / (2 * (n - 2))) * (dm[i].sum() - dm[j].sum())
    j_to_u = i_to_j - i_to_u

    if dissallow_negative_branch_length and i_to_u < 0:
        i_to_u = 0
    if dissallow_negative_branch_length and j_to_u < 0:
        j_to_u = 0

    return i_to_u, j_to_u
