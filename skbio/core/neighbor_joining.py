#!/usr/bin/env python
""" Neighbour joining phylogenetic tree estimation.

Note that by default,
negative branch lengths are reset to 0.0 during the calculations.

"""
from __future__ import division
import numpy
from skbio.core.tree import TreeNode


# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


def rnj(distmtx, no_negatives=True):
    """Computes a tree using the neighbor joining method
    
    Arguments:
        - dists: dict of (name1, name2): distance
        - no_negatives: negative branch lengths will be set to 0
    """

    d, names = distmtx.data.copy(), distmtx.ids 
    # need copy for safety -  we modify d in place (anc copy) below
    
    nodes = [TreeNode(name=name) for name in names]
    while len(nodes) > 2:
        # total distance from other tips / net divergence / r
        r = d.sum(0)

        # rate_corrected_d is used only to find neighbors i, j
        rate_corrected_d = d.copy()
        for i in range(len(nodes)):
            for j in range(len(nodes)):
                rate_corrected_d[i,j] -= (r[i] + r[j]) / (len(nodes)-2)
        diag_num = rate_corrected_d.max() + 1
        for i in range(len(nodes)):
            rate_corrected_d[i,i] = diag_num
            # will soon need to find minimum off diagonal
        i, j = numpy.unravel_index(rate_corrected_d.argmin(),
                                    rate_corrected_d.shape)

        num_nodes = len(nodes)
        n1 = nodes[i]
        n1.length = .5 * (d[i,j] + (r[i] - r[j])/(num_nodes-2))
        n2 = nodes[j]
        n2.length = .5 * (d[i,j] + (r[j] - r[i])/(num_nodes-2))


        # no negative branch lengths
        if no_negatives:
            nodes[i].length = max(0.0, nodes[i].length)
            nodes[j].length = max(0.0, nodes[j].length)
        

        new_node = TreeNode(children=[nodes[i], nodes[j]])
        
        # Store new_node at i
        new_dists = 0.5 * (d[i] + d[j] - d[i,j])
        d[:, i] = new_dists
        d[i, :] = new_dists
        d[i, i] = 0.0
        nodes[i] = new_node
        
        # rm j
        d[j, :] = d[num_nodes-1, :]
        d[:, j] = d[:, num_nodes-1]
        assert d[j, j] == 0.0, d
        d = d[0:num_nodes-1, 0:num_nodes-1]
        nodes[j] = nodes[num_nodes-1]
        nodes.pop()


    # join last two nodes and root halfway between the two nodes
    final_dist =  .5 * d[0,1]
    nodes[0].length = final_dist
    nodes[1].length = final_dist

    # arbitrary placement of root, nj doesn't know 
    # anything about rooting sans outgroups
    root_node = TreeNode(children=nodes)

    return root_node
