#!/usr/bin/env python
"""Relaxed Neighbour joining phylogenetic tree estimation.

Note that by default,
negative branch lengths are reset to 0.0 during the calculations.

This code is primarily based off of Gavin Huttley's nj.py code version 1.1,
however, the algorithm does only a local search for neighbors to join, thus
reducing (theoretically at least) computational time to typically O(N^2 log N)

See, for example:
Relaxed Neighbor Joining: A Fast Distance-Based Phylogenetic Tree
Construction Method, by Jason Evans, Luke Sheneman, James Foster

If this algorithm is the bottleneck of an expensive (in computation time) task,
it may be worthwhile to use clearcut, which is written in c.

"""
from __future__ import division
import numpy
from random import shuffle, seed
from skbio.core.tree import TreeNode
# from cogent.core.tree import TreeBuilder
# from cogent.phylo.util import distanceDictTo2D

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------


def rnj(distmtx, no_negatives=True, randomize=True):
    """Computes a tree using the relaxed neighbor joining method
    
    Arguments:
        - dists: dict of (name1, name2): distance
        - no_negatives: negative branch lengths will be set to 0
        - randomize: the algorithm will search nodes randomly until two
        neighbors are found.
    """
        
    # constructor = TreeBuilder(mutable=True).createEdge
    # (names, d) = distanceDictTo2D(dists)
    d, names = distmtx.data.copy(), distmtx.ids 
    # need copy because we modify d in place

    # nodes = [constructor([], name, {}) for name in names]
    
    nodes = [TreeNode(name=name) for name in names]
    # print nodes
    while len(nodes) > 2:
        r = d.sum(0) # rate / net divergence / r


        rate_corrected_d = d.copy()
        for i in range(len(nodes)):
            for j in range(len(nodes)):
                rate_corrected_d[i,j] -= (r[i] + r[j]) / (len(nodes)-2)
        diag_num = rate_corrected_d.max() + 1
        for i in range(len(nodes)):
            rate_corrected_d[i,i] = diag_num
            # will soon need to find minimum off diagonal
        # print rate_corrected_d
        i, j = numpy.unravel_index(rate_corrected_d.argmin(),rate_corrected_d.shape)


        num_nodes = len(nodes) # len will change below
        n1 = nodes[i]
        n1.length = .5 * (d[i,j] + (r[i] - r[j])/(num_nodes-2))

        n2 = nodes[j]
        n2.length = .5 * (d[i,j] + (r[j] - r[i])/(num_nodes-2))


        # no negative branch lengths
        if no_negatives:
            nodes[i].length = max(0.0, nodes[i].length)
            nodes[j].length = max(0.0, nodes[j].length)
        
        # Join i and k to make new node
        new_node = TreeNode(children=[nodes[i], nodes[j]])
        
        # Store new node at i
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
    
        # print 'd:', d, [node.name for node in nodes]

    # 2 left
    # if len(nodes[0].children) < len(nodes[1].children):
    #     nodes.reverse()


    final_dist =  .5 * d[0,1]
    nodes[0].length = final_dist
    nodes[1].length = final_dist

    root_node = TreeNode(children=nodes)


    return root_node
