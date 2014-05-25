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
        # Eliminate one node per iteration until 2 left
        num_nodes = len(nodes)
        
        # compute r (normalized), the sum of all pairwise distances
        # the normalization is over (num - 2), since later for a given i, j
        # distance(i, j) will be removed, and distance(i, i) = 0 always
        r = numpy.sum(d, 0) * 1./(num_nodes-2.)
        
        # find two nodes i, j that are minimize each other's 
        # transformed distance
        node_indices = range(num_nodes)
        if randomize == True:
            shuffle(node_indices)
        chose_pair = False
        
        # coefficient used calculating transformed distances
        coef = num_nodes * 1./(num_nodes - 2.)
        for i in node_indices:
        # find i's closest, call it j
        
            # xformed_dists is a list of T_i,j for all j
            xformed_dists = coef*d[i] - r - r[i]
        
            # give distance to self a bogus but nonminimum value
            xformed_dists[i] = numpy.abs(xformed_dists[0])*2. +\
                numpy.abs(xformed_dists[num_nodes - 1])*2.
        
            j = numpy.argmin(xformed_dists)
        
        
        # now find j's closest
            xformed_dists = coef*d[j] - r - r[j]
            xformed_dists[j] = numpy.abs(xformed_dists[0])*2. +\
                numpy.abs(xformed_dists[num_nodes - 1])*2.
            
            # if i and j are each other's minimum, choose this (i, j) pair
            if i == numpy.argmin(xformed_dists):
                # choose these i, j
                chose_pair = True
                break
        
        if not chose_pair:
            raise Exception("didn't choose a pair of nodes correctly")
        assert i != j, (i, j)
        
        # Branch lengths from i and j to new node
        nodes[i].length = 0.5 * (d[i,j] + r[i] - r[j])
        nodes[j].length = 0.5 * (d[i,j] + r[j] - r[i])
            
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
        
        # Eliminate j
        d[j, :] = d[num_nodes-1, :]
        d[:, j] = d[:, num_nodes-1]
        assert d[j, j] == 0.0, d
        d = d[0:num_nodes-1, 0:num_nodes-1]
        nodes[j] = nodes[num_nodes-1]
        nodes.pop()
    

    # 2 left
    # if len(nodes[0].children) < len(nodes[1].children):
    #     nodes.reverse()


    final_dist = 0.5 * d[0,1]
    nodes[0].length = final_dist
    nodes[1].length = final_dist

    root_node = TreeNode(children=nodes)


    return root_node
