#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division

import unittest, os
import numpy
# from cogent.phylo.distance import *
from skbio.core.rnj import rnj
# from cogent import LoadTree
import skbio.core.distance
from skbio.core.tree import TreeNode
import random

class RnjTests(unittest.TestCase):
    def setUp(self):
        self.mtx = numpy.array(
        [   [0,1.5,3],
            [1.5,0,2.5],
            [3,2.5,0]])
        self.names = ['A','B','C']
        
        # example from Yves Van de Peer, pg 150
        self.distmtx = skbio.core.distance.DistanceMatrix(self.mtx,self.names)
        self.mtx2 = numpy.array(
        [   [0,5,4,7,6,8],
            [5,0,7,10,9,11],
            [4,7,0,7,6,8],
            [7,10,7,0,5,9],
            [6,9,6,5,0,8],
            [8,11,8,9,8,0]] )
        self.names2 = [char for char in 'ABCDEF']
        self.distmtx2 = skbio.core.distance.DistanceMatrix(self.mtx2,self.names2)

    def test_rnj(self):
        """test nj with 2 hand-worked examples"""
        res = rnj(self.distmtx)
        dists = res.tip_tip_distances()
        self.assertEqual(dists, self.distmtx)

        res2 = rnj(self.distmtx2)
        # print res2.ascii_art(with_distances=True)
        # print res2.to_newick(with_distances=True)
        self.assertEqual(res2.tip_tip_distances(), self.distmtx2)

    def test_build_recover(self):
        nodes = [TreeNode(name=char) for char in 'abcde']
        while len(nodes) > 1:
            random.seed(0)
            random.shuffle(nodes)
            n1 = nodes.pop()
            n1.length = round(numpy.random.gamma(1,1),3)
            n2 = nodes.pop()
            n2.length = 1
            new_parent = TreeNode(children=[n1,n2])
            nodes.append(new_parent)
        tree = nodes[0]
        print tree.ascii_art(with_distances=True)
        dmtx = tree.tip_tip_distances()
        rnj_tree = rnj(dmtx)
        print 'rnj:'
        print rnj_tree.ascii_art(with_distances=True)
        print dmtx
        rnj_dmtx = rnj_tree.tip_tip_distances()
        print rnj_dmtx
        print 'eq to orig?', rnj_dmtx == dmtx

if __name__ == '__main__':
    unittest.main()
