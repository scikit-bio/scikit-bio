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


class RnjTests(unittest.TestCase):
    def setUp(self):
        self.mtx = numpy.array(
        [   [0,1.5,3],
            [1.5,0,2.5],
            [3,2.5,0]])
        self.names = ['A','B','C']
        
        self.distmtx = skbio.core.distance.DistanceMatrix(self.mtx,self.names)
        self.mtx2 = numpy.array(
        [   [0,5,4,7,6,8],
            [5,0,7,10,9,11],
            [4,7,0,7,6,8],
            [7,10,7,0,5,8],
            [6,9,6,5,0,8],
            [8,11,8,8,8,0]] )
        self.names2 = [char for char in 'ABCDEF']

    def test_rnj(self):
        """testing (well, exercising at least), rnj"""
        res = rnj(self.distmtx)
        print res.ascii_art(with_distances=True)
        print res.to_newick(with_distances=True)
        dists = res.tip_tip_distances()
        print dists
        print self.distmtx
        print dists == self.distmtx
        # res2 = rnj(self.mtx2, self.names2)
        # print res2.ascii_art()
        # print res2.to_newick(with_distances=True)
        # print res2.tip_tip_distances()

if __name__ == '__main__':
    unittest.main()
