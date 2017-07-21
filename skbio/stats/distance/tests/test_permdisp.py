# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from functools import partial
from unittest import TestCase, main

import numpy as np
import pandas as pd
from pandas.util.testing import assert_series_equal

from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permdisp
from skbio.stats.distance._permdisp import _compute_centroid_groups

class testPERMDISP(TestCase):
    
    def setUp(self):
        
        #test with 2 groups of equal size
        #when assigned different labels, results should be the same
        self.grouping_eq = ['foo', 'foo', 'foo', 'bar', 'bar', 'bar']
        self.grouping_eq_relab = ['pyt', 'pyt', 'pyt', 'hon', 'hon', 'hon']
        
        #test with 3 groups of different sizes
        #when assigned different labels results should be the same
        self.grouping_uneq = ['foo', 'foo', 'bar', 'bar', 'bar', 
                                    'qw', 'qw', 'qw', 'qw']

        self.grouping_uneq_relab = [12,12, 7, 7, 7, 23, 23, 23, 23]

        self.grouping_un_mixed = ['a', 'a', 7, 7, 7, 'b', 'b', 'b', 'b']
        
        eq_ids = ['s1', 's2', 's3', 's4', 's5', 's6']
        uneq_ids = ['s1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9']
        #matrix for equal grouping
        self.eq_mat = DistanceMatrix([[0, 4, 0, 0, 4, 2],
                                      [4, 0, 2, 0, 3, 1],
                                      [0, 2, 0, 5, 2, 5],
                                      [0, 0, 5, 0, 0, 2],
                                      [4, 3, 2, 0, 0, 2],
                                      [2, 1, 5, 2, 2, 0]], eq_ids)
        
        #matrix for unequal grouping
        self.uneq_mat = DistanceMatrix([[0, 0, 4, 0, 0, 3, 5, 3, 0],
                                        [0, 0, 0, 3, 4, 5, 3, 0, 3],
                                        [4, 0, 0, 4, 3, 1, 0, 5, 2],
                                        [0, 3, 4, 0, 0, 2, 1, 3, 5],
                                        [0, 4, 3, 0, 0, 1, 1, 5, 0],
                                        [3, 5, 1, 2, 1, 0, 2, 0, 5],
                                        [5, 3, 0, 1, 1, 2, 0, 4, 3],
                                        [3, 0, 5, 3, 5, 0, 4, 0, 4],
                                        [0, 3, 2, 5, 0, 5, 3, 4, 0]], uneq_ids)

        #null matrix for equal grouping
        self.null_mat = DistanceMatrix([[0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0],
                                        [0, 0, 0, 0, 0, 0]], eq_ids)
        

    
    def test_centroids_eq_groups(self):
        exp = [[1.2886811963240687, 1.890538910062923, 1.490527658097728],
                   [2.17349240061718, 2.3192679626679946, 2.028338553903792]]


        dm = pcoa(self.eq_mat)
        
        obs = _compute_centroid_groups(dm, self.grouping_eq)
        self.assertEqual(obs, exp)
        
        obs_relab = _compute_centroid_groups(dm, self.grouping_eq_relab)
        self.assertEqual(obs_relab, obs)

    def test_centroids_uneq_groups(self):
        exp = [[2.5847022428144935, 2.285624595858895, 
                                    1.7022431146340287],
                [1.724817266046108, 1.724817266046108],
                [2.4333280644972795, 2.389000390879655, 
                    2.8547180589306036, 3.218568759338847]]

        dm = pcoa(self.uneq_mat)

        obs = _compute_centroid_groups(dm, self.grouping_uneq)
        self.assertEqual(obs, exp)

        obs_relab = _compute_centroid_groups(dm, self.grouping_uneq_relab)
        self.assertEqual(obs, obs_relab)

    def test_centroids_mixedgroups(self):
        exp = [[2.5847022428144935, 2.285624595858895, 
                                    1.7022431146340287],
                [1.724817266046108, 1.724817266046108],
                [2.4333280644972795, 2.389000390879655, 
                    2.8547180589306036, 3.218568759338847]]
        dm = pcoa(self.uneq_mat)
        
        obs_mixed = _compute_centroid_groups(dm, self.grouping_un_mixed)
        self.assertEqual(exp, obs_mixed)

    def test_centroids_null(self):
        exp = [[0.00, 0.00, 0.00], [0.00, 0.00, 0.00]]
        dm = pcoa(self.null_mat)
        
        obs_null = _compute_centroid_groups(dm, self.grouping_eq)
        self.assertEqual(exp, obs_null)
    
    def test_no_permuations(self):
        stat, pval = permdisp(self.eq_mat, self.grouping_eq, permutations=0)
        print(pval, np.nan)
        np.isnan(pval)



if __name__ == '__main__':
    main()

