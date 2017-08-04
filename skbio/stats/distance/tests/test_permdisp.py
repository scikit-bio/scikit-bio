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
import numpy.testing as npt
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
        self.exp_index = ['method name', 'test statistic name', 'sample size',
                          'number of groups', 'test statistic', 'p-value',
                          'number of permutations']
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
        
        self.assert_series_equal = partial(assert_series_equal,
                                           check_index_type=True,
                                           check_series_type=True)
        

    
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
        obs = permdisp(self.eq_mat, self.grouping_eq, permutations=0)
        
        pval = obs['p-value']
        np.isnan(pval)

    def test_centroid_normal(self):
        exp = pd.Series(index=self.exp_index, 
                        data=['PERMDISP', 'F-value', 9, 2, 0.244501519876, .63, 99],
                        name='PERMDISP results')

        dm_ids = ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593', 'PC.607',
                'PC.634', 'PC.635', 'PC.636']
        
        dm = DistanceMatrix(
        [[0.0, 0.595483768391, 0.618074717633, 0.582763100909, 0.566949022108,
            0.714717232268, 0.772001731764, 0.690237118413, 0.740681707488],
	[0.595483768391, 0.0, 0.581427669668, 0.613726772383, 0.65945132763,
            0.745176523638, 0.733836123821, 0.720305073505, 0.680785600439],
	[0.618074717633, 0.581427669668, 0.0, 0.672149021573, 0.699416863323,
            0.71405573754, 0.759178215168, 0.689701276341, 0.725100672826],
	[0.582763100909, 0.613726772383, 0.672149021573, 0.0, 0.64756120797,
            0.666018240373, 0.66532968784, 0.650464714994, 0.632524644216],
        [0.566949022108, 0.65945132763, 0.699416863323, 0.64756120797, 0.0,
            0.703720200713, 0.748240937349, 0.73416971958, 0.727154987937],
	[0.714717232268, 0.745176523638, 0.71405573754, 0.666018240373,
            0.703720200713, 0.0, 0.707316869557, 0.636288883818, 0.699880573956],
	[0.772001731764, 0.733836123821, 0.759178215168, 0.66532968784,
            0.748240937349, 0.707316869557, 0.0, 0.565875193399, 0.560605525642],
	[0.690237118413, 0.720305073505, 0.689701276341, 0.650464714994,
            0.73416971958, 0.636288883818, 0.565875193399, 0.0, 0.575788039321],
	[0.740681707488, 0.680785600439, 0.725100672826, 0.632524644216,
            0.727154987937, 0.699880573956,0.560605525642, 0.575788039321, 0.0]],
        dm_ids)

        grouping = ['Control', 'Control', 'Control', 'Control', 'Control',
                'Fast', 'Fast', 'Fast', 'Fast']
        
        np.random.seed(0)
        obs = permdisp(dm, grouping, test='centroid', permutations=99)

        self.assert_series_equal(obs, exp)

    def test_median_normal(self):
        #so far I expect this one to fail, still dont have median working
        exp = pd.Series(index=self.exp_index,
                        data=['PERMDISP', 'F-value', 9, 2, 0.1395, 0.7199, 99],
                        name='PERMDISP results')

        #this distance matrix nd ids should really be established during setup
        #make sure to fix this plz
        dm_ids = ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593', 'PC.607',
            'PC.634', 'PC.635', 'PC.636']
            

        dm = DistanceMatrix(
        [[0.0, 0.595483768391, 0.618074717633, 0.582763100909, 0.566949022108,
            0.714717232268, 0.772001731764, 0.690237118413, 0.740681707488],
	[0.595483768391, 0.0, 0.581427669668, 0.613726772383, 0.65945132763,
            0.745176523638, 0.733836123821, 0.720305073505, 0.680785600439],
	[0.618074717633, 0.581427669668, 0.0, 0.672149021573, 0.699416863323,
            0.71405573754, 0.759178215168, 0.689701276341, 0.725100672826],
	[0.582763100909, 0.613726772383, 0.672149021573, 0.0, 0.64756120797,
            0.666018240373, 0.66532968784, 0.650464714994, 0.632524644216],
        [0.566949022108, 0.65945132763, 0.699416863323, 0.64756120797, 0.0,
            0.703720200713, 0.748240937349, 0.73416971958, 0.727154987937],
	[0.714717232268, 0.745176523638, 0.71405573754, 0.666018240373,
            0.703720200713, 0.0, 0.707316869557, 0.636288883818, 0.699880573956],
	[0.772001731764, 0.733836123821, 0.759178215168, 0.66532968784,
            0.748240937349, 0.707316869557, 0.0, 0.565875193399, 0.560605525642],
	[0.690237118413, 0.720305073505, 0.689701276341, 0.650464714994,
            0.73416971958, 0.636288883818, 0.565875193399, 0.0, 0.575788039321],
	[0.740681707488, 0.680785600439, 0.725100672826, 0.632524644216,
            0.727154987937, 0.699880573956, 0.560605525642, 0.575788039321, 0.0]],
        dm_ids)

        grouping = ['Control', 'Control', 'Control', 'Control', 'Control',
                'Fast', 'Fast', 'Fast', 'Fast']

        np.random.seed(0)
        obs = permdisp(dm, grouping, test='median', permutations=99)

        self.assertLess(obs['p-value'], exp['p-value'])
        self.assertEqual(obs['test statistic'], exp['test statistic'])

    def test_not_distance_matrix(self):
        dm = []
        grouping = ['Control', 'Control', 'Control', 'Control', 'Control',
                'Fast', 'Fast', 'Fast', 'Fast']

        npt.assert_raises(TypeError, permdisp, dm, grouping, permutations=0)

    def test_mismatched_group(self):
        dm_ids = ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593', 'PC.607',
            'PC.634', 'PC.635', 'PC.636']
        
        dm = DistanceMatrix(
        [[0.0, 0.595483768391, 0.618074717633, 0.582763100909, 0.566949022108,
            0.714717232268, 0.772001731764, 0.690237118413, 0.740681707488],
	[0.595483768391, 0.0, 0.581427669668, 0.613726772383, 0.65945132763,
            0.745176523638, 0.733836123821, 0.720305073505, 0.680785600439],
	[0.618074717633, 0.581427669668, 0.0, 0.672149021573, 0.699416863323,
            0.71405573754, 0.759178215168, 0.689701276341, 0.725100672826],
	[0.582763100909, 0.613726772383, 0.672149021573, 0.0, 0.64756120797,
            0.666018240373, 0.66532968784, 0.650464714994, 0.632524644216],
        [0.566949022108, 0.65945132763, 0.699416863323, 0.64756120797, 0.0,
            0.703720200713, 0.748240937349, 0.73416971958, 0.727154987937],
	[0.714717232268, 0.745176523638, 0.71405573754, 0.666018240373,
            0.703720200713, 0.0, 0.707316869557, 0.636288883818, 0.699880573956],
	[0.772001731764, 0.733836123821, 0.759178215168, 0.66532968784,
            0.748240937349, 0.707316869557, 0.0, 0.565875193399, 0.560605525642],
	[0.690237118413, 0.720305073505, 0.689701276341, 0.650464714994,
            0.73416971958, 0.636288883818, 0.565875193399, 0.0, 0.575788039321],
	[0.740681707488, 0.680785600439, 0.725100672826, 0.632524644216,
            0.727154987937, 0.699880573956, 0.560605525642, 0.575788039321, 0.0]],
        dm_ids)

        gr = ['foo', 'bar']

        npt.assert_raises(ValueError, permdisp, dm, gr)

    def test_single_group(self):
        dm_ids = ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593', 'PC.607',
            'PC.634', 'PC.635', 'PC.636']
        
        dm = DistanceMatrix(
        [[0.0, 0.595483768391, 0.618074717633, 0.582763100909, 0.566949022108,
            0.714717232268, 0.772001731764, 0.690237118413, 0.740681707488],
	[0.595483768391, 0.0, 0.581427669668, 0.613726772383, 0.65945132763,
            0.745176523638, 0.733836123821, 0.720305073505, 0.680785600439],
	[0.618074717633, 0.581427669668, 0.0, 0.672149021573, 0.699416863323,
            0.71405573754, 0.759178215168, 0.689701276341, 0.725100672826],
	[0.582763100909, 0.613726772383, 0.672149021573, 0.0, 0.64756120797,
            0.666018240373, 0.66532968784, 0.650464714994, 0.632524644216],
        [0.566949022108, 0.65945132763, 0.699416863323, 0.64756120797, 0.0,
            0.703720200713, 0.748240937349, 0.73416971958, 0.727154987937],
	[0.714717232268, 0.745176523638, 0.71405573754, 0.666018240373,
            0.703720200713, 0.0, 0.707316869557, 0.636288883818, 0.699880573956],
	[0.772001731764, 0.733836123821, 0.759178215168, 0.66532968784,
            0.748240937349, 0.707316869557, 0.0, 0.565875193399, 0.560605525642],
	[0.690237118413, 0.720305073505, 0.689701276341, 0.650464714994,
            0.73416971958, 0.636288883818, 0.565875193399, 0.0, 0.575788039321],
	[0.740681707488, 0.680785600439, 0.725100672826, 0.632524644216,
            0.727154987937, 0.699880573956, 0.560605525642, 0.575788039321, 0.0]],
        dm_ids)

        gr = ['f', 'f', 'f', 'f', 'f', 'f', 'f', 'f', 'f']

        npt.assert_raises(ValueError, permdisp, dm, gr)


if __name__ == '__main__':
    main()

