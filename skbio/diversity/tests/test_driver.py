# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

import pandas as pd
import numpy as np
import numpy.testing as npt

from skbio import DistanceMatrix, TreeNode
from skbio.util._testing import assert_series_almost_equal
from skbio.diversity import (alpha_diversity, beta_diversity,
                             partial_beta_diversity,
                             get_alpha_diversity_metrics,
                             get_beta_diversity_metrics)
from skbio.diversity.alpha import faith_pd, observed_otus
from skbio.diversity.beta import unweighted_unifrac, weighted_unifrac
from skbio.tree import DuplicateNodeError, MissingNodeError


class AlphaDiversityTests(TestCase):
    def setUp(self):
        self.table1 = np.array([[1, 3, 0, 1, 0],
                                [0, 2, 0, 4, 4],
                                [0, 0, 6, 2, 1],
                                [0, 0, 1, 1, 1]])
        self.sids1 = list('ABCD')
        self.oids1 = ['OTU%d' % i for i in range(1, 6)]
        self.tree1 = TreeNode.read(io.StringIO(
            '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):'
            '0.0,(OTU4:0.75,OTU5:0.75):1.25):0.0)root;'))

        self.table2 = np.array([[1, 3],
                                [0, 2],
                                [0, 0]])
        self.sids2 = list('xyz')
        self.oids2 = ['OTU1', 'OTU5']
        self.tree2 = TreeNode.read(io.StringIO(
            '(((((OTU1:42.5,OTU2:0.5):0.5,OTU3:1.0):1.0):'
            '0.0,(OTU4:0.75,OTU5:0.0001):1.25):0.0)root;'))

    def test_invalid_input(self):
        # number of ids doesn't match the number of samples
        self.assertRaises(ValueError, alpha_diversity, 'observed_otus',
                          self.table1, list('ABC'))

        # unknown metric provided
        self.assertRaises(ValueError, alpha_diversity, 'not-a-metric',
                          self.table1)

        # 3-D list provided as input
        self.assertRaises(ValueError, alpha_diversity, 'observed_otus',
                          [[[43]]])

        # negative counts
        self.assertRaises(ValueError, alpha_diversity, 'observed_otus',
                          [0, 3, -12, 42])

        # additional kwargs
        self.assertRaises(TypeError, alpha_diversity, 'observed_otus',
                          [0, 1], not_a_real_kwarg=42.0)
        self.assertRaises(TypeError, alpha_diversity, 'faith_pd',
                          [0, 1], tree=self.tree1, otu_ids=['OTU1', 'OTU2'],
                          not_a_real_kwarg=42.0)
        self.assertRaises(TypeError, alpha_diversity, faith_pd,
                          [0, 1], tree=self.tree1, otu_ids=['OTU1', 'OTU2'],
                          not_a_real_kwarg=42.0)

    def test_invalid_input_phylogenetic(self):
        # otu_ids not provided
        self.assertRaises(ValueError, alpha_diversity, 'faith_pd', self.table1,
                          list('ABC'), tree=self.tree1)
        # tree not provided
        self.assertRaises(ValueError, alpha_diversity, 'faith_pd', self.table1,
                          list('ABC'), otu_ids=self.oids1)

        # tree has duplicated tip ids
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU2:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(DuplicateNodeError, alpha_diversity, 'faith_pd',
                          counts, otu_ids=otu_ids, tree=t)

        # unrooted tree as input
        t = TreeNode.read(io.StringIO(
            '((OTU1:0.1, OTU2:0.2):0.3, OTU3:0.5,OTU4:0.7);'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, alpha_diversity, 'faith_pd',
                          counts, otu_ids=otu_ids, tree=t)

        # otu_ids has duplicated ids
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU2:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU2']
        self.assertRaises(ValueError, alpha_diversity, 'faith_pd',
                          counts, otu_ids=otu_ids, tree=t)

        # count and OTU vectors are not equal length
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU2:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2']
        self.assertRaises(ValueError, alpha_diversity, 'faith_pd',
                          counts, otu_ids=otu_ids, tree=t)
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU2:0.75):1.25):0.0)root;'))
        counts = [1, 2]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, alpha_diversity, 'faith_pd',
                          counts, otu_ids=otu_ids, tree=t)

        # tree with no branch lengths
        t = TreeNode.read(
            io.StringIO('((((OTU1,OTU2),OTU3)),(OTU4,OTU5));'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, alpha_diversity, 'faith_pd',
                          counts, otu_ids=otu_ids, tree=t)

        # tree missing some branch lengths
        t = TreeNode.read(
            io.StringIO('(((((OTU1,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                        '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, alpha_diversity, 'faith_pd',
                          counts, otu_ids=otu_ids, tree=t)

        # some otu_ids not present in tree
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU42']
        self.assertRaises(MissingNodeError, alpha_diversity, 'faith_pd',
                          counts, otu_ids=otu_ids, tree=t)

    def test_empty(self):
        # empty vector
        actual = alpha_diversity('observed_otus', np.array([], dtype=np.int64))
        expected = pd.Series([0])
        assert_series_almost_equal(actual, expected)

        # array of empty vector
        actual = alpha_diversity('observed_otus',
                                 np.array([[]], dtype=np.int64))
        expected = pd.Series([0])
        assert_series_almost_equal(actual, expected)

        # array of empty vectors
        actual = alpha_diversity('observed_otus',
                                 np.array([[], []], dtype=np.int64))
        expected = pd.Series([0, 0])
        assert_series_almost_equal(actual, expected)

        # empty vector
        actual = alpha_diversity('faith_pd', np.array([], dtype=np.int64),
                                 tree=self.tree1, otu_ids=[])
        expected = pd.Series([0.])
        assert_series_almost_equal(actual, expected)

        # array of empty vector
        actual = alpha_diversity('faith_pd',
                                 np.array([[]], dtype=np.int64),
                                 tree=self.tree1, otu_ids=[])
        expected = pd.Series([0.])
        assert_series_almost_equal(actual, expected)

        # array of empty vectors
        actual = alpha_diversity('faith_pd',
                                 np.array([[], []], dtype=np.int64),
                                 tree=self.tree1, otu_ids=[])
        expected = pd.Series([0., 0.])
        assert_series_almost_equal(actual, expected)

    def test_single_count_vector(self):
        actual = alpha_diversity('observed_otus', np.array([1, 0, 2]))
        expected = pd.Series([2])
        assert_series_almost_equal(actual, expected)

        actual = alpha_diversity('faith_pd', np.array([1, 3, 0, 1, 0]),
                                 tree=self.tree1, otu_ids=self.oids1)
        self.assertAlmostEqual(actual[0], 4.5)

    def test_input_types(self):
        list_result = alpha_diversity('observed_otus', [1, 3, 0, 1, 0])
        array_result = alpha_diversity('observed_otus',
                                       np.array([1, 3, 0, 1, 0]))
        self.assertAlmostEqual(list_result[0], 3)
        assert_series_almost_equal(list_result, array_result)

        list_result = alpha_diversity('faith_pd', [1, 3, 0, 1, 0],
                                      tree=self.tree1, otu_ids=self.oids1)
        array_result = alpha_diversity('faith_pd', np.array([1, 3, 0, 1, 0]),
                                       tree=self.tree1, otu_ids=self.oids1)
        self.assertAlmostEqual(list_result[0], 4.5)
        assert_series_almost_equal(list_result, array_result)

    def test_observed_otus(self):
        # expected values hand-calculated
        expected = pd.Series([3, 3, 3, 3], index=self.sids1)
        actual = alpha_diversity('observed_otus', self.table1, self.sids1)
        assert_series_almost_equal(actual, expected)
        # function passed instead of string
        actual = alpha_diversity(observed_otus, self.table1, self.sids1)
        assert_series_almost_equal(actual, expected)
        # alt input table
        expected = pd.Series([2, 1, 0], index=self.sids2)
        actual = alpha_diversity('observed_otus', self.table2, self.sids2)
        assert_series_almost_equal(actual, expected)

    def test_faith_pd(self):
        # calling faith_pd through alpha_diversity gives same results as
        # calling it directly
        expected = []
        for e in self.table1:
            expected.append(faith_pd(e, tree=self.tree1, otu_ids=self.oids1))
        expected = pd.Series(expected)
        actual = alpha_diversity('faith_pd', self.table1, tree=self.tree1,
                                 otu_ids=self.oids1)
        assert_series_almost_equal(actual, expected)

        # alt input table and tree
        expected = []
        for e in self.table2:
            expected.append(faith_pd(e, tree=self.tree2, otu_ids=self.oids2))
        expected = pd.Series(expected)
        actual = alpha_diversity('faith_pd', self.table2, tree=self.tree2,
                                 otu_ids=self.oids2)
        assert_series_almost_equal(actual, expected)

    def test_no_ids(self):
        # expected values hand-calculated
        expected = pd.Series([3, 3, 3, 3])
        actual = alpha_diversity('observed_otus', self.table1)
        assert_series_almost_equal(actual, expected)

    def test_optimized(self):
        # calling optimized faith_pd gives same results as calling unoptimized
        # version
        optimized = alpha_diversity('faith_pd', self.table1, tree=self.tree1,
                                    otu_ids=self.oids1)
        unoptimized = alpha_diversity(faith_pd, self.table1, tree=self.tree1,
                                      otu_ids=self.oids1)
        assert_series_almost_equal(optimized, unoptimized)


class BetaDiversityTests(TestCase):
    def setUp(self):
        self.table1 = [[1, 5],
                       [2, 3],
                       [0, 1]]
        self.sids1 = list('ABC')
        self.tree1 = TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        self.oids1 = ['O1', 'O2']

        self.table2 = [[23, 64, 14, 0, 0, 3, 1],
                       [0, 3, 35, 42, 0, 12, 1],
                       [0, 5, 5, 0, 40, 40, 0],
                       [44, 35, 9, 0, 1, 0, 0],
                       [0, 2, 8, 0, 35, 45, 1],
                       [0, 0, 25, 35, 0, 19, 0]]
        self.sids2 = list('ABCDEF')

    def test_qualitative_bug_issue_1549(self):
        mat = np.array([[42, 0, 37, 99, 1],
                        [12, 1, 22, 88, 0],
                        [25, 3, 23, 86, 0],
                        [0, 0, 87, 12, 0]])
        as_presence_absence = mat > 0
        obs_mat = beta_diversity('jaccard', mat)
        obs_presence_absence = beta_diversity('jaccard', as_presence_absence)
        self.assertEqual(obs_mat, obs_presence_absence)

    def test_invalid_input(self):
        # number of ids doesn't match the number of samples
        error_msg = (r"Number of rows")
        with self.assertRaisesRegex(ValueError, error_msg):
            beta_diversity(self.table1, list('AB'), 'euclidean')

        # unknown metric provided
        error_msg = r"not-a-metric"
        with self.assertRaisesRegex(ValueError, error_msg):
            beta_diversity('not-a-metric', self.table1)

        # 3-D list provided as input
        error_msg = (r"Only 1-D and 2-D")
        with self.assertRaisesRegex(ValueError, error_msg):
            beta_diversity('euclidean', [[[43]]])

        # negative counts
        error_msg = r"negative values."
        with self.assertRaisesRegex(ValueError, error_msg):
            beta_diversity('euclidean', [[0, 1, 3, 4], [0, 3, -12, 42]])
        with self.assertRaisesRegex(ValueError, error_msg):
            beta_diversity('euclidean', [[0, 1, 3, -4], [0, 3, 12, 42]])

        # additional kwargs
        error_msg = r"keyword argument"
        with self.assertRaisesRegex(TypeError, error_msg):
            beta_diversity('euclidean', [[0, 1, 3], [0, 3, 12]],
                           not_a_real_kwarg=42.0)
        with self.assertRaisesRegex(TypeError, error_msg):
            beta_diversity('unweighted_unifrac', [[0, 1, 3], [0, 3, 12]],
                           not_a_real_kwarg=42.0, tree=self.tree1,
                           otu_ids=['O1', 'O2', 'O3'])
        with self.assertRaisesRegex(TypeError, error_msg):
            beta_diversity('weighted_unifrac', [[0, 1, 3], [0, 3, 12]],
                           not_a_real_kwarg=42.0, tree=self.tree1,
                           otu_ids=['O1', 'O2', 'O3'])
        with self.assertRaisesRegex(TypeError, error_msg):
            beta_diversity(weighted_unifrac, [[0, 1, 3], [0, 3, 12]],
                           not_a_real_kwarg=42.0, tree=self.tree1,
                           otu_ids=['O1', 'O2', 'O3'])

    def test_invalid_input_phylogenetic(self):
        # otu_ids not provided
        self.assertRaises(ValueError, beta_diversity, 'weighted_unifrac',
                          self.table1, list('ABC'), tree=self.tree1)
        self.assertRaises(ValueError, beta_diversity, 'unweighted_unifrac',
                          self.table1, list('ABC'), tree=self.tree1)
        # tree not provided
        self.assertRaises(ValueError, beta_diversity, 'weighted_unifrac',
                          self.table1, list('ABC'), otu_ids=self.oids1)
        self.assertRaises(ValueError, beta_diversity, 'unweighted_unifrac',
                          self.table1, list('ABC'), otu_ids=self.oids1)

        # tree has duplicated tip ids
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU2:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(DuplicateNodeError, beta_diversity,
                          'weighted_unifrac', counts, otu_ids=otu_ids, tree=t)
        self.assertRaises(DuplicateNodeError, beta_diversity,
                          'unweighted_unifrac', counts, otu_ids=otu_ids,
                          tree=t)

        # unrooted tree as input
        t = TreeNode.read(io.StringIO('((OTU1:0.1, OTU2:0.2):0.3, OTU3:0.5,'
                                      'OTU4:0.7);'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, beta_diversity,
                          'weighted_unifrac', counts, otu_ids=otu_ids, tree=t)
        self.assertRaises(ValueError, beta_diversity,
                          'unweighted_unifrac', counts, otu_ids=otu_ids,
                          tree=t)

        # otu_ids has duplicated ids
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU2:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU2']
        self.assertRaises(ValueError, beta_diversity,
                          'weighted_unifrac', counts, otu_ids=otu_ids, tree=t)
        self.assertRaises(ValueError, beta_diversity,
                          'unweighted_unifrac', counts, otu_ids=otu_ids,
                          tree=t)

        # count and OTU vectors are not equal length
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU2:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2']
        self.assertRaises(ValueError, beta_diversity,
                          'weighted_unifrac', counts, otu_ids=otu_ids, tree=t)
        self.assertRaises(ValueError, beta_diversity,
                          'unweighted_unifrac', counts, otu_ids=otu_ids,
                          tree=t)
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU2:0.75):1.25):0.0)root;'))
        counts = [1, 2]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, beta_diversity,
                          'weighted_unifrac', counts, otu_ids=otu_ids, tree=t)
        self.assertRaises(ValueError, beta_diversity,
                          'unweighted_unifrac', counts, otu_ids=otu_ids,
                          tree=t)

        # tree with no branch lengths
        t = TreeNode.read(
            io.StringIO('((((OTU1,OTU2),OTU3)),(OTU4,OTU5));'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, beta_diversity,
                          'weighted_unifrac', counts, otu_ids=otu_ids, tree=t)
        self.assertRaises(ValueError, beta_diversity,
                          'unweighted_unifrac', counts, otu_ids=otu_ids,
                          tree=t)

        # tree missing some branch lengths
        t = TreeNode.read(
            io.StringIO('(((((OTU1,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                        '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, beta_diversity,
                          'weighted_unifrac', counts, otu_ids=otu_ids, tree=t)
        self.assertRaises(ValueError, beta_diversity,
                          'unweighted_unifrac', counts, otu_ids=otu_ids,
                          tree=t)

        # some otu_ids not present in tree
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU42']
        self.assertRaises(MissingNodeError, beta_diversity,
                          'weighted_unifrac', counts, otu_ids=otu_ids, tree=t)
        self.assertRaises(MissingNodeError, beta_diversity,
                          'unweighted_unifrac', counts, otu_ids=otu_ids,
                          tree=t)

    def test_empty(self):
        # array of empty vectors
        actual = beta_diversity('euclidean',
                                np.array([[], []], dtype=np.int64),
                                ids=['a', 'b'])
        expected_dm = DistanceMatrix([[0.0, 0.0], [0.0, 0.0]], ['a', 'b'])
        npt.assert_array_equal(actual, expected_dm)

        actual = beta_diversity('unweighted_unifrac',
                                np.array([[], []], dtype=np.int64),
                                ids=['a', 'b'], tree=self.tree1, otu_ids=[])
        expected_dm = DistanceMatrix([[0.0, 0.0], [0.0, 0.0]], ['a', 'b'])
        self.assertEqual(actual, expected_dm)

    def test_input_types(self):
        actual_array = beta_diversity('euclidean',
                                      np.array([[1, 5], [2, 3]]),
                                      ids=['a', 'b'])
        actual_list = beta_diversity('euclidean',
                                     [[1, 5], [2, 3]], ids=['a', 'b'])
        self.assertEqual(actual_array, actual_list)

    def test_euclidean(self):
        # TODO: update npt.assert_almost_equal calls to use DistanceMatrix
        # near-equality testing when that support is available
        actual_dm = beta_diversity('euclidean', self.table1, self.sids1)
        self.assertEqual(actual_dm.shape, (3, 3))
        npt.assert_almost_equal(actual_dm['A', 'A'], 0.0)
        npt.assert_almost_equal(actual_dm['B', 'B'], 0.0)
        npt.assert_almost_equal(actual_dm['C', 'C'], 0.0)
        npt.assert_almost_equal(actual_dm['A', 'B'], 2.23606798)
        npt.assert_almost_equal(actual_dm['B', 'A'], 2.23606798)
        npt.assert_almost_equal(actual_dm['A', 'C'], 4.12310563)
        npt.assert_almost_equal(actual_dm['C', 'A'], 4.12310563)
        npt.assert_almost_equal(actual_dm['B', 'C'], 2.82842712)
        npt.assert_almost_equal(actual_dm['C', 'B'], 2.82842712)

        actual_dm = beta_diversity('euclidean', self.table2, self.sids2)
        expected_data = [
            [0., 80.8455317, 84.0297566, 36.3042697, 86.0116271, 78.9176786],
            [80.8455317, 0., 71.0844568, 74.4714710, 69.3397433, 14.422205],
            [84.0297566, 71.0844568, 0., 77.2851861, 8.3066238, 60.7536007],
            [36.3042697, 74.4714710, 77.2851861, 0., 78.7908624, 70.7389567],
            [86.0116271, 69.3397433, 8.3066238, 78.7908624, 0., 58.4807660],
            [78.9176786, 14.422205, 60.7536007, 70.7389567, 58.4807660, 0.]]
        expected_dm = DistanceMatrix(expected_data, self.sids2)
        for id1 in self.sids2:
            for id2 in self.sids2:
                npt.assert_almost_equal(actual_dm[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_braycurtis(self):
        # TODO: update npt.assert_almost_equal calls to use DistanceMatrix
        # near-equality testing when that support is available
        actual_dm = beta_diversity('braycurtis', self.table1, self.sids1)
        self.assertEqual(actual_dm.shape, (3, 3))
        npt.assert_almost_equal(actual_dm['A', 'A'], 0.0)
        npt.assert_almost_equal(actual_dm['B', 'B'], 0.0)
        npt.assert_almost_equal(actual_dm['C', 'C'], 0.0)
        npt.assert_almost_equal(actual_dm['A', 'B'], 0.27272727)
        npt.assert_almost_equal(actual_dm['B', 'A'], 0.27272727)
        npt.assert_almost_equal(actual_dm['A', 'C'], 0.71428571)
        npt.assert_almost_equal(actual_dm['C', 'A'], 0.71428571)
        npt.assert_almost_equal(actual_dm['B', 'C'], 0.66666667)
        npt.assert_almost_equal(actual_dm['C', 'B'], 0.66666667)

        actual_dm = beta_diversity('braycurtis', self.table2, self.sids2)
        expected_data = [
            [0., 0.78787879, 0.86666667, 0.30927835, 0.85714286, 0.81521739],
            [0.78787879, 0., 0.78142077, 0.86813187, 0.75, 0.1627907],
            [0.86666667, 0.78142077, 0., 0.87709497, 0.09392265, 0.71597633],
            [0.30927835, 0.86813187, 0.87709497, 0., 0.87777778, 0.89285714],
            [0.85714286, 0.75, 0.09392265, 0.87777778, 0., 0.68235294],
            [0.81521739, 0.1627907, 0.71597633, 0.89285714, 0.68235294, 0.]]
        expected_dm = DistanceMatrix(expected_data, self.sids2)
        for id1 in self.sids2:
            for id2 in self.sids2:
                npt.assert_almost_equal(actual_dm[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_unweighted_unifrac(self):
        # TODO: update npt.assert_almost_equal calls to use DistanceMatrix
        # near-equality testing when that support is available
        # expected values calculated by hand
        dm1 = beta_diversity('unweighted_unifrac', self.table1, self.sids1,
                             otu_ids=self.oids1, tree=self.tree1)
        dm2 = beta_diversity(unweighted_unifrac, self.table1, self.sids1,
                             otu_ids=self.oids1, tree=self.tree1)
        self.assertEqual(dm1.shape, (3, 3))
        self.assertEqual(dm1, dm2)
        expected_data = [[0.0, 0.0, 0.25],
                         [0.0, 0.0, 0.25],
                         [0.25, 0.25, 0.0]]
        expected_dm = DistanceMatrix(expected_data, ids=self.sids1)
        for id1 in self.sids1:
            for id2 in self.sids1:
                npt.assert_almost_equal(dm1[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_weighted_unifrac(self):
        # TODO: update npt.assert_almost_equal calls to use DistanceMatrix
        # near-equality testing when that support is available
        # expected values calculated by hand
        dm1 = beta_diversity('weighted_unifrac', self.table1, self.sids1,
                             otu_ids=self.oids1, tree=self.tree1)
        dm2 = beta_diversity(weighted_unifrac, self.table1, self.sids1,
                             otu_ids=self.oids1, tree=self.tree1)
        self.assertEqual(dm1.shape, (3, 3))
        self.assertEqual(dm1, dm2)
        expected_data = [
            [0.0, 0.1750000, 0.12499999],
            [0.1750000, 0.0, 0.3000000],
            [0.12499999, 0.3000000, 0.0]]
        expected_dm = DistanceMatrix(expected_data, ids=self.sids1)
        for id1 in self.sids1:
            for id2 in self.sids1:
                npt.assert_almost_equal(dm1[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_weighted_unifrac_normalized(self):
        # TODO: update npt.assert_almost_equal calls to use DistanceMatrix
        # near-equality testing when that support is available
        # expected values calculated by hand
        dm1 = beta_diversity('weighted_unifrac', self.table1, self.sids1,
                             otu_ids=self.oids1, tree=self.tree1,
                             normalized=True)
        dm2 = beta_diversity(weighted_unifrac, self.table1, self.sids1,
                             otu_ids=self.oids1, tree=self.tree1,
                             normalized=True)
        self.assertEqual(dm1.shape, (3, 3))
        self.assertEqual(dm1, dm2)
        expected_data = [
            [0.0, 0.128834, 0.085714],
            [0.128834, 0.0, 0.2142857],
            [0.085714, 0.2142857, 0.0]]
        expected_dm = DistanceMatrix(expected_data, ids=self.sids1)
        for id1 in self.sids1:
            for id2 in self.sids1:
                npt.assert_almost_equal(dm1[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_scipy_kwargs(self):
        # confirm that p can be passed to SciPy's minkowski, and that it
        # gives a different result than not passing it (the off-diagonal
        # entries are not equal).
        dm1 = beta_diversity('minkowski', self.table1, self.sids1)
        dm2 = beta_diversity('minkowski', self.table1, self.sids1, p=42.0)

        for id1 in self.sids1:
            for id2 in self.sids1:
                if id1 != id2:
                    self.assertNotEqual(dm1[id1, id2], dm2[id1, id2])

    def test_alt_pairwise_func(self):
        # confirm that pairwise_func is actually being used
        def not_a_real_pdist(counts, metric):
            return [[0.0, 42.0], [42.0, 0.0]]
        dm1 = beta_diversity('unweighted_unifrac', self.table1,
                             otu_ids=self.oids1, tree=self.tree1,
                             pairwise_func=not_a_real_pdist)
        expected = DistanceMatrix([[0.0, 42.0], [42.0, 0.0]])
        self.assertEqual(dm1, expected)

        dm1 = beta_diversity('weighted_unifrac', self.table1,
                             otu_ids=self.oids1, tree=self.tree1,
                             pairwise_func=not_a_real_pdist)
        expected = DistanceMatrix([[0.0, 42.0], [42.0, 0.0]])
        self.assertEqual(dm1, expected)

        dm1 = beta_diversity(unweighted_unifrac, self.table1,
                             otu_ids=self.oids1, tree=self.tree1,
                             pairwise_func=not_a_real_pdist)
        expected = DistanceMatrix([[0.0, 42.0], [42.0, 0.0]])
        self.assertEqual(dm1, expected)

        dm1 = beta_diversity("euclidean", self.table1,
                             pairwise_func=not_a_real_pdist)
        expected = DistanceMatrix([[0.0, 42.0], [42.0, 0.0]])
        self.assertEqual(dm1, expected)


class MetricGetters(TestCase):

    def test_get_alpha_diversity_metrics(self):
        m = get_alpha_diversity_metrics()
        # basic sanity checks
        self.assertTrue('faith_pd' in m)
        self.assertTrue('chao1' in m)

    def test_get_alpha_diversity_metrics_sorted(self):
        m = get_alpha_diversity_metrics()
        n = sorted(list(m))
        self.assertEqual(m, n)

    def test_get_beta_diversity_metrics(self):
        m = get_beta_diversity_metrics()
        # basic sanity checks
        self.assertTrue('unweighted_unifrac' in m)
        self.assertTrue('weighted_unifrac' in m)

    def test_get_beta_diversity_metrics_sorted(self):
        m = get_beta_diversity_metrics()
        n = sorted(list(m))
        self.assertEqual(m, n)


class TestPartialBetaDiversity(TestCase):
    def setUp(self):
        self.table1 = [[1, 5],
                       [2, 3],
                       [0, 1]]
        self.sids1 = list('ABC')
        self.tree1 = TreeNode.read(io.StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        self.oids1 = ['O1', 'O2']

        self.table2 = [[23, 64, 14, 0, 0, 3, 1],
                       [0, 3, 35, 42, 0, 12, 1],
                       [0, 5, 5, 0, 40, 40, 0],
                       [44, 35, 9, 0, 1, 0, 0],
                       [0, 2, 8, 0, 35, 45, 1],
                       [0, 0, 25, 35, 0, 19, 0]]
        self.sids2 = list('ABCDEF')

    def test_id_pairs_as_iterable(self):
        id_pairs = iter([('B', 'C'), ])
        dm = partial_beta_diversity('unweighted_unifrac', self.table1,
                                    self.sids1, otu_ids=self.oids1,
                                    tree=self.tree1, id_pairs=id_pairs)
        self.assertEqual(dm.shape, (3, 3))
        expected_data = [[0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.25],
                         [0.0, 0.25, 0.0]]
        expected_dm = DistanceMatrix(expected_data, ids=self.sids1)
        for id1 in self.sids1:
            for id2 in self.sids1:
                npt.assert_almost_equal(dm[id1, id2],
                                        expected_dm[id1, id2], 6)

        # pass in iter(foo)

    def test_unweighted_unifrac_partial(self):
        # TODO: update npt.assert_almost_equal calls to use DistanceMatrix
        # near-equality testing when that support is available
        # expected values calculated by hand
        dm = partial_beta_diversity('unweighted_unifrac', self.table1,
                                    self.sids1, otu_ids=self.oids1,
                                    tree=self.tree1, id_pairs=[('B', 'C'), ])
        self.assertEqual(dm.shape, (3, 3))
        expected_data = [[0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.25],
                         [0.0, 0.25, 0.0]]
        expected_dm = DistanceMatrix(expected_data, ids=self.sids1)
        for id1 in self.sids1:
            for id2 in self.sids1:
                npt.assert_almost_equal(dm[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_weighted_unifrac_partial_full(self):
        # TODO: update npt.assert_almost_equal calls to use DistanceMatrix
        # near-equality testing when that support is available
        # expected values calculated by hand
        dm1 = partial_beta_diversity('weighted_unifrac', self.table1,
                                     self.sids1, otu_ids=self.oids1,
                                     tree=self.tree1, id_pairs=[('A', 'B'),
                                                                ('A', 'C'),
                                                                ('B', 'C')])
        dm2 = beta_diversity('weighted_unifrac', self.table1, self.sids1,
                             otu_ids=self.oids1, tree=self.tree1)

        self.assertEqual(dm1.shape, (3, 3))
        self.assertEqual(dm1, dm2)
        expected_data = [
            [0.0, 0.1750000, 0.12499999],
            [0.1750000, 0.0, 0.3000000],
            [0.12499999, 0.3000000, 0.0]]
        expected_dm = DistanceMatrix(expected_data, ids=self.sids1)
        for id1 in self.sids1:
            for id2 in self.sids1:
                npt.assert_almost_equal(dm1[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_self_self_pair(self):
        error_msg = (r"A duplicate or a self-self pair was observed.")
        with self.assertRaisesRegex(ValueError, error_msg):
            partial_beta_diversity((lambda x, y: x + y), self.table1,
                                   self.sids1, id_pairs=[('A', 'B'),
                                                         ('A', 'A')])

    def test_duplicate_pairs(self):
        # confirm that partial pairwise execution fails if duplicate pairs are
        # observed
        error_msg = (r"A duplicate or a self-self pair was observed.")
        with self.assertRaisesRegex(ValueError, error_msg):
            partial_beta_diversity((lambda x, y: x + y), self.table1,
                                   self.sids1, id_pairs=[('A', 'B'),
                                                         ('A', 'B')])

    def test_duplicate_transpose_pairs(self):
        # confirm that partial pairwise execution fails if a transpose
        # duplicate is observed
        error_msg = (r"A duplicate or a self-self pair was observed.")
        with self.assertRaisesRegex(ValueError, error_msg):
            partial_beta_diversity((lambda x, y: x + y), self.table1,
                                   self.sids1, id_pairs=[('A', 'B'),
                                                         ('A', 'B')])

    def test_pairs_not_subset(self):
        # confirm raise when pairs are not a subset of IDs
        error_msg = (r"`id_pairs` are not a subset of `ids`")
        with self.assertRaisesRegex(ValueError, error_msg):
            partial_beta_diversity((lambda x, y: x + y), self.table1,
                                   self.sids1, id_pairs=[('x', 'b'), ])

    def test_euclidean(self):
        # confirm that pw execution through partial is identical
        def euclidean(u, v, **kwargs):
            return np.sqrt(((u - v)**2).sum())

        id_pairs = [('A', 'B'), ('B', 'F'), ('D', 'E')]
        actual_dm = partial_beta_diversity(euclidean, self.table2, self.sids2,
                                           id_pairs=id_pairs)
        actual_dm = DistanceMatrix(actual_dm, self.sids2)

        expected_data = [
            [0., 80.8455317, 0., 0., 0., 0.],
            [80.8455317, 0., 0., 0., 0., 14.422205],
            [0., 0., 0., 0., 0., 0.],
            [0., 0., 0., 0., 78.7908624, 0.],
            [0., 0., 0., 78.7908624, 0., 0.],
            [0., 14.422205, 0., 0., 0., 0.]]

        expected_dm = DistanceMatrix(expected_data, self.sids2)
        for id1 in self.sids2:
            for id2 in self.sids2:
                npt.assert_almost_equal(actual_dm[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_unusable_metric(self):
        id_pairs = [('A', 'B'), ('B', 'F'), ('D', 'E')]
        error_msg = r"partial_beta_diversity is only compatible"
        with self.assertRaisesRegex(ValueError, error_msg):
            partial_beta_diversity('hamming', self.table2, self.sids2,
                                   id_pairs=id_pairs)


if __name__ == "__main__":
    main()
