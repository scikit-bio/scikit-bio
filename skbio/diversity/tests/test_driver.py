# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main

import pandas as pd
import numpy as np
import numpy.testing as npt

from skbio.io._fileobject import StringIO
from skbio import DistanceMatrix, TreeNode
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.diversity.alpha import faith_pd, observed_otus
from skbio.diversity.beta import unweighted_unifrac, weighted_unifrac
from skbio.diversity._driver import (_validate_counts_vector,
                                     _validate_counts_vectors,
                                     _validate_otu_ids_and_tree,
                                     _vectorize_counts_and_tree)
from skbio.tree import DuplicateNodeError, MissingNodeError


class AlphaDiversityTests(TestCase):
    def setUp(self):
        self.table1 = np.array([[1, 3, 0, 1, 0],
                                [0, 2, 0, 4, 4],
                                [0, 0, 6, 2, 1],
                                [0, 0, 1, 1, 1]])
        self.sids1 = list('ABCD')
        self.oids1 = ['OTU%d' % i for i in range(1, 6)]
        self.t1 = TreeNode.read(StringIO(
            u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):'
            u'0.0,(OTU4:0.75,OTU5:0.75):1.25):0.0)root;'))

    def test_invalid_input(self):
        # number of ids doesn't match the number of samples
        self.assertRaises(ValueError, alpha_diversity, 'observed_otus',
                          self.table1, list('ABC'))

        # otu_ids not provided for metric=faith_pd
        self.assertRaises(ValueError, alpha_diversity, 'faith_pd', self.table1,
                          list('ABC'), tree=self.t1)
        # tree not provided for metric=faith_pd
        self.assertRaises(ValueError, alpha_diversity, 'faith_pd', self.table1,
                          list('ABC'), otu_ids=self.oids1)

        # unknown metric provided
        self.assertRaises(ValueError, alpha_diversity, 'not-a-metric',
                          self.table1)

    def test_observed_otus(self):
        # expected values hand-calculated
        expected = pd.Series([3, 3, 3, 3], index=self.sids1)
        actual = alpha_diversity('observed_otus', self.table1, self.sids1)
        npt.assert_array_equal(actual, expected)
        # function passed instead of string
        actual = alpha_diversity(observed_otus, self.table1, self.sids1)
        npt.assert_array_equal(actual, expected)

    def test_no_ids(self):
        # expected values hand-calculated
        expected = pd.Series([3, 3, 3, 3])
        actual = alpha_diversity('observed_otus', self.table1)
        npt.assert_array_equal(actual, expected)

    def test_faith_pd(self):
        # calling faith_pd through alpha_diversity gives same results as
        # calling it directly
        expected = []
        for e in self.table1:
            expected.append(faith_pd(e, tree=self.t1, otu_ids=self.oids1))
        expected = pd.Series(expected)
        actual = alpha_diversity('faith_pd', self.table1, tree=self.t1,
                                 otu_ids=self.oids1)
        npt.assert_array_equal(actual, expected)

    def test_optimized(self):
        # calling optimized faith_pd gives same results as calling unoptimized
        # version
        optimized = alpha_diversity('faith_pd', self.table1, tree=self.t1,
                                    otu_ids=self.oids1)
        unoptimized = alpha_diversity(faith_pd, self.table1, tree=self.t1,
                                      otu_ids=self.oids1)
        npt.assert_array_equal(optimized, unoptimized)


class BetaDiversityTests(TestCase):
    def setUp(self):
        self.t1 = [[1, 5],
                   [2, 3],
                   [0, 1]]
        self.ids1 = list('ABC')
        self.tree1 = TreeNode.read(StringIO(
            '((O1:0.25, O2:0.50):0.25, O3:0.75)root;'))
        self.otu_ids1 = ['O1', 'O2']

        self.t2 = [[23, 64, 14, 0, 0, 3, 1],
                   [0, 3, 35, 42, 0, 12, 1],
                   [0, 5, 5, 0, 40, 40, 0],
                   [44, 35, 9, 0, 1, 0, 0],
                   [0, 2, 8, 0, 35, 45, 1],
                   [0, 0, 25, 35, 0, 19, 0]]
        self.ids2 = list('ABCDEF')

    def test_invalid_input(self):
        # number of ids doesn't match the number of samples
        self.assertRaises(ValueError, beta_diversity, self.t1, list('AB'),
                          'euclidean')

    def test_euclidean(self):
        actual_dm = beta_diversity('euclidean', self.t1, self.ids1)
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

        actual_dm = beta_diversity('euclidean', self.t2, self.ids2)
        expected_data = [
            [0., 80.8455317, 84.0297566, 36.3042697, 86.0116271, 78.9176786],
            [80.8455317, 0., 71.0844568, 74.4714710, 69.3397433, 14.422205],
            [84.0297566, 71.0844568, 0., 77.2851861, 8.3066238, 60.7536007],
            [36.3042697, 74.4714710, 77.2851861, 0., 78.7908624, 70.7389567],
            [86.0116271, 69.3397433, 8.3066238, 78.7908624, 0., 58.4807660],
            [78.9176786, 14.422205, 60.7536007, 70.7389567, 58.4807660, 0.]]
        expected_dm = DistanceMatrix(expected_data, self.ids2)
        for id1 in self.ids2:
            for id2 in self.ids2:
                npt.assert_almost_equal(actual_dm[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_braycurtis(self):
        actual_dm = beta_diversity('braycurtis', self.t1, self.ids1)
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

        actual_dm = beta_diversity('braycurtis', self.t2, self.ids2)
        expected_data = [
            [0., 0.78787879, 0.86666667, 0.30927835, 0.85714286, 0.81521739],
            [0.78787879, 0., 0.78142077, 0.86813187, 0.75, 0.1627907],
            [0.86666667, 0.78142077, 0., 0.87709497, 0.09392265, 0.71597633],
            [0.30927835, 0.86813187, 0.87709497, 0., 0.87777778, 0.89285714],
            [0.85714286, 0.75, 0.09392265, 0.87777778, 0., 0.68235294],
            [0.81521739, 0.1627907, 0.71597633, 0.89285714, 0.68235294, 0.]]
        expected_dm = DistanceMatrix(expected_data, self.ids2)
        for id1 in self.ids2:
            for id2 in self.ids2:
                npt.assert_almost_equal(actual_dm[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_unweighted_unifrac(self):
        # expected values calculated by hand
        dm1 = beta_diversity('unweighted_unifrac', self.t1, self.ids1,
                             otu_ids=self.otu_ids1, tree=self.tree1)
        dm2 = beta_diversity(unweighted_unifrac, self.t1, self.ids1,
                             otu_ids=self.otu_ids1, tree=self.tree1)
        self.assertEqual(dm1.shape, (3, 3))
        self.assertEqual(dm1, dm2)
        expected_data = [
            [0.0, 0.0, 0.25/1.0],
            [0.0, 0.0, 0.25/1.0],
            [0.25/1.0, 0.25/1.0, 0.0]]
        expected_dm = DistanceMatrix(expected_data, ids=self.ids1)
        for id1 in self.ids1:
            for id2 in self.ids1:
                npt.assert_almost_equal(dm1[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_weighted_unifrac(self):
        # expected values calculated by hand
        dm1 = beta_diversity('weighted_unifrac', self.t1, self.ids1,
                             otu_ids=self.otu_ids1, tree=self.tree1)
        dm2 = beta_diversity(weighted_unifrac, self.t1, self.ids1,
                             otu_ids=self.otu_ids1, tree=self.tree1)
        self.assertEqual(dm1.shape, (3, 3))
        self.assertEqual(dm1, dm2)
        expected_data = [
            [0.0, 0.1750000, 0.12499999],
            [0.1750000, 0.0, 0.3000000],
            [0.12499999, 0.3000000, 0.0]]
        expected_dm = DistanceMatrix(expected_data, ids=self.ids1)
        for id1 in self.ids1:
            for id2 in self.ids1:
                npt.assert_almost_equal(dm1[id1, id2],
                                        expected_dm[id1, id2], 6)

    def test_weighted_unifrac_normalized(self):
        # expected values calculated by hand
        dm1 = beta_diversity('weighted_unifrac', self.t1, self.ids1,
                             otu_ids=self.otu_ids1, tree=self.tree1,
                             normalized=True)
        dm2 = beta_diversity(weighted_unifrac, self.t1, self.ids1,
                             otu_ids=self.otu_ids1, tree=self.tree1,
                             normalized=True)
        self.assertEqual(dm1.shape, (3, 3))
        self.assertEqual(dm1, dm2)
        expected_data = [
            [0.0, 0.128834, 0.085714],
            [0.128834, 0.0, 0.2142857],
            [0.085714, 0.2142857, 0.0]]
        expected_dm = DistanceMatrix(expected_data, ids=self.ids1)
        for id1 in self.ids1:
            for id2 in self.ids1:
                npt.assert_almost_equal(dm1[id1, id2],
                                        expected_dm[id1, id2], 6)


class UtilityFunctionTests(TestCase):

    def test_validate_counts_vector(self):
        # python list
        obs = _validate_counts_vector([0, 2, 1, 3])
        npt.assert_array_equal(obs, np.array([0, 2, 1, 3]))
        self.assertEqual(obs.dtype, int)

        # numpy array (no copy made)
        data = np.array([0, 2, 1, 3])
        obs = _validate_counts_vector(data)
        npt.assert_array_equal(obs, data)
        self.assertEqual(obs.dtype, int)
        self.assertTrue(obs is data)

        # single element
        obs = _validate_counts_vector([42])
        npt.assert_array_equal(obs, np.array([42]))
        self.assertEqual(obs.dtype, int)
        self.assertEqual(obs.shape, (1,))

        # suppress casting to int
        obs = _validate_counts_vector([42.2, 42.1, 0], suppress_cast=True)
        npt.assert_array_equal(obs, np.array([42.2, 42.1, 0]))
        self.assertEqual(obs.dtype, float)

        # all zeros
        obs = _validate_counts_vector([0, 0, 0])
        npt.assert_array_equal(obs, np.array([0, 0, 0]))
        self.assertEqual(obs.dtype, int)

        # all zeros (single value)
        obs = _validate_counts_vector([0])
        npt.assert_array_equal(obs, np.array([0]))
        self.assertEqual(obs.dtype, int)

    def test_validate_counts_vector_invalid_input(self):
        # wrong dtype
        with self.assertRaises(TypeError):
            _validate_counts_vector([0, 2, 1.2, 3])

        # wrong number of dimensions (2-D)
        with self.assertRaises(ValueError):
            _validate_counts_vector([[0, 2, 1, 3], [4, 5, 6, 7]])

        # wrong number of dimensions (scalar)
        with self.assertRaises(ValueError):
            _validate_counts_vector(1)

        # negative values
        with self.assertRaises(ValueError):
            _validate_counts_vector([0, 0, 2, -1, 3])

    def test_validate_counts_vectors(self):
        # basic valid input (n=2)
        obs_u, obs_v = _validate_counts_vectors([0, 1, 1, 0, 2],
                                                [0, 0, 2, 1, 3])
        npt.assert_array_equal(obs_u, np.array([0, 1, 1, 0, 2]))
        npt.assert_array_equal(obs_v, np.array([0, 0, 2, 1, 3]))

        # basic valid input (n=3)
        actual = _validate_counts_vectors([0, 1, 1, 0, 2],
                                          [0, 0, 2, 1, 3],
                                          [1, 1, 1, 1, 1])
        npt.assert_array_equal(actual[0], np.array([0, 1, 1, 0, 2]))
        npt.assert_array_equal(actual[1], np.array([0, 0, 2, 1, 3]))
        npt.assert_array_equal(actual[2], np.array([1, 1, 1, 1, 1]))

        # empty counts vectors
        obs_u, obs_v = _validate_counts_vectors(np.array([], dtype=int),
                                                np.array([], dtype=int))
        npt.assert_array_equal(obs_u, np.array([]))
        npt.assert_array_equal(obs_v, np.array([]))

    def test_validate_counts_vectors_suppress_cast(self):
        # suppress_cast is passed through to _validate_counts_vector
        obs_u, obs_v = _validate_counts_vectors(
            [42.2, 42.1, 0], [42.2, 42.1, 1.0], suppress_cast=True)
        npt.assert_array_equal(obs_u, np.array([42.2, 42.1, 0]))
        npt.assert_array_equal(obs_v, np.array([42.2, 42.1, 1.0]))
        self.assertEqual(obs_u.dtype, float)
        self.assertEqual(obs_v.dtype, float)
        with self.assertRaises(TypeError):
            _validate_counts_vectors([0.0], [1], suppress_cast=False)

    def test_validate_counts_vectors_invalid_input(self):
        # checks that are caught by the calls to _validate_counts_vector
        with self.assertRaises(ValueError):
            _validate_counts_vectors([0, 1, 1, 0, 2], [0, 0, 2, -1, 3])
        with self.assertRaises(ValueError):
            _validate_counts_vectors([0, 0, 2, -1, 3], [0, 1, 1, 0, 2])

        # len of vectors not equal
        u_counts = [1, 2]
        v_counts = [1, 1, 1]
        self.assertRaises(ValueError, _validate_counts_vectors, u_counts,
                          v_counts)
        u_counts = [1, 2, 3]
        v_counts = [1, 1]
        self.assertRaises(ValueError, _validate_counts_vectors, u_counts,
                          v_counts)

    def test_validate_otu_ids_and_tree(self):
        # basic valid input
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

        # all tips observed
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 1, 1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5']
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

        # no tips observed
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = []
        otu_ids = []
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

        # all counts zero
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [0, 0, 0, 0, 0]
        otu_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5']
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

        # single node tree
        t = TreeNode.read(StringIO(u'root;'))
        counts = []
        otu_ids = []
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

    def test_validate_otu_ids_and_tree_invalid_input(self):
        # tree has duplicated tip ids
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU2:0.75):1.25):0.0)root;'))
        counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(DuplicateNodeError, _validate_otu_ids_and_tree,
                          counts, otu_ids, t)

        # unrooted tree as input
        t = TreeNode.read(StringIO(u'((OTU1:0.1, OTU2:0.2):0.3, OTU3:0.5,'
                                   u'OTU4:0.7);'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # otu_ids has duplicated ids
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU2']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # len of vectors not equal
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # tree with no branch lengths
        t = TreeNode.read(
            StringIO(u'((((OTU1,OTU2),OTU3)),(OTU4,OTU5));'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # tree missing some branch lengths
        t = TreeNode.read(
            StringIO(u'(((((OTU1,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # otu_ids not present in tree
        t = TreeNode.read(
            StringIO(u'(((((OTU1:0.25,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                     u'0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU32']
        self.assertRaises(MissingNodeError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

    def test_vectorize_counts_and_tree(self):
        t = TreeNode.read(StringIO(u"((a:1, b:2)c:3)root;"))
        counts = np.array([[0, 1], [1, 5], [10, 1]])
        count_array, indexed, branch_lengths = \
            _vectorize_counts_and_tree(counts, np.array(['a', 'b']), t)
        exp_counts = np.array([[0, 1, 10], [1, 5, 1], [1, 6, 11], [1, 6, 11]])
        npt.assert_equal(count_array, exp_counts.T)

    def test_vectorize_counts_and_tree_w_precomputed_index(self):
        t = TreeNode.read(StringIO(u"((a:1, b:2)c:3)root;"))
        counts = np.array([[0, 1], [1, 5], [10, 1]])
        indexed = t.to_array(nan_length_value=0.0)
        count_array, indexed, branch_lengths = _vectorize_counts_and_tree(
            counts, np.array(['a', 'b']), t, indexed)
        exp_counts = np.array([[0, 1, 10], [1, 5, 1], [1, 6, 11], [1, 6, 11]])
        npt.assert_equal(count_array, exp_counts.T)


if __name__ == "__main__":
    main()
