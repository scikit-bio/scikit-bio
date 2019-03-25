# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

import numpy as np
import pandas as pd
import numpy.testing as npt

from skbio import TreeNode
from skbio.diversity._util import (_validate_counts_vector,
                                   _validate_counts_matrix,
                                   _validate_otu_ids_and_tree,
                                   _vectorize_counts_and_tree)
from skbio.tree import DuplicateNodeError, MissingNodeError


class ValidationTests(TestCase):

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

        # wrong number of dimensions (2-D)
        with self.assertRaises(ValueError):
            _validate_counts_vector([[0, 2, 1, 3], [4, 5, 6, 7]])

        # wrong number of dimensions (scalar)
        with self.assertRaises(ValueError):
            _validate_counts_vector(1)

        # negative values
        with self.assertRaises(ValueError):
            _validate_counts_vector([0, 0, 2, -1, 3])

        # strings
        with self.assertRaises(ValueError):
            _validate_counts_vector([0, 0, 'a', -1, 3])

    def test_validate_counts_matrix(self):
        # basic valid input (n=2)
        obs = _validate_counts_matrix([[0, 1, 1, 0, 2],
                                       [0, 0, 2, 1, 3]])
        npt.assert_array_equal(obs[0], np.array([0, 1, 1, 0, 2]))
        npt.assert_array_equal(obs[1], np.array([0, 0, 2, 1, 3]))

        # basic valid input (n=3)
        obs = _validate_counts_matrix([[0, 1, 1, 0, 2],
                                       [0, 0, 2, 1, 3],
                                       [1, 1, 1, 1, 1]])
        npt.assert_array_equal(obs[0], np.array([0, 1, 1, 0, 2]))
        npt.assert_array_equal(obs[1], np.array([0, 0, 2, 1, 3]))
        npt.assert_array_equal(obs[2], np.array([1, 1, 1, 1, 1]))

        # empty counts vectors
        obs = _validate_counts_matrix(np.array([[], []], dtype=int))
        npt.assert_array_equal(obs[0], np.array([]))
        npt.assert_array_equal(obs[1], np.array([]))

    def test_validate_counts_matrix_pandas(self):
        obs = _validate_counts_matrix(pd.DataFrame([[0, 1, 1, 0, 2],
                                                    [0, 0, 2, 1, 3],
                                                    [1, 1, 1, 1, 1]]))
        npt.assert_array_equal(obs[0], np.array([0, 1, 1, 0, 2]))
        npt.assert_array_equal(obs[1], np.array([0, 0, 2, 1, 3]))
        npt.assert_array_equal(obs[2], np.array([1, 1, 1, 1, 1]))

    def test_validate_counts_matrix_suppress_cast(self):
        # suppress_cast is passed through to _validate_counts_vector
        obs = _validate_counts_matrix(
            [[42.2, 42.1, 0], [42.2, 42.1, 1.0]], suppress_cast=True)
        npt.assert_array_equal(obs[0], np.array([42.2, 42.1, 0]))
        npt.assert_array_equal(obs[1], np.array([42.2, 42.1, 1.0]))
        self.assertEqual(obs[0].dtype, float)
        self.assertEqual(obs[1].dtype, float)

    def test_validate_counts_matrix_negative_counts(self):
        with self.assertRaises(ValueError):
            _validate_counts_matrix([[0, 1, 1, 0, 2], [0, 0, 2, -1, 3]])
        with self.assertRaises(ValueError):
            _validate_counts_matrix([[0, 0, 2, -1, 3], [0, 1, 1, 0, 2]])

    def test_validate_counts_matrix_unequal_lengths(self):
        # len of vectors not equal
        with self.assertRaises(ValueError):
            _validate_counts_matrix([[0], [0, 0], [9, 8]])
        with self.assertRaises(ValueError):
            _validate_counts_matrix([[0, 0], [0, 0, 8], [9, 8]])
        with self.assertRaises(ValueError):
            _validate_counts_matrix([[0, 0, 75], [0, 0, 3], [9, 8, 22, 44]])

    def test_validate_otu_ids_and_tree(self):
        # basic valid input
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

        # all tips observed
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 1, 1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5']
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

        # no tips observed
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = []
        otu_ids = []
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

        # all counts zero
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [0, 0, 0, 0, 0]
        otu_ids = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5']
        self.assertTrue(_validate_otu_ids_and_tree(counts, otu_ids, t) is None)

    def test_validate_otu_ids_and_tree_invalid_input(self):
        # tree has duplicated tip ids
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU2:0.75):1.25):0.0)root;'))
        counts = [1, 1, 1]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(DuplicateNodeError, _validate_otu_ids_and_tree,
                          counts, otu_ids, t)

        # unrooted tree as input
        t = TreeNode.read(io.StringIO('((OTU1:0.1, OTU2:0.2):0.3, OTU3:0.5,'
                                      'OTU4:0.7);'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # otu_ids has duplicated ids
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU2']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # len of vectors not equal
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
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
            io.StringIO('((((OTU1,OTU2),OTU3)),(OTU4,OTU5));'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # tree missing some branch lengths
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # otu_ids not present in tree
        t = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.25,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        otu_ids = ['OTU1', 'OTU2', 'OTU32']
        self.assertRaises(MissingNodeError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

        # single node tree
        t = TreeNode.read(io.StringIO('root;'))
        counts = []
        otu_ids = []
        self.assertRaises(ValueError, _validate_otu_ids_and_tree, counts,
                          otu_ids, t)

    def test_vectorize_counts_and_tree(self):
        t = TreeNode.read(io.StringIO("((a:1, b:2)c:3)root;"))
        counts = np.array([[0, 1], [1, 5], [10, 1]])
        count_array, indexed, branch_lengths = \
            _vectorize_counts_and_tree(counts, np.array(['a', 'b']), t)
        exp_counts = np.array([[0, 1, 10], [1, 5, 1], [1, 6, 11], [1, 6, 11]])
        npt.assert_equal(count_array, exp_counts.T)


if __name__ == "__main__":
    main()
