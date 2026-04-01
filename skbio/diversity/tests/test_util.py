# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import pandas as pd
import scipy.sparse as sp
import numpy.testing as npt

from skbio import TreeNode
from skbio.diversity._util import (
    vectorize_counts_and_tree,
    _validate_counts_vector,
    _validate_counts_matrix,
    _qualify_counts,
)
from skbio.tree import DuplicateNodeError, MissingNodeError


class VectorizeTests(TestCase):

    def test_vectorize_counts_and_tree(self):
        tree = TreeNode.read(["((a:1, b:2)c:3)root;"])
        counts = np.array([[0, 1], [1, 5], [10, 1]])
        count_array, indexed, branch_lengths = \
            vectorize_counts_and_tree(counts, np.array(['a', 'b']), tree)
        exp_counts = np.array([[0, 1, 10], [1, 5, 1], [1, 6, 11], [1, 6, 11]])
        npt.assert_equal(count_array, exp_counts.T)

        # test with sparse matrix
        counts_sparse = sp.csr_matrix(counts)
        count_array_sparse, _, _ = vectorize_counts_and_tree(
            counts_sparse, np.array(['a', 'b']), tree)
        npt.assert_equal(count_array_sparse, exp_counts.T)

        # test with 1-D sparse matrix
        counts_1d = np.array([1, 5])
        counts_1d_sparse = sp.csr_matrix(counts_1d)
        count_array_1d, _, _ = vectorize_counts_and_tree(
            counts_1d_sparse, np.array(['a', 'b']), tree)
        exp_counts_1d = np.array([[1, 5, 6, 6]])
        npt.assert_equal(count_array_1d, exp_counts_1d)
        # Ensure the dtype matches the input dtype when using sparse matrices
        self.assertEqual(count_array_1d.dtype, counts_1d.dtype)

        # test with subset of taxa
        counts_sub = np.array([[0], [1], [10]])
        counts_sub_sparse = sp.csr_matrix(counts_sub)
        count_array_sub, _, _ = vectorize_counts_and_tree(
            counts_sub_sparse, np.array(['a']), tree)
        exp_counts_sub = np.array([[0, 0, 0, 0], [1, 0, 1, 1], [10, 0, 10, 10]])
        npt.assert_equal(count_array_sub, exp_counts_sub)

        # test missing taxa error
        with self.assertRaisesRegex(ValueError, "Taxon 'x' not found in tree."):
            vectorize_counts_and_tree(sp.csr_matrix([[1]]), np.array(['x']), tree)


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

        # keep float
        obs = _validate_counts_vector([42.2, 42.7, 0])
        npt.assert_array_equal(obs, np.array([42.2, 42.7, 0]))
        self.assertEqual(obs.dtype, float)

        # cast into int
        obs = _validate_counts_vector([42.2, 42.7, 0], cast_int=True)
        npt.assert_array_equal(obs, np.array([42, 42, 0]))
        self.assertEqual(obs.dtype, int)

        # all zeros
        obs = _validate_counts_vector([0, 0, 0])
        npt.assert_array_equal(obs, np.array([0, 0, 0]))
        self.assertEqual(obs.dtype, int)

        # all zeros (single value)
        obs = _validate_counts_vector([0])
        npt.assert_array_equal(obs, np.array([0]))
        self.assertEqual(obs.dtype, int)

    def test_validate_counts_vector_invalid_input(self):
        # wrong data type (strings)
        with self.assertRaises(ValueError):
            _validate_counts_vector(['a', 'b', 'c'])

        # wrong data type (complex numbers)
        with self.assertRaises(ValueError):
            _validate_counts_vector([1 + 2j, 3 + 4j])

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

    def test_validate_counts_matrix_cast_int(self):
        obs = _validate_counts_matrix(
            [[42.2, 42.1, 0], [42.2, 42.1, 1.0]], cast_int=True)
        npt.assert_array_equal(obs[0], np.array([42, 42, 0]))
        npt.assert_array_equal(obs[1], np.array([42, 42, 1]))
        self.assertEqual(obs[0].dtype, int)
        self.assertEqual(obs[1].dtype, int)

    def test_validate_counts_matrix_negative_counts(self):
        with self.assertRaises(ValueError):
            _validate_counts_matrix([[0, 1, 1, 0, 2], [0, 0, 2, -1, 3]])
        with self.assertRaises(ValueError):
            _validate_counts_matrix([[0, 0, 2, -1, 3], [0, 1, 1, 0, 2]])

    def test_validate_counts_matrix_unequal_lengths(self):
        with self.assertRaises(ValueError):
            _validate_counts_matrix([[0], [0, 0], [9, 8]])
        with self.assertRaises(ValueError):
            _validate_counts_matrix([[0, 0], [0, 0, 8], [9, 8]])
        with self.assertRaises(ValueError):
            _validate_counts_matrix([[0, 0, 75], [0, 0, 3], [9, 8, 22, 44]])

    def test_validate_counts_matrix_invalid_input(self):
        with self.assertRaises(ValueError):
            _validate_counts_matrix([['a', 'b', 'c']])
        with self.assertRaises(ValueError):
            _validate_counts_matrix([[1 + 2j, 3 + 4j]])

    def test_qualify_counts(self):
        counts = np.array([[0, 1], [1, 5], [10, 1]])
        exp = np.array([[False, True], [True, True], [True, True]])
        obs = _qualify_counts(counts)
        npt.assert_equal(obs, exp)

        counts = np.array([[0, 0, 0], [1, 0, 42]])
        exp = np.array([[False, False, False], [True, False, True]])
        obs = _qualify_counts(counts)
        npt.assert_equal(obs, exp)


if __name__ == "__main__":
    main()
