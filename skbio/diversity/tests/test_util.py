# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

import numpy as np
import pandas as pd
import numpy.testing as npt

from skbio import TreeNode
from skbio.table import example_table
from skbio.diversity._util import (_validate_counts_vector,
                                   _validate_counts_matrix,
                                   _validate_taxa_and_tree,
                                   _vectorize_counts_and_tree,
                                   _quantitative_to_qualitative_counts,
                                   _check_taxa_alias,
                                   _table_to_numpy,
                                   _validate_table)
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

    def test_validate_counts_matrix_unmatching_ids(self):
        with self.assertRaises(ValueError):
            _validate_counts_matrix([[0, 1, 1, 0, 2],
                                     [0, 0, 2, 1, 3],
                                     [1, 1, 1, 1, 1]], ids=['a', 'b'])
        with self.assertRaises(ValueError):
            obs = _validate_counts_matrix(pd.DataFrame(
                [[0, 1, 1, 0, 2],
                 [0, 0, 2, 1, 3],
                 [1, 1, 1, 1, 1]]), ids=['a', 'b'])

    def test_validate_counts_matrix_unequal_lengths(self):
        # len of vectors not equal
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

    def test_validate_taxa_and_tree(self):
        # basic valid input
        tree = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 1, 1]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertTrue(_validate_taxa_and_tree(counts, taxa, tree) is None)

        # all tips observed
        tree = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 1, 1, 1, 1]
        taxa = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5']
        self.assertTrue(_validate_taxa_and_tree(counts, taxa, tree) is None)

        # no tips observed
        tree = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = []
        taxa = []
        self.assertTrue(_validate_taxa_and_tree(counts, taxa, tree) is None)

        # all counts zero
        tree = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [0, 0, 0, 0, 0]
        taxa = ['OTU1', 'OTU2', 'OTU3', 'OTU4', 'OTU5']
        self.assertTrue(_validate_taxa_and_tree(counts, taxa, tree) is None)

    def test_validate_taxa_and_tree_invalid_input(self):
        # tree has duplicated tip ids
        tree = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU2:0.75):1.25):0.0)root;'))
        counts = [1, 1, 1]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(DuplicateNodeError, _validate_taxa_and_tree,
                          counts, taxa, tree)

        # unrooted tree as input
        tree = TreeNode.read(io.StringIO('((OTU1:0.1, OTU2:0.2):0.3, OTU3:0.5,'
                                      'OTU4:0.7);'))
        counts = [1, 2, 3]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_taxa_and_tree, counts, taxa, tree)

        # taxa has duplicated ids
        tree = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        taxa = ['OTU1', 'OTU2', 'OTU2']
        self.assertRaises(ValueError, _validate_taxa_and_tree, counts, taxa, tree)

        # len of vectors not equal
        tree = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_taxa_and_tree, counts, taxa, tree)
        counts = [1, 2, 3]
        taxa = ['OTU1', 'OTU2']
        self.assertRaises(ValueError, _validate_taxa_and_tree, counts, taxa, tree)

        # tree with no branch lengths
        tree = TreeNode.read(io.StringIO('((((OTU1,OTU2),OTU3)),(OTU4,OTU5));'))
        counts = [1, 2, 3]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_taxa_and_tree, counts, taxa, tree)

        # tree missing some branch lengths
        tree = TreeNode.read(
            io.StringIO(
                '(((((OTU1,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        taxa = ['OTU1', 'OTU2', 'OTU3']
        self.assertRaises(ValueError, _validate_taxa_and_tree, counts, taxa, tree)

        # taxa not present in tree
        tree = TreeNode.read(
            io.StringIO(
                '(((((OTU1:0.25,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:'
                '0.75,OTU5:0.75):1.25):0.0)root;'))
        counts = [1, 2, 3]
        taxa = ['OTU1', 'OTU2', 'OTU32']
        self.assertRaises(MissingNodeError, _validate_taxa_and_tree, counts, taxa, tree)

        # single node tree
        tree = TreeNode.read(io.StringIO('root;'))
        counts = []
        taxa = []
        self.assertRaises(ValueError, _validate_taxa_and_tree, counts, taxa, tree)

    def test_vectorize_counts_and_tree(self):
        tree = TreeNode.read(io.StringIO("((a:1, b:2)c:3)root;"))
        counts = np.array([[0, 1], [1, 5], [10, 1]])
        count_array, indexed, branch_lengths = \
            _vectorize_counts_and_tree(counts, np.array(['a', 'b']), tree)
        exp_counts = np.array([[0, 1, 10], [1, 5, 1], [1, 6, 11], [1, 6, 11]])
        npt.assert_equal(count_array, exp_counts.T)

    def test_quantitative_to_qualitative_counts(self):
        counts = np.array([[0, 1], [1, 5], [10, 1]])
        exp = np.array([[False, True], [True, True], [True, True]])
        obs = _quantitative_to_qualitative_counts(counts)
        npt.assert_equal(obs, exp)

        counts = np.array([[0, 0, 0], [1, 0, 42]])
        exp = np.array([[False, False, False], [True, False, True]])
        obs = _quantitative_to_qualitative_counts(counts)
        npt.assert_equal(obs, exp)

    def test_check_taxa_alias(self):
        # for backward compatibility; will be removed in the future
        msg = "A list of taxon IDs must be provided."
        with self.assertRaises(ValueError) as cm:
            _check_taxa_alias(None, None, None)
        self.assertEqual(str(cm.exception), msg)

        msg = "A phylogenetic tree must be provided."
        with self.assertRaises(ValueError) as cm:
            _check_taxa_alias([1], None, None)
        self.assertEqual(str(cm.exception), msg)

        obs = _check_taxa_alias([1], '1', None)
        self.assertListEqual(obs, [1])
        obs = _check_taxa_alias(None, '1', [1])
        self.assertListEqual(obs, [1])


class TableConversionTests(TestCase):
    def test_table_to_numpy(self):
        exp_data = np.array([[0, 1, 2], [3, 4, 5]]).T
        exp_ids = ['S1', 'S2', 'S3']
        exp_feat_ids = ['O1', 'O2']
        obs_data, obs_ids, obs_feat_ids = _table_to_numpy(example_table)
        npt.assert_equal(obs_data, exp_data)
        self.assertEqual(obs_ids, exp_ids)
        self.assertEqual(obs_feat_ids, exp_feat_ids)

    def test_validate_table(self):
        self.assertRaises(ValueError, _validate_table, example_table, ['foo', 'bar'], {})
        self.assertRaises(ValueError, _validate_table, example_table, None,
                          {'taxa': 'foo'})
        obs_data, obs_ids = _validate_table(example_table, None, {})
        exp_data = np.array([[0, 1, 2], [3, 4, 5]]).T
        exp_ids = ['S1', 'S2', 'S3']
        npt.assert_equal(obs_data, exp_data)
        self.assertEqual(obs_ids, exp_ids)


if __name__ == "__main__":
    main()
