# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd

from skbio.tree import TreeNode
from skbio.stats.composition import closure
from skbio.table._augment import (
    mixup,
    aitchison_mixup,
    compos_cutmix,
    phylomix,
    _validate_labels,
    _normalize_matrix,
    _all_pairs,
    _intra_class_pairs,
    _aitchison_add,
    _aitchison_multiply,
    _indices_under_nodes,
)


class AugmentationTests(TestCase):

    def setUp(self):
        # a normal example (10 x 10)
        self.matrix = np.array([
            [ 2,  0,  3,  0,  9, 31,  3,  1,  3,  7],
            [ 9,  0,  8,  4,  1, 24,  0,  0,  3,  6],
            [ 0,  0,  1,  0,  0,  1,  2,  1,  0,  1],
            [ 0,  5,  0,  1,  7,  0,  5,  2,  4,  0],
            [ 0,  0,  5, 22,  1, 10,  4, 16, 19,  0],
            [ 0,  0,  2,  0,  1,  2,  3,  2,  0,  1],
            [ 0,  0,  0,  2,  0,  3,  0,  3,  1,  0],
            [ 5,  2,  0, 14,  0,  3,  3,  5,  9, 51],
            [15,  3,  0,  8, 17,  0,  0,  5,  0,  0],
            [ 0,  3,  5,  0,  2,  2,  2,  0,  4,  1],
        ])
        self.labels_1 = np.zeros(10, dtype=int)  # one class
        self.labels_2 = np.array([0, 0, 0, 1, 0, 1, 1, 1, 0, 0])  # two classes
        self.labels_3 = np.array([1, 1, 0, 2, 0, 0, 0, 2, 1, 2])  # three classes

        self.taxa = [f"O{i}" for i in range(10)]
        self.tree = TreeNode.read([
            "((O6,((O1,O4),O7),((O5,O2),(O3,O8))),(O0,O9));"
        ])

    def test_validate_labels(self):
        # Create a simple matrix for testing
        matrix = np.array([[1, 2, 3], [4, 5, 6]])
        n = matrix.shape[0]

        # Test with valid 1D labels
        labels = np.array([0, 1])
        obs_idx, obs_1hot = _validate_labels(labels, n)
        npt.assert_array_equal(obs_idx, labels)
        npt.assert_array_equal(obs_1hot, np.array([[1, 0], [0, 1]]))

        # Test with a plain list
        obs_idx, obs_1hot = _validate_labels(labels.tolist(), n)
        npt.assert_array_equal(obs_idx, labels)
        npt.assert_array_equal(obs_1hot, np.array([[1, 0], [0, 1]]))

        # Test with integer-equivalent floats
        obs_idx, obs_1hot = _validate_labels(labels.astype(float), n)
        npt.assert_array_equal(obs_idx, labels)
        npt.assert_array_equal(obs_1hot, np.array([[1, 0], [0, 1]]))
        self.assertTrue(np.issubdtype(obs_idx.dtype, np.integer))

        # Test with valid 2D one-hot encoded labels (discrete)
        one_hot = np.array([[1, 0], [0, 1]])
        obs_idx, obs_1hot = _validate_labels(one_hot, n)
        npt.assert_array_equal(obs_idx, np.array([0, 1]))
        npt.assert_array_equal(obs_1hot, one_hot)

        # Test with non-integer but still valid one-hot labels
        one_hot = np.array([[0.6, 0.4], [0.2, 0.8]])
        obs_idx, obs_1hot = _validate_labels(one_hot, n)
        npt.assert_array_equal(obs_idx, np.array([0, 1]))
        npt.assert_array_equal(obs_1hot, one_hot)

    def test_validate_labels_errors(self):
        matrix = np.array([[1, 2, 3], [4, 5, 6]])
        n = matrix.shape[0]

        # Test with wrong dimensions
        msg = "Labels should be 1-D or 2-D, but got 3 dimensions instead."
        wrong_dim_labels = np.array([[[0], [1]]])
        with self.assertRaises(ValueError) as cm:
            _validate_labels(wrong_dim_labels, n)
        self.assertEqual(str(cm.exception), msg)

        # Test with wrong number of samples
        msg = "Number of labels (3) does not match number of samples in the data (2)."
        wrong_len_labels = np.array([0, 1, 2])
        with self.assertRaises(ValueError) as cm:
            _validate_labels(wrong_len_labels, n)
        self.assertEqual(str(cm.exception), msg)

        # Test with non-whole numbers
        msg = "Labels must only contain integer values."
        non_int_labels = np.array([1.5, 2.0])
        with self.assertRaises(ValueError) as cm:
            _validate_labels(non_int_labels, n)
        self.assertEqual(str(cm.exception), msg)

        # Test with non-zero indexed label
        msg = "Labels must be zero-indexed. Minimum value must be 0."
        non_zero_labels = np.array([1, 2])
        with self.assertRaises(ValueError) as cm:
            _validate_labels(non_zero_labels, n)
        self.assertEqual(str(cm.exception), msg)

        # Test with non-consecutive label
        msg = "Labels must be consecutive integers from 0 to n_classes - 1."
        non_consecutive = np.array([0, 2])
        with self.assertRaises(ValueError) as cm:
            _validate_labels(non_consecutive, n)
        self.assertEqual(str(cm.exception), msg)

        # Test with invalid one-hot encoding
        msg = "Labels are not properly one-hot encoded."
        invalid_one_hot = np.array([[1, 1], [0, 1]])
        with self.assertRaises(ValueError) as cm:
            _validate_labels(invalid_one_hot, n)
        self.assertEqual(str(cm.exception), msg)

    def test_normalize_matrix(self):
        matrix = np.array([[0.2, 0.5, 0.3], [0.3, 0.6, 0.1]])
        obs = _normalize_matrix(matrix)
        self.assertIs(obs, matrix)

        matrix = np.array([[2, 5, 3], [3, 6, 1]], dtype=float)
        obs = _normalize_matrix(matrix)
        exp = np.array([[0.2, 0.5, 0.3], [0.3, 0.6, 0.1]])
        npt.assert_allclose(obs, exp)

    def test_aitchison_add(self):
        """Test Aitchison addition."""
        # vector addition
        x1 = np.array([0.2, 0.3, 0.5])
        x2 = np.array([0.1, 0.4, 0.5])
        obs = _aitchison_add(x1, x2)
        exp = np.array([0.05128, 0.30769, 0.64103])
        npt.assert_array_equal(obs.round(5), exp)
        self.assertAlmostEqual(obs.sum(), 1.0)

        # matrix addition
        x1 = np.array([[0.2, 0.3, 0.5],
                       [0.1, 0.4, 0.5]])
        x2 = np.array([[0.3, 0.3, 0.4],
                       [0.5, 0.2, 0.3]])
        obs = _aitchison_add(x1, x2)
        exp = np. array([[0.17143, 0.25714, 0.57143],
                         [0.17857, 0.28571, 0.53571]])
        npt.assert_array_equal(obs.round(5), exp)
        self.assertTrue(np.allclose(obs.sum(axis=1), 1))

    def test_aitchison_multiply(self):
        """Test Aitchison scalar multiplication."""
        p = 2.0

        # vector multiplication
        x = np.array([0.2, 0.3, 0.5])
        obs = _aitchison_multiply(x, p)
        exp = np.array([0.10526, 0.23684, 0.65789])
        npt.assert_array_equal(obs.round(5), exp)
        self.assertAlmostEqual(obs.sum(), 1.0)

        # matrix multiplication
        x = np.array([[0.2, 0.3, 0.5],
                      [0.1, 0.4, 0.5]])
        obs = _aitchison_multiply(x, p)
        exp = np.array([[0.10526, 0.23684, 0.65789],
                        [0.02381, 0.38095, 0.59524]])
        npt.assert_array_equal(obs.round(5), exp)
        self.assertTrue(np.allclose(obs.sum(axis=1), 1))

    def test_all_pairs(self):
        obs = _all_pairs(2)
        exp = np.array([[0, 1]])
        npt.assert_array_equal(obs, exp)

        obs = _all_pairs(3)
        exp = np.array([[0, 1], [0, 2], [1, 2]])
        npt.assert_array_equal(obs, exp)

        obs = _all_pairs(5)
        exp = np.array([[0, 0, 0, 0, 1, 1, 1, 2, 2, 3],
                        [1, 2, 3, 4, 2, 3, 4, 3, 4, 4]]).T
        npt.assert_array_equal(obs, exp)

        # edge case: one sample -> no pair
        obs = _all_pairs(1)
        exp = np.empty((0, 2), dtype=int)
        npt.assert_array_equal(obs, exp)

    def test_intra_class_pairs(self):
        # one class only
        obs = _intra_class_pairs([0, 0, 0])
        exp = np.array([[0, 1], [0, 2], [1, 2]])
        npt.assert_array_equal(obs, exp)

        # two classes, one pair each
        obs = _intra_class_pairs([0, 0, 1, 1])
        exp = np.array([[0, 1], [2, 3]])
        npt.assert_array_equal(obs, exp)

        # two classes, interleaved
        obs = _intra_class_pairs([0, 1, 0, 1])
        exp = np.array([[0, 2], [1, 3]])
        npt.assert_array_equal(obs, exp)

        # two classes, three pairs each
        obs = _intra_class_pairs([0, 0, 1, 1, 1, 0])
        exp = np.array([[0, 0, 1, 2, 2, 3],
                        [1, 5, 5, 3, 4, 4]]).T
        npt.assert_array_equal(obs, exp)

        # three classes, mixed pair numbers
        obs = _intra_class_pairs([0, 0, 1, 1, 2, 0, 2, 1, 1, 1])
        exp = np.array([[0, 0, 1, 2, 2, 2, 2, 3, 3, 3, 7, 7, 8, 4],
                        [1, 5, 5, 3, 7, 8, 9, 7, 8, 9, 8, 9, 9, 6]]).T
        npt.assert_array_equal(obs, exp)

        # edge case: no pair
        obs = _intra_class_pairs([0, 1, 2])
        exp = np.empty((0, 2), dtype=int)
        npt.assert_array_equal(obs, exp)

        # edge case: no sample
        obs = _intra_class_pairs([])
        exp = np.empty((0, 2), dtype=int)
        npt.assert_array_equal(obs, exp)

    def test_indices_under_nodes(self):
        # a tree with 10 taxa
        tree = self.tree
        taxa = self.taxa
        obs = _indices_under_nodes(tree, taxa)
        # naive method (not optimized)
        taxon_map = {taxon: i for i, taxon in enumerate(taxa)}
        nodes = tree.non_tips(include_self=True)
        exp = [[taxon_map[tip.name] for tip in node.tips()] for node in nodes]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertListEqual(list(o.keys()), e)

        # taxon map doesn't cover entire tree (a common scenario)
        taxa = ["O3", "O7", "O1", "O4", "O9"]
        obs = _indices_under_nodes(tree, taxa)
        # manually constructed
        exp = [[2, 3], [2, 3, 1], [0], [0], [2, 3, 1, 0], [4], [2, 3, 1, 0, 4]]
        self.assertEqual(len(obs), len(exp))
        for o, e in zip(obs, exp):
            self.assertListEqual(list(o.keys()), e)

        # check if temporary attributes have been cleaned up
        for node in tree.traverse(include_self=True):
            self.assertFalse(hasattr(node, "_taxa"))

        # taxa contain duplicates
        msg = "All taxa must be unique."
        taxa = ["O1", "O2", "O1"]
        with self.assertRaises(ValueError) as cm:
            _indices_under_nodes(tree, taxa)
        self.assertEqual(str(cm.exception), msg)

        # taxa not found in the tree
        msg = "2 taxa are not present as tip names in the tree."
        taxa = ["O1", "X2", "O3", "X4", "O5"]
        with self.assertRaises(ValueError) as cm:
            _indices_under_nodes(tree, taxa)
        self.assertEqual(str(cm.exception), msg)

    def test_mixup(self):
        matrix, labels_2, labels_3 = self.matrix, self.labels_2, self.labels_3

        # no labels
        obs_mat, obs_lab = mixup(matrix, n=10, seed=42)
        exp_mat = np.array([
            [ 0.55,  0.  ,  2.27,  0.  ,  3.19,  9.95,  3.  ,  1.73,  0.82,  2.64],
            [ 0.  ,  2.36,  5.  ,  4.73,  1.79,  3.72,  2.43,  3.44,  7.22,  0.79],
            [ 0.  ,  4.64,  0.9 ,  0.82,  6.1 ,  0.36,  4.46,  1.64,  4.  ,  0.18],
            [ 0.  ,  0.  ,  1.69,  0.  ,  0.69,  1.69,  2.69,  1.69,  0.  ,  1.  ],
            [ 0.  ,  0.  ,  1.61,  0.  ,  0.61,  1.61,  2.61,  1.61,  0.  ,  1.  ],
            [ 0.  ,  1.91,  3.91,  0.  ,  1.64,  2.  ,  2.36,  0.73,  2.54,  1.  ],
            [ 1.31,  0.  ,  3.69,  7.55,  6.25, 23.79,  3.34,  6.15,  8.49,  4.6 ],
            [ 0.  ,  0.  ,  0.93,  5.72,  0.19,  4.3 ,  0.74,  5.42,  4.35,  0.  ],
            [ 2.08,  0.  ,  2.62,  0.93,  0.23,  6.33,  1.54,  0.77,  0.69,  2.16],
            [ 1.11,  0.  ,  2.55,  0.  ,  5.44, 18.08,  3.  ,  1.45,  1.66,  4.33],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)
        self.assertIsNone(obs_lab)

        # append mode
        obs_mat, obs_lab = mixup(matrix, n=10, seed=42, append=True)
        self.assertTupleEqual(obs_mat.shape, (20, 10))
        npt.assert_array_equal(obs_mat[:10], matrix)
        npt.assert_array_equal(obs_mat.round(2)[10:], exp_mat)
        self.assertIsNone(obs_lab)

        # two-class labels
        obs_mat, obs_lab = mixup(matrix, n=10, labels=labels_2, seed=42)
        npt.assert_array_equal(obs_mat.round(2), exp_mat)
        exp_lab = np.array([
            [0.27, 0.73],
            [1.  , 0.  ],
            [0.18, 0.82],
            [0.31, 0.69],
            [0.39, 0.61],
            [0.64, 0.36],
            [1.  , 0.  ],
            [0.19, 0.81],
            [1.  , 0.  ],
            [0.55, 0.45],
        ])
        npt.assert_array_equal(obs_lab.round(2), exp_lab)

        # append mode
        obs_mat, obs_lab = mixup(matrix, n=10, labels=labels_2, seed=42, append=True)
        self.assertTupleEqual(obs_mat.shape, (20, 10))
        npt.assert_array_equal(obs_mat[:10], matrix)
        npt.assert_array_equal(obs_mat.round(2)[10:], exp_mat)
        self.assertTupleEqual(obs_lab.shape, (20, 2))
        npt.assert_array_equal(obs_lab[:10], np.eye(2, dtype=int)[labels_2])
        npt.assert_array_equal(obs_lab.round(2)[10:], exp_lab)

        # alpha parameter
        obs_mat, obs_lab = mixup(matrix, n=10, alpha=1, seed=42)
        exp_mat = np.array([
            [ 1.72,  0.  ,  2.86,  0.  ,  7.88, 26.94,  3.  ,  1.14,  2.58,  6.16],
            [ 0.  ,  1.35,  5.  , 12.07,  1.45,  6.39,  3.1 ,  8.78, 12.23,  0.45],
            [ 0.  ,  3.58,  3.55,  0.29,  3.45,  1.42,  2.87,  0.58,  4.  ,  0.71],
            [ 0.  ,  0.  ,  1.93,  0.  ,  0.93,  1.93,  2.93,  1.93,  0.  ,  1.  ],
            [ 0.  ,  0.  ,  1.71,  0.  ,  0.71,  1.71,  2.71,  1.71,  0.  ,  1.  ],
            [ 0.  ,  2.34,  4.34,  0.  ,  1.78,  2.  ,  2.22,  0.44,  3.12,  1.  ],
            [ 0.88,  0.  ,  4.12, 12.3 ,  4.53, 19.26,  3.56,  9.38, 11.94,  3.09],
            [ 0.  ,  0.  ,  2.97, 13.86,  0.59,  7.15,  2.37, 10.71, 11.68,  0.  ],
            [ 6.09,  0.  ,  5.74,  2.71,  0.68, 16.57,  0.65,  0.32,  2.03,  4.39],
            [ 0.59,  0.  ,  2.3 ,  0.  ,  3.38, 10.61,  3.  ,  1.7 ,  0.89,  2.78],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)

        # two-class labels, intra-class mixing
        obs_mat, obs_lab = mixup(
            matrix, n=10, labels=labels_2, intra_class=True, seed=42
        )
        exp_mat = np.array([
            [ 0.55,  0.  ,  1.55,  0.  ,  2.47,  9.22,  2.27,  1.  ,  0.82,  2.64],
            [ 0.  ,  1.07,  0.  ,  1.79,  1.5 ,  2.36,  1.07,  2.79,  1.64,  0.  ],
            [ 0.  ,  0.54,  5.  , 18.05,  1.18,  8.56,  3.64, 13.13, 16.31,  0.18],
            [ 0.  ,  0.  ,  3.77, 15.25,  0.69,  7.24,  3.39, 11.4 , 13.17,  0.31],
            [ 0.  ,  0.  ,  3.45, 13.48,  0.61,  6.51,  3.23, 10.19, 11.64,  0.39],
            [ 0.  ,  0.  ,  0.73,  1.27,  0.36,  2.64,  1.09,  2.64,  0.64,  0.36],
            [ 1.31,  0.  ,  2.31,  0.  ,  5.91, 20.7 ,  2.66,  1.  ,  1.97,  4.94],
            [ 2.79,  3.  ,  4.07,  1.49,  4.79,  1.63,  1.63,  0.93,  3.26,  0.81],
            [ 0.46,  2.31,  4.54,  0.  ,  3.62,  8.72,  2.23,  0.23,  3.77,  2.39],
            [ 1.11,  0.  ,  2.11,  0.  ,  4.99, 17.64,  2.55,  1.  ,  1.66,  4.33],
        ])
        exp_lab = np.array([
            [1, 0],
            [0, 1],
            [1, 0],
            [1, 0],
            [1, 0],
            [0, 1],
            [1, 0],
            [1, 0],
            [1, 0],
            [1, 0],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)
        npt.assert_array_equal(obs_lab, exp_lab)

        # three-class labels
        obs_mat, obs_lab = mixup(matrix, n=10, labels=labels_3, seed=42)
        exp_mat = np.array([
            [ 0.55,  0.  ,  2.27,  0.  ,  3.19,  9.95,  3.  ,  1.73,  0.82,  2.64],
            [ 0.  ,  2.36,  5.  ,  4.73,  1.79,  3.72,  2.43,  3.44,  7.22,  0.79],
            [ 0.  ,  4.64,  0.9 ,  0.82,  6.1 ,  0.36,  4.46,  1.64,  4.  ,  0.18],
            [ 0.  ,  0.  ,  1.69,  0.  ,  0.69,  1.69,  2.69,  1.69,  0.  ,  1.  ],
            [ 0.  ,  0.  ,  1.61,  0.  ,  0.61,  1.61,  2.61,  1.61,  0.  ,  1.  ],
            [ 0.  ,  1.91,  3.91,  0.  ,  1.64,  2.  ,  2.36,  0.73,  2.54,  1.  ],
            [ 1.31,  0.  ,  3.69,  7.55,  6.25, 23.79,  3.34,  6.15,  8.49,  4.6 ],
            [ 0.  ,  0.  ,  0.93,  5.72,  0.19,  4.3 ,  0.74,  5.42,  4.35,  0.  ],
            [ 2.08,  0.  ,  2.62,  0.93,  0.23,  6.33,  1.54,  0.77,  0.69,  2.16],
            [ 1.11,  0.  ,  2.55,  0.  ,  5.44, 18.08,  3.  ,  1.45,  1.66,  4.33],
        ])
        exp_lab = np.array([
            [0.73, 0.27, 0.  ],
            [0.21, 0.  , 0.79],
            [0.  , 0.  , 1.  ],
            [1.  , 0.  , 0.  ],
            [1.  , 0.  , 0.  ],
            [0.36, 0.  , 0.64],
            [0.34, 0.66, 0.  ],
            [1.  , 0.  , 0.  ],
            [0.77, 0.23, 0.  ],
            [0.45, 0.55, 0.  ],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)
        npt.assert_array_equal(obs_lab.round(2), exp_lab)

    def test_mixup_errors(self):
        msg = "Cannot find a pair of samples to mix."

        # Each class has only one sample. In intra-class mode there is no pair to mix.
        labels = np.arange(10)
        with self.assertRaises(ValueError) as cm:
            mixup(self.matrix, n=10, labels=labels, intra_class=True, seed=42)
        self.assertEqual(str(cm.exception), msg)

        # Entire dataset has only one sample => no pair to mix.
        with self.assertRaises(ValueError) as cm:
            mixup(self.matrix[:1], n=10, seed=42)
        self.assertEqual(str(cm.exception), msg)

    def test_aitchison_mixup(self):
        matrix, labels_2, labels_3 = self.matrix, self.labels_2, self.labels_3

        # no labels
        obs_mat, obs_lab = aitchison_mixup(matrix, n=10, seed=42)
        exp_mat = np.array([
            [0.  , 0.  , 0.15, 0.  , 0.12, 0.29, 0.2 , 0.11, 0.  , 0.12],
            [0.  , 0.  , 0.29, 0.  , 0.1 , 0.16, 0.13, 0.  , 0.32, 0.  ],
            [0.  , 0.25, 0.  , 0.  , 0.3 , 0.  , 0.23, 0.  , 0.22, 0.  ],
            [0.  , 0.  , 0.19, 0.  , 0.  , 0.19, 0.31, 0.19, 0.  , 0.12],
            [0.  , 0.  , 0.19, 0.  , 0.  , 0.19, 0.31, 0.19, 0.  , 0.12],
            [0.  , 0.  , 0.34, 0.  , 0.15, 0.19, 0.22, 0.  , 0.  , 0.1 ],
            [0.  , 0.  , 0.09, 0.  , 0.1 , 0.52, 0.08, 0.06, 0.14, 0.  ],
            [0.  , 0.  , 0.  , 0.25, 0.  , 0.3 , 0.  , 0.32, 0.14, 0.  ],
            [0.  , 0.  , 0.31, 0.  , 0.  , 0.4 , 0.  , 0.  , 0.  , 0.29],
            [0.  , 0.  , 0.11, 0.  , 0.15, 0.41, 0.13, 0.06, 0.  , 0.13],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)
        self.assertIsNone(obs_lab)

        # check if the synthetic samples are also compositional
        self.assertTrue(np.allclose(obs_mat.sum(axis=1), 1))

        # append mode
        obs_mat, obs_lab = aitchison_mixup(matrix, n=10, seed=42, append=True)
        self.assertTupleEqual(obs_mat.shape, (20, 10))
        npt.assert_array_equal(obs_mat.round(2)[10:], exp_mat)

        # check if the original samples have been normalized
        npt.assert_array_equal(obs_mat[:10].round(2), closure(matrix).round(2))
        # feed pre-normalized data
        normed = closure(matrix)
        obs_mat, obs_lab = aitchison_mixup(normed, n=10, seed=42, normalize=False)
        npt.assert_array_equal(obs_mat.round(2), exp_mat)

        # skip normalization
        obs_mat, obs_lab = aitchison_mixup(matrix, n=10, seed=42, normalize=False)
        exp_mat = np.array([
            [0.  , 0.  , 0.15, 0.  , 0.12, 0.29, 0.2 , 0.11, 0.  , 0.12],
            [0.  , 0.  , 0.29, 0.  , 0.1 , 0.16, 0.13, 0.  , 0.32, 0.  ],
            [0.  , 0.25, 0.  , 0.  , 0.3 , 0.  , 0.23, 0.  , 0.22, 0.  ],
            [0.  , 0.  , 0.19, 0.  , 0.  , 0.19, 0.31, 0.19, 0.  , 0.12],
            [0.  , 0.  , 0.19, 0.  , 0.  , 0.19, 0.31, 0.19, 0.  , 0.12],
            [0.  , 0.  , 0.34, 0.  , 0.15, 0.19, 0.22, 0.  , 0.  , 0.1 ],
            [0.  , 0.  , 0.09, 0.  , 0.1 , 0.52, 0.08, 0.06, 0.14, 0.  ],
            [0.  , 0.  , 0.  , 0.25, 0.  , 0.3 , 0.  , 0.32, 0.14, 0.  ],
            [0.  , 0.  , 0.31, 0.  , 0.  , 0.4 , 0.  , 0.  , 0.  , 0.29],
            [0.  , 0.  , 0.11, 0.  , 0.15, 0.41, 0.13, 0.06, 0.  , 0.13],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)

        # two-class labels
        obs_mat, obs_lab = aitchison_mixup(matrix, n=10, labels=labels_2, seed=42)
        npt.assert_array_equal(obs_mat.round(2), exp_mat)
        # Note that the synthetic labels are identical to those of `mixup`. This is
        # because they are determined by sampling from beta distribution, which is
        # independent from the subsequent calculation.
        exp_lab = np.array([
            [0.27, 0.73],
            [1.  , 0.  ],
            [0.18, 0.82],
            [0.31, 0.69],
            [0.39, 0.61],
            [0.64, 0.36],
            [1.  , 0.  ],
            [0.19, 0.81],
            [1.  , 0.  ],
            [0.55, 0.45],
        ])
        npt.assert_array_equal(obs_lab.round(2), exp_lab)

        # alpha parameter
        obs_mat, obs_lab = aitchison_mixup(matrix, n=10, alpha=1, seed=42)
        exp_mat = np.array([
            [0.  , 0.  , 0.07, 0.  , 0.17, 0.53, 0.08, 0.03, 0.  , 0.13],
            [0.  , 0.  , 0.21, 0.  , 0.06, 0.21, 0.12, 0.  , 0.4 , 0.  ],
            [0.  , 0.27, 0.  , 0.  , 0.22, 0.  , 0.2 , 0.  , 0.31, 0.  ],
            [0.  , 0.  , 0.2 , 0.  , 0.  , 0.2 , 0.3 , 0.2 , 0.  , 0.1 ],
            [0.  , 0.  , 0.19, 0.  , 0.  , 0.19, 0.31, 0.19, 0.  , 0.12],
            [0.  , 0.  , 0.37, 0.  , 0.16, 0.18, 0.2 , 0.  , 0.  , 0.09],
            [0.  , 0.  , 0.1 , 0.  , 0.07, 0.41, 0.09, 0.12, 0.21, 0.  ],
            [0.  , 0.  , 0.  , 0.29, 0.  , 0.22, 0.  , 0.29, 0.2 , 0.  ],
            [0.  , 0.  , 0.25, 0.  , 0.  , 0.54, 0.  , 0.  , 0.  , 0.21],
            [0.  , 0.  , 0.15, 0.  , 0.13, 0.3 , 0.2 , 0.11, 0.  , 0.12],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)
        self.assertTrue(np.allclose(obs_mat.sum(axis=1), 1))

        # two-class labels, intra-class mixing
        obs_mat, obs_lab = aitchison_mixup(
            matrix, n=10, labels=labels_2, intra_class=True, seed=42
        )
        exp_mat = np.array([
            [0.  , 0.  , 0.15, 0.  , 0.  , 0.29, 0.25, 0.11, 0.  , 0.19],
            [0.  , 0.  , 0.  , 0.3 , 0.  , 0.  , 0.  , 0.47, 0.23, 0.  ],
            [0.  , 0.  , 0.16, 0.  , 0.04, 0.24, 0.11, 0.  , 0.46, 0.  ],
            [0.  , 0.  , 0.17, 0.  , 0.  , 0.27, 0.18, 0.38, 0.  , 0.  ],
            [0.  , 0.  , 0.18, 0.  , 0.  , 0.27, 0.2 , 0.36, 0.  , 0.  ],
            [0.  , 0.  , 0.  , 0.  , 0.  , 0.5 , 0.  , 0.5 , 0.  , 0.  ],
            [0.  , 0.  , 0.11, 0.  , 0.  , 0.51, 0.14, 0.05, 0.  , 0.19],
            [0.  , 0.5 , 0.  , 0.  , 0.5 , 0.  , 0.  , 0.  , 0.  , 0.  ],
            [0.  , 0.  , 0.24, 0.  , 0.15, 0.2 , 0.12, 0.  , 0.2 , 0.08],
            [0.  , 0.  , 0.12, 0.  , 0.  , 0.45, 0.17, 0.07, 0.  , 0.2 ],
        ])
        exp_lab = np.array([
            [1, 0],
            [0, 1],
            [1, 0],
            [1, 0],
            [1, 0],
            [0, 1],
            [1, 0],
            [1, 0],
            [1, 0],
            [1, 0],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)
        npt.assert_array_equal(obs_lab, exp_lab)
        self.assertTrue(np.allclose(obs_mat.sum(axis=1), 1))

        # three-class labels
        obs_mat, obs_lab = aitchison_mixup(matrix, n=10, labels=labels_3, seed=42)
        exp_mat = np.array([
            [0.  , 0.  , 0.15, 0.  , 0.12, 0.29, 0.2 , 0.11, 0.  , 0.12],
            [0.  , 0.  , 0.29, 0.  , 0.1 , 0.16, 0.13, 0.  , 0.32, 0.  ],
            [0.  , 0.25, 0.  , 0.  , 0.3 , 0.  , 0.23, 0.  , 0.22, 0.  ],
            [0.  , 0.  , 0.19, 0.  , 0.  , 0.19, 0.31, 0.19, 0.  , 0.12],
            [0.  , 0.  , 0.19, 0.  , 0.  , 0.19, 0.31, 0.19, 0.  , 0.12],
            [0.  , 0.  , 0.34, 0.  , 0.15, 0.19, 0.22, 0.  , 0.  , 0.1 ],
            [0.  , 0.  , 0.09, 0.  , 0.1 , 0.52, 0.08, 0.06, 0.14, 0.  ],
            [0.  , 0.  , 0.  , 0.25, 0.  , 0.3 , 0.  , 0.32, 0.14, 0.  ],
            [0.  , 0.  , 0.31, 0.  , 0.  , 0.4 , 0.  , 0.  , 0.  , 0.29],
            [0.  , 0.  , 0.11, 0.  , 0.15, 0.41, 0.13, 0.06, 0.  , 0.13],
        ])
        exp_lab = np.array([
            [0.73, 0.27, 0.  ],
            [0.21, 0.  , 0.79],
            [0.  , 0.  , 1.  ],
            [1.  , 0.  , 0.  ],
            [1.  , 0.  , 0.  ],
            [0.36, 0.  , 0.64],
            [0.34, 0.66, 0.  ],
            [1.  , 0.  , 0.  ],
            [0.77, 0.23, 0.  ],
            [0.45, 0.55, 0.  ],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)
        npt.assert_array_equal(obs_lab.round(2), exp_lab)
        self.assertTrue(np.allclose(obs_mat.sum(axis=1), 1))

    def test_compos_cutmix(self):
        matrix, labels_2, labels_3 = self.matrix, self.labels_2, self.labels_3

        # no labels
        obs_mat, obs_lab = compos_cutmix(matrix, n=10, seed=42)
        exp_mat = np.array([
            [0.03, 0.  , 0.05, 0.  , 0.15, 0.53, 0.05, 0.02, 0.05, 0.12],
            [0.  , 0.  , 0.06, 0.27, 0.01, 0.12, 0.1 , 0.2 , 0.23, 0.  ],
            [0.  , 0.22, 0.  , 0.04, 0.31, 0.  , 0.11, 0.09, 0.18, 0.06],
            [0.  , 0.  , 0.18, 0.  , 0.09, 0.18, 0.27, 0.18, 0.  , 0.09],
            [0.  , 0.  , 0.16, 0.  , 0.09, 0.16, 0.26, 0.16, 0.  , 0.16],
            [0.  , 0.13, 0.22, 0.  , 0.09, 0.09, 0.09, 0.15, 0.18, 0.04],
            [0.03, 0.  , 0.05, 0.  , 0.15, 0.53, 0.05, 0.02, 0.05, 0.12],
            [0.  , 0.  , 0.07, 0.3 , 0.01, 0.14, 0.  , 0.22, 0.26, 0.  ],
            [0.25, 0.  , 0.22, 0.11, 0.  , 0.25, 0.  , 0.  , 0.  , 0.17],
            [0.03, 0.  , 0.16, 0.  , 0.13, 0.16, 0.24, 0.16, 0.04, 0.08],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)
        self.assertIsNone(obs_lab)

        # check if the synthetic samples are also compositional
        self.assertTrue(np.allclose(obs_mat.sum(axis=1), 1))

        # append mode
        obs_mat, obs_lab = compos_cutmix(matrix, n=10, seed=42, append=True)
        self.assertTupleEqual(obs_mat.shape, (20, 10))
        npt.assert_array_equal(obs_mat.round(2)[10:], exp_mat)

        # check if the original samples have been normalized
        npt.assert_array_equal(obs_mat[:10].round(2), closure(matrix).round(2))

        # feed pre-normalized data
        normed = closure(matrix)
        obs_mat, obs_lab = compos_cutmix(normed, n=10, seed=42, normalize=False)
        npt.assert_array_equal(obs_mat.round(2), exp_mat)

        # skip normalization
        obs_mat, obs_lab = compos_cutmix(matrix, n=10, seed=42, normalize=False)
        exp_mat = np.array([
            [0.03, 0.  , 0.05, 0.  , 0.15, 0.53, 0.05, 0.02, 0.05, 0.12],
            [0.  , 0.  , 0.07, 0.29, 0.01, 0.13, 0.03, 0.21, 0.25, 0.  ],
            [0.  , 0.23, 0.  , 0.05, 0.32, 0.  , 0.09, 0.09, 0.18, 0.05],
            [0.  , 0.  , 0.18, 0.  , 0.09, 0.18, 0.27, 0.18, 0.  , 0.09],
            [0.  , 0.  , 0.12, 0.  , 0.12, 0.12, 0.38, 0.12, 0.  , 0.12],
            [0.  , 0.14, 0.24, 0.  , 0.1 , 0.1 , 0.1 , 0.1 , 0.19, 0.05],
            [0.03, 0.  , 0.05, 0.  , 0.15, 0.53, 0.05, 0.02, 0.05, 0.12],
            [0.  , 0.  , 0.07, 0.3 , 0.01, 0.14, 0.  , 0.22, 0.26, 0.  ],
            [0.32, 0.  , 0.29, 0.14, 0.  , 0.04, 0.  , 0.  , 0.  , 0.21],
            [0.08, 0.  , 0.08, 0.  , 0.38, 0.08, 0.12, 0.08, 0.12, 0.04],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)

        # although original samples are not compositional, synthetic samples are still
        # compositional.
        self.assertTrue(np.allclose(obs_mat.sum(axis=1), 1))

        # two-class labels
        # Note that compos_cutmix is always intra-class.
        obs_mat, obs_lab = compos_cutmix(matrix, n=10, labels=labels_2, seed=42)
        exp_mat = np.array([
            [0.03, 0.  , 0.05, 0.  , 0.15, 0.53, 0.05, 0.02, 0.05, 0.12],
            [0.  , 0.26, 0.  , 0.05, 0.37, 0.  , 0.  , 0.11, 0.21, 0.  ],
            [0.  , 0.  , 0.06, 0.26, 0.01, 0.12, 0.1 , 0.19, 0.22, 0.05],
            [0.  , 0.  , 0.06, 0.29, 0.01, 0.13, 0.05, 0.21, 0.25, 0.  ],
            [0.  , 0.  , 0.16, 0.28, 0.01, 0.16, 0.05, 0.16, 0.  , 0.16],
            [0.  , 0.  , 0.  , 0.26, 0.  , 0.39, 0.  , 0.21, 0.13, 0.  ],
            [0.03, 0.  , 0.05, 0.  , 0.15, 0.53, 0.05, 0.02, 0.05, 0.12],
            [0.26, 0.13, 0.  , 0.14, 0.29, 0.  , 0.09, 0.09, 0.  , 0.  ],
            [0.05, 0.  , 0.07, 0.  , 0.15, 0.15, 0.07, 0.02, 0.3 , 0.17],
            [0.03, 0.  , 0.13, 0.  , 0.12, 0.13, 0.27, 0.13, 0.04, 0.13],
        ])
        exp_lab = np.array([
            [1, 0],
            [0, 1],
            [1, 0],
            [1, 0],
            [1, 0],
            [0, 1],
            [1, 0],
            [1, 0],
            [1, 0],
            [1, 0],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)
        npt.assert_array_equal(obs_lab, exp_lab)
        self.assertTrue(np.allclose(obs_mat.sum(axis=1), 1))

        # three-class labels
        obs_mat, obs_lab = compos_cutmix(matrix, n=10, labels=labels_3, seed=42)
        exp_mat = np.array([
            [0.  , 0.  , 0.17, 0.  , 0.  , 0.17, 0.33, 0.17, 0.  , 0.17],
            [0.  , 0.25, 0.  , 0.05, 0.35, 0.  , 0.04, 0.1 , 0.2 , 0.  ],
            [0.04, 0.  , 0.06, 0.  , 0.18, 0.63, 0.  , 0.02, 0.06, 0.  ],
            [0.  , 0.  , 0.  , 0.22, 0.  , 0.33, 0.  , 0.33, 0.11, 0.  ],
            [0.  , 0.  , 0.21, 0.26, 0.  , 0.21, 0.  , 0.21, 0.  , 0.11],
            [0.  , 0.15, 0.24, 0.  , 0.1 , 0.1 , 0.1 , 0.08, 0.19, 0.05],
            [0.  , 0.  , 0.17, 0.  , 0.  , 0.17, 0.33, 0.17, 0.  , 0.17],
            [0.15, 0.06, 0.14, 0.07, 0.02, 0.41, 0.  , 0.  , 0.05, 0.1 ],
            [0.  , 0.  , 0.13, 0.  , 0.  , 0.26, 0.26, 0.13, 0.09, 0.13],
            [0.  , 0.  , 0.2 , 0.  , 0.  , 0.2 , 0.3 , 0.2 , 0.  , 0.1 ],
        ])
        exp_lab = np.array([
            [1, 0, 0],
            [0, 0, 1],
            [0, 1, 0],
            [1, 0, 0],
            [1, 0, 0],
            [0, 0, 1],
            [1, 0, 0],
            [0, 1, 0],
            [1, 0, 0],
            [1, 0, 0],
        ])
        npt.assert_array_equal(obs_mat.round(2), exp_mat)
        npt.assert_array_equal(obs_lab, exp_lab)
        self.assertTrue(np.allclose(obs_mat.sum(axis=1), 1))

    def test_phylomix(self):
        matrix, tree, taxa, = self.matrix, self.tree, self.taxa
        labels_2, labels_3 = self.labels_2, self.labels_3

        # no labels
        obs_mat, obs_lab = phylomix(matrix, n=10, tree=tree, taxa=taxa, seed=42)
        exp_mat = np.array([
            [ 0,  0,  0,  0,  1,  4,  3,  0,  0,  1],
            [ 0, 12, 21,  0,  8,  8,  8,  0, 17,  0],
            [ 0,  3,  5,  0,  2,  2,  2,  0,  4,  1],
            [ 0,  0,  2,  0,  0,  2,  3,  2,  0,  1],
            [ 0,  0,  2,  0,  0,  2,  3,  2,  0,  2],
            [ 0,  1,  2,  0,  1,  1,  1,  2,  2,  0],
            [ 2,  0,  3,  0,  9,  7,  3, 12, 14,  7],
            [ 0,  0,  0, 17,  0, 25,  0, 25,  8,  0],
            [ 0,  0, 13,  0,  0, 13,  0, 13,  0, 13],
            [ 0,  0,  0,  0,  1,  3,  3,  2,  0,  0],
        ])
        npt.assert_array_equal(obs_mat, exp_mat)
        self.assertIsNone(obs_lab)

        # append mode
        obs_mat, obs_lab = phylomix(
            matrix, n=10, tree=tree, taxa=taxa, seed=42, append=True
        )
        self.assertTupleEqual(obs_mat.shape, (20, 10))
        npt.assert_array_equal(obs_mat[:10], matrix)
        npt.assert_array_equal(obs_mat[10:], exp_mat)

        # two-class labels
        obs_mat, obs_lab = phylomix(
            matrix, n=10, tree=tree, taxa=taxa, labels=labels_2, seed=42
        )
        npt.assert_array_equal(obs_mat, exp_mat)
        # Note that the synthetic labels won't be the same as those of `mixup`. This is
        # because `phylomix` shuffles pairs of sample indices to ensure randomness.
        exp_lab = np.array([
            [0.73, 0.27],
            [1.  , 0.  ],
            [0.82, 0.18],
            [0.69, 0.31],
            [0.61, 0.39],
            [0.64, 0.36],
            [1.  , 0.  ],
            [0.19, 0.81],
            [1.  , 0.  ],
            [0.45, 0.55],
        ])
        npt.assert_array_equal(obs_lab.round(2), exp_lab)

        # alpha parameter
        obs_mat, obs_lab = phylomix(
            matrix, n=10, tree=tree, taxa=taxa, alpha=1, seed=42
        )
        exp_mat = np.array([
            [ 0,  0,  2,  0,  1,  2,  3,  2,  0,  1],
            [ 0,  0,  1,  5,  2,  2,  2,  0,  4,  1],
            [ 0,  4,  0,  0,  5,  0,  2,  1,  3,  1],
            [ 0,  0,  1,  0,  0,  1,  1,  1,  0,  0],
            [ 0,  0,  1,  0,  0,  1,  1,  1,  0,  0],
            [ 0,  2,  2,  0,  1,  1,  1,  0,  2,  0],
            [ 2,  0,  4, 17,  0,  8,  3,  1, 15,  7],
            [ 0,  0,  0,  1,  0,  3,  0,  3,  1,  0],
            [ 9,  0, 16,  0,  0, 24,  0,  0,  0,  6],
            [ 0,  0,  2,  0,  1,  5,  0,  0,  0,  1],
        ])
        npt.assert_array_equal(obs_mat, exp_mat)

        # two-class labels, intra-class mixing
        obs_mat, obs_lab = phylomix(
            matrix, n=10, tree=tree, taxa=taxa, labels=labels_2,
            intra_class=True, seed=42
        )
        exp_mat = np.array([
            [ 0,  0,  0,  0,  0,  2,  2,  0,  0,  0],
            [ 0,  0,  0,  5,  0,  8,  0,  8,  2,  0],
            [ 0,  3,  5,  0,  2,  2,  2,  0,  4,  1],
            [ 0,  0, 24,  0,  0, 24,  4, 24,  0,  0],
            [ 0,  0,  5,  0,  0, 22,  4, 22,  0, 22],
            [ 0,  0,  2,  2,  0,  3,  0,  2,  1,  0],
            [ 2,  0,  3,  0,  9, 17,  3, 17,  0,  7],
            [ 0,  8, 13,  0,  5,  5,  5,  0, 10,  0],
            [ 0,  0, 20,  0,  8,  8,  3,  0, 16,  4],
            [ 0,  0,  0,  0,  0,  2,  2,  1,  0,  0],
        ])
        exp_lab = np.array([
            [1, 0],
            [0, 1],
            [1, 0],
            [1, 0],
            [1, 0],
            [0, 1],
            [1, 0],
            [1, 0],
            [1, 0],
            [1, 0],
        ])
        npt.assert_array_equal(obs_mat, exp_mat)
        npt.assert_array_equal(obs_lab, exp_lab)

        # three-class labels
        obs_mat, obs_lab = phylomix(
            matrix, n=10, tree=tree, taxa=taxa, labels=labels_3, seed=42
        )
        exp_mat = np.array([
            [ 0,  0,  0,  0,  1,  4,  3,  0,  0,  1],
            [ 0, 12, 21,  0,  8,  8,  8,  0, 17,  0],
            [ 0,  3,  5,  0,  2,  2,  2,  0,  4,  1],
            [ 0,  0,  2,  0,  0,  2,  3,  2,  0,  1],
            [ 0,  0,  2,  0,  0,  2,  3,  2,  0,  2],
            [ 0,  1,  2,  0,  1,  1,  1,  2,  2,  0],
            [ 2,  0,  3,  0,  9,  7,  3, 12, 14,  7],
            [ 0,  0,  0, 17,  0, 25,  0, 25,  8,  0],
            [ 0,  0, 13,  0,  0, 13,  0, 13,  0, 13],
            [ 0,  0,  0,  0,  1,  3,  3,  2,  0,  0],
        ])
        exp_lab = np.array([
            [0.27, 0.73, 0.  ],
            [0.21, 0.  , 0.79],
            [0.  , 0.  , 1.  ],
            [1.  , 0.  , 0.  ],
            [1.  , 0.  , 0.  ],
            [0.36, 0.  , 0.64],
            [0.34, 0.66, 0.  ],
            [1.  , 0.  , 0.  ],
            [0.77, 0.23, 0.  ],
            [0.55, 0.45, 0.  ],
        ])
        npt.assert_array_equal(obs_mat, exp_mat)
        npt.assert_array_equal(obs_lab.round(2), exp_lab)

    def test_phylomix_errors(self):
        # a simple example (2 x 5)
        matrix = np.arange(10).reshape(2, 5)
        labels = np.array([0, 1])
        taxa = list("abcde")
        tree = TreeNode.read(["(((a,b),c),(d,e));"])

        # normal situation
        obs_mat, obs_lab = phylomix(
            matrix, n=5, tree=tree, taxa=taxa, labels=labels, seed=42
        )
        exp_mat = np.array([
            [ 0,  2,  2,  2,  3],
            [ 0,  1,  1,  3,  4],
            [ 0,  1,  2,  3,  3],
            [ 0,  6,  8, 12,  9],
            [ 5,  6,  6,  9,  9],
        ])
        exp_lab = np.array([
            [0.42, 0.58],
            [0.53, 0.47],
            [0.64, 0.36],
            [0.48, 0.52],
            [0.32, 0.68],
        ])
        npt.assert_array_equal(obs_mat, exp_mat)
        npt.assert_array_equal(obs_lab.round(2), exp_lab)

        # taxa not provided
        msg = "Taxa must be included in table or explicitly provided."
        with self.assertRaises(ValueError) as cm:
            phylomix(matrix, n=5, tree=tree, taxa=None)
        self.assertEqual(str(cm.exception), msg)

        # taxa included in table
        df = pd.DataFrame(matrix, columns=taxa)
        obs_mat, obs_lab = phylomix(df, n=5, tree=tree, seed=42)
        npt.assert_array_equal(obs_mat, exp_mat)

        # explicitly provided taxa override feature IDs in table
        df = pd.DataFrame(matrix, columns=list("uvwxy"))
        obs_mat, obs_lab = phylomix(df, n=5, tree=tree, taxa=taxa, seed=42)
        npt.assert_array_equal(obs_mat, exp_mat)


if __name__ == "__main__":
    main()
