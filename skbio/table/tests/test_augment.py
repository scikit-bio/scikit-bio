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
from skbio.table import Table
from skbio.tree import TreeNode
from skbio.table._augment import (
    mixup,
    aitchison_mixup,
    compositional_cutmix,
    phylomix,
    _validate_tree,
    _validate_label,
    _get_all_possible_pairs,
    _aitchison_addition,
    _aitchison_scalar_multiplication,
)


class TestAugmentFunctions(TestCase):
    def setUp(self):
        samples = 40
        n_features = 100
        self.data = np.arange(samples * n_features).reshape(samples, n_features)
        feature_metadata = [
            {"phylogeny": "a"},
            {"phylogeny": "b"},
            {"phylogeny": "x"},
            {"phylogeny": "y"},
            {"phylogeny": "w"},
            {"phylogeny": "z"},
            {"phylogeny": "c"},
            {"phylogeny": "d"},
            {"phylogeny": "e"},
            {"phylogeny": "f"},
            {"phylogeny": "a"},
            {"phylogeny": "b"},
            {"phylogeny": "x"},
            {"phylogeny": "y"},
            {"phylogeny": "w"},
            {"phylogeny": "a"},
            {"phylogeny": "b"},
            {"phylogeny": "x"},
            {"phylogeny": "a"},
            {"phylogeny": "b"},
            {"phylogeny": "x"},
            {"phylogeny": "y"},
            {"phylogeny": "w"},
            {"phylogeny": "z"},
            {"phylogeny": "c"},
            {"phylogeny": "d"},
            {"phylogeny": "e"},
            {"phylogeny": "f"},
            {"phylogeny": "y"},
            {"phylogeny": "a"},
            {"phylogeny": "b"},
            {"phylogeny": "x"},
            {"phylogeny": "y"},
            {"phylogeny": "w"},
            {"phylogeny": "z"},
            {"phylogeny": "c"},
            {"phylogeny": "d"},
            {"phylogeny": "e"},
            {"phylogeny": "a"},
            {"phylogeny": "b"},
            {"phylogeny": "x"},
            {"phylogeny": "y"},
            {"phylogeny": "w"},
            {"phylogeny": "z"},
            {"phylogeny": "c"},
            {"phylogeny": "d"},
            {"phylogeny": "e"},
            {"phylogeny": "f"},
            {"phylogeny": "f"},
            {"phylogeny": "w"},
            {"phylogeny": "z"},
            {"phylogeny": "c"},
            {"phylogeny": "a"},
            {"phylogeny": "b"},
            {"phylogeny": "x"},
            {"phylogeny": "y"},
            {"phylogeny": "w"},
            {"phylogeny": "z"},
            {"phylogeny": "c"},
            {"phylogeny": "d"},
            {"phylogeny": "e"},
            {"phylogeny": "f"},
            {"phylogeny": "d"},
            {"phylogeny": "e"},
            {"phylogeny": "f"},
            {"phylogeny": "z"},
            {"phylogeny": "c"},
            {"phylogeny": "d"},
            {"phylogeny": "a"},
            {"phylogeny": "b"},
            {"phylogeny": "x"},
            {"phylogeny": "y"},
            {"phylogeny": "w"},
            {"phylogeny": "z"},
            {"phylogeny": "c"},
            {"phylogeny": "a"},
            {"phylogeny": "b"},
            {"phylogeny": "x"},
            {"phylogeny": "y"},
            {"phylogeny": "w"},
            {"phylogeny": "z"},
            {"phylogeny": "c"},
            {"phylogeny": "d"},
            {"phylogeny": "e"},
            {"phylogeny": "f"},
            {"phylogeny": "d"},
            {"phylogeny": "a"},
            {"phylogeny": "b"},
            {"phylogeny": "x"},
            {"phylogeny": "y"},
            {"phylogeny": "w"},
            {"phylogeny": "z"},
            {"phylogeny": "c"},
            {"phylogeny": "d"},
            {"phylogeny": "e"},
            {"phylogeny": "f"},
            {"phylogeny": "e"},
            {"phylogeny": "f"},
            {"phylogeny": "e"},
            {"phylogeny": "f"},
        ]

        self.data_simple = np.arange(10).reshape(2, 5)
        feature_metadata_simple = [
            {"phylogeny": "a"},
            {"phylogeny": "b"},
            {"phylogeny": "c"},
            {"phylogeny": "x"},
            {"phylogeny": "y"},
        ]
        self.simple_tree = TreeNode.read(["(((a,b)int1,c)int2,(x,y)int3);"])
        self.complex_tree = TreeNode.read(
            ["(((a,b)int1,(x,y,(w,z)int2,(c,d)int3)int4),(e,f)int5);"]
        )
        self.complex_tree.bifurcate()
        self.simple_tree.bifurcate()
        tree_tips_simple = {tip.name for tip in self.simple_tree.tips()}
        tree_tips = {tip.name for tip in self.complex_tree.tips()}
        self.tips_to_feature_mapping = {}
        self.tips_to_feature_mapping_simple = {}
        for idx, metadata in enumerate(feature_metadata):
            if metadata and "phylogeny" in metadata:
                phylogeny_label = metadata["phylogeny"]
                if phylogeny_label in tree_tips:
                    self.tips_to_feature_mapping[phylogeny_label] = idx
        for idx, metadata in enumerate(feature_metadata_simple):
            if metadata and "phylogeny" in metadata:
                phylogeny_label = metadata["phylogeny"]
                if phylogeny_label in tree_tips_simple:
                    self.tips_to_feature_mapping_simple[phylogeny_label] = idx
        self.labels = np.random.randint(0, 2, size=self.data.shape[0])
        self.labels_simple = np.array([0, 1])
        self.one_hot_labels = np.eye(2)[self.labels]

    def test_validate_tree(self):
        # Test with valid tree
        _validate_tree(self.simple_tree)  # Should not raise

        # Test with None
        with self.assertRaisesRegex(
            TypeError, "`tree` must be a skbio.tree.TreeNode object."
        ):
            _validate_tree(None)

        # Test with invalid type
        with self.assertRaisesRegex(
            TypeError, "`tree` must be a skbio.tree.TreeNode object."
        ):
            _validate_tree("not_a_tree")

    def test_validate_label(self):
        # Create a simple matrix for testing
        matrix = np.array([[1, 2, 3], [4, 5, 6]])

        # Test with valid 1D label
        label = np.array([0, 1])
        obs_label, obs_one_hot = _validate_label(label, matrix)
        self.assertTrue(np.array_equal(obs_label, label))
        self.assertTrue(np.array_equal(obs_one_hot, np.array([[1, 0], [0, 1]])))

        # Test with valid 2D one-hot label
        one_hot = np.array([[1, 0], [0, 1]])
        obs_label, obs_one_hot = _validate_label(one_hot, matrix)
        self.assertTrue(np.array_equal(obs_label, np.array([0, 1])))
        self.assertTrue(np.array_equal(obs_one_hot, one_hot))

        # Test with wrong number of samples
        wrong_label = np.array([0, 1, 2])
        with self.assertRaisesRegex(ValueError, "Number of elements in label"):
            _validate_label(wrong_label, matrix)

        # Test with non-numpy array
        with self.assertRaisesRegex(ValueError, "label must be a numpy.ndarray"):
            _validate_label([0, 1], matrix)

        # Test with wrong dimensions
        wrong_dim_label = np.array([[[0], [1]]])
        with self.assertRaisesRegex(ValueError, "labels should have shape"):
            _validate_label(wrong_dim_label, matrix)

        # Test with non-zero indexed labels
        non_zero_label = np.array([1, 2])
        with self.assertRaisesRegex(ValueError, "Labels must be zero-indexed"):
            _validate_label(non_zero_label, matrix)

        # Test with non-consecutive labels
        non_consecutive = np.array([0, 2])
        with self.assertRaisesRegex(ValueError, "Labels must be consecutive integers"):
            _validate_label(non_consecutive, matrix)

        # Test with invalid one-hot encoding
        invalid_one_hot = np.array([[1, 1], [0, 1]])
        with self.assertRaisesRegex(
            ValueError, "label is not properly one hot encoded"
        ):
            _validate_label(invalid_one_hot, matrix)

    def test_aitchison_operations(self):
        # Test Aitchison addition
        x = np.array([0.2, 0.3, 0.5])
        v = np.array([0.1, 0.4, 0.5])
        res = _aitchison_addition(x, v)
        npt.assert_array_almost_equal(
            res, np.array([0.05128205, 0.30769231, 0.64102564])
        )
        self.assertAlmostEqual(np.sum(res), 1.0)

        # Test Aitchison scalar multiplication
        lam = 2.0
        x = np.array([0.2, 0.3, 0.5])
        res = _aitchison_scalar_multiplication(lam, x)
        npt.assert_array_almost_equal(
            res, np.array([0.10526316, 0.23684211, 0.65789474])
        )
        self.assertAlmostEqual(np.sum(res), 1.0)

    def test_get_all_possible_pairs(self):
        matrix = np.arange(40).reshape(4, 10)
        labels = np.array([0, 0, 1, 1])

        # Test all pairs
        all_pairs = _get_all_possible_pairs(matrix)
        self.assertEqual(len(all_pairs), 6)  # 4 choose 2

        # Test intra-class pairs
        intra_pairs = _get_all_possible_pairs(matrix, label=labels, intra_class=True)
        self.assertEqual(len(intra_pairs), 2)
        self.assertTrue(np.all(intra_pairs == np.array([(0, 1), (2, 3)])))

        # Test intra-class without label
        with self.assertRaisesRegex(
            ValueError, "Label is required for intra-class augmentation."
        ):
            _get_all_possible_pairs(matrix, intra_class=True)

    def test_mixup(self):
        augmented_matrix, augmented_label = mixup(
            self.data, samples=10, label=self.labels, alpha=2
        )

        self.assertEqual(augmented_matrix.shape[0], self.data.shape[0] + 10)
        self.assertEqual(augmented_matrix.shape[1], self.data.shape[1])
        self.assertEqual(augmented_label.shape[0], len(self.labels) + 10)
        self.assertEqual(augmented_label.shape[1], 2)
        self.assertTrue(np.allclose(np.sum(augmented_label, axis=1), 1.0))

    def test_mixup_no_label(self):
        augmented_matrix, augmented_label = mixup(self.data, samples=10, alpha=2)
        self.assertTrue(augmented_label.empty)

    def test_aitchison_mixup(self):
        augmented_matrix, augmented_label = aitchison_mixup(
            self.data, samples=20, label=self.labels, alpha=2
        )

        self.assertEqual(augmented_matrix.shape[0], self.data.shape[0] + 20)
        self.assertEqual(augmented_matrix.shape[1], self.data.shape[1])
        # Check if the augmented data is compositional
        self.assertTrue(np.allclose(np.sum(augmented_matrix, axis=1), 1.0))
        self.assertEqual(augmented_label.shape[0], len(self.labels) + 20)
        self.assertEqual(augmented_label.shape[1], 2)
        self.assertTrue(np.allclose(np.sum(augmented_label, axis=1), 1.0))

    # the previous test should be enough, but I'm not clear whether users should
    # even be allowed to turn of normalization for this function?
    # def test_aitchison_mixup_non_compositional(self):
    #     # Test with non-compositional data (should normalize)
    #     augmented_matrix, augmented_label = aitchison_mixup(
    #         self.table, samples=20, label=self.labels, alpha=2
    #     )

    #     self.assertEqual(augmented_matrix.shape[0], self.table.shape[1] + 20)
    #     self.assertEqual(augmented_matrix.shape[1], self.table.shape[0])
    #     # Check if the augmented data is compositional
    #     self.assertTrue(np.allclose(np.sum(augmented_matrix, axis=1), 1.0))
    #     self.assertEqual(augmented_label.shape[0], len(self.labels) + 20)
    #     self.assertEqual(augmented_label.shape[1], 2)
    #     self.assertTrue(np.allclose(np.sum(augmented_label, axis=1), 1.0))

    def test_aitchison_mixup_no_label(self):
        augmented_matrix, augmented_label = aitchison_mixup(
            self.data, samples=20, alpha=2
        )
        self.assertTrue(augmented_label.empty)

    def test_compositional_cutmix(self):
        augmented_matrix, augmented_label = compositional_cutmix(
            self.data, samples=20, label=self.labels, seed=42
        )

        self.assertEqual(augmented_matrix.shape[0], self.data.shape[0] + 20)
        self.assertEqual(augmented_matrix.shape[1], self.data.shape[1])
        # Check if the augmented data is compositional
        # this rtol is really low, the data seems to be close to compositionality, but
        # not quite, it is within 0.9 - 1.1
        self.assertTrue(np.allclose(np.sum(augmented_matrix, axis=1), 1.0, rtol=1e-01))
        self.assertEqual(augmented_label.shape[0], len(self.labels) + 20)
        # compositional_cutmix returns a 2D output for labels
        self.assertEqual(augmented_label.shape[1], 2)

    def test_phylomix_simple(self):
        augmented_matrix, augmented_label = phylomix(
            self.data_simple,
            tree=self.simple_tree,
            tip_to_obs_mapping=self.tips_to_feature_mapping_simple,
            samples=20,
            label=self.labels_simple,
        )

        self.assertEqual(augmented_matrix.shape[0], self.data_simple.shape[0] + 20)
        self.assertEqual(augmented_matrix.shape[1], self.data_simple.shape[1])

    def test_phylomix_no_tree(self):
        with self.assertRaisesRegex(
            TypeError, "`tree` must be a skbio.tree.TreeNode object."
        ):
            phylomix(
                self.data_simple,
                tree=None,
                tip_to_obs_mapping=self.tips_to_feature_mapping_simple,
                samples=20,
                label=self.labels_simple,
            )

    def test_phylomix_bad_tips(self):
        bad_mapping = {
            k: v for k, v in self.tips_to_feature_mapping_simple.items() if k != "a"
        }
        with self.assertRaisesRegex(
            ValueError, "tip_to_obs_mapping must contain all tips in the tree"
        ):
            phylomix(
                self.data_simple,
                tree=self.simple_tree,
                tip_to_obs_mapping=bad_mapping,
                samples=20,
                label=self.labels_simple,
            )

    def test_phylomix_no_label(self):
        augmented_matrix, augmented_label = phylomix(
            self.data_simple,
            tree=self.simple_tree,
            tip_to_obs_mapping=self.tips_to_feature_mapping_simple,
            samples=20,
        )

        self.assertEqual(augmented_matrix.shape[0], self.data_simple.shape[0] + 20)
        self.assertEqual(augmented_matrix.shape[1], self.data_simple.shape[1])
        self.assertTrue(augmented_label.empty)

    def test_phylomix(self):
        augmented_matrix, augmented_label = phylomix(
            self.data,
            tree=self.complex_tree,
            tip_to_obs_mapping=self.tips_to_feature_mapping,
            samples=20,
            label=self.labels,
        )

        self.assertEqual(augmented_matrix.shape[0], self.data.shape[0] + 20)
        self.assertEqual(augmented_matrix.shape[1], self.data.shape[1])
        self.assertEqual(augmented_label.shape[0], len(self.labels) + 20)
        self.assertEqual(augmented_label.shape[1], 2)
        self.assertTrue(np.allclose(np.sum(augmented_label, axis=1), 1.0))

    def test_multiclass_phylomix(self):
        labels_multiclass = np.random.randint(0, 3, size=self.data.shape[0])
        augmented_matrix, augmented_label = phylomix(
            self.data,
            tree=self.complex_tree,
            tip_to_obs_mapping=self.tips_to_feature_mapping,
            samples=20,
            label=labels_multiclass,
        )

        self.assertEqual(augmented_matrix.shape[0], self.data.shape[0] + 20)
        self.assertEqual(augmented_matrix.shape[1], self.data.shape[1])
        self.assertEqual(augmented_label.shape[0], len(labels_multiclass) + 20)
        self.assertEqual(augmented_label.shape[1], 3)  # 3 classes


if __name__ == "__main__":
    main(buffer=False)
