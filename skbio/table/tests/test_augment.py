# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import numpy as np
from skbio.table import Table
from skbio.table._augment import Augmentation
from skbio.tree import TreeNode

class TestAugmentation(TestCase):

    def setUp(self):
        n_samples = 40
        n_features = 100
        data = np.arange(n_samples * n_features).reshape(n_features, n_samples)
        sample_ids = ['S%d' % i for i in range(n_samples)]
        feature_ids = ['O%d' % i for i in range(n_features)]
        feature_metadata = [{'phylogeny': 'a'},
                             {'phylogeny': 'b'},
                             {'phylogeny': 'x'},
                             {'phylogeny': 'y'},
                             {'phylogeny': 'w'},
                             {'phylogeny': 'z'},
                             {'phylogeny': 'c'},
                             {'phylogeny': 'd'},
                             {'phylogeny': 'e'},
                             {'phylogeny': 'f'},
                             {'phylogeny': 'a'},
                             {'phylogeny': 'b'},
                             {'phylogeny': 'x'},
                             {'phylogeny': 'y'},
                             {'phylogeny': 'w'},
                             {'phylogeny': 'a'},
                             {'phylogeny': 'b'},
                             {'phylogeny': 'x'},
                             {'phylogeny': 'a'},
                             {'phylogeny': 'b'},
                             {'phylogeny': 'x'},
                             {'phylogeny': 'y'},
                             {'phylogeny': 'w'},
                             {'phylogeny': 'z'},
                             {'phylogeny': 'c'},
                             {'phylogeny': 'd'},
                             {'phylogeny': 'e'},
                             {'phylogeny': 'f'},
                             {'phylogeny': 'y'},
                             {'phylogeny': 'a'},
                             {'phylogeny': 'b'},
                             {'phylogeny': 'x'},
                             {'phylogeny': 'y'},
                             {'phylogeny': 'w'},
                             {'phylogeny': 'z'},
                             {'phylogeny': 'c'},
                             {'phylogeny': 'd'},
                             {'phylogeny': 'e'},
                             {'phylogeny': 'a'},
                             {'phylogeny': 'b'},
                             {'phylogeny': 'x'},
                             {'phylogeny': 'y'},
                             {'phylogeny': 'w'},
                             {'phylogeny': 'z'},
                             {'phylogeny': 'c'},
                             {'phylogeny': 'd'},
                             {'phylogeny': 'e'},
                             {'phylogeny': 'f'},
                             {'phylogeny': 'f'},
                             {'phylogeny': 'w'},
                             {'phylogeny': 'z'},
                             {'phylogeny': 'c'},
                             {'phylogeny': 'a'},
                             {'phylogeny': 'b'},
                             {'phylogeny': 'x'},
                             {'phylogeny': 'y'},
                             {'phylogeny': 'w'},
                             {'phylogeny': 'z'},
                             {'phylogeny': 'c'},
                             {'phylogeny': 'd'},
                             {'phylogeny': 'e'},
                             {'phylogeny': 'f'},
                             {'phylogeny': 'd'},
                             {'phylogeny': 'e'},
                             {'phylogeny': 'f'},
                             {'phylogeny': 'z'},
                             {'phylogeny': 'c'},
                             {'phylogeny': 'd'},
                             {'phylogeny': 'a'},
                             {'phylogeny': 'b'},
                             {'phylogeny': 'x'},
                             {'phylogeny': 'y'},
                             {'phylogeny': 'w'},
                             {'phylogeny': 'z'},
                             {'phylogeny': 'c'},
                             {'phylogeny': 'a'},
                             {'phylogeny': 'b'},
                             {'phylogeny': 'x'},
                             {'phylogeny': 'y'},
                             {'phylogeny': 'w'},
                             {'phylogeny': 'z'},
                             {'phylogeny': 'c'},
                             {'phylogeny': 'd'},
                             {'phylogeny': 'e'},
                             {'phylogeny': 'f'},
                             {'phylogeny': 'd'},
                             {'phylogeny': 'a'},
                             {'phylogeny': 'b'},
                             {'phylogeny': 'x'},
                             {'phylogeny': 'y'},
                             {'phylogeny': 'w'},
                             {'phylogeny': 'z'},
                             {'phylogeny': 'c'},
                             {'phylogeny': 'd'},
                             {'phylogeny': 'e'},
                             {'phylogeny': 'f'},
                             {'phylogeny': 'e'},
                             {'phylogeny': 'f'},
                             {'phylogeny': 'e'},
                             {'phylogeny': 'f'},
                             ]
        # a complex tree to test the phylomix method
        self.table = Table(data, feature_ids, sample_ids, feature_metadata)
        self.table_min = np.min(self.table.to_dataframe().values)
        self.table_max = np.max(self.table.to_dataframe().values)

        data_simple = np.arange(10).reshape(5, 2)
        sample_ids_simple = ['S%d' % i for i in range(2)]
        feature_ids_simple = ['O%d' % i for i in range(5)]
        feature_metadata_simple = [{'phylogeny': 'a'},
                                  {'phylogeny': 'b'},
                                  {'phylogeny': 'c'},
                                  {'phylogeny': 'x'},
                                  {'phylogeny': 'y'},
                                  ]
        self.table_simple = Table(data_simple, feature_ids_simple, sample_ids_simple, feature_metadata_simple)
        self.simple_tree = TreeNode.read([
            "(((a,b)int1,c)int2,(x,y)int3);"])
        #                              /-a
        #                    /int1----|
        #          /int2----|          \-b
        #         |         |
        #---------|          \-c
        #         |
        #         |          /-x
        #          \int3----|
        #                    \-y

        self.complex_tree = TreeNode.read([
            "(((a,b)int1,(x,y,(w,z)int2,(c,d)int3)int4),(e,f)int5);"])
        #                               /-a
        #                     /int1----|
        #                    |          \-b
        #                    |
        #           /--------|          /-x
        #          |         |         |
        #          |         |         |--y
        #          |         |         |
        #          |          \int4----|          /-w
        #          |                   |-int2----|
        # ---------|                   |          \-z
        #          |                   |
        #          |                   |          /-c
        #          |                    \int3----|
        #          |                              \-d
        #          |
        #          |          /-e
        #           \int5----|
        #                     \-f
        self.complex_tree.bifurcate()
        self.simple_tree.bifurcate()
        tree_tips_simple = {tip.name for tip in self.simple_tree.tips()}
        tree_tips = {tip.name for tip in self.complex_tree.tips()}
        self.tips_to_feature_mapping = {}
        self.tips_to_feature_mapping_simple = {}
        for idx, metadata in enumerate(self.table.metadata(axis="observation")):
            if metadata and "phylogeny" in metadata:
                phylogeny_label = metadata["phylogeny"]
                if phylogeny_label in tree_tips:
                    self.tips_to_feature_mapping[phylogeny_label] = idx
        for idx, metadata in enumerate(self.table_simple.metadata(axis="observation")):
            if metadata and "phylogeny" in metadata:
                phylogeny_label = metadata["phylogeny"]
                if phylogeny_label in tree_tips_simple:
                    self.tips_to_feature_mapping_simple[phylogeny_label] = idx
        self.labels = np.random.randint(0, 2, size=self.table.shape[1])
        self.labels_simple = np.array([0, 1])# np.random.randint(0, 2, size=self.table_simple.shape[1])
        self.one_hot_labels = np.eye(2)[self.labels]

    def test_init(self):
        augmentation = Augmentation(self.table, label=self.labels) #, num_classes=2)
        self.assertIsInstance(augmentation, Augmentation)
        self.assertEqual(augmentation.table, self.table)
        self.assertTrue(np.array_equal(augmentation.label, self.labels))
        
        augmentation_with_tree = Augmentation(
            self.table, label=self.labels, tree=self.complex_tree
            )
        self.assertEqual(augmentation_with_tree.tree, self.complex_tree)
        
        with self.assertRaisesRegex(ValueError, "table must be a skbio.table.Table"):
            Augmentation(np.array([[1, 2, 3], [4, 5, 6]]))
        
        with self.assertRaisesRegex(ValueError, "tree must be a skbio.tree.TreeNode"):
            Augmentation(self.table, self.labels, tree="not_a_tree")

        # with self.assertRaisesRegex(ValueError, "num_classes must be provided if the labels are 1D"):
        #     # when label is 1D but num_classes is not provided
        #     Augmentation(self.table, label=self.labels) #, num_classes=None)

        with self.assertRaises(ValueError):
            # when label is 2D, num_classes is provided 
            # but does not match the number of classes in the label
            Augmentation(self.table, self.one_hot_labels) # , num_classes=3)

        with self.assertRaises(ValueError):
            # when label has wrong dimension
            self.wrong_label = self.one_hot_labels[..., np.newaxis]
            Augmentation(self.table, self.wrong_label) #, num_classes=2)

        with self.assertRaisesRegex(
            ValueError,
            "label must be a numpy.ndarray, but got <class 'int'> instead."
        ):
            Augmentation(self.table, label=3)

    def test_str(self):
        exp = "Augmentation(shape=(40, 100))"
        obs = str(Augmentation(self.table))
        self.assertEqual(exp, obs)

    def test_get_all_possible_pairs(self):
        data_simple = np.arange(40).reshape(10, 4)
        sample_ids = ['S%d' % i for i in range(4)]
        feature_ids = ['O%d' % i for i in range(10)]
        labels_simple = np.array([0, 0, 1, 1])
        table_simple = Table(data_simple, feature_ids, sample_ids)
        augmentation = Augmentation(table_simple, label=labels_simple) #, num_classes=2)
        possible_pairs = augmentation._get_all_possible_pairs()
        possible_pairs_intra = augmentation._get_all_possible_pairs(intra_class=True)

        self.assertEqual(len(possible_pairs), 6)
        self.assertEqual(len(possible_pairs_intra), 2)
        self.assertTrue(np.all(possible_pairs_intra == np.array([(0, 1), (2, 3)])))

    def test_get_all_possible_pairs_no_label(self):
        data_simple = np.arange(40).reshape(10, 4)
        sample_ids = ['S%d' % i for i in range(4)]
        feature_ids = ['O%d' % i for i in range(10)]
        labels_simple = np.array([0, 0, 1, 1])
        table_simple = Table(data_simple, feature_ids, sample_ids)
        augmentation = Augmentation(table_simple, label=None) #, num_classes=2)
        with self.assertRaisesRegex(ValueError, "label is required for intra-class augmentation"):
            possible_pairs = augmentation._get_all_possible_pairs(intra_class=True)


    def test_mixup(self):
        augmentation = Augmentation(self.table, label=self.labels) #, num_classes=2)

        augmented_matrix, augmented_label = augmentation.mixup(alpha=2, n_samples=10)

        self.assertEqual(augmented_matrix.shape[1], self.table.shape[0])
        self.assertEqual(augmented_matrix.shape[0], self.table.shape[1] + 10)
        self.assertEqual(len(augmented_label), len(self.labels) + 10)
        self.assertEqual(len(augmented_label[0]), 2)
        self.assertTrue(np.allclose(np.sum(augmented_label, axis=1), 1.0))
        self.assertTrue(np.all(augmented_matrix >= self.table_min))
        self.assertTrue(np.all(augmented_matrix <= self.table_max))

    def test_mixup_no_label(self):
        augmentation = Augmentation(self.table, label=None)
        augmented_matrix, augmented_label = augmentation.mixup(alpha=2, n_samples=10)
        self.assertIsNone(augmented_label)

    def test_aitchison_mixup(self):
        table_compositional = self.table.norm(axis="sample")
        augmentation = Augmentation(table_compositional, label=self.labels) #, num_classes=2)

        augmented_matrix, augmented_label = augmentation.aitchison_mixup(alpha=2, n_samples=20)

        self.assertEqual(augmented_matrix.shape[1], self.table.shape[0])
        self.assertEqual(augmented_matrix.shape[0], self.table.shape[1] + 20)
        # check if the augmented data is compositional
        self.assertTrue(np.allclose(np.sum(augmented_matrix, axis=1), 1.0))
        self.assertEqual(len(augmented_label), len(self.labels) + 20)
        self.assertEqual(len(augmented_label[0]), 2)
        self.assertTrue(np.allclose(np.sum(augmented_label, axis=1), 1.0))
        self.assertTrue(np.all(augmented_matrix >= self.table_min))
        self.assertTrue(np.all(augmented_matrix <= self.table_max))

    def test_aitchison_mixup_non_compositional(self):
        augmentation = Augmentation(self.table, label=self.labels) #, num_classes=2)

        augmented_matrix, augmented_label = augmentation.aitchison_mixup(alpha=2, n_samples=20)

        self.assertEqual(augmented_matrix.shape[1], self.table.shape[0])
        self.assertEqual(augmented_matrix.shape[0], self.table.shape[1] + 20)
        # check if the augmented data is compositional
        self.assertTrue(np.allclose(np.sum(augmented_matrix, axis=1), 1.0))
        self.assertEqual(len(augmented_label), len(self.labels) + 20)
        self.assertEqual(len(augmented_label[0]), 2)
        self.assertTrue(np.allclose(np.sum(augmented_label, axis=1), 1.0))
        self.assertTrue(np.all(augmented_matrix >= self.table_min))
        self.assertTrue(np.all(augmented_matrix <= self.table_max))

    def test_aitchison_mixup_no_label(self):
        augmentation = Augmentation(self.table, label=None)
        augmented_matrix, augmented_label = augmentation.aitchison_mixup(alpha=2, n_samples=20)
        self.assertIsNone(augmented_label)

    def test_phylomix_simple(self):
        augmentation = Augmentation(self.table_simple, label=self.labels_simple, tree=self.simple_tree)

        augmented_matrix, _ = augmentation.phylomix(self.tips_to_feature_mapping_simple, n_samples=20)

        self.assertEqual(augmented_matrix.shape[1], self.table_simple.shape[0])
        self.assertEqual(augmented_matrix.shape[0], self.table_simple.shape[1] + 20)
        self.assertTrue(np.all(augmented_matrix >= self.table_min))
        self.assertTrue(np.all(augmented_matrix <= self.table_max))

    def test_phylomix_no_tree(self):
        augmentation = Augmentation(self.table_simple, label=self.labels_simple, tree=None)
        with self.assertRaisesRegex(ValueError, "tree is required for phylomix augmentation"):
            augmentation.phylomix(self.tips_to_feature_mapping_simple, n_samples=20)

    def test_phylomix_bad_tips(self):
        bad_mapping = {k: v for k, v in self.tips_to_feature_mapping_simple.items() if k != 'a'}
        augmentation = Augmentation(self.table_simple, label=self.labels_simple, tree=self.simple_tree)
        with self.assertRaisesRegex(ValueError, "tip_to_obs_mapping must contain all tips in the tree"):
            augmentation.phylomix(bad_mapping, n_samples=20)

    def test_phylomix_no_label(self):
        augmentation = Augmentation(self.table_simple, label=None, tree=self.simple_tree)

        augmented_matrix, augmented_label = augmentation.phylomix(self.tips_to_feature_mapping_simple, n_samples=20)

        self.assertEqual(augmented_matrix.shape[1], self.table_simple.shape[0])
        self.assertEqual(augmented_matrix.shape[0], self.table_simple.shape[1] + 20)
        self.assertTrue(np.all(augmented_matrix >= self.table_min))
        self.assertTrue(np.all(augmented_matrix <= self.table_max))

        self.assertIsNone(augmented_label)

    def test_phylomix(self):
        augmentation = Augmentation(self.table, label=self.labels, tree=self.complex_tree)

        augmented_matrix, augmented_label = augmentation.phylomix(self.tips_to_feature_mapping, n_samples=20)

        self.assertEqual(augmented_matrix.shape[1], self.table.shape[0])
        self.assertEqual(augmented_matrix.shape[0], self.table.shape[1] + 20)
        # check if the entry are all integers
        self.assertTrue(np.all(augmented_matrix == np.round(augmented_matrix)))
        self.assertEqual(len(augmented_label), len(self.labels) + 20)
        self.assertTrue(len(augmented_label[0]) == 2)
        self.assertTrue(np.allclose(np.sum(augmented_label, axis=1), 1.0))
        self.assertTrue(np.all(augmented_matrix >= self.table_min))
        self.assertTrue(np.all(augmented_matrix <= self.table_max))

    def test_compositional_cutmix(self):
        table_compositional = self.table.norm(axis="sample")
        augmentation = Augmentation(table_compositional, label=self.labels)

        augmented_matrix, augmented_label = augmentation.compositional_cutmix(n_samples=20)

        self.assertEqual(augmented_matrix.shape[1], self.table.shape[0])
        self.assertEqual(augmented_matrix.shape[0], self.table.shape[1] + 20)
        self.assertEqual(len(augmented_label), len(self.labels) + 20)
        self.assertTrue(np.all(augmented_matrix >= self.table_min))
        self.assertTrue(np.all(augmented_matrix <= self.table_max))

    def test_multiclass_phylomix(self):
        labels = np.random.randint(0, 3, size=self.table.shape[1])
        augmentation = Augmentation(self.table, label=labels, tree=self.complex_tree)

        augmented_matrix, augmented_label = augmentation.phylomix(self.tips_to_feature_mapping, n_samples=20)

        self.assertEqual(augmented_matrix.shape[1], self.table.shape[0])
        self.assertEqual(augmented_matrix.shape[0], self.table.shape[1] + 20)
        self.assertEqual(len(augmented_label), len(self.labels) + 20)
        self.assertTrue(np.all(augmented_matrix >= self.table_min))
        self.assertTrue(np.all(augmented_matrix <= self.table_max))


if __name__ == '__main__':
    main()
