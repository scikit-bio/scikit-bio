# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from skbio.table import Table
import numpy as np
from skbio.table._augment import Augmentation
from skbio.tree import TreeNode

class TestAugmentation(TestCase):

    def setUp(self):
        data = np.arange(40).reshape(10, 4)
        sample_ids = ['S%d' % i for i in range(4)]
        observ_ids = ['O%d' % i for i in range(10)]
        observ_metadata = [{'phylogeny': 'a'},
                             {'phylogeny': 'b'},
                             {'phylogeny': 'x'},
                             {'phylogeny': 'y'},
                             {'phylogeny': 'w'},
                             {'phylogeny': 'z'},
                             {'phylogeny': 'c'},
                             {'phylogeny': 'd'},
                             {'phylogeny': 'e'},
                             {'phylogeny': 'f'},
                             ]
        # a complex tree to test the phylomix method
        self.table = Table(data, observ_ids, sample_ids, observ_metadata)
        
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
        tree_tips = {tip.name for tip in self.complex_tree.tips()}
        self.tips_to_obs_mapping = {}
        for idx, metadata in enumerate(self.table.metadata(axis="observation")):
            if metadata and "phylogeny" in metadata:
                phylogeny_label = metadata["phylogeny"]
                if phylogeny_label in tree_tips:
                    self.tips_to_obs_mapping[phylogeny_label] = idx

        self.labels = np.random.randint(0, 2, size=self.table.shape[0])

    def test_mixup(self):
        augmentation = Augmentation(self.table, method="mixup", label=self.labels)

        augmented_matrix, augmented_label = augmentation.mixup(alpha=2, n_samples=10)

        self.assertEqual(augmented_matrix.shape[1], self.table.shape[1] + 10)
        self.assertEqual(augmented_matrix.shape[0], self.table.shape[0])
        self.assertEqual(len(augmented_label), len(self.labels) + 10)
        self.assertEqual(len(augmented_label[0]), 2)
        self.assertTrue(np.allclose(np.sum(augmented_label, axis=1), 1.0))

    def test_aitchison_mixup(self):
        table_compositional = self.table.norm(axis="sample")
        augmentation = Augmentation(table_compositional, method="aitchison_mixup", label=self.labels)

        augmented_matrix, augmented_label = augmentation.aitchison_mixup(alpha=2, n_samples=20)

        self.assertEqual(augmented_matrix.shape[1], self.table.shape[1] + 20)
        self.assertEqual(augmented_matrix.shape[0], self.table.shape[0])
        # check if the augmented data is compositional
        self.assertTrue(np.allclose(np.sum(augmented_matrix, axis=0), 1.0))
        self.assertEqual(len(augmented_label), len(self.labels) + 20)
        self.assertEqual(len(augmented_label[0]), 2)
        self.assertTrue(np.allclose(np.sum(augmented_label, axis=1), 1.0))


    def test_phylomix(self):
        augmentation = Augmentation(self.table, method="phylomix", label=self.labels, tree=self.complex_tree)

        augmented_matrix, augmented_label = augmentation.phylomix(self.tips_to_obs_mapping, n_samples=20)

        self.assertEqual(augmented_matrix.shape[1], self.table.shape[1] + 20)
        self.assertEqual(augmented_matrix.shape[0], self.table.shape[0])
        # check if the entry are all integers
        self.assertTrue(np.all(augmented_matrix == np.round(augmented_matrix)))
        self.assertEqual(len(augmented_label), len(self.labels) + 20)
        self.assertEqual(len(augmented_label[0]), 2)
        self.assertTrue(np.allclose(np.sum(augmented_label, axis=1), 1.0))

if __name__ == '__main__':
    main()
