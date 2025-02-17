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


class TestAugmentation(TestCase):
    def setUp(self):
        # Create a mock table
        data = np.array([[1, 2], [3, 4], [5, 6]])
        sample_ids = ['S1', 'S2']
        observ_ids = ['O1', 'O2', 'O3']
        self.table = Table(data, observ_ids, sample_ids)

        # Create mock labels
        self.labels = np.array([0, 1, 1])
        self.labels = self.labels.reshape(-1, 1)

    def test_mixup(self):
        augmentation = Augmentation(self.table, method="mixup", label=self.labels)

        augmented_matrix, augmented_label = augmentation.mixup(alpha=2, n_samples=2)

        self.assertEqual(augmented_matrix.shape[0], self.table.shape[0] + 2)
        self.assertEqual(augmented_matrix.shape[1], self.table.shape[1])
        self.assertEqual(len(augmented_label), len(self.labels) + 2)


if __name__ == '__main__':
    main()
