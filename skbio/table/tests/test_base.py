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

from skbio.table import Table, example_table
from skbio.table._base import _table_to_numpy


class TableTests(TestCase):
    def test_table(self):
        data = np.arange(40).reshape(10, 4)
        sample_ids = ['S%d' % i for i in range(4)]
        observ_ids = ['O%d' % i for i in range(10)]
        sample_metadata = [{'environment': 'A'},
                           {'environment': 'B'},
                           {'environment': 'A'},
                           {'environment': 'B'}]
        observ_metadata = [{'taxonomy': ['Bacteria', 'Firmicutes']},
                           {'taxonomy': ['Bacteria', 'Firmicutes']},
                           {'taxonomy': ['Bacteria', 'Proteobacteria']},
                           {'taxonomy': ['Bacteria', 'Proteobacteria']},
                           {'taxonomy': ['Bacteria', 'Proteobacteria']},
                           {'taxonomy': ['Bacteria', 'Bacteroidetes']},
                           {'taxonomy': ['Bacteria', 'Bacteroidetes']},
                           {'taxonomy': ['Bacteria', 'Firmicutes']},
                           {'taxonomy': ['Bacteria', 'Firmicutes']},
                           {'taxonomy': ['Bacteria', 'Firmicutes']}]

        table = Table(data, observ_ids, sample_ids,
                      observ_metadata, sample_metadata,
                      table_id='Example Table')

        self.assertEqual(list(table.ids()), ['S0', 'S1', 'S2', 'S3'])
        self.assertEqual(list(table.ids(axis='observation')),
                         ['O0', 'O1', 'O2', 'O3', 'O4',
                          'O5', 'O6', 'O7', 'O8', 'O9'])
        self.assertEqual(int(table.nnz), 39)


class TableUtilTests(TestCase):
    def test_table_to_numpy(self):
        exp_data = np.array([[0, 1, 2], [3, 4, 5]]).T
        exp_ids = ['S1', 'S2', 'S3']
        exp_feat_ids = ['O1', 'O2']
        obs_data, obs_ids, obs_feat_ids = _table_to_numpy(example_table)
        npt.assert_equal(obs_data, exp_data)
        self.assertEqual(obs_ids, exp_ids)
        self.assertEqual(obs_feat_ids, exp_feat_ids)


if __name__ == '__main__':
    main()
