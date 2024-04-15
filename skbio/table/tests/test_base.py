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


if __name__ == '__main__':
    main()
