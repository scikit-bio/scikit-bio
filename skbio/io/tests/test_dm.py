# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils.six import StringIO

from unittest import TestCase, main

import numpy as np

from skbio.io import DMFormatError
from skbio.io.dm import dm_to_DissimilarityMatrix
from skbio.stats.distance import DissimilarityMatrix
from skbio.util import get_data_path


class DMToDissimilarityMatrixTests(TestCase):
    def setUp(self):
        self.dm_1x1_data = [[0.0]]
        self.dm_1x1_f = StringIO(DM_1x1_F)

        self.dm_2x2_data = [[0.0, 0.123], [0.123, 0.0]]
        self.dm_2x2_f = StringIO(DM_2x2_F)

        self.dm_2x2_asym_data = [[0.0, 1.0], [-2.0, 0.0]]
        self.dm_2x2_asym_f = StringIO(DM_2x2_ASYM_F)

        self.dm_3x3_data = [[0.0, 0.01, 4.2], [0.01, 0.0, 12.0],
                            [4.2, 12.0, 0.0]]
        self.dm_3x3_f = StringIO(DM_3x3_F)
        self.dm_3x3_whitespace_f = StringIO('\n'.join(DM_3x3_WHITESPACE_F))

        self.bad_dm_f1 = StringIO(BAD_DM_F1)
        self.bad_dm_f2 = StringIO(BAD_DM_F2)
        self.bad_dm_f3 = StringIO(BAD_DM_F3)
        self.bad_dm_f4 = StringIO(BAD_DM_F4)
        self.bad_dm_f5 = StringIO(BAD_DM_F5)

        self.dm_1x1 = DissimilarityMatrix(self.dm_1x1_data, ['a'])
        self.dm_2x2 = DissimilarityMatrix(self.dm_2x2_data, ['a', 'b'])
        self.dm_2x2_asym = DissimilarityMatrix(self.dm_2x2_asym_data,
                                               ['a', 'b'])
        self.dm_3x3 = DissimilarityMatrix(self.dm_3x3_data, ['a', 'b', 'c'])

        self.dms = [self.dm_1x1, self.dm_2x2, self.dm_2x2_asym, self.dm_3x3,
                    self.dm_3x3]
        self.dm_fs = [self.dm_1x1_f, self.dm_2x2_f, self.dm_2x2_asym_f,
                      self.dm_3x3_f, self.dm_3x3_whitespace_f]
        self.dm_shapes = [(1, 1), (2, 2), (2, 2), (3, 3)]
        self.dm_sizes = [1, 4, 4, 9]
        self.dm_transposes = [
            self.dm_1x1, self.dm_2x2,
            DissimilarityMatrix([[0, -2], [1, 0]], ['a', 'b']), self.dm_3x3]
        self.dm_redundant_forms = [np.array(self.dm_1x1_data),
                                   np.array(self.dm_2x2_data),
                                   np.array(self.dm_2x2_asym_data),
                                   np.array(self.dm_3x3_data)]

    def test_valid_files(self):
        for dm_f, dm in zip(self.dm_fs, self.dms):
            obs = dm_to_DissimilarityMatrix(dm_f)
            self.assertEqual(obs, dm)
            self.assertIsInstance(obs, DissimilarityMatrix)

    def test_invalid_files(self):
        # Empty dm.
        with self.assertRaises(DMFormatError):
            dm_to_DissimilarityMatrix(StringIO())

        # Number of values don't match number of IDs.
        with self.assertRaises(DMFormatError):
            dm_to_DissimilarityMatrix(self.bad_dm_f1)

        # Mismatched IDs.
        with self.assertRaises(DMFormatError):
            dm_to_DissimilarityMatrix(self.bad_dm_f2)

        # Extra data at end.
        with self.assertRaises(DMFormatError):
            dm_to_DissimilarityMatrix(self.bad_dm_f3)

        # Missing data.
        with self.assertRaises(DMFormatError):
            dm_to_DissimilarityMatrix(self.bad_dm_f4)

        # Header, but no data.
        with self.assertRaises(DMFormatError):
            dm_to_DissimilarityMatrix(self.bad_dm_f5)


# 1x1:
#     0.0
DM_1x1_F = "\ta\na\t0.0\n"

# 2x2:
#       0.0  0.123
#     0.123    0.0
DM_2x2_F = "\ta\tb\na\t0.0\t0.123\nb\t0.123\t0.0\n"

DM_2x2_ASYM_F = "\ta\tb\na\t0.0\t1.0\nb\t-2.0\t0.0\n"

DM_3x3_F = "a\tb\tc\na\t0.0\t0.01\t4.2\nb\t0.01\t0.0\t12.0\nc\t4.2\t12.0\t0.0"

# Extra whitespace-only lines throughout. Also has comments before the header.
DM_3x3_WHITESPACE_F = ['# foo',
                       '      \t \t ',
                       ' #bar',
                       '',
                       '',
                       '\ta\t b \tc',
                       'a  \t0.0\t0.01\t4.2',
                       '     \t',
                       'b\t0.01\t0.0\t12.0',
                       '',
                       '\t     \t',
                       '',
                       'c\t4.2\t12.0\t0.0',
                       '',
                       '   \t ',
                       '\t\t\t',
                       ' ']

# missing data
BAD_DM_F1 = 'a\tb\na\t0\t1\nb\t1'

# mismatched IDs
BAD_DM_F2 = 'a\tb\nb\t0\t1\na\t1\t0'

# extra data lines
BAD_DM_F3 = '\ta\tb\na\t0\t1\nb\t1\t0\n  \nfoo\n\n\n'

# missing data lines
BAD_DM_F4 = '\ta\tb\na\t0\t1\n  \n'

# no data lines
BAD_DM_F5 = '\ta\tb\n'


if __name__ == '__main__':
    main()
