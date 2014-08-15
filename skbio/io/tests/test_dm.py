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
from skbio.io.dm import (dm_to_DissimilarityMatrix, dm_to_DistanceMatrix,
                         _dm_to_matrix)
from skbio.stats.distance import (DissimilarityMatrix, DistanceMatrix,
                                  DistanceMatrixError)
from skbio.util import get_data_path


class DMToMatrixTests(TestCase):
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

        self.invalid_fhs = [
            StringIO(),
            StringIO(BAD_DM_F1),
            StringIO(BAD_DM_F2),
            StringIO(BAD_DM_F3),
            StringIO(BAD_DM_F4),
            StringIO(BAD_DM_F5)
        ]

        self.dissim_objs = [
            DissimilarityMatrix(self.dm_1x1_data, ['a']),
            DissimilarityMatrix(self.dm_2x2_data, ['a', 'b']),
            DissimilarityMatrix(self.dm_2x2_asym_data, ['a', 'b']),
            DissimilarityMatrix(self.dm_3x3_data, ['a', 'b', 'c']),
            DissimilarityMatrix(self.dm_3x3_data, ['a', 'b', 'c'])
        ]

        self.dissim_fhs = [self.dm_1x1_f, self.dm_2x2_f, self.dm_2x2_asym_f,
                           self.dm_3x3_f, self.dm_3x3_whitespace_f]

        self.dist_objs = [
            DistanceMatrix(self.dm_1x1_data, ['a']),
            DistanceMatrix(self.dm_2x2_data, ['a', 'b']),
            DistanceMatrix(self.dm_3x3_data, ['a', 'b', 'c']),
            DistanceMatrix(self.dm_3x3_data, ['a', 'b', 'c'])
        ]

        self.dist_fhs = [self.dm_1x1_f, self.dm_2x2_f, self.dm_3x3_f,
                         self.dm_3x3_whitespace_f]

    def test_valid_files(self):
        for fn, cls, objs, fhs in (dm_to_DissimilarityMatrix, DissimilarityMatrix,
                self.dissim_objs, self.dissim_fhs), (dm_to_DistanceMatrix, DistanceMatrix, self.dist_objs, self.dist_fhs):
            for fh, obj in zip(fhs, objs):
                obs = fn(fh)
                self.assertEqual(obs, obj)
                self.assertIsInstance(obs, cls)

    def test_invalid_files(self):
        for fn in dm_to_DissimilarityMatrix, dm_to_DistanceMatrix:
            for invalid_fh in self.invalid_fhs:
                with self.assertRaises(DMFormatError):
                    fn(invalid_fh)

        # Asymmetric data only raises an error for DistanceMatrix.
        with self.assertRaises(DistanceMatrixError):
            dm_to_DistanceMatrix(self.dm_2x2_asym_f)


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
