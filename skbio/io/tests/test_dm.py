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

from skbio.io import DMFormatError
from skbio.io.dm import (dm_to_DissimilarityMatrix, dm_to_DistanceMatrix,
                         DissimilarityMatrix_to_dm, DistanceMatrix_to_dm,
                         dm_sniffer)
from skbio.stats.distance import (DissimilarityMatrix, DistanceMatrix,
                                  DistanceMatrixError)


class DissimilarityAndDistanceMatrixReaderWriterTests(TestCase):
    def setUp(self):
        self.dm_1x1_data = [[0.0]]
        self.dm_1x1_fh = StringIO(DM_1x1)

        self.dm_2x2_data = [[0.0, 0.123], [0.123, 0.0]]
        self.dm_2x2_fh = StringIO(DM_2x2)

        self.dm_2x2_asym_data = [[0.0, 1.0], [-2.0, 0.0]]
        self.dm_2x2_asym_fh = StringIO(DM_2x2_ASYM)

        self.dm_3x3_data = [[0.0, 0.01, 4.2], [0.01, 0.0, 12.0],
                            [4.2, 12.0, 0.0]]
        self.dm_3x3_fh = StringIO(DM_3x3)
        self.dm_3x3_whitespace_fh = StringIO(DM_3x3_WHITESPACE)

        self.invalid_fhs = [
            StringIO(),
            StringIO(INVALID_1),
            StringIO(INVALID_2),
            StringIO(INVALID_3),
            StringIO(INVALID_4),
            StringIO(INVALID_5)
        ]

        # We repeat the 3x3 example because there are two file format
        # representations of it, one that is messy and one that is not. Both
        # should be read into an equivalent object and written to an equivalent
        # format though, which is why we duplicate the 3x3 objects and strings.
        self.dissim_objs = [
            DissimilarityMatrix(self.dm_1x1_data, ['a']),
            DissimilarityMatrix(self.dm_2x2_data, ['a', 'b']),
            DissimilarityMatrix(self.dm_2x2_asym_data, ['a', 'b']),
            DissimilarityMatrix(self.dm_3x3_data, ['a', 'b', 'c']),
            DissimilarityMatrix(self.dm_3x3_data, ['a', 'b', 'c'])
        ]

        self.dissim_strs = [DM_1x1, DM_2x2, DM_2x2_ASYM, DM_3x3, DM_3x3]

        self.dissim_fhs = [self.dm_1x1_fh, self.dm_2x2_fh, self.dm_2x2_asym_fh,
                           self.dm_3x3_fh, self.dm_3x3_whitespace_fh]

        self.dist_objs = [
            DistanceMatrix(self.dm_1x1_data, ['a']),
            DistanceMatrix(self.dm_2x2_data, ['a', 'b']),
            DistanceMatrix(self.dm_3x3_data, ['a', 'b', 'c']),
            DistanceMatrix(self.dm_3x3_data, ['a', 'b', 'c'])
        ]

        self.dist_strs = [DM_1x1, DM_2x2, DM_3x3, DM_3x3]

        self.dist_fhs = [self.dm_1x1_fh, self.dm_2x2_fh, self.dm_3x3_fh,
                         self.dm_3x3_whitespace_fh]

    def test_read_valid_files(self):
        for fn, cls, objs, fhs in ((dm_to_DissimilarityMatrix,
                                    DissimilarityMatrix, self.dissim_objs,
                                    self.dissim_fhs),
                                   (dm_to_DistanceMatrix, DistanceMatrix,
                                    self.dist_objs, self.dist_fhs)):
            for fh, obj in zip(fhs, objs):
                obs = fn(fh)
                self.assertEqual(obs, obj)
                self.assertIsInstance(obs, cls)

    def test_read_invalid_files(self):
        for fn in dm_to_DissimilarityMatrix, dm_to_DistanceMatrix:
            for invalid_fh in self.invalid_fhs:
                with self.assertRaises(DMFormatError):
                    fn(invalid_fh)

        # Asymmetric data only raises an error for DistanceMatrix.
        with self.assertRaises(DistanceMatrixError):
            dm_to_DistanceMatrix(self.dm_2x2_asym_fh)

    def test_write(self):
        for fn, objs, strs in ((DissimilarityMatrix_to_dm, self.dissim_objs,
                                self.dissim_strs),
                               (DistanceMatrix_to_dm, self.dist_objs,
                                self.dist_strs)):
            for obj, str_ in zip(objs, strs):
                fh = StringIO()
                fn(obj, fh)
                obs = fh.getvalue()
                fh.close()
                self.assertEqual(obs, str_)

    def test_roundtrip_read_write(self):
        for reader_fn, writer_fn, fhs in ((dm_to_DissimilarityMatrix,
                                           DissimilarityMatrix_to_dm,
                                           self.dissim_fhs),
                                          (dm_to_DistanceMatrix,
                                           DistanceMatrix_to_dm,
                                           self.dist_fhs)):
            for fh in fhs:
                # Read.
                dm1 = reader_fn(fh)

                # Write.
                out_fh = StringIO()
                writer_fn(dm1, out_fh)
                out_fh.seek(0)

                # Read.
                dm2 = reader_fn(out_fh)

                self.assertEqual(dm1, dm2)


class SnifferTests(TestCase):
    def test_valid(self):
        obs = dm_sniffer(StringIO(DM_3x3))
        self.assertEqual(obs, (True, {'delimiter': '\t'}))


DM_1x1 = "\ta\na\t0.0\n"

DM_2x2 = "\ta\tb\na\t0.0\t0.123\nb\t0.123\t0.0\n"

DM_2x2_ASYM = "\ta\tb\na\t0.0\t1.0\nb\t-2.0\t0.0\n"

DM_3x3 = ("\ta\tb\tc\na\t0.0\t0.01\t4.2\nb\t0.01\t0.0\t12.0\nc\t4.2\t12.0\t"
          "0.0\n")

# Extra whitespace-only lines throughout. Also has comments before the header.
DM_3x3_WHITESPACE = '\n'.join(['# foo',
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
                               ' '])

# missing data
INVALID_1 = 'a\tb\na\t0\t1\nb\t1'

# mismatched IDs
INVALID_2 = 'a\tb\nb\t0\t1\na\t1\t0'

# extra data lines
INVALID_3 = '\ta\tb\na\t0\t1\nb\t1\t0\n  \nfoo\n\n\n'

# missing data lines
INVALID_4 = '\ta\tb\na\t0\t1\n  \n'

# no data lines
INVALID_5 = '\ta\tb\n'


if __name__ == '__main__':
    main()
