# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from unittest import TestCase, main

from skbio import DistanceMatrix
from skbio.io import LSMatFormatError
from skbio.io.format.lsmat import (
    _lsmat_to_dissimilarity_matrix, _lsmat_to_distance_matrix,
    _dissimilarity_matrix_to_lsmat, _distance_matrix_to_lsmat, _lsmat_sniffer)
from skbio.stats.distance import DissimilarityMatrix, DistanceMatrixError


class LSMatTestData(TestCase):
    def setUp(self):
        self.lsmat_1x1_fh = io.StringIO(LSMat_1x1)
        self.lsmat_2x2_fh = io.StringIO(LSMat_2x2)
        self.lsmat_2x2_asym_fh = io.StringIO(LSMat_2x2_ASYM)
        self.lsmat_3x3_fh = io.StringIO(LSMat_3x3)
        self.lsmat_3x3_whitespace_fh = io.StringIO(LSMat_3x3_WHITESPACE)
        self.lsmat_3x3_csv_fh = io.StringIO(LSMat_3x3_CSV)
        self.lsmat_3x3_fw_fh = io.StringIO(LSMat_3x3_FW)

        self.valid_fhs = [
            self.lsmat_1x1_fh,
            self.lsmat_2x2_fh,
            self.lsmat_2x2_asym_fh,
            self.lsmat_3x3_fh,
            self.lsmat_3x3_whitespace_fh
        ]

        self.empty_fh = io.StringIO()
        self.invalid_1_fh = io.StringIO(INVALID_1)
        self.invalid_2_fh = io.StringIO(INVALID_2)
        self.invalid_3_fh = io.StringIO(INVALID_3)
        self.invalid_4_fh = io.StringIO(INVALID_4)
        self.invalid_5_fh = io.StringIO(INVALID_5)
        self.invalid_6_fh = io.StringIO(INVALID_6)

        self.invalid_fhs = [
            (self.empty_fh, r'empty'),
            (self.invalid_1_fh, r'1 value\(s\).*2.*\(2\)'),
            (self.invalid_2_fh, r"'b'.*'a'"),
            (self.invalid_3_fh, r'extra row\(s\)'),
            (self.invalid_4_fh, r'2 row\(s\).*found 1'),
            (self.invalid_5_fh, r'2 row\(s\).*found 0'),
            (self.invalid_6_fh, r"delimiter '\\t'")
        ]


class DissimilarityAndDistanceMatrixReaderWriterTests(LSMatTestData):
    def setUp(self):
        super(DissimilarityAndDistanceMatrixReaderWriterTests, self).setUp()

        self.lsmat_1x1_data = [[0.0]]
        self.lsmat_2x2_data = [[0.0, 0.123], [0.123, 0.0]]
        self.lsmat_2x2_asym_data = [[0.0, 1.0], [-2.0, 0.0]]
        self.lsmat_3x3_data = [[0.0, 0.01, 4.2], [0.01, 0.0, 12.0],
                               [4.2, 12.0, 0.0]]

        # We repeat the 3x3 example because there are two file format
        # representations of it, one that is messy and one that is not. Both
        # should be read into an equivalent object and written to an equivalent
        # format though, which is why we duplicate the 3x3 objects and strings.
        self.dissim_objs = [
            DissimilarityMatrix(self.lsmat_1x1_data, ['a']),
            DissimilarityMatrix(self.lsmat_2x2_data, ['a', 'b']),
            DissimilarityMatrix(self.lsmat_2x2_asym_data, ['a', 'b']),
            DissimilarityMatrix(self.lsmat_3x3_data, ['a', 'b', 'c']),
            DissimilarityMatrix(self.lsmat_3x3_data, ['a', 'b', 'c'])
        ]

        self.dissim_strs = [LSMat_1x1, LSMat_2x2, LSMat_2x2_ASYM, LSMat_3x3,
                            LSMat_3x3]

        self.dissim_fhs = [self.lsmat_1x1_fh, self.lsmat_2x2_fh,
                           self.lsmat_2x2_asym_fh, self.lsmat_3x3_fh,
                           self.lsmat_3x3_whitespace_fh]

        self.dist_objs = [
            DistanceMatrix(self.lsmat_1x1_data, ['a']),
            DistanceMatrix(self.lsmat_2x2_data, ['a', 'b']),
            DistanceMatrix(self.lsmat_3x3_data, ['a', 'b', 'c']),
            DistanceMatrix(self.lsmat_3x3_data, ['a', 'b', 'c'])
        ]

        self.dist_strs = [LSMat_1x1, LSMat_2x2, LSMat_3x3, LSMat_3x3]

        self.dist_fhs = [self.lsmat_1x1_fh, self.lsmat_2x2_fh,
                         self.lsmat_3x3_fh, self.lsmat_3x3_whitespace_fh]

    def test_read_valid_files(self):
        for fn, cls, objs, fhs in ((_lsmat_to_dissimilarity_matrix,
                                    DissimilarityMatrix, self.dissim_objs,
                                    self.dissim_fhs),
                                   (_lsmat_to_distance_matrix, DistanceMatrix,
                                    self.dist_objs, self.dist_fhs)):
            for fh, obj in zip(fhs, objs):
                fh.seek(0)
                obs = fn(fh)
                self.assertEqual(obs, obj)
                self.assertIsInstance(obs, cls)

        # Above files are TSV (default delimiter). Test that CSV works too.
        for fn, cls in ((_lsmat_to_dissimilarity_matrix, DissimilarityMatrix),
                        (_lsmat_to_distance_matrix, DistanceMatrix)):
            exp = cls(self.lsmat_3x3_data, ['a', 'b', 'c'])
            self.lsmat_3x3_csv_fh.seek(0)
            obs = fn(self.lsmat_3x3_csv_fh, delimiter=',')
            self.assertEqual(obs, exp)
            self.assertIsInstance(obs, cls)

        # Test that fixed-width works too.
        for fn, cls in ((_lsmat_to_dissimilarity_matrix, DissimilarityMatrix),
                        (_lsmat_to_distance_matrix, DistanceMatrix)):
            exp = cls(self.lsmat_3x3_data, ['a', 'b', 'c'])
            self.lsmat_3x3_fw_fh.seek(0)
            obs = fn(self.lsmat_3x3_fw_fh, delimiter=None)
            self.assertEqual(obs, exp)
            self.assertIsInstance(obs, cls)

    def test_read_invalid_files(self):
        for fn in _lsmat_to_dissimilarity_matrix, _lsmat_to_distance_matrix:
            for invalid_fh, error_msg_regexp in self.invalid_fhs:
                with self.assertRaisesRegex(LSMatFormatError,
                                            error_msg_regexp):
                    invalid_fh.seek(0)
                    fn(invalid_fh)

        # Asymmetric data only raises an error for DistanceMatrix.
        with self.assertRaises(DistanceMatrixError):
            _lsmat_to_distance_matrix(self.lsmat_2x2_asym_fh)

    def test_write(self):
        for fn, objs, strs in ((_dissimilarity_matrix_to_lsmat,
                                self.dissim_objs, self.dissim_strs),
                               (_distance_matrix_to_lsmat, self.dist_objs,
                                self.dist_strs)):
            for obj, str_ in zip(objs, strs):
                fh = io.StringIO()
                fn(obj, fh)
                obs = fh.getvalue()
                fh.close()
                self.assertEqual(obs, str_)

        # Test writing CSV (TSV is written above).
        for fn, cls in ((_dissimilarity_matrix_to_lsmat, DissimilarityMatrix),
                        (_distance_matrix_to_lsmat, DistanceMatrix)):
            obj = cls(self.lsmat_3x3_data, ['a', 'b', 'c'])
            fh = io.StringIO()
            fn(obj, fh, delimiter=',')
            obs = fh.getvalue()
            fh.close()
            self.assertEqual(obs, LSMat_3x3_CSV)

    def test_roundtrip_read_write(self):
        for reader_fn, writer_fn, fhs in ((_lsmat_to_dissimilarity_matrix,
                                           _dissimilarity_matrix_to_lsmat,
                                           self.dissim_fhs),
                                          (_lsmat_to_distance_matrix,
                                           _distance_matrix_to_lsmat,
                                           self.dist_fhs)):
            for fh in fhs:
                # Read.
                fh.seek(0)
                lsmat1 = reader_fn(fh)

                # Write.
                out_fh = io.StringIO()
                writer_fn(lsmat1, out_fh)
                out_fh.seek(0)

                # Read.
                lsmat2 = reader_fn(out_fh)
                out_fh.close()

                self.assertEqual(lsmat1, lsmat2)


class SnifferTests(LSMatTestData):
    def setUp(self):
        super(SnifferTests, self).setUp()

    def test_match_tsv(self):
        # Sniffer should match all valid files, and will match some invalid
        # ones too because it doesn't exhaustively check the entire file.
        fhs = self.valid_fhs + [self.invalid_1_fh, self.invalid_3_fh,
                                self.invalid_4_fh]
        for fh in fhs:
            self.assertEqual(_lsmat_sniffer(fh), (True, {'delimiter': '\t'}))

    def test_match_csv(self):
        self.assertEqual(_lsmat_sniffer(self.lsmat_3x3_csv_fh),
                         (True, {'delimiter': ','}))

    def test_no_match(self):
        for fh in (self.empty_fh, self.invalid_2_fh, self.invalid_5_fh,
                   self.invalid_6_fh):
            self.assertEqual(_lsmat_sniffer(fh), (False, {}))


LSMat_1x1 = (
    '\ta\n'
    'a\t0.0\n')

LSMat_2x2 = (
    '\ta\tb\n'
    'a\t0.0\t0.123\n'
    'b\t0.123\t0.0\n')

LSMat_2x2_ASYM = (
    '\ta\tb\n'
    'a\t0.0\t1.0\n'
    'b\t-2.0\t0.0\n')

LSMat_3x3 = (
    '\ta\tb\tc\n'
    'a\t0.0\t0.01\t4.2\n'
    'b\t0.01\t0.0\t12.0\n'
    'c\t4.2\t12.0\t0.0\n')

# Extra whitespace-only lines throughout. Also has comments before the header.
LSMat_3x3_WHITESPACE = '\n'.join([
    '# foo',
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

# Same matrix as above, but delimited by commas instead of tabs.
LSMat_3x3_CSV = (
    ',a,b,c\n'
    'a,0.0,0.01,4.2\n'
    'b,0.01,0.0,12.0\n'
    'c,4.2,12.0,0.0\n')

# Same matrix as above, but delimited by whitespaces instead of tabs.
LSMat_3x3_FW = (
    '   a     b     c   \n'
    'a  0.0   0.01  4.2 \n'
    'b  0.01  0.0  12.0 \n'
    'c  4.2  12.0   0.0 \n')

# missing data
INVALID_1 = '\ta\tb\na\t0\t1\nb\t1'

# mismatched IDs
INVALID_2 = '\ta\tb\nb\t0\t1\na\t1\t0'

# extra data lines
INVALID_3 = '\ta\tb\na\t0\t1\nb\t1\t0\n  \nfoo\n\n\n'

# missing data lines
INVALID_4 = '\ta\tb\na\t0\t1\n  \n'

# no data lines
INVALID_5 = '\ta\tb\n'

# missing leading delimiter in header
INVALID_6 = "a\tb\na\t0.0\t0.123\nb\t0.123\t0.0\n"


if __name__ == '__main__':
    main()
