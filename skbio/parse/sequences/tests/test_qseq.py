#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function

import tempfile
from unittest import TestCase, main

import numpy.testing as npt

from skbio import parse_qseq
from skbio.core.exception import RecordError


QSEQ_PARSERS_DATA = {
    'one_seq': 'CRESSIA\t242\t1\t2204\t1491\t1930\t0\t1\t'
               '.TTGATAAGAATGTCTGTTGTGGCTTAAAA\t'
               'B[[[W][Y[Zcdccccccc\cccac_____\t1',
    'two_seq': 'CRESSIA\t242\t1\t2204\t1453\t1918\t0\t1\t'
               '.TTAATAAGAATGTCTGTTGTGGCTTAAAA\t'
               'B[[[W][Y[Zccccccccc\cccac_____\t1\n'
               'CRESSIA\t242\t1\t2204\t1490\t1921\t0\t2\t'
               '..GTAAAACCCATATATTGAAAACTACAAA\t'
               'BWUTWcXVXXcccc_cccccccccc_cccc\t0',
    'missing_items': 'CRESSIA\t242\t1\t1491\t0\t1\t'
                     '.TTGATAAGAATGTCTGTTGTGGCTTAAAA\t'
                     'B[[[W][Y[Zcdccccccc\cccac_____\t1',
    }


class IterableData(object):
    """Store qseq data as lists of strings."""
    def setUp(self):
        for attr, val in QSEQ_PARSERS_DATA.items():
            setattr(self, attr, val.split('\n'))


class FileData(object):
    """Store fasta data as file names pointing to the data."""
    def setUp(self):
        tmp_files = []
        for attr, val in QSEQ_PARSERS_DATA.items():
            tmp_file = tempfile.NamedTemporaryFile('r+')
            tmp_file.write(val)
            tmp_file.flush()
            tmp_file.seek(0)
            setattr(self, attr, tmp_file.name)
            tmp_files.append(tmp_file)
        self._tmp_files = tmp_files

    def tearDown(self):
        for tmp_file in self._tmp_files:
            tmp_file.close()


class ParseQseqTests(object):

    """Tests of parse_qseq: returns (seq_id, seq, qual, others) tuples."""

    def test_one_seq(self):
        """parse_qseq should return one record as tuple.
        """
        f = list(parse_qseq(self.one_seq))
        self.assertEqual(len(f), 1)
        a = f[0]
        # First record.
        self.assertEqual(a[0], 'CRESSIA_242:1:2204:1491:1930#0/1')
        self.assertEqual(a[1], '.TTGATAAGAATGTCTGTTGTGGCTTAAAA')
        npt.assert_equal(a[2], [2, 27, 27, 27, 23, 29, 27, 25, 27, 26, 35, 36,
                                35, 35, 35, 35, 35, 35, 35, 28, 35, 35, 35, 33,
                                35, 31, 31, 31, 31, 31])
        record = a[3]
        self.assertEqual(record.machine_name, 'CRESSIA')
        self.assertEqual(record.run, 242)
        self.assertEqual(record.lane, 1)
        self.assertEqual(record.tile, 2204)
        self.assertEqual(record.x, 1491)
        self.assertEqual(record.y, 1930)
        self.assertEqual(record.index, 0)
        self.assertEqual(record.read, 1)
        self.assertTrue(record.filtered)

    def test_two_seq(self):
        """parse_qseq should return two records as tuple.
        """
        f = list(parse_qseq(self.two_seq))
        self.assertEqual(len(f), 2)
        a, b = f
        # First record.
        self.assertEqual(a[0], 'CRESSIA_242:1:2204:1453:1918#0/1')
        self.assertEqual(a[1], '.TTAATAAGAATGTCTGTTGTGGCTTAAAA')
        npt.assert_equal(a[2], [2, 27, 27, 27, 23, 29, 27, 25, 27, 26, 35, 35,
                                35, 35, 35, 35, 35, 35, 35, 28, 35, 35, 35, 33,
                                35, 31, 31, 31, 31, 31])
        record = a[3]
        self.assertEqual(record.machine_name, 'CRESSIA')
        self.assertEqual(record.run, 242)
        self.assertEqual(record.lane, 1)
        self.assertEqual(record.tile, 2204)
        self.assertEqual(record.x, 1453)
        self.assertEqual(record.y, 1918)
        self.assertEqual(record.index, 0)
        self.assertEqual(record.read, 1)
        self.assertTrue(record.filtered)
        # Second record.
        self.assertEqual(b[0], 'CRESSIA_242:1:2204:1490:1921#0/2')
        self.assertEqual(b[1], '..GTAAAACCCATATATTGAAAACTACAAA')
        npt.assert_equal(b[2], [2, 23, 21, 20, 23, 35, 24, 22, 24, 24, 35, 35,
                                35, 35, 31, 35, 35, 35, 35, 35, 35, 35, 35, 35,
                                35, 31, 35, 35, 35, 35])
        record = b[3]
        self.assertEqual(record.machine_name, 'CRESSIA')
        self.assertEqual(record.run, 242)
        self.assertEqual(record.lane, 1)
        self.assertEqual(record.tile, 2204)
        self.assertEqual(record.x, 1490)
        self.assertEqual(record.y, 1921)
        self.assertEqual(record.index, 0)
        self.assertEqual(record.read, 2)
        self.assertFalse(record.filtered)

    def test_missing_items(self):
        """parse_qseq should raise RecordError.
        """
        self.assertRaises(RecordError, list, parse_qseq(self.missing_items))


class ParseQseqTestsInputIsIterable(IterableData, ParseQseqTests, TestCase):
    """Mixin: `parse_qseq` in ParseQseqTests gets lists
    of strings.
    """
    pass


class ParseFastaTestsInputIsFileNames(FileData, ParseQseqTests, TestCase):
    """Mixin: `parse_qseq` in ParseQseqTests gets a
    file name.
    """
    pass


if __name__ == "__main__":
    main()
