#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function

from unittest import TestCase, main
import tempfile

from numpy import array

from skbio.parse.sequences import parse_fastq
from skbio.core.exception import FastqParseError


class IterableData(object):
    def setUp(self):
        """ Initialize variables to be used by the tests as lists of strings"""
        self.FASTQ_EXAMPLE = FASTQ_EXAMPLE.split('\n')
        self.FASTQ_EXAMPLE_2 = FASTQ_EXAMPLE_2.split('\n')
        self.FASTQ_EXAMPLE_3 = FASTQ_EXAMPLE_3.split('\n')


class FileData(object):
    def setUp(self):
        """ Initialize variables to be used by the tests as file names"""
        tmp_files = []
        for attr, val in [('FASTQ_EXAMPLE', FASTQ_EXAMPLE),
                          ('FASTQ_EXAMPLE_2', FASTQ_EXAMPLE_2),
                          ('FASTQ_EXAMPLE_3', FASTQ_EXAMPLE_3)]:
            tmp_file = tempfile.NamedTemporaryFile('w')
            tmp_file.write(val)
            tmp_file.flush()
            tmp_file.seek(0)
            setattr(self, attr, tmp_file.name)
            tmp_files.append(tmp_file)
        self._tmp_files = tmp_files

    def tearDown(self):
        for tmp_file in self._tmp_files:
            tmp_file.close()


class ParseFastqTests(object):
    def test_parse(self):
        for label, seq, qual in parse_fastq(self.FASTQ_EXAMPLE,
                                            phred_offset=64):
            self.assertTrue(label in DATA)
            self.assertEqual(seq, DATA[label]["seq"])
            self.assertTrue((qual == DATA[label]["qual"]).all())

        # Make sure that enforce_qual_range set to False allows qual scores
        # to fall outside the typically acceptable range of 0-62
        for label, seq, qual in parse_fastq(self.FASTQ_EXAMPLE_2,
                                            phred_offset=33,
                                            enforce_qual_range=False):
            self.assertTrue(label in DATA_2)
            self.assertEqual(seq, DATA_2[label]["seq"])
            self.assertTrue((qual == DATA_2[label]["qual"]).all())

        # This should raise a FastqParseError since the qual scores are
        # intended to be interpreted with an offset of 64, and using 33 will
        # make the qual score fall outside the acceptable range of 0-62.
        with self.assertRaises(FastqParseError):
            list(parse_fastq(self.FASTQ_EXAMPLE, phred_offset=33))

    def test_parse_error(self):
        with self.assertRaises(FastqParseError):
            list(parse_fastq(self.FASTQ_EXAMPLE_2, strict=True))

        with self.assertRaises(FastqParseError):
            list(parse_fastq(self.FASTQ_EXAMPLE_3, phred_offset=64))

    def test_invalid_phred_offset(self):
        with self.assertRaises(ValueError):
            list(parse_fastq(self.FASTQ_EXAMPLE, phred_offset=42))


class ParseFastqTestsInputIsIterable(IterableData, ParseFastqTests, TestCase):
    pass


class ParseFastqTestsInputIsFileNames(FileData, ParseFastqTests, TestCase):
    pass


DATA = {
    "GAPC_0015:6:1:1259:10413#0/1":
    dict(seq='AACACCAAACTTCTCCACCACGTGAGCTACAAAAG',
         qual=array([32, 32, 32, 32, 25, 30, 20, 29, 32, 29, 35, 30, 35, 33,
                     34, 35, 33, 35, 35, 32, 30, 12, 34, 30, 35, 35, 25, 20,
                     28, 20, 28, 25, 28, 23, 6])),
    "GAPC_0015:6:1:1283:11957#0/1":
    dict(seq='TATGTATATATAACATATACATATATACATACATA',
         qual=array([29, 11, 26, 27, 16, 25, 29, 31, 27, 25, 25, 30, 32, 32,
                     32, 33, 35, 30, 28, 28, 32, 34, 20, 32, 32, 35, 32, 28,
                     33, 20, 32, 32, 34, 34, 34])),
    "GAPC_0015:6:1:1284:10484#0/1":
    dict(seq='TCAGTTTTCCTCGCCATATTTCACGTCCTAAAGCG',
         qual=array([21, 13, 31, 29, 29, 21, 31, 29, 26, 31, 25, 30, 28, 30,
                     30, 32, 32, 25, 29, 32, 30, 19, 26, 29, 28, 25, 34, 34,
                     32, 30, 31, 12, 34, 12, 31])),
    "GAPC_0015:6:1:1287:17135#0/1":
    dict(seq='TGTGCCTATGGAAGCAGTTCTAGGATCCCCTAGAA',
         qual=array([30, 33, 33, 35, 35, 35, 12, 28, 35, 35, 35, 28, 35, 28,
                     35, 20, 11, 20, 19, 29, 11, 26, 28, 29, 29,  9, 28, 27,
                     23, 33, 30, 20, 32, 30, 11])),
    "GAPC_0015:6:1:1293:3171#0/1":
    dict(seq="AAAGAAAGGAAGAAAAGAAAAAGAAACCCGAGTTA",
             qual=array([34, 32, 34, 34, 34, 21, 31, 27, 25, 25, 35, 33, 36,
                         35, 36, 33, 31, 12, 34, 33, 33, 33, 34, 23, 34, 33,
                         33, 35, 25, 35, 35, 32, 33, 30, 35])),
    "GAPC_0015:6:1:1297:10729#0/1":
    dict(seq="TAATGCCAAAGAAATATTTCCAAACTACATGCTTA",
             qual=array([20, 28, 35, 35, 12, 34, 34, 32, 32, 34, 33, 35, 35,
                         29, 31, 35, 33, 35, 35, 35, 35, 35, 12, 35, 35, 35,
                         28, 35, 35, 20, 35, 35, 25, 12, 30])),
    "GAPC_0015:6:1:1299:5940#0/1":
    dict(seq="AATCAAGAAATGAAGATTTATGTATGTGAAGAATA",
             qual=array([36, 35, 36, 36, 34, 35, 38, 38, 38, 36, 38, 38, 38,
                         36, 32, 36, 36, 32, 30, 32, 35, 32, 15, 35, 32, 25,
                         34, 34, 32, 30, 37, 37, 35, 36, 37])),
    "GAPC_0015:6:1:1308:6996#0/1":
    dict(seq="TGGGACACATGTCCATGCTGTGGTTTTAACCGGCA",
             qual=array([33, 29, 32, 33, 12, 25, 32, 25, 30, 30, 35, 35, 25,
                         33, 32, 30, 30, 20, 35, 35, 11, 31, 24, 29, 28, 35,
                         28, 35, 32, 35, 33, 20, 20, 20, 35])),
    "GAPC_0015:6:1:1314:13295#0/1":
    dict(seq="AATATTGCTTTGTCTGAACGATAGTGCTCTTTGAT",
         qual=array([35, 12, 35, 35, 28, 28, 36, 36, 36, 36, 36, 33, 33, 25,
                     36, 32, 20, 32, 32, 32, 34, 12, 25, 20, 28, 32, 33, 32,
                     32, 32, 34, 26, 35, 35, 35])),
    "GAPC_0015:6:1:1317:3403#0/1":
    dict(seq="TTGTTTCCACTTGGTTGATTTCACCCCTGAGTTTG",
         # had to add space in qual line
         qual=array([28, 28, 28, 26, 20, 25, 20, 19, 33, 12, 34, 34, 32, 32,
                     28, 31, 21, 26, 31, 34, 34, 35, 35, 32, 35, 35, 30, 27,
                     33, 35, 28, 33, 28, 20, 35]))

}


DATA_2 = {
    "GAPC_0017:6:1:1259:10413#0/1":
    dict(seq='AACACCAAACTTCTCCACCACGTGAGCTACAAAAG',
         qual=array([63, 63, 63, 63, 56, 61, 51, 60, 63, 60, 66, 61, 66, 64,
                     65, 66, 64, 66, 66, 63, 61, 43, 65, 61, 66, 66, 56, 51,
                     59, 51, 59, 56, 59, 54, 37])),
    "GAPC_0015:6:1:1283:11957#0/1":
    dict(seq='TATGTATATATAACATATACATATATACATACATA',
         qual=array([60, 42, 57, 58, 47, 56, 60, 62, 58, 56, 56, 61, 63, 63,
                     63, 64, 66, 61, 59, 59, 63, 65, 51, 63, 63, 66, 63, 59,
                     64, 51, 63, 63, 65, 65, 65]))
}


FASTQ_EXAMPLE = r"""@GAPC_0015:6:1:1259:10413#0/1
AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
+GAPC_0015:6:1:1259:10413#0/1
````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
@GAPC_0015:6:1:1283:11957#0/1
TATGTATATATAACATATACATATATACATACATA
+GAPC_0015:6:1:1283:11957#0/1
]KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb
@GAPC_0015:6:1:1284:10484#0/1
TCAGTTTTCCTCGCCATATTTCACGTCCTAAAGCG
+GAPC_0015:6:1:1284:10484#0/1
UM_]]U_]Z_Y^\^^``Y]`^SZ]\Ybb`^_LbL_
@GAPC_0015:6:1:1287:17135#0/1
TGTGCCTATGGAAGCAGTTCTAGGATCCCCTAGAA
+GAPC_0015:6:1:1287:17135#0/1
^aacccL\ccc\c\cTKTS]KZ\]]I\[Wa^T`^K
@GAPC_0015:6:1:1293:3171#0/1
AAAGAAAGGAAGAAAAGAAAAAGAAACCCGAGTTA
+GAPC_0015:6:1:1293:3171#0/1
b`bbbU_[YYcadcda_LbaaabWbaacYcc`a^c
@GAPC_0015:6:1:1297:10729#0/1
TAATGCCAAAGAAATATTTCCAAACTACATGCTTA
+GAPC_0015:6:1:1297:10729#0/1
T\ccLbb``bacc]_cacccccLccc\ccTccYL^
@GAPC_0015:6:1:1299:5940#0/1
AATCAAGAAATGAAGATTTATGTATGTGAAGAATA
+GAPC_0015:6:1:1299:5940#0/1
dcddbcfffdfffd`dd`^`c`Oc`Ybb`^eecde
@GAPC_0015:6:1:1308:6996#0/1
TGGGACACATGTCCATGCTGTGGTTTTAACCGGCA
+GAPC_0015:6:1:1308:6996#0/1
a]`aLY`Y^^ccYa`^^TccK_X]\c\c`caTTTc
@GAPC_0015:6:1:1314:13295#0/1
AATATTGCTTTGTCTGAACGATAGTGCTCTTTGAT
+GAPC_0015:6:1:1314:13295#0/1
cLcc\\dddddaaYd`T```bLYT\`a```bZccc
@GAPC_0015:6:1:1317:3403#0/1
TTGTTTCCACTTGGTTGATTTCACCCCTGAGTTTG
+GAPC_0015:6:1:1317:3403#0/1
\\\ZTYTSaLbb``\_UZ_bbcc`cc^[ac\a\Tc"""


FASTQ_EXAMPLE_2 = r"""@GAPC_0017:6:1:1259:10413#0/1
AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
+GAPC_0015:6:1:1259:10413#0/1
````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
@GAPC_0015:6:1:1283:11957#0/1
TATGTATATATAACATATACATATATACATACATA
+GAPC_0015:6:1:1283:11957#0/1
]KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb
"""


FASTQ_EXAMPLE_3 = r"""@GAPC_0017:6:1:1259:10413#0/1
AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
+GAPC_0015:6:1:1259:10413#0/1
````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
@GAPC_0015:6:1:1283:11957#0/1
"""


if __name__ == "__main__":
    main()
