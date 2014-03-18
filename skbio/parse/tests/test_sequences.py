#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------


from skbio.core.exception import FastqParseError, RecordError
from skbio.parse.sequences import parse_fastq
from skbio.parse.sequences import parse_fasta
from unittest import TestCase, main


class GenericFastaTest(TestCase):
    """Setup data for all the various FASTA parsers."""

    def setUp(self):
        """standard files"""
        self.labels = '>abc\n>def\n>ghi\n'.split('\n')
        self.oneseq = '>abc\nUCAG\n'.split('\n')
        self.multiline = '>xyz\nUUUU\nCC\nAAAAA\nG'.split('\n')
        self.threeseq = '>123\na\n> \t abc  \t \ncag\ngac\n>456\nc\ng'.split(
            '\n')
        self.twogood = '>123\n\n> \t abc  \t \ncag\ngac\n>456\nc\ng'.split(
            '\n')
        self.oneX = '>123\nX\n> \t abc  \t \ncag\ngac\n>456\nc\ng'.split('\n')
        self.nolabels = 'GJ>DSJGSJDF\nSFHKLDFS>jkfs\n'.split('\n')
        self.empty = []


class parse_fastaTests(GenericFastaTest):
    """Tests of parse_fasta: returns (label, seq) tuples."""

    def test_empty(self):
        """parse_fasta should return empty list from 'file' w/o labels
        """
        self.assertEqual(list(parse_fasta(self.empty)), [])
        self.assertEqual(list(parse_fasta(self.nolabels, strict=False)),
                         [])
        self.assertRaises(RecordError, list, parse_fasta(self.nolabels))

    def test_no_labels(self):
        """parse_fasta should return empty list from file w/o seqs"""
        # should fail if strict (the default)
        self.assertRaises(RecordError, list,
                          parse_fasta(self.labels, strict=True))
        # if not strict, should skip the records
        self.assertEqual(list(parse_fasta(self.labels, strict=False)),
                         [])

    def test_single(self):
        """parse_fasta should read single record as (label, seq) tuple
        """
        f = list(parse_fasta(self.oneseq))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('abc', 'UCAG'))

        f = list(parse_fasta(self.multiline))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('xyz', 'UUUUCCAAAAAG'))

    def test_gt_bracket_in_seq(self):
        """parse_fasta handles alternate finder function

            this test also illustrates how to use the parse_fasta
            to handle "sequences" that start with a > symbol, which can
            happen when we abuse the parse_fasta to parse
            fasta-like sequence quality files.
        """
        oneseq_w_gt = '>abc\n>CAG\n'.split('\n')

        def get_two_line_records(infile):
            line1 = None
            for line in infile:
                if line1 is None:
                    line1 = line
                else:
                    yield (line1, line)
                    line1 = None
        f = list(parse_fasta(oneseq_w_gt, finder=get_two_line_records))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('abc', '>CAG'))

    def test_multiple(self):
        """parse_fasta should read multiline records correctly"""
        f = list(parse_fasta(self.threeseq))
        self.assertEqual(len(f), 3)
        a, b, c = f
        self.assertEqual(a, ('123', 'a'))
        self.assertEqual(b, ('abc', 'caggac'))
        self.assertEqual(c, ('456', 'cg'))

    def test_multiple_bad(self):
        """parse_fasta should complain or skip bad records"""
        self.assertRaises(RecordError, list, parse_fasta(self.twogood))
        f = list(parse_fasta(self.twogood, strict=False))
        self.assertEqual(len(f), 2)
        a, b = f
        self.assertEqual(a, ('abc', 'caggac'))


class ParseFastqTests(TestCase):
    def setUp(self):
        """ Initialize variables to be used by the tests """
        self.fastq_example = fastq_example.split('\n')
        self.fastq_example_2 = fastq_example_2.split('\n')

    def test_parse(self):
        """sequence and info objects should correctly match"""
        for label, seq, qual in parse_fastq(self.fastq_example):
            self.assertTrue(label in data)
            self.assertEqual(seq, data[label]["seq"])
            self.assertEqual(qual, data[label]["qual"])

    def test_parse_error(self):
        """Does this raise a FastqParseError with incorrect input?"""
        with self.assertRaises(FastqParseError):
            list(parse_fastq(self.fastq_example_2, strict=True))

data = {
    "GAPC_0015:6:1:1259:10413#0/1":
    dict(seq='AACACCAAACTTCTCCACCACGTGAGCTACAAAAG',
         qual=r'````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF'),
    "GAPC_0015:6:1:1283:11957#0/1":
    dict(seq='TATGTATATATAACATATACATATATACATACATA',
         qual=r']KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb'),
    "GAPC_0015:6:1:1284:10484#0/1":
    dict(seq='TCAGTTTTCCTCGCCATATTTCACGTCCTAAAGCG',
         qual=r'UM_]]U_]Z_Y^\^^``Y]`^SZ]\Ybb`^_LbL_'),
    "GAPC_0015:6:1:1287:17135#0/1":
    dict(seq='TGTGCCTATGGAAGCAGTTCTAGGATCCCCTAGAA',
         qual=r'^aacccL\ccc\c\cTKTS]KZ\]]I\[Wa^T`^K'),
    "GAPC_0015:6:1:1293:3171#0/1":
    dict(seq="AAAGAAAGGAAGAAAAGAAAAAGAAACCCGAGTTA",
             qual=r"b`bbbU_[YYcadcda_LbaaabWbaacYcc`a^c"),
    "GAPC_0015:6:1:1297:10729#0/1":
    dict(seq="TAATGCCAAAGAAATATTTCCAAACTACATGCTTA",
             qual=r"T\ccLbb``bacc]_cacccccLccc\ccTccYL^"),
    "GAPC_0015:6:1:1299:5940#0/1":
    dict(seq="AATCAAGAAATGAAGATTTATGTATGTGAAGAATA",
             qual=r"dcddbcfffdfffd`dd`^`c`Oc`Ybb`^eecde"),
    "GAPC_0015:6:1:1308:6996#0/1":
    dict(seq="TGGGACACATGTCCATGCTGTGGTTTTAACCGGCA",
             qual=r"a]`aLY`Y^^ccYa`^^TccK_X]\c\c`caTTTc"),
    "GAPC_0015:6:1:1314:13295#0/1":
    dict(seq="AATATTGCTTTGTCTGAACGATAGTGCTCTTTGAT",
         qual=r"cLcc\\dddddaaYd`T```bLYT\`a```bZccc"),
    "GAPC_0015:6:1:1317:3403#0/1":
    dict(seq="TTGTTTCCACTTGGTTGATTTCACCCCTGAGTTTG",
         # had to add space in qual line
         qual=r"\\\ZTYTSaLbb``\_UZ_bbcc`cc^[ac\a\Tc ".strip())
}

fastq_example = r"""@GAPC_0015:6:1:1259:10413#0/1
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

fastq_example_2 = r"""@GAPC_0017:6:1:1259:10413#0/1
AACACCAAACTTCTCCACCACGTGAGCTACAAAAG
+GAPC_0015:6:1:1259:10413#0/1
````Y^T]`]c^cabcacc`^Lb^ccYT\T\Y\WF
@GAPC_0015:6:1:1283:11957#0/1
TATGTATATATAACATATACATATATACATACATA
+GAPC_0015:6:1:1283:11957#0/1
]KZ[PY]_[YY^```ac^\\`bT``c`\aT``bbb
@GAPC_0015:6:1:1284:10484#0/1
"""

if __name__ == "__main__":
    main()
