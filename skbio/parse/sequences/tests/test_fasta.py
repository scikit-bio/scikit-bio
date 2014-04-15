#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import division

from skbio.parse.sequences.fasta import parse_fasta, parse_qual
from skbio.core.exception import RecordError

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


class ParseFastaTests(GenericFastaTest):

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

    def test_parse_qual(self):
        """parse_qual should yield (id_, quals)"""
        scores = ['>x', '5 10 5', '12',
                  '>y', '30 40',
                  '>a', '5 10 5', '12',
                  '>b', '30 40']
        gen = list(parse_qual(scores))
        self.assertItemsEqual(gen[0][1], [5, 10, 5, 12])
        self.assertItemsEqual(gen[1][1], [30, 40])
        self.assertItemsEqual(gen[2][1], [5, 10, 5, 12])
        self.assertItemsEqual(gen[3][1], [30, 40])

if __name__ == "__main__":
    main()
