#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from bipy.parse.fasta import MinimalFastaParser
from bipy.parse.record import RecordError
from bipy.util.unit_test import TestCase, main

class GenericFastaTest(TestCase):
    """Setup data for all the various FASTA parsers."""
    def setUp(self):
        """standard files"""
        self.labels = '>abc\n>def\n>ghi\n'.split('\n')
        self.oneseq = '>abc\nUCAG\n'.split('\n')
        self.multiline = '>xyz\nUUUU\nCC\nAAAAA\nG'.split('\n')
        self.threeseq='>123\na\n> \t abc  \t \ncag\ngac\n>456\nc\ng'.split('\n')
        self.twogood='>123\n\n> \t abc  \t \ncag\ngac\n>456\nc\ng'.split('\n')
        self.oneX='>123\nX\n> \t abc  \t \ncag\ngac\n>456\nc\ng'.split('\n')
        self.nolabels = 'GJ>DSJGSJDF\nSFHKLDFS>jkfs\n'.split('\n')
        self.empty = []
 
class MinimalFastaParserTests(GenericFastaTest):
    """Tests of MinimalFastaParser: returns (label, seq) tuples."""
       
    def test_empty(self):
        """MinimalFastaParser should return empty list from 'file' w/o labels"""
        self.assertEqual(list(MinimalFastaParser(self.empty)), [])
        self.assertEqual(list(MinimalFastaParser(self.nolabels, strict=False)),
            [])
        self.assertRaises(RecordError, list, MinimalFastaParser(self.nolabels))

    def test_no_labels(self):
        """MinimalFastaParser should return empty list from file w/o seqs"""
        #should fail if strict (the default)
        self.assertRaises(RecordError, list, 
            MinimalFastaParser(self.labels,strict=True))
        #if not strict, should skip the records
        self.assertEqual(list(MinimalFastaParser(self.labels, strict=False)), 
            [])
        
    def test_single(self):
        """MinimalFastaParser should read single record as (label, seq) tuple"""
        f = list(MinimalFastaParser(self.oneseq))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('abc', 'UCAG'))

        f = list(MinimalFastaParser(self.multiline))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('xyz', 'UUUUCCAAAAAG'))
    
    def test_gt_bracket_in_seq(self):
        """MinimalFastaParser handles alternate finder function
            
            this test also illustrates how to use the MinimalFastaParser
            to handle "sequences" that start with a > symbol, which can
            happen when we abuse the MinimalFastaParser to parse
            fasta-like sequence quality files.
        """
        oneseq_w_gt = '>abc\n>CAG\n'.split('\n')
        def get_two_line_records(infile):
            line1 = None
            for line in infile:
                if line1 == None:
                    line1 = line
                else:
                    yield (line1, line)
                    line1 = None
        f = list(MinimalFastaParser(oneseq_w_gt,finder=get_two_line_records))
        self.assertEqual(len(f), 1)
        a = f[0]
        self.assertEqual(a, ('abc', '>CAG'))

    def test_multiple(self):
        """MinimalFastaParser should read multiline records correctly"""
        f = list(MinimalFastaParser(self.threeseq))
        self.assertEqual(len(f), 3)
        a, b, c = f
        self.assertEqual(a, ('123', 'a'))
        self.assertEqual(b, ('abc', 'caggac'))
        self.assertEqual(c, ('456', 'cg'))

    def test_multiple_bad(self):
        """MinimalFastaParser should complain or skip bad records"""
        self.assertRaises(RecordError, list, MinimalFastaParser(self.twogood))
        f = list(MinimalFastaParser(self.twogood, strict=False))
        self.assertEqual(len(f), 2)
        a, b = f
        self.assertEqual(a, ('abc', 'caggac'))


if __name__ == '__main__':
    main()
