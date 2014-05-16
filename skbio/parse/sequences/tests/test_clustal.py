#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import division

from skbio.parse.sequences.clustal import (_is_clustal_seq_line, last_space,
                                           _delete_trailing_number,
                                           parse_clustal)
from skbio.core.exception import RecordError

from unittest import TestCase, main


class ClustalTests(TestCase):

    """Tests of top-level functions."""

    def test_is_clustal_seq_line(self):
        """_is_clustal_seq_line should reject blanks and 'CLUSTAL'"""
        ic = _is_clustal_seq_line
        assert ic('abc')
        assert ic('abc  def')
        assert not ic('CLUSTAL')
        assert not ic('CLUSTAL W fsdhicjkjsdk')
        assert not ic('  *   *')
        assert not ic(' abc def')
        assert not ic('MUSCLE (3.41) multiple sequence alignment')

    def test_last_space(self):
        """last_space should split on last whitespace"""
        self.assertEqual(last_space('a\t\t\t  b    c'), ['a b', 'c'])
        self.assertEqual(last_space('xyz'), ['xyz'])
        self.assertEqual(last_space('  a b'), ['a', 'b'])

    def test_delete_trailing_number(self):
        """Should delete the trailing number if present"""
        dtn = _delete_trailing_number
        self.assertEqual(dtn('abc'), 'abc')
        self.assertEqual(dtn('a b c'), 'a b c')
        self.assertEqual(dtn('a \t  b  \t  c'), 'a \t  b  \t  c')
        self.assertEqual(dtn('a b 3'), 'a b')
        self.assertEqual(dtn('a b c \t 345'), 'a b c')


class ClustalParserTests(TestCase):

    """Tests of the parse_clustal function"""

    def test_null(self):
        """Should return empty dict and list on null input"""
        result = parse_clustal([])
        self.assertEqual(dict(result), {})

    def test_minimal(self):
        """Should handle single-line input correctly"""
        result = parse_clustal([MINIMAL])  # expects seq of lines
        self.assertEqual(dict(result), {'abc': 'ucag'})

    def test_two(self):
        """Should handle two-sequence input correctly"""
        result = parse_clustal(TWO)
        self.assertEqual(dict(result), {'abc': 'uuuaaa', 'def': 'cccggg'})

    def test_real(self):
        """Should handle real Clustal output"""
        data = parse_clustal(REAL)
        self.assertEqual(dict(data), {
            'abc':
            'GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA'
            'GUCGAUACGUACGUCAGUCAGUACGUCAGCAUGCAUACGUACGUCGUACGUACGU-CGAC'
            'UGACUAGUCAGCUAGCAUCGAUCAGU',
            'def':
            '------------------------------------------------------------'
            '-----------------------------------------CGCGAUGCAUGCAU-CGAU'
            'CGAUCAGUCAGUCGAU----------',
            'xyz':
            '------------------------------------------------------------'
            '-------------------------------------CAUGCAUCGUACGUACGCAUGAC'
            'UGCUGCAUCA----------------'
        })

    def test_bad(self):
        """Should reject bad data if strict"""
        result = parse_clustal(BAD, strict=False)
        self.assertEqual(dict(result), {})
        # should fail unless we turned strict processing off
        with self.assertRaises(RecordError):
            dict(parse_clustal(BAD))

    def test_space_labels(self):
        """Should tolerate spaces in labels"""
        result = parse_clustal(SPACE_LABELS)
        self.assertEqual(dict(result), {'abc': 'uca', 'def ggg': 'ccc'})


MINIMAL = 'abc\tucag'
TWO = 'abc\tuuu\ndef\tccc\n\n    ***\n\ndef ggg\nabc\taaa\n'.split('\n')

REAL = """CLUSTAL W (1.82) multiple sequence alignment


abc             GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA 60
def             ------------------------------------------------------------
xyz             ------------------------------------------------------------


abc             GUCGAUACGUACGUCAGUCAGUACGUCAGCAUGCAUACGUACGUCGUACGUACGU-CGAC 11
def             -----------------------------------------CGCGAUGCAUGCAU-CGAU 18
xyz             -------------------------------------CAUGCAUCGUACGUACGCAUGAC 23
                                                         *    * * * *    **

abc             UGACUAGUCAGCUAGCAUCGAUCAGU 145
def             CGAUCAGUCAGUCGAU---------- 34
xyz             UGCUGCAUCA---------------- 33
                *     ***""".split('\n')

BAD = ['dshfjsdfhdfsj', 'hfsdjksdfhjsdf']

SPACE_LABELS = ['abc uca', 'def ggg ccc']


if __name__ == '__main__':
    main()
