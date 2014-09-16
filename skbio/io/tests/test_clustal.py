#!/usr/bin/env python
# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import absolute_import, division, print_function
from future.utils.six import StringIO

from unittest import TestCase, main

from skbio.io.clustal import _clustal_to_generator, _generator_to_clustal
from skbio.io.clustal import (_is_clustal_seq_line, last_space,
                              _delete_trailing_number)
from skbio.io import RecordError


class ClustalHelperTests(TestCase):

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


class ClustalIOTests(TestCase):
    """Tests of the _clustal_to_generator function
    and _generator_to_clustal function"""

    def setUp(self):
        self.valid_clustal_out = [
            StringIO('abc\tucag'),
            StringIO('abc\tuuu\ndef\tccc\n\n    ***\n\ndef ggg\nabc\taaa\n'),
            StringIO('\n'.join(['abc uca', 'def ggg ccc'])),
            StringIO("""CLUSTAL W (1.82) multiple sequence alignment


abc             GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA
def             ------------------------------------------------------------
xyz             ------------------------------------------------------------


abc             GUCGAUACGUACGUCAGUCAGUACGUCAGCAUGCAUACGUACGUCGUACGUACGU-CGAC
def             -----------------------------------------CGCGAUGCAUGCAU-CGAU
xyz             -------------------------------------CAUGCAUCGUACGUACGCAUGAC


abc             UGACUAGUCAGCUAGCAUCGAUCAGU
def             CGAUCAGUCAGUCGAU----------
xyz             UGCUGCAUCA----------------""")

            ]
        self.invalid = [StringIO('\n'.join(['dshfjsdfhdfsj',
                                            'hfsdjksdfhjsdf']))]

    def test_generator_to_clustal_with_empty_input(self):
        """Should return empty dict and list on null input"""
        result = _clustal_to_generator(StringIO())
        self.assertEqual(dict(result), {})

    def test_real(self):
        """Should handle real Clustal output"""

        REAL = StringIO("""CLUSTAL W (1.82) multiple sequence alignment


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
                *     ***""")

        data = _clustal_to_generator(REAL)
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

    def test_valid(self):
        import os
        for valid_out in self.valid_clustal_out:
            fname = "test.aln"
            testfile = open(fname, 'w')
            result_before = _clustal_to_generator(valid_out)
            records = list(result_before)
            _generator_to_clustal(records, testfile)
            testfile.close()
            testfile = open(fname, 'r')
            result_after = _clustal_to_generator(testfile)
            self.assertEquals(set(records), set(result_after))
        os.remove(fname)

if __name__ == '__main__':
    main()
