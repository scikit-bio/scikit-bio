# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from six import StringIO

from unittest import TestCase, main

from skbio.io.clustal import (_clustal_to_alignment, _alignment_to_clustal,
                              _clustal_sniffer)
from skbio.io.clustal import (_is_clustal_seq_line, _delete_trailing_number,
                              _check_length, _label_line_parser)

from skbio.io import ClustalFormatError


class ClustalHelperTests(TestCase):
    def test_label_line_parser(self):
        self.assertEqual(_label_line_parser(StringIO('abc\tucag')),
                         ({"abc": ["ucag"]}, ['abc']))

        with self.assertRaises(ClustalFormatError):
            _label_line_parser(StringIO('abctucag'))

    def test_is_clustal_seq_line(self):
        ic = _is_clustal_seq_line
        self.assertTrue(ic('abc'))
        self.assertTrue(ic('abc  def'))
        self.assertFalse(ic('CLUSTAL'))
        self.assertFalse(ic('CLUSTAL W fsdhicjkjsdk'))
        self.assertFalse(ic('  *   *'))
        self.assertFalse(ic(' abc def'))
        self.assertFalse(ic('MUSCLE (3.41) multiple sequence alignment'))

    def test_delete_trailing_number(self):
        dtn = _delete_trailing_number
        self.assertEqual(dtn('abc'), 'abc')
        self.assertEqual(dtn('a b c'), 'a b c')
        self.assertEqual(dtn('a \t  b  \t  c'), 'a \t  b  \t  c')
        self.assertEqual(dtn('a b 3'), 'a b')
        self.assertEqual(dtn('a b c \t 345'), 'a b c')

    def test_check_lengh(self):
        self.assertEqual(False,
                         _check_length({'abc': ['adjfkadfjaksdlfadskfda'],
                                        'xyz': ['adjfkadfjaksdlfadsk']},
                                       ['abc', 'xyz'])),
        self.assertEqual(True,
                         _check_length({'abc': ['adjfkadfjaksdlfadskfda'],
                                        'xyz': ['adjfkadfjaksdlfadsksdf']},
                                       ['abc', 'xyz']))
        self.assertEqual(True,
                         _check_length({'abc': ['adjfkadfjaksdlfadskfda',
                                                'adjfkadfjaksdlfadskfda'],
                                        'xyz': ['adjfkadfjaksdlfadsksdf',
                                                'adjfkadfjaksdlfadsksdf']},
                                       ['abc', 'xyz']))
        self.assertEqual(False,
                         _check_length({'abc': ['adjfkadfjaksdlfadskfd',
                                                'adjfkadfjaksdlfadskfda'],
                                        'xyz': ['adjfkadfjaksdlfadsksdf',
                                                'adjfkadfjaksdlfadsksdf']},
                                       ['abc', 'xyz']))
        self.assertEqual(False,
                         _check_length({'abc': ['adjfkadfjaksdlfadskfda',
                                                'adjfkadfjaksdlfadskfda'],
                                        'xyz': ['adjfkadfjaksdlfadsksdf',
                                                'adjfkadfjaksdlfadsksd']},
                                       ['abc', 'xyz']))


class ClustalIOTests(TestCase):

    def setUp(self):
        self.valid_clustal_out = [
            StringIO('abc\tucag'),
            StringIO('abc\tuuu\ndef\tccc\n\n    ***\n\ndef ggg\nabc\taaa\n'),
            StringIO('\n'.join(['abc uca', 'def ggg ccc'])),
            StringIO('\n'.join(['abc uca ggg', 'def ggg ccc'])),
            StringIO("""CLUSTAL


abc             GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA
def             ------------------------------------------------------------
xyz             ------------------------------------------------------------


"""),
            StringIO("""CLUSTAL


abc             GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA
def             ------------------------------------------------------------
xyz             ------------------------------------------------------------


abc             GUCGAUACGUACGUCAGUCAGUACGUCAGCAUGCAUACGUACGUCGUACGUACGU-CGAC
def             -----------------------------------------CGCGAUGCAUGCAU-CGAU
xyz             -------------------------------------CAUGCAUCGUACGUACGCAUGAC
"""),
            StringIO("""CLUSTAL W (1.82) multiple sequence alignment


abc             GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA
def             ------------------------------------------------------------
xyz             ------------------------------------------------------------


abc             GUCGAUACGUACGUCAGUCAGUACGUCAGCAUGCAUACGUACGUCGUACGUACGU-CGAC
def             -----------------------------------------CGCGAUGCAUGCAU-CGAU
xyz             -------------------------------------CAUGCAUCGUACGUACGCAUGAC


abc             UGACUAGUCAGCUAGCAUCGAUCAGU
def             CGAUCAGUCAGUCGAU----------
xyz             UGCUGCAUCA----------------"""),
            StringIO("""CLUSTAL W (1.74) multiple sequence alignment


abc             GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA 60
def             ------------------------------------------------------------
xyz             ------------------------------------------------------------


abc             GUCGAUACGUACGUCAGUCAGUACGUCAGCAUGCAUACGUACGUCGUACGUACGU-CGAC 11
def             -----------------------------------------CGCGAUGCAUGCAU-CGAU 18
xyz             -------------------------------------CAUGCAUCGUACGUACGCAUGAC 23
                                                         :    * * * *    **

abc             UGACUAGUCAGCUAGCAUCGAUCAGU 145
def             CGAUCAGUCAGUCGAU---------- 34
xyz             UGCUGCAUCA---------------- 33
                *     ***""")
            ]
        self.invalid_clustal_out = [StringIO('\n'.join(['dshfjsdfhdfsj',
                                                        'hfsdjksdfhjsdf'])),
                                    StringIO('\n'.join(['hfsdjksdfhjsdf'])),
                                    StringIO('\n'.join(['dshfjsdfhdfsj',
                                                        'dshfjsdfhdfsj',
                                                        'hfsdjksdfhjsdf'])),
                                    StringIO('\n'.join(['dshfjsdfhdfsj',
                                                        '\t',
                                                        'hfsdjksdfhjsdf'])),
                                    StringIO('\n'.join(['dshfj\tdfhdfsj',
                                                        'hfsdjksdfhjsdf'])),
                                    StringIO('\n'.join(['dshfjsdfhdfsj',
                                                        'hfsdjk\tdfhjsdf'])),
                                    StringIO("""CLUSTAL W (1.74) multiple sequence alignment


adj GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA
------------------------------------------------------------
adk -----GGGGGGG------------------------------------------------
"""),
                                    StringIO("""CLUSTAL W (1.74) multiple sequence alignment


adj GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA
adk -----GGGGGGG------------------------------------------------


adj GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA
adk -----GGGGGGG---------------------------------------------
"""),
                                    StringIO("""CLUSTAL W (1.74) multiple sequence alignment


adj GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA
adk -----GGGGGGG---------------------------------------------


adj GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCA
adk -----GGGGGGG---------------------------------------------
"""),

                                    StringIO("""CLUSTAL W (1.74) multiple sequence alignment


adj GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA
------------------------------------------------------------
adk -----GGGGGGG------------------------------------------------
"""),

                                    StringIO("""CLUSTAL W (1.74) multiple sequence alignment


GCAUGCAUGCAUGAUCGUACGUCAGCAUGCUAGACUGCAUACGUACGUACGCAUGCAUCA
------------------------------------------------------------
------------------------------------------------------------


GUCGAUACGUACGUCAGUCAGUACGUCAGCAUGCAUACGUACGUCGUACGUACGU-CGAC
-----------------------------------------CGCGAUGCAUGCAU-CGAU
------------------------------------------------------------
                                         :    * * * *    **

UGACUAGUCAGCUAGCAUCGAUCAGU 145
CGAUCAGUCAGUCGAU---------- 34
UGCUGCAUCA---------------- 33
*     ***""")]

    def test_alignment_to_clustal_with_empty_input(self):
        result = _clustal_to_alignment(StringIO())
        self.assertEqual(dict(result), {})

    def test_alignment_to_clustal_with_bad_input(self):
        BAD = StringIO('\n'.join(['dshfjsdfhdfsj', 'hfsdjksdfhjsdf']))
        result = _clustal_to_alignment(BAD, strict=False)
        self.assertEqual(dict(result), {})
        # should fail unless we turned strict processing off
        with self.assertRaises(ClustalFormatError):
            BAD.seek(0)
            dict(_clustal_to_alignment(BAD))

    def test_valid_alignment_to_clustal_and_clustal_to_alignment(self):
        import os
        for valid_out in self.valid_clustal_out:
            fname = "test.aln"
            testfile = open(fname, 'w')
            result_before = _clustal_to_alignment(valid_out)
            _alignment_to_clustal(result_before, testfile)
            testfile.close()
            testfile = open(fname, 'r')
            result_after = _clustal_to_alignment(testfile)
            self.assertEqual(result_before, result_after)
        os.remove(fname)

    def test_invalid_alignment_to_clustal_and_clustal_to_alignment(self):
        for invalid_out in self.invalid_clustal_out:
            with self.assertRaises(ClustalFormatError):
                dict(_clustal_to_alignment(invalid_out, strict=True))

    def test_clustal_sniffer_valid_files(self):
        for valid_out in self.valid_clustal_out:
            self.assertEqual(_clustal_sniffer(valid_out), (True, {}))

    def test_clustal_sniffer_invalid_files(self):
        for invalid_out in self.invalid_clustal_out:
            self.assertEqual(_clustal_sniffer(invalid_out), (False, {}))
        # sniffer should return False on empty file (which isn't contained
        # in self.invalid_clustal_out since an empty file is a valid output)
        self.assertEqual(_clustal_sniffer(StringIO()), (False, {}))

if __name__ == '__main__':
    main()
