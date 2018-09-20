# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import string
from io import StringIO
from unittest import TestCase, main

from skbio import TabularMSA
from skbio.sequence import GrammaredSequence
from skbio.util import classproperty
from skbio.util._decorator import overrides
from skbio.io.format.clustal import (
    _clustal_to_tabular_msa, _tabular_msa_to_clustal, _clustal_sniffer,
    _is_clustal_seq_line, _delete_trailing_number, _check_length,
    _label_line_parser)

from skbio.io import ClustalFormatError


class CustomSequence(GrammaredSequence):
    @classproperty
    @overrides(GrammaredSequence)
    def gap_chars(cls):
        return set('-.')

    @classproperty
    @overrides(GrammaredSequence)
    def default_gap_char(cls):
        return '-'

    @classproperty
    @overrides(GrammaredSequence)
    def definite_chars(cls):
        return set(string.ascii_letters)

    @classproperty
    @overrides(GrammaredSequence)
    def degenerate_map(cls):
        return {}


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
            StringIO('CLUSTAL\n\nabc\tucag'),
            StringIO('CLUSTAL\n\nabc\tuuu\ndef\tccc\n\n    ***\n\ndef ggg\nab'
                     'c\taaa\n'),
            StringIO('\n'.join(['CLUSTAL\n', 'abc uca', 'def ggg ccc'])),
            StringIO('\n'.join(['CLUSTAL\n', 'abc uca ggg', 'def ggg ccc'])),
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

    def test_tabular_msa_to_clustal_with_empty_input(self):
        result = _clustal_to_tabular_msa(StringIO(),
                                         constructor=CustomSequence)
        self.assertEqual(dict(result), {})

    def test_tabular_msa_to_clustal_with_bad_input(self):
        BAD = StringIO('\n'.join(['dshfjsdfhdfsj', 'hfsdjksdfhjsdf']))

        with self.assertRaises(ClustalFormatError):
            dict(_clustal_to_tabular_msa(BAD, constructor=CustomSequence))

    def test_valid_tabular_msa_to_clustal_and_clustal_to_tabular_msa(self):
        for valid_out in self.valid_clustal_out:
            result_before = _clustal_to_tabular_msa(
                    valid_out, constructor=CustomSequence)
            with StringIO() as fh:
                _tabular_msa_to_clustal(result_before, fh)
                fh.seek(0)
                result_after = _clustal_to_tabular_msa(
                        fh, constructor=CustomSequence)
            self.assertEqual(result_before, result_after)

    def test_invalid_tabular_msa_to_clustal_and_clustal_to_tabular_msa(self):
        for invalid_out in self.invalid_clustal_out:
            with self.assertRaises(ClustalFormatError):
                dict(_clustal_to_tabular_msa(invalid_out,
                                             constructor=CustomSequence))

    def test_clustal_sniffer_valid_files(self):
        for valid_out in self.valid_clustal_out:
            self.assertEqual(_clustal_sniffer(valid_out), (True, {}))

    def test_clustal_sniffer_invalid_files(self):
        for invalid_out in self.invalid_clustal_out:
            self.assertEqual(_clustal_sniffer(invalid_out), (False, {}))
        # sniffer should return False on empty file (which isn't contained
        # in self.invalid_clustal_out since an empty file is a valid output)
        self.assertEqual(_clustal_sniffer(StringIO()), (False, {}))

    def test_no_constructor(self):
        with self.assertRaisesRegex(ValueError, r"`constructor`"):
            _clustal_to_tabular_msa(self.valid_clustal_out[0])

    def test_duplicate_labels(self):
        msa = TabularMSA([CustomSequence('foo'),
                          CustomSequence('bar')], index=['a', 'a'])

        with self.assertRaisesRegex(ClustalFormatError, r"index.*unique"):
            with StringIO() as fh:
                _tabular_msa_to_clustal(msa, fh)

    def test_invalid_lengths(self):
        fh = StringIO(
            "CLUSTAL\n"
            "\n\n"
            "abc             GCAU\n"
            "def             -----\n")

        with self.assertRaisesRegex(ClustalFormatError, r"not aligned"):
            _clustal_to_tabular_msa(fh, constructor=CustomSequence)


if __name__ == '__main__':
    main()
