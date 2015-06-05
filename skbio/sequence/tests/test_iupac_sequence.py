# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio.sequence import IUPACSequence
from skbio.util import classproperty


class ExampleIUPACSequence(IUPACSequence):
    @classproperty
    def degenerate_map(cls):
        return {"X": set("AB"), "Y": set("BC"), "Z": set("AC")}

    @classproperty
    def nondegenerate_chars(cls):
        return set("ABC")


class ExampleMotifsTester(ExampleIUPACSequence):
    @property
    def _motifs(self):
        # These aren't really motifs, just a way to excercise the code paths
        return {
            "name1": lambda x, _, __: str(x),
            "name2": lambda x, _, __: len(x)
        }


class TestIUPACSequence(TestCase):
    def test_instantiation_with_no_implementation(self):
        class IUPACSequenceSubclassNoImplementation(IUPACSequence):
            pass

        with self.assertRaises(TypeError) as cm:
            IUPACSequenceSubclassNoImplementation()

        self.assertIn("abstract class", str(cm.exception))
        self.assertIn("nondegenerate_chars", str(cm.exception))
        self.assertIn("degenerate_map", str(cm.exception))

    def test_init_default_parameters(self):
        seq = ExampleIUPACSequence('.-ABCXYZ')

        npt.assert_equal(seq.values, np.array('.-ABCXYZ', dtype='c'))
        self.assertEqual(seq.id, "")
        self.assertEqual(seq.description, "")
        self.assertIsNone(seq.quality)

    def test_init_nondefault_parameters(self):
        seq = ExampleIUPACSequence('.-ABCXYZ', id='foo', description='bar baz',
                                   quality=range(8))

        npt.assert_equal(seq.values, np.array('.-ABCXYZ', dtype='c'))
        self.assertEqual(seq.id, 'foo')
        self.assertEqual(seq.description, 'bar baz')
        npt.assert_equal(seq.quality, np.array(range(8), dtype='int'))

    def test_init_valid_empty_sequence(self):
        # just make sure we can instantiate an empty sequence regardless of
        # `validate` and `case_insensitive` parameters. more extensive tests
        # are performed in Sequence base class unit tests
        for validate in (True, False):
            for case_insensitive in (True, False):
                seq = ExampleIUPACSequence('', validate=validate,
                                           case_insensitive=case_insensitive)
                self.assertEqual(seq, ExampleIUPACSequence(''))

    def test_init_valid_single_character_sequence(self):
        for validate in (True, False):
            for case_insensitive in (True, False):
                seq = ExampleIUPACSequence('C', validate=validate,
                                           case_insensitive=case_insensitive)
                self.assertEqual(seq, ExampleIUPACSequence('C'))

    def test_init_valid_multiple_character_sequence(self):
        for validate in (True, False):
            for case_insensitive in (True, False):
                seq = ExampleIUPACSequence('BAACB.XYY-AZ', validate=validate,
                                           case_insensitive=case_insensitive)
                self.assertEqual(seq, ExampleIUPACSequence('BAACB.XYY-AZ'))

    def test_init_validate_parameter_single_character(self):
        seq = 'w'

        with self.assertRaisesRegexp(ValueError, "character.*'w'"):
            ExampleIUPACSequence(seq)

        # test that we can instantiate an invalid sequence. we don't guarantee
        # anything working beyond instantiation
        ExampleIUPACSequence(seq, validate=False)

    def test_init_validate_parameter_multiple_characters(self):
        # mix of valid and invalid characters with repeats and lowercased
        # alphabet characters
        seq = 'CBCBBbawCbbwBXYZ-.x'

        with self.assertRaisesRegexp(ValueError, "\['a', 'b', 'w', 'x'\]"):
            ExampleIUPACSequence(seq)

        ExampleIUPACSequence(seq, validate=False)

    def test_init_case_insensitive_lowercase(self):
        s = 'cbcbbbazcbbzbxyz-.x'

        with self.assertRaisesRegexp(ValueError,
                                     "\['a', 'b', 'c', 'x', 'y', 'z'\]"):
            ExampleIUPACSequence(s)

        seq = ExampleIUPACSequence(s, case_insensitive=True)
        self.assertEqual(seq, ExampleIUPACSequence('CBCBBBAZCBBZBXYZ-.X'))

    def test_init_case_insensitive_mixed_case(self):
        s = 'CBCBBbazCbbzBXYZ-.x'

        with self.assertRaisesRegexp(ValueError, "\['a', 'b', 'x', 'z'\]"):
            ExampleIUPACSequence(s)

        seq = ExampleIUPACSequence(s, case_insensitive=True)
        self.assertEqual(seq, ExampleIUPACSequence('CBCBBBAZCBBZBXYZ-.X'))

    def test_init_case_insensitive_no_validation(self):
        s = 'car'

        with self.assertRaisesRegexp(ValueError, "\['a', 'c', 'r'\]"):
            ExampleIUPACSequence(s)

        with self.assertRaisesRegexp(ValueError, "character.*'R'"):
            ExampleIUPACSequence(s, case_insensitive=True)

        ExampleIUPACSequence(s, case_insensitive=True, validate=False)

    def test_init_case_insensitive_byte_ownership(self):
        bytes = np.array([97, 98, 97], dtype=np.uint8)

        with self.assertRaisesRegexp(ValueError, "\['a', 'b'\]"):
            ExampleIUPACSequence(bytes)

        seq = ExampleIUPACSequence(bytes, case_insensitive=True)
        self.assertEqual(seq, ExampleIUPACSequence('ABA'))

        # should not share the same memory
        self.assertIsNot(seq._bytes, bytes)

        # we should have copied `bytes` before modifying in place to convert to
        # upper. make sure `bytes` hasn't been mutated
        npt.assert_equal(bytes, np.array([97, 98, 97], dtype=np.uint8))

    def test_degenerate_chars(self):
        expected = set("XYZ")
        self.assertIs(type(ExampleIUPACSequence.degenerate_chars), set)
        self.assertEqual(ExampleIUPACSequence.degenerate_chars, expected)

        ExampleIUPACSequence.degenerate_chars.add("W")
        self.assertEqual(ExampleIUPACSequence.degenerate_chars, expected)

        self.assertEqual(ExampleIUPACSequence('').degenerate_chars, expected)

        with self.assertRaises(AttributeError):
            ExampleIUPACSequence('').degenerate_chars = set("BAR")

    def test_nondegenerate_chars(self):
        expected = set("ABC")
        self.assertEqual(ExampleIUPACSequence.nondegenerate_chars, expected)

        ExampleIUPACSequence.degenerate_chars.add("D")
        self.assertEqual(ExampleIUPACSequence.nondegenerate_chars, expected)

        self.assertEqual(ExampleIUPACSequence('').nondegenerate_chars,
                         expected)

        with self.assertRaises(AttributeError):
            ExampleIUPACSequence('').nondegenerate_chars = set("BAR")

    def test_gap_chars(self):
        expected = set(".-")
        self.assertIs(type(ExampleIUPACSequence.gap_chars), set)
        self.assertEqual(ExampleIUPACSequence.gap_chars, expected)

        ExampleIUPACSequence.gap_chars.add("_")
        self.assertEqual(ExampleIUPACSequence.gap_chars, expected)

        self.assertEqual(ExampleIUPACSequence('').gap_chars, expected)

        with self.assertRaises(AttributeError):
            ExampleIUPACSequence('').gap_chars = set("_ =")

    def test_alphabet(self):
        expected = set("ABC.-XYZ")
        self.assertIs(type(ExampleIUPACSequence.alphabet), set)
        self.assertEqual(ExampleIUPACSequence.alphabet, expected)

        ExampleIUPACSequence.alphabet.add("DEF")
        self.assertEqual(ExampleIUPACSequence.alphabet, expected)

        self.assertEqual(ExampleIUPACSequence('').alphabet, expected)

        with self.assertRaises(AttributeError):
            ExampleIUPACSequence('').alphabet = set("ABCDEFG.-WXYZ")

    def test_degenerate_map(self):
        expected = {"X": set("AB"), "Y": set("BC"), "Z": set("AC")}
        self.assertEqual(ExampleIUPACSequence.degenerate_map, expected)

        ExampleIUPACSequence.degenerate_map['W'] = set("ABC")
        ExampleIUPACSequence.degenerate_map['X'] = set("CA")
        self.assertEqual(ExampleIUPACSequence.degenerate_map, expected)

        self.assertEqual(ExampleIUPACSequence('').degenerate_map, expected)

        with self.assertRaises(AttributeError):
            ExampleIUPACSequence('').degenerate_map = {'W': "ABC"}

    def test_gaps(self):
        self.assertIs(type(ExampleIUPACSequence("").gaps()), np.ndarray)
        self.assertIs(ExampleIUPACSequence("").gaps().dtype, np.dtype('bool'))
        npt.assert_equal(ExampleIUPACSequence("ABCXBZYABC").gaps(),
                         np.zeros(10).astype(bool))

        npt.assert_equal(ExampleIUPACSequence(".-.-.").gaps(),
                         np.ones(5).astype(bool))

        npt.assert_equal(ExampleIUPACSequence("A.B-C.X-Y.").gaps(),
                         np.array([0, 1] * 5, dtype=bool))

        npt.assert_equal(ExampleIUPACSequence("AB.AC.XY-").gaps(),
                         np.array([0, 0, 1] * 3, dtype=bool))

        npt.assert_equal(ExampleIUPACSequence("A.BC.-").gaps(),
                         np.array([0, 1, 0, 0, 1, 1], dtype=bool))

    def test_has_gaps(self):
        self.assertIs(type(ExampleIUPACSequence("").has_gaps()), bool)
        self.assertIs(type(ExampleIUPACSequence("-").has_gaps()), bool)

        self.assertFalse(ExampleIUPACSequence("").has_gaps())
        self.assertFalse(ExampleIUPACSequence("ABCXYZ").has_gaps())

        self.assertTrue(ExampleIUPACSequence("-").has_gaps())
        self.assertTrue(ExampleIUPACSequence("ABCXYZ-").has_gaps())

    def test_degenerates(self):
        self.assertIs(type(ExampleIUPACSequence("").degenerates()), np.ndarray)
        self.assertIs(ExampleIUPACSequence("").degenerates().dtype,
                      np.dtype('bool'))

        npt.assert_equal(ExampleIUPACSequence("ABCBC-.AB.").degenerates(),
                         np.zeros(10).astype(bool))

        npt.assert_equal(ExampleIUPACSequence("ZYZYZ").degenerates(),
                         np.ones(5).astype(bool))

        npt.assert_equal(ExampleIUPACSequence("AX.Y-ZBXCZ").degenerates(),
                         np.array([0, 1] * 5, dtype=bool))

        npt.assert_equal(ExampleIUPACSequence("ABXACY.-Z").degenerates(),
                         np.array([0, 0, 1] * 3, dtype=bool))

        npt.assert_equal(ExampleIUPACSequence("AZBCXY").degenerates(),
                         np.array([0, 1, 0, 0, 1, 1], dtype=bool))

    def test_has_degenerates(self):
        self.assertIs(type(ExampleIUPACSequence("").has_degenerates()), bool)
        self.assertIs(type(ExampleIUPACSequence("X").has_degenerates()), bool)

        self.assertFalse(ExampleIUPACSequence("").has_degenerates())
        self.assertFalse(ExampleIUPACSequence("A-.BC").has_degenerates())

        self.assertTrue(ExampleIUPACSequence("Z").has_degenerates())
        self.assertTrue(ExampleIUPACSequence("ABC.XYZ-").has_degenerates())

    def test_nondegenerates(self):
        self.assertIs(type(ExampleIUPACSequence("").nondegenerates()),
                      np.ndarray)
        self.assertIs(ExampleIUPACSequence("").nondegenerates().dtype,
                      np.dtype('bool'))

        npt.assert_equal(ExampleIUPACSequence("XYZYZ-.XY.").nondegenerates(),
                         np.zeros(10).astype(bool))

        npt.assert_equal(ExampleIUPACSequence("ABABA").nondegenerates(),
                         np.ones(5).astype(bool))

        npt.assert_equal(ExampleIUPACSequence("XA.B-AZCXA").nondegenerates(),
                         np.array([0, 1] * 5, dtype=bool))

        npt.assert_equal(ExampleIUPACSequence("XXAZZB.-C").nondegenerates(),
                         np.array([0, 0, 1] * 3, dtype=bool))

        npt.assert_equal(ExampleIUPACSequence("YB.-AC").nondegenerates(),
                         np.array([0, 1, 0, 0, 1, 1], dtype=bool))

    def test_has_nondegenerates(self):
        self.assertIs(type(ExampleIUPACSequence("").has_nondegenerates()),
                      bool)
        self.assertIs(type(ExampleIUPACSequence("A").has_nondegenerates()),
                      bool)

        self.assertFalse(ExampleIUPACSequence("").has_nondegenerates())
        self.assertFalse(ExampleIUPACSequence("X-.YZ").has_nondegenerates())

        self.assertTrue(ExampleIUPACSequence("C").has_nondegenerates())
        self.assertTrue(ExampleIUPACSequence(".XYZ-ABC").has_nondegenerates())

    def test_degap(self):
        kw = {
            'id': 'some_id',
            'description': 'some description',
        }

        self.assertEquals(ExampleIUPACSequence("", quality=[], **kw).degap(),
                          ExampleIUPACSequence("", quality=[], **kw))

        self.assertEquals(ExampleIUPACSequence("ABCXYZ", quality=np.arange(6),
                                               **kw).degap(),
                          ExampleIUPACSequence("ABCXYZ", quality=np.arange(6),
                                               **kw))

        self.assertEquals(ExampleIUPACSequence("ABC-XYZ", quality=np.arange(7),
                                               **kw).degap(),
                          ExampleIUPACSequence("ABCXYZ",
                                               quality=[0, 1, 2, 4, 5, 6],
                                               **kw))

        self.assertEquals(ExampleIUPACSequence(".-ABC-XYZ.",
                                               quality=np.arange(10), **kw
                                               ).degap(),
                          ExampleIUPACSequence("ABCXYZ",
                                               quality=[2, 3, 4, 6, 7, 8],
                                               **kw))

        self.assertEquals(ExampleIUPACSequence("---.-.-.-.-.",
                                               quality=np.arange(12), **kw
                                               ).degap(),
                          ExampleIUPACSequence("", quality=[], **kw))

    def test_expand_degenerates_no_degens(self):
        seq = ExampleIUPACSequence("ABCABCABC")
        self.assertEqual(list(seq.expand_degenerates()), [seq])

    def test_expand_degenerates_all_degens(self):
        exp = [ExampleIUPACSequence('ABA'), ExampleIUPACSequence('ABC'),
               ExampleIUPACSequence('ACA'), ExampleIUPACSequence('ACC'),
               ExampleIUPACSequence('BBA'), ExampleIUPACSequence('BBC'),
               ExampleIUPACSequence('BCA'), ExampleIUPACSequence('BCC')]
        # Sort based on sequence string, as order is not guaranteed.
        obs = sorted(ExampleIUPACSequence('XYZ').expand_degenerates(), key=str)
        self.assertEqual(obs, exp)

    def test_expand_degenerates_with_metadata(self):
        kw = {
            "quality": np.arange(3),
            "id": "some_id",
            "description": "some description"
        }
        exp = [ExampleIUPACSequence('ABA', **kw),
               ExampleIUPACSequence('ABC', **kw),
               ExampleIUPACSequence('BBA', **kw),
               ExampleIUPACSequence('BBC', **kw)]
        obs = sorted(ExampleIUPACSequence('XBZ', **kw).expand_degenerates(),
                     key=str)
        self.assertEqual(obs, exp)

    def test_find_motifs_no_motif(self):
        seq = ExampleMotifsTester("ABCABCABC")
        with self.assertRaises(ValueError) as cm:
            seq.find_motifs("doesn't-exist")
        self.assertIn("doesn't-exist", str(cm.exception))

        seq = ExampleIUPACSequence("ABCABCABC")
        with self.assertRaises(ValueError) as cm:
            seq.find_motifs("doesn't-exist")
        self.assertIn("doesn't-exist", str(cm.exception))

    def test_find_motifs(self):
        seq = ExampleMotifsTester("ABC")
        self.assertEqual(seq.find_motifs("name1"), "ABC")
        self.assertEqual(seq.find_motifs("name2"), 3)

if __name__ == "__main__":
    main()
