# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
import six

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio.sequence._iupac_sequence import IUPACSequence
from skbio.util._decorator import classproperty


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
    def setUp(self):
        self.lowercase_seq = ExampleIUPACSequence('AAAAaaaa', lowercase='key')

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
        self.assertFalse(seq.has_metadata())
        self.assertFalse(seq.has_positional_metadata())

    def test_init_nondefault_parameters(self):
        seq = ExampleIUPACSequence('.-ABCXYZ',
                                   metadata={'id': 'foo'},
                                   positional_metadata={'quality': range(8)})

        npt.assert_equal(seq.values, np.array('.-ABCXYZ', dtype='c'))
        self.assertTrue(seq.has_metadata())
        self.assertEqual(seq.metadata['id'], 'foo')
        self.assertTrue(seq.has_positional_metadata())
        npt.assert_equal(seq.positional_metadata['quality'], np.array(range(8),
                         dtype='int'))

    def test_init_valid_empty_sequence(self):
        # just make sure we can instantiate an empty sequence regardless of
        # `validate` and `lowercase` parameters. more extensive tests
        # are performed in Sequence base class unit tests
        for validate in (True, False):
            for lowercase in (True, False):
                seq = ExampleIUPACSequence('', validate=validate,
                                           lowercase=lowercase)
                self.assertEqual(seq, ExampleIUPACSequence(''))

    def test_init_valid_single_character_sequence(self):
        for validate in (True, False):
            for lowercase in (True, False):
                seq = ExampleIUPACSequence('C', validate=validate,
                                           lowercase=lowercase)
                self.assertEqual(seq, ExampleIUPACSequence('C'))

    def test_init_valid_multiple_character_sequence(self):
        for validate in (True, False):
            for lowercase in (True, False):
                seq = ExampleIUPACSequence('BAACB.XYY-AZ', validate=validate,
                                           lowercase=lowercase)
                self.assertEqual(seq, ExampleIUPACSequence('BAACB.XYY-AZ'))

    def test_init_validate_parameter_single_character(self):
        seq = 'w'

        with six.assertRaisesRegex(self, ValueError, "character.*'w'"):
            ExampleIUPACSequence(seq)

        # test that we can instantiate an invalid sequence. we don't guarantee
        # anything working beyond instantiation
        ExampleIUPACSequence(seq, validate=False)

    def test_init_validate_parameter_multiple_characters(self):
        # mix of valid and invalid characters with repeats and lowercased
        # alphabet characters
        seq = 'CBCBBbawCbbwBXYZ-.x'

        with six.assertRaisesRegex(self, ValueError, "\['a', 'b', 'w', 'x'\]"):
            ExampleIUPACSequence(seq)

        ExampleIUPACSequence(seq, validate=False)

    def test_init_lowercase_all_lowercase(self):
        s = 'cbcbbbazcbbzbxyz-.x'

        with six.assertRaisesRegex(self, ValueError,
                                   "\['a', 'b', 'c', 'x', 'y', 'z'\]"):
            ExampleIUPACSequence(s)

        seq = ExampleIUPACSequence(s, lowercase=True)
        self.assertEqual(seq, ExampleIUPACSequence('CBCBBBAZCBBZBXYZ-.X'))

    def test_init_lowercase_mixed_case(self):
        s = 'CBCBBbazCbbzBXYZ-.x'

        with six.assertRaisesRegex(self, ValueError, "\['a', 'b', 'x', 'z'\]"):
            ExampleIUPACSequence(s)

        seq = ExampleIUPACSequence(s, lowercase=True)
        self.assertEqual(seq, ExampleIUPACSequence('CBCBBBAZCBBZBXYZ-.X'))

    def test_init_lowercase_no_validation(self):
        s = 'car'

        with six.assertRaisesRegex(self, ValueError, "\['a', 'c', 'r'\]"):
            ExampleIUPACSequence(s)

        with six.assertRaisesRegex(self, ValueError, "character.*'R'"):
            ExampleIUPACSequence(s, lowercase=True)

        ExampleIUPACSequence(s, lowercase=True, validate=False)

    def test_init_lowercase_byte_ownership(self):
        bytes = np.array([97, 98, 97], dtype=np.uint8)

        with six.assertRaisesRegex(self, ValueError, "\['a', 'b'\]"):
            ExampleIUPACSequence(bytes)

        seq = ExampleIUPACSequence(bytes, lowercase=True)
        self.assertEqual(seq, ExampleIUPACSequence('ABA'))

        # should not share the same memory
        self.assertIsNot(seq._bytes, bytes)

        # we should have copied `bytes` before modifying in place to convert to
        # upper. make sure `bytes` hasn't been mutated
        npt.assert_equal(bytes, np.array([97, 98, 97], dtype=np.uint8))

    def test_init_lowercase_invalid_keys(self):
        for invalid_key in ((), [], 2):
            invalid_type = type(invalid_key)
            with six.assertRaisesRegex(self, TypeError,
                                       "lowercase keyword argument expected "
                                       "a bool or string, but got %s" %
                                       invalid_type):
                ExampleIUPACSequence('ACGTacgt', lowercase=invalid_key)

    def test_lowercase_mungeable_key(self):
        # NOTE: This test relies on Sequence._munge_to_index_array working
        # properly. If the internal implementation of the lowercase method
        # changes to no longer use _munge_to_index_array, this test may need
        # to be updated to cover cases currently covered by
        # _munge_to_index_array
        self.assertEqual('AAAAaaaa', self.lowercase_seq.lowercase('key'))

    def test_lowercase_array_key(self):
        # NOTE: This test relies on Sequence._munge_to_index_array working
        # properly. If the internal implementation of the lowercase method
        # changes to no longer use _munge_to_index_array, this test may need
        # to be updated to cover cases currently covered by
        # _munge_to_index_array
        self.assertEqual('aaAAaaaa',
                         self.lowercase_seq.lowercase(
                             np.array([True, True, False, False, True, True,
                                       True, True])))
        self.assertEqual('AaAAaAAA',
                         self.lowercase_seq.lowercase([1, 4]))

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
            'metadata': {
                'id': 'some_id',
                'description': 'some description',
            },
        }

        self.assertEqual(
            ExampleIUPACSequence("", positional_metadata={'qual': []},
                                 **kw).degap(),
            ExampleIUPACSequence("", positional_metadata={'qual': []},
                                 **kw))

        self.assertEqual(
            ExampleIUPACSequence(
                "ABCXYZ",
                positional_metadata={'qual': np.arange(6)},
                **kw).degap(),
            ExampleIUPACSequence(
                "ABCXYZ",
                positional_metadata={'qual': np.arange(6)},
                **kw))

        self.assertEqual(
            ExampleIUPACSequence(
                "ABC-XYZ",
                positional_metadata={'qual': np.arange(7)},
                **kw).degap(),
            ExampleIUPACSequence(
                "ABCXYZ",
                positional_metadata={'qual': [0, 1, 2, 4, 5, 6]},
                **kw))

        self.assertEqual(
            ExampleIUPACSequence(
                ".-ABC-XYZ.",
                positional_metadata={'qual': np.arange(10)},
                **kw).degap(),
            ExampleIUPACSequence(
                "ABCXYZ",
                positional_metadata={'qual': [2, 3, 4, 6, 7, 8]},
                **kw))

        self.assertEqual(
            ExampleIUPACSequence(
                "---.-.-.-.-.",
                positional_metadata={'quality': np.arange(12)},
                **kw).degap(),
            ExampleIUPACSequence(
                "",
                positional_metadata={'quality': np.array([], dtype=np.int64)},
                **kw))

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
            "metadata": {
                "id": "some_id",
                "description": "some description"
            },
            "positional_metadata": {
                "quality": np.arange(3),
            },
        }
        exp = [ExampleIUPACSequence('ABA', **kw),
               ExampleIUPACSequence('ABC', **kw),
               ExampleIUPACSequence('BBA', **kw),
               ExampleIUPACSequence('BBC', **kw)]
        obs = sorted(ExampleIUPACSequence('XBZ', **kw).expand_degenerates(),
                     key=str)
        self.assertEqual(obs, exp)

    def test_to_regex_no_degens(self):
        seq = ExampleIUPACSequence('ABC')
        regex = seq.to_regex()
        self.assertEquals(regex.pattern, str(seq))

    def test_to_regex_with_degens(self):
        seq = ExampleIUPACSequence('AYZ')
        regex = seq.to_regex()
        self.assertFalse(any(regex.match(s) is None
                             for s in 'ABA ABC ACA ACC'.split()))
        self.assertTrue(all(regex.match(s) is None
                            for s in 'CBA BBA ABB AAA'.split()))

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

    def test_repr(self):
        # basic sanity checks for custom repr stats. more extensive testing is
        # performed on Sequence.__repr__

        # minimal
        obs = repr(ExampleIUPACSequence(''))
        self.assertEqual(obs.count('\n'), 7)
        self.assertTrue(obs.startswith('ExampleIUPACSequence'))
        self.assertIn('length: 0', obs)
        self.assertIn('has gaps: False', obs)
        self.assertIn('has degenerates: False', obs)
        self.assertIn('has non-degenerates: False', obs)
        self.assertTrue(obs.endswith('-'))

        # no metadata, mix of gaps, degenerates, and non-degenerates
        obs = repr(ExampleIUPACSequence('AY-B'))
        self.assertEqual(obs.count('\n'), 8)
        self.assertTrue(obs.startswith('ExampleIUPACSequence'))
        self.assertIn('length: 4', obs)
        self.assertIn('has gaps: True', obs)
        self.assertIn('has degenerates: True', obs)
        self.assertIn('has non-degenerates: True', obs)
        self.assertTrue(obs.endswith('0 AY-B'))

        # metadata and positional metadata of mixed types
        obs = repr(
            ExampleIUPACSequence(
                'ABCA',
                metadata={'foo': 42, u'bar': 33.33, None: True, False: {},
                          (1, 2): 3, 'acb' * 100: "'"},
                positional_metadata={'foo': range(4),
                                     42: ['a', 'b', [], 'c']}))
        self.assertEqual(obs.count('\n'), 18)
        self.assertTrue(obs.startswith('ExampleIUPACSequence'))
        self.assertIn('None: True', obs)
        self.assertIn('\'foo\': 42', obs)
        self.assertIn('42: <dtype: object>', obs)
        self.assertIn('\'foo\': <dtype: int64>', obs)
        self.assertIn('length: 4', obs)
        self.assertIn('has gaps: False', obs)
        self.assertIn('has degenerates: False', obs)
        self.assertIn('has non-degenerates: True', obs)
        self.assertTrue(obs.endswith('0 ABCA'))

        # sequence spanning > 5 lines
        obs = repr(ExampleIUPACSequence('A' * 301))
        self.assertEqual(obs.count('\n'), 12)
        self.assertTrue(obs.startswith('ExampleIUPACSequence'))
        self.assertIn('length: 301', obs)
        self.assertIn('has gaps: False', obs)
        self.assertIn('has degenerates: False', obs)
        self.assertIn('has non-degenerates: True', obs)
        self.assertIn('...', obs)
        self.assertTrue(obs.endswith('300 A'))


if __name__ == "__main__":
    main()
