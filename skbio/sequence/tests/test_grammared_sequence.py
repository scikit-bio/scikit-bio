# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd

from skbio.sequence import GrammaredSequence
from skbio.util import classproperty
from skbio.util import assert_data_frame_almost_equal
from skbio.metadata import IntervalMetadata


class ExampleGrammaredSequence(GrammaredSequence):
    @classproperty
    def degenerate_map(cls):
        return {"X": set("AB"), "Y": set("BC"), "Z": set("AC"), "W": set("ABCQ")}

    @classproperty
    def definite_chars(cls):
        return set("ABCQ")

    @classproperty
    def default_gap_char(cls):
        return '-'

    @classproperty
    def gap_chars(cls):
        return set('-.')

    @classproperty
    def noncanonical_chars(cls):
        return "Q"

    @classproperty
    def wildcard_char(cls):
        return "W"


class ExampleMotifsTester(ExampleGrammaredSequence):
    @property
    def _motifs(self):
        # These aren't really motifs, just a way to excercise the code paths
        return {
            "name1": lambda x, _, __: str(x),
            "name2": lambda x, _, __: len(x)
        }


class TestGrammaredSequence(TestCase):
    def test_default_gap_must_be_in_gap_chars(self):
        with self.assertRaisesRegex(
                TypeError,
                r"default_gap_char must be in gap_chars for class "
                "GrammaredSequenceInvalidDefaultGap"):

            class GrammaredSequenceInvalidDefaultGap(ExampleGrammaredSequence):
                @classproperty
                def default_gap_char(cls):
                    return '*'

    def test_degenerates_must_expand_to_valid_definites(self):
        with self.assertRaisesRegex(
                TypeError,
                r"degenerate_map must expand only to characters included in "
                "definite_chars for class "
                "GrammaredSequenceInvalidDefaultGap"):

            class GrammaredSequenceInvalidDefaultGap(ExampleGrammaredSequence):
                @classproperty
                def degenerate_map(cls):
                    return {"X": set("B")}

                @classproperty
                def definite_chars(cls):
                    return set("A")

    def test_gap_chars_and_degenerates_share(self):
        with self.assertRaisesRegex(
                TypeError,
                r"gap_chars and degenerate_chars must not share any characters"
                " for class GrammaredSequenceGapInDegenerateMap"):

            class GrammaredSequenceGapInDegenerateMap(
                    ExampleGrammaredSequence):
                @classproperty
                def degenerate_map(cls):
                    return {"X": set("AB")}

                @classproperty
                def definite_chars(cls):
                    return set("ABC")

                @classproperty
                def gap_chars(cls):
                    return set(".-X")

    def test_gap_chars_and_definites_share(self):
        with self.assertRaisesRegex(
            TypeError,
            (r"gap_chars and definite_chars must not share any characters "
             "for class GrammaredSequenceGapInDefiniteMap")):

            class GrammaredSequenceGapInDefiniteMap(
                    ExampleGrammaredSequence):
                @classproperty
                def degenerate_map(cls):
                    return {"X": set("AB")}

                @classproperty
                def definite_chars(cls):
                    return set("ABC")

                @classproperty
                def gap_chars(cls):
                    return set(".-A")

    def test_degenerates_and_definites_share(self):
        with self.assertRaisesRegex(
            TypeError,
            (r"degenerate_chars and definite_chars must not share any "
             "characters for class GrammaredSequenceInvalid")):

            class GrammaredSequenceInvalid(ExampleGrammaredSequence):
                @classproperty
                def degenerate_map(cls):
                    return {"X": set("AB")}

                @classproperty
                def definite_chars(cls):
                    return set("ABCX")

    def test_instantiation_with_no_implementation(self):
        class GrammaredSequenceSubclassNoImplementation(GrammaredSequence):
            pass

        with self.assertRaises(TypeError) as cm:
            GrammaredSequenceSubclassNoImplementation()

        self.assertIn("abstract class", str(cm.exception))
        self.assertIn("definite_chars", str(cm.exception))
        self.assertIn("degenerate_map", str(cm.exception))

    def test_init_default_parameters(self):
        seq = ExampleGrammaredSequence('.-ABCXYZ')

        npt.assert_equal(seq.values, np.array('.-ABCXYZ', dtype='c'))
        self.assertEqual(seq.metadata, {})
        assert_data_frame_almost_equal(seq.positional_metadata,
                                       pd.DataFrame(index=range(8)))
        self.assertEqual(seq.interval_metadata,
                         IntervalMetadata(8))

    def test_init_nondefault_parameters(self):
        im = IntervalMetadata(8)
        im.add([(1, 8)], metadata={'gene': 'p53'})
        seq = ExampleGrammaredSequence(
            '.-ABCXYZ',
            metadata={'id': 'foo'},
            positional_metadata={'quality': range(8)},
            interval_metadata=im)

        npt.assert_equal(seq.values, np.array('.-ABCXYZ', dtype='c'))
        self.assertEqual(seq.metadata, {'id': 'foo'})
        assert_data_frame_almost_equal(seq.positional_metadata,
                                       pd.DataFrame({'quality': range(8)}))
        self.assertEqual(seq.interval_metadata, im)

    def test_init_valid_empty_sequence(self):
        # just make sure we can instantiate an empty sequence regardless of
        # `validate` and `lowercase` parameters. more extensive tests
        # are performed in Sequence base class unit tests
        for validate in (True, False):
            for lowercase in (True, False):
                seq = ExampleGrammaredSequence(
                    '', validate=validate, lowercase=lowercase)
                self.assertEqual(seq, ExampleGrammaredSequence(''))

    def test_init_valid_single_character_sequence(self):
        for validate in (True, False):
            for lowercase in (True, False):
                seq = ExampleGrammaredSequence(
                    'C', validate=validate, lowercase=lowercase)
                self.assertEqual(seq, ExampleGrammaredSequence('C'))

    def test_init_valid_multiple_character_sequence(self):
        for validate in (True, False):
            for lowercase in (True, False):
                seq = ExampleGrammaredSequence(
                    'BAACB.XYY-AZ', validate=validate, lowercase=lowercase)
                self.assertEqual(seq, ExampleGrammaredSequence('BAACB.XYY-AZ'))

    def test_init_validate_parameter_single_character(self):
        seq = 'w'

        with self.assertRaisesRegex(ValueError, r"character.*'w'"):
            ExampleGrammaredSequence(seq)

        # test that we can instantiate an invalid sequence. we don't guarantee
        # anything working beyond instantiation
        ExampleGrammaredSequence(seq, validate=False)

    def test_init_validate_parameter_multiple_characters(self):
        # mix of valid and invalid characters with repeats and lowercased
        # alphabet characters
        seq = 'CBCBBbawCbbwBXYZ-.x'

        with self.assertRaisesRegex(ValueError, r"\['a', 'b', 'w', 'x'\]"):
            ExampleGrammaredSequence(seq)

        ExampleGrammaredSequence(seq, validate=False)

    def test_init_lowercase_all_lowercase(self):
        s = 'cbcbbbazcbbzbxyz-.x'

        with self.assertRaisesRegex(ValueError,
                                    r"\['a', 'b', 'c', 'x', 'y', 'z'\]"):
            ExampleGrammaredSequence(s)

        seq = ExampleGrammaredSequence(s, lowercase=True)
        self.assertEqual(seq, ExampleGrammaredSequence('CBCBBBAZCBBZBXYZ-.X'))

    def test_init_lowercase_mixed_case(self):
        s = 'CBCBBbazCbbzBXYZ-.x'

        with self.assertRaisesRegex(ValueError, r"\['a', 'b', 'x', 'z'\]"):
            ExampleGrammaredSequence(s)

        seq = ExampleGrammaredSequence(s, lowercase=True)
        self.assertEqual(seq, ExampleGrammaredSequence('CBCBBBAZCBBZBXYZ-.X'))

    def test_init_lowercase_no_validation(self):
        s = 'car'

        with self.assertRaisesRegex(ValueError, r"\['a', 'c', 'r'\]"):
            ExampleGrammaredSequence(s)

        with self.assertRaisesRegex(ValueError, r"character.*'R'"):
            ExampleGrammaredSequence(s, lowercase=True)

        ExampleGrammaredSequence(s, lowercase=True, validate=False)

    def test_init_lowercase_byte_ownership(self):
        bytes = np.array([97, 98, 97], dtype=np.uint8)

        with self.assertRaisesRegex(ValueError, r"\['a', 'b'\]"):
            ExampleGrammaredSequence(bytes)

        seq = ExampleGrammaredSequence(bytes, lowercase=True)
        self.assertEqual(seq, ExampleGrammaredSequence('ABA'))

        # should not share the same memory
        self.assertIsNot(seq._bytes, bytes)

        # we should have copied `bytes` before modifying in place to convert to
        # upper. make sure `bytes` hasn't been mutated
        npt.assert_equal(bytes, np.array([97, 98, 97], dtype=np.uint8))

    def test_init_lowercase_invalid_keys(self):
        for invalid_key in ((), [], 2):
            invalid_type = type(invalid_key)
            with self.assertRaisesRegex(TypeError,
                                        r"lowercase keyword argument expected "
                                        "a bool or string, but got %s" %
                                        invalid_type):
                ExampleGrammaredSequence('ACGTacgt', lowercase=invalid_key)

    def test_definite_char_codes(self):
        definite_char_codes = set(ExampleGrammaredSequence._definite_char_codes)
        self.assertEqual(definite_char_codes, set([65, 66, 67, 81]))

    def test_gap_codes(self):
        gap_codes = set(ExampleGrammaredSequence._gap_codes)
        self.assertEqual(gap_codes, set([45, 46]))

    def test_noncanonical_codes(self):
        noncanonical_codes = set(ExampleGrammaredSequence._noncanonical_codes)
        self.assertEqual(noncanonical_codes, set([81]))

    def test_degenerate_chars(self):
        expected = set("WXYZ")
        self.assertIs(type(ExampleGrammaredSequence.degenerate_chars), set)
        self.assertEqual(ExampleGrammaredSequence.degenerate_chars, expected)

        ExampleGrammaredSequence.degenerate_chars.add("W")
        self.assertEqual(ExampleGrammaredSequence.degenerate_chars, expected)

        self.assertEqual(ExampleGrammaredSequence('').degenerate_chars,
                         expected)

        with self.assertRaises(AttributeError):
            ExampleGrammaredSequence('').degenerate_chars = set("BAR")

    # TODO: duplicate of test_definite_chars, remove when nondegenerate_chars,
    # is removed
    def test_nondegenerate_chars(self):
        expected = set("ABCQ")
        self.assertEqual(ExampleGrammaredSequence.nondegenerate_chars,
                         expected)

        ExampleGrammaredSequence.degenerate_chars.add("D")
        self.assertEqual(ExampleGrammaredSequence.nondegenerate_chars,
                         expected)

        self.assertEqual(ExampleGrammaredSequence('').nondegenerate_chars,
                         expected)

        with self.assertRaises(AttributeError):
            ExampleGrammaredSequence('').nondegenerate_chars = set("BAR")

    def test_definite_chars(self):
        expected = set("ABCQ")
        self.assertEqual(ExampleGrammaredSequence.definite_chars,
                         expected)

        ExampleGrammaredSequence.degenerate_chars.add("D")
        self.assertEqual(ExampleGrammaredSequence.definite_chars,
                         expected)

        self.assertEqual(ExampleGrammaredSequence('').definite_chars,
                         expected)

        with self.assertRaises(AttributeError):
            ExampleGrammaredSequence('').definite_chars = set("BAR")

    def test_gap_chars(self):
        expected = set(".-")
        self.assertIs(type(ExampleGrammaredSequence.gap_chars), set)
        self.assertEqual(ExampleGrammaredSequence.gap_chars, expected)

        ExampleGrammaredSequence.gap_chars.add("_")
        self.assertEqual(ExampleGrammaredSequence.gap_chars, expected)

        self.assertEqual(ExampleGrammaredSequence('').gap_chars, expected)

        with self.assertRaises(AttributeError):
            ExampleGrammaredSequence('').gap_chars = set("_ =")

    def test_default_gap_char(self):
        self.assertIs(type(ExampleGrammaredSequence.default_gap_char), str)
        self.assertEqual(ExampleGrammaredSequence.default_gap_char, '-')
        self.assertEqual(ExampleGrammaredSequence('').default_gap_char, '-')

        with self.assertRaises(AttributeError):
            ExampleGrammaredSequence('').default_gap_char = '.'

    def test_alphabet(self):
        expected = set("ABC.-XYZQW")
        self.assertIs(type(ExampleGrammaredSequence.alphabet), set)
        self.assertEqual(ExampleGrammaredSequence.alphabet, expected)

        ExampleGrammaredSequence.alphabet.add("DEF")
        self.assertEqual(ExampleGrammaredSequence.alphabet, expected)

        self.assertEqual(ExampleGrammaredSequence('').alphabet, expected)

        with self.assertRaises(AttributeError):
            ExampleGrammaredSequence('').alphabet = set("ABCDEFG.-WXYZ")

    def test_degenerate_map(self):
        expected = {"X": set("AB"), "Y": set("BC"), "Z": set("AC"), "W": set("ABCQ")}
        self.assertEqual(ExampleGrammaredSequence.degenerate_map, expected)

        ExampleGrammaredSequence.degenerate_map['W'] = set("ABC")
        ExampleGrammaredSequence.degenerate_map['X'] = set("CA")
        self.assertEqual(ExampleGrammaredSequence.degenerate_map, expected)

        self.assertEqual(ExampleGrammaredSequence('').degenerate_map, expected)

        with self.assertRaises(AttributeError):
            ExampleGrammaredSequence('').degenerate_map = {'W': "ABC"}

    def test_gaps(self):
        self.assertIs(type(ExampleGrammaredSequence("").gaps()), np.ndarray)
        self.assertIs(ExampleGrammaredSequence("").gaps().dtype,
                      np.dtype('bool'))
        npt.assert_equal(ExampleGrammaredSequence("ABCXBZYABC").gaps(),
                         np.zeros(10).astype(bool))

        npt.assert_equal(ExampleGrammaredSequence(".-.-.").gaps(),
                         np.ones(5).astype(bool))

        npt.assert_equal(ExampleGrammaredSequence("A.B-C.X-Y.").gaps(),
                         np.array([0, 1] * 5, dtype=bool))

        npt.assert_equal(ExampleGrammaredSequence("AB.AC.XY-").gaps(),
                         np.array([0, 0, 1] * 3, dtype=bool))

        npt.assert_equal(ExampleGrammaredSequence("A.BC.-").gaps(),
                         np.array([0, 1, 0, 0, 1, 1], dtype=bool))

    def test_has_gaps(self):
        self.assertIs(type(ExampleGrammaredSequence("").has_gaps()), bool)
        self.assertIs(type(ExampleGrammaredSequence("-").has_gaps()), bool)

        self.assertFalse(ExampleGrammaredSequence("").has_gaps())
        self.assertFalse(ExampleGrammaredSequence("ABCXYZ").has_gaps())

        self.assertTrue(ExampleGrammaredSequence("-").has_gaps())
        self.assertTrue(ExampleGrammaredSequence("ABCXYZ-").has_gaps())

    def test_degenerates(self):
        self.assertIs(type(ExampleGrammaredSequence("").degenerates()),
                      np.ndarray)
        self.assertIs(ExampleGrammaredSequence("").degenerates().dtype,
                      np.dtype('bool'))

        npt.assert_equal(ExampleGrammaredSequence("ABCBC-.AB.").degenerates(),
                         np.zeros(10).astype(bool))

        npt.assert_equal(ExampleGrammaredSequence("ZYZYZ").degenerates(),
                         np.ones(5).astype(bool))

        npt.assert_equal(ExampleGrammaredSequence("AX.Y-ZBXCZ").degenerates(),
                         np.array([0, 1] * 5, dtype=bool))

        npt.assert_equal(ExampleGrammaredSequence("ABXACY.-Z").degenerates(),
                         np.array([0, 0, 1] * 3, dtype=bool))

        npt.assert_equal(ExampleGrammaredSequence("AZBCXY").degenerates(),
                         np.array([0, 1, 0, 0, 1, 1], dtype=bool))

    def test_has_degenerates(self):
        self.assertIs(type(ExampleGrammaredSequence("").has_degenerates()),
                      bool)
        self.assertIs(type(ExampleGrammaredSequence("X").has_degenerates()),
                      bool)

        self.assertFalse(ExampleGrammaredSequence("").has_degenerates())
        self.assertFalse(ExampleGrammaredSequence("A-.BC").has_degenerates())

        self.assertTrue(ExampleGrammaredSequence("Z").has_degenerates())
        self.assertTrue(ExampleGrammaredSequence("ABC.XYZ-").has_degenerates())

    # TODO: duplicate of test_definites; remove when nondegenerates is removed
    def test_nondegenerates(self):
        self.assertIs(type(ExampleGrammaredSequence("").nondegenerates()),
                      np.ndarray)
        self.assertIs(ExampleGrammaredSequence("").nondegenerates().dtype,
                      np.dtype('bool'))

        npt.assert_equal(
            ExampleGrammaredSequence("XYZYZ-.XY.").nondegenerates(),
            np.zeros(10).astype(bool))

        npt.assert_equal(ExampleGrammaredSequence("ABABA").nondegenerates(),
                         np.ones(5).astype(bool))

        npt.assert_equal(
            ExampleGrammaredSequence("XA.B-AZCXA").nondegenerates(),
            np.array([0, 1] * 5, dtype=bool))

        npt.assert_equal(
            ExampleGrammaredSequence("XXAZZB.-C").nondegenerates(),
            np.array([0, 0, 1] * 3, dtype=bool))

        npt.assert_equal(ExampleGrammaredSequence("YB.-AC").nondegenerates(),
                         np.array([0, 1, 0, 0, 1, 1], dtype=bool))

    def test_definites(self):
        self.assertIs(type(ExampleGrammaredSequence("").definites()),
                      np.ndarray)
        self.assertIs(ExampleGrammaredSequence("").definites().dtype,
                      np.dtype('bool'))

        npt.assert_equal(
            ExampleGrammaredSequence("XYZYZ-.XY.").definites(),
            np.zeros(10).astype(bool))

        npt.assert_equal(ExampleGrammaredSequence("ABABA").definites(),
                         np.ones(5).astype(bool))

        npt.assert_equal(
            ExampleGrammaredSequence("XA.B-AZCXA").definites(),
            np.array([0, 1] * 5, dtype=bool))

        npt.assert_equal(
            ExampleGrammaredSequence("XXAZZB.-C").definites(),
            np.array([0, 0, 1] * 3, dtype=bool))

        npt.assert_equal(ExampleGrammaredSequence("YB.-AC").definites(),
                         np.array([0, 1, 0, 0, 1, 1], dtype=bool))

    # TODO: duplicate of test_has_definites; remove when has_nondegenerates is
    # removed.
    def test_has_nondegenerates(self):
        self.assertIs(type(ExampleGrammaredSequence("").has_nondegenerates()),
                      bool)
        self.assertIs(type(ExampleGrammaredSequence("A").has_nondegenerates()),
                      bool)

        self.assertFalse(ExampleGrammaredSequence("").has_nondegenerates())
        self.assertFalse(
            ExampleGrammaredSequence("X-.YZ").has_nondegenerates())

        self.assertTrue(ExampleGrammaredSequence("C").has_nondegenerates())
        self.assertTrue(
            ExampleGrammaredSequence(".XYZ-ABC").has_nondegenerates())

    def test_has_definites(self):
        self.assertIs(type(ExampleGrammaredSequence("").has_definites()),
                      bool)
        self.assertIs(type(ExampleGrammaredSequence("A").has_definites()),
                      bool)

        self.assertFalse(ExampleGrammaredSequence("").has_definites())
        self.assertFalse(
            ExampleGrammaredSequence("X-.YZ").has_definites())

        self.assertTrue(ExampleGrammaredSequence("C").has_definites())
        self.assertTrue(
            ExampleGrammaredSequence(".XYZ-ABC").has_definites())

    def test_degap(self):
        kw = {
            'metadata': {
                'id': 'some_id',
                'description': 'some description',
            },
        }

        self.assertEqual(
            ExampleGrammaredSequence(
                "", positional_metadata={'qual': []}, **kw).degap(),
            ExampleGrammaredSequence(
                "", positional_metadata={'qual': []}, **kw))

        self.assertEqual(
            ExampleGrammaredSequence(
                "ABCXYZ",
                positional_metadata={'qual': np.arange(6)},
                **kw).degap(),
            ExampleGrammaredSequence(
                "ABCXYZ",
                positional_metadata={'qual': np.arange(6)},
                **kw))

        self.assertEqual(
            ExampleGrammaredSequence(
                "ABC-XYZ",
                positional_metadata={'qual': np.arange(7, dtype=np.int64)},
                **kw).degap(),
            ExampleGrammaredSequence(
                "ABCXYZ",
                positional_metadata={'qual': np.asarray([0, 1, 2, 4, 5, 6],
                                                        dtype=np.int64)},
                **kw))

        self.assertEqual(
            ExampleGrammaredSequence(
                ".-ABC-XYZ.",
                positional_metadata={'qual': np.arange(10, dtype=np.int64)},
                **kw).degap(),
            ExampleGrammaredSequence(
                "ABCXYZ",
                positional_metadata={'qual': np.asarray([2, 3, 4, 6, 7, 8],
                                                        dtype=np.int64)},
                **kw))

        self.assertEqual(
            ExampleGrammaredSequence(
                "---.-.-.-.-.",
                positional_metadata={'quality': np.arange(12, dtype=np.int64)},
                **kw).degap(),
            ExampleGrammaredSequence(
                "",
                positional_metadata={'quality': np.array([], dtype=np.int64)},
                **kw))

    def test_expand_degenerates_no_degens(self):
        seq = ExampleGrammaredSequence("ABCABCABC")
        self.assertEqual(list(seq.expand_degenerates()), [seq])

    def test_expand_degenerates_all_degens(self):
        exp = [
            ExampleGrammaredSequence('ABA'), ExampleGrammaredSequence('ABC'),
            ExampleGrammaredSequence('ACA'), ExampleGrammaredSequence('ACC'),
            ExampleGrammaredSequence('BBA'), ExampleGrammaredSequence('BBC'),
            ExampleGrammaredSequence('BCA'), ExampleGrammaredSequence('BCC')
        ]
        # Sort based on sequence string, as order is not guaranteed.
        obs = sorted(ExampleGrammaredSequence('XYZ').expand_degenerates(),
                     key=str)
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
        exp = [ExampleGrammaredSequence('ABA', **kw),
               ExampleGrammaredSequence('ABC', **kw),
               ExampleGrammaredSequence('BBA', **kw),
               ExampleGrammaredSequence('BBC', **kw)]
        obs = sorted(
            ExampleGrammaredSequence('XBZ', **kw).expand_degenerates(),
            key=str)
        self.assertEqual(obs, exp)

    def test_to_regex_no_degens(self):
        seq = ExampleGrammaredSequence('ABC')
        regex = seq.to_regex()
        self.assertEqual(regex.pattern, str(seq))

    def test_to_regex_with_degens(self):
        seq = ExampleGrammaredSequence('AYZ')
        regex = seq.to_regex()
        self.assertFalse(any(regex.match(s) is None
                             for s in 'ABA ABC ACA ACC'.split()))
        self.assertTrue(all(regex.match(s) is None
                            for s in 'CBA BBA ABB AAA'.split()))

    def test_to_regex_within_capture(self):
        seq = ExampleGrammaredSequence('XYC')
        regex = seq.to_regex(within_capture=True)

        for ref in 'ABA BBB CCA'.split():
            self.assertFalse(any(len(match.groups()) == 1
                                 for match in regex.finditer(ref)))

        for ref in 'ABC BBC ACC'.split():
            self.assertTrue(all(len(match.groups()) == 1
                                for match in regex.finditer(ref)))

    def test_find_motifs_no_motif(self):
        seq = ExampleMotifsTester("ABCABCABC")
        with self.assertRaises(ValueError) as cm:
            seq.find_motifs("doesn't-exist")
        self.assertIn("doesn't-exist", str(cm.exception))

        seq = ExampleGrammaredSequence("ABCABCABC")
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
        obs = repr(ExampleGrammaredSequence(''))
        self.assertEqual(obs.count('\n'), 7)
        self.assertTrue(obs.startswith('ExampleGrammaredSequence'))
        self.assertIn('length: 0', obs)
        self.assertIn('has gaps: False', obs)
        self.assertIn('has degenerates: False', obs)
        self.assertIn('has definites: False', obs)
        self.assertTrue(obs.endswith('-'))

        # no metadata, mix of gaps, degenerates, and definites
        obs = repr(ExampleGrammaredSequence('AY-B'))
        self.assertEqual(obs.count('\n'), 8)
        self.assertTrue(obs.startswith('ExampleGrammaredSequence'))
        self.assertIn('length: 4', obs)
        self.assertIn('has gaps: True', obs)
        self.assertIn('has degenerates: True', obs)
        self.assertIn('has definites: True', obs)
        self.assertTrue(obs.endswith('0 AY-B'))

        # metadata and positional metadata of mixed types
        obs = repr(
            ExampleGrammaredSequence(
                'ABCA',
                metadata={'foo': 42, b'bar': 33.33, None: True, False: {},
                          (1, 2): 3, 'acb' * 100: "'"},
                positional_metadata={'foo': range(4),
                                     42: ['a', 'b', [], 'c']}))
        self.assertEqual(obs.count('\n'), 18)
        self.assertTrue(obs.startswith('ExampleGrammaredSequence'))
        self.assertIn('None: True', obs)
        self.assertIn('\'foo\': 42', obs)
        self.assertIn('42: <dtype: object>', obs)
        self.assertIn('\'foo\': <dtype: int64>', obs)
        self.assertIn('length: 4', obs)
        self.assertIn('has gaps: False', obs)
        self.assertIn('has degenerates: False', obs)
        self.assertIn('has definites: True', obs)
        self.assertTrue(obs.endswith('0 ABCA'))

        # sequence spanning > 5 lines
        obs = repr(ExampleGrammaredSequence('A' * 301))
        self.assertEqual(obs.count('\n'), 12)
        self.assertTrue(obs.startswith('ExampleGrammaredSequence'))
        self.assertIn('length: 301', obs)
        self.assertIn('has gaps: False', obs)
        self.assertIn('has degenerates: False', obs)
        self.assertIn('has definites: True', obs)
        self.assertIn('...', obs)
        self.assertTrue(obs.endswith('300 A'))

    def test_to_definites(self):
        seq = ExampleGrammaredSequence("ABCQXYZ")

        # default behavior, here I expect to see the sequence "ABCWWWW" returned
        obs = seq.to_definites()
        exp = ExampleGrammaredSequence("ABCWWWW")
        self.assertEqual(obs, exp)

        # noncanonical wildcard, expect to see "ABCQWWW" returned
        obs = seq.to_definites(noncanonical=False)
        exp = ExampleGrammaredSequence("ABCQWWW")

        # gap behavior, I expect to see the sequence "ABC----" returned
        obs = seq.to_definites(degenerate="gap")
        exp = ExampleGrammaredSequence("ABC----")
        self.assertEqual(obs, exp)

        # noncanonical gap
        obs = seq.to_definites(degenerate="gap", noncanonical=False)
        exp = ExampleGrammaredSequence("ABCQ---")
        self.assertEqual(obs, exp)

        # canonical trim
        obs = seq.to_definites(degenerate="del")
        exp = ExampleGrammaredSequence("ABC")
        self.assertEqual(obs, exp)

        # noncanonical trim
        obs = seq.to_definites(degenerate="del", noncanonical=False)
        exp = ExampleGrammaredSequence("ABCQ")
        self.assertEqual(obs, exp)

        # single char, acceptable input
        obs = seq.to_definites(degenerate="A")
        exp = ExampleGrammaredSequence("ABCAAAA")
        self.assertEqual(obs, exp)

        # noncanonical single char, acceptable input
        obs = seq.to_definites(degenerate="A", noncanonical=False)
        exp = ExampleGrammaredSequence("ABCQAAA")
        self.assertEqual(obs, exp)

        # test that single char outside of alphabet will throw error
        with self.assertRaises(ValueError):
            seq.to_definites("P")

        # test that an invalid wildcard (not a string) will throw an error
        ExampleGrammaredSequence.wildcard_char = 1
        with self.assertRaises(ValueError):
            seq.to_definites()
        ExampleGrammaredSequence.wildcard_char = 'W'

        # test that nonsense input for 'to' will throw error
        with self.assertRaises(ValueError):
            seq.to_definites(degenerate='nonsense')

    def test_noncanonical_chars(self):
        self.assertTrue(isinstance(GrammaredSequence.noncanonical_chars, set))
        self.assertEqual(len(GrammaredSequence.noncanonical_chars), 0)
    
    def test_wildcard_char(self):
        exp = None
        self.assertEqual(GrammaredSequence.wildcard_char, exp)

if __name__ == "__main__":
    main()
