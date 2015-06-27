# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

import numpy as np

from skbio.sequence import IUPACSequence, NucleotideMixin
from skbio.util import classproperty


class ExampleNucleotideMixin(IUPACSequence, NucleotideMixin):
    @classproperty
    def degenerate_map(cls):
        return {"X": set("AB"), "Y": set("BC"), "Z": set("AC")}

    @classproperty
    def nondegenerate_chars(cls):
        return set("ABC")

    @classproperty
    def complement_map(cls):
        comp_map = {
            'A': 'C', 'C': 'A',
            'B': 'B',
            'X': 'Y', 'Y': 'X',
            'Z': 'Z'
        }

        comp_map.update({c: c for c in cls.gap_chars})
        return comp_map


class ExampleNucleotideMixinSubclass(ExampleNucleotideMixin):
    pass


class TestNucelotideSequence(unittest.TestCase):
    def setUp(self):
        self.sequence_kinds = frozenset([
            str,
            ExampleNucleotideMixin,
            lambda s: np.fromstring(s, dtype='|S1'),
            lambda s: np.fromstring(s, dtype=np.uint8)])

    def test_instantiation_with_no_implementation(self):
        class NucleotideSequenceSubclassNoImplementation(NucleotideMixin):
            pass

        with self.assertRaises(TypeError) as cm:
            NucleotideSequenceSubclassNoImplementation()

        self.assertIn("abstract class", str(cm.exception))
        self.assertIn("complement_map", str(cm.exception))

    def test_complement_map(self):
        expected = {
            'A': 'C', 'C': 'A',
            'B': 'B',
            'X': 'Y', 'Y': 'X',
            'Z': 'Z',
            '.': '.',
            '-': '-'
        }

        self.assertEqual(ExampleNucleotideMixin.complement_map, expected)

        ExampleNucleotideMixin.complement_map['W'] = 'X'
        ExampleNucleotideMixin.complement_map['X'] = 'W'
        self.assertEqual(ExampleNucleotideMixin.complement_map, expected)

        self.assertEqual(ExampleNucleotideMixin('').complement_map,
                         expected)

        with self.assertRaises(AttributeError):
            ExampleNucleotideMixin('').complement_map = {'W': 'X'}

    def test_complement_without_reverse_empty(self):
        # without optional attributes
        comp = ExampleNucleotideMixin('').complement()
        self.assertEqual(comp, ExampleNucleotideMixin(''))

        # with optional attributes
        comp = ExampleNucleotideMixin(
            '',
            metadata={'id': 'foo', 'description': 'bar'},
            positional_metadata={'quality': []}).complement()
        self.assertEqual(
            comp,
            ExampleNucleotideMixin(
                '',
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': []}))

    def test_complement_without_reverse_non_empty(self):
        comp = ExampleNucleotideMixin('ABCXYZ.-BBZ').complement()
        self.assertEqual(comp, ExampleNucleotideMixin('CBAYXZ.-BBZ'))

        comp = ExampleNucleotideMixin(
            'ABCXYZ.-BBZ',
            metadata={'id': 'foo', 'description': 'bar'},
            positional_metadata={'quality': range(11)}).complement()
        self.assertEqual(
            comp,
            ExampleNucleotideMixin(
                'CBAYXZ.-BBZ',
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': range(11)}))

    def test_complement_with_reverse_empty(self):
        rc = ExampleNucleotideMixin('').complement(reverse=True)
        self.assertEqual(rc, ExampleNucleotideMixin(''))

        rc = ExampleNucleotideMixin(
            '',
            metadata={'id': 'foo', 'description': 'bar'},
            positional_metadata={'quality': []}).complement(reverse=True)
        self.assertEqual(
            rc,
            ExampleNucleotideMixin(
                '',
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': []}))

    def test_complement_with_reverse_non_empty(self):
        rc = ExampleNucleotideMixin('ABCXYZ.-BBZ').complement(reverse=True)
        self.assertEqual(rc, ExampleNucleotideMixin('ZBB-.ZXYABC'))

        rc = ExampleNucleotideMixin(
            'ABCXYZ.-BBZ',
            metadata={'id': 'foo', 'description': 'bar'},
            positional_metadata={'quality':
                                 range(11)}).complement(reverse=True)
        self.assertEqual(
            rc,
            ExampleNucleotideMixin(
                'ZBB-.ZXYABC',
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': list(range(11))[::-1]}))

    def test_reverse_complement(self):
        # light tests because this just calls
        # NucleotideSequence.complement(reverse=True), which is tested more
        # extensively
        rc = ExampleNucleotideMixin(
            'ABCXYZ.-BBZ',
            metadata={'id': 'foo', 'description': 'bar'},
            positional_metadata={'quality': range(11)}).reverse_complement()
        self.assertEqual(
            rc,
            ExampleNucleotideMixin(
                'ZBB-.ZXYABC',
                metadata={'id': 'foo', 'description': 'bar'},
                positional_metadata={'quality': list(range(11))[::-1]}))

    def test_is_reverse_complement_varied_types(self):
        tested = 0
        for constructor in self.sequence_kinds:
            tested += 1
            seq1 = ExampleNucleotideMixin('ABCXYZ.-BBZ')
            seq2 = constructor('ZBB-.ZXYABC')

            self.assertTrue(seq1.is_reverse_complement(seq2))

        self.assertEqual(tested, 4)

    def test_is_reverse_complement_empty(self):
        seq1 = ExampleNucleotideMixin('')
        self.assertTrue(seq1.is_reverse_complement(seq1))

        # optional attributes are ignored, only the sequence is compared
        seq2 = ExampleNucleotideMixin(
            '',
            metadata={'id': 'foo', 'description': 'bar'},
            positional_metadata={'quality':
                                 np.array([], dtype=np.int64)})
        self.assertTrue(seq2.is_reverse_complement(seq2))
        self.assertTrue(seq1.is_reverse_complement(seq2))
        self.assertTrue(seq2.is_reverse_complement(seq1))

    def test_is_reverse_complement_metadata_ignored(self):
        seq1 = ExampleNucleotideMixin('ABCXYZ.-BBZ')
        seq2 = ExampleNucleotideMixin(
            'ZBB-.ZXYABC',
            metadata={'id': 'foo', 'description': 'bar'},
            positional_metadata={'quality': range(11)})

        self.assertFalse(seq1.is_reverse_complement(seq1))
        self.assertFalse(seq2.is_reverse_complement(seq2))

        self.assertTrue(seq1.is_reverse_complement(seq2))
        self.assertTrue(seq2.is_reverse_complement(seq1))

    def test_is_reverse_complement_non_reverse_complements(self):
        # same length
        seq1 = ExampleNucleotideMixin('AABC')
        seq2 = ExampleNucleotideMixin('ABCX')

        self.assertFalse(seq1.is_reverse_complement(seq1))
        self.assertFalse(seq2.is_reverse_complement(seq2))

        self.assertFalse(seq1.is_reverse_complement(seq2))
        self.assertFalse(seq2.is_reverse_complement(seq1))

        # different length
        seq1 = ExampleNucleotideMixin('AABC')
        seq2 = ExampleNucleotideMixin('ABCXZ')

        self.assertFalse(seq1.is_reverse_complement(seq1))
        self.assertFalse(seq2.is_reverse_complement(seq2))

        self.assertFalse(seq1.is_reverse_complement(seq2))
        self.assertFalse(seq2.is_reverse_complement(seq1))

    def test_is_reverse_complement_type_mismatch(self):
        seq1 = ExampleNucleotideMixin('ABC')
        seq2 = ExampleNucleotideMixinSubclass('ABC')

        with self.assertRaises(TypeError):
            seq1.is_reverse_complement(seq2)


if __name__ == "__main__":
    unittest.main()
