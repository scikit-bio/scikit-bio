# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

from skbio.sequence._nucleotide_sequence import NucleotideSequence
from skbio.util import classproperty


class ExampleNucleotideSequence(NucleotideSequence):
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


class TestNucelotideSequence(unittest.TestCase):
    def test_instantiation_with_no_implementation(self):
        class NucleotideSequenceSubclassNoImplementation(NucleotideSequence):
            pass

        with self.assertRaises(TypeError) as cm:
            NucleotideSequenceSubclassNoImplementation()

        self.assertIn("abstract class", str(cm.exception))
        self.assertIn("nondegenerate_chars", str(cm.exception))
        self.assertIn("degenerate_map", str(cm.exception))
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

        self.assertEqual(ExampleNucleotideSequence.complement_map, expected)

        ExampleNucleotideSequence.complement_map['W'] = 'X'
        ExampleNucleotideSequence.complement_map['X'] = 'W'
        self.assertEqual(ExampleNucleotideSequence.complement_map, expected)

        self.assertEqual(ExampleNucleotideSequence('').complement_map, expected)

        with self.assertRaises(AttributeError):
            ExampleNucleotideSequence('').complement_map = {'W': 'X'}


# class NucelotideSequenceTests(TestCase):
#
#    def setUp(self):
#        self.empty = NucleotideSequence('')
#        self.b1 = NucleotideSequence('GATTACA')
#        self.b2 = NucleotideSequence(
#            'ACCGGUACC', id="test-seq-2",
#            description="A test sequence")
#        self.b3 = NucleotideSequence('G-AT-TG.AT.T')
#
#    def test_alphabet(self):
#        exp = {
#            'A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'U', 'T',
#            'W', 'V', 'Y', 'a', 'c', 'b', 'd', 'g', 'h', 'k', 'm', 'n', 's',
#            'r', 'u', 't', 'w', 'v', 'y'
#        }
#
#        # Test calling from an instance and purely static context.
#        self.assertEqual(self.b1.alphabet, exp)
#        self.assertEqual(NucleotideSequence.alphabet, exp)
#
#    def test_gap_chars(self):
#        self.assertEqual(self.b1.gap_chars, set('-.'))
#
#
#    def test_nondegenerate_chars(self):
#        exp = set("ACGTUacgtu")
#        self.assertEqual(self.b1.nondegenerate_chars, exp)
#        self.assertEqual(NucleotideSequence.nondegenerate_chars, exp)
#
#    def test_degenerate_map(self):
#        exp = {
#            # upper
#            'B': set(['C', 'U', 'T', 'G']), 'D': set(['A', 'U', 'T', 'G']),
#            'H': set(['A', 'C', 'U', 'T']), 'K': set(['U', 'T', 'G']),
#            'M': set(['A', 'C']), 'N': set(['A', 'C', 'U', 'T', 'G']),
#            'S': set(['C', 'G']), 'R': set(['A', 'G']),
#            'W': set(['A', 'U', 'T']), 'V': set(['A', 'C', 'G']),
#            'Y': set(['C', 'U', 'T']),
#            # lower
#            'b': set(['c', 'u', 't', 'g']), 'd': set(['a', 'u', 't', 'g']),
#            'h': set(['a', 'c', 'u', 't']), 'k': set(['u', 't', 'g']),
#            'm': set(['a', 'c']), 'n': set(['a', 'c', 'u', 't', 'g']),
#            's': set(['c', 'g']), 'r': set(['a', 'g']),
#            'w': set(['a', 'u', 't']), 'v': set(['a', 'c', 'g']),
#            'y': set(['c', 'u', 't'])
#        }
#        self.assertEqual(self.b1.degenerate_map, exp)
#        self.assertEqual(NucleotideSequence.degenerate_map, exp)
#
#        # Test that we can modify a copy of the mapping without altering the
#        # canonical representation.
#        degen = NucleotideSequence.degenerate_map
#        degen.update({'V': set("BRO"), 'Z': set("ZORRO")})
#        self.assertNotEqual(degen, exp)
#        self.assertEqual(NucleotideSequence.degenerate_map, exp)
#
#    def test_degenerate_chars(self):
#        exp = set(['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y',
#                   'b', 'd', 'h', 'k', 'm', 'n', 's', 'r', 'w', 'v', 'y'])
#        self.assertEqual(self.b1.degenerate_chars, exp)
#        self.assertEqual(NucleotideSequence.degenerate_chars, exp)
#
#    def test_complement(self):
#        self.assertRaises(SequenceError,
#                          self.b1.complement)
#
#    def test_reverse_complement(self):
#        self.assertRaises(SequenceError,
#                          self.b1.reverse_complement)
#
#    def test_is_reverse_complement(self):
#        self.assertRaises(SequenceError,
#                          self.b1.is_reverse_complement, self.b1)
#
#    def test_nondegenerates_invalid(self):
#        with self.assertRaises(SequenceError):
#            list(NucleotideSequence('AZA').nondegenerates())
#
#    def test_nondegenerates_empty(self):
#        self.assertEqual(list(self.empty.nondegenerates()), [self.empty])
#
#    def test_nondegenerates_no_degens(self):
#        self.assertEqual(list(self.b1.nondegenerates()), [self.b1])
#
#    def test_nondegenerates_all_degens(self):
#        # Same chars.
#        exp = [NucleotideSequence('CC'), NucleotideSequence('CG'),
#               NucleotideSequence('GC'), NucleotideSequence('GG')]
#        # Sort based on sequence string, as order is not guaranteed.
#        obs = sorted(NucleotideSequence('SS').nondegenerates(), key=str)
#        self.assertEqual(obs, exp)
#
#        # Different chars.
#        exp = [NucleotideSequence('AC'), NucleotideSequence('AG'),
#               NucleotideSequence('GC'), NucleotideSequence('GG')]
#        obs = sorted(NucleotideSequence('RS').nondegenerates(), key=str)
#        self.assertEqual(obs, exp)
#
#        # Odd number of chars.
#        obs = list(NucleotideSequence('NNN').nondegenerates())
#        self.assertEqual(len(obs), 5**3)
#
#    def test_nondegenerates_mixed_degens(self):
#        exp = [NucleotideSequence('AGC'), NucleotideSequence('AGT'),
#               NucleotideSequence('AGU'), NucleotideSequence('GGC'),
#               NucleotideSequence('GGT'), NucleotideSequence('GGU')]
#        obs = sorted(NucleotideSequence('RGY').nondegenerates(), key=str)
#        self.assertEqual(obs, exp)
#
#    def test_nondegenerates_gap_mixed_case(self):
#        exp = [NucleotideSequence('-A.a'), NucleotideSequence('-A.c'),
#               NucleotideSequence('-C.a'), NucleotideSequence('-C.c')]
#        obs = sorted(NucleotideSequence('-M.m').nondegenerates(), key=str)
#        self.assertEqual(obs, exp)
#
#    def test_find_features(self):
#        exp = [(0, 2, 'GA'), (4, 5, 'A'), (6, 7, 'A')]
#        obs = list(self.b1.find_features('purine_run'))
#        self.assertEqual(obs, exp)
#
#        exp = [(2, 4, 'TT'), (5, 6, 'C')]
#        obs = list(self.b1.find_features('pyrimidine_run'))
#        self.assertEqual(obs, exp)
#
#        exp = [(0, 1, 'A'), (3, 5, 'GG'), (6, 7, 'A')]
#        obs = list(self.b2.find_features('purine_run'))
#        self.assertEqual(obs, exp)
#
#        exp = [(1, 3, 'CC'), (5, 6, 'U'), (7, 9, 'CC')]
#        obs = list(self.b2.find_features('pyrimidine_run'))
#        self.assertEqual(obs, exp)
#
#    def test_find_features_min_length(self):
#        exp = [(0, 2, 'GA')]
#        obs = list(self.b1.find_features('purine_run', 2))
#        self.assertEqual(obs, exp)
#
#        exp = [(2, 4, 'TT')]
#        obs = list(self.b1.find_features('pyrimidine_run', 2))
#        self.assertEqual(obs, exp)
#
#        exp = [(3, 5, 'GG')]
#        obs = list(self.b2.find_features('purine_run', 2))
#        self.assertEqual(obs, exp)
#
#        exp = [(1, 3, 'CC'), (7, 9, 'CC')]
#        obs = list(self.b2.find_features('pyrimidine_run', 2))
#        self.assertEqual(obs, exp)
#
#    def test_find_features_no_feature_type(self):
#        with self.assertRaises(ValueError):
#            list(self.b1.find_features('nonexistent_feature_type'))
#
#    def test_find_features_allow_gaps(self):
#        exp = [(0, 3, 'G-A'), (6, 9, 'G.A')]
#        obs = list(self.b3.find_features('purine_run', 2, True))
#        self.assertEqual(obs, exp)
#
#        exp = [(3, 6, 'T-T'), (9, 12, 'T.T')]
#        obs = list(self.b3.find_features('pyrimidine_run', 2, True))
#        self.assertEqual(obs, exp)
#
#    def test_nondegenerates_propagate_optional_properties(self):
#        seq = NucleotideSequence('RS', id='foo', description='bar',
#                                 quality=[42, 999])
#
#        exp = [
#            NucleotideSequence('AC', id='foo', description='bar',
#                               quality=[42, 999]),
#            NucleotideSequence('AG', id='foo', description='bar',
#                               quality=[42, 999]),
#            NucleotideSequence('GC', id='foo', description='bar',
#                               quality=[42, 999]),
#            NucleotideSequence('GG', id='foo', description='bar',
#                               quality=[42, 999])
#        ]
#
#        obs = sorted(seq.nondegenerates(), key=str)
#
#        for o, e in zip(obs, exp):
#            # use equals method to ensure that id, description, and quality
#            # are correctly propagated to the resulting sequence
#            self.assertTrue(o.equals(e))


if __name__ == "__main__":
    unittest.main()
