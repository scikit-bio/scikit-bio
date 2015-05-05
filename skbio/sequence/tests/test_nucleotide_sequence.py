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

        self.assertEqual(ExampleNucleotideSequence('').complement_map,
                         expected)

        with self.assertRaises(AttributeError):
            ExampleNucleotideSequence('').complement_map = {'W': 'X'}

    def test_complement_without_reverse_empty(self):
        # without optional attributes
        comp = ExampleNucleotideSequence('').complement()
        self.assertEqual(comp, ExampleNucleotideSequence(''))

        # with optional attributes
        comp = ExampleNucleotideSequence('', id='foo', description='bar',
                                         quality=[]).complement()
        self.assertEqual(
            comp,
            ExampleNucleotideSequence('', id='foo', description='bar',
                                      quality=[]))

    def test_complement_without_reverse_non_empty(self):
        comp = ExampleNucleotideSequence('ABCXYZ.-BBZ').complement()
        self.assertEqual(comp, ExampleNucleotideSequence('CBAYXZ.-BBZ'))

        comp = ExampleNucleotideSequence(
            'ABCXYZ.-BBZ', id='foo', description='bar',
            quality=range(11)).complement()
        self.assertEqual(
            comp,
            ExampleNucleotideSequence('CBAYXZ.-BBZ', id='foo',
                                      description='bar', quality=range(11)))

    def test_complement_with_reverse_empty(self):
        rc = ExampleNucleotideSequence('').complement(reverse=True)
        self.assertEqual(rc, ExampleNucleotideSequence(''))

        rc = ExampleNucleotideSequence('', id='foo', description='bar',
                                       quality=[]).complement(reverse=True)
        self.assertEqual(
            rc,
            ExampleNucleotideSequence('', id='foo', description='bar',
                                      quality=[]))

    def test_complement_with_reverse_non_empty(self):
        rc = ExampleNucleotideSequence('ABCXYZ.-BBZ').complement(reverse=True)
        self.assertEqual(rc, ExampleNucleotideSequence('ZBB-.ZXYABC'))

        rc = ExampleNucleotideSequence(
            'ABCXYZ.-BBZ', id='foo', description='bar',
            quality=range(11)).complement(reverse=True)
        self.assertEqual(
            rc,
            ExampleNucleotideSequence('ZBB-.ZXYABC', id='foo',
                                      description='bar',
                                      quality=list(range(11))[::-1]))

    def test_reverse_complement(self):
        # light tests because this just calls
        # NucleotideSequence.complement(reverse=True), which is tested more
        # extensively
        rc = ExampleNucleotideSequence(
            'ABCXYZ.-BBZ', id='foo', description='bar',
            quality=range(11)).reverse_complement()
        self.assertEqual(
            rc,
            ExampleNucleotideSequence('ZBB-.ZXYABC', id='foo',
                                      description='bar',
                                      quality=list(range(11))[::-1]))


# class NucelotideSequenceTests(TestCase):
#
#    def setUp(self):
#        self.empty = NucleotideSequence('')
#        self.b1 = NucleotideSequence('GATTACA')
#        self.b2 = NucleotideSequence(
#            'ACCGGUACC', id="test-seq-2",
#            description="A test sequence")
#        self.b3 = NucleotideSequence('G-AT-TG.AT.T')
##
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
