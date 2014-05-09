#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division

from collections import Counter, defaultdict
from unittest import TestCase, main

from skbio.core.sequence import (
    BiologicalSequence, NucleotideSequence, DNASequence, RNASequence)
from skbio.core.exception import BiologicalSequenceError


class BiologicalSequenceTests(TestCase):
    """ Tests of the BiologicalSequence class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
        self.b1 = BiologicalSequence('GATTACA')
        self.b2 = BiologicalSequence(
            'ACCGGTACC', identifier="test-seq-2",
            description="A test sequence")
        self.b3 = BiologicalSequence(
            'GREG', identifier="test-seq-3", description="A protein sequence")
        self.b4 = BiologicalSequence(
            'PRTEIN', identifier="test-seq-4")
        self.b5 = BiologicalSequence(
            'LLPRTEIN', description="some description")
        self.b6 = BiologicalSequence('ACGTACGTACGT')
        self.b7 = BiologicalSequence('..--..')
        self.b8 = BiologicalSequence('HE..--..LLO')

    def test_init(self):
        """ Initialization functions as expected with varied input types
        """
        # init as string
        b = BiologicalSequence('ACCGGXZY')
        self.assertEqual(str(b), 'ACCGGXZY')
        self.assertEqual(b.identifier, "")
        self.assertEqual(b.description, "")

        # init as string with optional values
        b = BiologicalSequence(
            'ACCGGXZY', 'test-seq-1', 'The first test sequence')
        self.assertEqual(str(b), 'ACCGGXZY')
        self.assertEqual(b.identifier, "test-seq-1")
        self.assertEqual(b.description, "The first test sequence")

        # test init as a different string
        b = BiologicalSequence('WRRTY')
        self.assertEqual(str(b), 'WRRTY')

        # init as list
        b = BiologicalSequence(list('ACCGGXZY'))
        self.assertEqual(str(b), 'ACCGGXZY')
        self.assertEqual(b.identifier, "")
        self.assertEqual(b.description, "")

        # init as tuple
        b = BiologicalSequence(tuple('ACCGGXZY'))
        self.assertEqual(str(b), 'ACCGGXZY')
        self.assertEqual(b.identifier, "")
        self.assertEqual(b.description, "")

    def test_init_validate(self):
        """ initialization with validation functions as expected
        """
        self.assertRaises(BiologicalSequenceError, BiologicalSequence, "ACC",
                          validate=True)
        # no error raised when only allow characters are passed
        BiologicalSequence("..--..", validate=True)

    def test_contains(self):
        """ contains functions as expected
        """
        self.assertTrue('G' in self.b1)
        self.assertFalse('g' in self.b1)

    def test_eq(self):
        """ equality functions as expected
        """
        self.assertTrue(self.b1 == self.b1)
        self.assertTrue(self.b2 == self.b2)
        self.assertTrue(self.b3 == self.b3)

        self.assertTrue(self.b1 != self.b3)
        self.assertTrue(self.b1 != self.b2)
        self.assertTrue(self.b2 != self.b3)

        # identicial sequences of the same type are equal, even if they have
        # different identifiers and/or descriptions
        self.assertTrue(
            BiologicalSequence('ACGT') == BiologicalSequence('ACGT'))
        self.assertTrue(
            BiologicalSequence('ACGT', identifier='a') ==
            BiologicalSequence('ACGT', identifier='b'))
        self.assertTrue(
            BiologicalSequence('ACGT', description='c') ==
            BiologicalSequence('ACGT', description='d'))
        self.assertTrue(
            BiologicalSequence('ACGT', identifier='a', description='c') ==
            BiologicalSequence('ACGT', identifier='b', description='d'))

        # different type causes sequences to not be equal
        self.assertFalse(
            BiologicalSequence('ACGT') == NucleotideSequence('ACGT'))

    def test_getitem(self):
        """ getitem functions as expected
        """
        self.assertEqual(self.b1[0], BiologicalSequence('G'))
        self.assertEqual(self.b1[:], BiologicalSequence('GATTACA'))
        self.assertEqual(self.b1[::-1], BiologicalSequence('ACATTAG'))

    def test_iter(self):
        """ iter functions as expected
        """
        b1_iter = iter(self.b1)
        for actual, expected in zip(b1_iter, "GATTACA"):
            self.assertEqual(actual, expected)

        self.assertRaises(StopIteration, lambda: next(b1_iter))

    def test_k_words(self):
        """ k_words functions as expected
        """
        # overlapping = True
        self.assertEqual(list(self.b1.k_words(1, overlapping=True)),
                         ['G', 'A', 'T', 'T', 'A', 'C', 'A'])
        self.assertEqual(list(self.b1.k_words(2, overlapping=True)),
                         ['GA', 'AT', 'TT', 'TA', 'AC', 'CA'])
        self.assertEqual(list(self.b1.k_words(3, overlapping=True)),
                         ['GAT', 'ATT', 'TTA', 'TAC', 'ACA'])
        self.assertEqual(list(self.b1.k_words(7, overlapping=True)),
                         ['GATTACA'])
        self.assertEqual(list(self.b1.k_words(8, overlapping=True)),
                         [])

        # overlapping = False
        self.assertEqual(list(self.b1.k_words(1, overlapping=True)),
                         ['G', 'A', 'T', 'T', 'A', 'C', 'A'])
        self.assertEqual(list(self.b1.k_words(2, overlapping=False)),
                         ['GA', 'TT', 'AC'])
        self.assertEqual(list(self.b1.k_words(3, overlapping=False)),
                         ['GAT', 'TAC'])
        self.assertEqual(list(self.b1.k_words(7, overlapping=False)),
                         ['GATTACA'])
        self.assertEqual(list(self.b1.k_words(8, overlapping=False)),
                         [])

        # error on invalid k
        self.assertRaises(ValueError, list, self.b1.k_words(0))
        self.assertRaises(ValueError, list, self.b1.k_words(-42))

        # tests with different sequences
        self.assertEqual(list(self.b8.k_words(3, overlapping=False)),
                         ['HE.', '.--', '..L'])
        b = BiologicalSequence('')
        self.assertEqual(list(b.k_words(3)), [])

    def test_k_word_counts(self):
        """ k_word_counts functions as expected
        """
        # overlapping = True
        expected = Counter('GATTACA')
        self.assertEqual(self.b1.k_word_counts(1, overlapping=True),
                         expected)
        expected = Counter(['GAT', 'ATT', 'TTA', 'TAC', 'ACA'])
        self.assertEqual(self.b1.k_word_counts(3, overlapping=True),
                         expected)

        # overlapping = False
        expected = Counter(['GAT', 'TAC'])
        self.assertEqual(self.b1.k_word_counts(3, overlapping=False),
                         expected)
        expected = Counter(['GATTACA'])
        self.assertEqual(self.b1.k_word_counts(7, overlapping=False),
                         expected)

    def test_k_word_frequencies(self):
        """ k_word_frequencies functions as expected
        """
        # overlapping = True
        expected = defaultdict(int)
        expected['A'] = 3/7.
        expected['C'] = 1/7.
        expected['G'] = 1/7.
        expected['T'] = 2/7.
        self.assertEqual(self.b1.k_word_frequencies(1, overlapping=True),
                         expected)
        expected = defaultdict(int)
        expected['GAT'] = 1/5.
        expected['ATT'] = 1/5.
        expected['TTA'] = 1/5.
        expected['TAC'] = 1/5.
        expected['ACA'] = 1/5.
        self.assertEqual(self.b1.k_word_frequencies(3, overlapping=True),
                         expected)

        # overlapping = False
        expected = defaultdict(int)
        expected['GAT'] = 1/2.
        expected['TAC'] = 1/2.
        self.assertEqual(self.b1.k_word_frequencies(3, overlapping=False),
                         expected)
        expected = defaultdict(int)
        expected['GATTACA'] = 1.0
        self.assertEqual(self.b1.k_word_frequencies(7, overlapping=False),
                         expected)
        expected = defaultdict(int)
        empty = BiologicalSequence('')
        self.assertEqual(empty.k_word_frequencies(1, overlapping=False),
                         expected)

    def test_len(self):
        """ len functions as expected
        """
        self.assertEqual(len(self.b1), 7)
        self.assertEqual(len(self.b2), 9)
        self.assertEqual(len(self.b3), 4)

    def test_repr(self):
        """ repr functions as expected
        """
        self.assertEqual(repr(self.b1),
                         "<BiologicalSequence: GATTACA (length: 7)>")
        self.assertEqual(repr(self.b6),
                         "<BiologicalSequence: ACGTACGTAC... (length: 12)>")

    def test_reversed(self):
        """ reversed functions as expected
        """
        b1_reversed = reversed(self.b1)
        for actual, expected in zip(b1_reversed, "ACATTAG"):
            self.assertEqual(actual, expected)

        self.assertRaises(StopIteration, lambda: next(b1_reversed))

    def test_str(self):
        """ str functions as expected
        """
        self.assertEqual(str(self.b1), "GATTACA")
        self.assertEqual(str(self.b2), "ACCGGTACC")
        self.assertEqual(str(self.b3), "GREG")

    def test_alphabet(self):
        """ alphabet property functions as expected
        """
        self.assertEqual(self.b1.alphabet(), set())

    def test_description(self):
        """ description property functions as expected
        """
        self.assertEqual(self.b1.description, "")
        self.assertEqual(self.b2.description, "A test sequence")
        self.assertEqual(self.b3.description, "A protein sequence")

    def test_gap_alphabet(self):
        """ gap_alphabet property functions as expected
        """
        self.assertEqual(self.b1.gap_alphabet(), set('-.'))

    def test_identifier(self):
        """ identifier property functions as expected
        """
        self.assertEqual(self.b1.identifier, "")
        self.assertEqual(self.b2.identifier, "test-seq-2")
        self.assertEqual(self.b3.identifier, "test-seq-3")

    def test_count(self):
        """ count functions as expected
        """
        self.assertEqual(self.b1.count('A'), 3)
        self.assertEqual(self.b1.count('T'), 2)
        self.assertEqual(self.b1.count('TT'), 1)

    def test_degap(self):
        """ degap functions as expected
        """
        self.assertEqual(self.b1.degap(), self.b1)
        self.assertEqual(self.b7.degap(), BiologicalSequence(''))
        self.assertEqual(self.b8.degap(), BiologicalSequence('HELLO'))

    def test_distance(self):
        """ distance functions as expected
        """
        # note that test_hamming_distance covers default behavior more
        # extensively
        self.assertEqual(self.b1.distance(self.b1), 0.0)
        self.assertEqual(self.b1.distance(BiologicalSequence('GATTACC')), 1./7)

        def dumb_distance(x, y):
            return 42

        self.assertEqual(
            self.b1.distance(self.b1, distance_fn=dumb_distance), 42)

    def test_fraction_diff(self):
        """ fraction_diff functions as expected
        """
        self.assertEqual(self.b1.fraction_diff(self.b1), 0., 5)
        self.assertEqual(
            self.b1.fraction_diff(BiologicalSequence('GATTACC')), 1. / 7., 5)

    def test_fraction_same(self):
        """ fraction_same functions as expected
        """
        self.assertAlmostEqual(self.b1.fraction_same(self.b1), 1., 5)
        self.assertAlmostEqual(
            self.b1.fraction_same(BiologicalSequence('GATTACC')), 6. / 7., 5)

    def test_gap_maps(self):
        """ gap_maps functions as expected
        """
        # in sequence with no gaps, the gap_maps are identical
        self.assertEqual(self.b1.gap_maps(),
                         ([0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6]))
        # in sequence with all gaps, the map of degapped to gapped is the empty
        # list (bc its length is 0), and the map of gapped to degapped is all
        # None
        self.assertEqual(self.b7.gap_maps(),
                         ([], [None, None, None, None, None, None]))

        self.assertEqual(self.b8.gap_maps(),
                         ([0, 1, 8, 9, 10],
                          [0, 1, None, None, None, None, None, None, 2, 3, 4]))

        # example from the gap_maps doc string
        self.assertEqual(BiologicalSequence('-ACCGA-TA-').gap_maps(),
                         ([1, 2, 3, 4, 5, 7, 8],
                          [None, 0, 1, 2, 3, 4, None, 5, 6, None]))

    def test_gap_vector(self):
        """ gap_vector functions as expected
        """
        self.assertEqual(self.b1.gap_vector(),
                         [False] * len(self.b1))
        self.assertEqual(self.b7.gap_vector(),
                         [True] * len(self.b7))
        self.assertEqual(self.b8.gap_vector(),
                         [False, False, True, True, True, True,
                          True, True, False, False, False])

    def test_unsupported_characters(self):
        """ unsupported_characters functions as expected
        """
        self.assertEqual(self.b1.unsupported_characters(), set('GATC'))
        self.assertEqual(self.b7.unsupported_characters(), set())

    def test_has_unsupported_characters(self):
        """ has_unsupported_characters functions as expected
        """
        self.assertTrue(self.b1.has_unsupported_characters())
        self.assertFalse(self.b7.has_unsupported_characters())

    def test_index(self):
        """ index functions as expected """
        self.assertEqual(self.b1.index('G'), 0)
        self.assertEqual(self.b1.index('A'), 1)
        self.assertEqual(self.b1.index('AC'), 4)
        self.assertRaises(ValueError, self.b1.index, 'x')

    def test_is_gap(self):
        """ is_gap functions as expected """
        self.assertTrue(self.b1.is_gap('.'))
        self.assertTrue(self.b1.is_gap('-'))
        self.assertFalse(self.b1.is_gap('A'))
        self.assertFalse(self.b1.is_gap('x'))
        self.assertFalse(self.b1.is_gap(' '))
        self.assertFalse(self.b1.is_gap(''))

    def test_is_gapped(self):
        """ is_gapped functions as expected """
        self.assertFalse(self.b1.is_gapped())
        self.assertFalse(self.b2.is_gapped())
        self.assertTrue(self.b7.is_gapped())
        self.assertTrue(self.b8.is_gapped())

    def test_is_valid(self):
        """ is_valid functions as expected
        """
        self.assertFalse(self.b1.is_valid())
        self.assertTrue(self.b7.is_valid())

    def test_to_fasta(self):
        """ to_fasta functions as expected
        """
        self.assertEqual(self.b1.to_fasta(), ">\nGATTACA\n")
        self.assertEqual(self.b1.to_fasta(terminal_character=""), ">\nGATTACA")
        self.assertEqual(self.b2.to_fasta(),
                         ">test-seq-2 A test sequence\nACCGGTACC\n")
        self.assertEqual(self.b3.to_fasta(),
                         ">test-seq-3 A protein sequence\nGREG\n")
        self.assertEqual(self.b4.to_fasta(),
                         ">test-seq-4\nPRTEIN\n")
        self.assertEqual(self.b5.to_fasta(),
                         "> some description\nLLPRTEIN\n")

        # alt parameters
        self.assertEqual(self.b2.to_fasta(field_delimiter=":"),
                         ">test-seq-2:A test sequence\nACCGGTACC\n")
        self.assertEqual(self.b2.to_fasta(terminal_character="!"),
                         ">test-seq-2 A test sequence\nACCGGTACC!")
        self.assertEqual(
            self.b2.to_fasta(field_delimiter=":", terminal_character="!"),
            ">test-seq-2:A test sequence\nACCGGTACC!")

    def test_upper(self):
        """ upper functions as expected
        """
        b = NucleotideSequence('GAt.ACa-', identifier='x', description='42')
        expected = NucleotideSequence('GAT.ACA-', identifier='x',
                                      description='42')
        self.assertEqual(b.upper(), expected)

    def test_lower(self):
        """ lower functions as expected
        """
        b = NucleotideSequence('GAt.ACa-', identifier='x', description='42')
        expected = NucleotideSequence('gat.aca-', identifier='x',
                                      description='42')
        self.assertEqual(b.lower(), expected)


class NucelotideSequenceTests(TestCase):
    """ Tests of the BiologicalSequence class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
        self.empty = NucleotideSequence('')
        self.b1 = NucleotideSequence('GATTACA')
        self.b2 = NucleotideSequence(
            'ACCGGUACC', identifier="test-seq-2",
            description="A test sequence")

    def test_alphabet(self):
        """alphabet property functions as expected"""
        exp = {
            'A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'U', 'T',
            'W', 'V', 'Y', 'a', 'c', 'b', 'd', 'g', 'h', 'k', 'm', 'n', 's',
            'r', 'u', 't', 'w', 'v', 'y'
        }

        # Test calling from an instance and purely static context.
        self.assertEqual(self.b1.alphabet(), exp)
        self.assertEqual(NucleotideSequence.alphabet(), exp)

    def test_gap_alphabet(self):
        """gap_alphabet property functions as expected"""
        self.assertEqual(self.b1.gap_alphabet(), set('-.'))

    def test_complement_map(self):
        """complement_map property functions as expected"""
        exp = {}
        self.assertEqual(self.b1.complement_map(), exp)
        self.assertEqual(NucleotideSequence.complement_map(), exp)

    def test_iupac_standard_characters(self):
        """iupac_standard_characters property functions as expected"""
        exp = set("ACGTUacgtu")
        self.assertEqual(self.b1.iupac_standard_characters(), exp)
        self.assertEqual(NucleotideSequence.iupac_standard_characters(), exp)

    def test_iupac_degeneracies(self):
        """iupac_degeneracies property functions as expected"""
        exp = {
            # upper
            'B': set(['C', 'U', 'T', 'G']), 'D': set(['A', 'U', 'T', 'G']),
            'H': set(['A', 'C', 'U', 'T']), 'K': set(['U', 'T', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'U', 'T', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']),
            'W': set(['A', 'U', 'T']), 'V': set(['A', 'C', 'G']),
            'Y': set(['C', 'U', 'T']),
            # lower
            'b': set(['c', 'u', 't', 'g']), 'd': set(['a', 'u', 't', 'g']),
            'h': set(['a', 'c', 'u', 't']), 'k': set(['u', 't', 'g']),
            'm': set(['a', 'c']), 'n': set(['a', 'c', 'u', 't', 'g']),
            's': set(['c', 'g']), 'r': set(['a', 'g']),
            'w': set(['a', 'u', 't']), 'v': set(['a', 'c', 'g']),
            'y': set(['c', 'u', 't'])
        }
        self.assertEqual(self.b1.iupac_degeneracies(), exp)
        self.assertEqual(NucleotideSequence.iupac_degeneracies(), exp)

        # Test that we can modify a copy of the mapping without altering the
        # canonical representation.
        degen = NucleotideSequence.iupac_degeneracies()
        degen.update({'V': set("BRO"), 'Z': set("ZORRO")})
        self.assertNotEqual(degen, exp)
        self.assertEqual(NucleotideSequence.iupac_degeneracies(), exp)

    def test_iupac_degenerate_characters(self):
        """iupac_degenerate_characters property functions as expected"""
        exp = set(['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y',
                   'b', 'd', 'h', 'k', 'm', 'n', 's', 'r', 'w', 'v', 'y'])
        self.assertEqual(self.b1.iupac_degenerate_characters(), exp)
        self.assertEqual(NucleotideSequence.iupac_degenerate_characters(), exp)

    def test_iupac_characters(self):
        """iupac_characters property functions as expected"""
        exp = {
            'A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'U', 'T',
            'W', 'V', 'Y', 'a', 'c', 'b', 'd', 'g', 'h', 'k', 'm', 'n', 's',
            'r', 'u', 't', 'w', 'v', 'y'
        }

        self.assertEqual(self.b1.iupac_characters(), exp)
        self.assertEqual(NucleotideSequence.iupac_characters(), exp)

    def test_complement(self):
        """ complement fails (it's undefined for generic NucleotideSequence)
        """
        self.assertRaises(BiologicalSequenceError,
                          self.b1.complement)

    def test_reverse_complement(self):
        """ rev comp fails (it's undefined for generic NucleotideSequence)
        """
        self.assertRaises(BiologicalSequenceError,
                          self.b1.reverse_complement)

    def test_is_reverse_complement(self):
        """ is_reverse_complement fails (it's undefined)
        """
        self.assertRaises(BiologicalSequenceError,
                          self.b1.is_reverse_complement, self.b1)

    def test_nondegenerates_invalid(self):
        """Should raise error if seq has invalid characters."""
        with self.assertRaises(BiologicalSequenceError):
            _ = list(NucleotideSequence('AZA').nondegenerates())

    def test_nondegenerates_empty(self):
        """Should correctly handle an empty sequence."""
        self.assertEqual(list(self.empty.nondegenerates()), [self.empty])

    def test_nondegenerates_no_degens(self):
        """Should correctly handle a sequence without any degeneracies."""
        self.assertEqual(list(self.b1.nondegenerates()), [self.b1])

    def test_nondegenerates_all_degens(self):
        """Should correctly handle a purely-degenerate sequence."""
        # Same chars.
        exp = [NucleotideSequence('CC'), NucleotideSequence('CG'),
               NucleotideSequence('GC'), NucleotideSequence('GG')]
        # Sort based on sequence string, as order is not guaranteed.
        obs = sorted(NucleotideSequence('SS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Different chars.
        exp = [NucleotideSequence('AC'), NucleotideSequence('AG'),
               NucleotideSequence('GC'), NucleotideSequence('GG')]
        obs = sorted(NucleotideSequence('RS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Odd number of chars.
        obs = list(NucleotideSequence('NNN').nondegenerates())
        self.assertEqual(len(obs), 5**3)

    def test_nondegenerates_mixed_degens(self):
        """Should correctly handle a sequence with standard and degen chars."""
        exp = [NucleotideSequence('AGC'), NucleotideSequence('AGT'),
               NucleotideSequence('AGU'), NucleotideSequence('GGC'),
               NucleotideSequence('GGT'), NucleotideSequence('GGU')]
        obs = sorted(NucleotideSequence('RGY').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

    def test_nondegenerates_gap_mixed_case(self):
        """Should correctly handle a sequence with gap chars and mixed case."""
        exp = [NucleotideSequence('-A.a'), NucleotideSequence('-A.c'),
               NucleotideSequence('-C.a'), NucleotideSequence('-C.c')]
        obs = sorted(NucleotideSequence('-M.m').nondegenerates(), key=str)
        self.assertEqual(obs, exp)


class DNASequenceTests(TestCase):
    """ Tests of the DNASequence class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
        self.empty = DNASequence('')
        self.b1 = DNASequence('GATTACA')
        self.b2 = DNASequence(
            'ACCGGTACC', identifier="test-seq-2",
            description="A test sequence")
        self.b3 = DNASequence(
            'ACCGGUACC', identifier="bad-seq-1",
            description="Not a DNA sequence")
        self.b4 = DNASequence(
            'MRWSYKVHDBN', identifier="degen",
            description="All of the degenerate bases")
        self.b5 = DNASequence('.G--ATTAC-A...')

    def test_alphabet(self):
        """alphabet property functions as expected"""
        exp = {
            'A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'T', 'W',
            'V', 'Y', 'a', 'c', 'b', 'd', 'g', 'h', 'k', 'm', 'n', 's', 'r',
            't', 'w', 'v', 'y'
        }

        self.assertEqual(self.b1.alphabet(), exp)
        self.assertEqual(DNASequence.alphabet(), exp)

    def test_gap_alphabet(self):
        """gap_alphabet property functions as expected"""
        self.assertEqual(self.b1.gap_alphabet(), set('-.'))

    def test_complement_map(self):
        """complement_map property functions as expected"""
        exp = {
            '-': '-', '.': '.', 'A': 'T', 'C': 'G', 'B': 'V', 'D': 'H',
            'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', 'S': 'S',
            'R': 'Y', 'T': 'A', 'W': 'W', 'V': 'B', 'Y': 'R', 'a': 't',
            'c': 'g', 'b': 'v', 'd': 'h', 'g': 'c', 'h': 'd', 'k': 'm',
            'm': 'k', 'n': 'n', 's': 's', 'r': 'y', 't': 'a', 'w': 'w',
            'v': 'b', 'y': 'r'
        }
        self.assertEqual(self.b1.complement_map(), exp)
        self.assertEqual(DNASequence.complement_map(), exp)

    def test_iupac_standard_characters(self):
        """iupac_standard_characters property functions as expected"""
        exp = set("ACGTacgt")
        self.assertEqual(self.b1.iupac_standard_characters(), exp)
        self.assertEqual(DNASequence.iupac_standard_characters(), exp)

    def test_iupac_degeneracies(self):
        """iupac_degeneracies property functions as expected"""
        exp = {
            'B': set(['C', 'T', 'G']), 'D': set(['A', 'T', 'G']),
            'H': set(['A', 'C', 'T']), 'K': set(['T', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'T', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']), 'W': set(['A', 'T']),
            'V': set(['A', 'C', 'G']), 'Y': set(['C', 'T']),
            'b': set(['c', 't', 'g']), 'd': set(['a', 't', 'g']),
            'h': set(['a', 'c', 't']), 'k': set(['t', 'g']),
            'm': set(['a', 'c']), 'n': set(['a', 'c', 't', 'g']),
            's': set(['c', 'g']), 'r': set(['a', 'g']), 'w': set(['a', 't']),
            'v': set(['a', 'c', 'g']), 'y': set(['c', 't'])
        }
        self.assertEqual(self.b1.iupac_degeneracies(), exp)
        self.assertEqual(DNASequence.iupac_degeneracies(), exp)

    def test_iupac_degenerate_characters(self):
        """iupac_degenerate_characters property functions as expected"""
        exp = set(['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y',
                   'b', 'd', 'h', 'k', 'm', 'n', 's', 'r', 'w', 'v', 'y'])
        self.assertEqual(self.b1.iupac_degenerate_characters(), exp)
        self.assertEqual(DNASequence.iupac_degenerate_characters(), exp)

    def test_iupac_characters(self):
        """iupac_characters property functions as expected"""
        exp = {
            'A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'T', 'W',
            'V', 'Y', 'a', 'c', 'b', 'd', 'g', 'h', 'k', 'm', 'n', 's', 'r',
            't', 'w', 'v', 'y'
        }
        self.assertEqual(self.b1.iupac_characters(), exp)
        self.assertEqual(DNASequence.iupac_characters(), exp)

    def test_complement(self):
        """ complement functions as expected
        """
        self.assertEqual(self.b1.complement(), DNASequence("CTAATGT"))
        self.assertEqual(self.b2.complement(), DNASequence("TGGCCATGG"))
        self.assertRaises(BiologicalSequenceError, self.b3.complement)
        self.assertEqual(self.b4.complement(), DNASequence("KYWSRMBDHVN"))
        self.assertEqual(self.b5.complement(), DNASequence(".C--TAATG-T..."))

    def test_reverse_complement(self):
        """ reverse complement functions as expected
        """
        self.assertEqual(self.b1.reverse_complement(), DNASequence("TGTAATC"))
        self.assertEqual(self.b2.reverse_complement(),
                         DNASequence("GGTACCGGT"))
        self.assertRaises(BiologicalSequenceError, self.b3.reverse_complement)
        self.assertEqual(self.b4.reverse_complement(),
                         DNASequence("NVHDBMRSWYK"))

    def test_unsupported_characters(self):
        """ unsupported_characters functions as expected
        """
        self.assertEqual(self.b1.unsupported_characters(), set())
        self.assertEqual(self.b2.unsupported_characters(), set())
        self.assertEqual(self.b3.unsupported_characters(), set('U'))
        self.assertEqual(self.b4.unsupported_characters(), set())

    def test_has_unsupported_characters(self):
        """ has_unsupported_characters functions as expected
        """
        self.assertFalse(self.b1.has_unsupported_characters())
        self.assertFalse(self.b2.has_unsupported_characters())
        self.assertTrue(self.b3.has_unsupported_characters())
        self.assertFalse(self.b4.has_unsupported_characters())

    def test_is_reverse_complement(self):
        """ is_reverse_complement functions as expected
        """
        self.assertFalse(self.b1.is_reverse_complement(self.b1))
        self.assertTrue(
            self.b1.is_reverse_complement(DNASequence('TGTAATC')))
        self.assertTrue(
            self.b4.is_reverse_complement(DNASequence('NVHDBMRSWYK')))

    def test_nondegenerates_invalid(self):
        """Should raise error if seq has invalid characters."""
        with self.assertRaises(BiologicalSequenceError):
            _ = list(DNASequence('AZA').nondegenerates())

    def test_nondegenerates_empty(self):
        """Should correctly handle an empty sequence."""
        self.assertEqual(list(self.empty.nondegenerates()), [self.empty])

    def test_nondegenerates_no_degens(self):
        """Should correctly handle a sequence without any degeneracies."""
        self.assertEqual(list(self.b1.nondegenerates()), [self.b1])

    def test_nondegenerates_all_degens(self):
        """Should correctly handle a purely-degenerate sequence."""
        # Same chars.
        exp = [DNASequence('CC'), DNASequence('CG'), DNASequence('GC'),
               DNASequence('GG')]
        # Sort based on sequence string, as order is not guaranteed.
        obs = sorted(DNASequence('SS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Different chars.
        exp = [DNASequence('AC'), DNASequence('AG'), DNASequence('GC'),
               DNASequence('GG')]
        obs = sorted(DNASequence('RS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Odd number of chars.
        obs = list(DNASequence('NNN').nondegenerates())
        self.assertEqual(len(obs), 4**3)

    def test_nondegenerates_mixed_degens(self):
        """Should correctly handle a sequence with standard and degen chars."""
        exp = [DNASequence('AGC'), DNASequence('AGT'), DNASequence('GGC'),
               DNASequence('GGT')]
        obs = sorted(DNASequence('RGY').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

    def test_nondegenerates_gap_mixed_case(self):
        """Should correctly handle a sequence with gap chars and mixed case."""
        exp = [DNASequence('-A.a'), DNASequence('-A.c'),
               DNASequence('-C.a'), DNASequence('-C.c')]
        obs = sorted(DNASequence('-M.m').nondegenerates(), key=str)
        self.assertEqual(obs, exp)


class RNASequenceTests(TestCase):
    """ Tests of the RNASequence class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
        self.empty = RNASequence('')
        self.b1 = RNASequence('GAUUACA')
        self.b2 = RNASequence(
            'ACCGGUACC', identifier="test-seq-2",
            description="A test sequence")
        self.b3 = RNASequence(
            'ACCGGTACC', identifier="bad-seq-1",
            description="Not a RNA sequence")
        self.b4 = RNASequence(
            'MRWSYKVHDBN', identifier="degen",
            description="All of the degenerate bases")
        self.b5 = RNASequence('.G--AUUAC-A...')

    def test_alphabet(self):
        """alphabet property functions as expected"""
        exp = {
            'A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'U', 'W',
            'V', 'Y', 'a', 'c', 'b', 'd', 'g', 'h', 'k', 'm', 'n', 's', 'r',
            'u', 'w', 'v', 'y'
        }

        self.assertEqual(self.b1.alphabet(), exp)
        self.assertEqual(RNASequence.alphabet(), exp)

    def test_gap_alphabet(self):
        """gap_alphabet property functions as expected"""
        self.assertEqual(self.b1.gap_alphabet(), set('-.'))

    def test_complement_map(self):
        """complement_map property functions as expected"""
        exp = {
            '-': '-', '.': '.', 'A': 'U', 'C': 'G', 'B': 'V', 'D': 'H',
            'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', 'S': 'S',
            'R': 'Y', 'U': 'A', 'W': 'W', 'V': 'B', 'Y': 'R', 'a': 'u',
            'c': 'g', 'b': 'v', 'd': 'h', 'g': 'c', 'h': 'd', 'k': 'm',
            'm': 'k', 'n': 'n', 's': 's', 'r': 'y', 'u': 'a', 'w': 'w',
            'v': 'b', 'y': 'r'
        }
        self.assertEqual(self.b1.complement_map(), exp)
        self.assertEqual(RNASequence.complement_map(), exp)

    def test_iupac_standard_characters(self):
        """iupac_standard_characters property functions as expected"""
        exp = set("ACGUacgu")
        self.assertEqual(self.b1.iupac_standard_characters(), exp)
        self.assertEqual(RNASequence.iupac_standard_characters(), exp)

    def test_iupac_degeneracies(self):
        """iupac_degeneracies property functions as expected"""
        exp = {
            'B': set(['C', 'U', 'G']), 'D': set(['A', 'U', 'G']),
            'H': set(['A', 'C', 'U']), 'K': set(['U', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'U', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']), 'W': set(['A', 'U']),
            'V': set(['A', 'C', 'G']), 'Y': set(['C', 'U']),
            'b': set(['c', 'u', 'g']), 'd': set(['a', 'u', 'g']),
            'h': set(['a', 'c', 'u']), 'k': set(['u', 'g']),
            'm': set(['a', 'c']), 'n': set(['a', 'c', 'u', 'g']),
            's': set(['c', 'g']), 'r': set(['a', 'g']), 'w': set(['a', 'u']),
            'v': set(['a', 'c', 'g']), 'y': set(['c', 'u'])
        }
        self.assertEqual(self.b1.iupac_degeneracies(), exp)
        self.assertEqual(RNASequence.iupac_degeneracies(), exp)

    def test_iupac_degenerate_characters(self):
        """iupac_degenerate_characters property functions as expected"""
        exp = set(['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y',
                   'b', 'd', 'h', 'k', 'm', 'n', 's', 'r', 'w', 'v', 'y'])
        self.assertEqual(self.b1.iupac_degenerate_characters(), exp)
        self.assertEqual(RNASequence.iupac_degenerate_characters(), exp)

    def test_iupac_characters(self):
        """iupac_characters property functions as expected"""
        exp = {
            'A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'U', 'W',
            'V', 'Y', 'a', 'c', 'b', 'd', 'g', 'h', 'k', 'm', 'n', 's', 'r',
            'u', 'w', 'v', 'y'
        }
        self.assertEqual(self.b1.iupac_characters(), exp)
        self.assertEqual(RNASequence.iupac_characters(), exp)

    def test_complement(self):
        """ complement functions as expected
        """
        self.assertEqual(self.b1.complement(), RNASequence("CUAAUGU"))
        self.assertEqual(self.b2.complement(), RNASequence("UGGCCAUGG"))
        self.assertRaises(BiologicalSequenceError, self.b3.complement)
        self.assertEqual(self.b4.complement(), RNASequence("KYWSRMBDHVN"))
        self.assertEqual(self.b5.complement(), RNASequence(".C--UAAUG-U..."))

    def test_reverse_complement(self):
        """ reverse complement functions as expected
        """
        self.assertEqual(self.b1.reverse_complement(), RNASequence("UGUAAUC"))
        self.assertEqual(self.b2.reverse_complement(),
                         RNASequence("GGUACCGGU"))
        self.assertRaises(BiologicalSequenceError, self.b3.reverse_complement)
        self.assertEqual(self.b4.reverse_complement(),
                         RNASequence("NVHDBMRSWYK"))

    def test_unsupported_characters(self):
        """ unsupported_characters functions as expected
        """
        self.assertEqual(self.b1.unsupported_characters(), set())
        self.assertEqual(self.b2.unsupported_characters(), set())
        self.assertEqual(self.b3.unsupported_characters(), set('T'))
        self.assertEqual(self.b4.unsupported_characters(), set())

    def test_has_unsupported_characters(self):
        """ has_unsupported_characters functions as expected
        """
        self.assertFalse(self.b1.has_unsupported_characters())
        self.assertFalse(self.b2.has_unsupported_characters())
        self.assertTrue(self.b3.has_unsupported_characters())
        self.assertFalse(self.b4.has_unsupported_characters())

    def test_is_reverse_complement(self):
        """ is_reverse_complement functions as expected
        """
        self.assertFalse(self.b1.is_reverse_complement(self.b1))
        self.assertTrue(
            self.b1.is_reverse_complement(RNASequence('UGUAAUC')))
        self.assertTrue(
            self.b4.is_reverse_complement(RNASequence('NVHDBMRSWYK')))

    def test_nondegenerates_invalid(self):
        """Should raise error if seq has invalid characters."""
        with self.assertRaises(BiologicalSequenceError):
            _ = list(RNASequence('AZA').nondegenerates())

    def test_nondegenerates_empty(self):
        """Should correctly handle an empty sequence."""
        self.assertEqual(list(self.empty.nondegenerates()), [self.empty])

    def test_nondegenerates_no_degens(self):
        """Should correctly handle a sequence without any degeneracies."""
        self.assertEqual(list(self.b1.nondegenerates()), [self.b1])

    def test_nondegenerates_all_degens(self):
        """Should correctly handle a purely-degenerate sequence."""
        # Same chars.
        exp = [RNASequence('CC'), RNASequence('CG'), RNASequence('GC'),
               RNASequence('GG')]
        # Sort based on sequence string, as order is not guaranteed.
        obs = sorted(RNASequence('SS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Different chars.
        exp = [RNASequence('AC'), RNASequence('AG'), RNASequence('GC'),
               RNASequence('GG')]
        obs = sorted(RNASequence('RS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Odd number of chars.
        obs = list(RNASequence('NNN').nondegenerates())
        self.assertEqual(len(obs), 4**3)

    def test_nondegenerates_mixed_degens(self):
        """Should correctly handle a sequence with standard and degen chars."""
        exp = [RNASequence('AGC'), RNASequence('AGU'), RNASequence('GGC'),
               RNASequence('GGU')]
        obs = sorted(RNASequence('RGY').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

    def test_nondegenerates_gap_mixed_case(self):
        """Should correctly handle a sequence with gap chars and mixed case."""
        exp = [RNASequence('-A.a'), RNASequence('-A.c'),
               RNASequence('-C.a'), RNASequence('-C.c')]
        obs = sorted(RNASequence('-M.m').nondegenerates(), key=str)
        self.assertEqual(obs, exp)


if __name__ == "__main__":
    main()
