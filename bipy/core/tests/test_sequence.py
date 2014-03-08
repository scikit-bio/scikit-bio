#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from unittest import TestCase, main

from bipy.core.sequence import (
    BiologicalSequence, NucleotideSequence, DNASequence, RNASequence,
    DNA, RNA)
from bipy.core.exception import BiologicalSequenceError


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
        self.assertEqual(self.b1[0], 'G')
        self.assertEqual(self.b1[:], 'GATTACA')
        self.assertEqual(self.b1[::-1], 'ACATTAG')

    def test_iter(self):
        """ iter functions as expected
        """
        b1_iter = iter(self.b1)
        for actual, expected in zip(b1_iter, "GATTACA"):
            self.assertEqual(actual, expected)

        self.assertRaises(StopIteration, b1_iter.next)

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

        self.assertRaises(StopIteration, b1_reversed.next)

    def test_str(self):
        """ str functions as expected
        """
        self.assertEqual(str(self.b1), "GATTACA")
        self.assertEqual(str(self.b2), "ACCGGTACC")
        self.assertEqual(str(self.b3), "GREG")

    def test_alphabet(self):
        """ alphabet property functions as expected
        """
        self.assertEqual(self.b1.alphabet, set())

    def test_description(self):
        """ description property functions as expected
        """
        self.assertEqual(self.b1.description, "")
        self.assertEqual(self.b2.description, "A test sequence")
        self.assertEqual(self.b3.description, "A protein sequence")

    def test_gap_alphabet(self):
        """ gap_alphabet property functions as expected
        """
        self.assertEqual(self.b1.gap_alphabet, set('-.'))

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


class NucelotideSequenceTests(TestCase):
    """ Tests of the BiologicalSequence class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
        self.b1 = NucleotideSequence('GATTACA')
        self.b2 = NucleotideSequence(
            'ACCGGUACC', identifier="test-seq-2",
            description="A test sequence")

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


class DNASequenceTests(TestCase):
    """ Tests of the DNASequence class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
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

    def test_complement(self):
        """ complement functions as expected
        """
        self.assertEqual(self.b1.complement(), DNASequence("CTAATGT"))
        self.assertEqual(self.b2.complement(), DNASequence("TGGCCATGG"))
        self.assertRaises(BiologicalSequenceError, self.b3.complement)
        self.assertEqual(self.b4.complement(), DNASequence("KYWSRMBDHVN"))

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


class RNASequenceTests(TestCase):
    """ Tests of the DNASequence class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
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

    def test_complement(self):
        """ complement functions as expected
        """
        self.assertEqual(self.b1.complement(), RNASequence("CUAAUGU"))
        self.assertEqual(self.b2.complement(), RNASequence("UGGCCAUGG"))
        self.assertRaises(BiologicalSequenceError, self.b3.complement)
        self.assertEqual(self.b4.complement(), RNASequence("KYWSRMBDHVN"))

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

if __name__ == "__main__":
    main()
