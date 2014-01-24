#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from unittest import TestCase, main

from bipy.core.sequence import (
    BiologicalSequence, NucleotideSequence, DNASequence, RNASequence,
    DnaSequence, RnaSequence, BiologicalSequenceError)

class BiologicalSequenceTests(TestCase):
    """ Tests of the BiologicalSequence class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
        self.b1 = BiologicalSequence('GATTACA')
        self.b2 = BiologicalSequence(
         'ACCGGTACC', identifier="test-seq-2", description="A test sequence")
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
        self.assertEqual(str(b),'ACCGGXZY') 
        self.assertEqual(b.Identifier,"")
        self.assertEqual(b.Description,"")
        
        # init as string with optional values
        b = BiologicalSequence(
         'ACCGGXZY','test-seq-1','The first test sequence')
        self.assertEqual(str(b),'ACCGGXZY') 
        self.assertEqual(b.Identifier,"test-seq-1")
        self.assertEqual(b.Description,"The first test sequence")

        # test init as a different string
        b = BiologicalSequence('WRRTY')
        self.assertEqual(str(b),'WRRTY') 

        # init as list
        b = BiologicalSequence(list('ACCGGXZY'))
        self.assertEqual(str(b),'ACCGGXZY') 
        self.assertEqual(b.Identifier,"")
        self.assertEqual(b.Description,"")
        
        # init as tuple
        b = BiologicalSequence(tuple('ACCGGXZY'))
        self.assertEqual(str(b),'ACCGGXZY') 
        self.assertEqual(b.Identifier,"")
        self.assertEqual(b.Description,"")

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
         BiologicalSequence('ACGT',identifier='a') == 
         BiologicalSequence('ACGT',identifier='b'))
        self.assertTrue(
         BiologicalSequence('ACGT',description='c') == 
         BiologicalSequence('ACGT',description='d'))
        self.assertTrue(
         BiologicalSequence('ACGT',identifier='a',description='c') == 
         BiologicalSequence('ACGT',identifier='b',description='d'))
        
        # different type causes sequences to not be equal
        self.assertFalse(
         BiologicalSequence('ACGT') == NucleotideSequence('ACGT'))

    def test_getitem(self):
        """ getitem functions as expected
        """
        self.assertEqual(self.b1[0],'G')
        self.assertEqual(self.b1[:],'GATTACA')
        self.assertEqual(self.b1[::-1],'ACATTAG')

    def test_iter(self):
        """ iter functions as expected
        """
        b1_iter = iter(self.b1)
        for actual, expected in zip(b1_iter, "GATTACA"):
            self.assertEqual(actual,expected)

        self.assertRaises(StopIteration,b1_iter.next)

    def test_len(self):
        """ len functions as expected
        """
        self.assertEqual(len(self.b1),7)
        self.assertEqual(len(self.b2),9)
        self.assertEqual(len(self.b3),4)

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
            self.assertEqual(actual,expected)

        self.assertRaises(StopIteration,b1_reversed.next)

    def test_str(self):
        """ str functions as expected
        """
        self.assertEqual(str(self.b1),"GATTACA")
        self.assertEqual(str(self.b2),"ACCGGTACC")
        self.assertEqual(str(self.b3),"GREG")

    def test_minimum_edit_distance(self):
        """ _minimum_edit_distance functions as expected
        """
        s1 = BiologicalSequence('AAAAA')
        s2 = BiologicalSequence('AAAAA')
        s3 = BiologicalSequence('AAAAAA')
        s4 = BiologicalSequence('TAAAAA')
        s5 = BiologicalSequence('CCCCCCX')
        s6 = BiologicalSequence('CCCCCCC')
        
        self.assertEqual(s1._minimum_edit_distance(s2),0)
        self.assertEqual(s1._minimum_edit_distance(s3),0)
        self.assertEqual(s1._minimum_edit_distance(s4),1)
        self.assertEqual(s1._minimum_edit_distance(s5),5)
        self.assertEqual(s4._minimum_edit_distance(s5),6)
        self.assertEqual(s4._minimum_edit_distance(s6),6)

    def test_Alphabet(self):
        """ Alphabet property functions as expected
        """
        self.assertEqual(self.b1.Alphabet,set())

    def test_Description(self):
        """ Description property functions as expected
        """
        self.assertEqual(self.b1.Description,"")
        self.assertEqual(self.b2.Description,"A test sequence")
        self.assertEqual(self.b3.Description,"A protein sequence")

    def test_GapAlphabet(self):
        """ GapAlphabet property functions as expected
        """
        self.assertEqual(self.b1.GapAlphabet,set('-.'))

    def test_Identifier(self):
        """ Identifier property functions as expected
        """
        self.assertEqual(self.b1.Identifier,"")
        self.assertEqual(self.b2.Identifier,"test-seq-2")
        self.assertEqual(self.b3.Identifier,"test-seq-3")

    def test_count(self):
        """ count functions as expected
        """
        self.assertEqual(self.b1.count('A'),3)
        self.assertEqual(self.b1.count('T'),2)

    def test_degap(self):
        """ degap functions as expected
        """
        self.assertEqual(self.b1.degap(),self.b1)
        self.assertEqual(self.b7.degap(),BiologicalSequence(''))
        self.assertEqual(self.b8.degap(),BiologicalSequence('HELLO'))

    def test_distance(self):
        """ distance functions as expected
        """
        # note that test_minimum_edit_distance covers default behavior more
        # extensively
        self.assertEqual(self.b1.distance(self.b1),0)
        self.assertEqual(self.b1.distance(BiologicalSequence('GATTACC')),1)

        def dumb_distance(x,y):
            return 42

        self.assertEqual(
         self.b1.distance(self.b1,distance_fn=dumb_distance),42)

    def test_fractionDiff(self):
        """ fractionDiff functions as expected
        """
        self.assertEqual(self.b1.fractionDiff(self.b1),0.,5)
        self.assertEqual(
         self.b1.fractionDiff(BiologicalSequence('GATTACC')),1./7.,5)

    def test_fractionSame(self):
        """ fractionDiff functions as expected
        """
        self.assertAlmostEqual(self.b1.fractionSame(self.b1),1.,5)
        self.assertAlmostEqual(
         self.b1.fractionSame(BiologicalSequence('GATTACC')),6./7.,5)

    def test_gapVector(self):
        """ gapVector functions as expected
        """
        self.assertEqual(self.b1.gapVector(),
                         [False] * len(self.b1))
        self.assertEqual(self.b7.gapVector(),
                         [True] * len(self.b7))
        self.assertEqual(self.b8.gapVector(),
                        [False,False,True,True,True,True,
                         True,True,False,False,False])

    def test_toFasta(self):
        """ toFasta functions as expected
        """
        self.assertEqual(self.b1.toFasta(),">\nGATTACA\n")
        self.assertEqual(self.b1.toFasta(terminal_character=""),">\nGATTACA")
        self.assertEqual(self.b2.toFasta(),
                         ">test-seq-2 A test sequence\nACCGGTACC\n")
        self.assertEqual(self.b3.toFasta(),
                         ">test-seq-3 A protein sequence\nGREG\n")
        self.assertEqual(self.b4.toFasta(),
                         ">test-seq-4\nPRTEIN\n")
        self.assertEqual(self.b5.toFasta(),
                         "> some description\nLLPRTEIN\n")

        # alt parameters
        self.assertEqual(self.b2.toFasta(field_delimiter=":"),
                        ">test-seq-2:A test sequence\nACCGGTACC\n")
        self.assertEqual(self.b2.toFasta(terminal_character="!"),
                        ">test-seq-2 A test sequence\nACCGGTACC!")
        self.assertEqual(
         self.b2.toFasta(field_delimiter=":",terminal_character="!"),
         ">test-seq-2:A test sequence\nACCGGTACC!")

    def test_getUnsupportedCharacters(self):
        """ getUnsupportedCharacters functions as expected
        """
        self.assertEqual(self.b1.getUnsupportedCharacters(),set('GATC'))
        self.assertEqual(self.b7.getUnsupportedCharacters(),set())

    def test_hasUnsupportedCharacters(self):
        """ hasUnsupportedCharacters functions as expected
        """
        self.assertTrue(self.b1.hasUnsupportedCharacters())
        self.assertFalse(self.b7.hasUnsupportedCharacters())

class NucelotideSequenceTests(TestCase):
    """ Tests of the BiologicalSequence class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
        self.b1 = NucleotideSequence('GATTACA')
        self.b2 = NucleotideSequence(
         'ACCGGUACC', identifier="test-seq-2", description="A test sequence")

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

class DNASequenceTests(TestCase):
    """ Tests of the DNASequence class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
        self.b1 = DNASequence('GATTACA')
        self.b2 = DNASequence(
         'ACCGGTACC', identifier="test-seq-2", description="A test sequence")
        self.b3 = DNASequence(
         'ACCGGUACC', identifier="bad-seq-1", description="Not a DNA sequence")
        self.b4 = DNASequence(
         'MRWSYKVHDBN', identifier="degen", 
         description="All of the degenerate bases")

    def test_complement(self):
        """ complement functions as expected
        """
        self.assertEqual(self.b1.complement(),DNASequence("CTAATGT"))
        self.assertEqual(self.b2.complement(),DNASequence("TGGCCATGG"))
        self.assertRaises(BiologicalSequenceError,self.b3.complement)
        self.assertEqual(self.b4.complement(),DNASequence("KYWSRMBDHVN"))

    def test_reverse_complement(self):
        """ reverse complement functions as expected
        """
        self.assertEqual(self.b1.reverse_complement(),DNASequence("TGTAATC"))
        self.assertEqual(self.b2.reverse_complement(),DNASequence("GGTACCGGT"))
        self.assertRaises(BiologicalSequenceError,self.b3.reverse_complement)
        self.assertEqual(self.b4.reverse_complement(),
                         DNASequence("NVHDBMRSWYK"))

    def test_getUnsupportedCharacters(self):
        """ getUnsupportedCharacters functions as expected
        """
        self.assertEqual(self.b1.getUnsupportedCharacters(),set())
        self.assertEqual(self.b2.getUnsupportedCharacters(),set())
        self.assertEqual(self.b3.getUnsupportedCharacters(),set('U'))
        self.assertEqual(self.b4.getUnsupportedCharacters(),set())

    def test_hasUnsupportedCharacters(self):
        """ hasUnsupportedCharacters functions as expected
        """
        self.assertFalse(self.b1.hasUnsupportedCharacters())
        self.assertFalse(self.b2.hasUnsupportedCharacters())
        self.assertTrue(self.b3.hasUnsupportedCharacters())
        self.assertFalse(self.b4.hasUnsupportedCharacters())


class RNASequenceTests(TestCase):
    """ Tests of the DNASequence class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
        self.b1 = RNASequence('GAUUACA')
        self.b2 = RNASequence(
         'ACCGGUACC', identifier="test-seq-2", description="A test sequence")
        self.b3 = RNASequence(
         'ACCGGTACC', identifier="bad-seq-1", description="Not a RNA sequence")
        self.b4 = RNASequence(
         'MRWSYKVHDBN', identifier="degen", 
         description="All of the degenerate bases")

    def test_complement(self):
        """ complement functions as expected
        """
        self.assertEqual(self.b1.complement(),RNASequence("CUAAUGU"))
        self.assertEqual(self.b2.complement(),RNASequence("UGGCCAUGG"))
        self.assertRaises(BiologicalSequenceError,self.b3.complement)
        self.assertEqual(self.b4.complement(),RNASequence("KYWSRMBDHVN"))

    def test_reverse_complement(self):
        """ reverse complement functions as expected
        """
        self.assertEqual(self.b1.reverse_complement(),RNASequence("UGUAAUC"))
        self.assertEqual(self.b2.reverse_complement(),RNASequence("GGUACCGGU"))
        self.assertRaises(BiologicalSequenceError,self.b3.reverse_complement)
        self.assertEqual(self.b4.reverse_complement(),
                         RNASequence("NVHDBMRSWYK"))

    def test_getUnsupportedCharacters(self):
        """ getUnsupportedCharacters functions as expected
        """
        self.assertEqual(self.b1.getUnsupportedCharacters(),set())
        self.assertEqual(self.b2.getUnsupportedCharacters(),set())
        self.assertEqual(self.b3.getUnsupportedCharacters(),set('T'))
        self.assertEqual(self.b4.getUnsupportedCharacters(),set())

    def test_hasUnsupportedCharacters(self):
        """ hasUnsupportedCharacters functions as expected
        """
        self.assertFalse(self.b1.hasUnsupportedCharacters())
        self.assertFalse(self.b2.hasUnsupportedCharacters())
        self.assertTrue(self.b3.hasUnsupportedCharacters())
        self.assertFalse(self.b4.hasUnsupportedCharacters())


if __name__ == "__main__":
    main()
