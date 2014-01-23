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

    def test_getitem(self):
        """ getitem functions as expected
        """
        self.assertEqual(self.b1[0],'G')
        self.assertEqual(self.b1[:],'GATTACA')
        self.assertEqual(self.b1[::-1],'ACATTAG')

    def test_len(self):
        """ len functions as expected
        """
        self.assertEqual(len(self.b1),7)
        self.assertEqual(len(self.b2),9)
        self.assertEqual(len(self.b3),4)
        

if __name__ == "__main__":
    main()
