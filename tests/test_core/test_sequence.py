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

    def setup(self):
        """ Initialize values to be used in tests
        """
        pass

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


if __name__ == "__main__":
    main()
