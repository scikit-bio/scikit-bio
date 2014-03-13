#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from unittest import TestCase, main

from bipy.core.sequence import (BiologicalSequence, NucleotideSequence, 
        DNASequence, RNASequence)
from bipy.core.alignment import (SequenceCollection, Alignment)
from bipy.core.exception import SequenceCollectionError


class SequenceCollectionTests(TestCase):
    """ Tests of the SequenceCollection class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
        s1_seqs = [('s1', 'GATTACA'), ('s2', 'TTG')]
        self.s1 = SequenceCollection(s1_seqs, DNASequence)

    def test_init(self):
        """ Initialization functions as expected with varied input types
        """
        seqs = [('s1', 'GATTACA'), ('s2', 'TTG')]
        SequenceCollection(seqs, BiologicalSequence)
        SequenceCollection(seqs, NucleotideSequence)
        SequenceCollection(seqs, DNASequence)
        
        seqs = [('s1', 'GAUUACA'), ('s2', 'UUG')]
        SequenceCollection(seqs, BiologicalSequence)
        SequenceCollection(seqs, NucleotideSequence)
        SequenceCollection(seqs, RNASequence)

    def test_init_validate(self):
        """ initialization with validation functions as expected
        """
        seqs = [('s1', 'GATTACA'), ('s2', 'TTG')]
        SequenceCollection(seqs, NucleotideSequence, validate=True)
        SequenceCollection(seqs, DNASequence, validate=True)
        self.assertRaises(SequenceCollectionError, SequenceCollection, seqs,
                BiologicalSequence, validate=True)
        
        seqs = [('s1', 'GXTTACA'), ('s2', 'TTG')]
        self.assertRaises(SequenceCollectionError, SequenceCollection, seqs,
                DNASequence, validate=True)
        self.assertRaises(SequenceCollectionError, SequenceCollection, seqs,
                DNASequence, validate=True)
        
        seqs = [('s1', 'GATTACA'), ('s2', 'TXG')]
        self.assertRaises(SequenceCollectionError, SequenceCollection, seqs,
                DNASequence, validate=True)
        self.assertRaises(SequenceCollectionError, SequenceCollection, seqs,
                DNASequence, validate=True)

    def test_getitem(self):
        """ getitem functions as expected
        """
        raise NotImplementedError

    def test_iter(self):
        """ iter functions as expected
        """
        raise NotImplementedError

        self.assertRaises(StopIteration, b1_iter.next)

    def test_len(self):
        """ len functions as expected
        """
        raise NotImplementedError

    def test_is_valid(self):
        """ is_valid functions as expected
        """
        raise NotImplementedError

    def test_to_fasta(self):
        """ to_fasta functions as expected
        """
        raise NotImplementedError

if __name__ == "__main__":
    main()
