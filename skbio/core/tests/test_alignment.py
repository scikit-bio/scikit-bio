#!/usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from unittest import TestCase, main
from collections import Counter, defaultdict

import numpy as np

from skbio.core.sequence import (BiologicalSequence, NucleotideSequence, 
        DNASequence, RNASequence)
from skbio.core.alignment import (SequenceCollection, Alignment)
from skbio.core.exception import SequenceCollectionError
from skbio.core.distance import SymmetricDistanceMatrix

class SequenceCollectionTests(TestCase):
    """ Tests of the SequenceCollection class """

    def setUp(self):
        """ Initialize values to be used in tests
        """
        self.d1 = DNASequence('GATTACA', identifier="d1")
        self.d2 = DNASequence('TTG', identifier="d2")
        self.d1_lower = DNASequence('gattaca', identifier="d1")
        self.d2_lower = DNASequence('ttg', identifier="d2")
        self.r1 = RNASequence('GAUUACA', identifier="r1")
        self.r2 = RNASequence('UUG', identifier="r2")
        self.r3 = RNASequence('U-----UGCC--', identifier="r3")
        
        self.i1 = DNASequence('GATXACA', identifier="i1")

        self.seqs1 = [self.d1, self.d2]
        self.seqs1_lower = [self.d1_lower, self.d2_lower]
        self.seqs2 = [self.r1, self.r2, self.r3]
        self.seqs3 = self.seqs1 + self.seqs2

        self.seqs1_t = [('d1', 'GATTACA'), ('d2', 'TTG')]
        self.seqs2_t = [('r1', 'GAUUACA'), ('r2', 'UUG'),
                ('r3', 'U-----UGCC--')]
        self.seqs3_t = self.seqs1_t + self.seqs2_t

        self.s1 = SequenceCollection(self.seqs1)
        self.s1_lower = SequenceCollection(self.seqs1_lower)
        self.s2 = SequenceCollection(self.seqs2)
        self.s3 = SequenceCollection(self.seqs3)

        self.invalid_s1 = SequenceCollection([self.i1])

    def test_init(self):
        """ Initialization functions as expected with varied input types
        """
        SequenceCollection(self.seqs1)
        SequenceCollection(self.seqs1)
        SequenceCollection(self.seqs1)
        
        SequenceCollection(self.seqs2)
        SequenceCollection(self.seqs2)
        SequenceCollection(self.seqs2)
        
        SequenceCollection(self.seqs3)
        SequenceCollection(self.seqs3)

    def test_init_validate(self):
        """ initialization with validation functions as expected
        """
        SequenceCollection(self.seqs1, validate=True)
        SequenceCollection(self.seqs1, validate=True)
        # can't validate self.seqs2 as a DNASequence
        self.assertRaises(SequenceCollectionError, SequenceCollection,
                self.invalid_s1, validate=True)

    def test_from_fasta_records(self):
        """ Initialization from list of tuples functions as expected
        """
        SequenceCollection.from_fasta_records(self.seqs1_t, DNASequence)
        SequenceCollection.from_fasta_records(self.seqs2_t, RNASequence)
        SequenceCollection.from_fasta_records(self.seqs3_t, NucleotideSequence)

    def test_eq(self):
        """ equality operator functions as expected
        """
        self.assertTrue(self.s1 == self.s1)
        self.assertFalse(self.s1 == self.s2)
       
        # different objects can be equal
        self.assertTrue(self.s1 == 
                SequenceCollection([self.d1, self.d2]))
        self.assertTrue(SequenceCollection([self.d1, self.d2]) 
                == self.s1)

        # SequenceCollections with different number of sequences are not equal
        self.assertFalse(self.s1 == SequenceCollection([self.d1]))
        
        class FakeSequenceCollection(SequenceCollection):
            pass
        # SequenceCollections of different types are not equal
        self.assertFalse(self.s1 == FakeSequenceCollection([self.d1, self.d2]))
        self.assertFalse(self.s1 == Alignment([self.d1, self.d2]))

        # SequenceCollections with different sequences are not equal
        self.assertFalse(self.s1 == SequenceCollection([self.d1, self.r1]))

    def test_getitem(self):
        """ getitem functions as expected
        """
        self.assertEqual(self.s1[0], self.d1)
        self.assertEqual(self.s1[1], self.d2)
        self.assertEqual(self.s2[0], self.r1)
        self.assertEqual(self.s2[1], self.r2)
   
    def test_iter(self):
        """ iter functions as expected
        """
        s1_iter = iter(self.s1)
        count = 0
        for actual, expected in zip(s1_iter, self.seqs1):
            count += 1
            self.assertEqual(actual, expected)
        self.assertEqual(count, len(self.seqs1))
        self.assertRaises(StopIteration, s1_iter.next)

    def test_len(self):
        """ len functions as expected
        """
        self.assertEqual(len(self.s1),2)
        self.assertEqual(len(self.s2),3)
        self.assertEqual(len(self.s3),5)

    def test_ne(self):
        """ inequality operator functions as expected
        """
        self.assertFalse(self.s1 != self.s1)
        self.assertTrue(self.s1 != self.s2)
       
        # SequenceCollections with different number of sequences are not equal
        self.assertTrue(self.s1 != SequenceCollection([self.d1]))
        
        class FakeSequenceCollection(SequenceCollection):
            pass
        # SequenceCollections of different types are not equal
        self.assertTrue(self.s1 != FakeSequenceCollection([self.d1, self.d2]))
        self.assertTruee(self.s1 != Alignment([self.d1, self.d2]))

        # SequenceCollections with different sequences are not equal
        self.assertTrue(self.s1 != SequenceCollection([self.d1, self]))

    def test_repr(self):
        """
        """
        self.assertEqual(repr(self.s1), 
                "<SequenceCollection: n=2; mean +/- std length=5.00 +/- 2.00>")
        self.assertEqual(repr(self.s2), 
                "<SequenceCollection: n=3; mean +/- std length=7.33 +/- 3.68>")
        self.assertEqual(repr(self.s3), 
                "<SequenceCollection: n=5; mean +/- std length=6.40 +/- 3.32>")
   
    def test_count_center_spread(self):
        """ count_center_spread functions as expected
        """
        actual1 = self.s1.count_center_spread()
        self.assertEqual(actual1[0],2)
        self.assertAlmostEqual(actual1[1], 5.0, 3)
        self.assertAlmostEqual(actual1[2], 2.0, 3)

        actual2 = self.s2.count_center_spread()
        self.assertEqual(actual2[0],3)
        self.assertAlmostEqual(actual2[1], 7.333, 3)
        self.assertAlmostEqual(actual2[2], 3.682, 3)

        actual3 = self.s3.count_center_spread()
        self.assertEqual(actual3[0],5)
        self.assertAlmostEqual(actual3[1], 6.400, 3)
        self.assertAlmostEqual(actual3[2], 3.323, 3)

    def test_degap(self):
        """ degap functions as expected
        """
        expected = [(id_, seq.replace('.', '').replace('-', '')) 
                for id_, seq in self.seqs2_t]
        expected = SequenceCollection.from_fasta_records(expected, RNASequence)
        actual = self.s2.degap()
        self.assertEqual(actual, expected)

    def test_get_seq(self):
        """ getseq functions asexpected
        """
        self.assertEqual(self.s1.get_seq('d1'), self.d1)
        self.assertEqual(self.s1.get_seq('d2'), self.d2)

    def test_identifiers(self):
        """ identifiers functions as expected
        """
        self.assertEqual(sorted(self.s1.identifiers()), ['d1', 'd2']) 
        self.assertEqual(sorted(self.s2.identifiers()), ['r1', 'r2', 'r3']) 
        self.assertEqual(sorted(self.s3.identifiers()),
                ['d1', 'd2', 'r1', 'r2', 'r3']) 

    def test_int_map(self):
        """ int_map functions as expected
        """
        expected1 = {"0": self.d1, "1": self.d2}
        expected2 = {"0": "d1", "1": "d2"}
        self.assertEqual(self.s1.int_map(), (expected1, expected2))
        
        expected1 = {"h-0": self.d1, "h-1": self.d2}
        expected2 = {"h-0": "d1", "h-1": "d2"}
        self.assertEqual(self.s1.int_map(prefix='h-'), (expected1, expected2))

    def test_is_valid(self):
        """ is_valid functions as expected
        """
        self.assertTrue(self.s1.is_valid())
        self.assertTrue(self.s2.is_valid())
        self.assertTrue(self.s3.is_valid())

        self.assertFalse(self.invalid_s1.is_valid())

    def test_items(self):
        """ items functions as expected
        """
        self.assertEqual(list(self.s1.items()),
                [(s.identifier, s) for s in self.s1])

    def test_lower(self):
        """ lower functions as expected
        """
        self.assertEqual(self.s1.lower(), self.s1_lower)

    def test_sequence_count(self):
        """ num_seqs functions as expected
        """
        self.assertEqual(self.s1.sequence_count(), 2)
        self.assertEqual(self.s2.sequence_count(), 3)
        self.assertEqual(self.s3.sequence_count(), 5)

    def test_sequence_lengths(self):
        """ sequence_lengths functions as expected
        """
        self.assertEqual(self.s1.sequence_lengths(), [7, 3])
        self.assertEqual(self.s2.sequence_lengths(), [7, 3, 12])
        self.assertEqual(self.s3.sequence_lengths(), [7, 3, 7, 3, 12])

    def test_to_fasta(self):
        """ to_fasta functions as expected
        """
        exp1 = ">d1\nGATTACA\n>d2\nTTG\n"
        self.assertEqual(self.s1.to_fasta(),exp1)
        exp2 = ">r1\nGATTACA\n>r2\nTTG\n>r3'U-----UGCC--\n"
        self.assertEqual(self.s1.to_fasta(),exp1)

    def test_to_phylip(self):
        """ to_phylip functions as expected
        """
        raise NotImplementedError

    def test_upper(self):
        """ upper functions as expected
        """
        self.assertEqual(self.s1_lower.upper(), self.s1)

class AlignmentTests(TestCase):

    def setUp(self):
        self.d1 = DNASequence('..ACC-GTTGG..', identifier="d1")
        self.d2 = DNASequence('TTACCGGT-GGCC', identifier="d2")
        self.d3 = DNASequence('.-ACC-GTTGC--', identifier="d3")
        
        self.r1 = DNASequence('UUAU-', identifier="r1")
        self.r2 = DNASequence('ACGUU', identifier="r2")
        
        self.seqs1 = [self.d1, self.d2, self.d3]
        self.seqs2 = [self.r1, self.r2]

        self.seqs1_t = [('d1', '..ACC-GTTGG..'), ('d2', 'TTACCGGT-GGCC'),
                ('d3', '.-ACC-GTTGC--')]
        self.seqs2_t = [('r1', 'UUAU-'), ('r2', 'ACGUU')]

        self.a1 = Alignment(self.seqs1)
        self.a2 = Alignment(self.seqs2)

    def test_degap(self):
        """ degap functions as expected
        """
        expected = [(id_, seq.replace('.', '').replace('-', '')) 
                for id_, seq in self.seqs1_t]
        expected = SequenceCollection.from_fasta_records(expected, DNASequence)
        actual = self.a1.degap()
        self.assertEqual(actual, expected)
        
        expected = [(id_, seq.replace('.', '').replace('-', '')) 
                for id_, seq in self.seqs2_t]
        expected = SequenceCollection.from_fasta_records(expected, RNASequence)
        actual = self.a2.degap()
        print expected.to_fasta()
        print actual.to_fasta()
        self.assertEqual(actual, expected)

    def test_distances(self):
        """ distances functions as expected
        """
        expected = [[    0, 6./13, 4./13],
                    [6./13, 0,     7./13],
                    [4./13, 7./13, 0]]
        expected = SymmetricDistanceMatrix(expected, ['d1', 'd2', 'd3'])
        actual = self.a1.distances()
        self.assertEqual(actual, expected)

    def test_get_subalignment(self):
        """ get_sub_alignment functions as expected
        """
        raise NotImplementedError

    def test_init_validate(self):
        """ initialization with validation functions as expected
        """
        Alignment(self.seqs1, validate=True)

        # invalid DNA character
        invalid_seqs1 = [self.d1, self.d2, self.d3,
                DNASequence('.-ACC-GTXGC--', identifier="i1")]
        self.assertRaises(SequenceCollectionError, Alignment,
                invalid_seqs1, validate=True)

        # invalid lengths (they're not all equal)
        invalid_seqs2 = [self.d1, self.d2, self.d3,
                DNASequence('.-ACC-GTGC--', identifier="i2")]
        self.assertRaises(SequenceCollectionError, Alignment,
                invalid_seqs2, validate=True)

    def test_is_valid(self):
        """
        """
        raise NotImplementedError

    def test_iter_positions(self):
        """ iter_positions functions as expected
        """
        actual = list(self.a2.iter_positions())
        expected = [map(DNASequence,list('UA')), 
                    map(DNASequence,list('UC')),
                    map(DNASequence,list('AG')),
                    map(DNASequence,list('UU')),
                    map(DNASequence,list('-U'))]
        self.seqs2_t = [('r1', 'UUAU-'), ('r2', 'ACGUU')]
        self.assertEqual(actual, expected)

        actual = list(self.a2.iter_positions(constructor=str))
        expected = [list('UA'), 
                    list('UC'),
                    list('AG'),
                    list('UU'),
                    list('-U')]
        self.seqs2_t = [('r1', 'UUAU-'), ('r2', 'ACGUU')]
        self.assertEqual(actual, expected)

    def test_majority_consensus(self):
        """ majority_consensus functions as expected
        """
        d1 = DNASequence('TTT', identifier="d1")
        d2 = DNASequence('TT-', identifier="d2")
        d3 = DNASequence('TC-', identifier="d3")
        a1 = Alignment([d1, d2, d3])
        self.assertEqual(a1.majority_consensus(), DNASequence('TT-'))
        
        d1 = DNASequence('T', identifier="d1")
        d2 = DNASequence('A', identifier="d2")
        a1 = Alignment([d1, d2])
        self.assertTrue(a1.majority_consensus() in 
                [DNASequence('T'), DNASequence('A')])

    def test_omit_gap_positions(self):
        """ omitting gap positions functions as expected
        """
        raise NotImplementedError
   
    def test_omit_gap_sequences(self):
        """ omitting gap sequences functions as expected
        """
        raise NotImplementedError

    def test_position_counters(self):
        """ position_counters functions as expected
        """
        expected = [Counter({'U': 1, 'A': 1}),
                    Counter({'U': 1, 'C': 1}),
                    Counter({'A': 1, 'G': 1}),
                    Counter({'U': 2}),
                    Counter({'-': 1, 'U': 1})]
        self.assertEqual(self.a2.position_counters(), expected)

    def test_position_frequencies(self):
        """ computing position frequencies functions as expected
        """
        expected = [defaultdict(int, {'U': 0.5, 'A': 0.5}),
                    defaultdict(int, {'U': 0.5, 'C': 0.5}),
                    defaultdict(int, {'A': 0.5, 'G': 0.5}),
                    defaultdict(int, {'U': 1.0}),
                    defaultdict(int, {'-': 0.5, 'U': 0.5})]
        self.assertEqual(self.a2.position_frequencies(), expected)

    def test_position_entropies(self):
        """ computing positional uncertainties functions as expected
        """
        raise NotImplementedError

    def test_sequence_length(self):
        """
        """
        raise NotImplementedError

    def test_validate_lengths(self):
        """
        """
        raise NotImplementedError



if __name__ == "__main__":
    main()
