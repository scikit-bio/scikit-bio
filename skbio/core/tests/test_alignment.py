#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division

from unittest import TestCase, main
from collections import Counter, defaultdict

import numpy as np

from skbio.core.sequence import NucleotideSequence, DNASequence, RNASequence
from skbio.core.alignment import SequenceCollection, Alignment
from skbio.core.exception import SequenceCollectionError
from skbio.core.distance import DistanceMatrix


class SequenceCollectionTests(TestCase):
    """Tests of the SequenceCollection class """

    def setUp(self):
        """Initialize values to be used in tests
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
        self.empty = SequenceCollection([])

        self.invalid_s1 = SequenceCollection([self.i1])

    def test_init(self):
        """Initialization functions as expected with varied input types
        """
        SequenceCollection(self.seqs1)
        SequenceCollection(self.seqs2)
        SequenceCollection(self.seqs3)
        SequenceCollection([])

    def test_init_fail(self):
        """initialization with sequences with overlapping identifiers fails
        """
        s1 = [self.d1, self.d1]
        self.assertRaises(SequenceCollectionError, SequenceCollection, s1)

    def test_init_validate(self):
        """initialization with validation functions as expected
        """
        SequenceCollection(self.seqs1, validate=True)
        SequenceCollection(self.seqs1, validate=True)
        # can't validate self.seqs2 as a DNASequence
        self.assertRaises(SequenceCollectionError, SequenceCollection,
                          self.invalid_s1, validate=True)

    def test_from_fasta_records(self):
        """Initialization from list of tuples functions as expected
        """
        SequenceCollection.from_fasta_records(self.seqs1_t, DNASequence)
        SequenceCollection.from_fasta_records(self.seqs2_t, RNASequence)
        SequenceCollection.from_fasta_records(self.seqs3_t, NucleotideSequence)

    def test_contains(self):
        """in operator functions as expected
        """
        self.assertTrue('d1' in self.s1)
        self.assertTrue('r2' in self.s2)
        self.assertFalse('r2' in self.s1)

    def test_eq(self):
        """equality operator functions as expected
        """
        self.assertTrue(self.s1 == self.s1)
        self.assertFalse(self.s1 == self.s2)

        # different objects can be equal
        self.assertTrue(self.s1 == SequenceCollection([self.d1, self.d2]))
        self.assertTrue(SequenceCollection([self.d1, self.d2]) == self.s1)

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
        """getitem functions as expected
        """
        self.assertEqual(self.s1[0], self.d1)
        self.assertEqual(self.s1[1], self.d2)
        self.assertEqual(self.s2[0], self.r1)
        self.assertEqual(self.s2[1], self.r2)

        self.assertRaises(IndexError, self.empty.__getitem__, 0)
        self.assertRaises(KeyError, self.empty.__getitem__, '0')

    def test_iter(self):
        """iter functions as expected
        """
        s1_iter = iter(self.s1)
        count = 0
        for actual, expected in zip(s1_iter, self.seqs1):
            count += 1
            self.assertEqual(actual, expected)
        self.assertEqual(count, len(self.seqs1))
        self.assertRaises(StopIteration, next(s1_iter))

    def test_len(self):
        """len functions as expected
        """
        self.assertEqual(len(self.s1), 2)
        self.assertEqual(len(self.s2), 3)
        self.assertEqual(len(self.s3), 5)
        self.assertEqual(len(self.empty), 0)

    def test_ne(self):
        """inequality operator functions as expected
        """
        self.assertFalse(self.s1 != self.s1)
        self.assertTrue(self.s1 != self.s2)

        # SequenceCollections with different number of sequences are not equal
        self.assertTrue(self.s1 != SequenceCollection([self.d1]))

        class FakeSequenceCollection(SequenceCollection):
            pass
        # SequenceCollections of different types are not equal
        self.assertTrue(self.s1 != FakeSequenceCollection([self.d1, self.d2]))
        self.assertTrue(self.s1 != Alignment([self.d1, self.d2]))

        # SequenceCollections with different sequences are not equal
        self.assertTrue(self.s1 !=
                        SequenceCollection([self.d1, self.r1]))

    def test_repr(self):
        """repr functions as expected
        """
        self.assertEqual(repr(self.s1),
                         "<SequenceCollection: n=2; "
                         "mean +/- std length=5.00 +/- 2.00>")
        self.assertEqual(repr(self.s2),
                         "<SequenceCollection: n=3; "
                         "mean +/- std length=7.33 +/- 3.68>")
        self.assertEqual(repr(self.s3),
                         "<SequenceCollection: n=5; "
                         "mean +/- std length=6.40 +/- 3.32>")
        self.assertEqual(repr(self.empty),
                         "<SequenceCollection: n=0; "
                         "mean +/- std length=0.00 +/- 0.00>")

    def test_reversed(self):
        """reversed functions as expected
        """
        s1_iter = reversed(self.s1)
        count = 0
        for actual, expected in zip(s1_iter, self.seqs1[::-1]):
            count += 1
            self.assertEqual(actual, expected)
        self.assertEqual(count, len(self.seqs1))
        self.assertRaises(StopIteration, next(s1_iter))

    def test_k_word_frequencies(self):
        """k_word_frequencies functions as expected
        """
        expected1 = defaultdict(int)
        expected1['A'] = 3/7.
        expected1['C'] = 1/7.
        expected1['G'] = 1/7.
        expected1['T'] = 2/7.
        expected2 = defaultdict(int)
        expected2['G'] = 1/3.
        expected2['T'] = 2/3.
        self.assertEqual(self.s1.k_word_frequencies(k=1),
                         [expected1, expected2])

        expected1 = defaultdict(int)
        expected1['GAT'] = 1/2.
        expected1['TAC'] = 1/2.
        expected2 = defaultdict(int)
        expected2['TTG'] = 1/1.
        self.assertEqual(self.s1.k_word_frequencies(k=3, overlapping=False),
                         [expected1, expected2])

        self.assertEqual(self.empty.k_word_frequencies(k=1), [])

    def test_str(self):
        """str functions as expected
        """
        exp1 = ">d1\nGATTACA\n>d2\nTTG\n"
        self.assertEqual(str(self.s1), exp1)
        exp2 = ">r1\nGAUUACA\n>r2\nUUG\n>r3\nU-----UGCC--\n"
        self.assertEqual(str(self.s2), exp2)
        exp4 = ""
        self.assertEqual(str(self.empty), exp4)

    def test_distribution_stats(self):
        """distribution_stats functions as expected
        """
        actual1 = self.s1.distribution_stats()
        self.assertEqual(actual1[0], 2)
        self.assertAlmostEqual(actual1[1], 5.0, 3)
        self.assertAlmostEqual(actual1[2], 2.0, 3)

        actual2 = self.s2.distribution_stats()
        self.assertEqual(actual2[0], 3)
        self.assertAlmostEqual(actual2[1], 7.333, 3)
        self.assertAlmostEqual(actual2[2], 3.682, 3)

        actual3 = self.s3.distribution_stats()
        self.assertEqual(actual3[0], 5)
        self.assertAlmostEqual(actual3[1], 6.400, 3)
        self.assertAlmostEqual(actual3[2], 3.323, 3)

        actual4 = self.empty.distribution_stats()
        self.assertEqual(actual4[0], 0)
        self.assertEqual(actual4[1], 0.0)
        self.assertEqual(actual4[2], 0.0)

    def test_degap(self):
        """degap functions as expected
        """
        expected = [(id_, seq.replace('.', '').replace('-', ''))
                    for id_, seq in self.seqs2_t]
        expected = SequenceCollection.from_fasta_records(expected, RNASequence)
        actual = self.s2.degap()
        self.assertEqual(actual, expected)

    def test_get_seq(self):
        """getseq functions asexpected
        """
        self.assertEqual(self.s1.get_seq('d1'), self.d1)
        self.assertEqual(self.s1.get_seq('d2'), self.d2)

    def test_identifiers(self):
        """identifiers functions as expected
        """
        self.assertEqual(self.s1.identifiers(), ['d1', 'd2'])
        self.assertEqual(self.s2.identifiers(), ['r1', 'r2', 'r3'])
        self.assertEqual(self.s3.identifiers(),
                         ['d1', 'd2', 'r1', 'r2', 'r3'])
        self.assertEqual(self.empty.identifiers(), [])

    def test_int_map(self):
        """int_map functions as expected
        """
        expected1 = {"1": self.d1, "2": self.d2}
        expected2 = {"1": "d1", "2": "d2"}
        self.assertEqual(self.s1.int_map(), (expected1, expected2))

        expected1 = {"h-1": self.d1, "h-2": self.d2}
        expected2 = {"h-1": "d1", "h-2": "d2"}
        self.assertEqual(self.s1.int_map(prefix='h-'), (expected1, expected2))

    def test_is_empty(self):
        """is_empty functions as expected
        """
        self.assertFalse(self.s1.is_empty())
        self.assertFalse(self.s2.is_empty())
        self.assertFalse(self.s3.is_empty())

        self.assertTrue(self.empty.is_empty())

    def test_is_valid(self):
        """is_valid functions as expected
        """
        self.assertTrue(self.s1.is_valid())
        self.assertTrue(self.s2.is_valid())
        self.assertTrue(self.s3.is_valid())
        self.assertTrue(self.empty.is_valid())

        self.assertFalse(self.invalid_s1.is_valid())

    def test_iteritems(self):
        """iteritems functions as expected
        """
        self.assertEqual(list(self.s1.iteritems()),
                         [(s.identifier, s) for s in self.s1])

    def test_lower(self):
        """lower functions as expected
        """
        self.assertEqual(self.s1.lower(), self.s1_lower)

    def test_sequence_count(self):
        """num_seqs functions as expected
        """
        self.assertEqual(self.s1.sequence_count(), 2)
        self.assertEqual(self.s2.sequence_count(), 3)
        self.assertEqual(self.s3.sequence_count(), 5)
        self.assertEqual(self.empty.sequence_count(), 0)

    def test_sequence_lengths(self):
        """sequence_lengths functions as expected
        """
        self.assertEqual(self.s1.sequence_lengths(), [7, 3])
        self.assertEqual(self.s2.sequence_lengths(), [7, 3, 12])
        self.assertEqual(self.s3.sequence_lengths(), [7, 3, 7, 3, 12])
        self.assertEqual(self.empty.sequence_lengths(), [])

    def test_to_fasta(self):
        """to_fasta functions as expected
        """
        exp1 = ">d1\nGATTACA\n>d2\nTTG\n"
        self.assertEqual(self.s1.to_fasta(), exp1)
        exp2 = ">r1\nGAUUACA\n>r2\nUUG\n>r3\nU-----UGCC--\n"
        self.assertEqual(self.s2.to_fasta(), exp2)

    def test_upper(self):
        """upper functions as expected
        """
        self.assertEqual(self.s1_lower.upper(), self.s1)


class AlignmentTests(TestCase):

    def setUp(self):
        self.d1 = DNASequence('..ACC-GTTGG..', identifier="d1")
        self.d2 = DNASequence('TTACCGGT-GGCC', identifier="d2")
        self.d3 = DNASequence('.-ACC-GTTGC--', identifier="d3")

        self.r1 = RNASequence('UUAU-', identifier="r1")
        self.r2 = RNASequence('ACGUU', identifier="r2")

        self.seqs1 = [self.d1, self.d2, self.d3]
        self.seqs2 = [self.r1, self.r2]

        self.seqs1_t = [('d1', '..ACC-GTTGG..'), ('d2', 'TTACCGGT-GGCC'),
                        ('d3', '.-ACC-GTTGC--')]
        self.seqs2_t = [('r1', 'UUAU-'), ('r2', 'ACGUU')]

        self.a1 = Alignment(self.seqs1)
        self.a2 = Alignment(self.seqs2)
        self.empty = Alignment([])

    def test_degap(self):
        """degap functions as expected
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
        self.assertEqual(actual, expected)

    def test_distances(self):
        """distances functions as expected
        """
        expected = [[0, 6./13, 4./13],
                    [6./13, 0, 7./13],
                    [4./13, 7./13, 0]]
        expected = DistanceMatrix(expected, ['d1', 'd2', 'd3'])
        actual = self.a1.distances()
        self.assertEqual(actual, expected)

    def test_subalignment(self):
        """subalignment functions as expected
        """
        # keep seqs by identifiers
        actual = self.a1.subalignment(seqs_to_keep=['d1', 'd3'])
        expected = Alignment([self.d1, self.d3])
        self.assertEqual(actual, expected)

        # keep seqs by indices
        actual = self.a1.subalignment(seqs_to_keep=[0, 2])
        expected = Alignment([self.d1, self.d3])
        self.assertEqual(actual, expected)

        # keep seqs by identifiers (invert)
        actual = self.a1.subalignment(seqs_to_keep=['d1', 'd3'],
                                      invert_seqs_to_keep=True)
        expected = Alignment([self.d2])
        self.assertEqual(actual, expected)

        # keep seqs by indices (invert)
        actual = self.a1.subalignment(seqs_to_keep=[0, 2],
                                      invert_seqs_to_keep=True)
        expected = Alignment([self.d2])
        self.assertEqual(actual, expected)

        # keep positions
        actual = self.a1.subalignment(positions_to_keep=[0, 2, 3])
        d1 = DNASequence('.AC', identifier="d1")
        d2 = DNASequence('TAC', identifier="d2")
        d3 = DNASequence('.AC', identifier="d3")
        expected = Alignment([d1, d2, d3])
        self.assertEqual(actual, expected)

        # keep positions (invert)
        actual = self.a1.subalignment(positions_to_keep=[0, 2, 3],
                                      invert_positions_to_keep=True)
        d1 = DNASequence('.C-GTTGG..', identifier="d1")
        d2 = DNASequence('TCGGT-GGCC', identifier="d2")
        d3 = DNASequence('-C-GTTGC--', identifier="d3")
        expected = Alignment([d1, d2, d3])
        self.assertEqual(actual, expected)

        # keep seqs and positions
        actual = self.a1.subalignment(seqs_to_keep=[0, 2],
                                      positions_to_keep=[0, 2, 3])
        d1 = DNASequence('.AC', identifier="d1")
        d3 = DNASequence('.AC', identifier="d3")
        expected = Alignment([d1, d3])
        self.assertEqual(actual, expected)

        # keep seqs and positions (invert)
        actual = self.a1.subalignment(seqs_to_keep=[0, 2],
                                      positions_to_keep=[0, 2, 3],
                                      invert_seqs_to_keep=True,
                                      invert_positions_to_keep=True)
        d2 = DNASequence('TCGGT-GGCC', identifier="d2")
        expected = Alignment([d2])
        self.assertEqual(actual, expected)

    def test_init_validate(self):
        """initialization with validation functions as expected
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
        """is_valid functions as expected
        """
        self.assertTrue(self.a1.is_valid())
        self.assertTrue(self.a2.is_valid())
        self.assertTrue(self.empty.is_valid())

        # invalid because of length mismatch
        d1 = DNASequence('..ACC-GTTGG..', identifier="d1")
        d2 = DNASequence('TTACCGGT-GGC', identifier="d2")
        self.assertFalse(Alignment([d1, d2]).is_valid())

        # invalid because of invalid charaters
        d1 = DNASequence('..ACC-GTXGG..', identifier="d1")
        d2 = DNASequence('TTACCGGT-GGCC', identifier="d2")
        self.assertFalse(Alignment([d1, d2]).is_valid())

    def test_iter_positions(self):
        """iter_positions functions as expected
        """
        actual = list(self.a2.iter_positions())
        expected = [list(map(RNASequence, list(i))) for i in
                    ['UA', 'UC', 'AG', 'UU', '-U']]
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
        """majority_consensus functions as expected
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

        self.assertEqual(self.empty.majority_consensus(), '')

    def test_omit_gap_positions(self):
        """omitting gap positions functions as expected
        """
        expected = self.a2
        self.assertEqual(self.a2.omit_gap_positions(1.0), expected)
        self.assertEqual(self.a2.omit_gap_positions(0.51), expected)

        r1 = RNASequence('UUAU', identifier="r1")
        r2 = RNASequence('ACGU', identifier="r2")
        expected = Alignment([r1, r2])
        self.assertEqual(self.a2.omit_gap_positions(0.49), expected)

        r1 = RNASequence('UUAU', identifier="r1")
        r2 = RNASequence('ACGU', identifier="r2")
        expected = Alignment([r1, r2])
        self.assertEqual(self.a2.omit_gap_positions(0.0), expected)

        self.assertEqual(self.empty.omit_gap_positions(0.0), self.empty)
        self.assertEqual(self.empty.omit_gap_positions(0.49), self.empty)
        self.assertEqual(self.empty.omit_gap_positions(1.0), self.empty)

    def test_omit_gap_sequences(self):
        """omitting gap sequences functions as expected
        """
        expected = self.a2
        self.assertEqual(self.a2.omit_gap_sequences(1.0), expected)
        self.assertEqual(self.a2.omit_gap_sequences(0.20), expected)

        expected = Alignment([self.r2])
        self.assertEqual(self.a2.omit_gap_sequences(0.19), expected)

        self.assertEqual(self.empty.omit_gap_sequences(0.0), self.empty)
        self.assertEqual(self.empty.omit_gap_sequences(0.2), self.empty)
        self.assertEqual(self.empty.omit_gap_sequences(1.0), self.empty)

    def test_position_counters(self):
        """position_counters functions as expected
        """
        expected = [Counter({'U': 1, 'A': 1}),
                    Counter({'U': 1, 'C': 1}),
                    Counter({'A': 1, 'G': 1}),
                    Counter({'U': 2}),
                    Counter({'-': 1, 'U': 1})]
        self.assertEqual(self.a2.position_counters(), expected)

        self.assertEqual(self.empty.position_counters(), [])

    def test_position_frequencies(self):
        """computing position frequencies functions as expected
        """
        expected = [defaultdict(int, {'U': 0.5, 'A': 0.5}),
                    defaultdict(int, {'U': 0.5, 'C': 0.5}),
                    defaultdict(int, {'A': 0.5, 'G': 0.5}),
                    defaultdict(int, {'U': 1.0}),
                    defaultdict(int, {'-': 0.5, 'U': 0.5})]
        self.assertEqual(self.a2.position_frequencies(), expected)

        self.assertEqual(self.empty.position_frequencies(), [])

    def test_position_entropies(self):
        """computing positional uncertainties functions as expected

        tested by calculating values as described in this post:
         http://stackoverflow.com/a/15476958/3424666
        """
        expected = [0.69314, 0.69314, 0.69314, 0.0, np.nan]
        np.testing.assert_almost_equal(self.a2.position_entropies(),
                                       expected, 5)

        expected = [1.0, 1.0, 1.0, 0.0, np.nan]
        np.testing.assert_almost_equal(self.a2.position_entropies(base=2),
                                       expected, 5)

        np.testing.assert_almost_equal(self.empty.position_entropies(base=2),
                                       [])

    def test_k_word_frequencies(self):
        """k_word_frequencies functions as expected
        """
        expected = [defaultdict(int, {'U': 3/5, 'A': 1/5, '-': 1/5}),
                    defaultdict(int, {'A': 1/5, 'C': 1/5, 'G': 1/5, 'U': 2/5})]
        actual = self.a2.k_word_frequencies(k=1)
        for a, e in zip(actual, expected):
            self.assertEqual(sorted(a.keys()), sorted(e.keys()), 5)
            np.testing.assert_almost_equal(sorted(a.values()),
                                           sorted(e.values()), 5)

    def test_sequence_length(self):
        """sequence_length functions as expected
        """
        self.assertEqual(self.a1.sequence_length(), 13)
        self.assertEqual(self.a2.sequence_length(), 5)
        self.assertEqual(self.empty.sequence_length(), 0)

    def test_to_phylip(self):
        """to_phylip functions as expected
        """
        d1 = DNASequence('..ACC-GTTGG..', identifier="d1")
        d2 = DNASequence('TTACCGGT-GGCC', identifier="d2")
        d3 = DNASequence('.-ACC-GTTGC--', identifier="d3")
        a = Alignment([d1, d2, d3])

        phylip_str, id_map = a.to_phylip(map_labels=False)
        self.assertEqual(id_map, {'d1': 'd1',
                                  'd3': 'd3',
                                  'd2': 'd2'})
        expected = "\n".join(["3 13",
                              "d1 ..ACC-GTTGG..",
                              "d2 TTACCGGT-GGCC",
                              "d3 .-ACC-GTTGC--"])
        self.assertEqual(phylip_str, expected)

    def test_to_phylip_map_labels(self):
        """to_phylip functions as expected with label mapping
        """
        d1 = DNASequence('..ACC-GTTGG..', identifier="d1")
        d2 = DNASequence('TTACCGGT-GGCC', identifier="d2")
        d3 = DNASequence('.-ACC-GTTGC--', identifier="d3")
        a = Alignment([d1, d2, d3])

        phylip_str, id_map = a.to_phylip(map_labels=True, label_prefix="s")
        self.assertEqual(id_map, {'s1': 'd1',
                                  's3': 'd3',
                                  's2': 'd2'})
        expected = "\n".join(["3 13",
                              "s1 ..ACC-GTTGG..",
                              "s2 TTACCGGT-GGCC",
                              "s3 .-ACC-GTTGC--"])
        self.assertEqual(phylip_str, expected)

    def test_validate_lengths(self):
        """
        """
        self.assertTrue(self.a1._validate_lengths())
        self.assertTrue(self.a2._validate_lengths())
        self.assertTrue(self.empty._validate_lengths())

        self.assertTrue(Alignment([
            DNASequence('TTT', identifier="d1")])._validate_lengths())
        self.assertFalse(Alignment([
            DNASequence('TTT', identifier="d1"),
            DNASequence('TT', identifier="d2")])._validate_lengths())


if __name__ == "__main__":
    main()
