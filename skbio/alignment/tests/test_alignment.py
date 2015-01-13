# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main
from collections import Counter, defaultdict, OrderedDict
try:
    from StringIO import StringIO
except ImportError:  # python3 system
    from io import StringIO
import tempfile

import numpy as np
import numpy.testing as npt
from scipy.spatial.distance import hamming
import matplotlib as mpl
mpl.use('Agg')

from skbio import (NucleotideSequence, DNASequence, RNASequence, DNA, RNA,
                   DistanceMatrix, Alignment, SequenceCollection)
from skbio.alignment import (StockholmAlignment, SequenceCollectionError,
                             StockholmParseError, AlignmentError)


class SequenceCollectionTests(TestCase):
    def setUp(self):
        self.d1 = DNASequence('GATTACA', id="d1")
        self.d2 = DNASequence('TTG', id="d2")
        self.d3 = DNASequence('GTATACA', id="d3")
        self.d1_lower = DNASequence('gattaca', id="d1")
        self.d2_lower = DNASequence('ttg', id="d2")
        self.d3_lower = DNASequence('gtataca', id="d3")
        self.r1 = RNASequence('GAUUACA', id="r1")
        self.r2 = RNASequence('UUG', id="r2")
        self.r3 = RNASequence('U-----UGCC--', id="r3")

        self.i1 = DNASequence('GATXACA', id="i1")

        self.seqs1 = [self.d1, self.d2]
        self.seqs1_lower = [self.d1_lower, self.d2_lower]
        self.seqs2 = [self.r1, self.r2, self.r3]
        self.seqs3 = self.seqs1 + self.seqs2
        self.seqs4 = [self.d1, self.d3]

        self.seqs1_t = [('d1', 'GATTACA'), ('d2', 'TTG')]
        self.seqs2_t = [('r1', 'GAUUACA'), ('r2', 'UUG'),
                        ('r3', 'U-----UGCC--')]
        self.seqs3_t = self.seqs1_t + self.seqs2_t

        self.s1 = SequenceCollection(self.seqs1)
        self.s1_lower = SequenceCollection(self.seqs1_lower)
        self.s2 = SequenceCollection(self.seqs2)
        self.s3 = SequenceCollection(self.seqs3)
        self.s4 = SequenceCollection(self.seqs4)
        self.empty = SequenceCollection([])

        self.invalid_s1 = SequenceCollection([self.i1])

    def test_init(self):
        SequenceCollection(self.seqs1)
        SequenceCollection(self.seqs2)
        SequenceCollection(self.seqs3)
        SequenceCollection([])

    def test_init_fail(self):
        # sequences with overlapping ids
        s1 = [self.d1, self.d1]
        self.assertRaises(SequenceCollectionError, SequenceCollection, s1)

    def test_init_validate(self):
        SequenceCollection(self.seqs1, validate=True)
        SequenceCollection(self.seqs1, validate=True)
        # can't validate self.seqs2 as a DNASequence
        self.assertRaises(SequenceCollectionError, SequenceCollection,
                          self.invalid_s1, validate=True)

    def test_from_fasta_records(self):
        SequenceCollection.from_fasta_records(self.seqs1_t, DNASequence)
        SequenceCollection.from_fasta_records(self.seqs2_t, RNASequence)
        SequenceCollection.from_fasta_records(self.seqs3_t, NucleotideSequence)

    def test_contains(self):
        self.assertTrue('d1' in self.s1)
        self.assertTrue('r2' in self.s2)
        self.assertFalse('r2' in self.s1)

    def test_eq(self):
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
        self.assertFalse(self.s4 == FakeSequenceCollection([self.d1, self.d3]))
        self.assertFalse(self.s4 == Alignment([self.d1, self.d3]))

        # SequenceCollections with different sequences are not equal
        self.assertFalse(self.s1 == SequenceCollection([self.d1, self.r1]))

    def test_getitem(self):
        self.assertEqual(self.s1[0], self.d1)
        self.assertEqual(self.s1[1], self.d2)
        self.assertEqual(self.s2[0], self.r1)
        self.assertEqual(self.s2[1], self.r2)

        self.assertRaises(IndexError, self.empty.__getitem__, 0)
        self.assertRaises(KeyError, self.empty.__getitem__, '0')

    def test_iter(self):
        s1_iter = iter(self.s1)
        count = 0
        for actual, expected in zip(s1_iter, self.seqs1):
            count += 1
            self.assertEqual(actual, expected)
        self.assertEqual(count, len(self.seqs1))
        self.assertRaises(StopIteration, lambda: next(s1_iter))

    def test_len(self):
        self.assertEqual(len(self.s1), 2)
        self.assertEqual(len(self.s2), 3)
        self.assertEqual(len(self.s3), 5)
        self.assertEqual(len(self.empty), 0)

    def test_ne(self):
        self.assertFalse(self.s1 != self.s1)
        self.assertTrue(self.s1 != self.s2)

        # SequenceCollections with different number of sequences are not equal
        self.assertTrue(self.s1 != SequenceCollection([self.d1]))

        class FakeSequenceCollection(SequenceCollection):
            pass
        # SequenceCollections of different types are not equal
        self.assertTrue(self.s4 != FakeSequenceCollection([self.d1, self.d3]))
        self.assertTrue(self.s4 != Alignment([self.d1, self.d3]))

        # SequenceCollections with different sequences are not equal
        self.assertTrue(self.s1 !=
                        SequenceCollection([self.d1, self.r1]))

    def test_repr(self):
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
        s1_iter = reversed(self.s1)
        count = 0
        for actual, expected in zip(s1_iter, self.seqs1[::-1]):
            count += 1
            self.assertEqual(actual, expected)
        self.assertEqual(count, len(self.seqs1))
        self.assertRaises(StopIteration, lambda: next(s1_iter))

    def test_k_word_frequencies(self):
        expected1 = defaultdict(float)
        expected1['A'] = 3 / 7.
        expected1['C'] = 1 / 7.
        expected1['G'] = 1 / 7.
        expected1['T'] = 2 / 7.
        expected2 = defaultdict(float)
        expected2['G'] = 1 / 3.
        expected2['T'] = 2 / 3.
        self.assertEqual(self.s1.k_word_frequencies(k=1),
                         [expected1, expected2])

        expected1 = defaultdict(float)
        expected1['GAT'] = 1 / 2.
        expected1['TAC'] = 1 / 2.
        expected2 = defaultdict(float)
        expected2['TTG'] = 1 / 1.
        self.assertEqual(self.s1.k_word_frequencies(k=3, overlapping=False),
                         [expected1, expected2])

        self.assertEqual(self.empty.k_word_frequencies(k=1), [])

        # Test to ensure floating point precision bug isn't present. See the
        # tests for BiologicalSequence.k_word_frequencies for more details.
        sc = SequenceCollection([RNA('C' * 10, id='s1'),
                                 RNA('G' * 10, id='s2')])
        self.assertEqual(sc.k_word_frequencies(1),
                         [defaultdict(float, {'C': 1.0}),
                          defaultdict(float, {'G': 1.0})])

    def test_str(self):
        exp1 = ">d1\nGATTACA\n>d2\nTTG\n"
        self.assertEqual(str(self.s1), exp1)
        exp2 = ">r1\nGAUUACA\n>r2\nUUG\n>r3\nU-----UGCC--\n"
        self.assertEqual(str(self.s2), exp2)
        exp4 = ""
        self.assertEqual(str(self.empty), exp4)

    def test_distances(self):
        s1 = SequenceCollection([DNA("ACGT", "d1"), DNA("ACGG", "d2")])
        expected = [[0, 0.25],
                    [0.25, 0]]
        expected = DistanceMatrix(expected, ['d1', 'd2'])
        actual = s1.distances(hamming)
        self.assertEqual(actual, expected)

        # alt distance function provided
        def dumb_distance(s1, s2):
            return 42.
        expected = [[0, 42.],
                    [42., 0]]
        expected = DistanceMatrix(expected, ['d1', 'd2'])
        actual = s1.distances(dumb_distance)
        self.assertEqual(actual, expected)

    def test_distribution_stats(self):
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
        expected = [(id_, seq.replace('.', '').replace('-', ''))
                    for id_, seq in self.seqs2_t]
        expected = SequenceCollection.from_fasta_records(expected, RNASequence)
        actual = self.s2.degap()
        self.assertEqual(actual, expected)

    def test_get_seq(self):
        self.assertEqual(self.s1.get_seq('d1'), self.d1)
        self.assertEqual(self.s1.get_seq('d2'), self.d2)

    def test_ids(self):
        self.assertEqual(self.s1.ids(), ['d1', 'd2'])
        self.assertEqual(self.s2.ids(), ['r1', 'r2', 'r3'])
        self.assertEqual(self.s3.ids(),
                         ['d1', 'd2', 'r1', 'r2', 'r3'])
        self.assertEqual(self.empty.ids(), [])

    def _assert_sequence_collections_equal(self, observed, expected):
        """Compare SequenceCollections strictly."""
        # TODO remove this custom equality testing code when SequenceCollection
        # has an equals method (part of #656). We need this method to include
        # IDs in the comparison (not part of SequenceCollection.__eq__).
        self.assertEqual(observed, expected)
        for obs_seq, exp_seq in zip(observed, expected):
            self.assertTrue(obs_seq.equals(exp_seq))

    def test_update_ids_default_behavior(self):
        # 3 seqs
        exp_sc = SequenceCollection([
            RNA('GAUUACA', id="1"),
            RNA('UUG', id="2"),
            RNA('U-----UGCC--', id="3")
        ])
        exp_id_map = {'1': 'r1', '2': 'r2', '3': 'r3'}
        obs_sc, obs_id_map = self.s2.update_ids()
        self._assert_sequence_collections_equal(obs_sc, exp_sc)
        self.assertEqual(obs_id_map, exp_id_map)

        # empty
        obs_sc, obs_id_map = self.empty.update_ids()
        self._assert_sequence_collections_equal(obs_sc, self.empty)
        self.assertEqual(obs_id_map, {})

    def test_update_ids_prefix(self):
        # 3 seqs
        exp_sc = SequenceCollection([
            RNA('GAUUACA', id="abc1"),
            RNA('UUG', id="abc2"),
            RNA('U-----UGCC--', id="abc3")
        ])
        exp_id_map = {'abc1': 'r1', 'abc2': 'r2', 'abc3': 'r3'}
        obs_sc, obs_id_map = self.s2.update_ids(prefix='abc')
        self._assert_sequence_collections_equal(obs_sc, exp_sc)
        self.assertEqual(obs_id_map, exp_id_map)

        # empty
        obs_sc, obs_id_map = self.empty.update_ids(prefix='abc')
        self._assert_sequence_collections_equal(obs_sc, self.empty)
        self.assertEqual(obs_id_map, {})

    def test_update_ids_fn_parameter(self):
        def append_42(ids):
            return [id_ + '-42' for id_ in ids]

        # 3 seqs
        exp_sc = SequenceCollection([
            RNA('GAUUACA', id="r1-42"),
            RNA('UUG', id="r2-42"),
            RNA('U-----UGCC--', id="r3-42")
        ])
        exp_id_map = {'r1-42': 'r1', 'r2-42': 'r2', 'r3-42': 'r3'}
        obs_sc, obs_id_map = self.s2.update_ids(fn=append_42)
        self._assert_sequence_collections_equal(obs_sc, exp_sc)
        self.assertEqual(obs_id_map, exp_id_map)

        # empty
        obs_sc, obs_id_map = self.empty.update_ids(fn=append_42)
        self._assert_sequence_collections_equal(obs_sc, self.empty)
        self.assertEqual(obs_id_map, {})

    def test_update_ids_ids_parameter(self):
        # 3 seqs
        exp_sc = SequenceCollection([
            RNA('GAUUACA', id="abc"),
            RNA('UUG', id="def"),
            RNA('U-----UGCC--', id="ghi")
        ])
        exp_id_map = {'abc': 'r1', 'def': 'r2', 'ghi': 'r3'}
        obs_sc, obs_id_map = self.s2.update_ids(ids=('abc', 'def', 'ghi'))
        self._assert_sequence_collections_equal(obs_sc, exp_sc)
        self.assertEqual(obs_id_map, exp_id_map)

        # empty
        obs_sc, obs_id_map = self.empty.update_ids(ids=[])
        self._assert_sequence_collections_equal(obs_sc, self.empty)
        self.assertEqual(obs_id_map, {})

    def test_update_ids_sequence_attributes_propagated(self):
        # 1 seq
        exp_sc = Alignment([
            DNA('ACGT', id="abc", description='desc', quality=range(4))
        ])
        exp_id_map = {'abc': 'seq1'}

        obj = Alignment([
            DNA('ACGT', id="seq1", description='desc', quality=range(4))
        ])

        obs_sc, obs_id_map = obj.update_ids(ids=('abc',))
        self._assert_sequence_collections_equal(obs_sc, exp_sc)
        self.assertEqual(obs_id_map, exp_id_map)

        # 2 seqs
        exp_sc = Alignment([
            DNA('ACGT', id="abc", description='desc1', quality=range(4)),
            DNA('TGCA', id="def", description='desc2', quality=range(4)[::-1])
        ])
        exp_id_map = {'abc': 'seq1', 'def': 'seq2'}

        obj = Alignment([
            DNA('ACGT', id="seq1", description='desc1', quality=(0, 1, 2, 3)),
            DNA('TGCA', id="seq2", description='desc2', quality=(3, 2, 1, 0))
        ])

        obs_sc, obs_id_map = obj.update_ids(ids=('abc', 'def'))
        self._assert_sequence_collections_equal(obs_sc, exp_sc)
        self.assertEqual(obs_id_map, exp_id_map)

    def test_update_ids_invalid_parameter_combos(self):
        with self.assertRaisesRegexp(SequenceCollectionError, 'ids and fn'):
            self.s1.update_ids(fn=lambda e: e, ids=['foo', 'bar'])

        with self.assertRaisesRegexp(SequenceCollectionError, 'prefix'):
            self.s1.update_ids(ids=['foo', 'bar'], prefix='abc')

        with self.assertRaisesRegexp(SequenceCollectionError, 'prefix'):
            self.s1.update_ids(fn=lambda e: e, prefix='abc')

    def test_update_ids_invalid_ids(self):
        # incorrect number of new ids
        with self.assertRaisesRegexp(SequenceCollectionError, '3 != 2'):
            self.s1.update_ids(ids=['foo', 'bar', 'baz'])
        with self.assertRaisesRegexp(SequenceCollectionError, '4 != 2'):
            self.s1.update_ids(fn=lambda e: ['foo', 'bar', 'baz', 'abc'])

        # duplicates
        with self.assertRaisesRegexp(SequenceCollectionError, 'foo'):
            self.s2.update_ids(ids=['foo', 'bar', 'foo'])
        with self.assertRaisesRegexp(SequenceCollectionError, 'bar'):
            self.s2.update_ids(fn=lambda e: ['foo', 'bar', 'bar'])

    def test_int_map(self):
        expected1 = {"1": self.d1, "2": self.d2}
        expected2 = {"1": "d1", "2": "d2"}
        obs = npt.assert_warns(DeprecationWarning, self.s1.int_map)
        self.assertEqual(obs, (expected1, expected2))

        expected1 = {"h-1": self.d1, "h-2": self.d2}
        expected2 = {"h-1": "d1", "h-2": "d2"}
        obs = npt.assert_warns(DeprecationWarning, self.s1.int_map,
                               prefix='h-')
        self.assertEqual(obs, (expected1, expected2))

    def test_is_empty(self):
        self.assertFalse(self.s1.is_empty())
        self.assertFalse(self.s2.is_empty())
        self.assertFalse(self.s3.is_empty())

        self.assertTrue(self.empty.is_empty())

    def test_is_valid(self):
        self.assertTrue(self.s1.is_valid())
        self.assertTrue(self.s2.is_valid())
        self.assertTrue(self.s3.is_valid())
        self.assertTrue(self.empty.is_valid())

        self.assertFalse(self.invalid_s1.is_valid())

    def test_iteritems(self):
        self.assertEqual(list(self.s1.iteritems()),
                         [(s.id, s) for s in self.s1])

    def test_lower(self):
        self.assertEqual(self.s1.lower(), self.s1_lower)

    def test_sequence_count(self):
        self.assertEqual(self.s1.sequence_count(), 2)
        self.assertEqual(self.s2.sequence_count(), 3)
        self.assertEqual(self.s3.sequence_count(), 5)
        self.assertEqual(self.empty.sequence_count(), 0)

    def test_sequence_lengths(self):
        self.assertEqual(self.s1.sequence_lengths(), [7, 3])
        self.assertEqual(self.s2.sequence_lengths(), [7, 3, 12])
        self.assertEqual(self.s3.sequence_lengths(), [7, 3, 7, 3, 12])
        self.assertEqual(self.empty.sequence_lengths(), [])

    def test_to_fasta(self):
        exp1 = ">d1\nGATTACA\n>d2\nTTG\n"
        self.assertEqual(self.s1.to_fasta(), exp1)
        exp2 = ">r1\nGAUUACA\n>r2\nUUG\n>r3\nU-----UGCC--\n"
        self.assertEqual(self.s2.to_fasta(), exp2)

    def test_toFasta(self):
        exp = ">d1\nGATTACA\n>d2\nTTG\n"
        obs = npt.assert_warns(DeprecationWarning, self.s1.toFasta)
        self.assertEqual(obs, exp)

    def test_upper(self):
        self.assertEqual(self.s1_lower.upper(), self.s1)


class AlignmentTests(TestCase):

    def setUp(self):
        self.d1 = DNASequence('..ACC-GTTGG..', id="d1")
        self.d2 = DNASequence('TTACCGGT-GGCC', id="d2")
        self.d3 = DNASequence('.-ACC-GTTGC--', id="d3")

        self.r1 = RNASequence('UUAU-', id="r1")
        self.r2 = RNASequence('ACGUU', id="r2")

        self.seqs1 = [self.d1, self.d2, self.d3]
        self.seqs2 = [self.r1, self.r2]

        self.seqs1_t = [('d1', '..ACC-GTTGG..'), ('d2', 'TTACCGGT-GGCC'),
                        ('d3', '.-ACC-GTTGC--')]
        self.seqs2_t = [('r1', 'UUAU-'), ('r2', 'ACGUU')]

        self.a1 = Alignment(self.seqs1)
        self.a2 = Alignment(self.seqs2)
        self.a3 = Alignment(self.seqs2, score=42.0,
                            start_end_positions=[(0, 3), (5, 9)])
        self.a4 = Alignment(self.seqs2, score=-42.0,
                            start_end_positions=[(1, 4), (6, 10)])

        # no sequences
        self.empty = Alignment([])

        # sequences, but no positions
        self.no_positions = Alignment([RNA('', id='a'), RNA('', id='b')])

    def test_degap(self):
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
        expected = [[0, 6. / 13, 4. / 13],
                    [6. / 13, 0, 7. / 13],
                    [4. / 13, 7. / 13, 0]]
        expected = DistanceMatrix(expected, ['d1', 'd2', 'd3'])
        actual = self.a1.distances()
        self.assertEqual(actual, expected)

        # alt distance function provided
        def dumb_distance(s1, s2):
            return 42.
        expected = [[0, 42., 42.],
                    [42., 0, 42.],
                    [42., 42., 0]]
        expected = DistanceMatrix(expected, ['d1', 'd2', 'd3'])
        actual = self.a1.distances(dumb_distance)
        self.assertEqual(actual, expected)

    def test_score(self):
        self.assertEqual(self.a3.score(), 42.0)
        self.assertEqual(self.a4.score(), -42.0)

    def test_start_end_positions(self):
        self.assertEqual(self.a3.start_end_positions(), [(0, 3), (5, 9)])
        self.assertEqual(self.a4.start_end_positions(), [(1, 4), (6, 10)])

    def test_subalignment(self):
        # keep seqs by ids
        actual = self.a1.subalignment(seqs_to_keep=['d1', 'd3'])
        expected = Alignment([self.d1, self.d3])
        self.assertEqual(actual, expected)

        # keep seqs by indices
        actual = self.a1.subalignment(seqs_to_keep=[0, 2])
        expected = Alignment([self.d1, self.d3])
        self.assertEqual(actual, expected)

        # keep seqs by ids (invert)
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
        d1 = DNASequence('.AC', id="d1")
        d2 = DNASequence('TAC', id="d2")
        d3 = DNASequence('.AC', id="d3")
        expected = Alignment([d1, d2, d3])
        self.assertEqual(actual, expected)

        # keep positions (invert)
        actual = self.a1.subalignment(positions_to_keep=[0, 2, 3],
                                      invert_positions_to_keep=True)
        d1 = DNASequence('.C-GTTGG..', id="d1")
        d2 = DNASequence('TCGGT-GGCC', id="d2")
        d3 = DNASequence('-C-GTTGC--', id="d3")
        expected = Alignment([d1, d2, d3])
        self.assertEqual(actual, expected)

        # keep seqs and positions
        actual = self.a1.subalignment(seqs_to_keep=[0, 2],
                                      positions_to_keep=[0, 2, 3])
        d1 = DNASequence('.AC', id="d1")
        d3 = DNASequence('.AC', id="d3")
        expected = Alignment([d1, d3])
        self.assertEqual(actual, expected)

        # keep seqs and positions (invert)
        actual = self.a1.subalignment(seqs_to_keep=[0, 2],
                                      positions_to_keep=[0, 2, 3],
                                      invert_seqs_to_keep=True,
                                      invert_positions_to_keep=True)
        d2 = DNASequence('TCGGT-GGCC', id="d2")
        expected = Alignment([d2])
        self.assertEqual(actual, expected)

    def test_subalignment_filter_out_everything(self):
        exp = Alignment([])

        # no sequences
        obs = self.a1.subalignment(seqs_to_keep=None, invert_seqs_to_keep=True)
        self.assertEqual(obs, exp)

        # no positions
        obs = self.a1.subalignment(positions_to_keep=None,
                                   invert_positions_to_keep=True)
        self.assertEqual(obs, exp)

    def test_init_not_equal_lengths(self):
        invalid_seqs = [self.d1, self.d2, self.d3,
                        DNASequence('.-ACC-GTGC--', id="i2")]
        self.assertRaises(AlignmentError, Alignment,
                          invalid_seqs)

    def test_init_equal_lengths(self):
        seqs = [self.d1, self.d2, self.d3]
        Alignment(seqs)

    def test_init_validate(self):
        Alignment(self.seqs1, validate=True)

        # invalid DNA character
        invalid_seqs1 = [self.d1, self.d2, self.d3,
                         DNASequence('.-ACC-GTXGC--', id="i1")]
        self.assertRaises(SequenceCollectionError, Alignment,
                          invalid_seqs1, validate=True)

    def test_iter_positions(self):
        actual = list(self.a2.iter_positions())
        expected = [[RNASequence(j) for j in i] for i in
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
        d1 = DNASequence('TTT', id="d1")
        d2 = DNASequence('TT-', id="d2")
        d3 = DNASequence('TC-', id="d3")
        a1 = Alignment([d1, d2, d3])
        self.assertTrue(a1.majority_consensus().equals(DNASequence('TT-')))

        d1 = DNASequence('T', id="d1")
        d2 = DNASequence('A', id="d2")
        a1 = Alignment([d1, d2])
        self.assertTrue(a1.majority_consensus() in
                        [DNASequence('T'), DNASequence('A')])

        self.assertEqual(self.empty.majority_consensus(), '')

    def test_majority_consensus_constructor(self):
        d1 = DNASequence('TTT', id="d1")
        d2 = DNASequence('TT-', id="d2")
        d3 = DNASequence('TC-', id="d3")
        a1 = Alignment([d1, d2, d3])

        obs = npt.assert_warns(DeprecationWarning, a1.majority_consensus,
                               constructor=str)
        self.assertEqual(obs, 'TT-')

    def test_omit_gap_positions(self):
        expected = self.a2
        self.assertEqual(self.a2.omit_gap_positions(1.0), expected)
        self.assertEqual(self.a2.omit_gap_positions(0.51), expected)

        r1 = RNASequence('UUAU', id="r1")
        r2 = RNASequence('ACGU', id="r2")
        expected = Alignment([r1, r2])
        self.assertEqual(self.a2.omit_gap_positions(0.49), expected)

        r1 = RNASequence('UUAU', id="r1")
        r2 = RNASequence('ACGU', id="r2")
        expected = Alignment([r1, r2])
        self.assertEqual(self.a2.omit_gap_positions(0.0), expected)

        self.assertEqual(self.empty.omit_gap_positions(0.0), self.empty)
        self.assertEqual(self.empty.omit_gap_positions(0.49), self.empty)
        self.assertEqual(self.empty.omit_gap_positions(1.0), self.empty)

        # Test to ensure floating point precision bug isn't present. See the
        # tests for Alignment.position_frequencies for more details.
        seqs = []
        for i in range(33):
            seqs.append(DNA('-.', id=str(i)))
        aln = Alignment(seqs)
        self.assertEqual(aln.omit_gap_positions(1 - np.finfo(float).eps),
                         Alignment([DNA('', id=str(i)) for i in range(33)]))

    def test_omit_gap_sequences(self):
        expected = self.a2
        self.assertEqual(self.a2.omit_gap_sequences(1.0), expected)
        self.assertEqual(self.a2.omit_gap_sequences(0.20), expected)

        expected = Alignment([self.r2])
        self.assertEqual(self.a2.omit_gap_sequences(0.19), expected)

        self.assertEqual(self.empty.omit_gap_sequences(0.0), self.empty)
        self.assertEqual(self.empty.omit_gap_sequences(0.2), self.empty)
        self.assertEqual(self.empty.omit_gap_sequences(1.0), self.empty)

        # Test to ensure floating point precision bug isn't present. See the
        # tests for Alignment.position_frequencies for more details.
        aln = Alignment([DNA('.' * 33, id='abc'), DNA('-' * 33, id='def')])
        self.assertEqual(aln.omit_gap_sequences(1 - np.finfo(float).eps),
                         Alignment([]))

    def test_position_counters(self):
        self.assertEqual(self.empty.position_counters(), [])

        self.assertEqual(self.no_positions.position_counters(), [])

        expected = [Counter({'U': 1, 'A': 1}),
                    Counter({'U': 1, 'C': 1}),
                    Counter({'A': 1, 'G': 1}),
                    Counter({'U': 2}),
                    Counter({'-': 1, 'U': 1})]
        self.assertEqual(self.a2.position_counters(), expected)

    def test_position_frequencies(self):
        self.assertEqual(self.empty.position_frequencies(), [])

        self.assertEqual(self.no_positions.position_frequencies(), [])

        expected = [defaultdict(float, {'U': 0.5, 'A': 0.5}),
                    defaultdict(float, {'U': 0.5, 'C': 0.5}),
                    defaultdict(float, {'A': 0.5, 'G': 0.5}),
                    defaultdict(float, {'U': 1.0}),
                    defaultdict(float, {'-': 0.5, 'U': 0.5})]
        self.assertEqual(self.a2.position_frequencies(), expected)

    def test_position_frequencies_floating_point_precision(self):
        # Test that a position with no variation yields a frequency of exactly
        # 1.0. Note that it is important to use self.assertEqual here instead
        # of self.assertAlmostEqual because we want to test for exactly 1.0. A
        # previous implementation of Alignment.position_frequencies added
        # (1 / sequence_count) for each occurrence of a character in a position
        # to compute the frequencies (see
        # https://github.com/biocore/scikit-bio/issues/801). In certain cases,
        # this yielded a frequency slightly less than 1.0 due to roundoff
        # error. The test case here uses an alignment of 10 sequences with no
        # variation at a position. This test case exposes the roundoff error
        # present in the previous implementation because 1/10 added 10 times
        # yields a number slightly less than 1.0. This occurs because 1/10
        # cannot be represented exactly as a floating point number.
        seqs = []
        for i in range(10):
            seqs.append(DNA('A', id=str(i)))
        aln = Alignment(seqs)
        self.assertEqual(aln.position_frequencies(),
                         [defaultdict(float, {'A': 1.0})])

    def test_position_entropies(self):
        # tested by calculating values as described in this post:
        #  http://stackoverflow.com/a/15476958/3424666
        expected = [0.69314, 0.69314, 0.69314, 0.0, np.nan]
        np.testing.assert_almost_equal(self.a2.position_entropies(),
                                       expected, 5)

        expected = [1.0, 1.0, 1.0, 0.0, np.nan]
        np.testing.assert_almost_equal(self.a2.position_entropies(base=2),
                                       expected, 5)

        np.testing.assert_almost_equal(self.empty.position_entropies(base=2),
                                       [])

    def test_k_word_frequencies(self):
        expected = [defaultdict(float, {'U': 3 / 5, 'A': 1 / 5, '-': 1 / 5}),
                    defaultdict(float, {'A': 1 / 5, 'C': 1 / 5, 'G': 1 / 5,
                                        'U': 2 / 5})]
        actual = self.a2.k_word_frequencies(k=1)
        for a, e in zip(actual, expected):
            self.assertEqual(sorted(a), sorted(e), 5)
            np.testing.assert_almost_equal(sorted(a.values()),
                                           sorted(e.values()), 5)

    def test_sequence_length(self):
        self.assertEqual(self.a1.sequence_length(), 13)
        self.assertEqual(self.a2.sequence_length(), 5)
        self.assertEqual(self.empty.sequence_length(), 0)

    def test_to_phylip(self):
        d1 = DNASequence('..ACC-GTTGG..', id="d1")
        d2 = DNASequence('TTACCGGT-GGCC', id="d2")
        d3 = DNASequence('.-ACC-GTTGC--', id="d3")
        a = Alignment([d1, d2, d3])

        phylip_str, id_map = npt.assert_warns(DeprecationWarning, a.to_phylip,
                                              map_labels=False)
        self.assertEqual(id_map, {'d1': 'd1',
                                  'd3': 'd3',
                                  'd2': 'd2'})
        expected = "\n".join(["3 13",
                              "d1 ..ACC-GTTGG..",
                              "d2 TTACCGGT-GGCC",
                              "d3 .-ACC-GTTGC--"])
        self.assertEqual(phylip_str, expected)

    def test_to_phylip_map_labels(self):
        d1 = DNASequence('..ACC-GTTGG..', id="d1")
        d2 = DNASequence('TTACCGGT-GGCC', id="d2")
        d3 = DNASequence('.-ACC-GTTGC--', id="d3")
        a = Alignment([d1, d2, d3])

        phylip_str, id_map = npt.assert_warns(DeprecationWarning, a.to_phylip,
                                              map_labels=True,
                                              label_prefix="s")
        self.assertEqual(id_map, {'s1': 'd1',
                                  's3': 'd3',
                                  's2': 'd2'})
        expected = "\n".join(["3 13",
                              "s1 ..ACC-GTTGG..",
                              "s2 TTACCGGT-GGCC",
                              "s3 .-ACC-GTTGC--"])
        self.assertEqual(phylip_str, expected)

    def test_to_phylip_no_sequences(self):
        with self.assertRaises(SequenceCollectionError):
            npt.assert_warns(DeprecationWarning, Alignment([]).to_phylip)

    def test_to_phylip_no_positions(self):
        d1 = DNASequence('', id="d1")
        d2 = DNASequence('', id="d2")
        a = Alignment([d1, d2])

        with self.assertRaises(SequenceCollectionError):
            npt.assert_warns(DeprecationWarning, a.to_phylip)

    def test_heatmap_with_defaults(self):
        values, sequences, a1 = self.heatmap_set_values()
        fig = a1.heatmap(values)
        self.heatmap_basic_sanity(fig, ['A', 'A', 'C', 'C', 'C', 'G', 'T'],
                                  ['seq1', 'seq2'], ['Minimum', 'Median',
                                                     'Maximum'])

    def test_heatmap_with_custom(self):
        sequences = [DNA('AGTCGGT', id="seq1"),
                     DNA('CAACGGA', id="seq2"),
                     DNA('AACCTCT', id="seq3"),
                     DNA('TACTCGT', id="seq4")]
        a1 = Alignment(sequences)
        values = {'A': 0.61, 'C': 1.07, 'T': 0.05, 'G': 0.07}
        clabels = ['a', 'b', 'c']
        fig = a1.heatmap(values, fig_size=(15, 10), cmap='Blues',
                         legend_labels=clabels,
                         sequence_order=('seq4', 'seq3', 'seq2', 'seq1'))
        self.heatmap_basic_sanity(fig, ['A', 'A', 'C', 'C', 'G', 'G', 'T'],
                                  ['seq4', 'seq1'], clabels)
        self.assertEqual(fig.get_figwidth(), 15.0)
        self.assertEqual(fig.get_figheight(), 10.0)

    def test_heatmap_raises(self):
        values, sequences, a1 = self.heatmap_set_values()
        with self.assertRaises(ValueError):
            a1.heatmap(values, legend_labels=['a', 'b', 'c', 'd'])

    def test_heatmap_errors(self):
        sequences = [DNA('AGTCGGT', id="seq1"),
                     DNA('CAACGGA', id="seq2"),
                     DNA('AACCTCT', id="seq3"),
                     DNA('TACTCGT', id="seq4")]
        a1 = Alignment(sequences)
        values = {'A': 0.61, 'C': 1.07, 'T': 0.05, 'G': 0.07}
        clabels = ['a', 'b', 'c']
        with self.assertRaises(ValueError):
            a1.heatmap(values, fig_size=(15, 10), cmap='Blues',
                       legend_labels=clabels,
                       sequence_order=('seq1', 'seq2', 'seq3', 'seq3'))
        with self.assertRaises(ValueError):
            a1.heatmap(values, fig_size=(15, 10), cmap='Blues',
                       legend_labels=clabels,
                       sequence_order=('seq1', 'seq2', 'seq4', 'seq3',
                                       'seq5'))

    def heatmap_set_values(self):
        sequences = [DNA('AACCCGT', id="seq1"),
                     DNA('AACCGGT', id="seq2")]
        a1 = Alignment(sequences)
        values = {'A': 0.61, 'C': 1.07, 'T': 0.05, 'G': 0.07}
        return(values, sequences, a1)

    def heatmap_basic_sanity(self, fig, xt, yt, clabels):
        axes = fig.get_axes()
        self.assertEqual(len(axes), 2)
        ax = axes[0]
        axc = axes[1]
        xticks = []
        for tick in ax.get_xticklabels():
            xticks.append(tick.get_text())
        self.assertEqual(xticks, xt)
        yticks = []
        for tick in ax.get_yticklabels():
            yticks.append(tick.get_text())
        self.assertEqual(yticks, yt)
        cticks = []
        for tick in axc.get_xticklabels():
            cticks.append(tick.get_text())
        self.assertEqual(clabels, cticks)

    def test_validate_lengths(self):
        self.assertTrue(self.a1._validate_lengths())
        self.assertTrue(self.a2._validate_lengths())
        self.assertTrue(self.empty._validate_lengths())

        self.assertTrue(Alignment([
            DNASequence('TTT', id="d1")])._validate_lengths())


class StockholmAlignmentTests(TestCase):
    def setUp(self):
        self.seqs = [DNASequence("ACC-G-GGTA", id="seq1"),
                     DNASequence("TCC-G-GGCA", id="seq2")]
        self.GF = OrderedDict([
            ("AC", "RF00360"),
            ("BM", ["cmbuild  -F CM SEED",
                    "cmsearch  -Z 274931 -E 1000000"]),
            ("SQ", "9"),
            ("RT", ["TITLE1",  "TITLE2"]),
            ("RN", ["[1]", "[2]"]),
            ("RA", ["Auth1;", "Auth2;"]),
            ("RL", ["J Mol Biol", "Cell"]),
            ("RM", ["11469857", "12007400"]),
            ('RN', ['[1]', '[2]'])
        ])
        self.GS = {"AC": OrderedDict([("seq1", "111"), ("seq2", "222")])}
        self.GR = {"SS": OrderedDict([("seq1", "1110101111"),
                                      ("seq2", "0110101110")])}
        self.GC = {"SS_cons": "(((....)))"}
        self.st = StockholmAlignment(self.seqs, gc=self.GC, gf=self.GF,
                                     gs=self.GS, gr=self.GR)

    def test_retrieve_metadata(self):
        self.assertEqual(self.st.gc, self.GC)
        self.assertEqual(self.st.gf, self.GF)
        self.assertEqual(self.st.gs, self.GS)
        self.assertEqual(self.st.gr, self.GR)

    def test_from_file_alignment(self):
        # test that a basic stockholm file with interleaved alignment can be
        # parsed
        sto = StringIO("# STOCKHOLM 1.0\n"
                       "seq1      ACC-G\n"
                       "seq2      TCC-G\n\n"
                       "seq1      -GGTA\n"
                       "seq2      -GGCA\n//")
        obs_sto = next(StockholmAlignment.from_file(sto, DNA))
        exp_sto = StockholmAlignment(self.seqs)
        self.assertEqual(obs_sto, exp_sto)

    def test_from_file_GF(self):
        # remove rn line to make sure auto-added
        self.GF.pop("RN")
        sto = StringIO("# STOCKHOLM 1.0\n#=GF RN [1]\n#=GF RM 11469857\n"
                       "#=GF RT TITLE1\n#=GF RA Auth1;\n#=GF RL J Mol Biol\n"
                       "#=GF RN [2]\n#=GF RM 12007400\n#=GF RT TITLE2\n"
                       "#=GF RA Auth2;\n#=GF RL Cell\n#=GF AC RF00360\n"
                       "#=GF BM cmbuild  -F CM SEED\n"
                       "#=GF BM cmsearch  -Z 274931 -E 1000000\n#=GF SQ 9\n"
                       "seq1         ACC-G-GGTA\nseq2         TCC-G-GGCA\n//")
        obs_sto = next(StockholmAlignment.from_file(sto, DNA))
        exp_sto = StockholmAlignment(self.seqs, self.GF, {}, {}, {})
        self.assertEqual(obs_sto, exp_sto)

    def test_from_file_GC(self):
        sto = StringIO("# STOCKHOLM 1.0\n"
                       "seq1         ACC-G-GGTA\nseq2         TCC-G-GGCA\n"
                       "#=GC SS_cons (((....)))\n//")
        obs_sto = next(StockholmAlignment.from_file(sto, DNA))
        exp_sto = StockholmAlignment(self.seqs, {}, {}, {}, self.GC)
        self.assertEqual(obs_sto, exp_sto)

    def test_from_file_GS(self):
        sto = StringIO("# STOCKHOLM 1.0\n#=GS seq2 AC 222\n#=GS seq1 AC 111\n"
                       "seq1          ACC-G-GGTA\n"
                       "seq2          TCC-G-GGCA\n//")
        obs_sto = next(StockholmAlignment.from_file(sto, DNA))
        exp_sto = StockholmAlignment(self.seqs, {}, self.GS, {}, {})
        self.assertEqual(obs_sto, exp_sto)

    def test_from_file_GR(self):
        sto = StringIO("# STOCKHOLM 1.0\nseq1          ACC-G\n"
                       "#=GR seq1 SS  11101\nseq2          TCC-G\n"
                       "#=GR seq2 SS  01101\n\nseq1          -GGTA\n"
                       "#=GR seq1 SS  01111\nseq2          -GGCA\n"
                       "#=GR seq2 SS  01110\n//")
        obs_sto = next(StockholmAlignment.from_file(sto, DNA))
        exp_sto = StockholmAlignment(self.seqs, {}, {}, self.GR, {})
        self.assertEqual(obs_sto, exp_sto)

    def test_from_file_multi(self):
        sto = StringIO("# STOCKHOLM 1.0\n#=GS seq2 AC 222\n#=GS seq1 AC 111\n"
                       "seq1          ACC-G-GGTA\n"
                       "seq2          TCC-G-GGCA\n//\n"
                       "# STOCKHOLM 1.0\nseq1          ACC-G-GGTA\n"
                       "#=GR seq1 SS  1110101111\nseq2          TCC-G-GGCA\n"
                       "#=GR seq2 SS  0110101110\n//")
        obs_sto = StockholmAlignment.from_file(sto, DNA)
        count = 0
        for obs in obs_sto:
            if count == 0:
                exp_sto = StockholmAlignment(self.seqs, {}, self.GS, {}, {})
                self.assertEqual(obs, exp_sto)
            elif count == 1:
                exp_sto = StockholmAlignment(self.seqs, {}, {}, self.GR, {})
                self.assertEqual(obs, exp_sto)
            else:
                raise AssertionError("More than 2 sto alignments parsed!")
            count += 1

    def test_parse_gf_multiline_nh(self):
        sto = ["#=GF TN MULTILINE TREE",
               "#=GF NH THIS IS FIRST", "#=GF NH THIS IS SECOND",
               "#=GF AC 1283394"]
        exp = {'TN': 'MULTILINE TREE',
               'NH': 'THIS IS FIRST THIS IS SECOND',
               'AC': '1283394'}
        self.assertEqual(self.st._parse_gf_info(sto), exp)

    def test_parse_gf_multiline_cc(self):
        sto = ["#=GF CC THIS IS FIRST", "#=GF CC THIS IS SECOND"]
        exp = {'CC': 'THIS IS FIRST THIS IS SECOND'}
        self.assertEqual(self.st._parse_gf_info(sto), exp)

    def test_parse_gf_info_nongf(self):
        sto = ["#=GF AC BLAAAAAAAHHH", "#=GC HUH THIS SHOULD NOT BE HERE"]
        with self.assertRaises(StockholmParseError):
            self.st._parse_gf_info(sto)

    def test_parse_gf_info_malformed(self):
        # too short of a line
        sto = ["#=GF AC", "#=GF"]
        with self.assertRaises(StockholmParseError):
            self.st._parse_gf_info(sto)

    def test_parse_gc_info_nongf(self):
        sto = ["#=GC AC BLAAAAAAAHHH", "#=GF HUH THIS SHOULD NOT BE HERE"]
        with self.assertRaises(StockholmParseError):
            self.st._parse_gf_info(sto)

    def test_parse_gc_info_strict_len(self):
        sto = ["#=GC SS_cons (((..)))"]
        with self.assertRaises(StockholmParseError):
            self.st._parse_gc_info(sto, seqlen=20, strict=True)

    def test_parse_gc_info_strict_duplicate(self):
        sto = ["#=GC SS_cons (((..)))", "#=GC SS_cons (((..)))"]
        with self.assertRaises(StockholmParseError):
            self.st._parse_gc_info(sto, seqlen=8, strict=True)

    def test_parse_gc_info_malformed(self):
        # too short of a line
        sto = ["#=GC AC BLAAAAAAAHHH", "#=GC"]
        with self.assertRaises(StockholmParseError):
            self.st._parse_gc_info(sto)

    def test_parse_gs_gr_info_mixed(self):
        sto = ["#=GS seq1 AC BLAAA", "#=GR seq2 HUH THIS SHOULD NOT BE HERE"]
        with self.assertRaises(StockholmParseError):
            self.st._parse_gs_gr_info(sto)

    def test_parse_gs_gr_info_malformed(self):
        # too short of a line
        sto = ["#=GS AC BLAAAAAAAHHH", "#=GS"]
        with self.assertRaises(StockholmParseError):
            self.st._parse_gs_gr_info(sto)

    def test_parse_gs_gr_info_strict(self):
        sto = ["#=GR seq1 SS  10101111", "#=GR seq2 SS  01101"]
        with self.assertRaises(StockholmParseError):
            self.st._parse_gs_gr_info(sto, seqlen=20, strict=True)

    def test_str(self):
        st = StockholmAlignment(self.seqs, gc=self.GC, gf=self.GF, gs=self.GS,
                                gr=self.GR)
        obs = str(st)
        exp = ('# STOCKHOLM 1.0\n'
               '#=GF AC RF00360\n'
               '#=GF BM cmbuild  -F CM SEED\n'
               '#=GF BM cmsearch  -Z 274931 -E 1000000\n'
               '#=GF SQ 9\n'
               '#=GF RN [1]\n'
               '#=GF RM 11469857\n'
               '#=GF RT TITLE1\n'
               '#=GF RA Auth1;\n'
               '#=GF RL J Mol Biol\n'
               '#=GF RN [2]\n'
               '#=GF RM 12007400\n'
               '#=GF RT TITLE2\n'
               '#=GF RA Auth2;\n'
               '#=GF RL Cell\n'
               '#=GS seq1 AC 111\n'
               '#=GS seq2 AC 222\n'
               'seq1          ACC-G-GGTA\n'
               '#=GR seq1 SS  1110101111\n'
               'seq2          TCC-G-GGCA\n'
               '#=GR seq2 SS  0110101110\n'
               '#=GC SS_cons  (((....)))\n//')
        self.assertEqual(obs, exp)

    def test_to_file(self):
        st = StockholmAlignment(self.seqs, gc=self.GC, gf=self.GF, gs=self.GS,
                                gr=self.GR)

        with tempfile.NamedTemporaryFile('r+') as temp_file:
            st.to_file(temp_file)
            temp_file.flush()
            temp_file.seek(0)
            obs = temp_file.read()
            exp = ('# STOCKHOLM 1.0\n'
                   '#=GF AC RF00360\n'
                   '#=GF BM cmbuild  -F CM SEED\n'
                   '#=GF BM cmsearch  -Z 274931 -E 1000000\n'
                   '#=GF SQ 9\n'
                   '#=GF RN [1]\n'
                   '#=GF RM 11469857\n'
                   '#=GF RT TITLE1\n'
                   '#=GF RA Auth1;\n'
                   '#=GF RL J Mol Biol\n'
                   '#=GF RN [2]\n'
                   '#=GF RM 12007400\n'
                   '#=GF RT TITLE2\n'
                   '#=GF RA Auth2;\n'
                   '#=GF RL Cell\n'
                   '#=GS seq1 AC 111\n'
                   '#=GS seq2 AC 222\n'
                   'seq1          ACC-G-GGTA\n'
                   '#=GR seq1 SS  1110101111\n'
                   'seq2          TCC-G-GGCA\n'
                   '#=GR seq2 SS  0110101110\n'
                   '#=GC SS_cons  (((....)))\n//')
        self.assertEqual(obs, exp)

    def test_str_gc(self):
        st = StockholmAlignment(self.seqs, gc=self.GC, gf=None, gs=None,
                                gr=None)
        obs = str(st)
        exp = ("# STOCKHOLM 1.0\nseq1          ACC-G-GGTA\n"
               "seq2          TCC-G-GGCA\n"
               "#=GC SS_cons  (((....)))\n//")
        self.assertEqual(obs, exp)

    def test_str_gf(self):
        st = StockholmAlignment(self.seqs, gc=None, gf=self.GF, gs=None,
                                gr=None)
        obs = str(st)
        exp = ('# STOCKHOLM 1.0\n'
               '#=GF AC RF00360\n'
               '#=GF BM cmbuild  -F CM SEED\n'
               '#=GF BM cmsearch  -Z 274931 -E 1000000\n'
               '#=GF SQ 9\n'
               '#=GF RN [1]\n'
               '#=GF RM 11469857\n'
               '#=GF RT TITLE1\n'
               '#=GF RA Auth1;\n'
               '#=GF RL J Mol Biol\n'
               '#=GF RN [2]\n'
               '#=GF RM 12007400\n'
               '#=GF RT TITLE2\n'
               '#=GF RA Auth2;\n'
               '#=GF RL Cell\n'
               'seq1          ACC-G-GGTA\n'
               'seq2          TCC-G-GGCA\n//')
        self.assertEqual(obs, exp)

    def test_str_gs(self):
        st = StockholmAlignment(self.seqs, gc=None, gf=None, gs=self.GS,
                                gr=None)
        obs = str(st)
        exp = ('# STOCKHOLM 1.0\n'
               '#=GS seq1 AC 111\n'
               '#=GS seq2 AC 222\n'
               'seq1          ACC-G-GGTA\n'
               'seq2          TCC-G-GGCA\n//')
        self.assertEqual(obs, exp)

    def test_str_gr(self):
        st = StockholmAlignment(self.seqs, gc=None, gf=None, gs=None,
                                gr=self.GR)
        obs = str(st)
        exp = ("# STOCKHOLM 1.0\nseq1          ACC-G-GGTA\n"
               "#=GR seq1 SS  1110101111\nseq2          TCC-G-GGCA\n"
               "#=GR seq2 SS  0110101110\n//")
        self.assertEqual(obs, exp)

    def test_str_trees(self):
        GF = OrderedDict({"NH": ["IMATREE", "IMATREETOO"],
                          "TN": ["Tree2", "Tree1"]})
        st = StockholmAlignment(self.seqs, gc=None, gf=GF, gs=None,
                                gr=None)
        obs = str(st)
        exp = ("# STOCKHOLM 1.0\n#=GF TN Tree2\n#=GF NH IMATREE\n#=GF TN Tree1"
               "\n#=GF NH IMATREETOO\nseq1          ACC-G-GGTA\n"
               "seq2          TCC-G-GGCA\n//")

        self.assertEqual(obs, exp)


if __name__ == "__main__":
    main()
