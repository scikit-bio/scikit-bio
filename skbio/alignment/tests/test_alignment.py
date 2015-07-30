# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import six

from unittest import TestCase, main
from collections import Counter, defaultdict

import numpy as np
from scipy.spatial.distance import hamming

from skbio import (Sequence, DNA, RNA,
                   DistanceMatrix, Alignment, SequenceCollection)
from skbio.alignment import (SequenceCollectionError, AlignmentError)


class SequenceCollectionTests(TestCase):
    def setUp(self):
        self.d1 = DNA('GATTACA', metadata={'id': "d1"})
        self.d2 = DNA('TTG', metadata={'id': "d2"})
        self.d3 = DNA('GTATACA', metadata={'id': "d3"})
        self.r1 = RNA('GAUUACA', metadata={'id': "r1"})
        self.r2 = RNA('UUG', metadata={'id': "r2"})
        self.r3 = RNA('U-----UGCC--', metadata={'id': "r3"})

        self.seqs1 = [self.d1, self.d2]
        self.seqs2 = [self.r1, self.r2, self.r3]
        self.seqs3 = self.seqs1 + self.seqs2
        self.seqs4 = [self.d1, self.d3]

        self.s1 = SequenceCollection(self.seqs1)
        self.s2 = SequenceCollection(self.seqs2)
        self.s3 = SequenceCollection(self.seqs3)
        self.s4 = SequenceCollection(self.seqs4)
        self.empty = SequenceCollection([])

    def test_init(self):
        SequenceCollection(self.seqs1)
        SequenceCollection(self.seqs2)
        SequenceCollection(self.seqs3)
        SequenceCollection([])

    def test_init_fail(self):
        # sequences with overlapping ids
        s1 = [self.d1, self.d1]
        self.assertRaises(SequenceCollectionError, SequenceCollection, s1)

    def test_init_fail_no_id(self):
        seq = Sequence('ACGTACGT')
        with six.assertRaisesRegex(self, SequenceCollectionError,
                                   "'id' must be included in the sequence "
                                   "metadata"):
            SequenceCollection([seq])

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

    def test_kmer_frequencies(self):
        expected1 = Counter({'GAT': 1, 'TAC': 1})
        expected2 = Counter({'TTG': 1})
        self.assertEqual(
            self.s1.kmer_frequencies(k=3, overlap=False, relative=False),
            [expected1, expected2])

        expected1 = defaultdict(float)
        expected1['A'] = 3 / 7.
        expected1['C'] = 1 / 7.
        expected1['G'] = 1 / 7.
        expected1['T'] = 2 / 7.
        expected2 = defaultdict(float)
        expected2['G'] = 1 / 3.
        expected2['T'] = 2 / 3.
        self.assertEqual(self.s1.kmer_frequencies(k=1, relative=True),
                         [expected1, expected2])

        expected1 = defaultdict(float)
        expected1['GAT'] = 1 / 2.
        expected1['TAC'] = 1 / 2.
        expected2 = defaultdict(float)
        expected2['TTG'] = 1 / 1.
        self.assertEqual(
            self.s1.kmer_frequencies(k=3, overlap=False, relative=True),
            [expected1, expected2])

        self.assertEqual(self.empty.kmer_frequencies(k=1, relative=True), [])

        # Test to ensure floating point precision bug isn't present. See the
        # tests for Sequence.kmer_frequencies for more details.
        sc = SequenceCollection([RNA('C' * 10, metadata={'id': 's1'}),
                                 RNA('G' * 10, metadata={'id': 's2'})])
        self.assertEqual(sc.kmer_frequencies(1, relative=True),
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
        s1 = SequenceCollection([DNA("ACGT", metadata={'id': "d1"}),
                                 DNA("ACGG", metadata={'id': "d2"})])
        expected = [[0, 0.25],
                    [0.25, 0]]
        expected = DistanceMatrix(expected, ['d1', 'd2'])

        def h(s1, s2):
            return hamming(s1.values, s2.values)
        actual = s1.distances(h)
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
        expected = SequenceCollection([
            RNA('GAUUACA', metadata={'id': "r1"}),
            RNA('UUG', metadata={'id': "r2"}),
            RNA('UUGCC', metadata={'id': "r3"})])
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

    def test_update_ids_default_behavior(self):
        # 3 seqs
        exp_sc = SequenceCollection([
            RNA('GAUUACA', metadata={'id': "1"}),
            RNA('UUG', metadata={'id': "2"}),
            RNA('U-----UGCC--', metadata={'id': "3"})
        ])
        exp_id_map = {'1': 'r1', '2': 'r2', '3': 'r3'}
        obs_sc, obs_id_map = self.s2.update_ids()
        self.assertEqual(obs_sc, exp_sc)
        self.assertEqual(obs_id_map, exp_id_map)

        # empty
        obs_sc, obs_id_map = self.empty.update_ids()
        self.assertEqual(obs_sc, self.empty)
        self.assertEqual(obs_id_map, {})

    def test_update_ids_prefix(self):
        # 3 seqs
        exp_sc = SequenceCollection([
            RNA('GAUUACA', metadata={'id': "abc1"}),
            RNA('UUG', metadata={'id': "abc2"}),
            RNA('U-----UGCC--', metadata={'id': "abc3"})
        ])
        exp_id_map = {'abc1': 'r1', 'abc2': 'r2', 'abc3': 'r3'}
        obs_sc, obs_id_map = self.s2.update_ids(prefix='abc')
        self.assertEqual(obs_sc, exp_sc)
        self.assertEqual(obs_id_map, exp_id_map)

        # empty
        obs_sc, obs_id_map = self.empty.update_ids(prefix='abc')
        self.assertEqual(obs_sc, self.empty)
        self.assertEqual(obs_id_map, {})

    def test_update_ids_func_parameter(self):
        def append_42(ids):
            return [id_ + '-42' for id_ in ids]

        # 3 seqs
        exp_sc = SequenceCollection([
            RNA('GAUUACA', metadata={'id': "r1-42"}),
            RNA('UUG', metadata={'id': "r2-42"}),
            RNA('U-----UGCC--', metadata={'id': "r3-42"})
        ])
        exp_id_map = {'r1-42': 'r1', 'r2-42': 'r2', 'r3-42': 'r3'}
        obs_sc, obs_id_map = self.s2.update_ids(func=append_42)
        self.assertEqual(obs_sc, exp_sc)
        self.assertEqual(obs_id_map, exp_id_map)

        # empty
        obs_sc, obs_id_map = self.empty.update_ids(func=append_42)
        self.assertEqual(obs_sc, self.empty)
        self.assertEqual(obs_id_map, {})

    def test_update_ids_ids_parameter(self):
        # 3 seqs
        exp_sc = SequenceCollection([
            RNA('GAUUACA', metadata={'id': "abc"}),
            RNA('UUG', metadata={'id': "def"}),
            RNA('U-----UGCC--', metadata={'id': "ghi"})
        ])
        exp_id_map = {'abc': 'r1', 'def': 'r2', 'ghi': 'r3'}
        obs_sc, obs_id_map = self.s2.update_ids(ids=('abc', 'def', 'ghi'))
        self.assertEqual(obs_sc, exp_sc)
        self.assertEqual(obs_id_map, exp_id_map)

        # empty
        obs_sc, obs_id_map = self.empty.update_ids(ids=[])
        self.assertEqual(obs_sc, self.empty)
        self.assertEqual(obs_id_map, {})

    def test_update_ids_sequence_attributes_propagated(self):
        # 1 seq
        exp_sc = Alignment([
            DNA('ACGT', metadata={'id': "abc", 'description': 'desc'},
                positional_metadata={'quality': range(4)})
        ])
        exp_id_map = {'abc': 'seq1'}

        obj = Alignment([
            DNA('ACGT', metadata={'id': "seq1", 'description': 'desc'},
                positional_metadata={'quality': range(4)})
        ])

        obs_sc, obs_id_map = obj.update_ids(ids=('abc',))
        self.assertEqual(obs_sc, exp_sc)
        self.assertEqual(obs_id_map, exp_id_map)

        # 2 seqs
        exp_sc = Alignment([
            DNA('ACGT', metadata={'id': "abc", 'description': 'desc1'},
                positional_metadata={'quality': range(4)}),
            DNA('TGCA', metadata={'id': "def", 'description': 'desc2'},
                positional_metadata={'quality': range(4)[::-1]})
        ])
        exp_id_map = {'abc': 'seq1', 'def': 'seq2'}

        obj = Alignment([
            DNA('ACGT', metadata={'id': "seq1", 'description': 'desc1'},
                positional_metadata={'quality': (0, 1, 2, 3)}),
            DNA('TGCA', metadata={'id': "seq2", 'description': 'desc2'},
                positional_metadata={'quality': (3, 2, 1, 0)})
        ])

        obs_sc, obs_id_map = obj.update_ids(ids=('abc', 'def'))
        self.assertEqual(obs_sc, exp_sc)
        self.assertEqual(obs_id_map, exp_id_map)

    def test_update_ids_invalid_parameter_combos(self):
        with six.assertRaisesRegex(self, SequenceCollectionError,
                                   'ids and func'):
            self.s1.update_ids(func=lambda e: e, ids=['foo', 'bar'])

        with six.assertRaisesRegex(self, SequenceCollectionError, 'prefix'):
            self.s1.update_ids(ids=['foo', 'bar'], prefix='abc')

        with six.assertRaisesRegex(self, SequenceCollectionError, 'prefix'):
            self.s1.update_ids(func=lambda e: e, prefix='abc')

    def test_update_ids_invalid_ids(self):
        # incorrect number of new ids
        with six.assertRaisesRegex(self, SequenceCollectionError, '3 != 2'):
            self.s1.update_ids(ids=['foo', 'bar', 'baz'])
        with six.assertRaisesRegex(self, SequenceCollectionError, '4 != 2'):
            self.s1.update_ids(func=lambda e: ['foo', 'bar', 'baz', 'abc'])

        # duplicates
        with six.assertRaisesRegex(self, SequenceCollectionError, 'foo'):
            self.s2.update_ids(ids=['foo', 'bar', 'foo'])
        with six.assertRaisesRegex(self, SequenceCollectionError, 'bar'):
            self.s2.update_ids(func=lambda e: ['foo', 'bar', 'bar'])

    def test_is_empty(self):
        self.assertFalse(self.s1.is_empty())
        self.assertFalse(self.s2.is_empty())
        self.assertFalse(self.s3.is_empty())

        self.assertTrue(self.empty.is_empty())

    def test_iteritems(self):
        self.assertEqual(list(self.s1.iteritems()),
                         [(s.metadata['id'], s) for s in self.s1])

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


class AlignmentTests(TestCase):

    def setUp(self):
        self.d1 = DNA('..ACC-GTTGG..', metadata={'id': "d1"})
        self.d2 = DNA('TTACCGGT-GGCC', metadata={'id': "d2"})
        self.d3 = DNA('.-ACC-GTTGC--', metadata={'id': "d3"})

        self.r1 = RNA('UUAU-', metadata={'id': "r1"})
        self.r2 = RNA('ACGUU', metadata={'id': "r2"})

        self.seqs1 = [self.d1, self.d2, self.d3]
        self.seqs2 = [self.r1, self.r2]

        self.a1 = Alignment(self.seqs1)
        self.a2 = Alignment(self.seqs2)
        self.a3 = Alignment(self.seqs2, score=42.0,
                            start_end_positions=[(0, 3), (5, 9)])
        self.a4 = Alignment(self.seqs2, score=-42.0,
                            start_end_positions=[(1, 4), (6, 10)])

        # no sequences
        self.empty = Alignment([])

        # sequences, but no positions
        self.no_positions = Alignment([RNA('', metadata={'id': 'a'}),
                                       RNA('', metadata={'id': 'b'})])

    def test_degap(self):
        expected = SequenceCollection([
            DNA('ACCGTTGG', metadata={'id': "d1"}),
            DNA('TTACCGGTGGCC', metadata={'id': "d2"}),
            DNA('ACCGTTGC', metadata={'id': "d3"})])
        actual = self.a1.degap()
        self.assertEqual(actual, expected)

        expected = SequenceCollection([
            RNA('UUAU', metadata={'id': "r1"}),
            RNA('ACGUU', metadata={'id': "r2"})])
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

    def test_score_unset_returns_none(self):
        self.assertEqual(self.a1.score(), None)
        self.assertEqual(self.a2.score(), None)

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
        d1 = DNA('.AC', metadata={'id': "d1"})
        d2 = DNA('TAC', metadata={'id': "d2"})
        d3 = DNA('.AC', metadata={'id': "d3"})
        expected = Alignment([d1, d2, d3])
        self.assertEqual(actual, expected)

        # keep positions (invert)
        actual = self.a1.subalignment(positions_to_keep=[0, 2, 3],
                                      invert_positions_to_keep=True)
        d1 = DNA('.C-GTTGG..', metadata={'id': "d1"})
        d2 = DNA('TCGGT-GGCC', metadata={'id': "d2"})
        d3 = DNA('-C-GTTGC--', metadata={'id': "d3"})
        expected = Alignment([d1, d2, d3])
        self.assertEqual(actual, expected)

        # keep seqs and positions
        actual = self.a1.subalignment(seqs_to_keep=[0, 2],
                                      positions_to_keep=[0, 2, 3])
        d1 = DNA('.AC', metadata={'id': "d1"})
        d3 = DNA('.AC', metadata={'id': "d3"})
        expected = Alignment([d1, d3])
        self.assertEqual(actual, expected)

        # keep seqs and positions (invert)
        actual = self.a1.subalignment(seqs_to_keep=[0, 2],
                                      positions_to_keep=[0, 2, 3],
                                      invert_seqs_to_keep=True,
                                      invert_positions_to_keep=True)
        d2 = DNA('TCGGT-GGCC', metadata={'id': "d2"})
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
                        DNA('.-ACC-GTGC--', metadata={'id': "i2"})]
        self.assertRaises(AlignmentError, Alignment,
                          invalid_seqs)

    def test_init_equal_lengths(self):
        seqs = [self.d1, self.d2, self.d3]
        Alignment(seqs)

    def test_iter_positions(self):
        actual = list(self.a2.iter_positions())
        expected = [
            [RNA('U', metadata={'id': 'r1'}), RNA('A', metadata={'id': 'r2'})],
            [RNA('U', metadata={'id': 'r1'}), RNA('C', metadata={'id': 'r2'})],
            [RNA('A', metadata={'id': 'r1'}), RNA('G', metadata={'id': 'r2'})],
            [RNA('U', metadata={'id': 'r1'}), RNA('U', metadata={'id': 'r2'})],
            [RNA('-', metadata={'id': 'r1'}), RNA('U', metadata={'id': 'r2'})]
        ]
        self.assertEqual(actual, expected)

        actual = list(self.a2.iter_positions(constructor=str))
        expected = [list('UA'),
                    list('UC'),
                    list('AG'),
                    list('UU'),
                    list('-U')]
        self.assertEqual(actual, expected)

    def test_majority_consensus(self):
        # empty cases
        self.assertEqual(
            self.empty.majority_consensus(), Sequence(''))
        self.assertEqual(
            self.no_positions.majority_consensus(), RNA(''))

        # alignment where all sequences are the same
        aln = Alignment([DNA('AG', metadata={'id': 'a'}),
                         DNA('AG', metadata={'id': 'b'})])
        self.assertEqual(aln.majority_consensus(), DNA('AG'))

        # no ties
        d1 = DNA('TTT', metadata={'id': "d1"})
        d2 = DNA('TT-', metadata={'id': "d2"})
        d3 = DNA('TC-', metadata={'id': "d3"})
        a1 = Alignment([d1, d2, d3])
        self.assertEqual(a1.majority_consensus(), DNA('TT-'))

        # ties
        d1 = DNA('T', metadata={'id': "d1"})
        d2 = DNA('A', metadata={'id': "d2"})
        a1 = Alignment([d1, d2])
        self.assertTrue(a1.majority_consensus() in
                        [DNA('T'), DNA('A')])

    def test_omit_gap_positions(self):
        expected = self.a2
        self.assertEqual(self.a2.omit_gap_positions(1.0), expected)
        self.assertEqual(self.a2.omit_gap_positions(0.51), expected)

        r1 = RNA('UUAU', metadata={'id': "r1"})
        r2 = RNA('ACGU', metadata={'id': "r2"})
        expected = Alignment([r1, r2])
        self.assertEqual(self.a2.omit_gap_positions(0.49), expected)

        r1 = RNA('UUAU', metadata={'id': "r1"})
        r2 = RNA('ACGU', metadata={'id': "r2"})
        expected = Alignment([r1, r2])
        self.assertEqual(self.a2.omit_gap_positions(0.0), expected)

        self.assertEqual(self.empty.omit_gap_positions(0.0), self.empty)
        self.assertEqual(self.empty.omit_gap_positions(0.49), self.empty)
        self.assertEqual(self.empty.omit_gap_positions(1.0), self.empty)

        # Test to ensure floating point precision bug isn't present. See the
        # tests for Alignment.position_frequencies for more details.
        seqs = []
        for i in range(33):
            seqs.append(DNA('-.', metadata={'id': str(i)}))
        aln = Alignment(seqs)
        self.assertEqual(aln.omit_gap_positions(1 - np.finfo(float).eps),
                         Alignment([DNA('', metadata={'id': str(i)})
                                    for i in range(33)]))

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
        aln = Alignment([DNA('.' * 33, metadata={'id': 'abc'}),
                         DNA('-' * 33, metadata={'id': 'def'})])
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
            seqs.append(DNA('A', metadata={'id': str(i)}))
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

    def test_kmer_frequencies(self):
        expected = [defaultdict(float, {'U': 3 / 5, 'A': 1 / 5, '-': 1 / 5}),
                    defaultdict(float, {'A': 1 / 5, 'C': 1 / 5, 'G': 1 / 5,
                                        'U': 2 / 5})]
        actual = self.a2.kmer_frequencies(k=1, relative=True)
        for a, e in zip(actual, expected):
            self.assertEqual(sorted(a), sorted(e), 5)
            np.testing.assert_almost_equal(sorted(a.values()),
                                           sorted(e.values()), 5)

    def test_sequence_length(self):
        self.assertEqual(self.a1.sequence_length(), 13)
        self.assertEqual(self.a2.sequence_length(), 5)
        self.assertEqual(self.empty.sequence_length(), 0)

    def test_validate_lengths(self):
        self.assertTrue(self.a1._validate_lengths())
        self.assertTrue(self.a2._validate_lengths())
        self.assertTrue(self.empty._validate_lengths())

        self.assertTrue(Alignment([
            DNA('TTT', metadata={'id': "d1"})])._validate_lengths())

if __name__ == "__main__":
    main()
