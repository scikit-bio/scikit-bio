# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import itertools
import unittest

import numpy as np
import numpy.testing as npt

from skbio import Sequence, DNA
from skbio.sequence.distance import hamming, kmer_distance


class TestHamming(unittest.TestCase):
    def test_non_sequence(self):
        seq1 = Sequence('abc')
        seq2 = 'abc'

        with self.assertRaisesRegex(TypeError, r'seq1.*seq2.*Sequence.*str'):
            hamming(seq1, seq2)

        with self.assertRaisesRegex(TypeError, r'seq1.*seq2.*Sequence.*str'):
            hamming(seq2, seq1)

    def test_type_mismatch(self):
        seq1 = Sequence('ABC')
        seq2 = DNA('ACG')

        with self.assertRaisesRegex(TypeError,
                                    r'Sequence.*does not match.*DNA'):
            hamming(seq1, seq2)

    def test_length_mismatch(self):
        seq1 = Sequence('ABC')
        seq2 = Sequence('ABCD')

        with self.assertRaisesRegex(ValueError, r'equal length.*3 != 4'):
            hamming(seq1, seq2)

    def test_return_type(self):
        seq1 = Sequence('ABC')
        seq2 = Sequence('ABC')

        distance = hamming(seq1, seq2)

        self.assertIsInstance(distance, float)
        self.assertEqual(distance, 0.0)

    def test_minimum_distance(self):
        seq1 = Sequence('ABC')
        seq2 = Sequence('ABC')

        distance = hamming(seq1, seq2)

        self.assertEqual(distance, 0.0)

    def test_mid_range_distance(self):
        seq1 = Sequence("abcdefgh")
        seq2 = Sequence("1b23ef45")

        distance = hamming(seq1, seq2)

        self.assertEqual(distance, 5.0/8.0)

    def test_maximum_distance(self):
        seq1 = Sequence('ABC')
        seq2 = Sequence('CAB')

        distance = hamming(seq1, seq2)

        self.assertEqual(distance, 1.0)

    def test_empty_sequences(self):
        seq1 = Sequence('')
        seq2 = Sequence('')

        distance = hamming(seq1, seq2)

        npt.assert_equal(distance, np.nan)

    def test_single_character_sequences(self):
        seq1 = Sequence('a')
        seq2 = Sequence('b')

        self.assertEqual(hamming(seq1, seq1), 0.0)
        self.assertEqual(hamming(seq1, seq2), 1.0)

    def test_sequence_subclass(self):
        seq1 = DNA('ACG-T')
        seq2 = DNA('ACCTT')

        distance = hamming(seq1, seq2)

        self.assertEqual(distance, 2.0/5.0)

    def test_sequences_with_metadata(self):
        # test for #1254
        seqs1 = [
            Sequence("ACGT"),
            Sequence("ACGT", metadata={'id': 'abc'}),
            Sequence("ACGT", positional_metadata={'qual': range(4)})
        ]
        seqs2 = [
            Sequence("AAAA"),
            Sequence("AAAA", metadata={'id': 'def'}),
            Sequence("AAAA", positional_metadata={'qual': range(4, 8)})
        ]

        for seqs in seqs1, seqs2:
            for seq1, seq2 in itertools.product(seqs, repeat=2):
                distance = hamming(seq1, seq2)
                self.assertEqual(distance, 0.0)

        for seq1, seq2 in itertools.product(seqs1, seqs2):
            distance = hamming(seq1, seq2)
            self.assertEqual(distance, 0.75)


class TestKmerDistance(unittest.TestCase):
    def test_default_kwargs(self):
        seq1 = Sequence('AACCTAGCAATGGAT')
        seq2 = Sequence('CAGGCAGTTCTCACC')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 0.9130434782608695
        self.assertAlmostEqual(obs, exp)

    def test_nondefault_k(self):
        seq1 = Sequence('GCTTATGGAGAGAGA')
        seq2 = Sequence('CTCGAACTCCAGCCA')
        obs = kmer_distance(seq1, seq2, 2)
        exp = 0.7333333333333333
        self.assertAlmostEqual(obs, exp)
        seq1 = Sequence('EADDECAEECDEACD')
        seq2 = Sequence('DCBCBADADABCCDA')
        obs = kmer_distance(seq1, seq2, 1)
        exp = 0.4
        self.assertAlmostEqual(obs, exp)

    def test_overlap_false(self):
        seq1 = Sequence('CGTTATGTCTGTGAT')
        seq2 = Sequence('CTGAATCGGTAGTGT')
        obs = kmer_distance(seq1, seq2, 3, overlap=False)
        exp = 0.8888888888888888
        self.assertAlmostEqual(obs, exp)

    def test_entirely_different_sequences(self):
        seq1 = Sequence('CCGTGGTCGTATAAG')
        seq2 = Sequence('CGCCTTCCACATCAG')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 1.0
        self.assertEqual(obs, exp)

    def test_same_sequence(self):
        seq1 = Sequence('CTGCGACAGTTGGTA')
        seq2 = Sequence('CTGCGACAGTTGGTA')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 0.0
        self.assertEqual(obs, exp)

    def test_differing_length_seqs(self):
        seq1 = Sequence('AGAAATCTGAGCAAGGATCA')
        seq2 = Sequence('TTAGTGCGTAATCCG')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 0.9285714285714286
        self.assertAlmostEqual(obs, exp)

    def test_with_sequence_subclass(self):
        seq1 = DNA('GATGGTACTGTAGGT')
        seq2 = DNA('AGGGTGAAGGTATCA')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 0.8421052631578947
        self.assertAlmostEqual(obs, exp)

    def test_with_metadata_sanity(self):
        seq1 = Sequence('AACCTAGCAATGGAT',
                        metadata={'Name': 'Kestrel Gorlick'},
                        positional_metadata={'seq': list('ACTCAAGCTACGAAG')})
        seq2 = Sequence('CAGGCAGTTCTCACC')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 0.9130434782608695
        self.assertAlmostEqual(obs, exp)

    def test_return_type(self):
        seq1 = Sequence('ATCG')
        seq2 = Sequence('ATCG')
        obs = kmer_distance(seq1, seq2, 3)
        self.assertIsInstance(obs, float)
        self.assertEqual(obs, 0.0)

    def test_empty_sequences(self):
        seq1 = Sequence('')
        seq2 = Sequence('')
        obs = kmer_distance(seq1, seq2, 3)
        npt.assert_equal(obs, np.nan)

    def test_one_empty_sequence(self):
        seq1 = Sequence('')
        seq2 = Sequence('CGGGCAGCTCCTACCTGCTA')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 1.0
        self.assertAlmostEqual(obs, exp)

    def test_no_kmers_found(self):
        seq1 = Sequence('ATCG')
        seq2 = Sequence('ACGT')
        obs = kmer_distance(seq1, seq2, 5)
        npt.assert_equal(obs, np.nan)

    def test_k_less_than_one_error(self):
        seq1 = Sequence('ATCG')
        seq2 = Sequence('ACTG')
        with self.assertRaisesRegex(ValueError, r'k must be greater than 0.'):
            kmer_distance(seq1, seq2, 0)

    def test_type_mismatch_error(self):
        seq1 = Sequence('ABC')
        seq2 = DNA('ATC')
        with self.assertRaisesRegex(TypeError, r"Type 'Sequence'.*type 'DNA'"):
            kmer_distance(seq1, seq2, 3)

    def test_non_sequence_error(self):
        seq1 = Sequence('ATCG')
        seq2 = 'ATCG'
        with self.assertRaisesRegex(TypeError, r"not 'str'"):
            kmer_distance(seq1, seq2, 3)


if __name__ == "__main__":
    unittest.main()
