# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import itertools
import unittest

import six
import numpy as np
import numpy.testing as npt

from skbio import Sequence, DNA
from skbio.sequence.distance import hamming


class TestHamming(unittest.TestCase):
    def test_non_sequence(self):
        seq1 = Sequence('abc')
        seq2 = 'abc'

        with six.assertRaisesRegex(self, TypeError,
                                   'seq1.*seq2.*Sequence.*str'):
            hamming(seq1, seq2)

        with six.assertRaisesRegex(self, TypeError,
                                   'seq1.*seq2.*Sequence.*str'):
            hamming(seq2, seq1)

    def test_type_mismatch(self):
        seq1 = Sequence('ABC')
        seq2 = DNA('ACG')

        with six.assertRaisesRegex(self, TypeError,
                                   'Sequence.*does not match.*DNA'):
            hamming(seq1, seq2)

    def test_length_mismatch(self):
        seq1 = Sequence('ABC')
        seq2 = Sequence('ABCD')

        with six.assertRaisesRegex(self, ValueError, 'equal length.*3 != 4'):
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


if __name__ == "__main__":
    unittest.main()
