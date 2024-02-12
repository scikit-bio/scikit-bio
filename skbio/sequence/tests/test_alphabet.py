# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio.sequence._alphabet import (
    _encode_alphabet, _alphabet_to_hashes,
    _indices_in_alphabet, _indices_in_alphabet_ascii,
    _indices_in_observed)


class TestAlphabet(TestCase):

    def test_encode_alphabet(self):
        # ascii characters
        alpha = 'ACGT'
        exp = np.array([65, 67, 71, 84], dtype=np.uint8)
        npt.assert_equal(_encode_alphabet(alpha), exp)
        npt.assert_equal(_encode_alphabet(list(alpha)), exp)
        npt.assert_equal(_encode_alphabet(tuple(alpha)), exp)
        npt.assert_equal(_encode_alphabet(np.array(list(alpha))), exp)
        npt.assert_equal(_encode_alphabet(np.char.encode(list(alpha))), exp)

        # ascii code points
        codes = list(map(ord, alpha))
        npt.assert_equal(_encode_alphabet(codes), exp)
        npt.assert_equal(_encode_alphabet(np.array(codes)), exp)
        npt.assert_equal(_encode_alphabet(np.array(codes).astype(
            np.uint8)), exp)

        # wrong data types
        with self.assertRaises(TypeError):
            _encode_alphabet(123)
        with self.assertRaises(TypeError):
            _encode_alphabet(set(alpha))
        with self.assertRaises(TypeError):
            _encode_alphabet([1.0, 1.5, 2.0])
        with self.assertRaises(TypeError):
            _encode_alphabet([['a', 'b'], ['c', 'd']])

        # not single characters
        with self.assertRaises(ValueError):
            _encode_alphabet(['this', 'is', 'not'])

        # exceed ascii range
        with self.assertRaises(ValueError):
            _encode_alphabet([100, 200, 300])
        with self.assertRaises(UnicodeEncodeError):
            _encode_alphabet(chr(1234) + chr(5678))
        with self.assertRaises(UnicodeEncodeError):
            _encode_alphabet([chr(1234), chr(5678)])

    def test_alphabet_to_hashes(self):
        alpha = 'ATGCSWRYKMBVHDN'
        obs = _alphabet_to_hashes(alpha)
        exp = np.array([
            255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
            0, 10, 3, 13, 255, 255, 2, 12, 255, 255, 8, 255, 9,
            14, 255, 255, 255, 6, 4, 1, 255, 11, 5, 255, 7, 255,
            255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
            255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255],
            dtype=np.uint8)
        npt.assert_equal(obs, exp)

    def test_indices_in_alphabet(self):
        seq = 'GAGCTCA'  # DNA sequence without degenerate characters
        alpha = 'ACGTN'  # DNA alphabet
        exp = np.array([2, 0, 2, 1, 3, 1, 0])  # indices of characters

        # either alphabet or sequence may be string, list/tuple, or iterator
        npt.assert_equal(_indices_in_alphabet(seq, alpha), exp)
        npt.assert_equal(_indices_in_alphabet(list(seq), list(alpha)), exp)
        npt.assert_equal(_indices_in_alphabet(iter(seq), iter(alpha)), exp)

        # alphabet is a dictionary of character : index (most performant)
        npt.assert_equal(_indices_in_alphabet(seq, dict(zip(alpha, range(len(
            alpha))))), exp)

        # one character is absent from alphabet
        seq = 'GAGRCTCA'
        msg = ('One or multiple characters in the sequence are absent from '
               'the alphabet.')
        with self.assertRaisesRegex(ValueError, msg):
            _indices_in_alphabet(seq, alpha)

        # replace absent character with wildcard
        obs = _indices_in_alphabet(seq, alpha, wildcard='N')
        exp = np.array([2, 0, 2, 4, 1, 3, 1, 0])
        npt.assert_equal(obs, exp)

        # wildcard not in alphabet
        msg = 'Wildcard character "X" is not in the alphabet.'
        with self.assertRaisesRegex(ValueError, msg):
            _indices_in_alphabet(seq, alpha, wildcard='X')

        # amino acid
        seq = 'MEEPQSDPSV'
        alpha = 'ARNDCQEGHILKMFPSTWYVBZX'
        exp = np.array([12, 6, 6, 14, 5, 15, 3, 14, 15, 19])
        npt.assert_equal(_indices_in_alphabet(seq, alpha), exp)

        # natural language
        seq = 'The quick brown fox jumps over the lazy dog'.split()
        alpha = ['dog', 'fox', 'jumps', 'the']
        obs = _indices_in_alphabet(seq, alpha, wildcard='the')
        exp = np.array([3, 3, 3, 1, 2, 3, 3, 3, 0])
        npt.assert_equal(obs, exp)

        # empty sequence
        self.assertEqual(_indices_in_alphabet('', alpha).size, 0)

    def test_indices_in_alphabet_ascii(self):
        # convert a sequence into a vector of code points
        seq = 'GAGCTCA'
        seq = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)

        # convert an alphabet into a vector of indices
        alpha = 'ACGTN'
        idx = np.frombuffer(alpha.encode('ascii'), dtype=np.uint8)
        alpha = np.full(128, 255, dtype=np.uint8)
        alpha[idx] = np.arange(idx.size)

        # a typical case
        obs = _indices_in_alphabet_ascii(seq, alpha)
        exp = np.array([2, 0, 2, 1, 3, 1, 0])
        npt.assert_equal(obs, exp)
        self.assertTrue(obs.dtype.type is np.uint8)

        # one character is absent
        seq = np.insert(seq, 3, ord('R'))
        msg = ('One or multiple characters in the sequence are absent from '
               'the alphabet.')
        with self.assertRaisesRegex(ValueError, msg):
            _indices_in_alphabet_ascii(seq, alpha)

        # replace absent character
        obs = _indices_in_alphabet_ascii(seq, alpha, wildcard=ord('N'))
        exp = np.array([2, 0, 2, 4, 1, 3, 1, 0])
        npt.assert_equal(obs, exp)
        self.assertTrue(obs.dtype.type is np.uint8)

        # wildcard not in alphabet
        msg = 'Wildcard character "&" is not in the alphabet.'
        with self.assertRaisesRegex(ValueError, msg):
            _indices_in_alphabet_ascii(seq, alpha, wildcard=38)

    def test_indices_in_observed(self):
        # data from human TP53 protein (NP_000537.3)
        seqs = ('MEEPQSDPSVEPPLSQETFSDLWKLLPE',
                'NNVLSPLPSQAMDDLMLSP',
                'DDIEQWFTEDPGPDEAPRMPEAA')

        obs_idx, obs_alp = _indices_in_observed(seqs)
        exp_alp = np.array(tuple('ADEFGIKLMNPQRSTVW'))
        exp_idx = (
            np.array([8, 2, 2, 10, 11, 13, 1, 10, 13, 15, 2, 10, 10, 7, 13, 11,
                      2, 14, 3, 13, 1, 7, 16, 6, 7, 7, 10, 2]),
            np.array([9, 9, 15, 7, 13, 10, 7, 10, 13, 11, 0, 8, 1, 1, 7, 8, 7,
                      13, 10]),
            np.array([1, 1, 5, 2, 11, 16, 3, 14, 2, 1, 10, 4, 10, 1, 2, 0, 10,
                      12, 8, 10, 2, 0, 0]))
        npt.assert_equal(obs_alp, exp_alp)
        for obs, exp in zip(obs_idx, exp_idx):
            npt.assert_equal(obs, exp)

        # reconstruct original sequences
        for idx, seq in zip(obs_idx, seqs):
            self.assertEqual(''.join(obs_alp[idx]), seq)

        # sequences are numbers
        seqs = ([1, 4, 6, 7, 8],
                [3, 3, 4, 1, 0],
                [5, 2, 5, 8, 0])
        obs_idx, obs_alp = _indices_in_observed(seqs)
        npt.assert_equal(obs_alp, np.arange(9))
        for idx, seq in zip(obs_idx, seqs):
            npt.assert_equal(obs_alp[idx], np.array(seq))

        # sequences are natural language
        seqs = (['this', 'is', 'a', 'cat'],
                ['that', 'is', 'a', 'dog'],
                ['cat', 'is', 'not', 'dog'])
        obs_idx, obs_alp = _indices_in_observed(seqs)
        exp_alp = np.unique(np.concatenate(seqs))
        npt.assert_equal(obs_alp, exp_alp)
        for idx, seq in zip(obs_idx, seqs):
            npt.assert_equal(obs_alp[idx], np.array(seq))

        # sequences are individual characters
        obs_idx, obs_alp = _indices_in_observed(['hello'])
        npt.assert_equal(obs_alp, np.array(['e', 'h', 'l', 'o']))
        self.assertEqual(''.join(obs_alp[np.concatenate(obs_idx)]), 'hello')

        # empty sequence
        obs_idx, obs_alp = _indices_in_observed([[]])
        self.assertEqual(obs_alp.size, 0)
        self.assertEqual(len(obs_idx), 1)
        self.assertEqual(obs_idx[0].size, 0)


if __name__ == "__main__":
    main()
