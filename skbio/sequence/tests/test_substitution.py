# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
from numpy.testing import assert_array_equal

from skbio import SubstitutionMatrix


class TestSubstitutionMatrix(TestCase):
    def setUp(self):
        self.alphabet = 'ACGTN'
        self.scores = np.array([
            [1, -2, -2, -2, 0],
            [-2, 1, -2, -2, 0],
            [-2, -2, 1, -2, 0],
            [-2, -2, -2, 1, 0],
            [0, 0, 0, 0, 0]])

    def test_init(self):
        # typical usage
        # alphabet becomes tuple of characters
        alphabet = tuple(self.alphabet)
        obs = SubstitutionMatrix(self.alphabet, self.scores)
        self.assertTupleEqual(obs.alphabet, alphabet)

        # alphabet is an alias of ids
        self.assertTupleEqual(obs.alphabet, obs.ids)

        # matrix is ndarray (this is important for alignment efficiency)
        self.assertTrue(isinstance(obs.scores, np.ndarray))
        self.assertTupleEqual(obs.shape, (5, 5))
        assert_array_equal(obs.scores, self.scores)

        # data type becomes float
        self.assertEqual(obs.dtype, np.float64)

        # scores is an alias of data
        assert_array_equal(obs.scores, obs.data)

        # character to index mapping
        self.assertDictEqual(obs._char_map, dict(zip(
            alphabet, range(len(alphabet)))))

        # alphabet can be encoded as ASCII characters
        self.assertTrue(obs._is_ascii)

        # hash table of ASCII characters
        self.assertTrue(isinstance(obs._char_hash, np.ndarray))
        self.assertTrue(obs._char_hash.dtype.type is np.uint8)
        for i, char in enumerate(alphabet):
            self.assertEqual(i, obs._char_hash[ord(char)])

    def test_init_alt_alphabet(self):
        # alternative formats of alphabet: list, dictionary (only keys matter),
        # and iterator
        alphabet = tuple(self.alphabet)
        for alp in (list(alphabet),
                    dict.fromkeys(alphabet),
                    iter(alphabet)):
            obs = SubstitutionMatrix(alp, self.scores)
            self.assertTupleEqual(obs.alphabet, alphabet)

    def test_init_alt_scores(self):
        # alternative format of scores: nested list
        obs = SubstitutionMatrix(self.alphabet, self.scores.tolist())
        assert_array_equal(obs.scores, self.scores)

        # condensed matrix (less likely because diagonal is zero)
        obs = SubstitutionMatrix('ACGT', [-1] * 6)
        assert_array_equal(obs.scores, np.identity(4) - 1)

    def test_to_dict(self):
        mat = SubstitutionMatrix(self.alphabet, self.scores)
        obs = mat.to_dict()
        exp = {'A': {'A': 1., 'C': -2., 'G': -2., 'T': -2., 'N': 0.},
               'C': {'A': -2., 'C': 1., 'G': -2., 'T': -2., 'N': 0.},
               'G': {'A': -2., 'C': -2., 'G': 1., 'T': -2., 'N': 0.},
               'T': {'A': -2., 'C': -2., 'G': -2., 'T': 1., 'N': 0.},
               'N': {'A': 0., 'C': 0., 'G': 0., 'T': 0., 'N': 0.}}
        self.assertDictEqual(obs, exp)

    def test_from_dict(self):
        d = {'a': {'a': 1, 'b': 0, 'c': 0},
             'b': {'a': 0, 'b': 1, 'c': 0},
             'c': {'a': 0, 'b': 0, 'c': 1}}
        obs = SubstitutionMatrix.from_dict(d)
        self.assertTrue(isinstance(obs, SubstitutionMatrix))
        self.assertTupleEqual(obs.alphabet, tuple('abc'))
        exp = np.array([[1., 0., 0.],
                        [0., 1., 0.],
                        [0., 0., 1.]])
        assert_array_equal(obs.data, exp)

        # alphabet is inconsistent
        msg = ('The outer and inner layers of the dictionary must have the '
               'same set of keys.')
        d['d'] = {'a': 0, 'b': 0, 'c': 0}
        with self.assertRaisesRegex(ValueError, msg):
            SubstitutionMatrix.from_dict(d)
        del d['d']
        d['a']['d'] = 2
        with self.assertRaisesRegex(ValueError, msg):
            SubstitutionMatrix.from_dict(d)
        del d['a']['d']

        # scores are not numbers
        d['a']['b'] = 'hello'
        with self.assertRaises(ValueError):
            SubstitutionMatrix.from_dict(d)
        d['a']['b'] = None
        with self.assertRaises(TypeError):
            SubstitutionMatrix.from_dict(d)

    def test_identity(self):
        obs = SubstitutionMatrix.identity('ACGT', 1, -2)
        self.assertTrue(isinstance(obs, SubstitutionMatrix))
        self.assertTupleEqual(obs.alphabet, tuple('ACGT'))
        exp = np.array([[1., -2., -2., -2.],
                        [-2., 1., -2., -2.],
                        [-2., -2., 1., -2.],
                        [-2., -2., -2., 1.]])
        assert_array_equal(obs.scores, exp)

    def test_by_name(self):
        obs = SubstitutionMatrix.by_name('NUC.4.4')
        self.assertEqual(len(obs.alphabet), 15)
        self.assertEqual(obs['A', 'T'], -4)
        obs = SubstitutionMatrix.by_name('BLOSUM50')
        self.assertEqual(len(obs.alphabet), 24)
        self.assertEqual(obs['M', 'K'], -2)
        obs = SubstitutionMatrix.by_name('blosum50')
        self.assertEqual(len(obs.alphabet), 24)
        self.assertEqual(obs['M', 'K'], -2)
        msg = 'Substitution matrix "hello" does not exist.'
        with self.assertRaisesRegex(ValueError, msg):
            SubstitutionMatrix.by_name('hello')

    def test_get_names(self):
        obs = SubstitutionMatrix.get_names()
        self.assertTrue('NUC.4.4' in obs)
        self.assertTrue('PAM250' in obs)
        self.assertTrue('BLOSUM62' in obs)


if __name__ == "__main__":
    main()
