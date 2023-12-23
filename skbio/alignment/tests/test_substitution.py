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
        self.data = np.array([
            [1, -2, -2, -2, 0],
            [-2, 1, -2, -2, 0],
            [-2, -2, 1, -2, 0],
            [-2, -2, -2, 1, 0],
            [0, 0, 0, 0, 0]])

    def test_init_typical(self):
        # typical usage
        # alphabet becomes tuple of characters
        alphabet = tuple(self.alphabet)
        obs = SubstitutionMatrix(self.data, self.alphabet)
        self.assertTupleEqual(obs.ids, alphabet)

        # alphabet is an alias of ids
        self.assertTupleEqual(obs.alphabet, alphabet)

        # matrix is ndarray (this is important for alignment efficiency)
        self.assertTrue(isinstance(obs.data, np.ndarray))
        self.assertTupleEqual(obs.shape, (5, 5))
        assert_array_equal(obs.data, self.data)

        # data type becomes float
        self.assertEqual(obs.dtype, np.float64)

    def test_init_alt_alphabet(self):
        # alternative formats of alphabet: list, dictionary (only keys matter),
        # and iterator
        alphabet = tuple(self.alphabet)
        for alp in (list(self.alphabet),
                    dict.fromkeys(self.alphabet),
                    iter(self.alphabet)):
            obs = SubstitutionMatrix(self.data, alp)
            self.assertTupleEqual(obs.alphabet, alphabet)

    def test_to_dict(self):
        mat = SubstitutionMatrix(self.data, self.alphabet)
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
        msg = 'Scores must be integers or floating-point numbers.'
        d['a']['b'] = 'hello'
        with self.assertRaisesRegex(ValueError, msg):
            SubstitutionMatrix.from_dict(d)
        d['a']['b'] = None
        with self.assertRaisesRegex(ValueError, msg):
            SubstitutionMatrix.from_dict(d)

    def test_identity(self):
        obs = SubstitutionMatrix.identity(1, -2, 'ACGT')
        self.assertTrue(isinstance(obs, SubstitutionMatrix))
        self.assertTupleEqual(obs.ids, tuple('ACGT'))
        exp = np.array([[1., -2., -2., -2.],
                        [-2., 1., -2., -2.],
                        [-2., -2., 1., -2.],
                        [-2., -2., -2., 1.]])
        assert_array_equal(obs.data, exp)

        # # alphabet is a dictionary (only keys are considered)
        # submat = SubstitutionMatrix(dict.fromkeys(self.alphabet), self.matrix)
        # self.assertListEqual(submat.alphabet, alphabet)

        # # alphabet is an iterator
        # submat = SubstitutionMatrix(iter(self.alphabet), self.matrix)
        # self.assertListEqual(submat.alphabet, alphabet)

        # # alphabet is not an iterable
        # msg = ('Alphabet must be iterable, such as a string or a list of '
        #        'characters.')
        # with self.assertRaisesRegex(ValueError, msg):
        #     SubstitutionMatrix(0, self.matrix)
        # with self.assertRaisesRegex(ValueError, msg):
        #     SubstitutionMatrix(None, self.matrix)
        # with self.assertRaisesRegex(ValueError, msg):
        #     SubstitutionMatrix(None, self.matrix)

        # # alphabet is empty
        # msg = 'Alphabet must not be empty.'
        # with self.assertRaisesRegex(ValueError, msg):
        #     SubstitutionMatrix('', self.matrix)
        # with self.assertRaisesRegex(ValueError, msg):
        #     SubstitutionMatrix([], self.matrix)

        # # alphabet contains duplicates
        # msg = 'Alphabet must not contain duplicated elements.'
        # with self.assertRaisesRegex(ValueError, msg):
        #     SubstitutionMatrix('hello', self.matrix)

    # def test_init(self):
    #     # typical usage (alphabet becomes list, matrix becomes 2D ndarray)
    #     alphabet = list(self.alphabet)
    #     submat = SubstitutionMatrix(self.alphabet, self.matrix)
    #     self.assertListEqual(submat.alphabet, alphabet)

    #     # alphabet is a list
    #     submat = SubstitutionMatrix(alphabet, self.matrix)
    #     self.assertListEqual(submat.alphabet, alphabet)

    #     # alphabet is a dictionary (only keys are considered)
    #     submat = SubstitutionMatrix(dict.fromkeys(self.alphabet), self.matrix)
    #     self.assertListEqual(submat.alphabet, alphabet)

    #     # alphabet is an iterator
    #     submat = SubstitutionMatrix(iter(self.alphabet), self.matrix)
    #     self.assertListEqual(submat.alphabet, alphabet)

    #     # alphabet is not an iterable
    #     msg = ('Alphabet must be iterable, such as a string or a list of '
    #            'characters.')
    #     with self.assertRaisesRegex(ValueError, msg):
    #         SubstitutionMatrix(0, self.matrix)
    #     with self.assertRaisesRegex(ValueError, msg):
    #         SubstitutionMatrix(None, self.matrix)
    #     with self.assertRaisesRegex(ValueError, msg):
    #         SubstitutionMatrix(None, self.matrix)

    #     # alphabet is empty
    #     msg = 'Alphabet must not be empty.'
    #     with self.assertRaisesRegex(ValueError, msg):
    #         SubstitutionMatrix('', self.matrix)
    #     with self.assertRaisesRegex(ValueError, msg):
    #         SubstitutionMatrix([], self.matrix)

    #     # alphabet contains duplicates
    #     msg = 'Alphabet must not contain duplicated elements.'
    #     with self.assertRaisesRegex(ValueError, msg):
    #         SubstitutionMatrix('hello', self.matrix)


if __name__ == "__main__":
    main()
