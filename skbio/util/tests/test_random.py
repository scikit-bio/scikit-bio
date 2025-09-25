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

from skbio.util._random import get_rng


class RandomTests(TestCase):

    def test_get_rng(self):

        # returns random generator
        res = get_rng()
        self.assertIsInstance(res, np.random.Generator)

        # seed is Python integer
        res = get_rng(42)
        obs = np.array([res.integers(100) for _ in range(5)])
        exp = np.array([8, 77, 65, 43, 43])
        npt.assert_array_equal(obs, exp)

        # seed is NumPy integer
        res = get_rng(np.uint8(42))
        obs = np.array([res.integers(100) for _ in range(5)])
        npt.assert_array_equal(obs, exp)

        # seed is new Generator
        res = get_rng(res)
        obs = np.array([res.integers(100) for _ in range(5)])
        exp = np.array([85, 8, 69, 20, 9])
        npt.assert_array_equal(obs, exp)

        # test if integer seeds are disjoint
        obs = [get_rng(i).integers(1e6) for i in range(10)]
        exp = [850624, 473188, 837575, 811504, 726442,
               670790, 445045, 944904, 719549, 421547]
        self.assertListEqual(obs, exp)

        # no seed: use current random state
        np.random.seed(42)
        res = get_rng()
        obs = np.array([res.integers(100) for _ in range(5)])
        exp = np.array([90, 11, 93, 94, 34])
        npt.assert_array_equal(obs, exp)

        # reset random state to reproduce output
        np.random.seed(42)
        res = get_rng()
        obs = np.array([res.integers(100) for _ in range(5)])
        npt.assert_array_equal(obs, exp)

        # call also advances current random state
        np.random.seed(42)
        self.assertEqual(np.random.randint(100), 51)
        res = get_rng()
        self.assertEqual(np.random.randint(100), 14)

        # seed is legacy RandomState
        res = get_rng(np.random.RandomState(42))
        obs = np.array([res.integers(100) for _ in range(5)])
        npt.assert_array_equal(obs, exp)

        # test if legacy random states are disjoint
        obs = [get_rng(np.random.RandomState(i)).integers(1e6) for i in range(5)]
        exp = [368454, 346004, 189187, 324799, 924851]
        self.assertListEqual(obs, exp)

        # invalid seed
        msg = 'Invalid seed. It must be an integer or a random generator instance.'
        with self.assertRaises(ValueError) as cm:
            get_rng('hello')
        self.assertEqual(str(cm.exception), msg)

        # mimic legacy NumPy
        delattr(np.random, 'Generator')
        msg = ('The installed NumPy version does not support random.Generator. '
               'Please use NumPy >= 1.17.')
        with self.assertRaises(ValueError) as cm:
            get_rng()
        self.assertEqual(str(cm.exception), msg)
        with self.assertRaises(ValueError) as cm:
            get_rng('hello')
        self.assertEqual(str(cm.exception), msg)


if __name__ == '__main__':
    main()
