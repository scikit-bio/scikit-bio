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

from skbio.stats.distance._utils import (
    is_symmetric_and_hollow,
    is_symmetric,
    is_hollow,
    distmat_reorder,
    distmat_reorder_condensed,
)


class UtilsTests(TestCase):
    def test_is_symmetric_and_hollow(self):
        # a symmetric and hollow matrix
        data = np.array([
            [0, 1, 2, 3],
            [1, 0, 4, 5],
            [2, 4, 0, 6],
            [3, 5, 6, 0],
        ]).astype(float)
        obs = is_symmetric_and_hollow(data)
        exp = [True, True]
        self.assertListEqual(obs, exp)

        # convert matrix to F-contiguous order
        self.assertListEqual(is_symmetric_and_hollow(np.asfortranarray(data)), exp)

        # a symmetric but non-hollow matrix
        self.assertListEqual(is_symmetric_and_hollow(np.array([
            [1, 2, 3],
            [2, 4, 5],
            [3, 5, 6],
        ]).astype(float)), [True, False])

        # a hollow but asymmetric matrix
        self.assertListEqual(is_symmetric_and_hollow(np.array([
            [0, 1, 3],
            [2, 0, 5],
            [4, 6, 0],
        ]).astype(float)), [False, True])

        # an asymmetric and non-hollow matrix
        self.assertListEqual(is_symmetric_and_hollow(np.array([
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9],
        ]).astype(float)), [False, False])

    def test_is_symmetric(self):
        self.assertTrue(is_symmetric(np.array([
            [1, 2, 3], [2, 4, 5], [3, 5, 6],
        ]).astype(float)))

        self.assertFalse(is_symmetric(np.array([
            [0, 1, 3], [2, 0, 5], [4, 6, 0],
        ]).astype(float)))

    def test_is_hollow(self):
        self.assertTrue(is_hollow(np.array([
            [0, 1, 3], [2, 0, 5], [4, 6, 0],
        ]).astype(float)))

        self.assertFalse(is_hollow(np.array([
            [1, 2, 3], [2, 4, 5], [3, 5, 6],
        ]).astype(float)))

    def test_distmat_reorder(self):
        mat = np.array([
            [0, 1, 2, 3],
            [1, 0, 4, 5],
            [2, 4, 0, 6],
            [3, 5, 6, 0],
        ]).astype(float)
        vec = np.array([1, 0, 3, 2])
        obs = distmat_reorder(mat, vec)
        exp = np.array([
            [0, 1, 5, 4],
            [1, 0, 3, 2],
            [5, 3, 0, 6],
            [4, 2, 6, 0],
        ]).astype(float)
        npt.assert_array_equal(obs, exp)

        # validation
        obs = distmat_reorder(mat, vec, validate=True)
        npt.assert_array_equal(obs, exp)

        # F-contiguous order
        obs = distmat_reorder(np.asfortranarray(mat), vec)
        npt.assert_array_equal(obs, exp)

        # partial and repeated elements
        obs = distmat_reorder(mat, [2, 1, 2])
        exp = np.array([
            [0, 4, 0],
            [4, 0, 4],
            [0, 4, 0],
        ]).astype(float)
        npt.assert_array_equal(obs, exp)

        # invalid order vectors
        with self.assertRaises(ValueError):
            distmat_reorder(mat, np.array([1, 100]), validate=True)

        with self.assertRaises(ValueError):
            distmat_reorder(mat, np.array([-1, 2, 3]), validate=True)

    def test_distmat_reorder_condensed(self):
        mat = np.array([
            [0, 1, 2, 3],
            [1, 0, 4, 5],
            [2, 4, 0, 6],
            [3, 5, 6, 0],
        ]).astype(float)
        vec = np.array([1, 0, 3, 2])
        obs = distmat_reorder_condensed(mat, vec)
        exp = np.array([1, 5, 4, 3, 2, 6]).astype(float)
        npt.assert_array_equal(obs, exp)

        # validation
        obs = distmat_reorder_condensed(mat, vec, validate=True)
        npt.assert_array_equal(obs, exp)

        # F-contiguous order
        obs = distmat_reorder_condensed(np.asfortranarray(mat), vec)
        npt.assert_array_equal(obs, exp)

        # partial and repeated elements
        obs = distmat_reorder_condensed(mat, [2, 1, 2])
        exp = np.array([4, 0, 4]).astype(float)
        npt.assert_array_equal(obs, exp)


if __name__ == '__main__':
    main()
