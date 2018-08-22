# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import copy

import numpy as np
import numpy.testing as npt

from skbio.stats.ordination import corr, mean_and_std, e_matrix, f_matrix, \
    center_distance_matrix

from skbio.stats.ordination._utils import _e_matrix_inplace, _f_matrix_inplace


class TestUtils(TestCase):
    def setUp(self):
        self.x = np.array([[1, 2, 3], [4, 5, 6]])
        self.y = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

        self.matrix = np.arange(1, 7).reshape(2, 3)
        self.matrix2 = np.arange(1, 10).reshape(3, 3)

        self.small_mat = np.array([[7, 5, 5], [4, 4, 9], [7, 5, 3]])

    def test_mean_and_std(self):
        obs = mean_and_std(self.x)
        npt.assert_almost_equal((3.5, 1.707825127), obs)

        obs = mean_and_std(self.x, with_std=False)
        self.assertEqual((3.5, None), obs)

        obs = mean_and_std(self.x, ddof=2)
        npt.assert_almost_equal((3.5, 2.091650066), obs)

    def test_mean_and_std_no_mean_no_std(self):
        with npt.assert_raises(ValueError):
            mean_and_std(self.x, with_mean=False, with_std=False)

    def test_corr(self):
        obs = corr(self.small_mat)
        npt.assert_almost_equal(np.array([[1, 1, -0.94491118],
                                          [1, 1, -0.94491118],
                                          [-0.94491118, -0.94491118, 1]]),
                                obs)

    def test_corr_shape_mismatch(self):
        with npt.assert_raises(ValueError):
            corr(self.x, self.y)

    def test_e_matrix(self):
        E = e_matrix(self.matrix)
        expected_E = np.array([[-0.5, -2., -4.5],
                               [-8., -12.5, -18.]])
        npt.assert_almost_equal(E, expected_E)

    def test_f_matrix(self):
        F = f_matrix(self.matrix2)
        expected_F = np.zeros((3, 3))
        # Note that `test_make_F_matrix` in cogent is wrong
        npt.assert_almost_equal(F, expected_F)

    def test_e_matrix_inplace(self):
        E = _e_matrix_inplace(self.matrix)
        expected_E = np.array([[-0.5, -2., -4.5],
                               [-8., -12.5, -18.]])
        npt.assert_almost_equal(E, expected_E)

    def test_f_matrix_inplace(self):
        F = _f_matrix_inplace(self.matrix2)
        expected_F = np.zeros((3, 3))
        npt.assert_almost_equal(F, expected_F)

    def test_center_distance_matrix_inplace(self):
        dm_expected = f_matrix(e_matrix(self.small_mat))

        # make copy of matrix to test inplace centering
        matrix_copy = copy.deepcopy(self.small_mat)
        dm_centered = center_distance_matrix(matrix_copy, inplace=False)

        # ensure that matrix_copy was NOT modified inplace
        self.assertTrue(np.array_equal(matrix_copy, self.small_mat))

        # and ensure that the result of centering was correct
        npt.assert_almost_equal(dm_expected, dm_centered)

        # next, sort same matrix inplace
        matrix_copy2 = copy.deepcopy(self.small_mat)
        dm_centered_inp = center_distance_matrix(matrix_copy2, inplace=True)

        # and ensure that the result of inplace centering was correct
        npt.assert_almost_equal(dm_expected, dm_centered_inp)


if __name__ == '__main__':
    main()
