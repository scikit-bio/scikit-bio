# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio import DistanceMatrix

from skbio.stats.distance import mrm


class RegressionTests(TestCase):
    def setUp(self):
        self.x1 = DistanceMatrix([[0, 1, 2.1],
                                  [1, 0, 2.9],
                                  [2.1, 2.9, 0]])
        self.y1 = DistanceMatrix([[0, 2, 4],
                                  [2, 0, 6],
                                  [4, 6, 0]])

        self.x2 = DistanceMatrix([[0, 1, 2.1, 4.4, 5.3],
                                  [1, 0, 2.9, 6.2, 7.1],
                                  [2.1, 2.9, 0, 7.9, 8.8],
                                  [4.4, 6.2, 7.9, 0, 9.7],
                                  [5.3, 7.1, 8.8, 9.7, 0]])
        self.y2 = DistanceMatrix([[0, 2, 4, 9, 10],
                                  [2, 0, 6, 12, 14],
                                  [4, 6, 0, 16, 18],
                                  [9, 12, 16, 0, 20],
                                  [10, 14, 18, 20, 0]])
        self.x3 = DistanceMatrix(self.x2.data*0.5)

        # self.x2 with 4 and 5 swapped
        self.x4 = DistanceMatrix([[0, 1, 2.1, 5.3, 4.4],
                                  [1, 0, 2.9, 7.1, 6.2],
                                  [2.1, 2.9, 0, 8.8, 7.9],
                                  [5.3, 7.1, 8.8, 0, 9.7],
                                  [4.4, 6.2, 7.9, 9.7, 0]],
                                  ['1', '2', '3', '5', '4'])
        self.y3 = DistanceMatrix(self.y2.data,
                                 ['1', '2', '3', '4', '5'])
        self.x5 = DistanceMatrix(self.x3.data,
                                 ['1', '2', '3', '4', '5'])
        self.x6 = DistanceMatrix([[0, 1, 2.1, 4.4, 5.3, 1],
                                  [1, 0, 2.9, 6.2, 7.1, 1],
                                  [2.1, 2.9, 0, 7.9, 8.8, 1],
                                  [4.4, 6.2, 7.9, 0, 9.7, 1],
                                  [5.3, 7.1, 8.8, 9.7, 0, 1],
                                  [1, 1, 1, 1, 1, 0]],
                                  ['1', '2', '3', '5', '4', '6'])

    def test1(self):
        np.random.seed(0)
        B, T, pvals, F, model_pval, R2 = mrm(self.y1, self.x1,
                                                    permutations=10)
        npt.assert_allclose(B.values,
                            np.array([-0.17582418,  2.08791209]))
        npt.assert_allclose(T.values,
                            np.array([ -0.43039376,  10.96965511]))
        npt.assert_allclose(pvals.values,
                            np.array([ 0.45454545,  0.27272727]))
        self.assertAlmostEqual(F, 120.33333333340981)
        self.assertAlmostEqual(model_pval, 0.2727272727272727)
        self.assertAlmostEqual(R2, 0.9917582417582471)

        B, T, pvals, F, model_pval, R2 = mrm(self.y1, self.x1,
                                                    permutations=10,
                                                    random_state=0)
        npt.assert_allclose(B.values,
                            np.array([-0.17582418,  2.08791209]))
        npt.assert_allclose(T.values,
                            np.array([ -0.43039376,  10.96965511]))
        npt.assert_allclose(pvals.values,
                            np.array([ 0.45454545,  0.27272727]))
        self.assertAlmostEqual(F, 120.33333333340981)
        self.assertAlmostEqual(model_pval, 0.2727272727272727)
        self.assertAlmostEqual(R2, 0.9917582417582471)

        st = np.random.RandomState(seed=0)
        B, T, pvals, F, model_pval, R2 = mrm(self.y1, self.x1,
                                                    permutations=10,
                                                    random_state=st)
        npt.assert_allclose(B.values,
                            np.array([-0.17582418,  2.08791209]))
        npt.assert_allclose(T.values,
                            np.array([ -0.43039376,  10.96965511]))
        npt.assert_allclose(pvals.values,
                            np.array([ 0.45454545,  0.27272727]))
        self.assertAlmostEqual(F, 120.33333333340981)
        self.assertAlmostEqual(model_pval, 0.2727272727272727)
        self.assertAlmostEqual(R2, 0.9917582417582471)

    def test2(self):

        B, T, pvals, F, model_pval, R2 = mrm(self.y2, self.x2,
                                                    permutations=10,
                                                    random_state=0)
        npt.assert_allclose(B.values,
                            np.array([-0.25088147,  2.04889557]))
        npt.assert_allclose(T.values,
                            np.array([ -0.98072414,  49.6360995 ]))
        npt.assert_allclose(pvals.values,
                            np.array([ 0.90909091,  0.36363636]))
        self.assertAlmostEqual(F, 2463.7423732163857)
        self.assertAlmostEqual(model_pval, 0.36363636363636365)
        self.assertAlmostEqual(R2, 0.99676341673522)

        B, T, pvals, F, model_pval, R2 = mrm(self.y2, self.x2, self.x3,
                                                    permutations=100,
                                                    random_state=0)
        npt.assert_allclose(B.values,
                            np.array([-0.25088147, 1.63911646, 0.81955823]))
        npt.assert_allclose(T.values,
                            np.array([-0.91738343, 46.43031958,
                                      46.43031958]))
        npt.assert_allclose(pvals.values,
                            np.array([0.4950495, 0.14851485,
                                      0.14851485]))

        self.assertItemsEqual(B.index,
                              ['intercept'] + range(2))
        self.assertItemsEqual(T.index,
                              ['intercept'] + range(2))
        self.assertItemsEqual(pvals.index,
                              ['intercept'] + range(2))

        self.assertAlmostEqual(F, 1077.8872882816077)
        self.assertAlmostEqual(model_pval, 0.1485148514851485)
        self.assertAlmostEqual(R2, 0.9967634167352184)

    # Tests for reordering of ids in distance matrices
    def test3(self):
        B, T, pvals, F, model_pval, R2 = mrm(self.y3, self.x4, self.x5,
                                             permutations=100,
                                             random_state=0)
        npt.assert_allclose(B.values,
                            np.array([-0.25088147, 1.63911646, 0.81955823]))
        npt.assert_allclose(T.values,
                            np.array([-0.91738343, 46.43031958,
                                      46.43031958]))
        npt.assert_allclose(pvals.values,
                            np.array([0.4950495, 0.14851485,
                                      0.14851485]))
        self.assertItemsEqual(B.index,
                              ['intercept'] + range(2))
        self.assertItemsEqual(T.index,
                              ['intercept'] + range(2))
        self.assertItemsEqual(pvals.index,
                              ['intercept'] + range(2))

        self.assertAlmostEqual(F, 1077.8872882816077)
        self.assertAlmostEqual(model_pval, 0.1485148514851485)
        self.assertAlmostEqual(R2, 0.9967634167352184)

    # Tests for bad inputs (i.e. bad distance matrices)
    def test_bad(self):
        with self.assertRaises(ValueError):
            mrm(self.y3, self.x4, self.x2,
                       permutations=100,
                       random_state=0)

        with self.assertRaises(ValueError):
            mrm(self.y3, self.x4, self.x5,
                labels=['1'],
                permutations=100,
                random_state=0)

        with self.assertRaises(ValueError):
            mrm(self.y3, self.x4, self.x5,
                labels=['1', '1'],
                permutations=100,
                random_state=0)

    def test_distance_matrix_instances_with_lookup(self):
        self.x1.ids = ('a', 'b', 'c')
        self.y1.ids = ('d', 'e', 'f')
        lookup = {'a': 'A', 'b': 'B', 'c': 'C',
                  'd': 'A', 'e': 'B', 'f': 'C'}

        B, T, pvals, F, model_pval, R2 = mrm(self.y1, self.x1,
                                             permutations=10,
                                             lookup=lookup,
                                             random_state=0)
        npt.assert_allclose(B.values,
                            np.array([-0.17582418,  2.08791209]))
        npt.assert_allclose(T.values,
                            np.array([ -0.43039376,  10.96965511]))
        npt.assert_allclose(pvals.values,
                            np.array([ 0.45454545,  0.27272727]))

        self.assertItemsEqual(B.index,
                              ['intercept', 0])
        self.assertItemsEqual(T.index,
                              ['intercept', 0])
        self.assertItemsEqual(pvals.index,
                              ['intercept', 0])
        self.assertAlmostEqual(F, 120.33333333340981)
        self.assertAlmostEqual(model_pval, 0.2727272727272727)
        self.assertAlmostEqual(R2, 0.9917582417582471)

    def test_strict(self):
        B, T, pvals, F, model_pval, R2 = mrm(self.y3, self.x6, self.x5,
                                             permutations=100,
                                             strict=False,
                                             random_state=0)
        npt.assert_allclose(B.values,
                            np.array([-0.25088147, 1.63911646, 0.81955823]))
        npt.assert_allclose(T.values,
                            np.array([-0.91738343, 46.43031958,
                                      46.43031958]))
        npt.assert_allclose(pvals.values,
                            np.array([0.4950495, 0.14851485,
                                      0.14851485]))
        self.assertItemsEqual(B.index,
                              ['intercept'] + range(2))
        self.assertItemsEqual(T.index,
                              ['intercept'] + range(2))
        self.assertItemsEqual(pvals.index,
                              ['intercept'] + range(2))

        self.assertAlmostEqual(F, 1077.8872882816077)
        self.assertAlmostEqual(model_pval, 0.1485148514851485)
        self.assertAlmostEqual(R2, 0.9967634167352184)

if __name__ == "__main__":
    main()
