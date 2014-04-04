#!/usr/bin/env python
from __future__ import division

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from unittest import TestCase, main
import numpy as np
from numpy.testing import assert_almost_equal

from skbio.maths.gradient import (weight_by_vector, windowed_diff,
                                  AverageVectors, TrajectoryVectors,
                                  DifferenceVectors, WindowDifferenceVectors,
                                  VectorResults)


class GradientTests(TestCase):
    """"""
    def test_weight_by_vector_error(self):
        """Raises an error with erroneous inputs"""
        with self.assertRaises(ValueError):
            weight_by_vector([1, 2, 3, 4], [1, 2, 3])
        with self.assertRaises(TypeError):
            weight_by_vector(9, 1)
        with self.assertRaises(ValueError):
            weight_by_vector([1, 2, 3, 4], [1, 2, 3, 3])

    def test_weight_by_vector(self):
        """Test that the vector is weighted by the given vector"""
        # data for test_weight_by_vector
        in_vector_to_weight = np.array([1, 2, 3, 4, 5, 6, 7, 8])
        in_weighting_vector = np.array([1, 5, 8, 12, 45, 80, 85, 90])
        out_weighted_vector = np.array([1, 6.3571428571, 12.7142857142,
                                        12.7142857142, 1.9264069264,
                                        2.1795918367, 17.8, 20.3428571428])
        obs = weight_by_vector(in_vector_to_weight, in_weighting_vector)
        assert_almost_equal(obs, out_weighted_vector)

        in_flat_vector_to_weight = np.array([1, 2, 3, 4, 5, 6, 7, 8])
        in_flat_weighting_vector = np.array([1, 2, 3, 4, 5, 6, 7, 8])
        out_flat_weighted_vector = np.array([1, 2, 3, 4, 5, 6, 7, 8])
        obs = weight_by_vector(in_flat_vector_to_weight,
                               in_flat_weighting_vector)
        assert_almost_equal(obs, out_flat_weighted_vector)

        in_second_flat_vector_to_weight = np.array([2, 3, 4, 5, 6])
        in_second_flat_weighting_vector = np.array([25, 30, 35, 40, 45])
        out_second_flat_weighted_vector = np.array([2, 3, 4, 5, 6])
        obs = weight_by_vector(in_second_flat_vector_to_weight,
                               in_second_flat_weighting_vector)
        assert_almost_equal(obs, out_second_flat_weighted_vector)

        in_nd_array_vector_to_weight = np.array([[1, 2, 3], [2, 3, 4],
                                                 [5, 6, 7], [8, 9, 10]])
        in_nd_array_weighting_vector = np.array([1, 2, 3, 4])
        out_nd_array_flat_vector = np.array([[1, 2, 3], [2, 3, 4],
                                             [5, 6, 7], [8, 9, 10]])
        obs = weight_by_vector(in_nd_array_vector_to_weight,
                               in_nd_array_weighting_vector)
        assert_almost_equal(obs, out_nd_array_flat_vector)

    def test_windowed_diff_error(self):
        """Raises an error with erroneous inputs"""
        in_array_windowed_diff = np.array([24, 15, 28, 16, 28, 43, 12, 53])

        # check all the disallowed behaviors
        with self.assertRaises(ValueError):
            windowed_diff(in_array_windowed_diff, -10)
        with self.assertRaises(ValueError):
            windowed_diff(in_array_windowed_diff, 8)
        with self.assertRaises(ValueError):
            windowed_diff(in_array_windowed_diff, 50)

    def test_windowed_diff(self):
        """Test correct functioning of the modified first difference"""
        in_array_windowed_diff = np.array([24, 15, 28, 16, 28, 43, 12, 53])

        # checks for correct behavior with a 1-d vector
        obs = windowed_diff(in_array_windowed_diff, 3)
        exp = np.array([-4.3333333333, 9.0, 1.0, 11.6666666666, 8.0,
                        -3.6666666666, 41.0, 0.0])
        assert_almost_equal(obs, exp)

        obs = windowed_diff(in_array_windowed_diff, 1)
        exp = [-9.0, 13.0, -12.0, 12.0, 15.0, -31.0, 41.0, 0.0]
        assert_almost_equal(obs, exp)

        # checks for correct behaviour with an n-d vector
        in_nd_array_windowed_diff = np.array([[2, 5, 1, 6],
                                              [7, 1, 6, 3],
                                              [1, 8, 3, 5],
                                              [7, 3, 2, 6],
                                              [8, 1, 6, 1],
                                              [15, 5, 4, 1],
                                              [1, 5, 1, 2],
                                              [33, 5, 67, 12]])
        obs = windowed_diff(in_nd_array_windowed_diff, 2)
        exp = np.array([[2., -0.5,  3.5, -2.],
                       [-3., 4.5, -3.5, 2.5],
                       [6.5, -6.,  1., -1.5],
                       [4.5, 0., 3., -5.],
                       [0., 4., -3.5, 0.5],
                       [2., 0., 30., 6.],
                       [32., 0., 66., 10.],
                       [0., 0., 0., 0.]])
        assert_almost_equal(obs, exp)


class BaseTests(TestCase):
    """"""
    def setUp(self):
        """"""
        self.coord_header = ["Sample1", "Sample2", "Sample3"]
        self.coords = np.array([[-0.219044992, 0.079674486, 0.09233683],
                                [-0.042258081, 0.000204041, 0.024837603],
                                [0.080504323, -0.212014503, -0.088353435]])
        self.coord_dict = dict(zip(self.coord_header, self.coords))
        self.pct_var = np.array([25.00, 30.00, 35.00])
        self.eigvals = [i for i in reversed(self.pct_var)]
        self.ids = ["Sample1", "Sample2", "Sample3"]


class AverageVectorsTests(BaseTests):
    """"""
    def test_results(self):
        av = AverageVectors(self.coord_dict, self.eigvals, self.ids)
        obs = av.results()
        exp_vector = [6.9954747524, 1.5180408981, 7.4608959440]
        exp_calc = {'avg': 5.3248038648}

        self.assertEqual(type(obs), VectorResults)
        assert_almost_equal(obs.vector, exp_vector)
        self.assertEqual(obs.calc.keys(), exp_calc.keys())
        for key in obs.calc:
            assert_almost_equal(obs.calc[key], exp_calc[key])
        self.assertEqual(obs.message, None)


class TrajectoryVectorsTests(BaseTests):
    def test_results(self):
        tv = TrajectoryVectors(self.coord_dict, self.eigvals, self.ids)
        obs = tv.results()
        exp_vector = [14.3839779766, 6.8423140875]
        exp_calc = {'trajectory': 15.9284677387}

        self.assertEqual(type(obs), VectorResults)
        assert_almost_equal(obs.vector, exp_vector)
        self.assertEqual(obs.calc.keys(), exp_calc.keys())
        for key in obs.calc:
            assert_almost_equal(obs.calc[key], exp_calc[key])
        self.assertEqual(obs.message, None)


class DifferenceVectorsTests(BaseTests):
    def test_results(self):
        dv = DifferenceVectors(self.coord_dict, self.eigvals, self.ids)
        obs = dv.results()
        exp_vector = np.array([-7.54166389])
        exp_calc = {'mean': [-7.541663889], 'std': [0.0]}

        self.assertEqual(type(obs), VectorResults)
        assert_almost_equal(obs.vector, exp_vector)
        self.assertEqual(obs.calc.keys(), exp_calc.keys())
        for key in obs.calc:
            assert_almost_equal(obs.calc[key], exp_calc[key])
        self.assertEqual(obs.message, None)


class WindowDifferenceVectorsTests(BaseTests):
    def test_results(self):
        wdv = WindowDifferenceVectors(self.coord_dict, self.eigvals,
                                      self.ids, 1)
        obs = wdv.results()
        exp_vector = [-7.5416638890, 0.0]
        exp_calc = {'std': 3.7708319445,  'mean': -3.7708319445}

        self.assertEqual(type(obs), VectorResults)
        assert_almost_equal(obs.vector, exp_vector)
        self.assertEqual(obs.calc.keys(), exp_calc.keys())
        for key in obs.calc:
            assert_almost_equal(obs.calc[key], exp_calc[key])
        self.assertEqual(obs.message, None)

if __name__ == '__main__':
    main()
