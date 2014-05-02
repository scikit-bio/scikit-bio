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

from pandas import DataFrame

from skbio.maths.stats.ordination import OrdinationResults
from skbio.maths.gradient import (make_groups, weight_by_vector, windowed_diff,
                                  AverageVectors, TrajectoryVectors,
                                  DifferenceVectors, WindowDifferenceVectors,
                                  VectorResults)


class MakeGroupsTests(TestCase):
    def setUp(self):
        eigvals = np.array([0.512367260461, 0.300719094427,
                            0.267912066004, 0.208988681078, 0.19169895326,
                            0.16054234528, 0.15017695712, 0.122457748167,
                            0.0])
        site = np.array([[-0.258465461183, 0.173999546883, 0.0382875792552,
                          -0.19447750562, 0.0831176020844, 0.262430333201,
                          -0.0231636392235, -0.0184794039581, 0.0],
                         [-0.271001135391, -0.0185951319063, -0.0864841926349,
                          0.118064245315, -0.198808358437, -0.0211723599535,
                          -0.191024027565, 0.155646592377, 0.0],
                         [0.235077898175, 0.0962519254489, -0.345792726714,
                          -0.00320862577619, -0.0963777675519, 0.0457025386953,
                          0.185472813286, 0.0404093971793, 0.0],
                         [0.0261407664325, -0.0111459676533, 0.147660603015,
                          0.29087660853, 0.203945472801, 0.0619712384758,
                          0.101641328709, 0.105690998719, 0.0],
                         [0.285007552283, -0.0192549888483, 0.0623263375385,
                          0.138126799852, -0.104798602423, 0.0951720730628,
                          -0.129636097542, -0.220687170372, 0.0],
                         [0.204636326241, -0.139361150932, 0.291513819623,
                          -0.181566786821, -0.159580132715, -0.0246412130162,
                          0.0866252404441, 0.0996221476871, 0.0],
                         [0.233482403212, 0.225257974068, -0.0188623096268,
                          -0.107729981831, 0.177108999572, -0.192905835151,
                          -0.149819471408, 0.0383549037465, 0.0],
                         [-0.0949631911323, -0.420974802495, -0.154869454869,
                          -0.0898427509281, 0.152618194488, -0.0334232691501,
                          -0.0251224777303, -0.0508988536409, 0.0],
                         [-0.359915158638, 0.113822595435, 0.0662203444138,
                          0.0297579972788, -0.0572254078183, -0.193133506163,
                          0.145026331031, -0.149658611738, 0.0]])
        prop_expl = np.array([0.267573832777, 0.15704469605, 0.139911863774,
                              0.109140272454, 0.100111048503, 0.0838401161912,
                              0.0784269939011, 0.0639511763509, 0.0])
        site_ids = ['PC.636', 'PC.635', 'PC.356', 'PC.481', 'PC.354',
                    'PC.593', 'PC.355', 'PC.607', 'PC.634']
        self.ord_res = OrdinationResults(eigvals=eigvals, site=site,
                                         proportion_explained=prop_expl,
                                         site_ids=site_ids)

        mapdata = {'PC.354': {'BarcodeSequence': 'AGCACGAGCCTA',
                              'Treatment': 'Control',
                              'DOB': '20061218',
                              'Description': '354'},
                   'PC.355': {'BarcodeSequence': 'AACTCGTCGATG',
                              'Treatment': 'Control',
                              'DOB': '20061218',
                              'Description': '355'},
                   'PC.356': {'BarcodeSequence': 'ACAGACCACTCA',
                              'Treatment': 'Control',
                              'DOB': '20061126',
                              'Description': '356'},
                   'PC.481': {'BarcodeSequence': 'ACCAGCGACTAG',
                              'Treatment': 'Control',
                              'DOB': '20070314',
                              'Description': '481'},
                   'PC.593': {'BarcodeSequence': 'AGCAGCACTTGT',
                              'Treatment': 'Control',
                              'DOB': '20071210',
                              'Description': '593'},
                   'PC.607': {'BarcodeSequence': 'AACTGTGCGTAC',
                              'Treatment': 'Fast',
                              'DOB': '20071112',
                              'Description': '607'},
                   'PC.634': {'BarcodeSequence': 'ACAGAGTCGGCT',
                              'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Description': '634'},
                   'PC.635': {'BarcodeSequence': 'ACCGCAGAGTCA',
                              'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Description': '635'},
                   'PC.636': {'BarcodeSequence': 'ACGGTGAGTGTC',
                              'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Description': '636'}}
        self.metamap = DataFrame.from_dict(mapdata, orient='index')

    def test_make_groups(self):
        """Correctly groups samples"""
        obs = make_groups(self.ord_res, self.metamap, 'Treatment')
        exp = {'Control': [('PC.354', 'PC.354'),
                           ('PC.355', 'PC.355'),
                           ('PC.356', 'PC.356'),
                           ('PC.481', 'PC.481'),
                           ('PC.593', 'PC.593')],
               'Fast': [('PC.607', 'PC.607'),
                        ('PC.634', 'PC.634'),
                        ('PC.635', 'PC.635'),
                        ('PC.636', 'PC.636')]}
        self.assertEqual(obs, exp)

    def test_make_groups_sorted(self):
        """Correctly groups and sorts samples"""
        obs = make_groups(self.ord_res, self.metamap, 'Treatment', 'DOB')
        exp = {'Control': [('20061126', 'PC.356'),
                           ('20061218', 'PC.354'),
                           ('20061218', 'PC.355'),
                           ('20070314', 'PC.481'),
                           ('20071210', 'PC.593')],
               'Fast': [('20071112', 'PC.607'),
                        ('20080116', 'PC.634'),
                        ('20080116', 'PC.635'),
                        ('20080116', 'PC.636')]}
        self.assertEqual(obs, exp)


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
