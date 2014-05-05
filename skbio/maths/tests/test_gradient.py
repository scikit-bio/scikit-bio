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
import numpy.testing as npt

import pandas as pd
import pandas.util.testing as pdt

from skbio.maths.gradient import (BaseVectors, AverageVectors,
                                  TrajectoryVectors, DifferenceVectors,
                                  WindowDifferenceVectors, VectorResults)


class BaseTests(TestCase):
    """"""
    def setUp(self):
        """"""
        coord_data = {
            'PC.636': np.array([-0.212230626531, 0.216034194368, 0.03532727349,
                                -0.254450494129, -0.0687468542543,
                                0.231895596562, 0.00496549154314,
                                -0.0026246871695, 9.73837390723e-10]),
            'PC.635': np.array([-0.277487312135, -0.0295483215975,
                                -0.0744173437992, 0.0957182357964,
                                0.204714844022, -0.0055407341857,
                                -0.190287966833, 0.16307126638,
                                9.73837390723e-10]),
            'PC.356': np.array([0.220886492631, 0.0874848360559,
                                -0.351990132198, -0.00316535032886,
                                0.114635191853, -0.00019194106125,
                                0.188557853937, 0.030002427212,
                                9.73837390723e-10]),
            'PC.481': np.array([0.0308923744062, -0.0446295973489,
                                0.133996451689, 0.29318228566, -0.167812539312,
                                0.130996149793, 0.113551017379, 0.109987942454,
                                9.73837390723e-10]),
            'PC.354': np.array([0.27616778138, -0.0341866951102,
                                0.0633000238256, 0.100446653327,
                                0.123802521199, 0.1285839664, -0.132852841046,
                                -0.217514322505, 9.73837390723e-10]),
            'PC.593': np.array([0.202458130052, -0.115216120518,
                                0.301820871723, -0.18300251046, 0.136208248567,
                                -0.0989435556722, 0.0927738484879,
                                0.0909429797672, 9.73837390723e-10]),
            'PC.355': np.array([0.236467470907, 0.21863434374,
                                -0.0301637746424, -0.0225473129718,
                                -0.205287183891, -0.180224615141,
                                -0.165277751908, 0.0411933458557,
                                9.73837390723e-10]),
            'PC.607': np.array([-0.105517545144, -0.41405687433,
                                -0.150073017617, -0.116066751485,
                                -0.158763393475, -0.0223918378516,
                                -0.0263068046112, -0.0501209518091,
                                9.73837390723e-10]),
            'PC.634': np.array([-0.371636765565, 0.115484234741,
                                0.0721996475289, 0.0898852445906,
                                0.0212491652909, -0.184183028843,
                                0.114877153051, -0.164938000185,
                                9.73837390723e-10])
            }
        self.coords = pd.DataFrame.from_dict(coord_data, orient='index')

        coord_data = {
            'PC.636': np.array([-0.212230626531, 0.216034194368,
                                0.03532727349]),
            'PC.635': np.array([-0.277487312135, -0.0295483215975,
                                -0.0744173437992]),
            'PC.356': np.array([0.220886492631, 0.0874848360559,
                                -0.351990132198]),
            'PC.481': np.array([0.0308923744062, -0.0446295973489,
                                0.133996451689]),
            'PC.354': np.array([0.27616778138, -0.0341866951102,
                                0.0633000238256]),
            'PC.593': np.array([0.202458130052, -0.115216120518,
                                0.301820871723]),
            'PC.355': np.array([0.236467470907, 0.21863434374,
                                -0.0301637746424]),
            'PC.607': np.array([-0.105517545144, -0.41405687433,
                                -0.150073017617]),
            'PC.634': np.array([-0.371636765565, 0.115484234741,
                                0.0721996475289])
            }
        self.coords_3axes = pd.DataFrame.from_dict(coord_data, orient='index')

        metamap = {'PC.354': {'Treatment': 'Control',
                              'DOB': '20061218',
                              'Weight': '60',
                              'Description': 'Control_mouse_I.D._354'},
                   'PC.355': {'Treatment': 'Control',
                              'DOB': '20061218',
                              'Weight': '55',
                              'Description': 'Control_mouse_I.D._355'},
                   'PC.356': {'Treatment': 'Control',
                              'DOB': '20061126',
                              'Weight': '50',
                              'Description': 'Control_mouse_I.D._356'},
                   'PC.481': {'Treatment': 'Control',
                              'DOB': '20070314',
                              'Weight': '52',
                              'Description': 'Control_mouse_I.D._481'},
                   'PC.593': {'Treatment': 'Control',
                              'DOB': '20071210',
                              'Weight': '57',
                              'Description': 'Control_mouse_I.D._593'},
                   'PC.607': {'Treatment': 'Fast',
                              'DOB': '20071112',
                              'Weight': '65',
                              'Description': 'Fasting_mouse_I.D._607'},
                   'PC.634': {'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Weight': '68',
                              'Description': 'Fasting_mouse_I.D._634'},
                   'PC.635': {'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Weight': '70',
                              'Description': 'Fasting_mouse_I.D._635'},
                   'PC.636': {'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Weight': '72',
                              'Description': 'Fasting_mouse_I.D._636'}}
        self.metamap = pd.DataFrame.from_dict(metamap, orient='index')

        self.eigvals = np.array([0.494208869204, 0.304213438411,
                                 0.272385344185, 0.22551158501, 0.18961649413,
                                 0.164195653585, 0.152144472132,
                                 0.126593271547, 8.53523337213e-18])


class BaseVectorsTests(BaseTests):
    """"""
    def test_init(self):
        """Correctly initializes the class attributes"""
        # Note self._groups is tested on test_make_groups
        # so we are not testing it here

        # Test with weighted = False
        bv = BaseVectors(self.coords, self.eigvals, self.metamap)

        pdt.assert_frame_equal(bv._coords, self.coords_3axes)
        npt.assert_equal(bv._eigvals, self.eigvals)
        pdt.assert_frame_equal(bv._metamap, self.metamap)
        self.assertTrue(bv._weighting_vector is None)

        # Test with weighted = True
        bv = BaseVectors(self.coords, self.eigvals, self.metamap,
                         sort_category='Weight', weighted=True)

        pdt.assert_frame_equal(bv._coords, self.coords_3axes)
        npt.assert_equal(bv._eigvals, self.eigvals)
        pdt.assert_frame_equal(bv._metamap, self.metamap)
        exp_weighting_vector = pd.Series(
            np.array([60, 55, 50, 52, 57, 65, 68, 70, 72]),
            ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593', 'PC.607',
             'PC.634', 'PC.635', 'PC.636']
            ).astype(np.float64)
        pdt.assert_series_equal(bv._weighting_vector, exp_weighting_vector)

    def test_init_error(self):
        """Raises an error with bad parameters"""
        # Raises ValueError if any category in vector_categories is not
        # present in metamap
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.eigvals, self.metamap,
                        vector_categories=['foo'])
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.eigvals, self.metamap,
                        vector_categories=['Weight', 'Treatment', 'foo'])

        # Raises ValueError if sort_category is not present in metamap
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.eigvals, self.metamap,
                        sort_category='foo')

        # Raises ValueError if weighted == True and sort_category == None
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.eigvals, self.metamap,
                        weighted=True)

        # Raises ValueError if weighted == True and the values under
        # sort_category are not numerical
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.eigvals, self.metamap,
                        sort_category='Treatment', weighted=True)

        # Raises ValueError if axes > len(eigvals)
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.eigvals, self.metamap, axes=10)

        # Raises ValueError if axes < 0
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.eigvals, self.metamap, axes=-1)

    def test_normalize_samples(self):
        """Correctly normalizes the samples between coords and metamap"""
        coord_data = {
            'PC.636': np.array([-0.212230626531, 0.216034194368,
                                0.03532727349]),
            'PC.635': np.array([-0.277487312135, -0.0295483215975,
                                -0.0744173437992]),
            'PC.355': np.array([0.236467470907, 0.21863434374,
                                -0.0301637746424]),
            'PC.607': np.array([-0.105517545144, -0.41405687433,
                                -0.150073017617]),
            'PC.634': np.array([-0.371636765565, 0.115484234741,
                                0.0721996475289])
            }
        subset_coords = pd.DataFrame.from_dict(coord_data, orient='index')

        metamap = {'PC.355': {'Treatment': 'Control',
                              'DOB': '20061218',
                              'Weight': '55',
                              'Description': 'Control_mouse_I.D._355'},
                   'PC.607': {'Treatment': 'Fast',
                              'DOB': '20071112',
                              'Weight': '65',
                              'Description': 'Fasting_mouse_I.D._607'},
                   'PC.634': {'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Weight': '68',
                              'Description': 'Fasting_mouse_I.D._634'},
                   'PC.635': {'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Weight': '70',
                              'Description': 'Fasting_mouse_I.D._635'},
                   'PC.636': {'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Weight': '72',
                              'Description': 'Fasting_mouse_I.D._636'}}
        subset_metamap = pd.DataFrame.from_dict(metamap, orient='index')

        # Takes a subset from metamap
        bv = BaseVectors(subset_coords, self.eigvals, self.metamap)
        pdt.assert_frame_equal(bv._coords.sort(axis=0),
                               subset_coords.sort(axis=0))
        pdt.assert_frame_equal(bv._metamap.sort(axis=0),
                               subset_metamap.sort(axis=0))

        # Takes a subset from coords
        bv = BaseVectors(self.coords, self.eigvals, subset_metamap)
        pdt.assert_frame_equal(bv._coords.sort(axis=0),
                               subset_coords.sort(axis=0))
        pdt.assert_frame_equal(bv._metamap.sort(axis=0),
                               subset_metamap.sort(axis=0))

        # Takes a subset from metamap and coords at the same time
        coord_data = {
            'PC.636': np.array([-0.212230626531, 0.216034194368,
                                0.03532727349]),
            'PC.635': np.array([-0.277487312135, -0.0295483215975,
                                -0.0744173437992]),
            'PC.355': np.array([0.236467470907, 0.21863434374,
                                -0.0301637746424])
            }
        subset_coords = pd.DataFrame.from_dict(coord_data, orient='index')

        metamap = {'PC.355': {'Treatment': 'Control',
                              'DOB': '20061218',
                              'Weight': '55',
                              'Description': 'Control_mouse_I.D._355'},
                   'PC.607': {'Treatment': 'Fast',
                              'DOB': '20071112',
                              'Weight': '65',
                              'Description': 'Fasting_mouse_I.D._607'},
                   'PC.634': {'Treatment': 'Fast',
                              'DOB': '20080116',
                              'Weight': '68',
                              'Description': 'Fasting_mouse_I.D._634'}}
        subset_metamap = pd.DataFrame.from_dict(metamap, orient='index')

        bv = BaseVectors(subset_coords, self.eigvals, subset_metamap)
        exp_coords = pd.DataFrame.from_dict(
            {'PC.355': np.array([0.236467470907, 0.21863434374,
                                 -0.0301637746424])},
            orient='index')
        pdt.assert_frame_equal(bv._coords.sort(axis=0),
                               exp_coords.sort(axis=0))
        exp_metamap = pd.DataFrame.from_dict(
            {'PC.355': {'Treatment': 'Control',
                        'DOB': '20061218',
                        'Weight': '55',
                        'Description': 'Control_mouse_I.D._355'}},
            orient='index')
        pdt.assert_frame_equal(bv._metamap.sort(axis=0),
                               exp_metamap.sort(axis=0))

    def test_normalize_samples_error(self):
        """Raises an error if coords and metamap does not have samples in
        common"""
        error_metamap = pd.DataFrame.from_dict(
            {'Foo': {'Treatment': 'Control',
                     'DOB': '20061218',
                     'Weight': '55',
                     'Description': 'Control_mouse_I.D._355'},
             'Bar': {'Treatment': 'Fast',
                     'DOB': '20071112',
                     'Weight': '65',
                     'Description': 'Fasting_mouse_I.D._607'}},
            orient='index')
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.eigvals, error_metamap)

    def test_make_groups(self):
        """Correctly generates the groups for vector_categories"""
        # Test with all categories
        bv = BaseVectors(self.coords, self.eigvals, self.metamap)
        exp_groups = {'Treatment': {'Control': [('PC.354', 'PC.354'),
                                                ('PC.355', 'PC.355'),
                                                ('PC.356', 'PC.356'),
                                                ('PC.481', 'PC.481'),
                                                ('PC.593', 'PC.593')],
                                    'Fast': [('PC.607', 'PC.607'),
                                             ('PC.634', 'PC.634'),
                                             ('PC.635', 'PC.635'),
                                             ('PC.636', 'PC.636')]},
                      'DOB': {'20061218': [('PC.354', 'PC.354'),
                                           ('PC.355', 'PC.355')],
                              '20061126': [('PC.356', 'PC.356')],
                              '20070314': [('PC.481', 'PC.481')],
                              '20071210': [('PC.593', 'PC.593')],
                              '20071112': [('PC.607', 'PC.607')],
                              '20080116': [('PC.634', 'PC.634'),
                                           ('PC.635', 'PC.635'),
                                           ('PC.636', 'PC.636')]},
                      'Weight': {'60': [('PC.354', 'PC.354')],
                                 '55': [('PC.355', 'PC.355')],
                                 '50': [('PC.356', 'PC.356')],
                                 '52': [('PC.481', 'PC.481')],
                                 '57': [('PC.593', 'PC.593')],
                                 '65': [('PC.607', 'PC.607')],
                                 '68': [('PC.634', 'PC.634')],
                                 '70': [('PC.635', 'PC.635')],
                                 '72': [('PC.636', 'PC.636')]},
                      'Description': {'Control_mouse_I.D._354': [('PC.354',
                                                                  'PC.354')],
                                      'Control_mouse_I.D._355': [('PC.355',
                                                                  'PC.355')],
                                      'Control_mouse_I.D._356': [('PC.356',
                                                                  'PC.356')],
                                      'Control_mouse_I.D._481': [('PC.481',
                                                                  'PC.481')],
                                      'Control_mouse_I.D._593': [('PC.593',
                                                                  'PC.593')],
                                      'Fasting_mouse_I.D._607': [('PC.607',
                                                                  'PC.607')],
                                      'Fasting_mouse_I.D._634': [('PC.634',
                                                                  'PC.634')],
                                      'Fasting_mouse_I.D._635': [('PC.635',
                                                                  'PC.635')],
                                      'Fasting_mouse_I.D._636': [('PC.636',
                                                                  'PC.636')]}}
        self.assertEqual(bv._groups, exp_groups)

        # Test with user-defined categories
        bv = BaseVectors(self.coords, self.eigvals, self.metamap,
                         vector_categories=['Treatment', 'DOB'])
        exp_groups = {'Treatment': {'Control': [('PC.354', 'PC.354'),
                                                ('PC.355', 'PC.355'),
                                                ('PC.356', 'PC.356'),
                                                ('PC.481', 'PC.481'),
                                                ('PC.593', 'PC.593')],
                                    'Fast': [('PC.607', 'PC.607'),
                                             ('PC.634', 'PC.634'),
                                             ('PC.635', 'PC.635'),
                                             ('PC.636', 'PC.636')]},
                      'DOB': {'20061218': [('PC.354', 'PC.354'),
                                           ('PC.355', 'PC.355')],
                              '20061126': [('PC.356', 'PC.356')],
                              '20070314': [('PC.481', 'PC.481')],
                              '20071210': [('PC.593', 'PC.593')],
                              '20071112': [('PC.607', 'PC.607')],
                              '20080116': [('PC.634', 'PC.634'),
                                           ('PC.635', 'PC.635'),
                                           ('PC.636', 'PC.636')]}}
        self.assertEqual(bv._groups, exp_groups)

    def test_get_vectors(self):
        """"""
        pass

    def test_get_subgroup_vectors(self):
        """"""
        pass

    def test_compute_vector_results(self):
        """"""
        pass

    def test_weight_by_vector(self):
        """"""
        pass


# class AverageVectorsTests(BaseTests):
#     """"""
#     def test_results(self):
#         av = AverageVectors(self.coord_dict, self.eigvals, self.ids)
#         obs = av.results()
#         exp_vector = [6.9954747524, 1.5180408981, 7.4608959440]
#         exp_calc = {'avg': 5.3248038648}

#         self.assertEqual(type(obs), VectorResults)
#         assert_almost_equal(obs.vector, exp_vector)
#         self.assertEqual(obs.calc.keys(), exp_calc.keys())
#         for key in obs.calc:
#             assert_almost_equal(obs.calc[key], exp_calc[key])
#         self.assertEqual(obs.message, None)


# class TrajectoryVectorsTests(BaseTests):
#     def test_results(self):
#         tv = TrajectoryVectors(self.coord_dict, self.eigvals, self.ids)
#         obs = tv.results()
#         exp_vector = [14.3839779766, 6.8423140875]
#         exp_calc = {'trajectory': 15.9284677387}

#         self.assertEqual(type(obs), VectorResults)
#         assert_almost_equal(obs.vector, exp_vector)
#         self.assertEqual(obs.calc.keys(), exp_calc.keys())
#         for key in obs.calc:
#             assert_almost_equal(obs.calc[key], exp_calc[key])
#         self.assertEqual(obs.message, None)


# class DifferenceVectorsTests(BaseTests):
#     def test_results(self):
#         dv = DifferenceVectors(self.coord_dict, self.eigvals, self.ids)
#         obs = dv.results()
#         exp_vector = np.array([-7.54166389])
#         exp_calc = {'mean': [-7.541663889], 'std': [0.0]}

#         self.assertEqual(type(obs), VectorResults)
#         assert_almost_equal(obs.vector, exp_vector)
#         self.assertEqual(obs.calc.keys(), exp_calc.keys())
#         for key in obs.calc:
#             assert_almost_equal(obs.calc[key], exp_calc[key])
#         self.assertEqual(obs.message, None)


# class WindowDifferenceVectorsTests(BaseTests):
#     def test_results(self):
#         wdv = WindowDifferenceVectors(self.coord_dict, self.eigvals,
#                                       self.ids, 1)
#         obs = wdv.results()
#         exp_vector = [-7.5416638890, 0.0]
#         exp_calc = {'std': 3.7708319445,  'mean': -3.7708319445}

#         self.assertEqual(type(obs), VectorResults)
#         assert_almost_equal(obs.vector, exp_vector)
#         self.assertEqual(obs.calc.keys(), exp_calc.keys())
#         for key in obs.calc:
#             assert_almost_equal(obs.calc[key], exp_calc[key])
#         self.assertEqual(obs.message, None)

if __name__ == '__main__':
    main()
