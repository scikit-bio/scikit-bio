#!/usr/bin/env python

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from __future__ import division
from future.builtins import zip

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

import pandas as pd
import pandas.util.testing as pdt

from skbio.maths.gradient import (BaseVectors, AverageVectors,
                                  TrajectoryVectors, DifferenceVectors,
                                  WindowDifferenceVectors, GroupResults,
                                  CategoryResults, VectorsResults)


class BaseTests(TestCase):
    def setUp(self):
        """Initializes some data for testing"""
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

        self.prop_expl = np.array([25.6216900347, 15.7715955926,
                                   14.1215046787, 11.6913885817, 9.83044890697,
                                   8.51253468595, 7.88775505332, 6.56308246609,
                                   4.42499350906e-16])

    # This function makes the comparisons between the results classes easier
    def assert_group_results_almost_equal(self, obs, exp):
        """Tests that obs and exp are almost equal"""
        self.assertEqual(obs.name, exp.name)
        npt.assert_almost_equal(obs.vector, exp.vector)
        npt.assert_almost_equal(obs.mean, exp.mean)
        self.assertEqual(obs.info.keys(), exp.info.keys())
        for key in obs.info:
            npt.assert_almost_equal(obs.info[key], exp.info[key])
        self.assertEqual(obs.message, exp.message)

    def assert_category_results_almost_equal(self, obs, exp):
        """Tests that obs and exp are almost equal"""
        self.assertEqual(obs.category, exp.category)
        npt.assert_almost_equal(obs.probability, exp.probability)

        for o, e in zip(sorted(obs.groups), sorted(exp.groups)):
            self.assert_group_results_almost_equal(o, e)

    def assert_vectors_results_almost_equal(self, obs, exp):
        """Tests that obs and exp are almost equal"""
        self.assertEqual(obs.algorithm, exp.algorithm)
        self.assertEqual(obs.weighted, exp.weighted)

        for o, e in zip(sorted(obs.categories), sorted(exp.categories)):
            self.assert_category_results_almost_equal(o, e)


class BaseVectorsTests(BaseTests):
    def test_init(self):
        """Correctly initializes the class attributes"""
        # Note self._groups is tested on test_make_groups
        # so we are not testing it here

        # Test with weighted = False
        bv = BaseVectors(self.coords, self.prop_expl, self.metamap)

        pdt.assert_frame_equal(bv._coords, self.coords_3axes)
        exp_prop_expl = np.array([25.6216900347, 15.7715955926,
                                  14.1215046787])
        npt.assert_equal(bv._prop_expl, exp_prop_expl)
        pdt.assert_frame_equal(bv._metamap, self.metamap)
        self.assertTrue(bv._weighting_vector is None)
        self.assertFalse(bv._weighted)

        # Test with weighted = True
        bv = BaseVectors(self.coords, self.prop_expl, self.metamap,
                         sort_category='Weight', weighted=True)

        pdt.assert_frame_equal(bv._coords, self.coords_3axes)
        npt.assert_equal(bv._prop_expl, exp_prop_expl)
        pdt.assert_frame_equal(bv._metamap, self.metamap)
        exp_weighting_vector = pd.Series(
            np.array([60, 55, 50, 52, 57, 65, 68, 70, 72]),
            ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593', 'PC.607',
             'PC.634', 'PC.635', 'PC.636']
            ).astype(np.float64)
        pdt.assert_series_equal(bv._weighting_vector, exp_weighting_vector)
        self.assertTrue(bv._weighted)

    def test_init_error(self):
        """Raises an error with erroneous inputs"""
        # Raises ValueError if any category in vector_categories is not
        # present in metamap
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.prop_expl, self.metamap,
                        vector_categories=['foo'])
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.prop_expl, self.metamap,
                        vector_categories=['Weight', 'Treatment', 'foo'])

        # Raises ValueError if sort_category is not present in metamap
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.prop_expl, self.metamap,
                        sort_category='foo')

        # Raises ValueError if weighted == True and sort_category == None
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.prop_expl, self.metamap,
                        weighted=True)

        # Raises ValueError if weighted == True and the values under
        # sort_category are not numerical
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.prop_expl, self.metamap,
                        sort_category='Treatment', weighted=True)

        # Raises ValueError if axes > len(prop_expl)
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.prop_expl, self.metamap, axes=10)

        # Raises ValueError if axes < 0
        with self.assertRaises(ValueError):
            BaseVectors(self.coords, self.prop_expl, self.metamap, axes=-1)

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
        bv = BaseVectors(subset_coords, self.prop_expl, self.metamap)
        pdt.assert_frame_equal(bv._coords.sort(axis=0),
                               subset_coords.sort(axis=0))
        pdt.assert_frame_equal(bv._metamap.sort(axis=0),
                               subset_metamap.sort(axis=0))

        # Takes a subset from coords
        bv = BaseVectors(self.coords, self.prop_expl, subset_metamap)
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

        bv = BaseVectors(subset_coords, self.prop_expl, subset_metamap)
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
            BaseVectors(self.coords, self.prop_expl, error_metamap)

    def test_make_groups(self):
        """Correctly generates the groups for vector_categories"""
        # Test with all categories
        bv = BaseVectors(self.coords, self.prop_expl, self.metamap)
        exp_groups = {'Treatment': {'Control': ['PC.354', 'PC.355', 'PC.356',
                                                'PC.481', 'PC.593'],
                                    'Fast': ['PC.607', 'PC.634',
                                             'PC.635', 'PC.636']},
                      'DOB': {'20061218': ['PC.354', 'PC.355'],
                              '20061126': ['PC.356'],
                              '20070314': ['PC.481'],
                              '20071210': ['PC.593'],
                              '20071112': ['PC.607'],
                              '20080116': ['PC.634', 'PC.635', 'PC.636']},
                      'Weight': {'60': ['PC.354'],
                                 '55': ['PC.355'],
                                 '50': ['PC.356'],
                                 '52': ['PC.481'],
                                 '57': ['PC.593'],
                                 '65': ['PC.607'],
                                 '68': ['PC.634'],
                                 '70': ['PC.635'],
                                 '72': ['PC.636']},
                      'Description': {'Control_mouse_I.D._354': ['PC.354'],
                                      'Control_mouse_I.D._355': ['PC.355'],
                                      'Control_mouse_I.D._356': ['PC.356'],
                                      'Control_mouse_I.D._481': ['PC.481'],
                                      'Control_mouse_I.D._593': ['PC.593'],
                                      'Fasting_mouse_I.D._607': ['PC.607'],
                                      'Fasting_mouse_I.D._634': ['PC.634'],
                                      'Fasting_mouse_I.D._635': ['PC.635'],
                                      'Fasting_mouse_I.D._636': ['PC.636']}}
        self.assertEqual(bv._groups, exp_groups)

        # Test with user-defined categories
        bv = BaseVectors(self.coords, self.prop_expl, self.metamap,
                         vector_categories=['Treatment', 'DOB'])
        exp_groups = {'Treatment': {'Control': ['PC.354', 'PC.355', 'PC.356',
                                                'PC.481', 'PC.593'],
                                    'Fast': ['PC.607', 'PC.634',
                                             'PC.635', 'PC.636']},
                      'DOB': {'20061218': ['PC.354', 'PC.355'],
                              '20061126': ['PC.356'],
                              '20070314': ['PC.481'],
                              '20071210': ['PC.593'],
                              '20071112': ['PC.607'],
                              '20080116': ['PC.634', 'PC.635', 'PC.636']}}
        self.assertEqual(bv._groups, exp_groups)

    def test_get_vectors(self):
        """Should raise a NotImplementedError as this is a base class"""
        bv = BaseVectors(self.coords, self.prop_expl, self.metamap)
        with self.assertRaises(NotImplementedError):
            bv.get_vectors()

    def test_get_group_vectors(self):
        """Should raise a NotImplementedError in usual execution as this is
        a base class"""
        bv = BaseVectors(self.coords, self.prop_expl, self.metamap)
        with self.assertRaises(NotImplementedError):
            bv.get_vectors()

    def test_get_group_vectors_error(self):
        """Should raise a RuntimeError if the user call _get_group_vectors
        with erroneous inputs"""
        bv = BaseVectors(self.coords, self.prop_expl, self.metamap)
        with self.assertRaises(RuntimeError):
            bv._get_group_vectors("foo", ['foo'])
        with self.assertRaises(RuntimeError):
            bv._get_group_vectors("bar", [])

    def test_compute_vector_results(self):
        """Should raise a NotImplementedError as this is a base class"""
        bv = BaseVectors(self.coords, self.prop_expl, self.metamap)
        with self.assertRaises(NotImplementedError):
            bv._compute_vector_results("foo", [])

    def test_weight_by_vector(self):
        """Correctly weights the vectors"""
        bv = BaseVectors(self.coords, self.prop_expl, self.metamap)

        vector = np.array([1, 2, 3, 4, 5, 6, 7, 8])
        w_vector = np.array([1, 5, 8, 12, 45, 80, 85, 90])
        exp_vector = np.array([1, 6.3571428571, 12.7142857142,
                               12.7142857142, 1.9264069264,
                               2.1795918367, 17.8, 20.3428571428])
        obs = bv._weight_by_vector(vector, w_vector)
        npt.assert_almost_equal(obs, exp_vector)

        vector = np.array([1, 2, 3, 4, 5, 6, 7, 8])
        w_vector = np.array([1, 2, 3, 4, 5, 6, 7, 8])
        exp_vector = np.array([1, 2, 3, 4, 5, 6, 7, 8])
        obs = bv._weight_by_vector(vector, w_vector)
        npt.assert_almost_equal(obs, exp_vector)

        vector = np.array([2, 3, 4, 5, 6])
        w_vector = np.array([25, 30, 35, 40, 45])
        exp_vector = np.array([2, 3, 4, 5, 6])
        obs = bv._weight_by_vector(vector, w_vector)
        npt.assert_almost_equal(obs, exp_vector)

        vector = np.array([[1, 2, 3], [2, 3, 4], [5, 6, 7], [8, 9, 10]])
        w_vector = np.array([1, 2, 3, 4])
        exp_vector = np.array([[1, 2, 3], [2, 3, 4], [5, 6, 7], [8, 9, 10]])
        obs = bv._weight_by_vector(vector, w_vector)
        npt.assert_almost_equal(obs, exp_vector)

    def test_weight_by_vector_error(self):
        """Raises an error with erroneous inputs"""
        bv = BaseVectors(self.coords, self.prop_expl, self.metamap)
        # Different vector lengths
        with self.assertRaises(ValueError):
            bv._weight_by_vector([1, 2, 3, 4], [1, 2, 3])

        # Inputs are not iterables
        with self.assertRaises(TypeError):
            bv._weight_by_vector(9, 1)

        # Weighting vector is not a gradient
        with self.assertRaises(ValueError):
            bv._weight_by_vector([1, 2, 3, 4], [1, 2, 3, 3])


class AverageVectorsTests(BaseTests):
    def test_get_vectors_all(self):
        """get_vectors returns the results of all categories"""
        av = AverageVectors(self.coords, self.prop_expl, self.metamap)
        obs = av.get_vectors()

        exp_description = CategoryResults('Description', None, None,
                                          "This value can not be used.")
        exp_weight = CategoryResults('Weight', None, None,
                                     "This value can not be used.")
        exp_20070314_group = GroupResults('20070314',
                                          np.array([2.1685208937828686]),
                                          2.16852089378,
                                          {'avg': 2.1685208937828686}, None)
        exp_20071112_group = GroupResults('20071112',
                                          np.array([7.3787312682853168]),
                                          7.37873126829,
                                          {'avg': [7.3787312682853168]}, None)
        exp_20080116_group = GroupResults('20080116',
                                          np.array([2.3430994255305362,
                                                    2.3946056125630912,
                                                    2.666555024970267]),
                                          2.46808668769,
                                          {'avg': 2.4680866876879648}, None)
        exp_20061126_group = GroupResults('20061126',
                                          np.array([7.6577228425502213]),
                                          7.65772284255,
                                          {'avg': [7.6577228425502213]}, None)
        exp_20061218_group = GroupResults('20061218',
                                          np.array([2.1607848468202744,
                                                    2.1607848468202744]),
                                          2.16078484682,
                                          {'avg': 2.1607848468202744}, None)
        exp_20071210_group = GroupResults('20071210',
                                          np.array([6.9553100288582872]),
                                          6.95531002886,
                                          {'avg': [6.9553100288582872]}, None)
        exp_dob = CategoryResults('DOB', 0.000139, [exp_20070314_group,
                                                    exp_20071112_group,
                                                    exp_20080116_group,
                                                    exp_20061126_group,
                                                    exp_20061218_group,
                                                    exp_20071210_group], None)
        exp_control_group = GroupResults('Control',
                                         np.array([2.3694943596755276,
                                                   3.3716388181385781,
                                                   5.4452089176253367,
                                                   4.5704258453173559,
                                                   4.4972603724478377]),
                                         4.05080566264,
                                         {'avg': 4.0508056626409275}, None)
        exp_fast_group = GroupResults('Fast', np.array([7.2220488239279126,
                                                        4.2726021564374372,
                                                        1.1169097274372082,
                                                        4.02717600030876]),
                                      4.15968417703,
                                      {'avg': 4.1596841770278292}, None)
        exp_treatment = CategoryResults('Treatment', 0.93311555,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = VectorsResults('avg', False, [exp_description, exp_weight,
                                            exp_dob, exp_treatment])
        self.assert_vectors_results_almost_equal(obs, exp)

    def test_get_vectors_single(self):
        """get_vectors returns the results of the provided category"""
        av = AverageVectors(self.coords, self.prop_expl, self.metamap,
                            vector_categories=['Treatment'])
        obs = av.get_vectors()

        exp_control_group = GroupResults('Control',
                                         np.array([2.3694943596755276,
                                                   3.3716388181385781,
                                                   5.4452089176253367,
                                                   4.5704258453173559,
                                                   4.4972603724478377]),
                                         4.05080566264,
                                         {'avg': 4.0508056626409275}, None)
        exp_fast_group = GroupResults('Fast', np.array([7.2220488239279126,
                                                        4.2726021564374372,
                                                        1.1169097274372082,
                                                        4.02717600030876]),
                                      4.15968417703,
                                      {'avg': 4.1596841770278292}, None)
        exp_treatment = CategoryResults('Treatment', 0.93311555,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = VectorsResults('avg', False, [exp_treatment])

        self.assert_vectors_results_almost_equal(obs, exp)

    def test_get_vectors_weighted(self):
        """get_vectors returns the correct weighted results"""
        av = AverageVectors(self.coords, self.prop_expl, self.metamap,
                            vector_categories=['Treatment'],
                            sort_category='Weight', weighted=True)
        obs = av.get_vectors()

        exp_control_group = GroupResults('Control',
                                         np.array([3.8296365043700837,
                                                   1.8849504232819798,
                                                   3.133650469739389,
                                                   3.0818770261785962,
                                                   1.9743045285131615]),
                                         2.78088379042,
                                         {'avg': 2.7808837904166421}, None)
        exp_fast_group = GroupResults('Fast', np.array([7.2187223309514401,
                                                        2.5522161259650256,
                                                        2.2349795833225015,
                                                        4.527821517691037]),
                                      4.13343488948,
                                      {'avg': 4.1334348894825013}, None)
        exp_treatment = CategoryResults('Treatment', 0.255388,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = VectorsResults('avg', True, [exp_treatment])

        self.assert_vectors_results_almost_equal(obs, exp)

# class TrajectoryVectorsTests(BaseTests):
#     def test_results(self):
#         tv = TrajectoryVectors(self.coord_dict, self.prop_expl, self.ids)
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
#         dv = DifferenceVectors(self.coord_dict, self.prop_expl, self.ids)
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
#         wdv = WindowDifferenceVectors(self.coord_dict, self.prop_expl,
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
