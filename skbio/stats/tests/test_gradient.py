# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import io
from operator import attrgetter
from unittest import TestCase, main

import numpy as np
import pandas as pd
import numpy.testing as npt
import pandas.testing as pdt

from skbio.util import get_data_path, assert_data_frame_almost_equal
from skbio.stats.gradient import (GradientANOVA, AverageGradientANOVA,
                                  TrajectoryGradientANOVA,
                                  FirstDifferenceGradientANOVA,
                                  WindowDifferenceGradientANOVA, GroupResults,
                                  CategoryResults, GradientANOVAResults,
                                  _weight_by_vector, _ANOVA_trajectories)


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

        metadata_map = {'PC.354': {'Treatment': 'Control',
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
        self.metadata_map = pd.DataFrame.from_dict(metadata_map,
                                                   orient='index')

        self.prop_expl = np.array([25.6216900347, 15.7715955926,
                                   14.1215046787, 11.6913885817, 9.83044890697,
                                   8.51253468595, 7.88775505332, 6.56308246609,
                                   4.42499350906e-16])

        gr_wo_msg = GroupResults('Foo', np.array([-2.6750, -0.2510,
                                                  -2.8322, 0.]),
                                 -1.4398, {'mean': -1.4398, 'std': 1.3184},
                                 None)
        gr_w_msg = GroupResults('Bar', np.array([9.6823, 2.9511, 5.2434]),
                                5.9589, {'mean': 5.9589, 'std': 2.7942},
                                "Cannot calculate the first difference "
                                "with a window of size (3).")
        self.groups = [gr_wo_msg, gr_w_msg]

        cr_no_data = CategoryResults('foo', None, None,
                                     'This group can not be used. All groups '
                                     'should have more than 1 element.')
        cr_data = CategoryResults('bar', 0.0110, self.groups, None)
        self.categories = [cr_no_data, cr_data]

        vr = GradientANOVAResults('wdiff', True, self.categories)

        description = CategoryResults('Description', None, None,
                                      'This group can not be used. All groups '
                                      'should have more than 1 element.')
        weight = CategoryResults('Weight', None, None,
                                 'This group can not be used. All groups '
                                 'should have more than 1 element.')
        dob = CategoryResults('DOB', None, None,
                              'This group can not be used. All groups '
                              'should have more than 1 element.')
        control_group = GroupResults('Control', np.array([2.3694, 3.3716,
                                                          5.4452, 4.5704,
                                                          4.4972]),
                                     4.0508, {'avg': 4.0508}, None)
        fast_group = GroupResults('Fast', np.array([7.2220, 4.2726, 1.1169,
                                                    4.0271]),
                                  4.1596, {'avg': 4.1596}, None)
        treatment = CategoryResults('Treatment', 0.9331,
                                    [control_group, fast_group], None)
        vr_real = GradientANOVAResults('avg', False, [description, weight, dob,
                                                      treatment])

        self.vec_results = [vr, vr_real]

    # This function makes the comparisons between the results classes easier
    def assert_group_results_almost_equal(self, obs, exp):
        """Tests that obs and exp are almost equal"""
        self.assertEqual(obs.name, exp.name)
        npt.assert_almost_equal(obs.trajectory, exp.trajectory)
        npt.assert_almost_equal(obs.mean, exp.mean)
        self.assertEqual(obs.info.keys(), exp.info.keys())
        for key in obs.info:
            npt.assert_almost_equal(obs.info[key], exp.info[key])
        self.assertEqual(obs.message, exp.message)

    def assert_category_results_almost_equal(self, obs, exp):
        """Tests that obs and exp are almost equal"""
        self.assertEqual(obs.category, exp.category)

        if exp.probability is None:
            self.assertTrue(obs.probability is None)
            self.assertTrue(obs.groups is None)
        else:
            npt.assert_almost_equal(obs.probability, exp.probability)
            for o, e in zip(sorted(obs.groups, key=attrgetter('name')),
                            sorted(exp.groups, key=attrgetter('name'))):
                self.assert_group_results_almost_equal(o, e)

    def assert_gradientANOVA_results_almost_equal(self, obs, exp):
        """Tests that obs and exp are almost equal"""
        self.assertEqual(obs.algorithm, exp.algorithm)
        self.assertEqual(obs.weighted, exp.weighted)

        for o, e in zip(sorted(obs.categories, key=attrgetter('category')),
                        sorted(exp.categories, key=attrgetter('category'))):
            self.assert_category_results_almost_equal(o, e)


class GradientTests(BaseTests):
    def test_weight_by_vector(self):
        """Correctly weights the vectors"""
        trajectory = pd.DataFrame.from_dict({'s1': np.array([1]),
                                             's2': np.array([2]),
                                             's3': np.array([3]),
                                             's4': np.array([4]),
                                             's5': np.array([5]),
                                             's6': np.array([6]),
                                             's7': np.array([7]),
                                             's8': np.array([8])},
                                            orient='index')
        trajectory.sort_values(by=0, inplace=True)
        w_vector = pd.Series(np.array([1, 5, 8, 12, 45, 80, 85, 90]),
                             ['s1', 's2', 's3', 's4',
                              's5', 's6', 's7', 's8']).astype(np.float64)
        exp = pd.DataFrame.from_dict({'s1': np.array([1]),
                                      's2': np.array([6.3571428571]),
                                      's3': np.array([12.7142857142]),
                                      's4': np.array([12.7142857142]),
                                      's5': np.array([1.9264069264]),
                                      's6': np.array([2.1795918367]),
                                      's7': np.array([17.8]),
                                      's8': np.array([20.3428571428])},
                                     orient='index').astype(np.float64)
        obs = _weight_by_vector(trajectory, w_vector)
        assert_data_frame_almost_equal(obs.sort_index(), exp.sort_index())

        trajectory = pd.DataFrame.from_dict({'s1': np.array([1]),
                                             's2': np.array([2]),
                                             's3': np.array([3]),
                                             's4': np.array([4]),
                                             's5': np.array([5]),
                                             's6': np.array([6]),
                                             's7': np.array([7]),
                                             's8': np.array([8])},
                                            orient='index')
        trajectory.sort_values(by=0, inplace=True)
        w_vector = pd.Series(np.array([1, 2, 3, 4, 5, 6, 7, 8]),
                             ['s1', 's2', 's3', 's4',
                              's5', 's6', 's7', 's8']).astype(np.float64)
        exp = pd.DataFrame.from_dict({'s1': np.array([1]),
                                      's2': np.array([2]),
                                      's3': np.array([3]),
                                      's4': np.array([4]),
                                      's5': np.array([5]),
                                      's6': np.array([6]),
                                      's7': np.array([7]),
                                      's8': np.array([8])
                                      },
                                     orient='index').astype(np.float64)
        obs = _weight_by_vector(trajectory, w_vector)
        assert_data_frame_almost_equal(obs.sort_index(), exp.sort_index())

        trajectory = pd.DataFrame.from_dict({'s2': np.array([2]),
                                             's3': np.array([3]),
                                             's4': np.array([4]),
                                             's5': np.array([5]),
                                             's6': np.array([6])},
                                            orient='index')
        trajectory.sort_values(by=0, inplace=True)
        w_vector = pd.Series(np.array([25, 30, 35, 40, 45]),
                             ['s2', 's3', 's4', 's5', 's6']).astype(np.float64)
        exp = pd.DataFrame.from_dict({'s2': np.array([2]),
                                      's3': np.array([3]),
                                      's4': np.array([4]),
                                      's5': np.array([5]),
                                      's6': np.array([6])},
                                     orient='index').astype(np.float64)
        obs = _weight_by_vector(trajectory, w_vector)
        assert_data_frame_almost_equal(obs.sort_index(), exp.sort_index())

        trajectory = pd.DataFrame.from_dict({'s1': np.array([1, 2, 3]),
                                             's2': np.array([2, 3, 4]),
                                             's3': np.array([5, 6, 7]),
                                             's4': np.array([8, 9, 10])},
                                            orient='index')
        trajectory.sort_values(by=0, inplace=True)
        w_vector = pd.Series(np.array([1, 2, 3, 4]),
                             ['s1', 's2', 's3', 's4']).astype(np.float64)
        exp = pd.DataFrame.from_dict({'s1': np.array([1, 2, 3]),
                                      's2': np.array([2, 3, 4]),
                                      's3': np.array([5, 6, 7]),
                                      's4': np.array([8, 9, 10])},
                                     orient='index').astype(np.float64)
        obs = _weight_by_vector(trajectory, w_vector)
        assert_data_frame_almost_equal(obs.sort_index(), exp.sort_index())

        sample_ids = ['PC.356', 'PC.481', 'PC.355', 'PC.593', 'PC.354']
        trajectory = pd.DataFrame.from_dict({'PC.356': np.array([5.65948525,
                                                                 1.37977545,
                                                                 -4.9706303]),
                                             'PC.481': np.array([0.79151484,
                                                                 -0.70387996,
                                                                 1.89223152]),
                                             'PC.355': np.array([6.05869624,
                                                                 3.44821245,
                                                                 -0.42595788]),
                                             'PC.593': np.array([5.18731945,
                                                                 -1.81714206,
                                                                 4.26216485]),
                                             'PC.354': np.array([7.07588529,
                                                                 -0.53917873,
                                                                 0.89389158])
                                             }, orient='index')
        w_vector = pd.Series(np.array([50, 52, 55, 57, 60]),
                             sample_ids).astype(np.float64)
        exp = pd.DataFrame.from_dict({'PC.356': np.array([5.65948525,
                                                          1.37977545,
                                                          -4.9706303]),
                                      'PC.481': np.array([0.98939355,
                                                          -0.87984995,
                                                          2.3652894]),
                                      'PC.355': np.array([5.04891353,
                                                          2.87351038,
                                                          -0.3549649]),
                                      'PC.593': np.array([6.48414931,
                                                          -2.27142757,
                                                          5.32770606]),
                                      'PC.354': np.array([5.89657108,
                                                          -0.44931561,
                                                          0.74490965])
                                      }, orient='index')
        obs = _weight_by_vector(trajectory.loc[sample_ids],
                                w_vector[sample_ids])
        assert_data_frame_almost_equal(obs.sort_index(), exp.sort_index())

    def test_weight_by_vector_single_element(self):
        trajectory = pd.DataFrame.from_dict({'s1': np.array([42])},
                                            orient='index')
        w_vector = pd.Series(np.array([5]), ['s1']).astype(np.float64)

        obs = _weight_by_vector(trajectory, w_vector)
        assert_data_frame_almost_equal(obs, trajectory)

    def test_weight_by_vector_error(self):
        """Raises an error with erroneous inputs"""
        # Different vector lengths
        with self.assertRaises(ValueError):
            _weight_by_vector([1, 2, 3, 4], [1, 2, 3])

        # Inputs are not iterables
        with self.assertRaises(TypeError):
            _weight_by_vector(9, 1)

        # Weighting vector is not a gradient
        with self.assertRaises(ValueError):
            _weight_by_vector([1, 2, 3, 4], [1, 2, 3, 3])

    def test_ANOVA_trajectories(self):
        """Correctly performs the check before running ANOVA"""
        # Only one group in a given category
        group = GroupResults('Bar', np.array([2.3694943596755276,
                                              3.3716388181385781,
                                              5.4452089176253367,
                                              4.5704258453173559,
                                              4.4972603724478377]),
                             4.05080566264, {'avg': 4.0508056626409275}, None)
        obs = _ANOVA_trajectories('Foo', [group])
        exp = CategoryResults('Foo', None, None,
                              'Only one value in the group.')
        self.assert_category_results_almost_equal(obs, exp)

        # One element have only one element
        group2 = GroupResults('FooBar', np.array([4.05080566264]),
                              4.05080566264, {'avg': 4.05080566264}, None)
        obs = _ANOVA_trajectories('Foo', [group, group2])
        exp = CategoryResults('Foo', None, None,
                              'This group can not be used. All groups '
                              'should have more than 1 element.')
        self.assert_category_results_almost_equal(obs, exp)

        gr1 = GroupResults('Foo', np.array([-0.219044992, 0.079674486,
                                            0.09233683]),
                           -0.015677892, {'avg': -0.015677892}, None)
        gr2 = GroupResults('Bar', np.array([-0.042258081, 0.000204041,
                                            0.024837603]),
                           -0.0732878716, {'avg': -0.0732878716}, None)
        gr3 = GroupResults('FBF', np.array([0.080504323, -0.212014503,
                                            -0.088353435]),
                           -0.0057388123, {'avg': -0.0057388123}, None)
        obs = _ANOVA_trajectories('Cat', [gr1, gr2, gr3])
        exp = CategoryResults('Cat', 0.8067456876, [gr1, gr2, gr3], None)
        self.assert_category_results_almost_equal(obs, exp)


class GroupResultsTests(BaseTests):
    def test_to_file(self):
        out_paths = ['gr_wo_msg_out', 'gr_w_msg_out']
        raw_paths = ['gr_wo_msg_raw', 'gr_w_msg_raw']

        for gr, out_fp, raw_fp in zip(self.groups, out_paths, raw_paths):
            obs_out_f = io.StringIO()
            obs_raw_f = io.StringIO()
            gr.to_files(obs_out_f, obs_raw_f)
            obs_out = obs_out_f.getvalue()
            obs_raw = obs_raw_f.getvalue()
            obs_out_f.close()
            obs_raw_f.close()

            with open(get_data_path(out_fp)) as f:
                exp_out = f.read()

            with open(get_data_path(raw_fp)) as f:
                exp_raw = f.read()

            self.assertEqual(obs_out, exp_out)
            self.assertEqual(obs_raw, exp_raw)


class CategoryResultsTests(BaseTests):
    def test_to_file(self):
        out_paths = ['cr_no_data_out', 'cr_data_out']
        raw_paths = ['cr_no_data_raw', 'cr_data_raw']

        for cat, out_fp, raw_fp in zip(self.categories, out_paths, raw_paths):
            obs_out_f = io.StringIO()
            obs_raw_f = io.StringIO()
            cat.to_files(obs_out_f, obs_raw_f)
            obs_out = obs_out_f.getvalue()
            obs_raw = obs_raw_f.getvalue()
            obs_out_f.close()
            obs_raw_f.close()

            with open(get_data_path(out_fp)) as f:
                exp_out = f.read()

            with open(get_data_path(raw_fp)) as f:
                exp_raw = f.read()

            self.assertEqual(obs_out, exp_out)
            self.assertEqual(obs_raw, exp_raw)


class GradientANOVAResultsTests(BaseTests):
    def test_to_file(self):
        out_paths = ['vr_out']
        raw_paths = ['vr_raw']

        for vr, out_fp, raw_fp in zip(self.vec_results, out_paths, raw_paths):
            obs_out_f = io.StringIO()
            obs_raw_f = io.StringIO()
            vr.to_files(obs_out_f, obs_raw_f)
            obs_out = obs_out_f.getvalue()
            obs_raw = obs_raw_f.getvalue()
            obs_out_f.close()
            obs_raw_f.close()

            with open(get_data_path(out_fp)) as f:
                exp_out = f.read()

            with open(get_data_path(raw_fp)) as f:
                exp_raw = f.read()

            self.assertEqual(obs_out, exp_out)
            self.assertEqual(obs_raw, exp_raw)


class GradientANOVATests(BaseTests):
    def test_init(self):
        """Correctly initializes the class attributes"""
        # Note self._groups is tested on test_make_groups
        # so we are not testing it here

        # Test with weighted = False
        bv = GradientANOVA(self.coords, self.prop_expl, self.metadata_map)

        assert_data_frame_almost_equal(bv._coords.loc[self.coords_3axes.index],
                                       self.coords_3axes)
        exp_prop_expl = np.array([25.6216900347, 15.7715955926,
                                  14.1215046787])
        npt.assert_equal(bv._prop_expl, exp_prop_expl)
        assert_data_frame_almost_equal(bv._metadata_map.loc[self.metadata_map.index],  # noqa
                                       self.metadata_map)
        self.assertTrue(bv._weighting_vector is None)
        self.assertFalse(bv._weighted)

        # Test with weighted = True
        bv = GradientANOVA(self.coords, self.prop_expl, self.metadata_map,
                           sort_category='Weight', weighted=True)

        assert_data_frame_almost_equal(bv._coords.loc[self.coords_3axes.index],
                                       self.coords_3axes)
        npt.assert_equal(bv._prop_expl, exp_prop_expl)
        assert_data_frame_almost_equal(bv._metadata_map.loc[self.metadata_map.index],  # noqa
                                       self.metadata_map)
        exp_weighting_vector = pd.Series(
            np.array([60, 55, 50, 52, 57, 65, 68, 70, 72]),
            ['PC.354', 'PC.355', 'PC.356', 'PC.481', 'PC.593', 'PC.607',
             'PC.634', 'PC.635', 'PC.636'], name='Weight'
            ).astype(np.float64)
        pdt.assert_series_equal(bv._weighting_vector.loc[exp_weighting_vector.index],  # noqa
                                exp_weighting_vector)
        self.assertTrue(bv._weighted)

    def test_init_error(self):
        """Raises an error with erroneous inputs"""
        # Raises ValueError if any category in trajectory_categories is not
        # present in metadata_map
        with self.assertRaises(ValueError):
            GradientANOVA(self.coords, self.prop_expl, self.metadata_map,
                          trajectory_categories=['foo'])
        with self.assertRaises(ValueError):
            GradientANOVA(self.coords, self.prop_expl, self.metadata_map,
                          trajectory_categories=['Weight', 'Treatment', 'foo'])

        # Raises ValueError if sort_category is not present in metadata_map
        with self.assertRaises(ValueError):
            GradientANOVA(self.coords, self.prop_expl, self.metadata_map,
                          sort_category='foo')

        # Raises ValueError if weighted == True and sort_category == None
        with self.assertRaises(ValueError):
            GradientANOVA(self.coords, self.prop_expl, self.metadata_map,
                          weighted=True)

        # Raises ValueError if weighted == True and the values under
        # sort_category are not numerical
        with self.assertRaises(ValueError):
            GradientANOVA(self.coords, self.prop_expl, self.metadata_map,
                          sort_category='Treatment', weighted=True)

        # Raises ValueError if axes > len(prop_expl)
        with self.assertRaises(ValueError):
            GradientANOVA(self.coords, self.prop_expl, self.metadata_map,
                          axes=10)

        # Raises ValueError if axes < 0
        with self.assertRaises(ValueError):
            GradientANOVA(self.coords, self.prop_expl, self.metadata_map,
                          axes=-1)

    def test_normalize_samples(self):
        """Correctly normalizes the samples between coords and metadata_map"""
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

        metadata_map = {'PC.355': {'Treatment': 'Control',
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
        subset_metadata_map = pd.DataFrame.from_dict(metadata_map,
                                                     orient='index')

        # Takes a subset from metadata_map
        bv = GradientANOVA(subset_coords, self.prop_expl, self.metadata_map)
        assert_data_frame_almost_equal(
            bv._coords.sort_index(),
            subset_coords.sort_index())
        assert_data_frame_almost_equal(
            bv._metadata_map.sort_index(),
            subset_metadata_map.sort_index())

        # Takes a subset from coords
        bv = GradientANOVA(self.coords, self.prop_expl, subset_metadata_map)
        assert_data_frame_almost_equal(
            bv._coords.sort_index(),
            subset_coords.sort_index())
        assert_data_frame_almost_equal(
            bv._metadata_map.sort_index(),
            subset_metadata_map.sort_index())

        # Takes a subset from metadata_map and coords at the same time
        coord_data = {
            'PC.636': np.array([-0.212230626531, 0.216034194368,
                                0.03532727349]),
            'PC.635': np.array([-0.277487312135, -0.0295483215975,
                                -0.0744173437992]),
            'PC.355': np.array([0.236467470907, 0.21863434374,
                                -0.0301637746424])
            }
        subset_coords = pd.DataFrame.from_dict(coord_data, orient='index')

        metadata_map = {'PC.355': {'Treatment': 'Control',
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
        subset_metadata_map = pd.DataFrame.from_dict(metadata_map,
                                                     orient='index')

        bv = GradientANOVA(subset_coords, self.prop_expl, subset_metadata_map)
        exp_coords = pd.DataFrame.from_dict(
            {'PC.355': np.array([0.236467470907, 0.21863434374,
                                 -0.0301637746424])},
            orient='index')
        assert_data_frame_almost_equal(
            bv._coords.sort_index(),
            exp_coords.sort_index())
        exp_metadata_map = pd.DataFrame.from_dict(
            {'PC.355': {'Treatment': 'Control',
                        'DOB': '20061218',
                        'Weight': '55',
                        'Description': 'Control_mouse_I.D._355'}},
            orient='index')
        assert_data_frame_almost_equal(
            bv._metadata_map.sort_index(),
            exp_metadata_map.sort_index())

    def test_normalize_samples_error(self):
        """Raises an error if coords and metadata_map does not have samples in
        common"""
        error_metadata_map = pd.DataFrame.from_dict(
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
            GradientANOVA(self.coords, self.prop_expl, error_metadata_map)

    def test_make_groups(self):
        """Correctly generates the groups for trajectory_categories"""
        # Test with all categories
        bv = GradientANOVA(self.coords, self.prop_expl, self.metadata_map)
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
        bv = GradientANOVA(self.coords, self.prop_expl, self.metadata_map,
                           trajectory_categories=['Treatment', 'DOB'])
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

    def test_make_groups_natural_sorting(self):
        # Ensure sample IDs are sorted using a natural sorting algorithm.
        df = pd.DataFrame.from_dict({
            'a2': {'Col1': 'foo', 'Col2': '1.0'},
            'a1': {'Col1': 'bar', 'Col2': '-42.0'},
            'a11.0': {'Col1': 'foo', 'Col2': '2e-5'},
            'a-10': {'Col1': 'foo', 'Col2': '5'},
            'a10': {'Col1': 'bar', 'Col2': '5'}},
            orient='index')

        coords = pd.DataFrame.from_dict({
            'a10': np.array([-0.212230626531, 0.216034194368, 0.03532727349]),
            'a11.0': np.array([-0.277487312135, -0.0295483215975,
                               -0.0744173437992]),
            'a1': np.array([0.220886492631, 0.0874848360559,
                            -0.351990132198]),
            'a2': np.array([0.0308923744062, -0.0446295973489,
                            0.133996451689]),
            'a-10': np.array([0.27616778138, -0.0341866951102,
                              0.0633000238256])},
            orient='index')

        prop_expl = np.array([25.6216900347, 15.7715955926, 14.1215046787,
                              11.6913885817, 9.83044890697])

        # Sort by sample IDs.
        ga = GradientANOVA(coords, prop_expl, df)

        exp_groups = {
            'Col1': {
                'foo': ['a-10', 'a2', 'a11.0'],
                'bar': ['a1', 'a10']
            },
            'Col2': {
                '1.0': ['a2'],
                '-42.0': ['a1'],
                '2e-5': ['a11.0'],
                '5': ['a-10', 'a10']
            }
        }

        self.assertEqual(ga._groups, exp_groups)

        # Sort sample IDs by Col2.
        ga = GradientANOVA(coords, prop_expl, df,
                           trajectory_categories=['Col1'],
                           sort_category='Col2')

        exp_groups = {
            'Col1': {
                'foo': ['a11.0', 'a2', 'a-10'],
                'bar': ['a1', 'a10']
            }
        }

        self.assertEqual(ga._groups, exp_groups)

    def test_get_trajectories(self):
        """Should raise a NotImplementedError as this is a base class"""
        bv = GradientANOVA(self.coords, self.prop_expl, self.metadata_map)
        with self.assertRaises(NotImplementedError):
            bv.get_trajectories()

    def test_get_group_trajectories(self):
        """Should raise a NotImplementedError in usual execution as this is
        a base class"""
        bv = GradientANOVA(self.coords, self.prop_expl, self.metadata_map)
        with self.assertRaises(NotImplementedError):
            bv.get_trajectories()

    def test_get_group_trajectories_error(self):
        """Should raise a RuntimeError if the user call _get_group_trajectories
        with erroneous inputs"""
        bv = GradientANOVA(self.coords, self.prop_expl, self.metadata_map)
        with self.assertRaises(KeyError):
            bv._get_group_trajectories("foo", ['foo'])
        with self.assertRaises(RuntimeError):
            bv._get_group_trajectories("bar", [])

    def test_compute_trajectories_results(self):
        """Should raise a NotImplementedError as this is a base class"""
        bv = GradientANOVA(self.coords, self.prop_expl, self.metadata_map)
        with self.assertRaises(NotImplementedError):
            bv._compute_trajectories_results("foo", [])


class AverageGradientANOVATests(BaseTests):
    def test_get_trajectories_all(self):
        """get_trajectories returns the results of all categories"""
        av = AverageGradientANOVA(self.coords, self.prop_expl,
                                  self.metadata_map)
        obs = av.get_trajectories()

        exp_description = CategoryResults('Description', None, None,
                                          'This group can not be used. All '
                                          'groups should have more than 1 '
                                          'element.')
        exp_weight = CategoryResults('Weight', None, None,
                                     'This group can not be used. All groups '
                                     'should have more than 1 element.')
        exp_dob = CategoryResults('DOB', None, None,
                                  'This group can not be used. All groups '
                                  'should have more than 1 element.')
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
        exp = GradientANOVAResults('avg', False, [exp_description, exp_weight,
                                                  exp_dob, exp_treatment])
        self.assert_gradientANOVA_results_almost_equal(obs, exp)

    def test_get_trajectories_single(self):
        """get_trajectories returns the results of the provided category"""
        av = AverageGradientANOVA(self.coords, self.prop_expl,
                                  self.metadata_map,
                                  trajectory_categories=['Treatment'])
        obs = av.get_trajectories()

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
        exp = GradientANOVAResults('avg', False, [exp_treatment])

        self.assert_gradientANOVA_results_almost_equal(obs, exp)

    def test_get_trajectories_weighted(self):
        """get_trajectories returns the correct weighted results"""
        av = AverageGradientANOVA(self.coords, self.prop_expl,
                                  self.metadata_map,
                                  trajectory_categories=['Treatment'],
                                  sort_category='Weight', weighted=True)
        obs = av.get_trajectories()
        exp_control_group = GroupResults('Control', np.array([5.7926887872,
                                                              4.3242308936,
                                                              2.9212403501,
                                                              5.5400792151,
                                                              1.2326804315]),
                                         3.9621839355,
                                         {'avg': 3.9621839355}, None)
        exp_fast_group = GroupResults('Fast', np.array([7.2187223286,
                                                        2.5522161282,
                                                        2.2349795861,
                                                        4.5278215248]),
                                      4.1334348919,
                                      {'avg': 4.1334348919}, None)
        exp_treatment = CategoryResults('Treatment', 0.9057666800,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = GradientANOVAResults('avg', True, [exp_treatment])
        self.assert_gradientANOVA_results_almost_equal(obs, exp)


class TrajectoryGradientANOVATests(BaseTests):

    def test_get_trajectories(self):
        tv = TrajectoryGradientANOVA(self.coords, self.prop_expl,
                                     self.metadata_map,
                                     trajectory_categories=['Treatment'],
                                     sort_category='Weight')
        obs = tv.get_trajectories()
        exp_control_group = GroupResults('Control', np.array([8.6681963576,
                                                              7.0962717982,
                                                              7.1036434615,
                                                              4.0675712674]),
                                         6.73392072123,
                                         {'2-norm': 13.874494152}, None)
        exp_fast_group = GroupResults('Fast', np.array([11.2291654905,
                                                        3.9163741156,
                                                        4.4943507388]),
                                      6.5466301150,
                                      {'2-norm': 12.713431181}, None)
        exp_treatment = CategoryResults('Treatment', 0.9374500147,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = GradientANOVAResults('trajectory', False, [exp_treatment])
        self.assert_gradientANOVA_results_almost_equal(obs, exp)

    def test_get_trajectories_weighted(self):
        tv = TrajectoryGradientANOVA(self.coords, self.prop_expl,
                                     self.metadata_map,
                                     trajectory_categories=['Treatment'],
                                     sort_category='Weight', weighted=True)
        obs = tv.get_trajectories()
        exp_control_group = GroupResults('Control', np.array([8.9850643421,
                                                              6.1617529749,
                                                              7.7989125908,
                                                              4.9666249268]),
                                         6.9780887086,
                                         {'2-norm': 14.2894710091}, None)
        exp_fast_group = GroupResults('Fast', np.array([9.6823682852,
                                                        2.9511115209,
                                                        5.2434091953]),
                                      5.9589630005,
                                      {'2-norm': 11.3995901159}, None)
        exp_treatment = CategoryResults('Treatment', 0.6248157720,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = GradientANOVAResults('trajectory', True, [exp_treatment])
        self.assert_gradientANOVA_results_almost_equal(obs, exp)


class FirstDifferenceGradientANOVATests(BaseTests):
    def test_get_trajectories(self):
        dv = FirstDifferenceGradientANOVA(self.coords, self.prop_expl,
                                          self.metadata_map,
                                          trajectory_categories=['Treatment'],
                                          sort_category='Weight')
        obs = dv.get_trajectories()
        exp_control_group = GroupResults('Control', np.array([-1.5719245594,
                                                              0.0073716633,
                                                              -3.0360721941]),
                                         -1.5335416967,
                                         {'mean': -1.5335416967,
                                          'std': 1.2427771485}, None)
        exp_fast_group = GroupResults('Fast', np.array([-7.3127913749,
                                                        0.5779766231]),
                                      -3.3674073758,
                                      {'mean': -3.3674073758,
                                       'std': 3.9453839990}, None)
        exp_treatment = CategoryResults('Treatment', 0.6015260608,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = GradientANOVAResults('diff', False, [exp_treatment])
        self.assert_gradientANOVA_results_almost_equal(obs, exp)

    def test_get_trajectories_weighted(self):
        dv = FirstDifferenceGradientANOVA(self.coords, self.prop_expl,
                                          self.metadata_map,
                                          trajectory_categories=['Treatment'],
                                          sort_category='Weight',
                                          weighted=True)
        obs = dv.get_trajectories()
        exp_control_group = GroupResults('Control', np.array([-2.8233113671,
                                                              1.6371596158,
                                                              -2.8322876639]),
                                         -1.3394798050,
                                         {'mean': -1.3394798050,
                                          'std': 2.1048051097}, None)
        exp_fast_group = GroupResults('Fast', np.array([-6.7312567642,
                                                        2.2922976743]),
                                      -2.2194795449,
                                      {'mean': -2.2194795449,
                                       'std': 4.5117772193}, None)
        exp_treatment = CategoryResults('Treatment', 0.8348644420,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = GradientANOVAResults('diff', True, [exp_treatment])
        self.assert_gradientANOVA_results_almost_equal(obs, exp)


class WindowDifferenceGradientANOVATests(BaseTests):
    def test_get_trajectories(self):
        wdv = WindowDifferenceGradientANOVA(
            self.coords, self.prop_expl, self.metadata_map, 3,
            trajectory_categories=['Treatment'], sort_category='Weight')
        obs = wdv.get_trajectories()
        exp_control_group = GroupResults('Control', np.array([-2.5790341819,
                                                              -2.0166764661,
                                                              -3.0360721941,
                                                              0.]),
                                         -1.9079457105,
                                         {'mean': -1.9079457105,
                                          'std': 1.1592139913}, None)
        exp_fast_group = GroupResults('Fast', np.array([11.2291654905,
                                                        3.9163741156,
                                                        4.4943507388]),
                                      6.5466301150,
                                      {'mean': 6.5466301150,
                                       'std': 3.3194494926},
                                      "Cannot calculate the first difference "
                                      "with a window of size (3).")
        exp_treatment = CategoryResults('Treatment', 0.0103976830,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = GradientANOVAResults('wdiff', False, [exp_treatment])
        self.assert_gradientANOVA_results_almost_equal(obs, exp)

    def test_get_trajectories_weighted(self):
        wdv = WindowDifferenceGradientANOVA(
            self.coords, self.prop_expl, self.metadata_map, 3,
            trajectory_categories=['Treatment'], sort_category='Weight',
            weighted=True)
        obs = wdv.get_trajectories()
        exp_control_group = GroupResults('Control', np.array([-2.6759675112,
                                                              -0.2510321601,
                                                              -2.8322876639,
                                                              0.]),
                                         -1.4398218338,
                                         {'mean': -1.4398218338,
                                          'std': 1.31845790844}, None)
        exp_fast_group = GroupResults('Fast', np.array([9.6823682852,
                                                        2.9511115209,
                                                        5.2434091953]),
                                      5.9589630005,
                                      {'mean': 5.9589630005,
                                       'std': 2.7942163293},
                                      "Cannot calculate the first difference "
                                      "with a window of size (3).")
        exp_treatment = CategoryResults('Treatment', 0.0110675605,
                                        [exp_control_group, exp_fast_group],
                                        None)
        exp = GradientANOVAResults('wdiff', True, [exp_treatment])
        self.assert_gradientANOVA_results_almost_equal(obs, exp)


if __name__ == '__main__':
    main()
