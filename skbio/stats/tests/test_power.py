
# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------


from __future__ import absolute_import, division, print_function
from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd
from scipy.stats import kruskal

from skbio.stats.power import (subsample_power,
                               subsample_paired_power,
                               _check_strs,
                               confidence_bound,
                               _calculate_power,
                               _compare_distributions,
                               _calculate_power_curve,
                               bootstrap_power_curve,
                               paired_subsamples)


class PowerAnalysisTest(TestCase):

    def setUp(self):
        # Defines a testing function
        def test_meta(ids, meta, cat, div):
            """Checks thhe div metric with a kruskal wallis"""
            out = [meta.loc[id_, div] for id_ in ids]
            return kruskal(*out)[1]
        self.test_meta = test_meta
        # Sets the random seed
        np.random.seed(5)
        # Sets up the distributions of data for use
        self.s1 = np.arange(0, 10, 1)
        # Sets up two distributions which will never be equal by a rank-sum
        # test.
        self.samps = [np.ones((10))/10., np.ones((10))]
        self.pop = [np.arange(0, 10, 0.1), np.arange(0, 20, 0.2)]
        # Sets up a vector of alpha values
        self.alpha = np.power(10, np.array([-1, -1.301, -2, -3])).round(3)
        # Sets up a vector of samples
        self.num_samps = np.arange(10, 100, 10)
        # Sets up the test function, a rank-sum test
        self.f = lambda x: kruskal(*x)[1]
        # Sets up a mapping file
        meta = {'GW': {'INT': 'N', 'ABX': np.nan, 'DIV': 19.5, 'AGE': '30s',
                       'SEX': 'M'},
                'CB': {'INT': 'Y', 'ABX': np.nan, 'DIV': 42.7, 'AGE': '30s',
                       'SEX': 'M'},
                'WM': {'INT': 'N', 'ABX': 'N', 'DIV': 27.5, 'AGE': '20s',
                       'SEX': 'F'},
                'MH': {'INT': 'Y', 'ABX': 'N', 'DIV': 62.3, 'AGE': '30s',
                       'SEX': 'F'},
                'CD': {'INT': 'Y', 'ABX': 'Y', 'DIV': 36.4, 'AGE': '40s',
                       'SEX': 'F'},
                'LF': {'INT': 'Y', 'ABX': 'N', 'DIV': 50.2, 'AGE': '20s',
                       'SEX': 'M'},
                'PP': {'INT': 'N', 'ABX': 'Y', 'DIV': 10.8, 'AGE': '30s',
                       'SEX': 'F'},
                'MM': {'INT': 'N', 'ABX': 'N', 'DIV': 55.6, 'AGE': '40s',
                       'SEX': 'F'},
                'SR': {'INT': 'N', 'ABX': 'Y', 'DIV': 2.2, 'AGE': '20s',
                       'SEX': 'M'},
                'TS': {'INT': 'N', 'ABX': 'Y', 'DIV': 16.1, 'AGE': '40s',
                       'SEX': 'M'},
                'PC': {'INT': 'Y', 'ABX': 'N', 'DIV': 82.6, 'AGE': '40s',
                       'SEX': 'M'},
                'NR': {'INT': 'Y', 'ABX': 'Y', 'DIV': 15.7, 'AGE': '20s',
                       'SEX': 'F'}}
        self.meta = pd.DataFrame.from_dict(meta, orient='index')
        self.meta_f = lambda x: test_meta(x, self.meta, 'INT', 'DIV')
        self.counts = np.array([5, 15, 25, 35, 45])
        self.powers = [np.array([[0.105, 0.137, 0.174, 0.208, 0.280],
                                 [0.115, 0.135, 0.196, 0.204, 0.281],
                                 [0.096, 0.170, 0.165, 0.232, 0.256],
                                 [0.122, 0.157, 0.202, 0.250, 0.279],
                                 [0.132, 0.135, 0.173, 0.203, 0.279]]),
                       np.array([[0.157, 0.345, 0.522, 0.639, 0.739],
                                 [0.159, 0.374, 0.519, 0.646, 0.757],
                                 [0.161, 0.339, 0.532, 0.634, 0.745],
                                 [0.169, 0.372, 0.541, 0.646, 0.762],
                                 [0.163, 0.371, 0.522, 0.648, 0.746]]),
                       np.array([[0.276, 0.626, 0.865, 0.927, 0.992],
                                 [0.267, 0.667, 0.848, 0.937, 0.978],
                                 [0.236, 0.642, 0.850, 0.935, 0.977],
                                 [0.249, 0.633, 0.828, 0.955, 0.986],
                                 [0.249, 0.663, 0.869, 0.951, 0.985]])]
        self.power_alpha = 0.1
        self.effects = np.array([0.15245, 0.34877, 0.55830])
        self.bounds = np.array([0.01049, 0.00299, 0.007492])
        self.labels = np.array(['Age', 'Intervenption', 'Antibiotics'])
        self.cats = np.array(['AGE', 'INT', 'ABX'])
        self.cat = "AGE"
        self.control_cats = ['INT', 'ABX']

    def test_subsample_power_matched_relationship_error(self):
        with self.assertRaises(ValueError):
            subsample_power(self.f,
                            samples=[np.ones((2)), np.ones((5))],
                            draw_mode="matched")

    def test_subsample_power_min_observations_error(self):
        with self.assertRaises(ValueError):
            subsample_power(self.f,
                            samples=[np.ones((2)), np.ones((5))])

    def test_subsample_power_interval_error(self):
        with self.assertRaises(ValueError):
            subsample_power(self.f,
                            samples=[np.ones((2)), np.ones((5))],
                            min_observations=2,
                            min_counts=5,
                            max_counts=7)

    def test_subsample_power(self):
        test_p, test_c = subsample_power(self.f,
                                         samples=self.pop,
                                         num_iter=10,
                                         num_runs=2,
                                         min_counts=5)
        self.assertEqual(test_p.shape, (2, 5))
        npt.assert_array_equal(np.arange(5, 50, 10), test_c)

    def test_subsample_paired_power_min_observations_error(self):
        with self.assertRaises(ValueError):
            subsample_paired_power(self.f,
                                   self.meta,
                                   cat=self.cat,
                                   control_cats=self.control_cats)

    def test_subsample_paired_power_interval_error(self):
        with self.assertRaises(ValueError):
            subsample_paired_power(self.f,
                                   self.meta,
                                   cat=self.cat,
                                   control_cats=self.control_cats,
                                   min_observations=2,
                                   min_counts=5,
                                   max_counts=7)

    def test_subsample_paired_power(self):
        known_c = np.array([1, 2, 3, 4, 5])
        # Sets up the handling values
        cat = 'INT'
        control_cats = ['SEX']
        # Tests for the control cats
        test_p, test_c = subsample_paired_power(self.meta_f,
                                                meta=self.meta,
                                                cat=cat,
                                                control_cats=control_cats,
                                                min_observations=1,
                                                counts_interval=1,
                                                num_iter=10,
                                                num_runs=2)
        # Test the output shapes are sane
        npt.assert_array_equal(test_p.shape, (2, 5))
        npt.assert_array_equal(known_c, test_c)

    def test__check_strs_str(self):
        self.assertTrue(_check_strs('string'))

    def test__check_strs_num(self):
        self.assertTrue(_check_strs(4.2))

    def test__check_str_nan(self):
        self.assertFalse(_check_strs(np.nan))

    def test__check_str_error(self):
        with self.assertRaises(TypeError):
            _check_strs(self.f)

    def test_confidence_bound_default(self):
        # Sets the know confidence bound
        known = 2.2830070
        test = confidence_bound(self.s1)
        npt.assert_almost_equal(test, known, 3)

    def test_confidence_bound_df(self):
        known = 2.15109
        test = confidence_bound(self.s1, df=15)
        npt.assert_almost_equal(known, test, 3)

    def test_confidence_bound_alpha(self):
        known = 3.2797886
        test = confidence_bound(self.s1, alpha=0.01)
        npt.assert_almost_equal(known, test, 3)

    def test_confidence_bound_nan(self):
        # Sets the value to test
        samples = np.array([[4, 3.2, 3.05],
                            [2, 2.8, 2.95],
                            [5, 2.9, 3.07],
                            [1, 3.1, 2.93],
                            [3, np.nan, 3.00]])
        # Sets the know value
        known = np.array([2.2284, 0.2573, 0.08573])
        # Tests the function
        test = confidence_bound(samples, axis=0)
        npt.assert_almost_equal(known, test, 3)

    def test_confidence_bound_axis_none(self):
        # Sets the value to test
        samples = np.array([[4, 3.2, 3.05],
                            [2, 2.8, 2.95],
                            [5, 2.9, 3.07],
                            [1, 3.1, 2.93],
                            [3, np.nan, 3.00]])
        # Sest the known value
        known = 0.52852
        # Tests the output
        test = confidence_bound(samples, axis=None)
        npt.assert_almost_equal(known, test, 3)

    def test__calculate_power(self):
        # Sets up the values to test
        crit = 0.025
        # Sets the known value
        known = 0.5
        # Calculates the test value
        test = _calculate_power(self.alpha, crit)
        # Checks the test value
        npt.assert_almost_equal(known, test)

    def test__compare_distributions_mode_error(self):
        with self.assertRaises(ValueError):
            _compare_distributions(self.f, self.samps, mode='fig')

    def test__compare_distributions_count_error(self):
        with self.assertRaises(ValueError):
            _compare_distributions(self.f, self.samps, counts=[1, 2, 3],
                                   num_iter=100)

    def test__compare_distributions_all_mode(self):
        known = np.ones((100))*0.0026998
        test = _compare_distributions(self.f, self.samps, num_iter=100)
        npt.assert_allclose(known, test, 5)

    def test__compare_distributions_matched_mode(self):
        # Sets the known value
        known_mean = 0.162195
        known_std = 0.121887
        known_shape = (100,)
        # Sets the sample value
        # Tests the sample value
        test = _compare_distributions(self.f, self.pop, mode='matched',
                                      num_iter=100)
        npt.assert_allclose(known_mean, test.mean(), rtol=0.1, atol=0.02)
        npt.assert_allclose(known_std, test.std(), rtol=0.1, atol=0.02)
        self.assertEqual(known_shape, test.shape)

    def test__calculate_power_curve_ratio_error(self):
        with self.assertRaises(ValueError):
            _calculate_power_curve(self.f, self.pop, self.num_samps,
                                   ratio=np.array([0.1, 0.2, 0.3]),
                                   num_iter=100)

    def test__calculate_power_curve_default(self):
        # Sets the know output
        known = np.array([0.509, 0.822, 0.962, 0.997, 1.000, 1.000, 1.000,
                          1.000,  1.000])

        # Generates the test values.
        test = _calculate_power_curve(self.f,
                                      self.pop,
                                      self.num_samps,
                                      num_iter=100)
        # Checks the samples returned sanely
        npt.assert_allclose(test, known, rtol=0.1, atol=0.1)

    def test__calculate_power_curve_alpha(self):
        # Sets the know output
        known = np.array([0.31, 0.568, 0.842, 0.954, 0.995, 1.000, 1.000,
                          1.000, 1.000])

        # Generates the test values
        test = _calculate_power_curve(self.f,
                                      self.pop,
                                      self.num_samps,
                                      alpha=0.01,
                                      num_iter=100)

        # Checks the samples returned sanely
        npt.assert_allclose(test, known, rtol=0.1, atol=0.1)

    def test__calculate_power_curve_ratio(self):
        # Sets the know output
        known = np.array([0.096, 0.333, 0.493, 0.743, 0.824, 0.937, 0.969,
                          0.996, 0.998])

        # Generates the test values
        test = _calculate_power_curve(self.f,
                                      self.pop,
                                      self.num_samps,
                                      ratio=np.array([0.25, 0.75]),
                                      num_iter=100)

        # Checks the samples returned sanely
        npt.assert_allclose(test, known, rtol=0.1, atol=0.1)

    def test_bootstrap_power_curve(self):
        # Sets the known values
        known_mean = np.array([0.500, 0.82, 0.965, 0.995, 1.000, 1.000,
                               1.000, 1.000,  1.000])
        known_bound = np.array([0.03, 0.02, 0.01, 0.01, 0.00, 0.00, 0.00, 0.00,
                                0.00])
        # Generates the test values
        test_mean, test_bound = bootstrap_power_curve(self.f,
                                                      self.pop,
                                                      self.num_samps,
                                                      num_iter=100)
        # Checks the function returned sanely
        npt.assert_allclose(test_mean, known_mean, rtol=0.05, atol=0.05)
        npt.assert_allclose(test_bound, known_bound, rtol=0.1, atol=0.01)

    def test_paired_subsamples_default(self):
        # Sets the known np.array set
        known_array = [sorted(['MM', 'SR', 'TS', 'GW', 'PP', 'WM']),
                       sorted(['CD', 'LF', 'PC', 'CB', 'MH', 'NR'])]

        # Gets the test value
        cat = 'INT'
        control_cats = ['SEX', 'AGE']
        test_array = paired_subsamples(self.meta, cat, control_cats)
        test_array[0] = sorted(test_array[0])
        test_array[1] = sorted(test_array[1])
        npt.assert_array_equal(known_array, test_array)

    def test_paired_subsamples_break(self):
        # Sets known np.array set
        known_array = [np.array([]), np.array([])]
        # Gets the test value
        cat = 'ABX'
        control_cats = ['SEX', 'AGE', 'INT']
        test_array = paired_subsamples(self.meta, cat, control_cats)
        npt.assert_array_equal(known_array, test_array)

    def test_paired_subsample_fewer(self):
        # Set known value
        known_array = {'PP', 'MH', 'CD', 'PC', 'TS', 'MM'}
        # Sets up test values
        cat = 'AGE'
        order = ['30s', '40s']
        control_cats = ['ABX']
        test_array = paired_subsamples(self.meta, cat, control_cats,
                                       order=order)
        for v in test_array[1]:
            self.assertTrue(v in known_array)
        for v in test_array[1]:
            self.assertTrue(v in known_array)

    def test_paired_subsamples_not_strict(self):
        known_array = [sorted(['WM', 'MM', 'GW', 'SR', 'TS']),
                       sorted(['LF', 'PC', 'CB', 'NR', 'CD'])]

        # Gets the test values
        cat = 'INT'
        control_cats = ['ABX', 'AGE']
        test_array = paired_subsamples(self.meta, cat, control_cats,
                                       strict_match=False)
        test_array[0] = sorted(test_array[0])
        test_array[1] = sorted(test_array[1])
        npt.assert_array_equal(known_array, test_array)

if __name__ == '__main__':
    main()
