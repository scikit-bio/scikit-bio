# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd
from scipy.stats import kruskal

from skbio.stats.power import (subsample_power,
                               subsample_paired_power,
                               _check_nans,
                               confidence_bound,
                               _calculate_power,
                               _compare_distributions,
                               _calculate_power_curve,
                               _check_subsample_power_inputs,
                               _identify_sample_groups,
                               _draw_paired_samples,
                               _get_min_size,
                               paired_subsamples
                               )


class PowerAnalysisTest(TestCase):

    def setUp(self):

        # Defines a testing functions
        def test_meta(ids, meta, cat, div):
            """Checks thhe div metric with a kruskal wallis"""
            out = [meta.loc[id_, div] for id_ in ids]
            return kruskal(*out)[1]

        def meta_f(x):
            """Applies `test_meta` to a result"""
            return test_meta(x, self.meta, 'INT', 'DIV')

        def f(x):
            """returns the p value of a kruskal wallis test"""
            return kruskal(*x)[1]

        self.test_meta = test_meta
        self.f = f
        self.meta_f = meta_f
        self.num_p = 1

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
        self.meta_pairs = {0: [['GW', 'SR', 'TS'], ['CB', 'LF', 'PC']],
                           1: [['MM', 'PP', 'WM'], ['CD', 'MH', 'NR']]}
        self.pair_index = np.array([0, 0, 0, 1, 1, 1])
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

    def test_subsample_power_defaults(self):
        test_p, test_c = subsample_power(self.f, self.pop,
                                         num_iter=10, num_runs=5)
        self.assertEqual(test_p.shape, (5, 4))
        npt.assert_array_equal(np.array([10, 20, 30, 40]), test_c)

    def test_subsample_power_counts(self):
        test_p, test_c = subsample_power(self.f,
                                         samples=self.pop,
                                         num_iter=10,
                                         num_runs=2,
                                         min_counts=5)
        self.assertEqual(test_p.shape, (2, 5))
        npt.assert_array_equal(np.arange(5, 50, 10), test_c)

    def test_subsample_power_matches(self):
        test_p, test_c = subsample_power(self.f,
                                         samples=self.pop,
                                         num_iter=10,
                                         num_runs=5,
                                         draw_mode="matched")
        self.assertEqual(test_p.shape, (5, 4))
        npt.assert_array_equal(np.array([10, 20, 30, 40]), test_c)

    def test_subsample_power_multi_p(self):
        test_p, test_c = subsample_power(lambda x: np.array([0.5, 0.5]),
                                         samples=self.pop,
                                         num_iter=10,
                                         num_runs=5)
        self.assertEqual(test_p.shape, (5, 4, 2))
        npt.assert_array_equal(np.array([10, 20, 30, 40]), test_c)

    def test_subsample_paired_power(self):
        known_c = np.array([1, 2, 3, 4])
        # Sets up the handling values
        cat = 'INT'
        control_cats = ['SEX']

        # Tests for the control cats
        test_p, test_c = subsample_paired_power(self.meta_f,
                                                meta=self.meta,
                                                cat=cat,
                                                control_cats=control_cats,
                                                counts_interval=1,
                                                num_iter=10,
                                                num_runs=2)
        # Test the output shapes are sane
        self.assertEqual(test_p.shape, (2, 4))
        npt.assert_array_equal(known_c, test_c)

    def test_subsample_paired_power_multi_p(self):
        def f(x):
            return np.array([0.5, 0.5, 0.005])
        cat = 'INT'
        control_cats = ['SEX']
        # Tests for the control cats
        test_p, test_c = subsample_paired_power(f,
                                                meta=self.meta,
                                                cat=cat,
                                                control_cats=control_cats,
                                                counts_interval=1,
                                                num_iter=10,
                                                num_runs=2)
        self.assertEqual(test_p.shape, (2, 4, 3))

    def test_check_nans_str(self):
        self.assertTrue(_check_nans('string'))

    def test_check_nans_num(self):
        self.assertTrue(_check_nans(4.2))

    def test__check_nans_nan(self):
        self.assertFalse(_check_nans(np.nan))

    def test__check_nans_clean_list(self):
        self.assertTrue(_check_nans(['foo', 'bar'], switch=True))

    def test__check_nans_list_nan(self):
        self.assertFalse(_check_nans(['foo', np.nan], switch=True))

    def test__check_str_error(self):
        with self.assertRaises(TypeError):
            _check_nans(self.f)

    def test__get_min_size_strict(self):
        known = 5
        test = _get_min_size(self.meta, 'INT', ['ABX', 'SEX'], ['Y', 'N'],
                             True)
        self.assertEqual(test, known)

    def test__get_min_size_relaxed(self):
        known = 5
        test = _get_min_size(self.meta, 'INT', ['ABX', 'SEX'], ['Y', 'N'],
                             False)
        self.assertEqual(known, test)

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

    def test__calculate_power_n(self):
        crit = 0.025
        known = np.array([0.5, 0.5])
        alpha = np.vstack((self.alpha, self.alpha))
        test = _calculate_power(alpha, crit)
        npt.assert_almost_equal(known, test)

    def test__compare_distributions_sample_counts_error(self):
        with self.assertRaises(ValueError):
            _compare_distributions(self.f, [self.pop[0][:5], self.pop[1]], 1,
                                   counts=25)

    def test__compare_distributions_all_mode(self):
        known = np.ones((100))*0.0026998
        test = _compare_distributions(self.f, self.samps, 1, num_iter=100)
        npt.assert_allclose(known, test, 5)

    def test__compare_distributions_matched_mode(self):
        # Sets the known value
        known_mean = 0.162195
        known_std = 0.121887
        known_shape = (100,)
        # Tests the sample value
        test = _compare_distributions(self.f, self.pop, self.num_p,
                                      mode='matched', num_iter=100)
        npt.assert_allclose(known_mean, test.mean(), rtol=0.1, atol=0.02)
        npt.assert_allclose(known_std, test.std(), rtol=0.1, atol=0.02)
        self.assertEqual(known_shape, test.shape)

    def test__compare_distributions_draw_mode(self):
        draw_mode = 'Ultron'
        with self.assertRaises(ValueError):
            _check_subsample_power_inputs(self.f, self.pop, draw_mode,
                                          self.num_p)

    def test__compare_distributions_multiple_returns(self):
        known = np.array([[1, 1, 1], [2, 2, 2], [3, 3, 3]])

        def f(x):
            return np.array([1, 2, 3])

        test = _compare_distributions(f, self.pop, 3, mode='matched',
                                      num_iter=3)
        npt.assert_array_equal(known, test)

    def test_check_subsample_power_inputs_matched_mode(self):
        with self.assertRaises(ValueError):
            _check_subsample_power_inputs(self.f,
                                          samples=[np.ones((2)), np.ones((5))],
                                          draw_mode="matched")

    def test_check_subsample_power_inputs_counts(self):
        with self.assertRaises(ValueError):
            _check_subsample_power_inputs(self.f,
                                          samples=[np.ones((3)), np.ones((5))],
                                          min_counts=5,
                                          counts_interval=1000,
                                          max_counts=7)

    def test_check_subsample_power_inputs_ratio(self):
        with self.assertRaises(ValueError):
            _check_subsample_power_inputs(self.f,
                                          self.samps,
                                          ratio=np.array([1, 2, 3]))

    def test_check_subsample_power_inputs_test(self):
        # Defines a test function
        def test(x):
            return 'Hello World!'
        with self.assertRaises(TypeError):
            _check_subsample_power_inputs(test, self.samps)

    def test_check_sample_power_inputs(self):
        # Defines the know returns
        known_num_p = 1
        known_ratio = np.ones((2))
        known_counts = np.arange(2, 10, 2)
        # Runs the code for the returns
        test_ratio, test_num_p, test_counts = \
            _check_subsample_power_inputs(self.f,
                                          self.samps,
                                          counts_interval=2,
                                          max_counts=10)
        # Checks the returns are sane
        self.assertEqual(known_num_p, test_num_p)
        npt.assert_array_equal(known_ratio, test_ratio)
        npt.assert_array_equal(known_counts, test_counts)

    def test__calculate_power_curve_ratio_error(self):
        with self.assertRaises(ValueError):
            _calculate_power_curve(self.f, self.pop, self.num_samps,
                                   ratio=np.array([0.1, 0.2, 0.3]),
                                   num_iter=100)

    def test__calculate_power_curve_default(self):
        # Sets the known output
        known = np.array([0.509, 0.822, 0.962, 0.997, 1.000, 1.000, 1.000,
                          1.000, 1.000])
        # Generates the test values
        test = _calculate_power_curve(self.f,
                                      self.pop,
                                      self.num_samps,
                                      num_iter=100)
        # Checks the samples returned sanely
        npt.assert_allclose(test, known, rtol=0.1, atol=0.01)

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

    def test_paired_subsamples_default(self):
        # Sets the known np.array set
        known_array = [{'MM', 'SR', 'TS', 'GW', 'PP', 'WM'},
                       {'CD', 'LF', 'PC', 'CB', 'MH', 'NR'}]

        # Gets the test value
        cat = 'INT'
        control_cats = ['SEX', 'AGE']
        test_array = paired_subsamples(self.meta, cat, control_cats)
        self.assertEqual(known_array[0], set(test_array[0]))
        self.assertEqual(known_array[1], set(test_array[1]))

    def test_paired_subsamples_break(self):
        # Sets known np.array set
        known_array = [np.array([]), np.array([])]
        # Gets the test value
        cat = 'ABX'
        control_cats = ['SEX', 'AGE', 'INT']
        test_array = paired_subsamples(self.meta, cat, control_cats)
        npt.assert_array_equal(known_array, test_array)

    def test_paired_subsample_undefined(self):
        known_array = np.zeros((2, 0))
        cat = 'INT'
        order = ['Y', 'N']
        control_cats = ['AGE', 'ABX', 'SEX']
        test_array = paired_subsamples(self.meta, cat, control_cats,
                                       order=order)
        npt.assert_array_equal(test_array, known_array)

    def test_paired_subsample_fewer(self):
        # Set known value
        known_array = {'PP', 'MH', 'CD', 'PC', 'TS', 'MM'}
        # Sets up test values
        cat = 'AGE'
        order = ['30s', '40s']
        control_cats = ['ABX']
        test_array = paired_subsamples(self.meta, cat, control_cats,
                                       order=order)
        for v in test_array[0]:
            self.assertTrue(v in known_array)
        for v in test_array[1]:
            self.assertTrue(v in known_array)

    def test_paired_subsamples_not_strict(self):
        known_array = [{'WM', 'MM', 'GW', 'SR', 'TS'},
                       {'LF', 'PC', 'CB', 'NR', 'CD'}]

        # Gets the test values
        cat = 'INT'
        control_cats = ['ABX', 'AGE']
        test_array = paired_subsamples(self.meta, cat, control_cats,
                                       strict_match=False)
        self.assertEqual(set(test_array[0]), known_array[0])
        self.assertEqual(set(test_array[1]), known_array[1])

    def test__identify_sample_groups(self):
        # Defines the know values
        known_pairs = {0: [['MM'], ['CD']],
                       1: [['SR'], ['LF']],
                       2: [['TS'], ['PC']],
                       3: [['GW'], ['CB']],
                       4: [['PP'], ['MH']],
                       5: [['WM'], ['NR']]}
        known_index = np.array([0, 1, 2, 3, 4, 5])
        test_pairs, test_index = _identify_sample_groups(self.meta,
                                                         'INT',
                                                         ['SEX', 'AGE'],
                                                         order=['N', 'Y'],
                                                         strict_match=True)
        self.assertEqual(known_pairs.keys(), test_pairs.keys())
        self.assertEqual(sorted(known_pairs.values()),
                         sorted(test_pairs.values()))
        npt.assert_array_equal(known_index, test_index)

    def test__identify_sample_groups_not_strict(self):
        # Defines the know values
        known_pairs = {1: [np.array(['PP'], dtype=object),
                           np.array(['CD', 'NR'], dtype=object)],
                       0: [np.array(['MM', 'WM'], dtype=object),
                           np.array(['MH'], dtype=object)],
                       2: [np.array(['GW'], dtype=object),
                           np.array(['CB'], dtype=object)]}
        known_index = np.array([0, 1, 2])
        test_pairs, test_index = _identify_sample_groups(self.meta,
                                                         'INT',
                                                         ['SEX', 'ABX'],
                                                         order=['N', 'Y'],
                                                         strict_match=False)
        self.assertEqual(known_pairs.keys(), test_pairs.keys())

        for k in known_pairs:
            for i in range(2):
                npt.assert_array_equal(known_pairs[k][i], test_pairs[k][i])
        npt.assert_array_equal(known_index, test_index)

    def test__draw_paired_samples(self):
        num_samps = 3
        known_sets = [{'GW', 'SR', 'TS', 'MM', 'PP', 'WM'},
                      {'CB', 'LF', 'PC', 'CD', 'MH', 'NR'}]
        test_samps = _draw_paired_samples(self.meta_pairs, self.pair_index,
                                          num_samps)
        for i, t in enumerate(test_samps):
            self.assertTrue(set(t).issubset(known_sets[i]))


if __name__ == '__main__':
    main()
