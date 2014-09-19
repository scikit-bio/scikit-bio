#!/usr/bin/env python
# test_TaxPlot.py

from __future__ import division
from unittest import TestCase, main
import numpy as np
import numpy.testing as npt
import pandas as pd
from scipy.stats import kruskal
from skbio.stats.power import (get_subsampled_power,
                               _check_strs,
                               confidence_bound,
                               _calculate_power,
                               _compare_distributions,
                               _calculate_power_curve,
                               bootstrap_power_curve,
                               get_significant_subsample,
                               get_paired_subsamples)


class PowerAnalysisTest(TestCase):

    def setUp(self):
        """Initializes data for each test instance"""
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
        self.labels = np.array(['Age', 'INTervenption', 'Anptibiotics'])
        self.cats = np.array(['AGE', 'INT', 'ABX'])
        self.cat = "AGE"

    def test_get_subsampled_power_paired_error_meta_none(self):
        """Checks get_subsampled_power errors in PAIRED when missing meta"""
        with self.assertRaises(ValueError):
            get_subsampled_power("PAIRED", self.test_meta, cat=self.cat,
                                 control_cats=self.cats)

    def test_get_subsampled_power_paired_error_cat_none(self):
        """Checks get_subsampled_power errors in PAIRED when missing cat"""
        with self.assertRaises(ValueError):
            get_subsampled_power("PAIRED", self.test_meta, meta=self.meta,
                                 control_cats=self.cats)

    def test_get_subsampled_power_paired_error_ctrl_cats_none(self):
        """Checks get_subsampled_power errors in PAIRED when missing cat"""
        with self.assertRaises(ValueError):
            get_subsampled_power("PAIRED", self.test_meta, meta=self.meta,
                                 cat=self.cat)

    def test_get_subsampled_power_sig_error_no_samples(self):
        """Checks get_subsampled_power error when no samples in SIG mode"""
        with self.assertRaises(ValueError):
            get_subsampled_power("SIG", self.f)

    def test_get_subsampled_power_all_error_no_samples(self):
        """Checks get_subsampled_power error when no samples in ALL mode"""
        with self.assertRaises(ValueError):
            get_subsampled_power("ALL", self.f)

    def test_get_subsampled_power_bad_mode(self):
        """Checks get_subsampled_power errors when the mode is not supported"""
        with self.assertRaises(ValueError):
            get_subsampled_power("foo", self.f)

    def test_get_subsampled_power_min_counts_error(self):
        """Checks get_subsampled_power errors when there are not enough samples
        """
        with self.assertRaises(RuntimeError):
            get_subsampled_power('ALL', self.f, samples=[np.ones((2)),
                                 np.ones((5))])

    def test_get_subsampled_power_interval_error(self):
        """Checks get_subsampled_power errors when counts_interval is too big
        """
        with self.assertRaises(RuntimeError):
            get_subsampled_power('ALL', self.f, samples=[np.ones((2)),
                                 np.ones((5))], counts_start=5, max_counts=7)

    def test_get_subsampled_power_paired(self):
        """Checks get_subsampled_power generates a reasonable subsample"""
        known_c = np.array([1, 2, 3, 4, 5])
        # Sets up the handling values
        cat = 'INT'
        control_cats = ['SEX']
        # Tests for the control cats
        test_p, test_c = get_subsampled_power("PAIRED", self.meta_f,
                                              meta=self.meta,
                                              cat=cat,
                                              control_cats=control_cats,
                                              min_counts=1,
                                              counts_interval=1,
                                              num_iter=10,
                                              num_runs=2)
        # Test the output shapes are sane
        npt.assert_array_equal(test_p.shape, (2, 5))
        npt.assert_array_equal(known_c, test_c)

    def test_get_subsampled_power_all_samples(self):
        """Checks get_unpaired_power handles all samples correctly"""
        test_p, test_c = get_subsampled_power('ALL', self.f, samples=self.pop,
                                              num_iter=10, num_runs=2,
                                              counts_start=5)
        self.assertEqual(test_p.shape, (2, 5))
        npt.assert_array_equal(np.arange(5, 50, 10), test_c)

    def test_get_get_subsampled_power_significant_samples(self):
        """Checks get_unpaired_power handles all samples correctly"""
        test_p, test_c = get_subsampled_power("SIG", self.f, samples=self.pop,
                                              num_iter=10, num_runs=2)
        self.assertEqual(test_p.shape, (2, 4))
        npt.assert_array_equal(np.arange(10, 50, 10), test_c)

    def test__check_strs_str(self):
        """Test check_strs returns sanely when passed a string"""
        self.assertTrue(_check_strs('string'))

    def test__check_strs_num(self):
        """Tests check_strs returns sanely when passed a number"""
        self.assertTrue(_check_strs(4.2))

    def test__check_str_nan(self):
        """Tests check_strs retruns sanely when passed a np.nan"""
        self.assertFalse(_check_strs(np.nan))

    def test__check_str_error(self):
        """Tests check_strs errors when not passed a string or number"""
        with self.assertRaises(TypeError):
            _check_strs(self.f)

    def test_confidence_bound_default(self):
        """Checks confidence_bound correctly determines an INTerval"""
        # Sets the know confidence bound
        known = 2.2830070
        test = confidence_bound(self.s1)
        npt.assert_almost_equal(test, known, 3)

    def test_confidence_bound_df(self):
        """Checks a custom df for confidence_bound"""
        known = 2.15109
        test = confidence_bound(self.s1, df=15)
        npt.assert_almost_equal(known, test, 3)

    def test_confidence_bound_alpha(self):
        """Checks a custom df for confidence_bound"""
        known = 3.2797886
        test = confidence_bound(self.s1, alpha=0.01)
        npt.assert_almost_equal(known, test, 3)

    def test_confidence_bound_nan(self):
        """Tests confidence_bound can handle np.nans"""
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
        """Tests confidence_bound can handle a None axis"""
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
        """Tests calculate_power is sane"""
        # Sets up the values to test
        crit = 0.025
        # Sets the known value
        known = 0.5
        # Calculates the test value
        test = _calculate_power(self.alpha, crit)
        # Checks the test value
        npt.assert_almost_equal(known, test)

    def test__compare_distributions_counpt_error(self):
        """Checks error is raised when there is not a counpt for each group"""
        with self.assertRaises(ValueError):
            _compare_distributions(self.f, self.samps, counts=[1, 2, 3],
                                   num_iter=100)

    def test__compare_distributions(self):
        """Checks _compare_distributions is sane"""
        known = np.ones((100))*0.0026998
        test = _compare_distributions(self.f, self.samps, num_iter=100)
        npt.assert_allclose(known, test, 5)

    def test__calculate_power_curve_ratio_error(self):
        """Checks the function errors correctly for a non-sane ratio argumenpt

        the ratio argumenpt must be a none-type (flat distribution), or an
        array with the same shape as population.

        """
        with self.assertRaises(ValueError):
            _calculate_power_curve(self.f, self.pop, self.num_samps,
                                   ratio=np.array([0.1, 0.2, 0.3]),
                                   num_iter=100)

    def test__calculate_power_curve_default(self):
        """Checks the power np.array is within a sane range for default values
        """
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
        """Checks the power np.array is in a sane range when alpha is varied"""
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
        """Checks the power np.array is in a sane range when ratio is varied"""
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
        """Checks the power estimate is handled in a sane way"""
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

    def test_get_significant_subsample_no_tests(self):
        """Checks get_significant_subsample errors when inputs are too similar
        """
        with self.assertRaises(RuntimeError):
            get_significant_subsample([None], self.samps, num_rounds=100)

    def test_get_significant_subsample_no_results(self):
        """Checks get_significant_subsample errors when inputs are too similar
        """
        # Sets up a function which will fail testing
        def test_f(x):
            if len(x[0]) == 100:
                return 0.001
            else:
                return 0.5
        # Tests a value error is raised
        with self.assertRaises(RuntimeError):
            get_significant_subsample([test_f], self.samps, sub_size=10,
                                      num_rounds=5)

    def test_get_signfigianpt_subsample_no_iteration(self):
        """Checks get_significant_subsample errors when iteration is not found
        """
        # Sets up a function which will fail testing
        def test_f(x):
            if len(x[0]) > 5:
                return 0.001
            else:
                return 0.5
        # Tests if a RuntimeError is raised
        with self.assertRaises(RuntimeError):
            get_significant_subsample([test_f], self.samps, sub_size=5,
                                      num_rounds=5)

    def test_get_significant_subsample_default(self):
        """Checks get_significant_subsample functions sanely under defaults"""
        pop = [np.arange(0, 10, 1), np.arange(0, 20, 0.2)]
        # Checks the overall data meets the parameters
        self.assertNotEqual(len(pop[0]), len(pop[1]))
        # Generates subsamples
        test_ids = get_significant_subsample([self.f], pop)
        # Checks the results
        self.assertEqual(len(test_ids[0]), len(test_ids[1]))
        self.assertTrue(self.f(test_ids) < 0.05)

    def test_get_paired_subsamples_default(self):
        """Checks controlled subsets can be generated sanely"""
        # Sets the known np.array set
        known_array = [sorted(['MM', 'SR', 'TS', 'GW', 'PP', 'WM']),
                       sorted(['CD', 'LF', 'PC', 'CB', 'MH', 'NR'])]

        # Gets the test value
        cat = 'INT'
        control_cats = ['SEX', 'AGE']
        test_array = get_paired_subsamples(self.meta, cat, control_cats)
        test_array[0] = sorted(test_array[0])
        test_array[1] = sorted(test_array[1])
        npt.assert_array_equal(known_array, test_array)

    def test_get_paired_subsamples_break(self):
        """Checks controlled susbets can skip sanely when there are no matches
        """
        # Sets known np.array set
        known_array = [np.array([]), np.array([])]
        # Gets the test value
        cat = 'ABX'
        control_cats = ['SEX', 'AGE', 'INT']
        test_array = get_paired_subsamples(self.meta, cat, control_cats)
        npt.assert_array_equal(known_array, test_array)

    def test_get_paired_subsample_fewer(self):
        """Checks controlled subsets can handle fewer samples sanely"""
        # Set known value
        known_array1 = {'PP', 'MH'}
        known_array2 = {'CD', 'PC', 'TS', 'MM'}
        # Sets up test values
        cat = 'AGE'
        order = ['30s', '40s']
        control_cats = ['ABX']
        test_array = get_paired_subsamples(self.meta, cat, control_cats,
                                           order=order)
        for v in test_array[1]:
            self.assertTrue(v in known_array1)
        for v in test_array[1]:
            self.assertTrue(v in known_array2)

    def test_get_paired_subsamples_not_strict(self):
        """Checks controlled subsets can be generated with missing values"""
        known_array = [sorted(['WM', 'MM', 'GW', 'SR', 'TS']),
                       sorted(['LF', 'PC', 'CB', 'NR', 'CD'])]

        # Gets the test values
        cat = 'INT'
        control_cats = ['ABX', 'AGE']
        test_array = get_paired_subsamples(self.meta, cat, control_cats,
                                           strict=False)
        test_array[0] = sorted(test_array[0])
        test_array[1] = sorted(test_array[1])
        npt.assert_array_equal(known_array, test_array)

if __name__ == '__main__':
    main()
