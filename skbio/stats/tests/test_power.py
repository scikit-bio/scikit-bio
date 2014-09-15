#!/usr/bin/env python
# test_TaxPlot.py

from __future__ import division
from unittest import TestCase, main
from numpy import ones, power, array, arange, nan, isnan
import numpy.random
from numpy.testing import (assert_almost_equal,
                           assert_array_equal,
                           assert_allclose)
from pandas import DataFrame
from scipy.stats import kruskal
from matplotlib.figure import Figure
from skbio.stats.power import (_check_strs,
                               _confidence_bound,
                               _calculate_power,
                               compare_distributions,
                               calculate_power_curve,
                               bootstrap_power_curve,
                               get_signifigant_subsample,
                               get_paired_subsamples,
                               plot_effects,
                               collate_effect_size,
                               get_paired_effect,
                               get_unpaired_effect)


def test_meta(ids, meta, cat, div):
    """Checks thhe div metric with a kruskal wallis"""
    out = [meta.loc[id_, div] for id_ in ids]
    return kruskal(*out)[1]


class PowerAnalysisTest(TestCase):

    def setUp(self):
        """Initializes data for each test instance"""
        # Sets the random seed
        numpy.random.seed(5)
        # Sets up the distributions of data for use
        self.s1 = arange(0, 10, 1)
        # Sets up two distributions which will never be equal by a rank-sum
        # test.
        self.samps = [ones((10))/10., ones((10))]
        self.pop = [arange(0, 10, 0.1), arange(0, 20, 0.2)]
        # Sets up a vector of alpha values
        self.alpha = power(10, array([-1, -1.301, -2, -3])).round(3)
        # Sets up a vector of samples
        self.num_samps = arange(10, 100, 10)
        # Sets up the test function, a rank-sum test
        self.f = lambda x: kruskal(*x)[1]
        # Sets up a mapping file
        meta = {'NR': {'RANGE': 'M', 'SEX': 'F', 'AGE': nan,   'ABX': 'Y', 'DAMAGE': 150},
                'MH': {'RANGE': 'L', 'SEX': 'F', 'AGE': '30s', 'ABX': 'Y', 'DAMAGE': 100},
                'PP': {'RANGE': 'M', 'SEX': 'F', 'AGE': '30s', 'ABX': 'N', 'DAMAGE': 250},
                'CD': {'RANGE': 'L', 'SEX': 'F', 'AGE': '30s', 'ABX': 'Y', 'DAMAGE': 100},
                'MM': {'RANGE': 'C', 'SEX': 'F', 'AGE': '30s', 'ABX': 'Y', 'DAMAGE': 500},
                'SW': {'RANGE': 'M', 'SEX': 'M', 'AGE': '30s', 'ABX': 'N', 'DAMAGE': 050},
                'TS': {'RANGE': 'M', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y', 'DAMAGE': 750},
                'CB': {'RANGE': 'L', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y', 'DAMAGE': 250},
                'BB': {'RANGE': 'C', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y', 'DAMAGE': 900},
                'WS': {'RANGE': 'M', 'SEX': 'M', 'AGE': nan,   'ABX': 'N', 'DAMAGE': 200},
                'SR': {'RANGE': 'C', 'SEX': 'M', 'AGE': nan,   'ABX': 'N', 'DAMAGE': 300},
                'NF': {'RANGE': 'L', 'SEX': 'M', 'AGE': '60s', 'ABX': 'N', 'DAMAGE': 400}}
        self.meta = DataFrame.from_dict(meta, orient='index')
        self.meta_f = lambda x: test_meta(x, self.meta, 'RANGE', 'DAMAGE')
        self.counts = array([5, 15, 25, 35, 45])
        self.powers = [array([[0.105, 0.137, 0.174, 0.208, 0.280],
                              [0.115, 0.135, 0.196, 0.204, 0.281],
                              [0.096, 0.170, 0.165, 0.232, 0.256],
                              [0.122, 0.157, 0.202, 0.250, 0.279],
                              [0.132, 0.135, 0.173, 0.203, 0.279]]),
                       array([[0.157, 0.345, 0.522, 0.639, 0.739],
                              [0.159, 0.374, 0.519, 0.646, 0.757],
                              [0.161, 0.339, 0.532, 0.634, 0.745],
                              [0.169, 0.372, 0.541, 0.646, 0.762],
                              [0.163, 0.371, 0.522, 0.648, 0.746]]),
                       array([[0.276, 0.626, 0.865, 0.927, 0.992],
                              [0.267, 0.667, 0.848, 0.937, 0.978],
                              [0.236, 0.642, 0.850, 0.935, 0.977],
                              [0.249, 0.633, 0.828, 0.955, 0.986],
                              [0.249, 0.663, 0.869, 0.951, 0.985]])]
        self.power_alpha = 0.1
        self.effects = array([0.15245, 0.34877, 0.55830])
        self.bounds = array([0.01027, 0.00293, 0.00734])
        self.labels = array(['Range', 'Sex', 'Abx'])

    def test_get_paired_effect_runtime_error(self):
        """Checks get_paired_effect generates a runtime error correctly"""
        # Sets up values for handling the data
        cat = 'RANGE'
        control_cats = ['SEX', 'AGE', 'ABX']
        # Checks the error is raised
        with self.assertRaises(RuntimeError):
            get_paired_effect(self.f, self.meta, cat, control_cats)

    def test_get_paired_effects(self):
        """Checks get_paired_effect generates a reasonably sized subsample"""
        # Sets up the handling values
        cat = 'RANGE'
        control_cats = ['SEX']
        # Tests for the control cats
        test_p, test_c = get_paired_effect(self.meta_f, self.meta, cat,
                                           control_cats, min_counts=1,
                                           counts_interval=1, num_iter=10,
                                           num_runs=2)
        # Test the output shapes are sane
        self.assertEqual(test_p.shape, (2, 2))
        assert_array_equal(test_c, array([1, 2]))

    def test_get_unpaired_effect_runtime_error(self):
        """Checks get_unpaired_effect generates a runtime error correctly"""
        with self.assertRaises(RuntimeError):
            get_unpaired_effect('ALL', self.f, [ones((2)), ones((5))])

    def test_get_unpaired_effects_all_samples(self):
        """Checks get_unpaired_effect handles all samples correctly"""
        test_p, test_c = get_unpaired_effect('ALL', self.f, self.pop)
        self.assertEqual(test_p.shape, (10, 4))
        assert_array_equal(arange(10, 50, 10), test_c)

    def test_get_unpaired_effects_signifigant_samples(self):
        """Checks get_unpaired_effect handles all samples correctly"""
        test_p, test_c = get_unpaired_effect('SIGNIFICANT', self.f, self.pop)
        self.assertEqual(test_p.shape, (10, 4))
        assert_array_equal(arange(10, 50, 10), test_c)

    def test__check_strs_str(self):
        """Test check_strs returns sanely when passed a string"""
        self.assertTrue(_check_strs('string'))

    def test__check_strs_num(self):
        """Tests check_strs returns sanely when passed a number"""
        self.assertTrue(_check_strs(4))

    def test__check_str_nan(self):
        """Tests check_strs retruns sanely when passed a nan"""
        self.assertFalse(_check_strs(nan))

    def test__check_str_error(self):
        """Tests check_strs errors when not passed a string or number"""
        with self.assertRaises(TypeError):
            _check_strs(self.f)

    def test__confidence_bound_default(self):
        """Checks confidence_bound correctly determines an interval"""
        # Sets the know confidence bound
        known = 2.1658506
        test = _confidence_bound(self.s1)
        assert_almost_equal(test, known, 3)

    def test__confidence_bound_df(self):
        """Checks a custom df for confidence_bound"""
        known = 2.0407076
        test = _confidence_bound(self.s1, df=15)
        assert_almost_equal(known, test, 3)

    def test__confidence_bound_alpha(self):
        """Checks a custom df for confidence_bound"""
        known = 3.111481
        test = _confidence_bound(self.s1, alpha=0.01)
        assert_almost_equal(known, test, 3)

    def test__confidence_bound_nan(self):
        """Tests _confidence_bound can handle nans"""
        # Sets the value to test
        samples = array([[4, 3.2, 3.05],
                         [2, 2.8, 2.95],
                         [5, 2.9, 3.07],
                         [1, 3.1, 2.93],
                         [3, nan, 3.00]])
        # Sets the know value
        known = array([1.9931, 0.2228, 0.07668])
        # Tests the function
        test = _confidence_bound(samples, axis=0)
        assert_almost_equal(known, test, 3)

    def test__confidence_bound_axis_none(self):
        """Tests _confidence_bound can handle a None axis"""
        # Sets the value to test
        samples = array([[4, 3.2, 3.05],
                         [2, 2.8, 2.95],
                         [5, 2.9, 3.07],
                         [1, 3.1, 2.93],
                         [3, nan, 3.00]])
        # Sest the known value
        known = 0.50930
        # Tests the output
        test = _confidence_bound(samples, axis=None)
        assert_almost_equal(known, test, 3)

    def test__calculate_power(self):
        """Tests calculate_power is sane"""
        # Sets up the values to test
        crit = 0.025
        # Sets the known value
        known = 0.5
        # Calculates the test value
        test = _calculate_power(self.alpha, crit)
        # Checks the test value
        assert_almost_equal(known, test)

    def test_compare_distributions_count_error(self):
        """Checks error is raised when there is not a count for each group"""
        with self.assertRaises(ValueError):
            compare_distributions(self.f, self.samps, counts=[1, 2, 3])

    def test_compare_distributions(self):
        """Checks compare_distributions is sane"""
        known = ones((100))*0.0026998
        test = compare_distributions(self.f, self.samps, num_iter=100)
        assert_allclose(known, test, 5)

    def test_calculate_power_curve_ratio_error(self):
        """Checks the function errors correctly for a non-sane ratio argument

        the ratio argument must be a none-type (flat distribution), or an array
        with the same shape as population.

        """
        self.assertRaises(ValueError, calculate_power_curve, self.f,
                          self.pop, self.num_samps,
                          ratio=array([0.1, 0.2, 0.3]))

    def test_calculate_power_curve_default(self):
        """Checks the power array is within a sane range for default values"""
        # Sets the know output
        known = array([0.509, 0.822, 0.962, 0.997, 1.000, 1.000, 1.000,
                       1.000,  1.000])

        # Generates the test values.
        test = calculate_power_curve(self.f,
                                     self.pop,
                                     self.num_samps,
                                     num_iter=100)
        # Checks the samples returned sanely
        assert_allclose(test, known, rtol=0.1, atol=0.1)

    def test_calculate_power_curve_alpha(self):
        """Checks the power array is in a sane range when alpha is varied"""
        # Sets the know output
        known = array([0.31, 0.568, 0.842, 0.954, 0.995, 1.000, 1.000, 1.000,
                       1.000])

        # Generates the test values
        test = calculate_power_curve(self.f,
                                     self.pop,
                                     self.num_samps,
                                     alpha=0.01,
                                     num_iter=100)

        # Checks the samples returned sanely
        assert_allclose(test, known, rtol=0.1, atol=0.1)

    def test_calculate_power_curve_ratio(self):
        """Checks the power array is in a sane range when ratio is varied"""
        # Sets the know output
        known = array([0.096, 0.333, 0.493, 0.743, 0.824, 0.937, 0.969,
                       0.996, 0.998])

        # Generates the test values
        test = calculate_power_curve(self.f,
                                     self.pop,
                                     self.num_samps,
                                     ratio=array([0.25, 0.75]),
                                     num_iter=100)

        # Checks the samples returned sanely
        assert_allclose(test, known, rtol=0.1, atol=0.1)

    def test_bootstrap_power_curve(self):
        """Checks the power estimate is handled in a sane way"""
        # Sets the known values
        known_mean = array([0.500, 0.82, 0.965, 0.995, 1.000, 1.000,
                            1.000, 1.000,  1.000])
        known_bound = array([0.03, 0.02, 0.01, 0.01, 0.00, 0.00, 0.00, 0.00,
                             0.00])
        # Generates the test values
        test_mean, test_bound = bootstrap_power_curve(self.f,
                                                      self.pop,
                                                      self.num_samps,
                                                      num_iter=100)
        # Checks the function returned sanely
        assert_allclose(test_mean, known_mean, rtol=0.05, atol=0.05)
        assert_allclose(test_bound, known_bound, rtol=0.1, atol=0.01)

    def test_get_signifigant_subsample_no_tests(self):
        """Checks get_signifigant_subsample errors when inputs are too similar
        """
        with self.assertRaises(RuntimeError):
            get_signifigant_subsample([None], self.samps)

    def test_get_signifigant_subsample_no_results(self):
        """Checks get_signifigant_subsample errors when inputs are too similar
        """
        # Sets up a function which will fail testing
        def test_f(x):
            if len(x[0]) == 100:
                return 0.001
            else:
                return 0.5
        # Tests a value error is raised
        with self.assertRaises(RuntimeError):
            get_signifigant_subsample([test_f], self.samps, sub_size=10,
                                      num_rounds=5)

    def test_get_signfigiant_subsample_no_iteration(self):
        """Checks get_signifigant_subsample errors when iteration is not found
        """
        # Sets up a function which will fail testing
        def test_f(x):
            if len(x[0]) > 5:
                return 0.001
            else:
                return 0.5
        # Tests if a RuntimeError is raised
        with self.assertRaises(RuntimeError):
            get_signifigant_subsample([test_f], self.samps, sub_size=5,
                                      num_rounds=5)

    def test_get_signifigant_subsample_default(self):
        """Checks get_signifigant_subsample functions sanely under defaults"""
        pop = [arange(0, 10, 1), arange(0, 20, 0.2)]
        # Checks the overall data meets the parameters
        self.assertNotEqual(len(pop[0]), len(pop[1]))
        # Generates subsamples
        test_ids = get_signifigant_subsample([self.f], pop)
        # Checks the results
        self.assertEqual(len(test_ids[0]), len(test_ids[1]))
        self.assertTrue(self.f(test_ids) < 0.05)

    def test_get_paired_subsamples_default(self):
        """Checks controlled subsets can be generated sanely"""
        # Sets the known array set
        known_array = [array(['BB']), array(['CB']), array(['TS'])]
        # Gets the test value
        cat = 'RANGE'
        control_cats = ['SEX', 'AGE', 'ABX']
        test_array = get_paired_subsamples(self.meta, cat, control_cats)
        assert_array_equal(known_array, test_array)

    def test_get_paired_subsamples_skips(self):
        """Checks controlled susbets can skip sanely"""
        # Sets known array set
        known_array = [array(['PP']), array(['SW'])]
        # Gets the test value
        cat = 'SEX'
        control_cats = ['AGE', 'RANGE']
        test_array = get_paired_subsamples(self.meta, cat, control_cats)
        assert_array_equal(known_array, test_array)

    def test_get_paired_subsamples_not_strict(self):
        """Checks controlled subsets can be generated with missing values"""
        known_array = [array(['WS']), array(['NR'])]
        # Gets the test values
        cat = 'ABX'
        control_cats = ['AGE', 'RANGE']
        order = ['N', 'Y']
        test_array = get_paired_subsamples(self.meta, cat, control_cats,
                                           order, strict=False)
        assert_array_equal(known_array, test_array)

    def test_plot_effects_kwargs_error(self):
        """Tests that plot_effects errors when a bad keyword is passed"""
        self.assertRaises(ValueError, plot_effects, self.effects, self.bounds,
                          self.labels, arange(5, 50, 5), doctor=10)

    def test_plot_effects_2d_effects_error(self):
        """Tests plot_effects errors when effects_mean is a 2d array"""
        self.assertRaises(ValueError, plot_effects, ones((2, 2)), self.bounds,
                          self.labels, arange(5, 50, 5))

    def test_plot_effects_effects_and_bounds_or_labels_different_error(self):
        """Tests plot_effects erros when each mean doesnt have a bound"""
        # Checks the error is triggered when bounds is a different shape
        self.assertRaises(ValueError, plot_effects, self.effects, ones((2, 2)),
                          self.labels, arange(5, 50, 5))
        # Checks error is triggered when labels is a different shape
        self.assertRaises(ValueError, plot_effects, self.effects, self.bounds,
                          ones((2, 2)), arange(5, 50, 5))

    def test_plot_effects_size_returns(self):
        """Tests plot_effect_size returns a figure"""
        # Generates an effect figure
        fig = plot_effects(self.effects, self.bounds, self.labels,
                           arange(5, 50, 5), grid=False)
        # Test the output is a figure
        self.assertTrue(isinstance(fig, Figure))

    def test_collate_effect_size_counts_shape_error(self):
        """Checks collate_effect_size errors when counts is not a 1d array"""
        self.assertRaises(TypeError, collate_effect_size, self.powers[0],
                          self.powers, self.power_alpha)

    def test_collate_effect_size_different_power_shapes(self):
        """Checks collate_effect_size errors when power shapes are unequal"""
        powers = [ones((5)), ones((5, 2))]
        self.assertRaises(ValueError, collate_effect_size, self.counts,
                          powers, self.power_alpha)

    def test_collate_effect_size_power_and_counts_not_same_lenght(self):
        """Checks error is thrown when there is not a power for each count"""
        powers = [ones(4)]
        self.assertRaises(ValueError, collate_effect_size, self.counts,
                          powers, self.power_alpha)

    def test_collate_effect_size_list(self):
        """Checks collate_effect_size returns sanely for a list of powers"""
        # Calulates the test value
        test_m, test_b = collate_effect_size(self.counts,
                                             self.powers,
                                             self.power_alpha)
        # Checks the test values
        assert_almost_equal(self.effects, test_m, 4)
        assert_almost_equal(self.bounds, test_b, 4)

    def test_collate_effect_size_array(self):
        """Checks collate_effect_size returns sanely for a power array."""
        # Sets the know values
        known_m = self.effects[0]
        known_b = self.bounds[0]
        # Calulates the test value
        test_m, test_b = collate_effect_size(self.counts,
                                             self.powers[0],
                                             self.power_alpha)
        # Checks the test values
        assert_almost_equal(known_m, test_m, 4)
        assert_almost_equal(known_b, test_b, 4)

    def test_collate_effect_size_array_nans_2d(self):
        """Checks colate_effect_size handles an array with nans sanely"""
        self.powers[0][2, 1] = nan
        test_m, test_b = collate_effect_size(self.counts,
                                             self.powers,
                                             self.power_alpha)
        # Checks the test values
        assert_almost_equal(self.effects, test_m, 3)
        assert_almost_equal(self.bounds, test_b, 3)

    def test_collate_effect_size_under_power_1d(self):
        """Checks collate_effect_size handles undefined power sanely"""
        self.powers[0][0, 1] = 0.01
        powers = [power[0, :] for power in self.powers]
        known_m = array([0.1320, 0.3419, 0.5691])
        known_b = array([0.0308, 0.0076, 0.0334])

        test_m, test_b = collate_effect_size(self.counts,
                                             powers,
                                             self.power_alpha)
        # Checks the test values
        assert_almost_equal(known_m, test_m, 3)
        assert_almost_equal(known_b, test_b, 3)

    def test_collate_effect_size_nans_1d(self):
        """Checks collate_effect_size responds correctly to nans"""
        powers = [ones((5))*nan]
        test_m, test_b = collate_effect_size(self.counts,
                                             powers,
                                             self.power_alpha)
        self.assertTrue(isnan(test_m).all())
        self.assertTrue(isnan(test_b).all())

    def test_collate_effect_size_1d_power(self):
        """Tests collate_effect_size can handle power vectors"""
        powers = [power[0, :] for power in self.powers]
        known_m = array([0.1320, 0.3419, 0.5691])
        known_b = array([0.0217, 0.0076, 0.0334])
        test_m, test_b = collate_effect_size(self.counts,
                                             powers,
                                             self.power_alpha)
        # Checks the test values
        assert_almost_equal(known_m, test_m, 3)
        assert_almost_equal(known_b, test_b, 3)

if __name__ == '__main__':
    main()
