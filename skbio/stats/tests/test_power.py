#!/usr/bin/env python
# test_TaxPlot.py

from __future__ import division
from unittest import TestCase, main
from numpy import ones, power, array, arange, nan
import numpy.random
from numpy.testing import (assert_almost_equal,
                           assert_array_equal,
                           assert_allclose)
from pandas import DataFrame
from scipy.stats import kruskal
from skbio.stats.power import (_check_strs,
                               _confidence_bound,
                               _calculate_power,
                               compare_distributions,
                               calculate_power_curve,
                               bootstrap_power_curve,
                               get_signifigant_subsample,
                               get_paired_subsamples)


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
        meta = {'NR': {'RANGE': 'M', 'SEX': 'F', 'AGE': nan, 'ABX': 'Y'},
                'MH': {'RANGE': 'L', 'SEX': 'F', 'AGE': '30s', 'ABX': 'Y'},
                'PP': {'RANGE': 'M', 'SEX': 'F', 'AGE': '30s', 'ABX': 'N'},
                'CD': {'RANGE': 'L', 'SEX': 'F', 'AGE': '30s', 'ABX': 'Y'},
                'MM': {'RANGE': 'C', 'SEX': 'F', 'AGE': '30s', 'ABX': 'Y'},
                'SW': {'RANGE': 'M', 'SEX': 'M', 'AGE': nan, 'ABX': 'N'},
                'TS': {'RANGE': 'M', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'},
                'CB': {'RANGE': 'L', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'},
                'BB': {'RANGE': 'C', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'}}

        self.meta = DataFrame.from_dict(meta, orient='index')

    def test_get_paired_subsamples_default(self):
        """Checks controlled subsets can be generated sanely"""
        # Sets the known array set
        known_array = [array(['BB']), array(['CB']), array(['TS'])]
        # Gets the test value
        cat = 'RANGE'
        control_cats = ['SEX', 'AGE', 'ABX']
        test_array = get_paired_subsamples(self.meta, cat, control_cats)
        assert_array_equal(known_array, test_array)

    def test_get_paired_subsamples_not_strict(self):
        """Checks controlled subsets can be generated with missing values"""
        known_array = [array(['SW']), array(['NR'])]
        # Gets the test values
        cat = 'ABX'
        control_cats = ['AGE', 'RANGE']
        order = ['N', 'Y']
        test_array = get_paired_subsamples(self.meta, cat, control_cats,
                                           order, strict=False)
        assert_array_equal(known_array, test_array)

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
        self.assertRaises(TypeError, _check_strs, {1, 2, 3})

    def test__confidence_bound_default(self):
        """Checks confidence_bound correctly determines an interval"""
        # Sets the know confidence bound
        known = 2.1658506
        test = _confidence_bound(self.s1)
        assert_almost_equal(known, test, 5)

    def test__confidence_bound_df(self):
        """Checks a custom df for confidence_bound"""
        known = 2.0407076
        test = _confidence_bound(self.s1, df=15)
        assert_almost_equal(known, test, 5)

    def test__confidence_bound_alpha(self):
        """Checks a custom df for confidence_bound"""
        known = 3.111481
        test = _confidence_bound(self.s1, alpha=0.01)
        assert_almost_equal(known, test, 5)

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
        self.assertRaises(ValueError, get_signifigant_subsample, [None],
                          self.samps)

    def test_get_signifigant_subsample_no_results(self):
        """Checks get_signifigant_subsample errors when inputs are too similar
        """
        # Sets up a function which will all testing
        def test_f(x):
            if len(x[0]) == 100:
                return 0.001
            else:
                return 0.5
        # Tests a value error is raised
        self.assertRaises(ValueError, get_signifigant_subsample, [test_f],
                          self.samps, sub_size=10, num_rounds=5)

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

if __name__ == '__main__':
    main()
