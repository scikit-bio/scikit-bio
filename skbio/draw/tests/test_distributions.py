#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

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
from matplotlib import use
use('Agg', warn=False)
import matplotlib.pyplot as plt

from skbio.draw.distributions import (boxplots, grouped_distributions,
                                      _calc_data_point_locations,
                                      _calc_data_point_ticks, _color_box_plot,
                                      _create_legend,
                                      _get_distribution_markers,
                                      _is_single_matplotlib_color,
                                      _plot_bar_data, _plot_box_data,
                                      _plot_scatter_data, _set_axes_options,
                                      _set_figure_size, _validate_input,
                                      _validate_x_values)


class DistributionsTests(TestCase):

    def setUp(self):
        # Test null data list.
        self.Null = None

        # Test empty data list.
        self.Empty = []

        # Test nested empty data list.
        self.EmptyNested = [[]]

        # Test nested empty data list (for bar/scatter plots).
        self.EmptyDeeplyNested = [[[]]]

        # Test invalid number of samples in data list (for bar/scatter plots).
        self.InvalidNumSamples = [[[1, 2, 3, 4, 5]],
                                  [[4, 5, 6, 7, 8], [2, 3, 2]],
                                  [[4, 7, 10, 33, 32, 6, 7, 8]]]

        # Test valid data with three samples and four data points
        # (for bar/scatter plots).
        self.ValidTypicalData = [[[1.0, 2, 3.5, 5], [2, 3, 5, 6], [2, 3, 8]],
                                 [[4, 7, 8], [8, 9, 10, 11], [9.0, 4, 1, 1]],
                                 [[4, 33, 32, 6, 8], [5, 4, 8, 13], [1, 1, 2]],
                                 [[2, 2, 2, 2], [3, 9, 8], [2, 1, 6, 7, 4, 5]]]

        # Test valid data with one sample (for bar/scatter plots).
        self.ValidSingleSampleData = [[[1, 2, 3, 4, 5]],
                                      [[4, 5, 6, 7, 8]],
                                      [[4, 7, 10, 33, 32, 6, 7, 8]]]

        # Test typical data to be plotted by the boxplot function.
        self.ValidTypicalBoxData = [[3.4, 10, 11.67, 12.0, 2, 2, 99.99],
                                    [2.3, 4, 5, 88, 9, 10, 11, 1, 0, 3, -8],
                                    [2, 9, 7, 5, 6]]

    def tearDown(self):
        # We get a warning from mpl if we don't clean up our figures.
        plt.close('all')

    def test_validate_input_null(self):
        """_validate_input() should raise a ValueError if null data is passed
        to it."""
        with npt.assert_raises(ValueError):
            _validate_input(self.Null, None, None, None)

    def test_validate_input_empty(self):
        """_validate_input() should raise a ValueError if empty data is passed
        to it."""
        with npt.assert_raises(ValueError):
            _validate_input(self.Empty, None, None, None)

    def test_validate_input_empty_nested(self):
        """_validate_input() should raise a ValueError if empty nested data is
        passed to it."""
        with npt.assert_raises(ValueError):
            _validate_input(self.EmptyNested, None, None, None)

    def test_validate_input_empty_deeply_nested(self):
        """_validate_input() should pass for deeply nested empty data."""
        num_points, num_samples = _validate_input(self.EmptyDeeplyNested,
                                                  None, None, None)
        self.assertEqual(num_points, 1)
        self.assertEqual(num_samples, 1)

    def test_validate_input_empty_point(self):
        with npt.assert_raises(ValueError):
            _validate_input([[[1, 2, 3], [4, 5]], []], None, None, None)

    def test_validate_input_invalid_num_samples(self):
        """_validate_input() should raise a ValueError if an inconsistent
        number of samples in included in the data."""
        with npt.assert_raises(ValueError):
            _validate_input(self.InvalidNumSamples, None, None, None)

    def test_validate_input_invalid_data_point_names(self):
        """_validate_input() should raise a ValueError on data_point_names that
        are an invalid length."""
        with npt.assert_raises(ValueError):
            _validate_input(self.ValidSingleSampleData, None, ["T0", "T1"],
                            None)

    def test_validate_input_invalid_sample_names(self):
        """_validate_input() should raise a ValueError on sample_names that are
        an invalid length."""
        with npt.assert_raises(ValueError):
            _validate_input(self.ValidSingleSampleData, None, None,
                            ["Men", "Women"])

    def test_validate_input_all_valid_input(self):
        """_validate_input() should return valid information about the data
        without throwing an exception."""
        self.assertEqual(_validate_input(self.ValidTypicalData, [1, 3, 4, 8],
                                         ["T0", "T1", "T2", "T3"],
                                         ["Infants", "Children", "Teens"]),
                         (4, 3))

    def test_validate_x_values_invalid_x_values(self):
        """_validate_x_values() should raise a ValueError on an invalid number
        of x_values."""
        with npt.assert_raises(ValueError):
            _validate_x_values([1, 2, 3, 4], ["T0", "T1", "T2"],
                               len(self.ValidSingleSampleData))

    def test_validate_x_values_invalid_x_tick_labels(self):
        """_validate_x_values() should raise a ValueError on an invalid number
        of x_tick_labels."""
        with npt.assert_raises(ValueError):
            _validate_x_values(None, ["T0"], len(self.ValidSingleSampleData))

    def test_validate_x_values_nonnumber_x_values(self):
        """_validate_x_values() should raise a ValueError on x_values that
        aren't numbers."""
        with npt.assert_raises(ValueError):
            _validate_x_values(["foo", 2, 3], None,
                               len(self.ValidSingleSampleData))

    def test_validate_x_values_valid_x_values(self):
        """_validate_x_values() should not throw an exception."""
        _validate_x_values([1, 2.0, 3], None, 3)

    def test_get_distribution_markers_null_marker_list(self):
        """_get_distribution_markers() should return a list of predefined
        matplotlib markers."""
        self.assertEqual(_get_distribution_markers('colors', None, 5),
                         ['b', 'g', 'r', 'c', 'm'])

    def test_get_distribution_markers_empty_marker_list(self):
        """_get_distribution_markers() should return a list of predefined
        matplotlib markers."""
        self.assertEqual(_get_distribution_markers('colors', None, 4),
                         ['b', 'g', 'r', 'c'])

    def test_get_distribution_markers_insufficient_markers(self):
        """_get_distribution_markers() should return a wrapped list of
        predefined markers."""
        self.assertEqual(npt.assert_warns(RuntimeWarning,
                                          _get_distribution_markers,
                                          'colors', None, 10),
                         ['b', 'g', 'r', 'c', 'm', 'y', 'w', 'b', 'g', 'r'])
        self.assertEqual(npt.assert_warns(RuntimeWarning,
                                          _get_distribution_markers,
                                          'symbols', ['^', '>', '<'], 5),
                         ['^', '>', '<', '^', '>'])

    def test_get_distribution_markers_bad_marker_type(self):
        """_get_distribution_markers() should raise a ValueError."""
        with npt.assert_raises(ValueError):
            _get_distribution_markers('shapes', [], 3)

    def test_get_distribution_markers_zero_markers(self):
        """_get_distribution_markers() should return an empty list."""
        self.assertEqual(_get_distribution_markers('symbols', None, 0), [])
        self.assertEqual(_get_distribution_markers('symbols', ['^'], 0), [])

    def test_get_distribution_markers_negative_num_markers(self):
        with npt.assert_raises(ValueError):
            _get_distribution_markers('symbols', [], -1)

    def test_plot_bar_data(self):
        """_plot_bar_data() should return a list of Rectangle objects."""
        fig, ax = plt.subplots()
        result = _plot_bar_data(ax, [1, 2, 3], 'red', 0.5, 3.75, 1.5, 'stdv')
        self.assertEqual(result[0].__class__.__name__, "Rectangle")
        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(result[0].get_width(), 0.5)
        self.assertAlmostEqual(result[0].get_facecolor(), (1.0, 0.0, 0.0, 1.0))
        self.assertAlmostEqual(result[0].get_height(), 2.0)

        fig, ax = plt.subplots()
        result = _plot_bar_data(ax, [1, 2, 3], 'red', 0.5, 3.75, 1.5, 'sem')
        self.assertEqual(result[0].__class__.__name__, "Rectangle")
        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(result[0].get_width(), 0.5)
        self.assertAlmostEqual(result[0].get_facecolor(), (1.0, 0.0, 0.0, 1.0))
        self.assertAlmostEqual(result[0].get_height(), 2.0)

    def test_plot_bar_data_bad_error_bar_type(self):
        """_plot_bar_data() should raise an exception on bad error bar type."""
        fig, ax = plt.subplots()
        with npt.assert_raises(ValueError):
            _plot_bar_data(ax, [1, 2, 3], 'red', 0.5, 3.75, 1.5, 'var')

    def test_plot_bar_data_empty(self):
        """_plot_bar_data() should not error when given empty list of data,
        but should not plot anything."""
        fig, ax = plt.subplots()
        result = _plot_bar_data(ax, [], 'red', 0.5, 3.75, 1.5, 'stdv')
        self.assertTrue(result is None)

        fig, ax = plt.subplots()
        result = _plot_bar_data(ax, [], 'red', 0.5, 3.75, 1.5, 'sem')
        self.assertTrue(result is None)

    def test_plot_scatter_data(self):
        """_plot_scatter_data() should return a Collection instance."""
        fig, ax = plt.subplots()
        result = _plot_scatter_data(ax, [1, 2, 3], '^', 0.77, 1, 1.5, 'stdv')
        self.assertEqual(result.get_sizes(), 20)

    def test_plot_scatter_data_empty(self):
        """_plot_scatter_data() should not error when given empty list of data,
        but should not plot anything."""
        fig, ax = plt.subplots()
        result = _plot_scatter_data(ax, [], '^', 0.77, 1, 1.5, 'stdv')
        self.assertTrue(result is None)

    def test_plot_box_data(self):
        """_plot_box_data() should return a dictionary for Line2D's."""
        fig, ax = plt.subplots()
        result = _plot_box_data(ax, [0, 0, 7, 8, -3, 44], 'blue', 0.33, 55,
                                1.5, 'stdv')
        self.assertEqual(result.__class__.__name__, "dict")
        self.assertEqual(len(result['boxes']), 1)
        self.assertEqual(len(result['medians']), 1)
        self.assertEqual(len(result['whiskers']), 2)
        self.assertEqual(len(result['fliers']), 2)
        self.assertEqual(len(result['caps']), 2)

    def test_plot_box_data_empty(self):
        """Should ignore empty distribution."""
        fig, ax = plt.subplots()
        result = _plot_box_data(ax, [], 'blue', 0.33, 55, 1.5, 'stdv')
        self.assertTrue(result is None)

    def test_calc_data_point_locations_invalid_x_values(self):
        """Should raise error when invalid x_values are encountered."""
        with npt.assert_raises(ValueError):
            _calc_data_point_locations(3, [1, 10.5])

    def test_calc_data_point_locations_default_spacing(self):
        """Should return evenly-spaced x-axis locations."""
        locs = _calc_data_point_locations(4)
        np.testing.assert_allclose(locs, [1, 2, 3, 4])

    def test_calc_data_point_locations_custom_spacing(self):
        """Should return non-evenly-spaced x-axis locations."""
        # Scaling down from 3..12 to 1..4.
        locs = _calc_data_point_locations(4, [3, 4, 10, 12])
        np.testing.assert_allclose(locs,
                                   np.array([1, 1.33333333, 3.33333333, 4]))

        # Sorted order shouldn't affect scaling.
        locs = _calc_data_point_locations(4, [4, 3, 12, 10])
        np.testing.assert_allclose(locs,
                                   np.array([1.33333333, 1, 4, 3.33333333]))

        # Scaling up from 0.001..0.87 to 1..3.
        locs = _calc_data_point_locations(3, [0.001, 0.2543, 0.87])
        np.testing.assert_allclose(locs,
                                   np.array([1, 1.58296893, 3]))

    def test_calc_data_point_ticks(self):
        """_calc_data_point_ticks() should return an array containing the
        x-axis locations for each data point tick."""
        ticks = _calc_data_point_ticks(np.array([1, 5, 9, 11]), 1, 0.5, False)
        np.testing.assert_allclose(ticks, [1.25, 5.25, 9.25, 11.25])

        ticks = _calc_data_point_ticks(np.array([0]), 3, 0.5, False)
        np.testing.assert_allclose(ticks, [0.75])

    def test_set_axes_options(self):
        """_set_axes_options() should set the labels on the axes and not raise
        any exceptions."""
        fig, ax = plt.subplots()
        _set_axes_options(ax, "Plot Title", "x-axis label", "y-axis label",
                          x_tick_labels=["T0", "T1"])
        self.assertEqual(ax.get_title(), "Plot Title")
        self.assertEqual(ax.get_ylabel(), "y-axis label")
        self.assertEqual(ax.get_xticklabels()[0].get_text(), "T0")
        self.assertEqual(ax.get_xticklabels()[1].get_text(), "T1")

    def test_set_axes_options_ylim(self):
        """_set_axes_options() should set the y-axis limits."""
        fig, ax = plt.subplots()
        _set_axes_options(ax, "Plot Title", "x-axis label", "y-axis label",
                          x_tick_labels=["T0", "T1", "T2"], y_min=0, y_max=1)
        self.assertEqual(ax.get_title(), "Plot Title")
        self.assertEqual(ax.get_ylabel(), "y-axis label")
        self.assertEqual(ax.get_xticklabels()[0].get_text(), "T0")
        self.assertEqual(ax.get_xticklabels()[1].get_text(), "T1")
        self.assertEqual(ax.get_ylim(), (0.0, 1.0))

    def test_set_axes_options_x_values_as_tick_labels(self):
        fig, ax = plt.subplots()
        _set_axes_options(ax, "Plot Title", "x-axis label", "y-axis label",
                          x_values=[42, 45, 800])

        self.assertEqual(ax.get_title(), "Plot Title")
        self.assertEqual(ax.get_ylabel(), "y-axis label")
        self.assertEqual(ax.get_xticklabels()[0].get_text(), '42')
        self.assertEqual(ax.get_xticklabels()[1].get_text(), '45')
        self.assertEqual(ax.get_xticklabels()[2].get_text(), '800')

    def test_set_axes_options_bad_ylim(self):
        """_set_axes_options() should raise an exception when given non-numeric
        y limits."""
        fig, ax = plt.subplots()
        with npt.assert_raises(ValueError):
            _set_axes_options(ax, "Plot Title", "x-axis label", "y-axis label",
                              x_tick_labels=["T0", "T1", "T2"], y_min='car',
                              y_max=30)

    def test_set_axes_options_invalid_x_tick_labels_orientation(self):
        fig, ax = plt.subplots()
        with npt.assert_raises(ValueError):
            _set_axes_options(ax, "Plot Title", "x-axis label", "y-axis label",
                              x_tick_labels=["T0", "T1"],
                              x_tick_labels_orientation='brofist')

    def test_create_legend(self):
        """_create_box_plot_legend() should create a legend on valid input."""
        fig, ax = plt.subplots()
        _create_legend(ax, ['b', 'r'], ['dist1', 'dist2'], 'colors')
        self.assertEqual(len(ax.get_legend().get_texts()), 2)

        fig, ax = plt.subplots()
        _create_legend(ax, ['^', '<', '>'], ['dist1', 'dist2', 'dist3'],
                       'symbols')
        self.assertEqual(len(ax.get_legend().get_texts()), 3)

    def test_create_legend_invalid_input(self):
        """Test raises error on bad input."""
        fig, ax = plt.subplots()
        with npt.assert_raises(ValueError):
            _create_legend(ax, ['^', '<', '>'], ['dist1', 'dist2'], 'symbols')
        with npt.assert_raises(ValueError):
            _create_legend(ax, ['^', '<', '>'], ['dist1', 'dist2', 'dist3'],
                           'foo')

    def test_grouped_distributions_bar(self):
        """Should return a valid barchart Figure object."""
        fig = grouped_distributions('bar', self.ValidTypicalData,
                                    [1, 4, 10, 11], ["T0", "T1", "T2", "T3"],
                                    ["Infants", "Children", "Teens"],
                                    ['b', 'r', 'g'], "x-axis label",
                                    "y-axis label", "Test")
        ax = fig.get_axes()[0]
        self.assertEqual(ax.get_title(), "Test")
        self.assertEqual(ax.get_xlabel(), "x-axis label")
        self.assertEqual(ax.get_ylabel(), "y-axis label")
        self.assertEqual(len(ax.get_xticklabels()), 4)
        np.testing.assert_allclose(ax.get_xticks(),
                                   [1.1125, 2.0125, 3.8125, 4.1125])

    def test_grouped_distributions_insufficient_colors(self):
        """grouped_distributions() should work even when there aren't
        enough colors. We should capture a warning."""
        args = ('bar', self.ValidTypicalData, [1, 4, 10, 11],
                ["T0", "T1", "T2", "T3"], ["Infants", "Children", "Teens"],
                ['b', 'r'], "x-axis label", "y-axis label", "Test")

        npt.assert_warns(RuntimeWarning,
                         grouped_distributions,
                         *args)

    def test_grouped_distributions_scatter(self):
        """Should return a valid scatterplot Figure object."""
        fig = grouped_distributions('scatter', self.ValidTypicalData,
                                    [1, 4, 10, 11], ["T0", "T1", "T2", "T3"],
                                    ["Infants", "Children", "Teens"],
                                    ['^', '>', '<'], "x-axis label",
                                    "y-axis label", "Test")
        ax = fig.get_axes()[0]
        self.assertEqual(ax.get_title(), "Test")
        self.assertEqual(ax.get_xlabel(), "x-axis label")
        self.assertEqual(ax.get_ylabel(), "y-axis label")
        self.assertEqual(len(ax.get_xticklabels()), 4)
        np.testing.assert_allclose(ax.get_xticks(),
                                   [1.075, 1.975, 3.775, 4.075])

    def test_grouped_distributions_insufficient_symbols(self):
        """grouped_distributions() should work even when there aren't
        enough symbols. We should capture a warning."""
        args = ('scatter', self.ValidTypicalData, [1, 4, 10, 11],
                ["T0", "T1", "T2", "T3"], ["Infants", "Children", "Teens"],
                ['^'], "x-axis label", "y-axis label", "Test")

        npt.assert_warns(RuntimeWarning, grouped_distributions, *args)

    def test_grouped_distributions_empty_marker_list(self):
        """grouped_distributions() should use the predefined list of
        markers if an empty list is provided by the user."""
        grouped_distributions('scatter', self.ValidTypicalData,
                              [1, 4, 10, 11], ["T0", "T1", "T2", "T3"],
                              ["Infants", "Children", "Teens"], [],
                              "x-axis label", "y-axis label", "Test")

    def test_grouped_distributions_box(self):
        """Should return a valid boxplot Figure object."""
        fig = grouped_distributions('box', self.ValidTypicalData,
                                    [1, 4, 10, 11], ["T0", "T1", "T2", "T3"],
                                    ["Infants", "Children", "Teens"],
                                    ['b', 'g', 'y'], "x-axis label",
                                    "y-axis label", "Test")
        ax = fig.get_axes()[0]
        self.assertEqual(ax.get_title(), "Test")
        self.assertEqual(ax.get_xlabel(), "x-axis label")
        self.assertEqual(ax.get_ylabel(), "y-axis label")
        self.assertEqual(len(ax.get_xticklabels()), 4)
        np.testing.assert_allclose(ax.get_xticks(),
                                   [1.075, 1.975, 3.775, 4.075])

    def test_grouped_distributions_error(self):
        """grouped_distributions() should raise a ValueError for an
        invalid plot type."""
        with npt.assert_raises(ValueError):
            grouped_distributions('pie', self.ValidTypicalData,
                                  [1, 4, 10, 11], ["T0", "T1", "T2", "T3"],
                                  ["Infants", "Children", "Teens"],
                                  ['b', 'g', 'y'],
                                  "x-axis label", "y-axis label", "Test")

    def test_grouped_distributions_negative_distribution_width(self):
        args = ('box', self.ValidTypicalData, [1, 4, 10, 11],
                ["T0", "T1", "T2", "T3"], ["Infants", "Children", "Teens"],
                ['b', 'g', 'y'], "x-axis label", "y-axis label", "Test")

        with self.assertRaises(ValueError):
            grouped_distributions(*args, distribution_width=0)

        with self.assertRaises(ValueError):
            grouped_distributions(*args, distribution_width=-42)

    def test_boxplots(self):
        """boxplots() should return a valid Figure object."""
        fig = boxplots(self.ValidTypicalBoxData, [1, 4, 10],
                       ["Data 1", "Data 2", "Data 3"], "Test", "x-axis label",
                       "y-axis label",
                       legend=(('blue', 'red'), ('foo', 'bar')))
        ax = fig.get_axes()[0]
        self.assertEqual(ax.get_title(), "Test")
        self.assertEqual(ax.get_xlabel(), "x-axis label")
        self.assertEqual(ax.get_ylabel(), "y-axis label")
        self.assertEqual(len(ax.get_xticklabels()), 3)
        self.assertTrue(np.array_equal(ax.get_xticks(), [1, 4, 10]))

    def test_boxplots_empty_distributions(self):
        """Test functions correctly with empty distributions."""
        fig = boxplots([[1, 2, 3], [], [4, 5, 6]], [1, 4, 10],
                       ["Data 1", "Data 2", "Data 3"], "Test", "x-axis label",
                       "y-axis label")
        ax = fig.get_axes()[0]
        self.assertEqual(ax.get_title(), "Test")
        self.assertEqual(ax.get_xlabel(), "x-axis label")
        self.assertEqual(ax.get_ylabel(), "y-axis label")
        self.assertEqual(len(ax.get_xticklabels()), 3)
        self.assertTrue(np.array_equal(ax.get_xticks(), [1, 4, 10]))

        # All distributions are empty.
        fig = boxplots([[], [], []], [1, 4, 10],
                       ["Data 1", "Data 2", "Data 3"], "Test", "x-axis label",
                       "y-axis label")
        ax = fig.get_axes()[0]
        self.assertEqual(ax.get_title(), "Test")
        self.assertEqual(ax.get_xlabel(), "x-axis label")
        self.assertEqual(ax.get_ylabel(), "y-axis label")
        self.assertEqual(len(ax.get_xticklabels()), 3)
        self.assertTrue(np.array_equal(ax.get_xticks(), [1, 4, 10]))

    def test_boxplots_box_colors(self):
        """Test correctly handles coloring of box plots."""
        # Coloring works with all empty distributions.
        fig = boxplots([[], [], []], box_colors=['blue', 'red', 'yellow'])
        ax = fig.get_axes()[0]
        self.assertEqual(len(ax.get_xticklabels()), 3)

        fig = boxplots([[], [], []], box_colors='pink')
        ax = fig.get_axes()[0]
        self.assertEqual(len(ax.get_xticklabels()), 3)

        # Coloring works with some empty distributions.
        fig = boxplots([[], [1, 2, 3.5], []],
                       box_colors=['blue', 'red', 'yellow'])
        ax = fig.get_axes()[0]
        self.assertEqual(len(ax.get_xticklabels()), 3)

    def test_boxplots_invalid_input(self):
        """Test correctly throws error on invalid input."""
        # Non-numeric entries in distribution.
        with npt.assert_raises(ValueError):
            boxplots([[1, 'foo', 3]])

        # Number of colors doesn't match number of distributions.
        with npt.assert_raises(ValueError):
            boxplots([[1, 2, 3], [], [4, 5, 6]], box_colors=['blue', 'red'])

        # Invalid legend.
        with npt.assert_raises(ValueError):
            boxplots([[1, 2, 3]], legend=('foo', 'bar', 'baz'))

    def test_color_box_plot(self):
        """Should not throw an exception when passed the proper input."""
        fig, ax = plt.subplots()
        box_plot = plt.boxplot(self.ValidTypicalBoxData)
        _color_box_plot(ax, box_plot, ['blue', 'w', (1, 1, 0.9)])

        # Some colors are None.
        fig, ax = plt.subplots()
        box_plot = plt.boxplot(self.ValidTypicalBoxData)
        _color_box_plot(ax, box_plot, ['blue', None, (1, 1, 0.9)])

        # All colors are None.
        fig, ax = plt.subplots()
        box_plot = plt.boxplot(self.ValidTypicalBoxData)
        _color_box_plot(ax, box_plot, [None, None, None])

    def test_color_box_plot_invalid_input(self):
        """Should throw an exception on invalid input."""
        # Invalid color.
        fig, ax = plt.subplots()
        box_plot = plt.boxplot(self.ValidTypicalBoxData)
        with npt.assert_raises(ValueError):
            _color_box_plot(ax, box_plot, ['red', 'foobarbaz', 'blue'])

        # Wrong number of colors.
        fig, ax = plt.subplots()
        box_plot = plt.boxplot(self.ValidTypicalBoxData)
        with npt.assert_raises(ValueError):
            _color_box_plot(ax, box_plot, ['blue', (1, 1, 0.9)])

    def test_is_single_matplotlib_color(self):
        """Test correct identification of single versus multiple mpl colors."""
        self.assertTrue(_is_single_matplotlib_color('w'))
        self.assertTrue(_is_single_matplotlib_color('white'))
        self.assertTrue(_is_single_matplotlib_color([1, 1, 1]))
        self.assertTrue(_is_single_matplotlib_color([1, 1, 1, 1]))
        self.assertTrue(_is_single_matplotlib_color((1, 1, 1)))
        self.assertTrue(_is_single_matplotlib_color((1, 1, 1, 1)))
        self.assertTrue(_is_single_matplotlib_color((1.0, 1.0, 1.0, 1.0)))
        self.assertTrue(_is_single_matplotlib_color((1.0, 1, 1.0)))
        self.assertTrue(_is_single_matplotlib_color((2.0, 1, 1.0)))

        self.assertFalse(_is_single_matplotlib_color(['w', 'r']))
        self.assertFalse(_is_single_matplotlib_color(['w']))
        self.assertFalse(_is_single_matplotlib_color(('w',)))
        self.assertFalse(_is_single_matplotlib_color(((1.0, 1.0, 1),)))
        self.assertFalse(_is_single_matplotlib_color(((1.0, 1.0, 1),
                                                      (0.9, 0.9))))

    def test_set_figure_size(self):
        """Test setting a valid figure size."""
        fig, ax = plt.subplots()
        _set_axes_options(ax, 'foo', 'x_foo', 'y_foo',
                          x_tick_labels=['foofoofoo', 'barbarbar'],
                          x_tick_labels_orientation='vertical')
        _set_figure_size(fig, 3, 4)
        self.assertTrue(np.array_equal(fig.get_size_inches(), (3, 4)))

    def test_set_figure_size_defaults(self):
        """Test setting a figure size using matplotlib defaults."""
        fig, ax = plt.subplots()
        _set_axes_options(ax, 'foo', 'x_foo', 'y_foo',
                          x_tick_labels=['foofoofoo', 'barbarbar'],
                          x_tick_labels_orientation='vertical')
        orig_fig_size = fig.get_size_inches()
        _set_figure_size(fig)
        self.assertTrue(np.array_equal(fig.get_size_inches(), orig_fig_size))

    def test_set_figure_size_invalid(self):
        """Test setting a figure size using invalid dimensions."""
        fig, ax = plt.subplots()
        _set_axes_options(ax, 'foo', 'x_foo', 'y_foo',
                          x_tick_labels=['foofoofoo', 'barbarbar'],
                          x_tick_labels_orientation='vertical')
        orig_fig_size = fig.get_size_inches()
        _set_figure_size(fig, -1, 0)
        self.assertTrue(np.array_equal(fig.get_size_inches(), orig_fig_size))

    def test_set_figure_size_long_labels(self):
        """Test setting a figure size that has really long labels."""
        fig, ax = plt.subplots()
        _set_axes_options(ax, 'foo', 'x_foo', 'y_foo',
                          x_tick_labels=['foofoofooooooooooooooooooooooooo'
                                         'oooooooooooooooooooooooooooooooo'
                                         'oooooooooooooooooooooooooooooooo'
                                         'oooo', 'barbarbar'],
                          x_tick_labels_orientation='vertical')
        npt.assert_warns(RuntimeWarning, _set_figure_size, fig, 3, 3)
        npt.assert_array_equal(fig.get_size_inches(), (3, 3))


if __name__ == '__main__':
    main()
