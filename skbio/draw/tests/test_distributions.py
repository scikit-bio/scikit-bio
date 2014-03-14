#!/usr/bin/env python
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import sys
from StringIO import StringIO
from unittest import TestCase, main

import numpy as np
import matplotlib.pyplot as plt

from skbio.draw.distributions import (boxplots, _validate_x_values,
                                      _create_plot,
                                      _is_single_matplotlib_color,
                                      _color_box_plot, _set_axes_options,
                                      _create_legend, _set_figure_size)


class DistributionsTests(TestCase):

    def setUp(self):
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

    def test_validate_x_values_invalid_x_values(self):
        """_validate_x_values() should raise a ValueError on an invalid number
        of x_values."""
        self.assertRaises(ValueError, _validate_x_values,
                          [1, 2, 3, 4], ["T0", "T1", "T2"],
                          len(self.ValidSingleSampleData))

    def test_validate_x_values_invalid_x_tick_labels(self):
        """_validate_x_values() should raise a ValueError on an invalid number
        of x_tick_labels."""
        self.assertRaises(ValueError, _validate_x_values,
                          None, ["T0"], len(self.ValidSingleSampleData))

    def test_validate_x_values_nonnumber_x_values(self):
        """_validate_x_values() should raise a ValueError on x_values that
        aren't numbers."""
        self.assertRaises(ValueError, _validate_x_values,
                          ["foo", 2, 3], None, len(self.ValidSingleSampleData))

    def test_validate_x_values_valid_x_values(self):
        """_validate_x_values() should not throw an exception."""
        _validate_x_values([1, 2.0, 3], None, 3)

    def test_create_plot(self):
        """_create_plot() should return a tuple containing a Figure and
        Axes."""
        fig, ax = _create_plot()
        self.assertEqual(fig.__class__.__name__, "Figure")
        self.assertEqual(ax.__class__.__name__, "AxesSubplot")

    def test_set_axes_options(self):
        """_set_axes_options() should set the labels on the axes and not raise
        any exceptions."""
        fig, ax = _create_plot()
        _set_axes_options(ax, "Plot Title", "x-axis label", "y-axis label",
                          x_tick_labels=["T0", "T1"])
        self.assertEqual(ax.get_title(), "Plot Title")
        self.assertEqual(ax.get_ylabel(), "y-axis label")
        self.assertEqual(ax.get_xticklabels()[0].get_text(), "T0")
        self.assertEqual(ax.get_xticklabels()[1].get_text(), "T1")

    def test_set_axes_options_ylim(self):
        """_set_axes_options() should set the y-axis limits."""
        fig, ax = _create_plot()
        _set_axes_options(ax, "Plot Title", "x-axis label", "y-axis label",
                          x_tick_labels=["T0", "T1", "T2"], y_min=0, y_max=1)
        self.assertEqual(ax.get_title(), "Plot Title")
        self.assertEqual(ax.get_ylabel(), "y-axis label")
        self.assertEqual(ax.get_xticklabels()[0].get_text(), "T0")
        self.assertEqual(ax.get_xticklabels()[1].get_text(), "T1")
        self.assertEqual(ax.get_ylim(), (0.0, 1.0))

    def test_set_axes_options_bad_ylim(self):
        """_set_axes_options() should raise an exception when given non-numeric
        y limits."""
        fig, ax = _create_plot()
        self.assertRaises(ValueError, _set_axes_options, ax, "Plot Title",
                          "x-axis label", "y-axis label",
                          x_tick_labels=["T0", "T1", "T2"], y_min='car',
                          y_max=30)

    def test_create_legend(self):
        """_create_box_plot_legend() should create a legend on valid input."""
        fig, ax = _create_plot()
        _create_legend(ax, ['b', 'r'], ['dist1', 'dist2'], 'colors')
        self.assertEqual(len(ax.get_legend().get_texts()), 2)

        fig, ax = _create_plot()
        _create_legend(ax, ['^', '<', '>'], ['dist1', 'dist2', 'dist3'],
                       'symbols')
        self.assertEqual(len(ax.get_legend().get_texts()), 3)

    def test_create_legend_invalid_input(self):
        """Test raises error on bad input."""
        fig, ax = _create_plot()
        self.assertRaises(ValueError, _create_legend, ax,
                          ['^', '<', '>'], ['dist1', 'dist2'], 'symbols')
        self.assertRaises(ValueError, _create_legend, ax, ['^', '<', '>'],
                          ['dist1', 'dist2', 'dist3'], 'foo')

    def test_boxplots(self):
        """boxplots() should return a valid Figure object."""
        fig = boxplots(self.ValidTypicalBoxData, [1, 4, 10],
                       ["Data 1", "Data 2", "Data 3"], "Test", "x-axis label",
                       "y-axis label")
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
        self.assertRaises(ValueError, boxplots, [[1, 'foo', 3]])

        # Number of colors doesn't match number of distributions.
        self.assertRaises(ValueError, boxplots, [[1, 2, 3], [],
                          [4, 5, 6]], box_colors=['blue', 'red'])

        # Invalid legend.
        self.assertRaises(ValueError, boxplots, [[1, 2, 3]],
                          legend=('foo', 'bar', 'baz'))

    def test_color_box_plot(self):
        """Should not throw an exception when passed the proper input."""
        fig, ax = _create_plot()
        box_plot = plt.boxplot(self.ValidTypicalBoxData)
        _color_box_plot(ax, box_plot, ['blue', 'w', (1, 1, 0.9)])

        # Some colors are None.
        fig, ax = _create_plot()
        box_plot = plt.boxplot(self.ValidTypicalBoxData)
        _color_box_plot(ax, box_plot, ['blue', None, (1, 1, 0.9)])

        # All colors are None.
        fig, ax = _create_plot()
        box_plot = plt.boxplot(self.ValidTypicalBoxData)
        _color_box_plot(ax, box_plot, [None, None, None])

    def test_color_box_plot_invalid_input(self):
        """Should throw an exception on invalid input."""
        # Invalid color.
        fig, ax = _create_plot()
        box_plot = plt.boxplot(self.ValidTypicalBoxData)
        self.assertRaises(ValueError, _color_box_plot, ax, box_plot,
                          ['red', 'foobarbaz', 'blue'])

        # Wrong number of colors.
        fig, ax = _create_plot()
        box_plot = plt.boxplot(self.ValidTypicalBoxData)
        self.assertRaises(ValueError, _color_box_plot, ax, box_plot,
                          ['blue', (1, 1, 0.9)])

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
        fig, ax = _create_plot()
        _set_axes_options(ax, 'foo', 'x_foo', 'y_foo',
                          x_tick_labels=['foofoofoo', 'barbarbar'],
                          x_tick_labels_orientation='vertical')
        _set_figure_size(fig, 3, 4)
        self.assertTrue(np.array_equal(fig.get_size_inches(), (3, 4)))

    def test_set_figure_size_defaults(self):
        """Test setting a figure size using matplotlib defaults."""
        fig, ax = _create_plot()
        _set_axes_options(ax, 'foo', 'x_foo', 'y_foo',
                          x_tick_labels=['foofoofoo', 'barbarbar'],
                          x_tick_labels_orientation='vertical')
        orig_fig_size = fig.get_size_inches()
        _set_figure_size(fig)
        self.assertTrue(np.array_equal(fig.get_size_inches(), orig_fig_size))

    def test_set_figure_size_invalid(self):
        """Test setting a figure size using invalid dimensions."""
        fig, ax = _create_plot()
        _set_axes_options(ax, 'foo', 'x_foo', 'y_foo',
                          x_tick_labels=['foofoofoo', 'barbarbar'],
                          x_tick_labels_orientation='vertical')
        orig_fig_size = fig.get_size_inches()
        _set_figure_size(fig, -1, 0)
        self.assertTrue(np.array_equal(fig.get_size_inches(), orig_fig_size))

    def test_set_figure_size_long_labels(self):
        """Test setting a figure size that has really long labels."""
        saved_stdout = sys.stdout
        try:
            out = StringIO()
            sys.stdout = out

            fig, ax = _create_plot()
            _set_axes_options(ax, 'foo', 'x_foo', 'y_foo',
                              x_tick_labels=['foofoofooooooooooooooooooooooooo'
                                             'oooooooooooooooooooooooooooooooo'
                                             'oooooooooooooooooooooooooooooooo'
                                             'oooo', 'barbarbar'],
                              x_tick_labels_orientation='vertical')
            _set_figure_size(fig, 3, 3)
            self.assertTrue(np.array_equal(fig.get_size_inches(),
                            (3, 3)))
            output = out.getvalue().strip()
            self.assertEqual(
                output,
                "Warning: could not automatically resize plot to make room "
                "for axes labels and plot title. This can happen if the "
                "labels or title are extremely long and the plot size is too "
                "small. Your plot may have its labels and/or title cut-off. "
                "To fix this, try increasing the plot's size (in inches) and "
                "try again.")
        finally:
            sys.stdout = saved_stdout


if __name__ == '__main__':
    main()
