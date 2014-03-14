#!/usr/bin/env python
"""
Distribution visualizations (:mod:`skbio.draw.distributions`)
=============================================================

.. currentmodule:: skbio.draw.distributions

This module provides plotting functionality for visualizing distributions.

Functions
---------

.. autosummary::
   :toctree: generated/

   boxplots

"""
from __future__ import division

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from matplotlib import use
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle
use('Agg', warn=False)
import matplotlib.pyplot as plt


def boxplots(distributions, x_values=None, x_tick_labels=None, title=None,
             x_label=None, y_label=None, x_tick_labels_orientation='vertical',
             y_min=None, y_max=None, whisker_length=1.5, box_width=0.5,
             box_colors=None, figure_width=None, figure_height=None,
             legend=None):
    """Generate a figure with a boxplot for each distribution.

    Parameters
    ----------
    distributions: list of lists
        List of distributions.
    x_values : list of numbers, optional
        List indicating where each boxplot should be placed. Must be the same
        length as `distributions` if provided.
    x_tick_labels : list of str, optional
        List of x-axis tick labels.
    title : str, optional
        Title of the plot.
    x_label : str, optional
        x-axis label.
    y_label : str, optional
        y-axis label.
    x_tick_labels_orientation : {'vertical', 'horizontal'}
        Orientation of the x-axis labels.
    y_min : scalar, optional
        Minimum value of the y-axis. If ``None``, uses matplotlib's autoscale.
    y_max : scalar, optional
        Maximum value of the y-axis. If ``None``, uses matplotlib's autoscale.
    whisker_length : scalar, optional
        Length of the whiskers as a function of the IQR. For example, if 1.5,
        the whiskers extend to ``1.5 * IQR``. Anything outside of that range is
        treated as an outlier.
    box_width : scalar, optional
        Width of each box in plot units.
    box_colors : str, tuple, or list of colors, optional
        Either a matplotlib-compatible string or tuple that indicates the color
        to be used for every boxplot, or a list of colors to color each boxplot
        individually. If ``None``, boxes will be the same color as the plot
        background. If a list of colors is provided, a color must be provided
        for each boxplot. Can also supply ``None`` instead of a color, which
        will color the box the same color as the plot background.
    figure_width : scalar, optional
        Width of the plot figure in inches. If not provided, will default to
        matplotlib's default figure width.
    figure_height : scalar, optional
        Height of the plot figure in inches. If not provided, will default to
        matplotlib's default figure height.
    legend : tuple or list, optional
        Two-element tuple or list that contains a list of valid matplotlib
        colors as the first element and a list of labels (strings) as the
        second element. The lengths of the first and second elements must be
        the same. If ``None``, a legend will not be plotted.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing a boxplot for each distribution.

    See Also
    --------
    matplotlib.pyplot.boxplot

    Notes
    -----
    This is a convenience wrapper around matplotlib's ``boxplot`` function that
    allows for coloring of boxplots and legend generation.

    Examples
    --------
    Create a plot with two boxplots:

    .. plot::

       >>> from skbio.draw.distributions import boxplots
       >>> fig = boxplots([[2, 2, 1, 3, 4, 4.2, 7], [0, -1, 4, 5, 6, 7]])

    Plot three distributions with custom colors and labels:

    .. plot::

       >>> from skbio.draw.distributions import boxplots
       >>> fig = boxplots(
       ...     [[2, 2, 1, 3], [0, -1, 0, 0.1, 0.3], [4, 5, 6, 3]],
       ...     x_tick_labels=('Control', 'Treatment 1', 'Treatment 2'),
       ...     box_colors=('green', 'blue', 'red'))

    """
    # Make sure our input makes sense.
    for distribution in distributions:
        try:
            map(float, distribution)
        except:
            raise ValueError("Each value in each distribution must be a "
                             "number.")

    _validate_x_values(x_values, x_tick_labels, len(distributions))

    # Create a new figure to plot our data on, and then plot the distributions.
    result, plot_axes = _create_plot()
    box_plot = plt.boxplot(distributions, positions=x_values,
                           whis=whisker_length, widths=box_width)

    if box_colors is not None:
        if _is_single_matplotlib_color(box_colors):
            box_colors = [box_colors] * len(box_plot['boxes'])
        else:
            # We check against the number of input distributions because mpl
            # will only return non-empty boxplots from the boxplot() call
            # above.
            if len(box_colors) != len(distributions):
                raise ValueError("Not enough colors were supplied to color "
                                 "each boxplot.")

            # Filter out colors corresponding to empty distributions.
            box_colors = [color for distribution, color in zip(distributions,
                                                               box_colors)
                          if distribution]

        _color_box_plot(plot_axes, box_plot, box_colors)

    # Set up the various plotting options, such as x- and y-axis labels, plot
    # title, and x-axis values if they have been supplied.
    _set_axes_options(plot_axes, title, x_label, y_label,
                      x_tick_labels=x_tick_labels,
                      x_tick_labels_orientation=x_tick_labels_orientation,
                      y_min=y_min, y_max=y_max)

    if legend is not None:
        if len(legend) != 2:
            raise ValueError("Invalid legend was provided. The legend must be "
                             "a two-element tuple/list where the first "
                             "element is a list of colors and the second "
                             "element is a list of labels.")
        _create_legend(plot_axes, legend[0], legend[1], 'colors')

    _set_figure_size(result, figure_width, figure_height)
    return result


def _validate_x_values(x_values, x_tick_labels, num_expected_values):
    """Validates the x values provided by the user, making sure they are the
    correct length and are all numbers.

    Also validates the number of x-axis tick labels.

    Raises a ValueError if these conditions are not met.
    """
    if x_values is not None:
        if len(x_values) != num_expected_values:
            raise ValueError("The number of x values must match the number "
                             "of data points.")
        try:
            map(float, x_values)
        except:
            raise ValueError("Each x value must be a number.")

    if x_tick_labels is not None:
        if len(x_tick_labels) != num_expected_values:
            raise ValueError("The number of x-axis tick labels must match the "
                             "number of data points.")


def _create_plot():
    """Creates a plot and returns the associated Figure and Axes objects."""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    return fig, ax


def _is_single_matplotlib_color(color):
    """Returns True if color is a single (not a list) mpl color."""
    single_color = False

    if (isinstance(color, str)):
        single_color = True
    elif len(color) == 3 or len(color) == 4:
        single_color = True

        for e in color:
            if not (isinstance(e, float) or isinstance(e, int)):
                single_color = False

    return single_color


def _color_box_plot(plot_axes, box_plot, colors):
    """Color boxes in the box plot with the specified colors.

    If any of the colors are None, the box will not be colored.

    The box_plot argument must be the dictionary returned by the call to
    matplotlib's boxplot function, and the colors argument must consist of
    valid matplotlib colors.
    """
    # Note: the following code is largely taken from this matplotlib boxplot
    # example:
    # http://matplotlib.sourceforge.net/examples/pylab_examples/
    #     boxplot_demo2.html
    if len(colors) != len(box_plot['boxes']):
        raise ValueError("Not enough colors were supplied to color each "
                         "boxplot.")

    for box, median, color in zip(box_plot['boxes'],
                                  box_plot['medians'],
                                  colors):
        if color is not None:
            box_x = []
            box_y = []

            # There are five points in the box. The first is the same as
            # the last.
            for i in range(5):
                box_x.append(box.get_xdata()[i])
                box_y.append(box.get_ydata()[i])

            box_coords = zip(box_x, box_y)
            box_polygon = Polygon(box_coords, facecolor=color)
            plot_axes.add_patch(box_polygon)

            # Draw the median lines back over what we just filled in with
            # color.
            median_x = []
            median_y = []
            for i in range(2):
                median_x.append(median.get_xdata()[i])
                median_y.append(median.get_ydata()[i])
                plot_axes.plot(median_x, median_y, 'black')


def _set_axes_options(plot_axes, title=None, x_label=None, y_label=None,
                      x_values=None, x_tick_labels=None,
                      x_tick_labels_orientation='vertical', y_min=None,
                      y_max=None):
    """Applies various labelling options to the plot axes."""
    if title is not None:
        plot_axes.set_title(title)
    if x_label is not None:
        plot_axes.set_xlabel(x_label)
    if y_label is not None:
        plot_axes.set_ylabel(y_label)

    if (x_tick_labels_orientation != 'vertical' and
            x_tick_labels_orientation != 'horizontal'):
        raise ValueError("Invalid orientation for x-axis tick labels: %s. "
                         "Valid orientations are 'vertical' or 'horizontal'."
                         % x_tick_labels_rotation)

    # If labels are provided, always use them. If they aren't, use the x_values
    # that denote the spacing between data points as labels. If that isn't
    # available, simply label the data points in an incremental fashion,
    # i.e. 1, 2, 3, ..., n, where n is the number of data points on the plot.
    if x_tick_labels is not None:
        labels = plot_axes.set_xticklabels(x_tick_labels,
                                           rotation=x_tick_labels_orientation)
    elif x_tick_labels is None and x_values is not None:
        labels = plot_axes.set_xticklabels(x_values,
                                           rotation=x_tick_labels_orientation)
    else:
        labels = plot_axes.set_xticklabels(
            range(1, len(plot_axes.get_xticklabels()) + 1),
            rotation=x_tick_labels_orientation)

    # Set the y-axis range if specified.
    if y_min is not None:
        plot_axes.set_ylim(bottom=float(y_min))
    if y_max is not None:
        plot_axes.set_ylim(top=float(y_max))


def _create_legend(plot_axes, distribution_markers, distribution_labels,
                   marker_type):
    """Creates a legend on the supplied axes."""
    # We have to use a proxy artist for the legend because box plots currently
    # don't have a very useful legend in matplotlib, and using the default
    # legend for bar/scatterplots chokes on empty/null distributions.
    #
    # Note: This code is based on the following examples:
    #   http://matplotlib.sourceforge.net/users/legend_guide.html
    #   http://stackoverflow.com/a/11423554
    if len(distribution_markers) != len(distribution_labels):
        raise ValueError("The number of distribution markers does not match "
                         "the number of distribution labels.")
    if marker_type == 'colors':
        legend_proxy = [Rectangle((0, 0), 1, 1, fc=marker)
                        for marker in distribution_markers]
        plot_axes.legend(legend_proxy, distribution_labels, loc='best')
    elif marker_type == 'symbols':
        legend_proxy = [Line2D(range(1), range(1), color='white',
                        markerfacecolor='black', marker=marker)
                        for marker in distribution_markers]
        plot_axes.legend(legend_proxy, distribution_labels, numpoints=3,
                         scatterpoints=3, loc='best')
    else:
        raise ValueError("Invalid marker_type: '%s'. marker_type must be "
                         "either 'colors' or 'symbols'." % marker_type)


def _set_figure_size(fig, width=None, height=None):
    """Sets the plot figure size and makes room for axis labels, titles, etc.

    If both width and height are not provided, will use matplotlib defaults.

    Making room for labels will not always work, and if it fails, the user will
    be warned that their plot may have cut-off labels.
    """
    # Set the size of the plot figure, then make room for the labels so they
    # don't get cut off. Must be done in this order.
    if width is not None and height is not None and width > 0 and height > 0:
        fig.set_size_inches(width, height)
    try:
        fig.tight_layout()
    except ValueError:
        print ("Warning: could not automatically resize plot to make room for "
               "axes labels and plot title. This can happen if the labels or "
               "title are extremely long and the plot size is too small. Your "
               "plot may have its labels and/or title cut-off. To fix this, "
               "try increasing the plot's size (in inches) and try again.")
