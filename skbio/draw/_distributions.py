#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils.six import string_types

from itertools import cycle
import warnings

import numpy as np

from matplotlib import use
use('Agg', warn=False)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon, Rectangle


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
    scipy.stats.ttest_ind

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
            list(map(float, distribution))
        except:
            raise ValueError("Each value in each distribution must be a "
                             "number.")

    _validate_x_values(x_values, x_tick_labels, len(distributions))

    # Create a new figure to plot our data on, and then plot the distributions.
    result, plot_axes = plt.subplots()
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


def grouped_distributions(plot_type, data, x_values=None,
                          data_point_labels=None, distribution_labels=None,
                          distribution_markers=None, x_label=None,
                          y_label=None, title=None,
                          x_tick_labels_orientation='vertical', y_min=None,
                          y_max=None, whisker_length=1.5,
                          error_bar_type='stdv', distribution_width=None,
                          figure_width=None, figure_height=None):
    """Generate a figure with distributions grouped at points along the x-axis.

    Parameters
    ----------
    plot_type : {'bar', 'scatter', 'box'}
        Type of plot to visualize distributions with.
    data : list of lists of lists
        Each inner list represents a data point along the x-axis. Each data
        point contains lists of data for each distribution in the group at that
        point. This nesting allows for the grouping of distributions at each
        data point.
    x_values : list of scalars, optional
        Spacing of data points along the x-axis. Must be the same length as the
        number of data points and be in ascending sorted order. If not
        provided, plots will be spaced evenly.
    data_point_labels : list of str, optional
        Labels for data points.
    distribution_labels : list of str, optional
        Labels for each distribution in a data point grouping.
    distribution_markers : list of str or list of tuple, optional
        Matplotlib-compatible strings or tuples that indicate the color or
        symbol to be used to distinguish each distribution in a data point
        grouping. Colors will be used for bar charts or box plots, while
        symbols will be used for scatter plots.
    x_label : str, optional
        x-axis label.
    y_label : str, optional
        y-axis label.
    title : str, optional
        Plot title.
    x_tick_labels_orientation : {'vertical', 'horizontal'}
        Orientation of x-axis labels.
    y_min : scalar, optional
        Minimum value of the y-axis. If ``None``, uses matplotlib's autoscale.
    y_max : scalar, optional
        Maximum value of the y-axis. If ``None``, uses matplotlib's autoscale.
    whisker_length : scalar, optional
        If `plot_type` is ``'box'``, determines the length of the whiskers as a
        function of the IQR. For example, if 1.5, the whiskers extend to
        ``1.5 * IQR``. Anything outside of that range is seen as an outlier.
        If `plot_type` is not ``'box'``, this parameter is ignored.
    error_bar_type : {'stdv', 'sem'}
        Type of error bars to use if `plot_type` is ``'bar'``. Can be either
        ``'stdv'`` (for standard deviation) or ``'sem'`` for the standard error
        of the mean. If `plot_type` is not ``'bar'``, this parameter is
        ignored.
    distribution_width : scalar, optional
        Width in plot units of each individual distribution (e.g. each bar if
        the plot type is a bar chart, or the width of each box if the plot type
        is a boxplot). If None, will be automatically determined.
    figure_width : scalar, optional
        Width of the plot figure in inches. If not provided, will default to
        matplotlib's default figure width.
    figure_height : scalar, optional
        Height of the plot figure in inches. If not provided, will default to
        matplotlib's default figure height.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing distributions grouped at points along the x-axis.

    Examples
    --------
    Create a plot with two distributions grouped at three points:

    .. plot::

       >>> from skbio.draw.distributions import grouped_distributions
       >>> fig = grouped_distributions('bar',
       ...                             [[[2, 2, 1,], [0, 1, 4]],
       ...                             [[1, 1, 1], [4, 4.5]],
       ...                             [[2.2, 2.4, 2.7, 1.0], [0, 0.2]]],
       ...                             distribution_labels=['Treatment 1',
       ...                                                  'Treatment 2'])

    """
    # Set up different behavior based on the plot type.
    if plot_type == 'bar':
        plotting_function = _plot_bar_data
        distribution_centered = False
        marker_type = 'colors'
    elif plot_type == 'scatter':
        plotting_function = _plot_scatter_data
        distribution_centered = True
        marker_type = 'symbols'
    elif plot_type == 'box':
        plotting_function = _plot_box_data
        distribution_centered = True
        marker_type = 'colors'
    else:
        raise ValueError("Invalid plot type '%s'. Supported plot types are "
                         "'bar', 'scatter', or 'box'." % plot_type)

    num_points, num_distributions = _validate_input(data, x_values,
                                                    data_point_labels,
                                                    distribution_labels)

    # Create a list of matplotlib markers (colors or symbols) that can be used
    # to distinguish each of the distributions. If the user provided a list of
    # markers, use it and loop around to the beginning if there aren't enough
    # markers. If they didn't provide a list, or it was empty, use our own
    # predefined list of markers (again, loop around to the beginning if we
    # need more markers).
    distribution_markers = _get_distribution_markers(marker_type,
                                                     distribution_markers,
                                                     num_distributions)

    # Now calculate where each of the data points will start on the x-axis.
    x_locations = _calc_data_point_locations(num_points, x_values)
    assert (len(x_locations) == num_points), "The number of x_locations " +\
        "does not match the number of data points."

    if distribution_width is None:
        # Find the smallest gap between consecutive data points and divide this
        # by the number of distributions + 1 for some extra spacing between
        # data points.
        min_gap = max(x_locations)
        for i in range(len(x_locations) - 1):
            curr_gap = x_locations[i + 1] - x_locations[i]
            if curr_gap < min_gap:
                min_gap = curr_gap

        distribution_width = min_gap / float(num_distributions + 1)
    else:
        if distribution_width <= 0:
            raise ValueError("The width of a distribution cannot be less than "
                             "or equal to zero.")

    result, plot_axes = plt.subplots()

    # Iterate over each data point, and plot each of the distributions at that
    # data point. Increase the offset after each distribution is plotted,
    # so that the grouped distributions don't overlap.
    for point, x_pos in zip(data, x_locations):
        dist_offset = 0
        for dist_index, dist, dist_marker in zip(range(num_distributions),
                                                 point, distribution_markers):
            dist_location = x_pos + dist_offset
            plotting_function(plot_axes, dist, dist_marker, distribution_width,
                              dist_location, whisker_length, error_bar_type)
            dist_offset += distribution_width

    # Set up various plot options that are best set after the plotting is done.
    # The x-axis tick marks (one per data point) are centered on each group of
    # distributions.
    plot_axes.set_xticks(_calc_data_point_ticks(x_locations,
                                                num_distributions,
                                                distribution_width,
                                                distribution_centered))
    _set_axes_options(plot_axes, title, x_label, y_label, x_values,
                      data_point_labels, x_tick_labels_orientation, y_min,
                      y_max)

    if distribution_labels is not None:
        _create_legend(plot_axes, distribution_markers, distribution_labels,
                       marker_type)

    _set_figure_size(result, figure_width, figure_height)

    # matplotlib seems to sometimes plot points on the rightmost edge of the
    # plot without adding padding, so we need to add our own to both sides of
    # the plot. For some reason this has to go after the call to draw(),
    # otherwise matplotlib throws an exception saying it doesn't have a
    # renderer. Boxplots need extra padding on the left.
    if plot_type == 'box':
        left_pad = 2 * distribution_width
    else:
        left_pad = distribution_width
    plot_axes.set_xlim(plot_axes.get_xlim()[0] - left_pad,
                       plot_axes.get_xlim()[1] + distribution_width)

    return result


def _validate_input(data, x_values, data_point_labels, distribution_labels):
    """Returns a tuple containing the number of data points and distributions
    in the data.

    Validates plotting options to make sure they are valid with the supplied
    data.
    """
    if data is None or not data or isinstance(data, string_types):
        raise ValueError("The data must be a list type, and it cannot be "
                         "None or empty.")

    num_points = len(data)
    num_distributions = len(data[0])

    empty_data_error_msg = ("The data must contain at least one data "
                            "point, and each data point must contain at "
                            "least one distribution to plot.")
    if num_points == 0 or num_distributions == 0:
        raise ValueError(empty_data_error_msg)

    for point in data:
        if len(point) == 0:
            raise ValueError(empty_data_error_msg)
        if len(point) != num_distributions:
            raise ValueError("The number of distributions in each data point "
                             "grouping must be the same for all data points.")

    # Make sure we have the right number of x values (one for each data point),
    # and make sure they are numbers.
    _validate_x_values(x_values, data_point_labels, num_points)

    if (distribution_labels is not None and
            len(distribution_labels) != num_distributions):
        raise ValueError("The number of distribution labels must be equal "
                         "to the number of distributions.")
    return num_points, num_distributions


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
            list(map(float, x_values))
        except:
            raise ValueError("Each x value must be a number.")

    if x_tick_labels is not None:
        if len(x_tick_labels) != num_expected_values:
            raise ValueError("The number of x-axis tick labels must match the "
                             "number of data points.")


def _get_distribution_markers(marker_type, marker_choices, num_markers):
    """Returns a list of length num_markers of valid matplotlib colors or
    symbols.

    The markers will be comprised of those found in marker_choices (if not None
    and not empty) or a list of predefined markers (determined by marker_type,
    which can be either 'colors' or 'symbols'). If there are not enough
    markers, the list of markers will be reused from the beginning again (as
    many times as are necessary).
    """
    if num_markers < 0:
        raise ValueError("num_markers must be greater than or equal to zero.")
    if marker_choices is None or len(marker_choices) == 0:
        if marker_type == 'colors':
            marker_choices = ['b', 'g', 'r', 'c', 'm', 'y', 'w']
        elif marker_type == 'symbols':
            marker_choices = \
                ['s', 'o', '^', '>', 'v', '<', 'd', 'p', 'h', '8', '+', 'x']
        else:
            raise ValueError("Invalid marker_type: '%s'. marker_type must be "
                             "either 'colors' or 'symbols'." % marker_type)
    if len(marker_choices) < num_markers:
        # We don't have enough markers to represent each distribution uniquely,
        # so let the user know. We'll add as many markers (starting from the
        # beginning of the list again) until we have enough, but the user
        # should still know because they may want to provide a new list of
        # markers.
        warnings.warn(
            "There are not enough markers to uniquely represent each "
            "distribution in your dataset. You may want to provide a list "
            "of markers that is at least as large as the number of "
            "distributions in your dataset.",
            RuntimeWarning)
        marker_cycle = cycle(marker_choices[:])
        while len(marker_choices) < num_markers:
            marker_choices.append(next(marker_cycle))
    return marker_choices[:num_markers]


def _calc_data_point_locations(num_points, x_values=None):
    """Returns the x-axis location for each of the data points to start at.

    Note: A numpy array is returned so that the overloaded "+" operator can be
    used on the array.

    The x-axis locations are scaled by x_values if it is provided, or else the
    x-axis locations are evenly spaced. In either case, the x-axis locations
    will always be in the range [1, num_points].
    """
    if x_values is None:
        # Evenly space the x-axis locations.
        x_locs = np.arange(1, num_points + 1)
    else:
        if len(x_values) != num_points:
            raise ValueError("The number of x-axis values must match the "
                             "number of data points.")

        # Scale to the range [1, num_points]. Taken from
        # http://www.heatonresearch.com/wiki/Range_Normalization
        x_min = min(x_values)
        x_max = max(x_values)
        x_range = x_max - x_min
        n_range = num_points - 1
        x_locs = np.array([(((x_val - x_min) * n_range) / float(x_range)) + 1
                           for x_val in x_values])

    return x_locs


def _calc_data_point_ticks(x_locations, num_distributions, distribution_width,
                           distribution_centered):
    """Returns a 1D numpy array of x-axis tick positions.

    These positions will be centered on each data point.

    Set distribution_centered to True for scatter and box plots because their
    plot types naturally center over a given horizontal position. Bar charts
    should use distribution_centered = False because the leftmost edge of a bar
    starts at a given horizontal position and extends to the right for the
    width of the bar.
    """
    dist_size = num_distributions - 1 if distribution_centered else\
        num_distributions
    return x_locations + ((dist_size * distribution_width) / 2)


def _plot_bar_data(plot_axes, distribution, distribution_color,
                   distribution_width, x_position, whisker_length,
                   error_bar_type):
    """Returns the result of plotting a single bar in matplotlib."""
    result = None

    # We do not want to plot empty distributions because matplotlib will not be
    # able to render them as PDFs.
    if len(distribution) > 0:
        avg = np.mean(distribution)
        if error_bar_type == 'stdv':
            error_bar = np.std(distribution)
        elif error_bar_type == 'sem':
            error_bar = np.std(distribution) / np.sqrt(len(distribution))
        else:
            raise ValueError(
                "Invalid error bar type '%s'. Supported error bar types are "
                "'stdv' and 'sem'." % error_bar_type)
        result = plot_axes.bar(x_position, avg, distribution_width,
                               yerr=error_bar, ecolor='black',
                               facecolor=distribution_color)
    return result


def _plot_scatter_data(plot_axes, distribution, distribution_symbol,
                       distribution_width, x_position, whisker_length,
                       error_bar_type):
    """Returns the result of plotting a single scatterplot in matplotlib."""
    result = None
    x_vals = [x_position] * len(distribution)

    # matplotlib's scatter function doesn't like plotting empty data.
    if len(x_vals) > 0 and len(distribution) > 0:
        result = plot_axes.scatter(x_vals, distribution,
                                   marker=distribution_symbol, c='k')
    return result


def _plot_box_data(plot_axes, distribution, distribution_color,
                   distribution_width, x_position, whisker_length,
                   error_bar_type):
    """Returns the result of plotting a single boxplot in matplotlib."""
    result = None

    if len(distribution) > 0:
        result = plot_axes.boxplot([distribution], positions=[x_position],
                                   widths=distribution_width,
                                   whis=whisker_length)
        _color_box_plot(plot_axes, result, [distribution_color])

    return result


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

            box_coords = list(zip(box_x, box_y))
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
        raise ValueError("Invalid orientation for x-axis tick labels: '%s'. "
                         "Valid orientations are 'vertical' or 'horizontal'."
                         % x_tick_labels_orientation)

    # If labels are provided, always use them. If they aren't, use the x_values
    # that denote the spacing between data points as labels. If that isn't
    # available, simply label the data points in an incremental fashion,
    # i.e. 1, 2, 3, ..., n, where n is the number of data points on the plot.
    if x_tick_labels is not None:
        plot_axes.set_xticklabels(x_tick_labels,
                                  rotation=x_tick_labels_orientation)
    elif x_tick_labels is None and x_values is not None:
        plot_axes.set_xticklabels(x_values, rotation=x_tick_labels_orientation)
    else:
        plot_axes.set_xticklabels(
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
        warnings.warn(
            "Could not automatically resize plot to make room for "
            "axes labels and plot title. This can happen if the labels or "
            "title are extremely long and the plot size is too small. Your "
            "plot may have its labels and/or title cut-off. To fix this, "
            "try increasing the plot's size (in inches) and try again.",
            RuntimeWarning)
