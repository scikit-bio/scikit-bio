r"""Empirical Power Estimation (:mod: `skbio.stats.power`)
=========================================================

.. currentmodule:: skbio.stats.power

The purpose of this module is to provide empirical, post-hoc power estimation
of microbiome data. It also provides support to subsample data to faciliate
this analysis.

The underlying principle is based on subsampling and monte carlo simulation.
Assume that there is some set of populations, $K_{1}, K_{2}, ... K_{n}$ which
have some property, $\mu$ such that $\mu_{1} \neq \mu_{2} \neq ... \neq
\mu_{n}$. For each of the populations, a sample, $S$ can be drawn, with a
parameter, $x$ where $x \aeq \mu$ and for the samples, we can use a test, f,
to show that $x_{1} \neq x_{2} \neq ... \neq x_{n}$.

Since we known that $\mu_{1} \neq \mu_{2} \neq ... \neq \mu_{n}$, we know we
should reject the null hypothesis. If we fail to reject the null hypothesis,
we have comitted a Type II error and our result is a False negative. We can
estimate the frequency of Type II errors at various sampling depths by
repeatedly subsampling the populations and observing how often we see a False
negative. If we repeat this several times for each subsampling depth, and vary
the depths we use, we can start to approximate a relationship between the
number of samples we use and the rate of false negatives, also called the
statistical power of the test.

We can then use the rate of false negatives and use the `statsmodels.power`
package to solve for an effect size. This can be used to extrapolate a power
curve for the data.

The general format for functions in this module is to define a statistical
test function which will take a list of ids, and return a p value. The test is
then evaluated over a series of subsample sizes.

With microbiome data, there are three ways we can approach selecting our
sample. We may choose to simply draw $n$ observations at random from the two
underlying samples. Alternatively, we can draw subsamples which are
significantly different. Finally, we can try to match samples based on a set
of control categories.

Example
-------
 Suppose we have 100 samples randomly drawn from two normal distribitions,
the first with mean 0 and standard devation 1, and the second with mean 3
and standard deviation 1.5

>>> import numpy as np
>>> samples_1 = np.random.randn(100)
>>> samples_2 = 1.5*np.random.randn(100) + 1

We want to test the statistical power of a kruskal-wallis test comparing
the two populations. We can define a test function, f, to perform the
comparison. The test function will take a list of value vectors and
return a p value.

>>> from scipy.stats import ttest_ind
>>> f = lambda x: ttest_ind(x[0], x[1])[1]

Now, we can determine the statitical power, or the probability that do not
have a false positive given that we do not have a false negative by varying
a number of subsamples.

>>> sample_counts = np.arange(5, 80, 5)
>>> power_mean, power_bound = bootstrap_power_curve(f,
...                                                 [samples_1, samples_2],
...                                                 sample_counts)
>>> print power_mean
    [ 0.2772  0.569   0.7744  0.9052  0.969   0.9898  0.9984  0.9998  1.
      1.      1.      1.      1.      1.      1.    ]
>>> print power_bound
    [ 0.0178  0.0124  0.0145  0.0097  0.0053  0.0027  0.0013  0.0004  0.
      0.      0.      0.      0.      0.      0.    ]

"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------

from __future__ import division
from future.utils import viewitems
from copy import deepcopy
from numpy import (array, zeros, ones, round as nround, hstack, isnan, ndarray,
                   nan, sqrt, arange, delete, where)
from numpy.random import choice
from scipy.stats import t, nanmean, nanstd
from statsmodels.stats.power import FTestAnovaPower
from matplotlib import rcParams

# Gets a power solving instance
ft = FTestAnovaPower()

# Sets up plotting parameters so that the default setting is use to Helvetica
# in plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Arial']
rcParams['text.usetex'] = True


def make_power_curves(mode, tests, cats, samples=None, meta=None, **kwargs):
    r"""Plots a power curve for a specified set of data and categories

    Parameters
    ----------
    mode : {"SIGNIFICANT", "ALL", "PAIRED"}
        how random observations should be drawn.
        "SIGNIFICANT" indicates that observations should be drawn from two
        subsets of samples which are significantly different at the level
        indictated by `alpha_pwr` / `scaling`.
        "ALL" indicates the observations should be drawn from the sample, and
        not subsamples.
        "PAIRED" indicates samples should be drawn using metadata to control
        for variable categories. The samples drawn are limited by the `cats`
        being varied and the `control_cats` for the samples given in `order`.
        Samples will differ in the `cat`, and matched in `control_cats` so that
        for a set of observations, all the values in the `control_cats` will
        be the same, and there will be an observation representing each
        sample specified in `order`.
    test : function
        the statistical test which accepts an array-like of sample ids
        (list of lists) and returns a p-value.
    cats : array
        a list of the category names used for calculating the effect sizes.
    samples : {None, array}, optional
        Default is None. `samples` can be a list of lists or an array where
        each sublist or row in the array corresponds to a sampled group.
        `samples` must be defined if `mode` is "ALL" or "SIGNIFICANT".
    meta : {None, dataframe}, optional
        Default is None.the metadata associated with the samples. `meta` must
        be defined if `mode` is "PAIRED".
    control_cats : {array-like, None}, optional
        Default is None. `control_cats` cannot be None if `mode` is "PAIRED".
        The metadata categories to be used as controls. For example, if you
        wanted to control age (`cat` = "AGE"), you might want to control for
        gender and health status (i.e. `control_cats` = ["SEX", "HEALTHY"]).
        `control_cats` can intersect with `cats`. If this is the case, the
        experimental cat will be removed from the control_categories during
        that analysis.
    order : {None, list}, optional
        Default is None. The order of groups in the category. This can be used
        to limit the groups selected. For example, if there's a category with
        groups 'A', 'B' and 'C', and you only want to look at A vs B, `order`
        would be set to ['A', 'B'].
    sub_size : {None, int}, optional
        the maximum number of samples to select from a group. If no value is
        provided, this will be the same as the size of the smallest group.
        Otherwise, this will be compared to the size of the smallest group, and
        which ever is lower will be used.
    counts : array, optional
        the sample counts to be plotted in the power curve.

    Returns
    -------
    eff_fig : figure
        figure showing the effect_size curves
    eff_means : array
        the mean effect size for each category
    eff_bounds : array
        the bound for each category, where the confidence interval around the
        effect size is defined as [`eff_means` - `eff_bounds`, `eff_means`
        + `eff_bounds`]
    labels : array
        the cleaned up category names

    Raises
    ------
    ValueError
        if mode is not "PAIRED", "SIGNIFICANT", or "ALL"
    ValueError
        if mode is "PAIRED" and meta or control_cats are not defined
    ValueError
        if mode is "ALL" or "SIGNIFICANT" and samples is not defined

    Also See
    --------
    get_paired_effect
    get_unpaired_effect
    plot_effects

    """

    # Handles keyword arguments
    eff_kwds = {'control_cats': None,
                'order': None,
                'sub_size': None,
                'alpha_pwr': 0.05,
                'min_counts': 20,
                'max_counts': 50,
                'count_int': 10,
                'num_iter': 500,
                'num_runs': 10,
                'strict': True,
                'scaling': 5,
                'labels': None,
                'sort_plot': True,
                'counts': hstack(((2, 5), arange(10, 255, 10)))}

    plot_kwds = {'alpha': 0.05,
                 'colormap': None,
                 'grid': True,
                 'title': '',
                 'show_bound': True,
                 'leg_offset': None,
                 'tick_size': 12,
                 'label_size': 15,
                 'title_size': 18,
                 'legend_size': 11}

    # Loops through keyword arguments
    for key, value in viewitems(kwargs):
        if key in eff_kwds:
            eff_kwds[key] = value
        elif key in plot_kwds:
            plot_kwds[key] = value
        else:
            raise ValueError('%s is not a supported keyword argument.' % key)

    # Checks the mode argument is sane
    if mode not in {'PAIRED', 'SIGNIFICANT', 'ALL'}:
        raise ValueError('%s is not a supported mode.' % mode)
    elif mode == 'PAIRED' and (meta is None or
                               eff_kwds['control_cats'] is None):
        raise ValueError('Metadata and control_cats must be supplied for '
                         'PAIRED data.')
    elif mode in {'SIGNIFICANT', 'ALL'} and samples is None:
        raise ValueError('samples must be supplied for SIGNIFICANT or ALL '
                         'mode.')

    # Prealocates a holding objects
    counts = []
    powers = []

    num_cats = len(cats)

    # Gets the power matrix
    for idx, test in enumerate(tests):
        # Handles paired data
        if mode == 'PAIRED':
            # Makes sure the category is not in the control cats
            cat = cats[idx]
            ctrl_cats = deepcopy(eff_kwds['control_cats'])
            if cat in ctrl_cats:
                control_cats = delete(ctrl_cats, where(ctrl_cats == cat))
            # Calcualtes the power
            try:
                pwr, cnt = get_paired_effect(test, meta, cat, control_cats,
                                             eff_kwds['order'],
                                             eff_kwds['alpha_pwr'],
                                             eff_kwds['min_counts'],
                                             eff_kwds['max_counts'],
                                             eff_kwds['count_int'],
                                             eff_kwds['num_iter'],
                                             eff_kwds['num_runs'],
                                             eff_kwds['strict'])
            except:
                pwr = array([nan])
                cnt = array([nan])
        # Handles unpaired data
        else:
            pwr, cnt = get_unpaired_effect(mode, test, samples,
                                           eff_kwds['sub_size'],
                                           eff_kwds['alpha_pwr'],
                                           eff_kwds['min_counts'],
                                           eff_kwds['max_counts'],
                                           eff_kwds['count_int'],
                                           eff_kwds['num_iter'],
                                           eff_kwds['num_runs'],
                                           eff_kwds['scaling'])

        # Adds the power and count data to the holding object
        powers.append(pwr)
        counts.append(cnt)

    # Gets the effect sizes
    eff_means, eff_bounds = collate_effect_size(counts, powers,
                                                eff_kwds['alpha_pwr'])

    # Checks the label
    if eff_kwds['labels'] is not None:
        if not len(eff_kwds['labels']) == num_cats:
            raise ValueError('Each category must have a label.')
        labels = array(eff_kwds['labels'])
    else:
        labels = array([cat_.replace('_', ' ').title() for cat_ in cats])

    # Sorts the data
    removes = isnan(eff_means)
    if eff_kwds['sort_plot']:
        eff_means = eff_means[True - removes]
        eff_order = eff_means.argsort()[::-1]
        eff_means = eff_means[eff_order]
        eff_bounds = eff_bounds[eff_order]
        labels = labels[eff_order]
        if plot_kwds['colormap'] is not None:
            if len(plot_kwds['colormap'].shape) == 1:
                plot_kwds['colormap'] = plot_kwds['colormap'][eff_order]
            else:
                plot_kwds['colormap'] = plot_kwds['colormap'][eff_order, :]

    # Plots the effect sizes
    if not isnan(eff_means).all():
        eff_fig = plot_effects(eff_means, eff_bounds, labels,
                               eff_kwds['counts'], **plot_kwds)
    else:
        eff_fig = None

    return eff_fig, eff_means, eff_bounds, labels


def get_paired_effect(test, meta, cat, control_cats, order=None,
                      alpha_pwr=0.05, min_counts=20, max_counts=50,
                      counts_interval=10, num_iter=500, num_runs=10,
                      strict=True):
    """Calculates the effect size using paired subsampling

    Parameters
    ----------
    test : function
        the statistical test which accepts an array-like of sample ids
        (list of lists) and returns a p-value.

    meta : dataframe
        the metadata associated with the samples

    cat : str
        the metadata categories for comparison

    control_cats : list
        the metadata categories to be used as controls. For example, if you
        wanted to control age (`cat` = "AGE"), you might want to control for
        gender and health status (i.e. `control_cats` = ["SEX", "HEALTHY"])

    order : {None, list}, optional
        Default is None. The order of groups in the category. This can be used
        to limit the groups selected. For example, if there's a category with
        groups 'A', 'B' and 'C', and you only want to look at A vs B, `order`
        would be set to ['A', 'B'].

    alpha_pwr : float, optional
        Default is 0.05. The alpha value used to calculate the power.

    min_counts : unsigned int, optional
        Default is 20. The minimum number of paired samples which must exist
        for a category and set of control categories to be able to subsample
        and make power calculations.

    max_counts : unsigned int, optional
        Default is 50. The maximum number of samples per group to draw for
        effect size calculation.

    counts_interval : unsigned int, optional
        Default is 10.

    num_iter : unsigned int
        Default is 1000. The number of p-values to generate for each point
        on the curve.

    num_runs : unsigned int
        Default is 10. The number of times to calculate each curve.

    strict: bool, optional
        Default is True. If a missing value (nan) is encountered, the group
        will be skipped when 'strict' is True.

    Returns
    -------
    power : array
        power calculated for each subsample at each count

    sample_counts : array
        the number of samples drawn at each power calculation

    Raises
    ------
    RuntimeError
        if the paired samples contains less than the minimum numebr of samples.

    """

    # Gets a paired sample population to check the number of pairs generated
    paired_ids = get_paired_subsamples(meta, cat, control_cats, order, strict)
    num_paired = paired_ids[0].shape[0]

    # Checks there are enough paired ids to subsample
    if num_paired <= min_counts:
        raise RuntimeError('There are not enough samples for subsampling.')

    # Gets the sampling array
    sample_counts = arange(counts_interval,
                           min(max_counts, num_paired),
                           counts_interval)

    # Prealocates a power array
    power = zeros((num_runs, len(sample_counts)))

    # Calculates power or hte first curve
    power[0, :] = calculate_power_curve(test, paired_ids, sample_counts,
                                        num_iter=num_iter, alpha=alpha_pwr)
    # Gets iteraitons and subsequent power
    for id1 in arange(1, num_runs):
        paired_ids = get_paired_subsamples(meta, cat, control_cats, order,
                                           strict)
        power[id1, :] = calculate_power_curve(test, paired_ids, sample_counts,
                                              num_iter=num_iter,
                                              alpha=alpha_pwr)
    return power, sample_counts


def get_unpaired_effect(mode, test, samples, sub_size=None, alpha_pwr=0.05,
                        min_counts=20, max_counts=50, counts_interval=10,
                        num_iter=500, num_runs=10, scaling=5):
    """Calculates the effect size for unpaired random sampling methods

    Parameters
    ----------
    mode : {"SIGNIFICANT", "ALL"}
        how random observations should be drawn. "SIGNIFICANT" indicates that
        observations should be drawn from two subsets of samples which are
        significantly different at the level indictated by `alpha_pwr` /
        `scaling`, while "ALL" indicates the observations should be drawn from
        the sample, and not subsamples.

    test : function
        the statistical test which accepts an array-like of sample ids
        (list of lists) and returns a p-value.

    samples : array-like
        samples can be a list of lists or an array where each sublist or row in
        the array corresponds to a sampled group.

    sub_size : {None, int}, optional
        the maximum number of samples to select from a group. If no value is
        provided, this will be the same as the size of the smallest group.
        Otherwise, this will be compared to the size of the smallest group, and
        which ever is lower will be used.

    alpha_pwr : float, optional
        default is 0.05. The critical value for the power calculation.

    min_counts : unsigned int, optional
        Default is 20. The minimum number of paired samples which must exist
        for a category and set of control categories to be able to subsample
        and make power calculations.

    max_counts : unsigned int, optional
        Default is 50. The maximum number of samples per group to draw for
        effect size calculation.

    counts_interval : unsigned int, optional
        Default is 10.

    num_iter : unsigned int, optional
        Default is 1000. The number of p-values to generate for each point
        on the curve.

    num_runs : unsigned int, optional
        Default is 10. The number of times to calculate each curve.

    scaling : int, optional
        a penalty scale on `alpha_pwr`, so the probability that two
        distributions are different in "SIGNIFICANT" mode is less than
        `alpha_pwr` / `scaling`.

    labels : 1d array
        a list of formatted strings describing the effects, to be used in the
        legend.

    counts : 1d array
        the counts where power should be calculated.

    Returns
    -------
    power : array
        power calculated for each subsample at each count

    sample_counts : array
        the number of samples drawn at each power calculation

    Raises
    ------
    RuntimeError
        if the paired samples contains less than the minimum numebr of samples.

    """
    # Gets a population of sample ids to check the number of subsamples
    # generated
    if mode == 'SIGNIFICANT':
            sub_ids = get_significant_subsample([test], samples, sub_size,
                                                alpha_pwr, num_iter, scaling)
    else:
        sub_ids = samples
    num_ids = len(sub_ids[0])

    # Checks there are enough samples to subsample
    if num_ids <= min_counts:
        raise RuntimeError('There are not enough samples for subsampling.')

    # Calculates the effect size vector
    sample_counts = arange(counts_interval,
                           min(max_counts, num_ids),
                           counts_interval)

    # Prealocates the power array
    power = zeros((num_runs, len(sample_counts)))

    # Calculates the first power curve instance
    power[0, :] = calculate_power_curve(test, sub_ids, sample_counts,
                                        num_iter=num_iter, alpha=alpha_pwr)

    # Calculates the power instances
    for id1 in arange(1, num_runs):
        # Gets the subsample
        if mode == 'SIGNIFICANT':
            sub_ids = get_significant_subsample([test], samples, sub_size,
                                                alpha_pwr, num_iter, scaling)
        else:
            sub_ids = samples
        # Calculates the power curve
        power[id1, :] = calculate_power_curve(test, sub_ids, sample_counts,
                                              num_iter=num_iter,
                                              alpha=alpha_pwr)

    return power, sample_counts


def plot_effects(effect_means, effect_bounds, labels, sample_counts, **kwargs):
    """Makes a power curve plot

    Parameters
    ----------
    effect_means: 1d array
        the mean effect sizes to plots.

    effect_bounds : {None, 1d array}
        the range used for the confidence interval. If there is no effect to
        show, this should be None.

    labels : 1d array
        a list of formatted strings describing the effects, to be used in the
        legend.

    sample_counts : 1d array
        the counts where power should be calculated.

    alpha : int, optional
        Default is 0.05. The critical value for the power curves.

    colormap : {None, array}, optional
        Default is None. A colormap to use for the lines. Each color
        designation must appear in a new row. If no colormap is supplied, the
        defualt colormap will be used.

    grid : bool, optional
        Default is True. Show grid.

    show_bound : bool
        Default is True. Shows the confidence bounds on the effect size. If
        `effect_bounds` is None, no bounds will be shown.

    Returns
    -------
    fig : figure
        a figure with the power curves plotted.

    Other parameters
    ----------------
    leg_offset : tuple
        Changes the legend position.

    tick_size : usigned int
        sets the font size for tick labels

    label_size : unsigned int
        sets the font size for the axis labels

    title_size : unsigned int
        sets the font size for the title

    legend_size : unsigned int
        sets the font size for enteries in the legend

    """
    # Sets the keyword properties
    kwds = {'alpha': 0.05,
            'colormap': None,
            'grid': True,
            'title': '',
            'show_bound': True,
            'leg_offset': None,
            'tick_size': 12,
            'label_size': 15,
            'title_size': 18,
            'legend_size': 11}
    for key, value in viewitems(kwargs):
        if key in kwds:
            kwds[key] = value
        else:
            raise ValueError('%s is not a property of plot_effects.' % key)

    # Checks the effect, bound, and mean argument is sane
    mean_shape = effect_means.shape
    if effect_bounds is None:
        kwds['show_bound'] = False
        effect_bounds = zeros(mean_shape)
    bound_shape = effect_bounds.shape
    label_shape = labels.shape

    if not len(mean_shape) == 1:
        raise ValueError('Effect Mean must be a 1d numpy array')
    elif mean_shape != bound_shape or mean_shape != label_shape:
        raise ValueError('There must be a label and bound for each effect.')

    # Plots the the lower bound data
    fig = ft.plot_power(dep_var='nobs',
                        nobs=sample_counts,
                        effect_size=effect_means - effect_bounds,
                        alpha=kwds['alpha'])
    # Gets the axis of the first plot and its position
    lax = fig.axes[0]
    # Makes the lower bound lines dashed and thin, and changes the color if
    # desired
    for idx, l in enumerate(lax.get_lines()):
        l.set_linestyle(':')
        l.set_linewidth(1.5)
        if kwds['colormap'] is not None and len(kwds['colormap'].shape) == 1:
            l.set_color(kwds['colormap'][idx])
        elif kwds['colormap'] is not None:
            l.set_color(kwds['colormap'][idx, :])
    # Hides the x ticks and labels
    lax.set_title('')
    lax.set_xticklabels('')
    lax.set_yticklabels('')
    lax.set_xlabel('')
    # Hides the legend
    lax.get_legend().set_visible(False)

    # Plots the upper bound data
    uax = fig.add_axes(lax.get_position())
    fig = ft.plot_power('nobs', sample_counts, effect_means + effect_bounds,
                        alpha=kwds['alpha'], ax=uax)
    # Makes the lower bound axes visable, if desired
    if kwds['show_bound']:
        uax.set_axis_bgcolor('none')
    # Makes the lower bound lines dashed and thin, and changes the color if
    # desired
    for idx, l in enumerate(uax.get_lines()):
        l.set_linestyle(':')
        l.set_linewidth(1.5)
        if kwds['colormap'] is not None and len(kwds['colormap'].shape) == 1:
            l.set_color(kwds['colormap'][idx])
        elif kwds['colormap'] is not None:
            l.set_color(kwds['colormap'][idx, :])
    # Hides the x ticks and labels
    uax.set_title('')
    uax.set_xticklabels('')
    uax.set_yticklabels('')
    uax.set_xlabel('')
    # Hides the legend
    uax.get_legend().set_visible(False)

    # Plots the mean data
    axm = fig.add_axes(lax.get_position())
    fig = ft.plot_power('nobs', sample_counts, effect_means, ax=axm,
                        alpha=kwds['alpha'])

    # Shows the confidence bounds, if desired
    if kwds['show_bound']:
        axm.set_axis_bgcolor('none')

    # Recolors the lines, if desired
    if kwds['colormap'] is not None and len(kwds['colormap'].shape) == 1:
        for idx, l in enumerate(axm.get_lines()):
            l.set_color(kwds['colormap'][idx])
    elif kwds['colormap'] is not None:
        for idx, l in enumerate(axm.get_lines()):
            l.set_color(kwds['colormap'][idx, :])

    # Sets up the labels
    axm.set_xticklabels(map(int, axm.get_xticks()), size=kwds['tick_size'])
    axm.set_yticklabels(axm.get_yticks(), size=kwds['tick_size'])
    axm.set_xlabel('Number of Observations', size=kwds['label_size'])
    axm.set_ylabel('Power of the Test', size=kwds['label_size'])
    axm.set_title(kwds['title'], size=kwds['title_size'])

    # Adds the grid, if desired
    if kwds['grid']:
        axm.grid()

    leg = axm.get_legend()
    # Sets the legend position
    if kwds['leg_offset'] is not None:
        leg.set_bbox_to_anchor(kwds['leg_offset'])
    # Sets up the legend text
    for idx, txt in enumerate(leg.get_texts()):
        txt.set_text(labels[idx])
        txt.set_size(kwds['legend_size'])

    # Returns the figure
    return fig


def collate_effect_size(counts, powers, alpha):
    """Calculates the effects for power values

    Parameters
    ----------
    counts : array
        the number of samples used to calculate the power

    powers : {list, ndarray}
        list of arrays of power values. If there are multiple power arrays,
        each power array should have the same dimensions. If samples are
        missing, these should be deonted by a nan.

    alpha : float
        the critical value used to calculate the power.

    Returns
    -------
    effect_means : 1d array
    effect_bounds : 1d array

    Raises
    ------
    TypeError
        if counts is not a one-dimensional array
    ValueError
        if the arrays in powers have different shapes
    ValueError
        if the length of the power arrays and the length of the count arrays
        are different
    """

    # Checks the power and counts iteratibility
    if isinstance(powers, ndarray):
        powers = [powers]

    num_powers = len(powers)
    if isinstance(counts, ndarray):
        counts = [counts]*num_powers

    # Checks there is a count for each power
    if not len(counts) == len(powers):
        raise ValueError('There must be a counts array for each power array.')

    # Checks the shape array
    for idx in xrange(num_powers):
        count_shape = counts[idx].shape
        power_shape = powers[idx].shape
        # Checks the count array is 1d
        if not len(count_shape) == 1:
            raise TypeError('Each count array must be a 1d array.')

        if len(power_shape) == 1:
            if not count_shape[0] == power_shape[0]:
                raise ValueError('There must be a sample count for each '
                                 'power.')
        elif not count_shape[0] == power_shape[1]:
            raise ValueError('There must be a sample count for each power.')

    # Prealocates the output arrays
    effect_means = zeros((num_powers))
    effect_bounds = zeros((num_powers))

    # Iterates through the powers and calculates the effect sizes
    for idp, pwr in enumerate(powers):
        count = counts[idp]
        pwr_shape = pwr.shape
        # Calculates the effect size for the power array
        eff = zeros(pwr_shape)
        for id2, cnt in enumerate(count):
            if len(pwr_shape) == 1:
                if isnan(pwr[id2]):
                    eff[id2] = nan
                    continue
                try:
                    eff[id2] = ft.solve_power(effect_size=None,
                                              nobs=cnt,
                                              alpha=alpha,
                                              power=pwr[id2])
                except:
                    eff[id2] = nan
            else:
                for id1 in xrange(pwr_shape[0]):
                    if isnan(pwr[id1, id2]):
                        eff[id1, id2] = nan
                        continue
                    try:
                        eff[id1, id2] = ft.solve_power(effect_size=None,
                                                       nobs=cnt,
                                                       alpha=alpha,
                                                       power=pwr[id1, id2])
                    except:
                        eff[id1, id2] = nan
        # Caluclates the mean and bound
        if isnan(eff).all():
            effect_means[idp] = nan
            effect_bounds[idp] = nan
        else:
            effect_means[idp] = nanmean(eff, None)
            effect_bounds[idp] = _confidence_bound(eff, alpha, None)

    return effect_means, effect_bounds


def _check_strs(x):
    r"""Determines if x is a string, number or nan"""

    if isinstance(x, str):
        return True
    elif isnan(x):
        return False
    elif isinstance(x, (float, int)):
        return True
    else:
        raise TypeError('input must be a string, float or a nan')


def _confidence_bound(vec, alpha=0.05, df=None, axis=None):
    r"""Calculates a confidence bound assuming a normal distribution

    Parameters
    ----------
    vec : array

    alpha : {0.05, float}
        the critical value

    df : {None, float}, optional
        the degrees of freedom associated with the distribution. If None is
        given, df is assumed to be the number elements in specified axis.

    axis : {None, unsigned int}, optional
        Default is None. The axis over which to take the devation.

    Return
    ------
    bound : float
        the confidence bound around the mean. The confidence interval is
        [mean - bound, mean + bound].

    """

    # Determines the number of non-nan counts
    vec_shape = vec.shape
    if axis is None and len(vec_shape) == 1:
        num_counts = vec_shape[0] - isnan(vec).sum()
        axis = None
    elif axis is None:
        num_counts = vec_shape[0] * vec_shape[1] - isnan(vec).sum()
    else:
        num_counts = vec_shape[axis] - isnan(vec).sum() / \
            (vec_shape[0] * vec_shape[1])

    # Gets the df if not supplied
    if df is None:
        df = num_counts - 1

    # Calculates the bound
    bound = nanstd(vec, axis=axis) / sqrt(num_counts - 1) * \
        t.ppf(1 - alpha / 2, df)

    return bound


def _calculate_power(p_values, alpha=0.05):
    r"""Calculates statical power empirically

    Parameters
    ----------
    p_values : 1d array

    alpha : float
        the critical value for the power calculation

    Returns
    -------
    power : float
        the emperical power, or the fraction of observed p values below the
        critical value

    """

    w = (p_values < float(alpha)).sum()/float(p_values.shape[0])

    return w


def compare_distributions(test, samples, counts=5, num_iter=1000):
    r"""Compares two distribution arrays iteratively

    Parameters
    ----------
    test : function
        the statistical test which accepts an array-like of sample ids
        (list of lists) and returns a p-value.

    samples : array-like
        samples can be a list of lists or an array where each sublist or row in
        the array corresponds to a sampled group.

    counts : {unsigned int, 1d array}, optional
        Default is 5. The number of samples to draw from each distribution.
        If this is a 1d array, the length must correspond to the number of
        samples.

    num_iter : int, optional
        Default 1000. The number of p-values to generate for each point on the
        curve.

    Returns
    -------
    p_values : array
        the bootstrapped p-values

    Raises
    ------
    ValueError
        if counts is a 1d array and counts and samples are different lengths

    """
    # Determines the number of groups
    num_groups = len(samples)

    # Handles the number of samples for later instances
    if isinstance(counts, int):
        counts = array([counts]*num_groups)
    elif not len(counts) == num_groups:
        raise ValueError('If counts is a 1d array, there must be a count to '
                         'draw for each group.')

    # Prealocates the pvalue matrix
    p_values = zeros((num_iter))

    for idx in range(num_iter):
        subs = [choice(array(pop), counts[i], replace=False)
                for i, pop in enumerate(samples)]
        p_values[idx] = test(subs)

    return p_values


def calculate_power_curve(test, samples, sample_counts, ratio=None,
                          num_iter=1000, alpha=0.05):
    """Generates an empirical power curve for the samples.

    ... ....

    Parameters
    ----------
    test : function
        the statistical test which accepts an array-like of sample ids
        (list of lists) and returns a p-value.

    samples : array-like
        samples can be a list of lists or an array where each sublist or row in
        the array corresponds to a sampled group.

    sample_counts : 1d array
        A vector of the number of samples which should be sampled in each curve

    ratio : {None, 1d array}
        The fraction of the sample counts which should be assigned to each
        group. This must be a none-type object, or the same length as samples.

    num_iter : int
        The default is 1000. The number of p-values to generate for each point
        on the curve.

    Returns
    -------
    p_values : array
        the bootstrapped p-values

    Raises
    ------
    ValueError
        if ratio is an array and ratio is not the same length as samples

    """

    # Determines the number of groups
    num_groups = len(samples)
    num_samps = len(sample_counts)
    if isinstance(alpha, float):
        vec = True
        pwr = zeros((num_samps))
        alpha = array([alpha])
    else:
        vec = False
        num_crit = alpha.shape[0]
        pwr = zeros((num_crit, num_samps))

    # Checks the ratio argument
    if ratio is None:
        ratio = ones((num_groups))
    elif not ratio.shape == (num_groups,):
        raise ValueError('There must be a ratio for each group.')

    # Loops through the sample sizes
    for id2, s in enumerate(sample_counts):
        count = nround(s*ratio, 0).astype(int)
        for id1, a in enumerate(alpha):
            ps = compare_distributions(test, samples, count, num_iter)
            if vec:
                pwr[id2] = _calculate_power(ps, a)
            else:
                pwr[id1, id2] = _calculate_power(ps, a)

    return pwr


def bootstrap_power_curve(test, samples, sample_counts, ratio=None,
                          alpha=0.05, num_iter=500, num_runs=10):
    r"""Repeatedly calculates the power curve for a specified alpha level

    Parameters
    ----------
    test : function
        the statistical test which accepts an array-like of sample ids
        (list of lists) and returns a p-value.

    samples : array-like
        samples can be a list of lists or an array where each sublist or row in
        the array corresponds to a sampled group.

    sample_counts : 1d array
        A vector of the number of samples which should be sampled in each curve

    ratio : {None, 1d array}
        The fraction of the sample counts which should be assigned to each
        group. This must be a none-type object, or the same length as samples.

    num_iter : unsigned int
        Default is 1000. The number of p-values to generate for each point
        on the curve.

    num_runs : unsigned int
        Default is 5. The number of times to calculate each curve.

    Returns
    -------
    p_mean : 1d array
        the mean p-values from the iterations

    p_std : vector
        the variance in the p-values

    Example
    -------
    Suppose we have 100 samples randomly drawn from two normal distribitions,
    the first with mean 0 and standard devation 1, and the second with mean 3
    and standard deviation 1.5

    >>> import numpy as np
    >>> samples_1 = np.random.randn(100)
    >>> samples_2 = 1.5*np.random.randn(100) + 1

    We want to test the statistical power of a kruskal-wallis test comparing
    the two populations. We can define a test function, f, to perform the
    comparison. The test function will take a list of value vectors and
    return a p value.

    >>> from scipy.stats import ttest_ind
    >>> f = lambda x: ttest_ind(x[0], x[1])[1]

    Now, we can determine the statitical power, or the probability that do not
    have a false positive given that we do not have a false negative by varying
    a number of subsamples.

    >>> sample_counts = np.arange(5, 80, 5)
    >>> power_mean, power_bound = bootstrap_power_curve(f,
    ...                                                 [samples_1, samples_2],
    ...                                                 sample_counts)
    >>> print power_mean
        [ 0.2772  0.569   0.7744  0.9052  0.969   0.9898  0.9984  0.9998  1.
          1.      1.      1.      1.      1.      1.    ]
    >>> print power_bound
        [ 0.0178  0.0124  0.0145  0.0097  0.0053  0.0027  0.0013  0.0004  0.
          0.      0.      0.      0.      0.      0.    ]

    """

    # Corrects the alpha value into a matrix
    alpha = ones((num_runs))*alpha

    # Boot straps the power curve
    power = calculate_power_curve(test,
                                  samples,
                                  sample_counts,
                                  ratio,
                                  num_iter,
                                  alpha)

    # Calculates two summary statitics
    power_mean = power.mean(0)
    power_bound = _confidence_bound(power, alpha=alpha[0], axis=0)

    # Calculates summary statitics
    return power_mean, power_bound


def get_significant_subsample(tests, samples, sub_size=None, p_crit=0.05,
                              num_rounds=500, p_scaling=5):
    """
    Subsamples data to an even sample number for all groups

    Parameters
    ----------
    tests : list
        the statistical tests to performed on the data. These tests should
        take a list of integers or sample ids, and return a p-value.

    samples : array-like
        samples can be a list of lists or an array where each sublist or row in
        the array corresponds to a sampled group.

    sub_size : {None, int}
        the maximum number of samples to select from a group. If no value is
        provided, this will be the same as the size of the smallest group.
        Otherwise, this will be compared to the size of the smallest group, and
        which ever is lower will be used.

    p_crit : {float, list}
        The critical p value or p values for the function.

    num_rounds : {500, int}
        the number of times the code should attempt to subsample the data
        before determining it has tried too many times and should quit.

    p_scaling : {5, int}
        a penalty scale on p_crit, so that the total distribution must be less
        than p_crit/p_scaling.

    Returns
    -------
    ids : array
        All sample ids in the subset

    sub_size : float
        the number of samples selected from each group

    Raises
    ------
    RuntimeError
        if all the tests are None, or no signfiiant difference can be found
        between samples

    RuntimeError
        if not iteration can be found that satisfies the signfigiant difference
        between groups

    """

    # Determines the size of the groups
    check_size = array([len(i) for i in samples])
    if sub_size is None:
        sub_size = check_size.min()
    else:
        sub_size = min([sub_size, check_size.min()])

    # Checks the critical value is the same length as the tests
    if isinstance(p_crit, float):
        p_crit = p_crit*ones((len(tests)))

    # Verifies testing is reasonable for the
    for idx, f in enumerate(tests):
        if f is not None and p_crit[idx]/p_scaling < f(samples):
            tests[idx] = None
    # Checks the functions are defined
    if (tests == array([None]*len(tests))).all():
        raise RuntimeError('There is no test defined')

    # Loops through to get a signfigant difference
    for i in xrange(num_rounds+1):
        # Subsamples the larger dataset
        sub_samps = []
        for ids in samples:
            sub_samps.append(choice(ids, size=sub_size, replace=False))

        # Tests the subsample
        test_res = ones((len(tests)))
        for idx, f in enumerate(tests):
            test_res[idx] = f(sub_samps)

        # Checks the critical values have been satisifed
        if (test_res < p_crit).all():
            return sub_samps

        # If no iteration has been found, this is supplied
        elif i == num_rounds:
            raise RuntimeError('There is no iteration which satisfies your '
                               'requirements.')


def get_paired_subsamples(meta, cat, control_cats, order=None, strict=True):
    r"""Gets a set samples to serve as controls

    Parameters
    ----------
    meta : dataframe

    cat : str
        the metadata categories for comparison

    control_cats : list
        the metadata categories to be used as controls. For example, if you
        wanted to control age (`cat` = "AGE"), you might want to control for
        gender and health status (i.e. `control_cats` = ["SEX", "HEALTHY"])

    order : {None, list}
        Default is None. The order of groups in the category. This can be used
        to limit the groups selected. For example, if there's a category with
        groups 'A', 'B' and 'C', and you only want to look at A vs B, `order`
        would be set to ['A', 'B'].

    strict: bool
        Default is True. If a missing value (nan) is encountered, the group
        will be skipped when 'strict' is True.

    Returns
    -------
    ids : array
        a set of arrays which satisfy the criteria. These are not grouped by
        cat. An empty array indicates there are no sample ids which satisfy
        the requirements.

    Example
    -------
    If we have a mapping file for a set of random samples looking at housing,
    sex, age and antibiotic use.

    >>> import pandas as pd
    >>> meta = {'NR': {'HOUSING': '2', 'SEX': 'F', 'AGE': nan, 'ABX': 'Y'},
    ...         'MH': {'HOUSING': '3', 'SEX': 'F', 'AGE': '30s', 'ABX': 'Y'},
    ...         'PP': {'HOUSING': '2', 'SEX': 'F', 'AGE': '30s', 'ABX': 'N'},
    ...         'CD': {'HOUSING': '3', 'SEX': 'F', 'AGE': '30s', 'ABX': 'Y'},
    ...         'MM': {'HOUSING': '1', 'SEX': 'F', 'AGE': '30s', 'ABX': 'Y'},
    ...         'SW': {'HOUSING': '2', 'SEX': 'M', 'AGE': nan, 'ABX': 'N'},
    ...         'TS': {'HOUSING': '2', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'},
    ...         'CB': {'HOUSING': '3', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'},
    ...         'BB': {'HOUSING': '1', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'}}
    >>>  meta = pd.DataFrame.from_dict(meta)

    Let's say we want to vary housing, controlling for sex, age, antibiotics
    and sex.
    >>> ids = get_paired_subsamples(meta, 'HOUSING', ['SEX', 'AGE', 'ABX'])
    >>> print ids
        [['BB'], ['CB'], ['TS']]

    """
    # Groups meta by category
    cat_groups = meta.groupby(cat).groups

    # Handles the order argument
    if order is None:
        order = sorted(cat_groups.keys())
    order = array(order)
    num_groups = len(order)

    # Determines the number of samples, and the experimental and control group
    group_size = array([len(cat_groups[o]) for o in order])
    ctrl_name = order[group_size == group_size.min()][0]
    order = order[order != ctrl_name]

    # Gets a control group table
    ctrl_group = meta.loc[cat_groups[ctrl_name]
                          ].groupby(list(control_cats)).groups
    exp_groups = [meta.loc[cat_groups[o]
                           ].groupby(list(control_cats)).groups for o in order]

    ids = [array([])]*num_groups
    # Loops through samples in the experimental group to match for controls
    for check_group, ctrl_ids in viewitems(ctrl_group):
        # Checks the categories have been defined
        undefed_check = array([_check_strs(p) for p in check_group])
        if not undefed_check.all() and strict:
            continue
        num_samps = len(ctrl_ids)
        exp_ids = []
        # Loops through the other groups
        for exp_group in exp_groups:
            # Checks group to be considered is included in the grouping
            if check_group not in exp_group:
                break
            # Gets the id associated with the group
            pos_ids = exp_group[check_group]
            # Randomly subsamples the possible ids
            num_draw = min([len(pos_ids), num_samps])
            exp_ids.append(choice(ctrl_ids, num_draw, replace=False))
            exp_ids.append(choice(pos_ids, num_draw, replace=False))

        if len(exp_ids) == num_groups:
            for idx in xrange(num_groups):
                ids[idx] = hstack((ids[idx], exp_ids[idx]))

    return ids
