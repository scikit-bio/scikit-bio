r"""Empirical Power Estimation (:mod: `skbio.stats.power`)
=========================================================

.. currentmodule:: skbio.stats.power

The purpose of this module is to provide empirical, post-hoc power estimation
of microbiome data. It also provides support to subsample data to faciliate
this analysis.

The power estimates here can use be used in conjunction with the statsmodel
power module to estimate effect size.

The general format is to define a test for the data which takes a list of
values (i.e. sample ids) and returns a p-value.

"""

# -----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# -----------------------------------------------------------------------------
from future.utils import viewitems
from numpy import (array, zeros, ones, round as nround, hstack, isnan, ndarray,
                   nan, nanmean, nanstd, sqrt)
from numpy.random import choice
from scipy.stats import t
from statsmodels.stats.power import FTestAnovaPower
from matplotlib import rcParams

# Gets a power solving instance
ft = FTestAnovaPower()

# Sets up plotting parameters so that the default setting is use to Helvetica
# in plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Helvetica', 'Arial']
rcParams['text.usetex'] = True


def collate_effect_size(counts, powers, alpha):
    """ """

    if isinstance(powers, ndarray):
        powers = [powers]
    num_powers = len(powers)
    # Prealocates the output arrays
    effects = array(num_powers, counts.shape[0])

    # Iterates through the powers and calculates the effect sizes
    for id2, count in enumerate(counts):
        for id1, pwr in enumerate(powers):
            try:
                effects[id1, id2] = ft.solve_power(effect_size=None,
                                                   nobs=count,
                                                   alpha=alpha,
                                                   power=pwr[id2])
            except:
                effects[id1, id2] = nan

    # Calculates the mean effect for each group
    effect_means = zeros((num_powers))
    effect_bounds = zeros((num_powers))
    for id1 in xrange(num_powers):
        if isnan(effects[id1, :]).all():
            effect_means[id1] = nan
            effect_bounds[id1] = nan
        else:
            effect_means[id1] = nanmean(effects[id1, :])
            effect_bounds[id1] = _confidence_bound(effects[id1, :],
                                                   alpha=alpha)
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


def _confidence_bound(vec, alpha=0.05, df=None):
    r"""Calculates a confidence bound assuming a normal distribution

    Parameters
    ----------
    vec : 1d array

    alpha : {0.05, float}
        the critical value

    df : {None, float}
        the degrees of freedom associated with the distribution. If None is
        given, df is assumed to be the number elements in the vector - 1.

    Return
    ------
    bound : float
        the confidence bound around the mean

    """

    # Determines the number of non-nan counts
    num_counts = len(vec) - isnan(vec).sum()

    # Gets the df if not supplied
    if df is None:
        df = num_counts - 1

    # Calculates the bound
    bound = nanstd(vec) / sqrt(num_counts - 1) * t.ppf((1 - alpha / 2), df)

    return bound


def _calculate_power(p_values, alpha=0.05):
    r"""Calculates statical power empirically

    Parameters
    ----------
    p_values : 1d array-like
    alpha : float
        the critical value for the power calculation

    Returns
    -------
    power : float
        the emperical power, or the fraction of observed p values below the
        critical value

    """
    if isinstance(p_values, list):
        p_values = array(p_values)

    w = (p_values < float(alpha)).sum()/float(len(p_values))

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
    power_bound = _confidence_bound(power, alpha=alpha[0])

    # Calculates summary statitics
    return power_mean, power_bound


def get_signifigant_subsample(tests, samples, sub_size=None, p_crit=0.05,
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

    """

    # Determines the size of the groups
    check_size = array([len(i) for i in samples])
    if sub_size is None:
        sub_size = check_size.min()
    else:
        sub_size = min([sub_size, check_size])

    # Checks the critical value is the same length as the tests
    if isinstance(p_crit, float) and isinstance(tests, list):
        p_crit = p_crit*ones((len(tests)))
    elif isinstance(p_crit, list):
        p_crit = array(p_crit)

    # Verifies testing is reasonable for the
    for idx, f in enumerate(tests):
        if f is not None and p_crit[idx]/p_scaling < f(samples):
            tests[idx] = None
            print f(samples)
    # Checks the functions are defined
    if (tests == array([None]*len(tests))).all():
        raise ValueError('There is no test defined')

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
            break
        # If no iteration has been found, this is supplied
        elif i == num_rounds:
            raise ValueError('There is no iteration which satisfies your '
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
    ctrl_group = order[group_size == group_size.min()][0]
    order = order[order != ctrl_group]

    # Gets a control group table
    ctrl_group = meta.loc[cat_groups[ctrl_group]].groupby(control_cats).groups
    exp_groups = [meta.loc[cat_groups[o]].groupby(control_cats).groups
                  for o in order]

    ids = [array([])]*num_groups
    # Loops through samples in the experimental group to match for controls
    for check_group, ctrl_ids in viewitems(ctrl_group):
        # Checks the categories have been defined
        undefed_check = array([_check_strs(p) for p in check_group])
        if not undefed_check.all() and strict:
            continue
        num_samps = len(ctrl_ids)
        exp_ids = [array(ctrl_ids)]
        # Loops through the other groups
        for exp_group in exp_groups:
            # Checks group to be considered is included in the grouping
            if check_group not in exp_group:
                break
            # Gets the id associated with the group
            pos_ids = exp_group[check_group]
            # Randomly subsamples the possible ids
            if len(pos_ids) < num_samps:
                continue
            exp_ids.append(choice(pos_ids, num_samps, replace=False))

        if len(exp_ids) == num_groups:
            for idx in xrange(num_groups):
                ids[idx] = hstack((ids[idx], exp_ids[idx]))

    return ids


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
    fig = ft.plot_power('nobs', sample_counts, effect_means - effect_bounds,
                        alpha=kwds['alpha'])
    # Gets the axis of the first plot and its position
    lax = fig.axes[0]
    # Makes the lower bound lines dashed and thin, and changes the color if
    # desired
    for idx, l in enumerate(lax.get_lines()):
        l.set_linestyle(':')
        l.set_linewidth(1.5)
        if kwds['colormap'] is not None:
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
        if kwds['colormap'] is not None:
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
    if kwds['colormap'] is not None:
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
