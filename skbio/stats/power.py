r"""
Empirical Power Estimation (:mod:`skbio.stats.power`)
=====================================================

.. currentmodule:: skbio.stats.power

The purpose of this module is to provide empirical, post-hoc power estimation
of normally and non-normally distributed data. It also provides support to
subsample data to facilitate this analysis.

The underlying principle is based on subsampling and Monte Carlo simulation.
Assume that there is some set of populations, :math:`K_{1}, K_{2}, ... K_{n}`
which have some property, :math:`\mu` such that :math:`\mu_{1} \neq \mu_{2}
\neq ... \neq \mu_{n}`. For each of the populations, a sample, :math:`S` can be
drawn, with a parameter, :math:`x` where :math:`x \approx \mu` and for the
samples, we can use a test, :math:`f`, to show that :math:`x_{1} \neq x_{2}
\neq ... \neq x_{n}`.

Since we know that :math:`\mu_{1} \neq \mu_{2} \neq ... \neq \mu_{n}`,
we know we should reject the null hypothesis. If we fail to reject the null
hypothesis, we have committed a Type II error and our result is a false
negative. We can estimate the frequency of Type II errors at various sampling
depths by repeatedly subsampling the populations and observing how often we
see a false negative. If we repeat this several times for each subsampling
depth, and vary the depths we use, we can start to approximate a relationship
between the number of samples we use and the rate of false negatives, also
called the statistical power of the test.

To generate complete power curves from data which appears underpowered, the
`statsmodels.stats.power` package can be used to solve for an effect size. The
effect size can be used to extrapolate a power curve for the data.

Most functions in this module accept a statistical test function which takes a
list of samples and returns a p value. The test is then evaluated over a series
of subsamples.

Sampling may be handled in two ways. For any set of samples, we may simply
choose to draw :math:`n` observations at random for each sample. Alternatively,
if metadata is avalaible, samples can be matched based on a set of control
categories so that paired samples are drawn at random from the set of avaliable
matches.

Functions
---------

.. autosummary::
    :toctree: generated/

    subsample_power
    subsample_paired_power
    confidence_bound
    paired_subsamples
    bootstrap_power_curve

Examples
--------
Suppose we wanted to test that there's a relationship between two random
variables, `ind` and `dep`. Let's use random subsampling to estimate the
statistical power of our test with an alpha of 0.1, 0.01, and 0.001.

To control for the pseudo-random number generation, we will use a seed.
When using these functions with your own data, you don't need to include the
step.

>>> import numpy as np
>>> np.random.seed(20)
>>> ind = np.random.randint(0, 20, 15)
>>> ind
array([ 3, 15,  9, 11,  7,  2,  0,  8, 19, 16,  6,  6, 16,  9,  5])
>>> dep = (3 * ind + 5 + np.random.randn(15) * 5).round(3)
>>> dep
array([ 15.617,  47.533,  28.04 ,  33.788,  19.602,  12.229,   4.779,
        36.838,  67.256,  55.032,  22.157,   7.051,  58.601,  38.664,
        18.783])

Let's define a test that will draw a list of sample pairs and determine
if they're correlated. We'll use `scipy.stats.pearsonr` which takes two arrays
and returns a correlation coefficient and a p-value representing the
probability the two distributions are correlated.

>>> from scipy.stats import pearsonr
>>> f = lambda x: pearsonr(x[0], x[1])[1]

Now, let's use random sampling to estimate the power of our test on
the first distribution.

>>> samples = [ind, dep]
>>> f(samples)
3.6459452596563003e-08

In `subsample_power`, we can maintain a paired relationship between samples
by setting `draw_mode` to "matched". We can also set our critical value, so
that we estimate power for a critical value of :math:`\alpha = 0.05`, an
estimate for the critical value of 0.01, and a critical value of 0.001.

>>> from skbio.stats.power import subsample_power
>>> pwr_100, counts_100 = subsample_power(test=f,
...                                       samples=samples,
...                                       min_observations=3,
...                                       max_counts=10,
...                                       min_counts=3,
...                                       counts_interval=1,
...                                       draw_mode="matched",
...                                       alpha_pwr=0.1)
>>> pwr_010, counts_010 = subsample_power(test=f,
...                                       samples=samples,
...                                       min_observations=3,
...                                       max_counts=10,
...                                       min_counts=3,
...                                       counts_interval=1,
...                                       draw_mode="matched",
...                                       alpha_pwr=0.01)
>>> pwr_001, counts_001 = subsample_power(test=f,
...                                       samples=samples,
...                                       min_observations=3,
...                                       max_counts=10,
...                                       min_counts=3,
...                                       counts_interval=1,
...                                       draw_mode="matched",
...                                       alpha_pwr=0.001)
>>> counts_100
array([3, 4, 5, 6, 7, 8, 9])
>>> pwr_100.mean(0)
array([ 0.4716,  0.8226,  0.9424,  0.986 ,  0.9988,  1.    ,  1.    ])
>>> pwr_010.mean(0)
array([ 0.0492,  0.2368,  0.5462,  0.823 ,  0.9474,  0.9828,  0.9982])
>>> pwr_001.mean(0)
array([ 0.0028,  0.0174,  0.1262,  0.342 ,  0.5928,  0.8256,  0.9594])

Based on this power estimate, as we increase our confidence that we have not
committed a type I error and identified a false positive, the number of samples
we need to be confident that we have not committed a type II error increases.

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils import viewitems
from future.builtins import range

from warnings import warn

import numpy as np
import scipy.stats


def subsample_power(test, samples, draw_mode='ind', alpha_pwr=0.05, ratio=None,
                    min_observations=20, max_counts=50,
                    counts_interval=10, min_counts=None, num_iter=500,
                    num_runs=10):
    r"""Subsamples data to iteratively calculate power

    Parameters
    ----------
    test : function
        The statistical test which accepts a list of arrays of values
        (sample ids or numeric values) and returns a p value or one-dimensional
        array of p values.
    samples : array_like
        `samples` can be a list of lists or a list of arrays where each
        sublist or row in the array corresponds to a sampled group.
    draw_mode : {"ind", "matched"}, optional
        "matched" samples should be used when observations in
        samples have corresponding observations in other groups. For instance,
        this may be useful when working with regression data where
        :math:`x_{1}, x_{2}, ..., x_{n}` maps to
        :math:`y_{1}, y_{2}, ..., y_{n}`. Sample vectors must be the same
        length in "matched" mode.
        If there is no reciprocal relationship between samples, then
        "ind" mode should be used.
    alpha_pwr : float, optional
        The critical value used to calculate the power.
    ratio : 1-D array, optional
        The fraction of the sample counts which should be
        assigned to each group. If this is a 1-D array, it must be the same
        length as `samples`. If no value is supplied (`ratio` is None),
        then an equal number of observations will be drawn for each sample. In
        `matched` mode, this will be set to one.
    min_observations : positive int, optional
        The minimum number of observations in any sample to perform power
        analysis. Note that this is not the same as the minimum number of
        samples drawn per group.
    max_counts : positive int, optional
        The maximum number of samples per group to draw for effect size
        calculation.
    counts_interval : positive int, optional
        The difference between each subsampling count.
    min_counts : positive int, optional
        How many samples should be drawn for the smallest
        subsample. If this is None, the `counts_interval` will be used.
    num_iter : positive int, optional
        The number of p-values to generate for each point
        on the curve.
    num_runs : positive int, optional
        The number of times to calculate each curve.

    Returns
    -------
    power : array
        The power calculated for each subsample at each count.
    sample_counts : array
        The number of samples drawn at each power calculation.

    Raises
    ------
    ValueError
        If the `mode` is "matched", an error will occur if the arrays in
        `samples` are not the same length.
    ValueError
        There is a ValueError if there are fewer samples than the minimum
        count.
    ValueError
        If the `counts_interval` is greater than the difference between the
        sample start and the max value, the function raises a ValueError.
    ValueError
        There are not an equal number of groups in `samples` and in `ratios`.
    TypeError
        `test` does not return a float or a 1-dimensional numpy array.


    Examples
    --------
    Let's say we wanted to look at the relationship between the presence of a
    specific bacteria and the probability of a pre or post menopausal woman
    experiencing a health outcome. Healthy women were enrolled in the study
    either before or after menopause, and followed for five years. They
    submitted fecal samples at regular intervals during that period, and were
    assessed for a particular irreversible health outcome over that period.

    16S sequencing and available literature suggest a set of candidate taxa
    may be associated with the health outcome. Assume there are 100 samples
    (50 premenopausal samples and 50 postmenopausal samples) where the taxa
    of interest was identified by 16S sequencing and the taxonomic abundance
    was confirmed in a certain fraction of samples at a minimum level.

    We can simulate the probability that a woman positive for this taxa
    experiences the health outcome using a binomial distribution.

    >>> import numpy as np
    >>> np.random.seed(25)
    >>> pre_rate = np.random.binomial(1, 0.75, size=(50,))
    >>> pre_rate
    array([0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1,
           0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0,
           1, 1, 1, 1])
    >>> pos_rate = np.random.binomial(1, 0.25, size=(50,))
    >>> pos_rate
    array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
           0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0,
           0, 1, 0, 0])

    Let's set up a test function, so we can test the probability of
    finding a difference in frequency between the two groups. We'll use
    `scipy.stats.chisquare` to look for the difference in frequency between
    groups.

    >>> from scipy.stats import chisquare, nanmean
    >>> test = lambda x: chisquare(np.array([x[i].sum() for i in
    ...     xrange(len(x))]))[1]

    Let's make sure that our two distributions are different.

    >>> round(test([pre_rate, pos_rate]), 5)
    9e-05

    Since there are an even number of samples, and we don't have enough
    information to try controlling the data, we'll use
    `skbio.stats.power.subsample_power` to compare the two groups. If we had
    metadata about other risk factors, like a family history, BMI, tobacco use,
    we might want to use `skbio.stats.power.subsample_paired_power`.
    We'll also use "ind" `draw_mode`, since there is no linkage between the
    two groups of samples.

    >>> from skbio.stats.power import subsample_power
    >>> pwr_est, counts = subsample_power(test=test,
    ...                                   samples=[pre_rate, pos_rate],
    ...                                   counts_interval=5)
    >>> counts
    array([ 5, 10, 15, 20, 25, 30, 35, 40, 45])
    >>> nanmean(pwr_est, 0)
    array([ 0.178 ,  0.3354,  0.658 ,  0.8992,  0.9818,  0.9984,  1.    ,
            1.    ,  1.    ])

    So, we can estimate that we will see a significant difference between
    the two groups (:math:`\alpha \leq 0.05`) at least 80% of the time if we
    use 20 observations per group.

    If we wanted to test the relationship of a second candidate taxa which is
    more rare in the population, but may have a similar effect, based on
    available literature, we might also start by trying to identify 20
    samples per group where the second candidate taxa is present.

    """

    # Checks the inputs
    ratio, num_p, sample_counts = \
        _check_subsample_power_inputs(test=test,
                                      samples=samples,
                                      draw_mode=draw_mode,
                                      ratio=ratio,
                                      min_observations=min_observations,
                                      min_counts=min_counts,
                                      max_counts=max_counts,
                                      counts_interval=counts_interval)

    # Prealocates the power array
    power = np.zeros((num_runs, len(sample_counts), num_p))

    # Calculates the power instances
    for id2, c in enumerate(sample_counts):
        count = np.round(c * ratio, 0).astype(int)
        for id1 in range(num_runs):
            ps = _compare_distributions(test=test,
                                        samples=samples,
                                        num_p=num_p,
                                        counts=count,
                                        num_iter=num_iter,
                                        mode=draw_mode)
            power[id1, id2, :] = _calculate_power(ps, alpha_pwr)

    if num_p == 1:
        power = power[:, :, 0]

    return power, sample_counts


def subsample_paired_power(test, meta, cat, control_cats, order=None,
                           strict_match=True, alpha_pwr=0.05,
                           min_observations=20, max_counts=50,
                           counts_interval=10, min_counts=None,
                           num_iter=500, num_runs=10):
    r"""Estimates power iteratively using samples with matching metadata

    Parameters
    ----------
    test : function
        The statistical test which accepts a list of arrays sample ids and
        returns a p value.
    meta : pandas.DataFrame
        The metadata associated with the samples.
    cat : str
        The metadata category being varied between samples.
    control_cats : list
        The metadata categories to be used as controls. For example, if
        you wanted to vary age (`cat` = "AGE"), you might want to control
        for gender and health status (i.e. `control_cats` = ["SEX",
        "HEALTHY"]).
    order : list, optional
        The order of groups in the category. This can be used
        to limit the groups selected. For example, if there's a category with
        groups 'A', 'B' and 'C', and you only want to look at A vs B, `order`
        would be set to ['A', 'B'].
    strict_match : bool, optional
        This determines how data is grouped using
        `control_cats`. If a sample within `meta` has an undefined value (NaN)
        for any of the columns in `control_cats`, the sample will not be
        considered as having a match and will be ignored when `strict_match`
        is True. If `strict_match` is False, missing values (NaN) in the
        `control_cats` can be considered matches.
    alpha_pwr : float, optional
        The critical value used to calculate the power.
    min_observations : positive int, optional
        The minimum number of paired samples which must exist
        for a category and set of control categories to be able to subsample
        and make power calculations. This is not the same as the minimum
        number of observations to draw during subsampling.
    max_counts : positive int, optional
        The maximum number of observations per sample to draw
        for effect size calculation.
    counts_interval : positive int, optional
        The difference between each subsampling count.
    min_counts : positive int, optional
        How many samples should be drawn for the smallest
        subsample. If this is None, the `counts_interval` will be used.
    num_iter : positive int, optional
        The number of p-values to generate for each point on the curve.
    num_runs : positive int, optional
        The number of times to calculate each curve.

    Returns
    -------
    power : array
        The power calculated for each subsample at each count.
    sample_counts : array
        The number of samples drawn at each power calculation.

    Raises
    ------
    ValueError
        There is a ValueError if there are fewer samples than the minimum
        count.
    ValueError
        If the `counts_interval` is greater than the difference between the
        sample start and the max value, the function raises a ValueError.
    TypeError
        `test` does not return a float or a 1-dimensional numpy array.


    Examples
    --------
    Assume you are interested in the role of a specific cytokine of protein
    translocation in myloid-lineage cells. You are able to culture two
    macrophage lineages (bone marrow derived phagocytes and
    peritoneally-derived macrophages). Due to unfortunate circumstances, your
    growth media must be acquired from multiple sources (lab, company A,
    company B). Also unfortunate, you must use labor-intense low throughput
    assays. You have some preliminary measurements, and you'd like to
    predict how many (more) cells you need to analyze for 80% power.

    You have information about 60 cells, which we'll simulate below. Note
    that we are setting a random seed value for consistency.

    >>> import numpy as np
    >>> import pandas as pd
    >>> np.random.seed(25)
    >>> data = pd.DataFrame.from_dict({
    ...     'CELL_LINE': np.random.binomial(1, 0.5, size=(60,)),
    ...     'SOURCE': np.random.binomial(2, 0.33, size=(60,)),
    ...     'TREATMENT': np.hstack((np.zeros((30)), np.ones((30)))),
    ...     'INCUBATOR': np.random.binomial(1, 0.2, size=(60,))})
    >>> data['OUTCOME'] = (0.25 + data.TREATMENT * 0.25) + \
    ...     np.random.randn(60) * (0.1 + data.SOURCE/10 + data.CELL_LINE/5)
    >>> data.loc[data.OUTCOME < 0, 'OUTCOME'] = 0
    >>> data.loc[data.OUTCOME > 1, 'OUTCOME'] = 1

    We will approach this by assuming that the distribution of our outcome is
    not normally distributed, and apply a kruskal-wallis test to compare
    between the cytokine treated and untreated cells.

    >>> from scipy.stats import kruskal
    >>> f = lambda x: kruskal(*[data.loc[i, 'OUTCOME'] for i in x])[1]

    Let's check that cytokine treatment has a signifigant effect across all
    the cells.

    >>> treatment_stat = [g for g in data.groupby('TREATMENT').groups.values()]
    >>> f(treatment_stat)
    0.0019386336266250209

    Now, let's pick the control categories. It seems reasonable to assume there
    may be an effect of cell line on the treatment outcome, which may be
    attributed to differences in receptor expression. It may also be possible
    that there are differences due cytokine source. Incubators were maintained
    under the same conditions throughout the experiment, within one degree of
    temperature difference at any given time, and the same level of CO2.
    So, at least initially, let's ignore differences due to the incubator.

    It's recommended that as a first pass analysis, control variables be
    selected based on an idea of what may be biologically relevant to the
    system, although further iteration might encourage the consideration of
    variable with effect sizes similar, or larger than the variable of
    interest.

    >>> control_cats = ['SOURCE', 'CELL_LINE']
    >>> from skbio.stats.power import subsample_paired_power
    >>> pwr, cnt = subsample_paired_power(test=f,
    ...                                   meta=data,
    ...                                   cat='TREATMENT',
    ...                                   control_cats=control_cats,
    ...                                   min_observations=5,
    ...                                   counts_interval=5,
    ...                                   num_iter=100,
    ...                                   num_runs=5)
    >>> cnt
    array([ 5, 10, 15, 20])
    >>> pwr.mean(0)
    array([ 0.19 ,  0.358,  0.708,  0.772])
    >>> pwr.std(0).round(3)
    array([ 0.061,  0.075,  0.114,  0.181])

    Estimating off the power curve, it looks like 20 cells per group may
    provide addiquite power for this experiment, although the large variance
    in power might suggest extending the curves or increasing the number of
    samples per group.

    """

    # Checks for the number of sampling pairs avaliable
    sub_ids = paired_subsamples(meta, cat, control_cats, order, strict_match)

    ratio, num_p, sample_counts = \
        _check_subsample_power_inputs(test=test,
                                      samples=sub_ids,
                                      draw_mode='matched',
                                      min_observations=min_observations,
                                      min_counts=min_counts,
                                      max_counts=max_counts,
                                      counts_interval=counts_interval)

    # Prealocates the power array
    power = np.zeros((num_runs, len(sample_counts), num_p))

    # Calculates power instances
    for id2, c in enumerate(sample_counts):
        count = np.round(c * ratio, 0).astype(int)
        for id1 in range(num_runs):
            sub_ids = paired_subsamples(meta, cat, control_cats, order,
                                        strict_match)
            ps = _compare_distributions(test=test,
                                        samples=sub_ids,
                                        num_p=num_p,
                                        counts=count,
                                        num_iter=num_iter,
                                        mode="matched")
            power[id1, id2, :] = _calculate_power(ps, alpha_pwr)

    if num_p == 1:
        power = power[:, :, 0]

    return power, sample_counts


def confidence_bound(vec, alpha=0.05, df=None, axis=None):
    r"""Calculates a confidence bound assuming a normal distribution

    Parameters
    ----------
    vec : array_like
        The array of values to use in the bound calculation.
    alpha : float, optional
        The critical value, used for the confidence bound calculation.
    df : float, optional
        The degrees of freedom associated with the
        distribution. If None is given, df is assumed to be the number of
        elements in specified axis.
    axis : positive int, optional
        The axis over which to take the deviation. When axis
        is None, a single value will be calculated for the whole matrix.

    Returns
    -------
    bound : float
        The confidence bound around the mean. The confidence interval is
        [mean - bound, mean + bound].

    """

    # Determines the number of non-nan counts
    vec = np.asarray(vec)
    vec_shape = vec.shape
    if axis is None and len(vec_shape) == 1:
        num_counts = vec_shape[0] - np.isnan(vec).sum()
    elif axis is None:
        num_counts = vec_shape[0] * vec_shape[1] - np.isnan(vec).sum()
    else:
        num_counts = vec_shape[axis] - np.isnan(vec).sum() / \
            (vec_shape[0] * vec_shape[1])

    # Gets the df if not supplied
    if df is None:
        df = num_counts - 1

    # Calculates the bound
    bound = scipy.stats.nanstd(vec, axis=axis) / np.sqrt(num_counts - 1) * \
        scipy.stats.t.ppf(1 - alpha / 2, df)

    return bound


def bootstrap_power_curve(test, samples, sample_counts, ratio=None,
                          alpha=0.05, mode='ind', num_iter=500, num_runs=10):
    r"""Repeatedly calculates the power curve for a specified alpha level

    .. note:: Deprecated in scikit-bio 0.2.3-dev
       ``bootstrap_power_curve`` will be removed in scikit-bio 0.3.1. It is
       Deprecated in favor of using ``subsample_power`` or
       ``sample_paired_power`` to calculate a power array, and then using
       ``confidence_bound`` to perform bootstrapping.

    Parameters
    ----------
    test : function
        The statistical test which accepts an array_like of sample ids
        (list of lists or arrays) and returns a p-value.
    samples : array_like
        samples can be a list of lists or an array where each sublist or row in
        the array corresponds to a sampled group.
    sample_counts : 1-D array_like
        A vector of the number of samples which should be sampled in each curve
    ratio : 1-D array_like, optional
        The fraction of the sample counts which should be
        assigned to each
        group. This must be a none-type object, or the same length as samples.
        If Ratio is None, the same number of observations are drawn from
        each sample.
    alpha : float, optional
        The default is 0.05. The critical value for calculating power.
    mode : {"ind", "matched"}, optional
        "matched" samples should be used when observations in
        samples have corresponding observations in other groups. For instance,
        this may be useful when working with regression data where
        :math:`x_{1}, x_{2}, ..., x_{n}` maps to :math:`y_{1}, y_{2}, ... ,
        y_{n}`.
    num_iter : positive int, optional
        The number of p-values to generate for each point on the curve.
    num_runs : positive int, optional
        The number of times to calculate each curve.

    Returns
    -------
    power_mean : 1-D array
        The mean p-values from the iterations.
    power_bound : vector
        The variance in the p-values.

    Examples
    --------
    Suppose we have 100 samples randomly drawn from two normal distribitions,
    the first with mean 0 and standard devation 1, and the second with mean 3
    and standard deviation 1.5

    >>> import numpy as np
    >>> np.random.seed(20)
    >>> samples_1 = np.random.randn(100)
    >>> samples_2 = 1.5 * np.random.randn(100) + 1

    We want to test the statistical power of a independent two sample t-test
    comparing the two populations. We can define an anonymous function, `f`,
    to wrap the scipy function for independent t tests,
    `scipy.stats.ttest_ind`. The test function will take a list of value
    vectors and return a p value.

    >>> from scipy.stats import ttest_ind
    >>> f = lambda x: ttest_ind(x[0], x[1])[1]

    Now, we can determine the statistical power, or the probability that we do
    not have a false negative given that we do not have a false positive, by
    varying a number of subsamples.

    >>> from skbio.stats.power import bootstrap_power_curve
    >>> sample_counts = np.arange(5, 80, 5)
    >>> power_mean, power_bound = bootstrap_power_curve(f,
    ...                                                 [samples_1, samples_2],
    ...                                                 sample_counts)
    >>> sample_counts[power_mean - power_bound.round(3) > .80].min()
    20

    Based on this analysis, it looks like we need at least 20 observations
    from each distribution to avoid committing a type II error more than 20%
    of the time.

    """

    warn("skbio.stats.power.bootstrap_power_curve is deprecated. Please "
         "use skbio.stats.power.subsample_power or "
         "skbio.stats.power.subsample_paired_power followed by "
         "confidence_bound.", DeprecationWarning)

    # Corrects the alpha value into a matrix
    alpha = np.ones((num_runs)) * alpha

    # Boot straps the power curve
    power = _calculate_power_curve(test=test,
                                   samples=samples,
                                   sample_counts=sample_counts,
                                   ratio=ratio,
                                   num_iter=num_iter,
                                   alpha=alpha,
                                   mode=mode)
    # Calculates two summary statitics
    power_mean = power.mean(0)
    power_bound = confidence_bound(power, alpha=alpha[0], axis=0)

    # Calculates summary statitics
    return power_mean, power_bound


def paired_subsamples(meta, cat, control_cats, order=None, strict_match=True):
    r"""Gets a set of samples to serve as controls

    This function is designed to provide controlled samples, based on a
    metadata category. For example, one could control for age, sex, education
    level, and diet type while measuring exercise frequency. No outcome
    value is considered in this subsampling process.

    Parameters
    ----------
    meta : pandas.DataFrame
        The metadata associated with the samples.
    cat : str, list
        The metadata category (or a list of categories) for comparison.
    control_cats : list
        The metadata categories to be used as controls. For example, if you
        wanted to vary age (`cat` = "AGE"), you might want to control for
        gender and health status (i.e. `control_cats` = ["SEX", "HEALTHY"])
    order : list, optional
        The order of groups in the category. This can be used
        to limit the groups selected. For example, if there's a category with
        groups 'A', 'B' and 'C', and you only want to look at A vs B, `order`
        would be set to ['A', 'B'].
    strict_match: bool, optional
        This determines how data is grouped using
        `control_cats`. If a sample within `meta` has an undefined value (NaN)
        for any of the columns in `control_cats`, the sample will not be
        considered as having a match and will be ignored when `strict_match`
        is True. If `strict_match` is False, missing values (NaN) in the
        `control_cats` can be considered matches.

    Returns
    -------
    ids : array
        a set of ids which satisfy the criteria. These are not grouped by
        `cat`. An empty array indicates there are no sample ids which satisfy
        the requirements.

    Examples
    --------
    If we have a mapping file for a set of random individuals looking at
    housing, sex, age and antibiotic use.

    >>> import pandas as pd
    >>> import numpy as np
    >>> meta = {'SW': {'HOUSING': '2', 'SEX': 'M', 'AGE': np.nan, 'ABX': 'Y'},
    ...         'TS': {'HOUSING': '2', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'},
    ...         'CB': {'HOUSING': '3', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'},
    ...         'BB': {'HOUSING': '1', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'}}
    >>> meta = pd.DataFrame.from_dict(meta, orient="index")
    >>> meta #doctest: +SKIP
       ABX HOUSING  AGE SEX
    BB   Y       1  40s   M
    CB   Y       3  40s   M
    SW   Y       2  NaN   M
    TS   Y       2  40s   M

    We may want to vary an individual's housing situation, while holding
    constant their age, sex and antibiotic use so we can estimate the effect
    size for housing, and later compare it to the effects of other variables.

    >>> from skbio.stats.power import paired_subsamples
    >>> ids = paired_subsamples(meta, 'HOUSING', ['SEX', 'AGE', 'ABX'])
    >>> np.hstack(ids) #doctest: +ELLIPSIS
    array(['BB', 'TS', 'CB']...

    So, for this set of data, we can match TS, CB, and BB based on their age,
    sex, and antibiotic use. SW cannot be matched in either group becuase
    `strict_match` was true, and there is missing AGE data for this sample.

    """

    # Sets the index data
    # Groups meta by category
    cat_groups = meta.groupby(cat).groups

    # Handles the order argument
    if order is None:
        order = sorted(cat_groups.keys())
    order = np.array(order)
    num_groups = len(order)

    # Determines the number of samples, and the experimental and control group
    group_size = np.array([len(cat_groups[o]) for o in order])
    ctrl_name = order[group_size == group_size.min()][0]
    order = order[order != ctrl_name]

    # Gets a control group table
    ctrl_match_groups = meta.groupby(control_cats).groups
    ctrl_group = meta.loc[cat_groups[ctrl_name]
                          ].groupby(list(control_cats)).groups

    ids = [np.array([])] * num_groups
    # Loops through samples in the experimental group to match for controls
    for check_group, ctrl_ids in viewitems(ctrl_group):
        # Checks the categories have been defined
        undefed_check = np.array([_check_strs(p) for p in check_group])
        if not undefed_check.all() and strict_match:
            continue
        # Removes the matched ids from order
        matched_ids = ctrl_match_groups[check_group]
        for id_ in ctrl_ids:
            matched_ids.remove(id_)
        pos_ids = []
        num_ids = [len(ctrl_ids)]
        # Gets the matrix of the matched ids and groups them
        exp_group = meta.loc[matched_ids].groupby(cat).groups
        for grp in order:
            # Checks group to be considered is included in the grouping
            if grp not in exp_group:
                break
            # Gets the id associated with the group
            pos_ids.append(exp_group[grp])
            num_ids.append(len(exp_group[grp]))
        # Determines the minimum number of samples
        num_draw = np.array(num_ids).min()
        # Draws samples from possible ids
        exp_ids = [np.random.choice(ctrl_ids, num_draw, replace=False)]
        exp_ids.extend([np.random.choice(id_, num_draw, replace=False)
                        for id_ in pos_ids])

        if len(exp_ids) == num_groups:
            for idx in range(num_groups):
                ids[idx] = np.hstack((ids[idx], exp_ids[idx]))

    return ids


def _check_strs(x):
    r"""Returns False if x is a nan and True is x is a string or number"""

    if isinstance(x, str):
        return True
    elif isinstance(x, (float, int)):
        return not np.isnan(x)
    else:
        raise TypeError('input must be a string, float or a nan')


def _calculate_power(p_values, alpha=0.05):
    r"""Calculates statical power empirically

    Parameters
    ----------
    p_values : 1-D array
        A 1-D numpy array with the test results.

    alpha : float
        The critical value for the power calculation.

    Returns
    -------
    power : float
        The emperical power, or the fraction of observed p values below the
        critical value.

    """

    if p_values.ndim == 1:
        p_values = p_values * np.array([[1]])

    w = (p_values < float(alpha)).sum(axis=1)/float(p_values.shape[1])

    return w


def _compare_distributions(test, samples, num_p, counts=5, mode="ind",
                           num_iter=1000):
    r"""Compares two distribution arrays iteratively

    Parameters
    ----------
    test : function
        The statistical test which accepts an array_like of sample ids
        (list of lists) and returns a p-value. This can be a one-dimesnional
        array, or a float.
    samples : list of arrays
        A list where each 1-d array represents a sample. If `mode` is
        "matched", there must be an equal number of observations in each
        sample.
    num_p : positive int, optional
        The number of p-values returned by the test.
    counts : positive int or 1-D array, optional
        The number of samples to draw from each distribution.
        If this is a 1-D array, the length must correspond to the number of
        samples. The function will not draw more observations than are in a
        sample. In "matched" `mode`, the same number of observations will be
        drawn from each group.
    mode : {"ind", "matched"}, optional
        "matched" samples should be used when observations in
        samples have corresponding observations in other groups. For instance,
        this may be useful when working with regression data where
        :math:`x_{1}, x_{2}, ..., x_{n}` maps to :math:`y_{1}, y_{2}, ... ,
        y_{n}`.
    num_iter : positive int, optional
        Default 1000. The number of p-values to generate for each point on the
        curve.

    Returns
    -------
    p_values : array
        The p-values for `n_iter` subsampled tests.

    Raises
    ------
    ValueError
        If mode is not "ind" or "matched".
    ValueError
        If the arrays in samples are not the same length in "matched" mode.
    ValueError
        If counts is a 1-D array and counts and samples are different lengths.

    """

    # Determines the number of groups
    num_groups = len(samples)

    # Checks the mode
    if mode not in {'ind', 'matched'}:
        raise ValueError('Supported sample modes are "ind" and "matched".')

    # Handles the number of samples for later instances
    if isinstance(counts, int):
        counts = np.array([counts] * num_groups)

    if not len(counts) == num_groups:
        raise ValueError('If counts is a 1-D array, there must be a count to'
                         ' draw for each group.')

    # Checks the group length
    samp_lens = [len(sample) for sample in samples]
    # Checks the group length
    if mode == 'matched' and np.array([samp_lens[i] != samp_lens[i+1] for i in
                                       range(num_groups-1)]).all():
        raise ValueError('In "matched" mode, each sample must have the same'
                         ' number of observations.')
    if np.array([samp_lens[i] < counts[i] for i in range(num_groups)]).any():
        raise ValueError('You cannot choose more observations that exist '
                         'in a sample.')

    # Prealocates the pvalue matrix
    p_values = np.zeros((num_p, num_iter))

    for idx in range(num_iter):
        if mode == "matched":
            pos = np.random.choice(np.arange(0, samp_lens[0]), counts[0],
                                   replace=False)
            subs = [sample[pos] for sample in samples]
        else:
            subs = [np.random.choice(np.array(pop), counts[i], replace=False)
                    for i, pop in enumerate(samples)]

        p_values[:, idx] = test(subs)

    if num_p == 1:
        p_values = np.hstack(p_values)

    return p_values


def _check_subsample_power_inputs(test, samples, draw_mode='ind', ratio=None,
                                  min_observations=20, max_counts=50,
                                  counts_interval=10, min_counts=None):
    r"""Makes sure that everything is sane before power calculations

    Parameters
    ----------
    test : function
        The statistical test which accepts a list of arrays of values
        (sample ids or numeric values) and returns a p value or one-dimensional
        array of p values.
    samples : array_like
        `samples` can be a list of lists or a list of arrays where each
        sublist or row in the array corresponds to a sampled group.
    draw_mode : {"ind", "matched"}, optional
        "matched" samples should be used when observations in
        samples have corresponding observations in other groups. For instance,
        this may be useful when working with regression data where
        :math:`x_{1}, x_{2}, ..., x_{n}` maps to
        :math:`y_{1}, y_{2}, ..., y_{n}`. Sample vectors must be the same
        length in "matched" mode.
        If there is no reciprocal relationship between samples, then
        "ind" mode should be used.
    ratio : 1-D array, optional
        The fraction of the sample counts which should be
        assigned to each group. If this is a 1-D array, it must be the same
        length as `samples`. If no value is supplied (`ratio` is None),
        then an equal number of observations will be drawn for each sample. In
        `matched` mode, this will be set to one.
    min_observations : positive int, optional
        The minimum number of observations in any sample to perform power
        analysis. Note that this is not the same as the minimum number of
        samples drawn per group.
        max_counts : positive int, optional
        The maximum number of samples per group to draw for effect size
        calculation.
    counts_interval : positive int, optional
        The difference between each subsampling count.
    min_counts : positive int, optional
        How many samples should be drawn for the smallest
        subsample. If this is None, the `counts_interval` will be used.

    Returns
    -------
    ratio : 1-D array
        The fraction of the sample counts which should be assigned to each
        group.
    num_p : positive integer
        The number of p values returned by `test`.
    sample_counts : array
        The number of samples drawn at each power calculation.

    Raises
    ------
    ValueError
        If the `mode` is "matched", an error will occur if the arrays in
        `samples` are not the same length.
    ValueError
        There is a ValueError if there are fewer samples than the minimum
        count.
    ValueError
        If the `counts_interval` is greater than the difference between the
        sample start and the max value, the function raises a ValueError.
    ValueError
        There are not an equal number of groups in `samples` and in `ratios`.
    TypeError
        `test` does not return a float or a 1-dimensional numpy array.

    """

    # Determines the minimum number of ids in a category
    num_ids = np.array([len(id_) for id_ in samples]).min()
    # Determines the number of groups
    num_groups = len(samples)

    # Checks there are enough samples to subsample
    if num_ids <= min_observations:
        raise ValueError('There are not enough samples for subsampling.')

    # Checks that "matched" mode is handled appropriately
    if draw_mode == "matched":
        for id_ in samples:
            if not len(id_) == num_ids:
                raise ValueError('Each vector in samples must be the same '
                                 'length in "matched" draw_mode.')

    # Checks the number of counts is appropriate
    if min_counts is None:
        min_counts = counts_interval
    if (max_counts - min_counts) < counts_interval:
        raise ValueError("No subsamples of the specified size can be drawn.")

    # Checks the ratio argument is sane
    if ratio is None or draw_mode == 'matched':
        ratio = np.ones((num_groups))
    else:
        ratio = np.asarray(ratio)
    if not ratio.shape == (num_groups,):
        raise ValueError('There must be a ratio for each group.')

    # Determines the number of p values returned by the test
    p_return = test(samples)
    if isinstance(p_return, float):
        num_p = 1
    elif isinstance(p_return, np.ndarray) and len(p_return.shape) == 1:
        num_p = p_return.shape[0]
    else:
        raise TypeError('test must return a float or one-dimensional array.')

    # Calculates the same counts
    sample_counts = np.arange(min_counts,
                              min(max_counts, num_ids),
                              counts_interval)

    return ratio, num_p, sample_counts


def _calculate_power_curve(test, samples, sample_counts, ratio=None,
                           mode='ind', num_iter=1000, alpha=0.05):
    r"""Generates an empirical power curve for the samples.

    Parameters
    ----------
    test : function
        The statistical test which accepts an list of arrays of values and
        returns a p value.
    samples : array_like
        `samples` can be a list of lists or an array where each sublist or row
        in the array corresponds to a sampled group.
    sample_counts : 1-D array
        A vector of the number of samples which should be sampled in each
        curve.
    mode : {"ind", "matched"}, optional
        "matched" samples should be used when observations in
        samples have corresponding observations in other groups. For instance,
        this may be useful when working with regression data where
        :math:`x_{1}, x_{2}, ..., x_{n}` maps to :math:`y_{1}, y_{2}, ... ,
        y_{n}`.
    ratio : 1-D array, optional
        The fraction of the sample counts which should be
        assigned to each group. If this is a 1-D array, it must be the same
        length as `samples`. If no value is supplied (`ratio` is None),
        then an equal number of observations will be drawn for each sample.
    num_iter : int
        The default is 1000. The number of p-values to generate for each point
        on the curve.

    Returns
    -------
    p_values : array
        The p-values associated with the input sample counts.
    Raises
    ------
    ValueError
        If ratio is an array and ratio is not the same length as samples
    """
    # Casts array-likes to arrays
    sample_counts = np.asarray(sample_counts)
    # Determines the number of groups
    num_groups = len(samples)
    num_samps = len(sample_counts)
    if isinstance(alpha, float):
        vec = True
        pwr = np.zeros((num_samps))
        alpha = np.array([alpha])
    else:
        vec = False
        num_crit = alpha.shape[0]
        pwr = np.zeros((num_crit, num_samps))
    # Checks the ratio argument
    if ratio is None:
        ratio = np.ones((num_groups))
    ratio = np.asarray(ratio)
    if not ratio.shape == (num_groups,):
        raise ValueError('There must be a ratio for each group.')

    # Loops through the sample sizes
    for id2, s in enumerate(sample_counts):
        count = np.round(s * ratio, 0).astype(int)
        for id1, a in enumerate(alpha):
            ps = _compare_distributions(test=test,
                                        samples=samples,
                                        counts=count,
                                        num_p=1,
                                        num_iter=num_iter,
                                        mode=mode)
            if vec:
                pwr[id2] = _calculate_power(ps, a)
            else:
                pwr[id1, id2] = _calculate_power(ps, a)

    return pwr
