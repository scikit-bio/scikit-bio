r"""Empirical Power Estimation (:mod:`skbio.stats.power`)
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
if metadata is available, samples can be matched based on a set of control
categories so that paired samples are drawn at random from the set of available
matches.

Functions
---------

.. autosummary::
    :toctree:

    subsample_power
    subsample_paired_power
    confidence_bound
    paired_subsamples

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
>>> ind # doctest: +ELLIPSIS
array([ 3, 15,  9, 11,  7,  2,  0,  8, 19, 16,  6,  6, 16,  9,  5]...
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
>>> print("%.3e" % f(samples))
3.646e-08

In `subsample_power`, we can maintain a paired relationship between samples
by setting `draw_mode` to "matched". We can also set our critical value, so
that we estimate power for a critical value of :math:`\alpha = 0.05`, an
estimate for the critical value of 0.01, and a critical value of 0.001.

>>> from skbio.stats.power import subsample_power
>>> pwr_100, counts_100 = subsample_power(test=f,
...                                       samples=samples,
...                                       max_counts=10,
...                                       min_counts=3,
...                                       counts_interval=1,
...                                       draw_mode="matched",
...                                       alpha_pwr=0.1,
...                                       num_iter=25)
>>> pwr_010, counts_010 = subsample_power(test=f,
...                                       samples=samples,
...                                       max_counts=10,
...                                       min_counts=3,
...                                       counts_interval=1,
...                                       draw_mode="matched",
...                                       alpha_pwr=0.01,
...                                       num_iter=25)
>>> pwr_001, counts_001 = subsample_power(test=f,
...                                       samples=samples,
...                                       max_counts=10,
...                                       min_counts=3,
...                                       counts_interval=1,
...                                       draw_mode="matched",
...                                       alpha_pwr=0.001,
...                                       num_iter=25)
>>> counts_100
array([3, 4, 5, 6, 7, 8, 9])
>>> pwr_100.mean(0)
array([ 0.484,  0.844,  0.932,  0.984,  1.   ,  1.   ,  1.   ])
>>> pwr_010.mean(0)
array([ 0.044,  0.224,  0.572,  0.836,  0.928,  0.996,  1.   ])
>>> pwr_001.mean(0)
array([ 0.   ,  0.016,  0.108,  0.332,  0.572,  0.848,  0.956])

Based on this power estimate, as we increase our confidence that we have not
committed a type I error and identified a false positive, the number of samples
we need to be confident that we have not committed a type II error increases.


"""  # noqa: D205, D415

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections
import copy

import numpy as np
import scipy.stats


def subsample_power(
    test,
    samples,
    draw_mode="ind",
    alpha_pwr=0.05,
    ratio=None,
    max_counts=50,
    counts_interval=10,
    min_counts=None,
    num_iter=500,
    num_runs=10,
):
    r"""Subsample data to iteratively calculate power.

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
        The power calculated for each subsample at each count. The array has
        `num_runs` rows, a length with the same number of elements as
        `sample_counts` and a depth equal to the number of p values returned by
        `test`. If `test` returns a float, the returned array will be
        two-dimensional instead of three.
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
    specific bacteria, *Gardnerella vaginalis* in the vaginal community, and
    the probability of a pre or post menopausal woman experiencing a urinary
    tract infection (UTI). Healthy women were enrolled in the study either
    before or after menopause, and followed for eight weeks. Participants
    submitted fecal samples at the beginning of the study, and were then
    followed for clinical symptoms of a UTI. A confirmed UTI was an endpoint
    in the study.

    Using available literature and 16S sequencing, a set of candidate taxa were
    identified as correlated with UTIs, including *G. vaginalis*. In the 100
    women (50 premenopausal and 50 postmenopausal samples) who had UTIs, the
    presence or absence of *G. vaginalis* was confirmed with quantitative PCR.

    We can model the probability that detectable *G. vaginalis* was found in
    these samples using a binomial model. (*Note that this is a simulation.*)

    >>> import numpy as np
    >>> np.random.seed(25)
    >>> pre_rate = np.random.binomial(1, 0.85, size=(50,))
    >>> pre_rate.sum()
    45
    >>> pos_rate = np.random.binomial(1, 0.40, size=(50,))
    >>> pos_rate.sum()
    21

    Let's set up a test function, so we can test the probability of
    finding a difference in frequency between the two groups. We'll use
    `scipy.stats.chisquare` to look for the difference in frequency between
    groups.

    >>> from scipy.stats import chisquare
    >>> test = lambda x: chisquare(np.array([x[i].sum() for i in
    ...     range(len(x))]))[1]

    Let's make sure that our two distributions are different.

    >>> print(round(test([pre_rate, pos_rate]), 3))
    0.003

    Since there are an even number of samples, and we don't have enough
    information to try controlling the data, we'll use
    `skbio.stats.power.subsample_power` to compare the two groups. If we had
    metadata about other risk factors, like a reproductive history, BMI,
    tobacco use, we might want to use
    `skbio.stats.power.subsample_paired_power`.
    We'll also use "ind" `draw_mode`, since there is no linkage between the
    two groups of samples.

    >>> from skbio.stats.power import subsample_power
    >>> pwr_est, counts = subsample_power(test=test,
    ...                                   samples=[pre_rate, pos_rate],
    ...                                   num_iter=100,
    ...                                   num_runs=5,
    ...                                   counts_interval=5)
    >>> counts
    array([ 5, 10, 15, 20, 25, 30, 35, 40, 45])
    >>> np.nanmean(pwr_est, axis=0) # doctest: +NORMALIZE_WHITESPACE
    array([ 0.056,  0.074,  0.226,  0.46 ,  0.61 ,  0.806,  0.952,  1.   ,
            1.   ])
    >>> counts[np.nanmean(pwr_est, axis=0) > 0.8].min()
    30

    So, we can estimate that we will see a significant difference in the
    presence of *G. vaginalis* in the stool of pre and post women with UTIs if
    we have at least 30 samples per group.

    If we wanted to test the relationship of a second candidate taxa which is
    more rare in the population, but may have a similar effect, based on
    available literature, we might also start by trying to identify 30
    samples per group where the second candidate taxa is present.

    Suppose, now, that we want to test that a secondary metabolite seen only in
    the presence of *G vaginalis* to see if it is also correlated with UTIs. We
    can model the abundance of the metabolite as a normal distribution.

    >>> met_pos = (np.random.randn(pre_rate.sum() + pos_rate.sum()) * 2000 +
    ...     2500)
    >>> met_pos[met_pos < 0] = 0
    >>> met_neg = met_neg = (np.random.randn(100 - (pre_rate.sum() +
    ...     pos_rate.sum())) * 2000 + 500)
    >>> met_neg[met_neg < 0] = 0

    Let's compare the populations with a kruskal-wallis test. Physically, there
    cannot be a negative concentration of a chemical, so we've set the lower
    bound at 0. This means that we can no longer assume our distribution is
    normal.

    >>> from scipy.stats import kruskal
    >>> def metabolite_test(x):
    ...     return kruskal(x[0], x[1])[1]
    >>> print(round(metabolite_test([met_pos, met_neg]), 3))
    0.005

    When we go to perform the statistical test on all the data, you might
    notice that there are twice as many samples from women with *G. vaginalis*
    than those without. It might make sense to account for this difference when
    we're testing power. So, we're going to set the `ratio` parameter, which
    lets us draw twice as many samples from women with *G. vaginalis*.

    >>> pwr_est2, counts2 = subsample_power(test=metabolite_test,
    ...                                     samples=[met_pos, met_neg],
    ...                                     counts_interval=5,
    ...                                     num_iter=100,
    ...                                     num_runs=5,
    ...                                     ratio=[2, 1])
    >>> counts2
    array([  5.,  10.,  15.,  20.,  25.,  30.])
    >>> np.nanmean(pwr_est2, axis=0)
    array([ 0.14 ,  0.272,  0.426,  0.646,  0.824,  0.996])
    >>> counts2[np.nanmean(pwr_est2, axis=0) > 0.8].min()
    25.0

    When we consider the number of samples per group needed in the power
    analysis, we need to look at the ratio. The analysis says that we need 25
    samples in the smallest group, in this case, the group of women without
    *G. vaginalis* and 50 samples from women with *G. vaginalis* to see a
    significant difference in the abundance of our secondary metabolite at 80%
    power.

    """
    # Checks the inputs
    ratio, num_p, sample_counts = _check_subsample_power_inputs(
        test=test,
        samples=samples,
        draw_mode=draw_mode,
        ratio=ratio,
        min_counts=min_counts,
        max_counts=max_counts,
        counts_interval=counts_interval,
    )

    # Prealocates the power array
    power = np.zeros((num_runs, len(sample_counts), num_p))

    # Calculates the power instances
    for id2, c in enumerate(sample_counts):
        count = np.round(c * ratio, 0).astype(int)
        for id1 in range(num_runs):
            ps = _compare_distributions(
                test=test,
                samples=samples,
                num_p=num_p,
                counts=count,
                num_iter=num_iter,
                mode=draw_mode,
            )
            power[id1, id2, :] = _calculate_power(ps, alpha_pwr)

    power = power.squeeze()

    return power, sample_counts


def subsample_paired_power(
    test,
    meta,
    cat,
    control_cats,
    order=None,
    strict_match=True,
    alpha_pwr=0.05,
    max_counts=50,
    counts_interval=10,
    min_counts=None,
    num_iter=500,
    num_runs=10,
):
    r"""Estimate power iteratively using samples with matching metadata.

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
        The power calculated for each subsample at each count. The array is
        `num_runs` rows, a length with the same number of elements as
        `sample_counts` and a depth equal to the number of p values returned by
        `test`. If `test` returns a float, the returned array will be
        two-dimensional instead of three.
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
    translocation in myeloid-lineage cells. You are able to culture two
    macrophage lineages (bone marrow derived phagocytes and
    peritoneally-derived macrophages). Due to unfortunate circumstances, your
    growth media must be acquired from multiple sources (lab, company A,
    company B). Also unfortunate, you must use labor-intensive low throughput
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

    Let's check that cytokine treatment has a significant effect across all
    the cells.

    >>> treatment_stat = [g for g in data.groupby('TREATMENT').groups.values()]
    >>> round(f(treatment_stat), 17)
    0.00193863362662502

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
    ...                                   counts_interval=5,
    ...                                   num_iter=25,
    ...                                   num_runs=5)
    >>> cnt
    array([  5.,  10.,  15.,  20.])
    >>> pwr.mean(0)
    array([ 0.24 ,  0.528,  0.68 ,  0.88 ])
    >>> pwr.std(0).round(3)
    array([ 0.088,  0.127,  0.168,  0.08 ])

    Estimating off the power curve, it looks like 20 cells per group may
    provide adequate power for this experiment, although the large variance
    in power might suggest extending the curves or increasing the number of
    samples per group.

    """
    # Handles the order argument
    if order is None:
        order = sorted(meta.groupby(cat).groups.keys())
    order = np.array(order)

    # Checks for the number of sampling pairs available
    meta_pairs, index = _identify_sample_groups(
        meta, cat, control_cats, order, strict_match
    )
    min_obs = min(
        [
            _get_min_size(meta, cat, control_cats, order, strict_match),
            np.floor(len(index) * 0.9),
        ]
    )
    sub_ids = _draw_paired_samples(meta_pairs, index, min_obs)

    ratio, num_p, sample_counts = _check_subsample_power_inputs(
        test=test,
        samples=sub_ids,
        draw_mode="matched",
        min_counts=min_counts,
        max_counts=max_counts,
        counts_interval=counts_interval,
    )

    # Prealocates the power array
    power = np.zeros((num_runs, len(sample_counts), num_p))

    # Calculates power instances
    for id2, c in enumerate(sample_counts):
        for id1 in range(num_runs):
            ps = np.zeros((num_p, num_iter))
            for id3 in range(num_iter):
                subs = _draw_paired_samples(meta_pairs, index, c)
                ps[:, id3] = test(subs)
            power[id1, id2, :] = _calculate_power(ps, alpha_pwr)

    power = power.squeeze()

    return power, sample_counts


def confidence_bound(vec, alpha=0.05, df=None, axis=None):
    r"""Calculate a confidence bound assuming a normal distribution.

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
        num_counts = vec_shape[axis] - np.isnan(vec).sum() / (
            vec_shape[0] * vec_shape[1]
        )

    # Gets the df if not supplied
    if df is None:
        df = num_counts - 1

    # Calculates the bound
    # In the conversion from scipy.stats.nanstd -> np.nanstd `ddof=1` had to be
    # added to match the scipy default of `bias=False`.
    bound = (
        np.nanstd(vec, axis=axis, ddof=1)
        / np.sqrt(num_counts - 1)
        * scipy.stats.t.ppf(1 - alpha / 2, df)
    )

    return bound


def paired_subsamples(meta, cat, control_cats, order=None, strict_match=True):
    r"""Draw a list of samples varied by `cat` and matched for `control_cats`.

    This function is designed to provide controlled samples, based on a
    metadata category. For example, one could control for age, sex, education
    level, and diet type while measuring exercise frequency.

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
        This determines how data is grouped using `control_cats`. If a sample
        within `meta` has an undefined value (`NaN`) for any of the columns in
        `control_cats`, the sample will not be considered as having a match and
        will be ignored when `strict_match` is True. If `strict_match` is
        False, missing values (NaN) in the `control_cats` can be considered
        matches.


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
    array(['BB', 'TS', 'CB']...)

    So, for this set of data, we can match TS, CB, and BB based on their age,
    sex, and antibiotic use. SW cannot be matched in either group because
    `strict_match` was true, and there is missing AGE data for this sample.

    """
    # Handles the order argument
    if order is None:
        order = sorted(meta.groupby(cat).groups.keys())
    order = np.array(order)

    # Checks the groups in the category
    min_obs = _get_min_size(meta, cat, control_cats, order, strict_match)

    # Identifies all possible subsamples
    meta_pairs, index = _identify_sample_groups(
        meta=meta,
        cat=cat,
        control_cats=control_cats,
        order=order,
        strict_match=strict_match,
    )

    # Draws paired ids
    ids = _draw_paired_samples(meta_pairs=meta_pairs, index=index, num_samps=min_obs)

    return ids


def _get_min_size(meta, cat, control_cats, order, strict_match):
    """Determine the smallest group represented."""
    if strict_match:
        all_cats = copy.deepcopy(control_cats)
        all_cats.append(cat)
        meta = meta[all_cats].dropna()

    return meta.groupby(cat).count().loc[order, control_cats[0]].min()


def _check_nans(x, switch=False):
    r"""Return False if x is a nan and True is x is a string or number."""
    if isinstance(x, str):
        return True
    elif isinstance(x, (float, int)):
        return not np.isnan(x)
    elif switch and isinstance(x, (list, tuple)) and np.nan in x:
        return False
    elif switch and isinstance(x, (list, tuple)):
        return True
    else:
        raise TypeError("input must be a string, float or a nan")


def _calculate_power(p_values, alpha=0.05):
    r"""Calculate statistical power empirically.

    Parameters
    ----------
    p_values : 1-D array
        A 1-D numpy array with the test results.

    alpha : float
        The critical value for the power calculation.

    Returns
    -------
    power : float
        The empirical power, or the fraction of observed p values below the
        critical value.

    """
    p_values = np.atleast_2d(p_values)

    w = (p_values < alpha).sum(axis=1) / p_values.shape[1]

    return w


def _compare_distributions(test, samples, num_p, counts=5, mode="ind", num_iter=100):
    r"""Compare two distribution arrays iteratively.

    Parameters
    ----------
    test : function
        The statistical test which accepts an array_like of sample ids
        (list of lists) and returns a p-value. This can be a one-dimensional
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
    mode : {"ind", "matched", "paired"}, optional
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
        The p-values for the subsampled tests. If `test` returned a single
        p value, p_values is a one-dimensional array. If `test` returned an
        array, `p_values` has dimensions `num_iter` x `num_p`

    Raises
    ------
    ValueError
        If mode is not "ind" or "matched".
    ValueError
        If the arrays in samples are not the same length in "matched" mode.
    ValueError
        If counts is a 1-D array and counts and samples are different lengths.

    """
    # Prealocates the pvalue matrix
    p_values = np.zeros((num_p, num_iter))

    # Determines the number of samples per group
    num_groups = len(samples)
    samp_lens = [len(sample) for sample in samples]

    if isinstance(counts, int):
        counts = np.array([counts] * num_groups)

    for idx in range(num_iter):
        if mode == "matched":
            pos = np.random.choice(np.arange(0, samp_lens[0]), counts[0], replace=False)
            subs = [sample[pos] for sample in samples]
        else:
            subs = [
                np.random.choice(np.array(pop), counts[i], replace=False)
                for i, pop in enumerate(samples)
            ]

        p_values[:, idx] = test(subs)

    if num_p == 1:
        p_values = p_values.squeeze()

    return p_values


def _check_subsample_power_inputs(
    test,
    samples,
    draw_mode="ind",
    ratio=None,
    max_counts=50,
    counts_interval=10,
    min_counts=None,
):
    r"""Make sure that everything is sane before power calculations.

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
    if draw_mode not in {"ind", "matched"}:
        raise ValueError('mode must be "matched" or "ind".')

    # Determines the minimum number of ids in a category
    id_counts = np.array([len(id_) for id_ in samples])
    num_ids = id_counts.min()
    # Determines the number of groups
    num_groups = len(samples)

    # Checks that "matched" mode is handled appropriately
    if draw_mode == "matched":
        for id_ in samples:
            if not len(id_) == num_ids:
                raise ValueError(
                    "Each vector in samples must be the same "
                    'length in "matched" draw_mode.'
                )

    # Checks the number of counts is appropriate
    if min_counts is None:
        min_counts = counts_interval
    if (max_counts - min_counts) < counts_interval:
        raise ValueError("No subsamples of the specified size can be drawn.")

    # Checks the ratio argument is sane
    if ratio is None or draw_mode == "matched":
        ratio = np.ones((num_groups))
    else:
        ratio = np.asarray(ratio)
    if not ratio.shape == (num_groups,):
        raise ValueError("There must be a ratio for each group.")

    ratio_counts = np.array([id_counts[i] / ratio[i] for i in range(num_groups)])
    largest = ratio_counts.min()

    # Determines the number of p values returned by the test
    p_return = test(samples)
    if isinstance(p_return, float):
        num_p = 1
    elif isinstance(p_return, np.ndarray) and len(p_return.shape) == 1:
        num_p = p_return.shape[0]
    else:
        raise TypeError("test must return a float or one-dimensional array.")

    # Calculates the same counts
    sample_counts = np.arange(min_counts, min(max_counts, largest), counts_interval)

    return ratio, num_p, sample_counts


def _identify_sample_groups(meta, cat, control_cats, order, strict_match):
    """Aggregate samples matches for `control_cats` that vary by `cat`.

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
    order : list
        The order of groups in the category. This can be used
        to limit the groups selected. For example, if there's a category with
        groups 'A', 'B' and 'C', and you only want to look at A vs B, `order`
        would be set to ['A', 'B'].
    ctrl_pos : int
        The location of the smallest group in `order`.
    strict_match: bool, optional
        This determines how data is grouped using `control_cats`. If a sample
        within `meta` has an undefined value (`NaN`) for any of the columns in
        `control_cats`, the sample will not be considered as having a match and
        will be ignored when `strict_match` is True. If `strict_match` is
        False, missing values (NaN) in the `control_cats` can be considered
        matches.

    Returns
    -------
    meta_pairs : dict
        Describes the categories matched for metadata. The
        `control_cat`-grouped samples are numbered, corresponding to the
        second list in `index`. The group is keyed to the list of sample arrays
        with the same length of `order`.
    index : list
        A list of numpy arrays describing the positions of samples to be drawn.
        The first array is an index array. The second gives an integer
        corresponding to the `control_cat`-group, and the third lists the
        position of the reference group sample in the list of samples.

    """
    # Sets up variables to be filled
    meta_pairs = {}
    index = []
    i1 = 0

    # Groups the data by the control groups
    ctrl_groups = meta.groupby(control_cats).groups
    # Identifies the samples that satisfy the control pairs. Keys are iterated
    # in sorted order so that results don't change with different dictionary
    # ordering.
    for g in sorted(ctrl_groups, key=lambda k: str(k)):
        ids = ctrl_groups[g]
        # If strict_match, Skips over data that has nans
        if not _check_nans(g, switch=True) and strict_match:
            continue
        # Draws the samples that are matched for control cats
        m_ids = meta.loc[ids].groupby(cat).groups
        # Checks if samples from the cat groups are represented in those
        # Samples
        id_vecs = [sorted(m_ids[o]) for o in order if o in m_ids]
        # If all groups are represented, the index and results are retained
        if len(id_vecs) == len(order):
            min_vec = np.array([len(v) for v in id_vecs])
            loc_vec = np.arange(0, min_vec.min())
            meta_pairs[i1] = id_vecs
            index.append(np.zeros(loc_vec.shape) + i1)
            i1 += 1
        # If the groups are not represented, an empty array gets passed
        else:
            index.append(np.array([]))

    # Converts index to a 1d array
    index = np.hstack(index)

    # If index is empty, sets up meta_paris with a no key.
    if not meta_pairs:
        meta_pairs["no"] = order

    return meta_pairs, index


def _draw_paired_samples(meta_pairs, index, num_samps):
    """Draw a random set of ids from a matched list.

    Parameters
    ----------
    meta_pairs : dict
        Describes the categories matched for metadata. The
        `control_cat`-grouped samples are numbered, corresponding to the
        second list in `index`. The group is keyed to the list of sample arrays
        with the same length of `order`.
    index : list
        A list of numpy arrays describing the positions of samples to be drawn.
        The first array is an index array. The second gives an integer
        corresponding to the `control_cat`-group, and the third lists the
        position of the reference group sample in the list of samples.
    num_samps : int
        The number of samples.

    Returns
    -------
    ids : list
        A set of randomly selected ids groups from each group.

    """
    # Handles an empty paired vector
    if "no" in meta_pairs:
        return [np.array([]) for o in meta_pairs["no"]]

    # Identifies the absolute positions of the control group being drawn
    set_pos = np.random.choice(index, int(num_samps), replace=False).astype(int)

    subs = []

    # Draws the other groups. Get a collection.Counter object for simplicity
    counter = collections.Counter(set_pos)

    # counter.items() order isn't guaranteed in python 3.6 and then the random
    # choice isn't reproducible between python version, even specifying seed;
    # so we access such values through sets.
    set_list = set(set_pos)

    # then, as stated by @RNAer, since we can't assure that items in sets are
    # ordered, we choose to order set_list before accessing values
    set_list = sorted(set_list)

    # now set_list is ordered and we can iterate over it to get counter obj
    for set_ in set_list:
        num_ = counter[set_]
        r2 = [np.random.choice(col, num_, replace=False) for col in meta_pairs[set_]]
        subs.append(r2)

    ids = [np.hstack(ids) for ids in zip(*subs)]

    return ids


def _calculate_power_curve(
    test, samples, sample_counts, ratio=None, mode="ind", num_iter=1000, alpha=0.05
):
    r"""Generate an empirical power curve for the samples.

    Parameters
    ----------
    test : function
        The statistical test which accepts a list of arrays of values and
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
    alpha : float, optional
        The significance level for the statistical test. Defaults to 0.05.

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
        raise ValueError("There must be a ratio for each group.")

    # Loops through the sample sizes
    for id2, s in enumerate(sample_counts):
        count = np.round(s * ratio, 0).astype(int)
        for id1, a in enumerate(alpha):
            ps = _compare_distributions(
                test=test,
                samples=samples,
                counts=count,
                num_p=1,
                num_iter=num_iter,
                mode=mode,
            )
            if vec:
                pwr[id2] = (_calculate_power(ps, a)).item()

            else:
                pwr[id1, id2] = _calculate_power(ps, a)

    return pwr
