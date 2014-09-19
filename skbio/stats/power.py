r"""Empirical Power Estimation (:mod: `skbio.stats.power`)
=========================================================

.. currentmodule:: skbio.stats.power

The purpose of this module is to provide empirical, post-hoc power estimation
of microbiome data. It also provides support to subsample data to faciliate
this analysis.

The underlying principle is based on subsampling and monte carlo simulation.
Assume that there is some set of populations, :math: `K_{1}, K_{2}, ... K_{n}`
which have some property, u such that :math: `\mu_{1} \neq \mu_{2} \neq ...
\neq \mu_{n}`. For each of the populations, a sample, S can be drawn, with a
parameter, x where :math: `x \approx \mu` and for the samples, we can use a
test, f, to show that :math: `x_{1} \neq x_{2} \neq ... \neq x_{n}`.

Since we known that :math: `\mu_{1} \neq \mu_{2} \neq ... \neq \mu_{n}`,
we know we should reject the null hypothesis. If we fail to reject the null
hypothesis, we have comitted a Type II error and our result is a False
negative. We can estimate the frequency of Type II errors at various sampling
depths by repeatedly subsampling the populations and observing how often we
see a False negative. If we repeat this several times for each subsampling
depth, and vary the depths we use, we can start to approximate a relationship
between the number of samples we use and the rate of false negatives, also
called the statistical power of the test.

We can then use the rate of false negatives and use the `statsmodels.power`
package to solve for an effect size. This can be used to extrapolate a power
curve for the data.

The general format for functions in this module is to define a statistical
test function which will take a list of ids, and return a p value. The test is
then evaluated over a series of subsample sizes.

With microbiome data, there are three ways we can approach selecting our
sample. We may choose to simply draw n observations at random from the two
underlying samples. Alternatively, we can draw subsamples which are
significantly different. Finally, we can try to match samples based on a set
of control categories.

Functions
---------

.. autosummary::
    :toctree: generated/

    get_subsampled_power
    confidence_bound
    bootstrap_power_curve
    get_significant_subsample
    get_paired_subsamples

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
from numpy import (array, zeros, ones, round as nround, hstack, isnan,
                   sqrt, arange)
from numpy.random import choice
from scipy.stats import t, nanstd


def get_subsampled_power(mode, test, meta=None, cat=None, control_cats=None,
                         order=None, strict=True, samples=None, sub_size=None,
                         scaling=5, alpha_pwr=0.05, min_counts=20,
                         max_counts=50, counts_interval=10, counts_start=None,
                         num_iter=500, num_runs=10):
    r"""Subsamples data to iterative calculate power

      Parameters
    ----------
    test : function
        the statistical test which accepts an array-like of sample ids
        (list of lists) and returns a p-value.
    meta : {None, dataframe}
        the metadata associated with the samples. Required for "PAIRED" mode.
    cat : {None, str}
        the metadata categories for comparison. Required for "PAIRED" mode.
    control_cats : {None, list}
        the metadata categories to be used as controls. For example, if you
        wanted to control age (`cat` = "AGE"), you might want to control for
        gender and health status (i.e. `control_cats` = ["SEX", "HEALTHY"]).
        Required for "PAIRED" mode.
    order : {None, list}, optional
        Default is None. The order of groups in the category. This can be used
        to limit the groups selected. For example, if there's a category with
        groups 'A', 'B' and 'C', and you only want to look at A vs B, `order`
        would be set to ['A', 'B'].
    strict: bool, optional
        Default is True. If a missing value (nan) is encountered, the group
        will be skipped when 'strict' is True.
    samples : {None, array-like}
        samples can be a list of lists or an array where each sublist or row in
        the array corresponds to a sampled group. Required for "ALL" and "SIG"
        mode.
    sub_size : {None, int}, optional
        the maximum number of samples to select from a group. If no value is
        provided, this will be the same as the size of the smallest group.
        Otherwise, this will be compared to the size of the smallest group, and
        which ever is lower will be used.
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
    counts_start : {None, unsigned int}, optional
        Defualt is None. How many samples should be drawn for the smallest
        subsample. If this is None, the `counts_interval` will be used.
    num_iter : unsigned int
        Default is 1000. The number of p-values to generate for each point
        on the curve.
    num_runs : unsigned int
        Default is 10. The number of times to calculate each curve.

    Returns
    -------
    power : array
        power calculated for each subsample at each count
    sample_counts : array
        the number of samples drawn at each power calculation

    Raises
    ------
    ValueError
        if mode is PAIRED and meta, cat or control_cats is None
    ValueError
        if mode is ALL or SIG and samples is None
    RuntimeError
        if there are fewer samples than the minimum count
    RuntimeError
        if the `counts_interval` is greater than the difference between the
        sample start and the max value.

    Examples
    --------
    Suppose we wanted to look at the power curve for two varaibles, `ind` and
    `dep`, using completely random subsampling. To control for the pseudo
    random number generation, we will use a seed.

    >>> import numpy as np
    >>> np.random.seed(20)
    >>> ind = np.random.randint(0, 20, 15)
    >>> print ind
    [ 3 15  9 11  7  2  0  8 19 16  6  6 16  9  5]
    >>> dep = (3 * ind + 5 + np.random.randn(15)*5).round(3)
    >>> print dep
    [ 15.617  47.533  28.04   33.788  19.602  12.229   4.779  36.838  67.256
      55.032  22.157   7.051  58.601  38.664  18.783]


    Let's define a test that will draw a list of sample pairs and determine
    if they're correlated. We'll use the `pearsonr` function from scipy, which
    returns the pearson correlation coeffectient and the probability value
    that the data is not correlated. The function takes two vectors as its
    input.

    >>> from scipy.stats import pearsonr
    >>> f = lambda x: pearsonr(ind[x[0]], dep[x[0]])[1]

    Now, let's use random sampling to estimate the power of our test on
    the first distribution. Since our test picks pairs, the "samples" vector
    will just be a list of positions in each array.

    >>> samples = [np.arange(0, 15, 1)]
    >>> print f(samples)
    3.64594525966e-08

    Since we know the two samples are correlated overall, let's try using
    completely random subsampling. This is recommended when sample populations
    are of simillar size, giving each sampled contained in the population a
    simillar probability of being drawn. If one sample is much larger than the
    other, signifigant subsampling can help decrease some of the noise.

    >>> from skbio.stats.power import get_subsampled_power
    >>> pwr_ests, counts = get_subsampled_power(mode="ALL",
    ...                                         test=f,
    ...                                         samples=samples,
    ...                                         min_counts=3,
    ...                                         max_counts=10,
    ...                                         counts_start=3,
    ...                                         counts_interval=1)
    >>> print counts
    [3 4 5 6 7 8 9]
    >>> print pwr_ests
    [[ 0.22   0.652  0.89   0.958  0.992  1.     1.   ]
     [ 0.234  0.642  0.876  0.96   0.99   1.     1.   ]
     [ 0.242  0.654  0.848  0.946  0.998  1.     1.   ]
     [ 0.244  0.664  0.884  0.946  0.988  1.     1.   ]
     [ 0.248  0.666  0.866  0.948  0.986  1.     1.   ]
     [ 0.242  0.658  0.9    0.94   0.99   1.     1.   ]
     [ 0.242  0.638  0.874  0.952  0.992  1.     1.   ]
     [ 0.24   0.66   0.904  0.95   0.988  1.     1.   ]
     [ 0.232  0.64   0.912  0.972  0.988  1.     1.   ]
     [ 0.256  0.646  0.854  0.952  0.992  1.     1.   ]]

    The power_est can then be used to fit an effect_size using the power module
    if the statsmodel library, or can be average and plotted.

    """

    # Checks the mode arguments
    if mode == "PAIRED":
        if meta is None or cat is None or control_cats is None:
            raise ValueError("PAIRED mode requires a meta dataframe, a "
                             "cat to vary and a set of control_cats.")
        else:
            sub_ids = get_paired_subsamples(meta, cat, control_cats, order,
                                            strict)
    elif mode == 'SIG':
        if samples is None:
            raise ValueError("SIG mode requires samples be defined.")
        else:
            sub_ids = get_significant_subsample([test], samples, sub_size,
                                                alpha_pwr, num_iter, scaling)
    elif mode == "ALL":
        if samples is None:
            raise ValueError("ALL mode requires samples be defined.")
        else:
            sub_ids = samples
    else:
        raise ValueError('%s is not a supported mode. Modes are "ALL", "SIG", '
                         'and "PAIRED".' % mode)

    # Determines the minium number of ids in a category
    num_ids = array([len(id_) for id_ in sub_ids]).min()

    # Checks there are enough samples to subsample
    if num_ids <= min_counts:
        raise RuntimeError('There are not enough samples for subsampling.')

    # Calculates the effect size vector
    if counts_start is None:
        counts_start = counts_interval

    if (max_counts - counts_start) < counts_interval:
        raise RuntimeError("No subsamples of the specified size can be drawn.")

    sample_counts = arange(counts_start,
                           min(max_counts, num_ids),
                           counts_interval)

    # Prealocates the power array
    power = zeros((num_runs, len(sample_counts)))

    # Calculates the first power curve instance
    power[0, :] = _calculate_power_curve(test, sub_ids, sample_counts,
                                         num_iter=num_iter, alpha=alpha_pwr)

    # Calculates the power instances
    for id1 in arange(1, num_runs):
        # Gets the subsample
        if mode == "PAIRED":
            sub_ids = get_paired_subsamples(meta, cat, control_cats, order,
                                            strict)
        elif mode == 'SIGNIFICANT':
            sub_ids = get_significant_subsample([test], samples, sub_size,
                                                alpha_pwr, num_iter, scaling)
        else:
            sub_ids = samples
        # Calculates the power curve
        power[id1, :] = _calculate_power_curve(test, sub_ids, sample_counts,
                                               num_iter=num_iter,
                                               alpha=alpha_pwr)

    return power, sample_counts


def _check_strs(x):
    r"""Returns False if x is a nan and True is x is a string or number"""

    if isinstance(x, str):
        return True
    elif isnan(x):
        return False
    elif isinstance(x, (float, int)):
        return True
    else:
        raise TypeError('input must be a string, float or a nan')


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


def _compare_distributions(test, samples, counts=5, num_iter=1000):
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


def confidence_bound(vec, alpha=0.05, df=None, axis=None):
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


def _calculate_power_curve(test, samples, sample_counts, ratio=None,
                           num_iter=1000, alpha=0.05):
    r"""Generates an empirical power curve for the samples.

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
            ps = _compare_distributions(test, samples, count, num_iter)
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
    >>> np.random.seed(20)
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

    >>> from skbio.stats.power import bootstrap_power_curve
    >>> sample_counts = np.arange(5, 80, 5)
    >>> power_mean, power_bound = bootstrap_power_curve(f,
    ...                                                 [samples_1, samples_2],
    ...                                                 sample_counts)
    >>> print power_mean
    [ 0.2546  0.4736  0.6732  0.821   0.9084  0.9602  0.9846  0.9956  0.9996
      1.      1.      1.      1.      1.      1.    ]
    >>> print power_bound.round(3)
    [ 0.011  0.012  0.014  0.015  0.01   0.008  0.004  0.002  0.001  0.     0.
      0.     0.     0.     0.   ]

    """

    # Corrects the alpha value into a matrix
    alpha = ones((num_runs))*alpha

    # Boot straps the power curve
    power = _calculate_power_curve(test,
                                   samples,
                                   sample_counts,
                                   ratio,
                                   num_iter,
                                   alpha)

    # Calculates two summary statitics
    power_mean = power.mean(0)
    power_bound = confidence_bound(power, alpha=alpha[0], axis=0)

    # Calculates summary statitics
    return power_mean, power_bound


def get_significant_subsample(tests, samples, sub_size=None, p_crit=0.05,
                              num_rounds=500, p_scaling=5):
    r"""Subsamples data to an even sample number with a signficiant difference

    This function is recommended for use when sample sizes are severely
    skewed. For example, comparing a sample with 10 observations to a sample
    with 100 observations, it's likely the observed range, and the observed
    varaince of the larger sample will be greater. To control for a
    difference in sample size and limit uneven weighting due to this disparity,
    `get_significant_subsample` allows the user to select a group of subsamples
    which are signfigantly different at some critical value (a required
    assumption for this iterative post-hoc power analysis). The function will
    terminate if a signifigantly different subsample cannot be identified in
    a specified number of iterations, typically 500.

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

    Example
    -------
    Let's assume we have two samples drawn from two populations. The first
    sample has 25 observations, and the second has 200.

    >>> import numpy as np
    >>> np.random.seed(25)
    >>> y1 =  0.15*np.random.randint(20, 30, 25)
    >>> y1 = y1 + np.random.randint(25, 35, 25) + 5*np.random.randn(25)
    >>> print y1.mean().round(3)
    33.65
    >>> print np.array([y1.min(), y1.max()]).round(3)
    [ 17.954  48.704]
    >>> y2 = 0.15*np.random.randint(50, 60, 200) + 30 + 7*np.random.randn(200)
    >>> print y2.mean().round(3)
    37.841
    >>> print np.array([y2.min(), y2.max()]).round(3)
    [ 20.229  56.001]

    We can compare the two sample populations using a kruskal-wallis test.

    >>> from scipy.stats import kruskal
    >>> f = lambda x: kruskal(*x)[1]
    >>> print f([y1, y2]).round(3)
    0.001

    However, looking at the overlap in the sample range, it's possible that
    if we draw a random sample to compare the two distributions, we could get
    samples which intersect and do not reflect the true distibution and
    difference in the samples. So, insteaad, we use get_significant_subsample
    to get a subsample of y2 the same size as y1. Let's compare the
    signifigantly subsampled population to a random subsample.

    >>> from skbio.stats.power import (get_significant_subsample,
    ...                                bootstrap_power_curve)
    >>> vals = get_significant_subsample([f], [y1, y2])
    >>> all_mean, all_std = bootstrap_power_curve(f,
    ...                                           [y1, y2],
    ...                                           np.arange(5, 25, 5))
    >>> sub_mean, sub_std = bootstrap_power_curve(f,
    ...                                           vals,
    ...                                           np.arange(5, 25, 5))
    >>> print all_mean
    [ 0.1894  0.32    0.4562  0.6056]
    >>> print sub_mean
    [ 0.2212  0.4658  0.7008  0.965 ]

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
    for i in range(num_rounds+1):
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

    This function is designed to provide controlled samples, based on a
    metadata category. For example, one could control for age, sex, education
    level, and diet type while measuring exercise frequency. No outcome
    value is considered in this subsampling process.

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
    >>> import numpy as np
    >>> meta = {'SW': {'HOUSING': '2', 'SEX': 'M', 'AGE': np.nan, 'ABX': 'N'},
    ...         'TS': {'HOUSING': '2', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'},
    ...         'CB': {'HOUSING': '3', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'},
    ...         'BB': {'HOUSING': '1', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'}}
    >>> meta = pd.DataFrame.from_dict(meta, orient="index")


    Let's say we want to vary housing, controlling for sex, age, antibiotics
    and sex.

    >>> from skbio.stats.power import get_paired_subsamples
    >>> ids = get_paired_subsamples(meta, 'HOUSING', ['SEX', 'AGE', 'ABX'])
    >>> ids
    [array(['BB'], 
          dtype='|S2'), array(['TS'], 
          dtype='|S2'), array(['CB'], 
          dtype='|S2')]

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
    ctrl_match_groups = meta.groupby(control_cats).groups
    ctrl_group = meta.loc[cat_groups[ctrl_name]
                          ].groupby(list(control_cats)).groups

    ids = [array([])]*num_groups
    # Loops through samples in the experimental group to match for controls
    for check_group, ctrl_ids in viewitems(ctrl_group):
        # Checks the categories have been defined
        undefed_check = array([_check_strs(p) for p in check_group])
        if not undefed_check.all() and strict:
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
        num_draw = array(num_ids).min()
        # Draws samples from possible ids
        exp_ids = [choice(ctrl_ids, num_draw, replace=False)]
        exp_ids.extend([choice(id_, num_draw, replace=False) for id_ in
                        pos_ids])

        if len(exp_ids) == num_groups:
            for idx in range(num_groups):
                ids[idx] = hstack((ids[idx], exp_ids[idx]))

    return ids
