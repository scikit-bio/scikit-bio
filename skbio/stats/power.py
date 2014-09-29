r"""
Empirical Power Estimation (:mod:`skbio.stats.power`)
=====================================================

.. currentmodule:: skbio.stats.power

The purpose of this module is to provide empirical, post-hoc power estimation
of normally and non-normally distributed data. It also provides support to
subsample data to faciliate this analysis.

The underlying principle is based on subsampling and monte carlo simulation.
Assume that there is some set of populations, :math:`K_{1}, K_{2}, ... K_{n}`
which have some property, :math:`\mu` such that :math:`\mu_{1} \neq \mu_{2}
\neq ... \neq \mu_{n}`. For each of the populations, a sample, S can be drawn,
with a parameter, :math:`x` where :math:`x \approx \mu` and for the samples,
we can use a test, f, to show that :math:`x_{1} \neq x_{2} \neq ...
\neq x_{n}`.

Since we known that :math:`\mu_{1} \neq \mu_{2} \neq ... \neq \mu_{n}`,
we know we should reject the null hypothesis. If we fail to reject the null
hypothesis, we have comitted a Type II error and our result is a false
negative. We can estimate the frequency of Type II errors at various sampling
depths by repeatedly subsampling the populations and observing how often we
see a false negative. If we repeat this several times for each subsampling
depth, and vary the depths we use, we can start to approximate a relationship
between the number of samples we use and the rate of false negatives, also
called the statistical power of the test.

To generate complete power curves from data which appears underpowered, the
`statsmodels.stats.power` package can be used to solve for an effect size. The
effect size can be used to extrapolate a power curve for the data.

The general format for functions in this module is to define a statistical
test function which will take a list of ids, and return a p value. The test is
then evaluated over a series of subsample sizes.

If metadata is avaliable, there are three ways we can approach selecting our
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

Examples
--------
Suppose we wanted to look at the power curve for two variables, `ind` and
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
if they're correlated. We'll use the `scipy.stats.pearsonr` function, which
returns the pearson correlation coefficient and the probability value
that the data is not correlated. The function takes two vectors as its
input.

>>> from scipy.stats import pearsonr
>>> f = lambda x: pearsonr(x[0], x[1])[1]

Now, let's use random sampling to estimate the power of our test on
the first distribution.

>>> samples = [ind, dep]
>>> print f(samples)
3.64594525966e-08

Because we have relatively equal sample sizes, we should try completely
random sampling, or set `mode` to "all". This is recommended when sample
population are of simillar size, giving each sampled contained in the
population a simillar probability of being drawn. If one sample is much larger
than the other, signifigant subsampling can help decrease some of the noise.
In `get_subsampled_power`, we can maintain a paired relationship between
samples by setting `draw_mode` to "matched".

>>> from skbio.stats.power import get_subsampled_power
>>> pwr_ests, counts = get_subsampled_power(mode="all",
...                                         test=f,
...                                         samples=samples,
...                                         min_counts=3,
...                                         max_counts=10,
...                                         counts_start=3,
...                                         counts_interval=1,
...                                         draw_mode="matched")
>>> print counts
[3 4 5 6 7 8 9]
>>> print pwr_ests
[[ 0.234  0.642  0.876  0.96   0.99   1.     1.   ]
 [ 0.242  0.654  0.848  0.946  0.998  1.     1.   ]
 [ 0.244  0.664  0.884  0.946  0.988  1.     1.   ]
 [ 0.248  0.666  0.866  0.948  0.986  1.     1.   ]
 [ 0.242  0.658  0.9    0.94   0.99   1.     1.   ]
 [ 0.242  0.638  0.874  0.952  0.992  1.     1.   ]
 [ 0.24   0.66   0.904  0.95   0.988  1.     1.   ]
 [ 0.232  0.64   0.912  0.972  0.988  1.     1.   ]
 [ 0.256  0.646  0.854  0.952  0.992  1.     1.   ]
 [ 0.216  0.646  0.882  0.962  0.998  1.     1.   ]]

The `pwr_est` can then be used to fit an effect size using
`statsmodel.stats.power` module or the results can be average and plotted.

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
import numpy as np
from scipy.stats import t, nanstd


def get_subsampled_power(mode, test, meta=None, cat=None, control_cats=None,
                         order=None, strict=True, samples=None, sub_size=None,
                         draw_mode='ind', scaling=5, alpha_pwr=0.05,
                         min_counts=20, max_counts=50, counts_interval=10,
                         counts_start=None, num_iter=500, num_runs=10):
    r"""Subsamples data to iterative calculate power

    Parameters
    ----------
    mode : {"all", "sig", "paired"}
        This designates the way subsamples will be selected.
        In "all" mode, samples from `samples` are randomly drawn restrictions
        on sampling. In "all" mode, the way samples are handled in power curve
        calculation can be modified by setting `draw_mode`.
        In `sig` mode, there must be a signifigant difference between the
        subsamples generated between the groups. The overall difference between
        the groups must be as large as `alpha_pwr` / `scaling`. The difference
        between the subsamples must be as large a `alpha_pwr`. The
        signifigantly different subsamples are then subsampled to generate a
        curve. The process is repeated over multiple subsamples. This method is
        recommended for datasets where there is a large size difference between
        heterogenous groups.
        In "paired" mode, observations are matched using `meta`, where the
        groups in `cat`, specified by `order`, are varied and the
        `control_cats` are held constant between paired observations. If more
        than one observations satisfies the sampling criteria, the matched
        observation will be randomly selected.
        "paired" mode uses the metadata in `meta` to vary `cat` while keeping
        the metadata fields in `control_cats` constant between paired samples.
    test : function
        The statistical test which accepts a list of arrays of values
        (sample ids or numeric values) and returns a p value.
    meta : pandas.dataframe
        Default is None. The metadata associated with the samples.
        Required for "paired" mode.
    cat : str
        Default is None. The metadata categories for comparison.
        Required for "paired" mode.
    control_cats : list
        Default is None. The metadata categories to be used as controls. For
        example, if you wanted to control age (`cat` = "AGE"), you might want
        to control for gender and health status (i.e. `control_cats` = ["SEX",
        "HEALTHY"]).
        Required for "paired" mode.
    order : list, optional
        Default is None. The order of groups in the category. This can be used
        to limit the groups selected. For example, if there's a category with
        groups 'A', 'B' and 'C', and you only want to look at A vs B, `order`
        would be set to ['A', 'B']. `order` is only used in "paired" mode.
    strict: bool, optional
        Default is True. If a missing value (nan) is encountered, the group
        will be skipped when `strict` is True.
    samples : array-like
        Default is None. `samples` can be a list of lists or a list of arrays
        where each sublist or row in the array corresponds to a sampled group.
        Required for "all" and "sig" mode.
    sub_size : {None, int}, optional
        The maximum number of samples to select from a group. If no value is
        provided, this will be the same as the size of the smallest group.
        Otherwise, this will be compared to the size of the smallest group, and
        which ever is lower will be used. This is only used in "all" and "sig"
        modes.
    draw_mode : {"ind", "matched"}, optional
        This value is can only be modified in "all" mode.
        Default is "ind". "matched" samples should be used when observations in
        samples have corresponding observations in other groups. For instance,
        this may be useful when working with regression data where
        :math:`x_{1}, x_{2}, ..., x_{n}` maps to :math:`y_{1}, y_{2}, ... ,
        y_{n}`.
        In "sig" mode,`draw_mode` will be set to "ind", since in "sig" mode
        should only be used when observations are of uneven sizes. In "paired"
        mode, this is set to "matched" because this retains the case/control
        relationship between samples.
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
        Default is 10. The difference between each subsampling count.
    counts_start : unsigned int, optional
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
        The power calculated for each subsample at each count.
    sample_counts : array
        The number of samples drawn at each power calculation.

    Raises
    ------
    ValueError
        If mode is "paired" and meta, cat or control_cats is None, a
        ValueError is raised.
    ValueError
        If mode is "all" or "sig" and samples is None, a ValueError is raised.
    RuntimeError
        There is a runtime error if there are fewer samples than the minimum
        count.
    RuntimeError
        If the `counts_interval` is greater than the difference between the
        sample start and the max value, the function raises a RuntimeError.

    Examples
    --------
    Let's say we wanted to look at the presence of a spectific genus,
    :math:`\textit{Gardnerella}`, is pre and post menopasual women.
    We've collected samples from 100 women: 50 in each group. Let's start by
    simulating the probability that women in each group have the vaginosis.
    We'll set a random seed in numpy to maintain consistent results.

    >>> import numpy as np
    >>> np.random.seed(25)
    >>> pre_rate = np.random.binomial(1, 0.75, size=(50,))
    >>> print pre_rate
    [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1 1 0 1 1 1 1 1 1 1 0 1 1 1 1 1
     1 0 0 1 1 1 1 0 0 1 1 1 1]
    >>> pos_rate = np.random.binomial(1, 0.25, size=(50,))
    >>> print pos_rate
    [0 1 0 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 0 0 0 0 0 0 1 1 0 1 0 0 0 1 1 0 0 0 0
     1 0 0 0 1 0 1 0 0 0 1 0 0]

    Let's set up a test function, now, so we can test the probability of
    finding a difference in frequency between the two groups. We'll use
    `scipy.stats.chisquare` to look for the difference in frequency between
    groups.

    >>> from scipy.stats import chisquare, nanmean
    >>> test = lambda x: chisquare(np.array([x[i].sum() for i in
    ...     xrange(len(x))]))[1]

    Let's make sure that our two distributions are different.

    >>> print round(test([pre_rate, pos_rate]), 5)
    9e-05

    Since there are an even number of samples, and we don't have enough
    information to try controlling the data, let's subsample for power using
    "all" mode. We'll also use "ind" draw mode, since there is no linkage
    between the two groups of samples.

    >>> from skbio.stats.power import get_subsampled_power
    >>> pwr_est, counts = get_subsampled_power(mode="all",
    ...                                        test=test,
    ...                                        samples=[pre_rate, pos_rate],
    ...                                        counts_interval=5)
    >>> print counts
    [ 5 10 15 20 25 30 35 40 45]
    >>> print nanmean(pwr_est, 0)
    [ 0.176   0.3376  0.6582  0.884   0.9796  0.9986  1.      1.      1.    ]

    So, we can estimate the difference between the two populations is powered
    at 80% with between 15 and 20 samples.
    """

    # Checks the mode arguments
    if mode == "paired":
        if meta is None or cat is None or control_cats is None:
            raise ValueError("paired mode requires a meta dataframe, a "
                             "cat to vary and a set of control_cats.")
        else:
            sub_ids = get_paired_subsamples(meta, cat, control_cats, order,
                                            strict)
            draw_mode = 'matched'
    elif mode == 'sig':
        if samples is None:
            raise ValueError("sig mode requires samples be defined.")
        else:
            sub_ids = get_significant_subsample([test], samples, sub_size,
                                                alpha_pwr, num_iter, scaling)
            draw_mode = 'ind'
    elif mode == "all":
        if samples is None:
            raise ValueError("all mode requires samples be defined.")
        else:
            sub_ids = samples
    else:
        raise ValueError('%s is not a supported mode. Modes are "all", "sig", '
                         'and "paired".' % mode)

    # Determines the minium number of ids in a category
    num_ids = np.array([len(id_) for id_ in sub_ids]).min()

    # Checks there are enough samples to subsample
    if num_ids <= min_counts:
        raise RuntimeError('There are not enough samples for subsampling.')

    # Calculates the effect size vector
    if counts_start is None:
        counts_start = counts_interval

    if (max_counts - counts_start) < counts_interval:
        raise RuntimeError("No subsamples of the specified size can be drawn.")

    sample_counts = np.arange(counts_start,
                              min(max_counts, num_ids),
                              counts_interval)

    # Prealocates the power array
    power = np.zeros((num_runs, len(sample_counts)))

    power[0, :] = _calculate_power_curve(test,
                                         sub_ids,
                                         sample_counts,
                                         mode=draw_mode,
                                         num_iter=num_iter,
                                         alpha=alpha_pwr)

    # Calculates the power instances
    for id1 in range(num_runs):
        # Gets the subsample
        if mode == "paired":
            sub_ids = get_paired_subsamples(meta, cat, control_cats, order,
                                            strict)

        elif mode == 'sig':
            sub_ids = get_significant_subsample([test], samples, sub_size,
                                                alpha_pwr, num_iter,
                                                scaling)
        else:
            sub_ids = samples
            # Calculates the power curve
            power[id1, :] = _calculate_power_curve(test,
                                                   sub_ids,
                                                   sample_counts,
                                                   num_iter=num_iter,
                                                   alpha=alpha_pwr,
                                                   mode=draw_mode)

    return power, sample_counts


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
    p_values : 1d array
        A 1d numpy array with the test results.

    alpha : float
        The critical value for the power calculation.

    Returns
    -------
    power : float
        The emperical power, or the fraction of observed p values below the
        critical value.

    """

    w = (p_values < float(alpha)).sum()/float(p_values.shape[0])

    return w


def _compare_distributions(test, samples, counts=5, mode="ind", num_iter=1000):
    r"""Compares two distribution arrays iteratively

    Parameters
    ----------
    test : function
        The statistical test which accepts an array-like of sample ids
        (list of lists) and returns a p-value.
    samples : list of arrays
        A list where each 1-d array represents a sample. If `mode` is
        "matched", there must be an equal number of observations in each
        sample.
    counts : unsigned int or 1d array, optional
        Default is 5. The number of samples to draw from each distribution.
        If this is a 1d array, the length must correspond to the number of
        samples. The function will not draw more observations than are in a
        sample. In "matched" `mode`, the same number of observations will be
        drawn from each group.
    mode : {"ind", "matched"}, optional
        Default is "ind". "matched" samples should be used when observations in
        samples have corresponding observations in other groups. For instance,
        this may be useful when working with regression data where
        :math:`x_{1}, x_{2}, ..., x_{n}` maps to :math:`y_{1}, y_{2}, ... ,
        y_{n}`.
    num_iter : int, optional
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
    RuntimeError
        If the arrays in samples are not the same length in "matched" mode.
    RuntimeError
        If counts is a 1d array and counts and samples are different lengths.

    """

    # Determines the number of groups
    num_groups = len(samples)

    # Checks the mode
    if mode not in {'ind', 'matched'}:
        raise ValueError('Supported sample modes are "ind" and "matched".')

    # Handles the number of samples for later instances
    if isinstance(counts, int):
        counts = np.array([counts]*num_groups)
    elif not len(counts) == num_groups:
        raise RuntimeError('If counts is a 1d array, there must be a count to '
                           'draw for each group.')

    # Checks the group length
    samp_lens = [len(sample) for sample in samples]
    # Checks the group length
    if mode == 'matched' and np.array([samp_lens[i] != samp_lens[i+1] for i in
                                       range(num_groups-1)]).all():
        raise RuntimeError('In "matched" mode, each sample must have the same'
                           ' number of observations.')
    if np.array([samp_lens[i] < counts[i] for i in range(num_groups)]).any():
        raise RuntimeError('You cannot choose more observations that exist '
                           'in a sample.')

    # Prealocates the pvalue matrix
    p_values = np.zeros((num_iter))

    for idx in range(num_iter):
        if mode == "matched":
            pos = np.random.choice(np.arange(0, samp_lens[0]), counts[0],
                                   replace=False)
            subs = [sample[pos] for sample in samples]
        else:
            subs = [np.random.choice(np.array(pop), counts[i], replace=False)
                    for i, pop in enumerate(samples)]

        p_values[idx] = test(subs)

    return p_values


def confidence_bound(vec, alpha=0.05, df=None, axis=None):
    r"""Calculates a confidence bound assuming a normal distribution

    Parameters
    ----------
    vec : array
        A 1d numpy array of the values to use in the bound calculation.
    alpha : {0.05, float}
        The critical value, used for the confidence bound calculation.
    df : {None, float}, optional
        The degrees of freedom associated with the distribution. If None is
        given, df is assumed to be the number elements in specified axis.
    axis : {None, unsigned int}, optional
        Default is None. The axis over which to take the devation.

    Return
    ------
    bound : float
        The confidence bound around the mean. The confidence interval is
        [mean - bound, mean + bound].

    """

    # Determines the number of non-nan counts
    vec_shape = vec.shape
    if axis is None and len(vec_shape) == 1:
        num_counts = vec_shape[0] - np.isnan(vec).sum()
        axis = None
    elif axis is None:
        num_counts = vec_shape[0] * vec_shape[1] - np.isnan(vec).sum()
    else:
        num_counts = vec_shape[axis] - np.isnan(vec).sum() / \
            (vec_shape[0] * vec_shape[1])

    # Gets the df if not supplied
    if df is None:
        df = num_counts - 1

    # Calculates the bound
    bound = nanstd(vec, axis=axis) / np.sqrt(num_counts - 1) * \
        t.ppf(1 - alpha / 2, df)

    return bound


def _calculate_power_curve(test, samples, sample_counts, ratio=None,
                           mode='ind', num_iter=1000, alpha=0.05):
    r"""Generates an empirical power curve for the samples.

    Parameters
    ----------
    test : function
        The statistical test which accepts an list of arrays of values and
        returns a p value.
    samples : array-like
        `samples` can be a list of lists or an array where each sublist or row
        in the array corresponds to a sampled group.
    sample_counts : 1d array
        A vector of the number of samples which should be sampled in each
        curve.
    mode : {"ind", "matched"}, optional
        Default is "ind". "matched" samples should be used when observations in
        samples have corresponding observations in other groups. For instance,
        this may be useful when working with regression data where
        :math:`x_{1}, x_{2}, ..., x_{n}` maps to :math:`y_{1}, y_{2}, ... ,
        y_{n}`.
    ratio : 1d array, optional
        Default is None. The fraction of the sample counts which should be
        assigned to each group. If this is a 1d array, it must be the same
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
    elif not ratio.shape == (num_groups,):
        raise ValueError('There must be a ratio for each group.')

    # Loops through the sample sizes
    for id2, s in enumerate(sample_counts):
        count = np.round(s*ratio, 0).astype(int)
        for id1, a in enumerate(alpha):
            ps = _compare_distributions(test=test,
                                        samples=samples,
                                        counts=count,
                                        num_iter=num_iter,
                                        mode=mode)
            if vec:
                pwr[id2] = _calculate_power(ps, a)
            else:
                pwr[id1, id2] = _calculate_power(ps, a)

    return pwr


def bootstrap_power_curve(test, samples, sample_counts, ratio=None,
                          alpha=0.05, mode='ind', num_iter=500, num_runs=10):
    r"""Repeatedly calculates the power curve for a specified alpha level

    Parameters
    ----------
    test : function
        The statistical test which accepts an array-like of sample ids
        (list of lists or list ) and returns a p-value.
    samples : array-like
        samples can be a list of lists or an array where each sublist or row in
        the array corresponds to a sampled group.
    sample_counts : 1d array
        A vector of the number of samples which should be sampled in each curve
    ratio : {None, 1d array}
        The fraction of the sample counts which should be assigned to each
        group. This must be a none-type object, or the same length as samples.
    alpha : float, optional
        The critical value for calculating power. The default is 0.05.
    mode : {"ind", "matched"}, optional
        Default is "ind". "matched" samples should be used when observations in
        samples have corresponding observations in other groups. For instance,
        this may be useful when working with regression data where
        :math:`x_{1}, x_{2}, ..., x_{n}` maps to :math:`y_{1}, y_{2}, ... ,
        y_{n}`.
    num_iter : unsigned int
        Default is 1000. The number of p-values to generate for each point
        on the curve.
    num_runs : unsigned int
        Default is 5. The number of times to calculate each curve.

    Returns
    -------
    p_mean : 1d array
        The mean p-values from the iterations.
    p_std : vector
        The variance in the p-values.

    Examples
    --------
    Suppose we have 100 samples randomly drawn from two normal distribitions,
    the first with mean 0 and standard devation 1, and the second with mean 3
    and standard deviation 1.5

    >>> import numpy as np
    >>> np.random.seed(20)
    >>> samples_1 = np.random.randn(100)
    >>> samples_2 = 1.5*np.random.randn(100) + 1

    We want to test the statistical power of a independent two sample t-test
    comparing the two populations. We can define a test function, f, to perform
    the comparison. The test function will take a list of value vectors and
    return a p value.

    >>> from scipy.stats import ttest_ind
    >>> f = lambda x: ttest_ind(x[0], x[1])[1]

    Now, we can determine the statitical power, or the probability that do not
    have a false negative given that we do not have a false positive, by
    varying a number of subsamples.

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
    alpha = np.ones((num_runs))*alpha

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
    test : function
        The statistical test which accepts a list of arrays of values
        (sample ids or numeric values) and returns a p value.
    samples : array-like
        Default is None. `samples` can be a list of lists or a list of arrays
        where each sublist or row in the array corresponds to a sampled group.
    sub_size : int
        Default is None. The maximum number of samples to select from a group.
        If no value is provided, this will be the same as the size of the
        smallest group. Otherwise, this will be compared to the size of the
        smallest group, and which ever is lower will be used.
    p_crit : {float, list}
        The critical p value or p values for the function.
    num_rounds : {500, int}
        The number of times the code should attempt to subsample the data
        before determining it has tried too many times and should quit.
    p_scaling : {5, int}
        A penalty scale on p_crit, so that the total distribution must be less
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

    Examples
    --------
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
    check_size = np.array([len(i) for i in samples])
    if sub_size is None:
        sub_size = check_size.min()
    else:
        sub_size = min([sub_size, check_size.min()])

    # Checks the critical value is the same length as the tests
    if isinstance(p_crit, float):
        p_crit = p_crit*np.ones((len(tests)))

    # Verifies testing is reasonable for the
    for idx, f in enumerate(tests):
        if f is not None and p_crit[idx]/p_scaling < f(samples):
            tests[idx] = None
    # Checks the functions are defined
    if (tests == np.array([None]*len(tests))).all():
        raise RuntimeError('There is no test defined')

    # Loops through to get a signfigant difference
    for i in range(num_rounds+1):
        # Subsamples the larger dataset
        sub_samps = []
        for ids in samples:
            sub_samps.append(np.random.choice(ids, size=sub_size,
                                              replace=False))

        # Tests the subsample
        test_res = np.ones((len(tests)))
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
    meta : pandas.dataframe
        The metadata associated with the samples.
    cat : str, list
        The metadata category (or a list of categories) for comparison.
    control_cats : list
        The metadata categories to be used as controls. For example, if you
        wanted to control age (`cat` = "AGE"), you might want to control for
        gender and health status (i.e. `control_cats` = ["SEX", "HEALTHY"])
    order : list, optional
        Default is None. The order of groups in the category. This can be used
        to limit the groups selected. For example, if there's a category with
        groups 'A', 'B' and 'C', and you only want to look at A vs B, `order`
        would be set to ['A', 'B'].
    strict: bool
        Default is True. If a missing value (nan) is encountered, the group
        will be skipped when `strict` is True.

    Returns
    -------
    ids : array
        a set of arrays which satisfy the criteria. These are not grouped by
        `cat`. An empty array indicates there are no sample ids which satisfy
        the requirements.

    Examples
    --------
    If we have a mapping file for a set of random samples looking at housing,
    sex, age and antibiotic use.

    >>> import pandas as pd
    >>> import numpy as np
    >>> meta = {'SW': {'HOUSING': '2', 'SEX': 'M', 'AGE': np.nan, 'ABX': 'N'},
    ...         'TS': {'HOUSING': '2', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'},
    ...         'CB': {'HOUSING': '3', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'},
    ...         'BB': {'HOUSING': '1', 'SEX': 'M', 'AGE': '40s', 'ABX': 'Y'}}
    >>> meta = pd.DataFrame.from_dict(meta, orient="index")
    >>> meta
       ABX HOUSING  AGE SEX
    BB   Y       1  40s   M
    CB   Y       3  40s   M
    SW   N       2  NaN   M
    TS   Y       2  40s   M
    <BLANKLINE>
    [4 rows x 4 columns]

    Let's say we want to vary housing, controlling for sex, age, antibiotics
    and sex.

    >>> from skbio.stats.power import get_paired_subsamples
    >>> ids = get_paired_subsamples(meta, 'HOUSING', ['SEX', 'AGE', 'ABX'])
    >>> ids #doctest: +NORMALIZE_WHITESPACE
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

    ids = [np.array([])]*num_groups
    # Loops through samples in the experimental group to match for controls
    for check_group, ctrl_ids in viewitems(ctrl_group):
        # Checks the categories have been defined
        undefed_check = np.array([_check_strs(p) for p in check_group])
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
        num_draw = np.array(num_ids).min()
        # Draws samples from possible ids
        exp_ids = [np.random.choice(ctrl_ids, num_draw, replace=False)]
        exp_ids.extend([np.random.choice(id_, num_draw, replace=False)
                        for id_ in pos_ids])

        if len(exp_ids) == num_groups:
            for idx in range(num_groups):
                ids[idx] = np.hstack((ids[idx], exp_ids[idx]))

    return ids
