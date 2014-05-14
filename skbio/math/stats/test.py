#!/usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division

import numpy as np
from scipy.stats import spearmanr, kruskal, mannwhitneyu, kendalltau

from skbio.math.stats.special import MACHEP, ndtri
from skbio.math.stats.distribution import (chi_high, zprob, f_high, t_high,
                                           t_low, tprob)

np.seterr(divide='raise')


class ZeroExpectedError(ValueError):

    """Class for handling tests where an expected value was zero."""
    pass


def G_2_by_2(a, b, c, d, williams=1, directional=1):
    """G test for independence in a 2 x 2 table.

    Usage: G, prob = G_2_by_2(a, b, c, d, willliams, directional)

    Cells are in the order:

        a b
        c d

    a, b, c, and d can be int, float, or long.
    williams is a boolean stating whether to do the Williams correction.
    directional is a boolean stating whether the test is 1-tailed.

    Briefly, computes sum(f ln f) for cells - sum(f ln f) for
    rows and columns + f ln f for the table.

    Always has 1 degree of freedom

    To generalize the test to r x c, use the same protocol:
    2*(cells - rows/cols + table), then with (r-1)(c-1) df.

    Note that G is always positive: to get a directional test,
    the appropriate ratio (e.g. a/b > c/d) must be tested
    as a separate procedure. Find the probability for the
    observed G, and then either halve or halve and subtract from
    one depending on whether the directional prediction was
    upheld.

    The default test is now one-tailed (Rob Knight 4/21/03).

    See Sokal & Rohlf (1995), ch. 17. Specifically, see box 17.6 (p731).
    """
    cells = [a, b, c, d]
    n = np.sum(cells)
    # return 0 if table was empty
    if not n:
        return (0, 1)
    # raise error if any counts were negative
    if min(cells) < 0:
        raise ValueError(
            "G_2_by_2 got negative cell counts(s): must all be >= 0.")

    G = 0
    # Add x ln x for items, adding zero for items whose counts are zero
    for i in filter(None, cells):
        G += i * np.log(i)
    # Find totals for rows and cols
    ab = a + b
    cd = c + d
    ac = a + c
    bd = b + d
    rows_cols = [ab, cd, ac, bd]
    # exit if we are missing a row or column entirely: result counts as
    # never significant
    if min(rows_cols) == 0:
        return (0, 1)
    # Subtract x ln x for rows and cols
    for i in filter(None, rows_cols):
        G -= i * np.log(i)
    # Add x ln x for table
    G += n * np.log(n)
    # Result needs to be multiplied by 2
    G *= 2

    # apply Williams correction
    if williams:
        q = 1 + \
            ((((n / ab) + (n / cd)) - 1) * (((n / ac) + (n / bd)) - 1)) / \
            (6 * n)
        G /= q

    p = chi_high(max(G, 0), 1)

    # find which tail we were in if the test was directional
    if directional:
        is_high = ((b == 0) or (d != 0 and (a / b > c / d)))
        p = tail(p, is_high)
        if not is_high:
            G = -1 * G
    return G, p


def safe_sum_p_log_p(a, base=None):
    """Calculates p * log(p) safely for an array that may contain zeros."""
    flat = np.ravel(a)
    nz = np.take(flat, np.nonzero(flat)[0])
    logs = np.log(nz)
    if base:
        logs /= np.log(base)
    return np.sum(nz * logs, 0)


def g_stat(data):
    """Calculate the G statistic for data.

    Parameters
    ----------
    data : iterable of 1-D array_like
        Iterable of 1D array_like, each array is 1D with any length. Each array
        represents the observed frequencies of a given OTU in one of the sample
        classes.

    Returns
    -------
    G : float
        The G statistic that is additive between all groups.

    Notes
    -----
    For discussion read [1]_ pg. 695-699. The G test
    is normally applied to data when you have only one observation of any given
    sample class (e.g. you observe 90 wildtype and 30 mutants). In microbial
    ecology it is normal to have multiple samples which contain a given feature
    where those samples share a metadata class (eg. you observe OTUX at certain
    frequencies in 12 samples, 6 of which are treatment samples, 6 of which are
    control samples). To reconcile these approaches this function averages the
    frequency of the given feature (OTU) across all samples in the metadata
    class (e.g. in the 6 treatment samples, the value for OTUX is averaged, and
    this forms the average frequency which represents all treatment samples in
    aggregate). This means that this version of the G stat cannot detect sample
    heterogeneity as a replicated goodness of fit test would be able to. In
    addition, this function assumes the extrinsic hypothesis is that the
    mean frequency in all the samples groups is the same.

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    """
    # G = 2*sum(f_i*ln(f_i/f_i_hat)) over all i phenotypes/sample classes
    # calculate the total number of observations under the consideration that
    # multiple observations in a given group are averaged.
    avgs = np.array([np.array(i).mean() for i in data])
    n = avgs.sum()
    a = len(data)  # a is number of phenotypes or sample classes
    obs_freqs = avgs  # f_i values
    exp_freqs = np.zeros(a) + (n / a)  # f_i_hat vals
    G = 2. * (obs_freqs * np.log(obs_freqs / exp_freqs)).sum()
    return G


def g_fit(data, williams=True):
    """Calculate G statistic and compare to one tailed chi-squared distribution.

    Parameters
    ----------
    data : list of arrays
        List of arrays, each array is 1D with any length. Each array
        represents the observed frequencies of a given OTU in one of the sample
        classes.
    williams : boolean
        Whether or not to apply the Williams correction before comparing to the
        chi-squared distribution.

    Returns
    -------
    G : float
        The G statistic that is additive between all groups.
    pval : float
        The pvalue associated with the given G statistic.

    Notes
    -----
    For discussion read [1]_. This function
    compares the calculated G statistic (with Williams correction by default)
    to the chi-squared distribution with the appropriate number of degrees of
    freedom. If the data do not pass sanity checks for basic assumptions then
    this function will return nans.

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    """
    # sanity checks on the data to return nans if conditions are not met
    if not all([(i >= 0).all() for i in data]):
        # g_fit: data contains negative values. G test would be
        # undefined. Ignoring this OTU.
        return np.nan, np.nan
    if not all([i.sum() > 0 for i in data]):
        # g_fit: data contains sample group with zero only values. This
        # means that the given OTU was never observed in this sample class
        # The G test fails in this case because we would be forced to
        # take log(0). Ignoring this OTU.
        return np.nan, np.nan
    G = g_stat(data)
    a = len(data)  # a is number of phenotypes or sample classes
    if williams:
        # calculate the total number of observations under the consideration
        # that multiple observations in a given group are averaged.
        n = sum([arr.mean() for arr in data])  # total observations
        G = williams_correction(n, a, G)
    return G, chi_high(G, a - 1)  # a-1 degrees of freedom bc of sum constraint


def williams_correction(n, a, G):
    """Return the Williams corrected G statistic for G goodness of fit test.

    For discussion read [1]_ pg 698-699.

    Parameters
    ----------
    n : int
        Sum of observed frequencies.
    a : int
        Number of groups that are being compared.
    G : float
        Uncorrected G statistic

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    """
    # q = 1. + (a**2 - 1)/(6.*n*a - 6.*n) == 1. + (a+1.)/(6.*n)
    q = 1. + (a + 1.) / (6. * n)
    return G / q


def t_paired(a, b, tails=None, exp_diff=0):
    """Returns t and prob for TWO RELATED samples of scores a and b.

    From Sokal and Rohlf (1995), p. 354.
    Calculates the vector of differences and compares it to exp_diff
    using the 1-sample t test.

    Usage:   t, prob = t_paired(a, b, tails, exp_diff)

    t is a float; prob is a probability.
    a and b should be equal-length lists of paired observations (numbers).
    tails should be None (default), 'high', or 'low'.
    exp_diff should be the expected difference in means (a-b); 0 by default.
    """
    n = len(a)
    if n != len(b):
        raise ValueError('Unequal length lists in ttest_paired.')
    try:
        diffs = np.array(a) - np.array(b)
        return t_one_sample(diffs, popmean=exp_diff, tails=tails)
    except (ZeroDivisionError, ValueError, AttributeError, TypeError,
            FloatingPointError):
        return (None, None)


def t_one_sample(a, popmean=0, tails=None):
    """Returns t for ONE group of scores a, given a population mean.

    Usage:   t, prob = t_one_sample(a, popmean, tails)

    t is a float; prob is a probability.
    a should support Mean, StandardDeviation, and Count.
    popmean should be the expected mean; 0 by default.
    tails should be None (default), 'high', or 'low'.
    """
    try:
        n = len(a)
        t = (np.mean(a) - popmean) / (np.std(a, ddof=1) / np.sqrt(n))
    except (ZeroDivisionError, ValueError, AttributeError, TypeError,
            FloatingPointError):
        return None, None
    if np.isnan(t) or np.isinf(t):
        return None, None

    prob = t_tailed_prob(t, n - 1, tails)
    return t, prob


def t_two_sample(a, b, tails=None, exp_diff=0, none_on_zero_variance=True):
    """Returns t, prob for two INDEPENDENT samples of scores a, and b.

    From Sokal and Rohlf, p 223.

    Usage:   t, prob = t_two_sample(a,b, tails, exp_diff)

    t is a float; prob is a probability.
    a and b should be sequences of observations (numbers). Need not be equal
        lengths.
    tails should be None (default), 'high', or 'low'.
    exp_diff should be the expected difference in means (a-b); 0 by default.
    none_on_zero_variance: if True, will return (None,None) if both a and b
        have zero variance (e.g. a=[1,1,1] and b=[2,2,2]). If False, the
        following values will be returned:

            Two-tailed test (tails=None):
                a < b: (-inf,0.0)
                a > b: (+inf,0.0)

            One-tailed 'high':
                a < b: (-inf,1.0)
                a > b: (+inf,0.0)

            One-tailed 'low':
                a < b: (-inf,0.0)
                a > b: (+inf,1.0)

        If a and b both have no variance and have the same single value (e.g.
        a=[1,1,1] and b=[1,1,1]), (None,None) will always be returned.
    """
    if tails is not None and tails != 'high' and tails != 'low':
        raise ValueError("Invalid tail type '%s'. Must be either None, "
                         "'high', or 'low'." % tails)

    try:
        # see if we need to back off to the single-observation for single-item
        # groups
        n1 = len(a)
        if n1 < 2:
            t, prob = \
                t_one_observation(np.sum(a), b, tails, exp_diff,
                                  none_on_zero_variance=none_on_zero_variance)
            return t, prob

        n2 = len(b)
        if n2 < 2:
            t, prob = \
                t_one_observation(np.sum(b), a, reverse_tails(tails),
                                  exp_diff,
                                  none_on_zero_variance=none_on_zero_variance)

            # Negate the t-statistic because we swapped the order of the inputs
            # in the t_one_observation call, as well as tails.
            if t != 0:
                t = -1 * t

            return (t, prob)

        # otherwise, calculate things properly
        x1 = np.mean(a)
        x2 = np.mean(b)

        # pass ddof=1 to estimate the unbiased variance
        var1 = np.var(a, ddof=1)
        var2 = np.var(b, ddof=1)

        if var1 == 0 and var2 == 0:
            # Both lists do not vary.
            if x1 == x2 or none_on_zero_variance:
                result = (None, None)
            else:
                result = _t_test_no_variance(x1, x2, tails)
        else:
            # At least one list varies.
            df = n1 + n2 - 2
            svar = ((n1 - 1) * var1 + (n2 - 1) * var2) / df
            t = (x1 - x2 - exp_diff) / np.sqrt(svar * (1 / n1 + 1 / n2))

            if np.isnan(t) or np.isinf(t):
                result = (None, None)
            else:
                prob = t_tailed_prob(t, df, tails)
                result = (t, prob)
    except (ZeroDivisionError, ValueError, AttributeError, TypeError,
            FloatingPointError) as e:
        # invalidate if the sample sizes are wrong, the values
        # aren't numeric or aren't present, etc.
        result = (None, None)

    return result


def _t_test_no_variance(mean1, mean2, tails):
    """Handles case where two distributions have no variance."""
    if tails is not None and tails != 'high' and tails != 'low':
        raise ValueError("Invalid tail type '%s'. Must be either None, "
                         "'high', or 'low'." % tails)

    if tails is None:
        if mean1 < mean2:
            result = (float('-inf'), 0.0)
        else:
            result = (float('inf'), 0.0)
    elif tails == 'high':
        if mean1 < mean2:
            result = (float('-inf'), 1.0)
        else:
            result = (float('inf'), 0.0)
    else:
        if mean1 < mean2:
            result = (float('-inf'), 0.0)
        else:
            result = (float('inf'), 1.0)

    return result


def reverse_tails(tails):
    """Swaps high for low or vice versa, leaving other values alone."""
    if tails == 'high':
        return 'low'
    elif tails == 'low':
        return 'high'
    else:
        return tails


def mc_t_two_sample(x_items, y_items, tails=None, permutations=999,
                    exp_diff=0):
    """Performs a two-sample t-test with Monte Carlo permutations.

    x_items and y_items must be INDEPENDENT observations (sequences of
    numbers). They do not need to be of equal length.

    Returns the observed t statistic, the parametric p-value, a list of t
    statistics obtained through Monte Carlo permutations, and the nonparametric
    p-value obtained from the Monte Carlo permutations test.

    This code is partially based on Jeremy Widmann's
    qiime.make_distance_histograms.monte_carlo_group_distances code.

    Arguments:
        x_items - the first list of observations
        y_items - the second list of observations
        tails - if None (the default), a two-sided test is performed. 'high'
            or 'low' for one-tailed tests
        permutations - the number of permutations to use in calculating the
            nonparametric p-value. Must be a number greater than or equal to 0.
            If 0, the nonparametric test will not be performed. In this case,
            the list of t statistics obtained from permutations will be empty,
            and the nonparametric p-value will be None
        exp_diff - the expected difference in means (x_items - y_items)
    """
    if tails is not None and tails != 'high' and tails != 'low':
        raise ValueError("Invalid tail type '%s'. Must be either None, "
                         "'high', or 'low'." % tails)
    if permutations < 0:
        raise ValueError("Invalid number of permutations: %d. Must be greater "
                         "than or equal to zero." % permutations)

    if (len(x_items) == 1 and len(y_items) == 1) or \
       (len(x_items) < 1 or len(y_items) < 1):
        raise ValueError("At least one of the sequences of observations is "
                         "empty, or the sequences each contain only a single "
                         "observation. Cannot perform the t-test.")

    # Perform t-test using original observations.
    obs_t, param_p_val = t_two_sample(x_items, y_items, tails=tails,
                                      exp_diff=exp_diff,
                                      none_on_zero_variance=False)

    # Only perform the Monte Carlo test if we got a sane answer back from the
    # initial t-test and we have been specified permutations.
    nonparam_p_val = np.nan
    perm_t_stats = []
    if permutations > 0 and obs_t is not None and param_p_val is not None:
        # Permute observations between x_items and y_items the specified number
        # of times.
        perm_x_items, perm_y_items = _permute_observations(x_items, y_items,
                                                           permutations)
        perm_t_stats = [t_two_sample(perm_x_items[n], perm_y_items[n],
                                     tails=tails, exp_diff=exp_diff,
                                     none_on_zero_variance=False)[0]
                        for n in range(permutations)]

        # Compute nonparametric p-value based on the permuted t-test results.
        if tails is None:
            better = (np.absolute(np.array(perm_t_stats)) >=
                      np.absolute(obs_t)).sum()
        elif tails == 'low':
            better = (np.array(perm_t_stats) <= obs_t).sum()
        elif tails == 'high':
            better = (np.array(perm_t_stats) >= obs_t).sum()
        nonparam_p_val = (better + 1) / (permutations + 1)

    return obs_t, param_p_val, perm_t_stats, nonparam_p_val


def _permute_observations(x, y, num_perms):
    """Return num_perms pairs of permuted vectors x,y.

    Parameters
    ----------
    x : 1-D array-like
        Lists or arrays of values to be permuted.
    y : 1-D array-like
        Lists or arrays of values to be permuted.

    Returns
    -------
    xs
        Permuted vectors x
    ys
        Permuted vectors y
    """
    vals = np.hstack([np.array(x), np.array(y)])
    lenx = len(x)
    # sorting step is unnecessary for this code, but it ensure that test code
    # which relies on seeding the prng works (if we dont do this then different
    # observation orders in x and y for eg. the mc_t_two_sample test will fail
    # to produce the same results)
    vals.sort()
    inds = np.arange(vals.size)
    xs, ys = [], []
    for i in range(num_perms):
        np.random.shuffle(inds)
        xs.append(vals[inds[:lenx]])
        ys.append(vals[inds[lenx:]])
    return xs, ys


def t_one_observation(x, sample, tails=None, exp_diff=0,
                      none_on_zero_variance=True):
    """Returns t-test for significance of single observation versus a sample.

    Equation for 1-observation t (Sokal and Rohlf 1995 p 228):
    t = obs - mean - exp_diff / (var * sqrt((n+1)/n))
    df = n - 1

    none_on_zero_variance: see t_two_sample for details. If sample has no
        variance and its single value is the same as x (e.g. x=1 and
        sample=[1,1,1]), (None,None) will always be returned
    """
    try:
        sample_mean = np.mean(sample)
        sample_std = np.std(sample, ddof=1)

        if sample_std == 0:
            # The list does not vary.
            if sample_mean == x or none_on_zero_variance:
                result = (None, None)
            else:
                result = _t_test_no_variance(x, sample_mean, tails)
        else:
            # The list varies.
            n = len(sample)
            t = ((x - sample_mean - exp_diff) / sample_std / np.sqrt((n + 1) /
                                                                     n))
            prob = t_tailed_prob(t, n - 1, tails)
            result = (t, prob)
    except (ZeroDivisionError, ValueError, AttributeError, TypeError,
            FloatingPointError):
        result = (None, None)

    return result


def pearson(v1, v2):
    '''Pearson correlation using numpy.corrcoef.

    Parameters
    ----------
    v1 : array-like
        List or array of ints or floats to be correlated.
    v2 : array-like
        List or array of ints or floats to be correlated.

    Returns
    -------
    corrcoef : float
        Pearson correlation between the vectors.

    Raises
    ------
    ValueError
        If the vectors are not equally sized or if they are only a single
        element a ValueError will be returned.

    Examples
    --------
    >>> from skbio.math.stats.test import pearson
    >>> v1 = [.1, .2, .5, .3, .4]
    >>> v2 = [.9, .01, .5, .6, .7]
    >>> pearson(v1, v2)
    -0.052364331421504685
    '''
    v1, v2 = np.array(v1), np.array(v2)
    if not (v1.size == v2.size > 1):
        raise ValueError('One or more vectors isn\'t long enough to correlate '
                         ' or they have unequal lengths.')
    return np.corrcoef(v1, v2)[0][1]  # 2x2 symmetric unit matrix


def spearman(v1, v2):
    """Returns Spearman's rho.

    Parameters
    ----------
    v1 : array-like
        List or array of ints or floats to be correlated.
    v2 : array-like
        List or array of ints or floats to be correlated.

    Returns
    -------
    rho : float
        Spearman correlation between the vectors.

    Raises
    ------
    ValueError
        If the vectors are not equally sized or if they are only a single
        element a ValueError will be returned.

    See Also
    --------
    scipy.stats.spearmanr

    Notes
    -----
    This will always be a value between -1.0 and +1.0. v1 and v2 must
    be the same length, and cannot have fewer than 2 elements each. If one or
    both of the input vectors do not have any variation, the return value will
    be nan.
    """
    v1, v2 = np.array(v1), np.array(v2)
    if not (v1.size == v2.size > 1):
        raise ValueError('One or more vectors isn\'t long enough to correlate '
                         ' or they have unequal lengths.')
    return spearmanr(v1, v2)[0]  # return only the rho-value


def kendall(v1, v2):
    """Compute Kendall's Tau between v1 and v2 using scipy.stats.kendalltau

    Parameters
    ----------
    v1 : array-like
        List or array of ints or floats to be correlated.
    v2 : array-like
        List or array of ints or floats to be correlated.

    Returns
    -------
    rho : float
        Spearman correlation between the vectors.

    Raises
    ------
    ValueError
        If the vectors are not equally sized or if they are only a single
        element a ValueError will be returned.
    """
    v1, v2 = np.array(v1), np.array(v2)
    if not (v1.size == v2.size > 1):
        raise ValueError('One or more vectors isn\'t long enough to correlate '
                         ' or they have unequal lengths.')
    return kendalltau(v1, v2)[0]  # return only the tau correlation coeff


def kendall_pval(tau, n):
    '''Calculate the p-value for the passed tau and vector length n.'''
    test_stat = tau / ((2 * (2 * n + 5)) / float(9 * n * (n - 1))) ** .5
    return zprob(test_stat)


def assign_correlation_pval(corr, n, method, permutations=None,
                            perm_test_fn=None, v1=None, v2=None):
    """Assign pval to a correlation score with given method.

    This function will assign significance to the correlation score passed
    given the method that is passed. Some of the methods are appropriate only
    for certain types of data and there is no way for this test to determine
    the appropriateness, thus you must use this function only when with the
    proper prior knowledge. The 'parametric_t_distribution' method is described
    in [1]_ pg. 576, the 'fisher_z_transform' method
    is described on pg 576 and 577. The 'bootstrap' method calculates the given
    correlation permutations number of times using perm_test_fn.
    Also note, this does *not* take the place of FDR correction.

    Paramters
    ---------
    corr : float
        Correlation score from Kendall's Tau, Spearman's Rho, or Pearson.
    n : int
        Length of the vectors that were correlated.
    method : str
        One of ['parametric_t_distribution', 'fisher_z_transform',
        'bootstrapped', 'kendall'].
    permutations : int
        Number of permutations to use if bootstrapped selected.
    perm_test_fn : function
        Used to use to calculate correlation if permuation test desired.
    v1 : array-like or None
        List or array of ints or floats to be correlated. Passed if
        method='bootstrapped'.
    v2 : array-like or None
        List or array of ints or floats to be correlated. Passed if
        method='bootstrapped'

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    """
    if method == 'parametric_t_distribution':
        df = n - 2
        if df <= 1:
            raise ValueError("Must have more than 1 degree of freedom. "
                             "Can't Continue.")
        try:
            ts = corr * ((df / (1. - corr ** 2)) ** .5)
            return tprob(ts, df)  # two tailed test because H0 is corr=0
        except (ValueError, FloatingPointError, ZeroDivisionError):
            # something unpleasant happened, most likely r or rho where +- 1
            # which means the parametric p val should be 1 or 0 or nan
            return np.nan
    elif method == 'fisher_z_transform':
        # Sokal and Rohlf indicate that for n<50, the Fisher Z transform for
        # assigning correlation probabilities is not accurate. Currently no
        # check is in place
        z = fisher_z_transform(corr)
        # the z transform pval compares against a t distribution with inf
        # degrees of freedom which is equal to a z distribution.
        return z_transform_pval(z, n)
    elif method == 'bootstrapped':
        if any([i is None for i in [v1, v2, permutations, perm_test_fn]]):
            raise ValueError('You must specify vectors, permutation '
                             'function, and number of permutations to calc '
                             'bootstrapped pvalues. Cant continue.')
        r = np.empty(permutations)
        for i in range(permutations):
            r[i] = perm_test_fn(v1, np.random.permutation(v2))
        return (abs(r) >= abs(corr)).sum() / float(permutations)
    elif method == 'kendall':
        return kendall_pval(corr, n)
    else:
        raise ValueError("'%s' method is unknown." % method)


def correlation_t(x_items, y_items, method='pearson', tails=None,
                  permutations=999, confidence_level=0.95):
    """Computes the correlation between two vectors and its significance.

    Computes a parametric p-value by using Student's t-distribution with df=n-2
    to perform the test of significance, as well as a nonparametric p-value
    obtained by permuting one of the input vectors the specified number of
    times given by the permutations parameter. A confidence interval is also
    computed using Fisher's Z transform if the number of observations is
    greater than 3. Please see Sokal and Rohlf pp. 575-580 and pg. 598-601 for
    more details regarding these techniques.

    Warning: the parametric p-value is unreliable when the method is spearman
    and there are less than 11 observations in each vector.

    Returns the correlation coefficient (r or rho), the parametric p-value, a
    list of the r or rho values obtained from permuting the input, the
    nonparametric p-value, and a tuple for the confidence interval, with the
    first element being the lower bound of the confidence interval and the
    second element being the upper bound for the confidence interval. The
    confidence interval will be (None, None) if the number of observations is
    not greater than 3.

    x_items and y_items must be the same length, and cannot have fewer than 2
    elements each. If one or both of the input vectors do not have any
    variation, r or rho will be 0.0.

    Note: the parametric portion of this function is based on the correlation
    function in this module.

    Arguments:
        x_items - the first list of observations
        y_items - the second list of observations
        method - 'pearson' or 'spearman'
        tails - if None (the default), a two-sided test is performed. 'high'
            for a one-tailed test for positive association, or 'low' for a
            one-tailed test for negative association. This parameter affects
            both the parametric and nonparametric tests, but the confidence
            interval will always be two-sided
        permutations - the number of permutations to use in the nonparametric
            test. Must be a number greater than or equal to 0. If 0, the
            nonparametric test will not be performed. In this case, the list of
            correlation coefficients obtained from permutations will be empty,
            and the nonparametric p-value will be None
        confidence_level - the confidence level to use when constructing the
            confidence interval. Must be between 0 and 1 (exclusive)
    """
    # Perform some initial error checking.
    if method == 'pearson':
        corr_fn = pearson
    elif method == 'spearman':
        corr_fn = spearman
    else:
        raise ValueError("Invalid method '%s'. Must be either 'pearson' or "
                         "'spearman'." % method)
    if tails is not None and tails != 'high' and tails != 'low':
        raise ValueError("Invalid tail type '%s'. Must be either None, "
                         "'high', or 'low'." % tails)
    if permutations < 0:
        raise ValueError("Invalid number of permutations: %d. Must be greater "
                         "than or equal to zero." % permutations)
    if confidence_level <= 0 or confidence_level >= 1:
        raise ValueError("Invalid confidence level: %.4f. Must be between "
                         "zero and one." % confidence_level)

    # Calculate the correlation coefficient.
    corr_coeff = corr_fn(x_items, y_items)

    # Perform the parametric test first.
    x_items, y_items = np.array(x_items), np.array(y_items)
    n = len(x_items)
    df = n - 2
    if n < 3:
        parametric_p_val = 1
    else:
        try:
            t = corr_coeff / np.sqrt((1 - (corr_coeff * corr_coeff)) / df)
            parametric_p_val = t_tailed_prob(t, df, tails)
        except (ZeroDivisionError, FloatingPointError):
            # r/rho was presumably 1.
            parametric_p_val = 0

    # Perform the nonparametric test.
    permuted_corr_coeffs = []
    nonparametric_p_val = None
    better = 0
    for i in range(permutations):
        permuted_y_items = y_items[np.random.permutation(n)]
        permuted_corr_coeff = corr_fn(x_items, permuted_y_items)
        permuted_corr_coeffs.append(permuted_corr_coeff)

        if tails is None:
            if abs(round(permuted_corr_coeff, 10)) >= \
               abs(round(corr_coeff, 10)):
                better += 1
        elif tails == 'high':
            if round(permuted_corr_coeff, 10) >= round(corr_coeff, 10):
                better += 1
        elif tails == 'low':
            if round(permuted_corr_coeff, 10) <= round(corr_coeff, 10):
                better += 1
        else:
            # Not strictly necessary since this was checked above, but included
            # for safety in case the above check gets removed or messed up. We
            # don't want to return a p-value of 0 if someone passes in a bogus
            # tail type somehow.
            raise ValueError("Invalid tail type '%s'. Must be either None, "
                             "'high', or 'low'." % tails)
    if permutations > 0:
        nonparametric_p_val = (better + 1) / (permutations + 1)

    # Compute the confidence interval for corr_coeff using Fisher's Z
    # transform.
    z_crit = abs(ndtri((1 - confidence_level) / 2))
    ci_low, ci_high = None, None

    if n > 3:
        try:
            ci_low = np.tanh(np.arctanh(corr_coeff) - (z_crit /
                                                       np.sqrt(n - 3)))
            ci_high = np.tanh(np.arctanh(corr_coeff) + (z_crit /
                                                        np.sqrt(n - 3)))
        except (ZeroDivisionError, FloatingPointError):
            # r/rho was presumably 1 or -1. Match what R does in this case.
            ci_low, ci_high = corr_coeff, corr_coeff

    return (corr_coeff, parametric_p_val, permuted_corr_coeffs,
            nonparametric_p_val, (ci_low, ci_high))


def fisher(probs):
    """Uses Fisher's method to combine multiple tests of a hypothesis.

    -2 * SUM(ln(P)) gives chi-squared distribution with 2n degrees of freedom.
    """
    try:
        return chi_high(-2 * np.sum(np.log(probs)), 2 * len(probs))
    except OverflowError:
        return 0.0


def ANOVA_one_way(a):
    """Performs a one way analysis of variance

    a is a list of lists of observed values. Each list is the values
    within a category. The analysis must include 2 or more categories(lists).
    Each category of the list, and overall list, is converted to a numpy array.

    An F value is first calculated as the variance of the group means
    divided by the mean of the within-group variances.
    """
    group_means = []
    group_variances = []
    num_cases = 0  # total observations in all groups
    all_vals = []
    for i in a:
        num_cases += len(i)
        group_means.append(np.mean(i))
        group_variances.append(i.var(ddof=1) * (len(i) - 1))
        all_vals.extend(i)

    # Get within Group variances (denominator)
    dfd = num_cases - len(group_means)
    # need to add a check -- if the sum of the group variances is zero it will
    # error, but only if the between_Groups value is not zero
    within_Groups = np.sum(group_variances) / dfd
    if within_Groups == 0.:
        return np.nan, np.nan
    # Get between Group variances (numerator)
    all_vals = np.array(all_vals)
    grand_mean = all_vals.mean()
    between_Groups = 0
    for i in a:
        diff = i.mean() - grand_mean
        diff_sq = diff * diff
        x = diff_sq * len(i)
        between_Groups += x

    dfn = len(group_means) - 1
    between_Groups = between_Groups / dfn
    F = between_Groups / within_Groups
    return F, f_high(dfn, dfd, F)


def _average_rank(start_rank, end_rank):
    ave_rank = np.sum(range(start_rank, end_rank + 1)) / \
        (1 + end_rank - start_rank)
    return ave_rank


def _get_bootstrap_sample(x, y, num_reps):
    """yields num_reps random samples drawn with replacement from x and y"""
    combined = np.array(list(x) + list(y))
    total_obs = len(combined)
    num_x = len(x)
    for i in range(num_reps):
        # sampling with replacement
        indices = np.random.randint(0, total_obs, total_obs)
        sampled = combined.take(indices)
        # split into the two populations
        sampled_x = sampled[:num_x]
        sampled_y = sampled[num_x:]
        yield sampled_x, sampled_y


def mw_t(x, y, continuity=True, two_sided=True):
    '''Compute the Mann Whitney U statistic using scipy.stats.mannwhitneyu

    This wrapper controls whether the continuity correction will be applied
    and whether or not a two sided hypothesis is specified.

    Parameters
    ----------
    x : array-like
        List or array of numeric values to be tested.
    y : array-like
        List or array of numeric values to be tested.
    continuity : boolean
        Whether or not to use the continuity correction.
    two_sided: boolean
        Whether or not to use a two sided test. See Notes.

    Returns
    -------
    U stat : float
        The MWU U statistic.
    p-value : float
        The pvalue associated with the given U statistic assuming a normal
        probability distribution.

    See Also
    --------
    scipy.stats.mannwhitneyu

    Notes
    -----
    Two tails is appropriate because we do not know which of our groups has a
    higher mean, thus our alternate hypothesis is that the distributions from
    which the two samples come are not the same (FA!=FB) and we must account
    for E[FA] > E[FB] and E[FB] < E[FA]. See [1]_ pgs 427-431.

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    '''
    u, pval = mannwhitneyu(x, y, continuity)
    if two_sided:
        return u, 2. * pval
    else:
        return u, pval


def mw_boot(x, y, num_reps=999):
    """Bootstrapped version of Mann-Whitney-U test

    Parameters
    ----------
    x : array-like
        List or array of numeric values to be tested.
    y : array-like
        List or array of numeric values to be tested.
    num_reps : int
        Number of permutations tests to do.

    Returns
    -------
    observed_stat : float
        Value of the U statistic for the comparison of x and y.
    pval : float
        Number of times a U statistic as small or smaller than the observed U
        statistic was found.

    Notes
    -----
    The u statistic must be smaller than the observed u statistic to count as
    more extreme according to [1]_. Only a two tailed test is allowed through
    this function.

    Examples
    --------
    >>> from skbio.math.stats.test import mw_boot
    >>> x = [1.5, 4.6, 7.8, 10.2, 23.4]
    >>> y = [3.4, 10.1, 100.3, 45.6, 45.6, 78.9]
    >>> mw_boot(x, y, num_reps = 999)
    (6.0, 0.079)

    References
    ----------
    .. [1] http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/scipy.sta
       ts.mannwhitneyu.html
    """
    tol = MACHEP * 100
    combined = np.array(list(x) + list(y))
    np.random.seed(0)
    observed_stat, obs_p = mw_t(x, y)
    # from qiime.pycogent_backports.test import mw_test
    # observed_stat, obs_p = mw_test(x, y)
    total_obs = len(combined)
    num_x = len(x)
    u_stats_as_or_more_extreme = 0
    for sampled_x, sampled_y in _get_bootstrap_sample(x, y, num_reps):
        sample_stat, sample_p = mw_t(sampled_x, sampled_y)
        if sample_stat <= (observed_stat - tol):
            # the u statistic must be smaller than the observed u statistic to
            # count as more extreme. see [1]
            u_stats_as_or_more_extreme += 1
    return observed_stat, (u_stats_as_or_more_extreme + 1) / (num_reps + 1)


def kruskal_wallis(data):
    '''Calculate Kruskal Wallis U stat and pval using scipy.stats.kruskal

    Parameters
    ----------
    data : list of array-likes
        data is a nested list whose elements are arrays of float data. The
        different lists correpsond to the groups being tested.

    Returns
    -------
    U stat : float
        The Kruskal Wallis U statistic.
    pval : float
        The pvalue associated with the given U statistic assuming a chi-squared
        probability distribution.

    Examples
    --------
    >>> from skbio.math.stats.test import kruskal_wallis
    >>> data = [[1, 4.5, 67, 100, 2], [145, 100, 3, 14.5, -19], [2, 1.1, 5.5,
    ...         3.3, 16.7, 18, 100.3]]
    >>> rho, pval = kruskal_wallis(data)
    >>> print rho == 0.16848016848016789
    True
    >>> print pval == 0.91921054163678728
    True
    '''
    return kruskal(*data)


def permute_2d(m, p):
    """Performs 2D permutation of matrix m according to p."""
    return m[p][:, p]
    # unused below
    m_t = np.transpose(m)
    r_t = np.take(m_t, p, axis=0)
    return np.take(np.transpose(r_t), p, axis=0)


def mantel(m1, m2, n):
    """Compares two distance matrices. Reports P-value for correlation.

    The p-value is based on a two-sided test.

    WARNING: The two distance matrices must be symmetric, hollow distance
    matrices, as only the lower triangle (excluding the diagonal) will be used
    in the calculations (matching R's vegan::mantel function).

    This function is retained for backwards-compatibility. Please use
    mantel_t() for more control over how the test is performed.
    """
    return mantel_t(m1, m2, n)[0]


def mantel_t(m1, m2, n, alt="two sided",
             suppress_symmetry_and_hollowness_check=False):
    """Runs a Mantel test on two distance matrices.

    Returns the p-value, Mantel correlation statistic, and a list of Mantel
    correlation statistics for each permutation test.

    WARNING: The two distance matrices must be symmetric, hollow distance
    matrices, as only the lower triangle (excluding the diagonal) will be used
    in the calculations (matching R's vegan::mantel function).

    Arguments:
        m1  - the first distance matrix to use in the test (should be a numpy
            array or convertible to a numpy array)
        m2  - the second distance matrix to use in the test (should be a numpy
            array or convertible to a numpy array)
        n   - the number of permutations to test when calculating the p-value
        alt - the type of alternative hypothesis to test (can be either
            'two sided' for a two-sided test, 'greater' or 'less' for one-sided
            tests)
        suppress_symmetry_and_hollowness_check - by default, the input distance
            matrices will be checked for symmetry and hollowness. It is
            recommended to leave this check in place for safety, as the check
            is fairly fast. However, if you *know* you have symmetric and
            hollow distance matrices, you can disable this check for small
            performance gains on extremely large distance matrices
    """
    # Perform some sanity checks on our input.
    if alt not in ("two sided", "greater", "less"):
        raise ValueError("Unrecognized alternative hypothesis. Must be either "
                         "'two sided', 'greater', or 'less'.")
    m1, m2 = np.asarray(m1), np.asarray(m2)
    if m1.shape != m2.shape:
        raise ValueError("Both distance matrices must be the same size.")
    if n < 1:
        raise ValueError("The number of permutations must be greater than or "
                         "equal to one.")
    if not suppress_symmetry_and_hollowness_check:
        if not (is_symmetric_and_hollow(m1) and is_symmetric_and_hollow(m2)):
            raise ValueError("Both distance matrices must be symmetric and "
                             "hollow.")

    # Get a flattened list of lower-triangular matrix elements (excluding the
    # diagonal) in column-major order. Use these values to calculate the
    # correlation statistic.
    m1_flat, m2_flat = (_flatten_lower_triangle(m1),
                        _flatten_lower_triangle(m2))
    orig_stat = pearson(m1_flat, m2_flat)

    # Run our permutation tests so we can calculate a p-value for the test.
    size = len(m1)
    better = 0
    perm_stats = []
    for i in range(n):
        perm = permute_2d(m1, np.random.permutation(size))
        perm_flat = _flatten_lower_triangle(perm)
        r = pearson(perm_flat, m2_flat)

        if alt == 'two sided':
            if abs(r) >= abs(orig_stat):
                better += 1
        else:
            if ((alt == 'greater' and r >= orig_stat) or
                    (alt == 'less' and r <= orig_stat)):
                better += 1
        perm_stats.append(r)
    return (better + 1) / (n + 1), orig_stat, perm_stats


def t_tailed_prob(t, df, tails):
    """Return appropriate p-value for given t and df, depending on tails."""
    if tails == 'high':
        return t_high(t, df)
    elif tails == 'low':
        return t_low(t, df)
    else:
        return tprob(t, df)


def is_symmetric_and_hollow(matrix):
    """Return True if matrix is symmetric and hollow, otherwise False."""
    return (matrix.T == matrix).all() and (np.trace(matrix) == 0)


def _flatten_lower_triangle(matrix):
    """Returns a list containing the flattened lower triangle of the matrix.

    The returned list will contain the elements in column-major order. The
    diagonal will be excluded.

    Arguments:
        matrix - numpy array containing the matrix data
    """
    matrix = np.asarray(matrix)
    flattened = []
    for col_num in range(matrix.shape[1]):
        for row_num in range(matrix.shape[0]):
            if col_num < row_num:
                flattened.append(matrix[row_num][col_num])
    return flattened


def tail(prob, test):
    """If test is true, returns prob/2. Otherwise returns 1-(prob/2).
    """
    prob /= 2
    if test:
        return prob
    else:
        return 1 - prob


def bonferroni_correction(pvals):
    """Adjust pvalues for multiple tests using the Bonferroni method.

    In short: multiply all pvals by the number of comparisons.

    Parameters
    ----------
    pvals : list or array
        List or array of floats.

    Returns
    -------
    list of pvals
        Returns the list of pvals multiplied by their length. Pvals are
        still unsorted (i.e. order has not changed).

    See Also
    --------
    benjamini_hochberg_step_down

    Examples
    --------
    >>> from skbio.math.stats.test import bonferroni_correction
    >>> bonferroni_correction([0.1, 0.21, 0.5, 0.2, 0.6])
    array([ 0.5 ,  1.05,  2.5 ,  1.  ,  3.  ])
    """
    return np.array(pvals, dtype=float) * len(pvals)  # float conv: Nones->nans


def fdr_correction(pvals):
    """Adjust pvalues for multiple tests using the false discovery rate method.

    Parameters
    ----------
    pvals : list or array
        List or array of floats.

    Returns
    -------
    list of pvals
        Returns the list of pvals properly adjusted based on the FDR. Pvals are
        still unsorted (i.e. order has not changed).

    See Also
    --------
    benjamini_hochberg_step_down

    Notes
    -----
    In short: ranks the p-values in ascending order and multiplies each p-value
    by the number of comparisons divided by the rank of the p-value in the
    sorted list. Input is list of floats.  Does *not* assume pvals is sorted.

    Examples
    --------
    >>> from skbio.math.stats.test import fdr_correction
    >>> fdr_correction([.01, .2, .5, .1, .3])
    array([ 0.05      ,  0.33333333,  0.5       ,  0.25      ,  0.375     ])
    """
    tmp = np.array(pvals).astype(float)  # this converts Nones to nans
    return tmp * tmp.size / (1. + np.argsort(np.argsort(tmp)).astype(float))


def benjamini_hochberg_step_down(pvals):
    """Perform Benjamini and Hochberg's 1995 FDR step down procedure.

    Parameters
    ----------
    pvals : list or array
        List or array of floats.

    Returns
    -------
    list of pvals
        Returns the list of pvals properly adjusted based on the and then
        adjusted according to the BH rules.

    See Also
    --------
    fdr_correction

    Notes
    -----
    In short, computes the fdr adjusted pvals (ap_i's), and working from
    the largest to smallest, compare ap_i to ap_i-1. If ap_i < ap_i-1 set
    ap_i-1 equal to ap_i. Does *not* assume pvals is sorted. Described in [1]_.

    Examples
    --------
    >>> from skbio.math.stats.test import fdr_correction
    >>> benjamini_hochberg_step_down([0.1, 0.21, 0.5, 0.2, 0.6])
    array([ 0.35,  0.35,  0.6 ,  0.35,  0.6 ])

    References
    ----------
    .. [1] Controlling the False Discovery Rate: A Practical and Powerful
       Approach to Multiple Testing' Yoav Benjamini and Yosef Hochberg. Journal
       of the Royal Statistical Society. Series B (Methodological), Vol. 57,
       No. 1 (1995) 289-300.
    """
    tmp = fdr_correction(pvals)
    corrected_vals = np.empty(len(pvals))
    max_pval = 1.
    for i in np.argsort(pvals)[::-1]:
        if tmp[i] < max_pval:
            corrected_vals[i] = tmp[i]
            max_pval = tmp[i]
        else:
            corrected_vals[i] = max_pval
    return corrected_vals


def fisher_z_transform(r):
    """Calculate the Fisher Z transform of a correlation coefficient.

    Relies on formulation in [1_] pg 575.

    Parameters
    ----------
    r : float
        Correlation coefficient to transform.

    Returns
    -------
    z value of r

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    """
    if abs(r) >= 1:  # fisher z transform is undefined, have to return nan
        return np.nan
    return .5 * np.log((1. + r) / (1. - r))


def inverse_fisher_z_transform(z):
    """Calculate the inverse of the Fisher Z transform on a z value.

    Relies on formulation in [1_] pg 576.

    Parameters
    ----------
    z : float
        z value of a correlation coefficient that has undergone transformation.

    Returns
    -------
    r : float
        Rho or correlation coefficient that would produce given z score.

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    """
    return ((np.e ** (2 * z)) - 1.) / ((np.e ** (2 * z)) + 1.)


def z_transform_pval(z, n):
    '''Calculate two tailed probability of value as or more extreme than z.

    Relies on formulation in [1_] pg. 576.

    Parameters
    ----------
    z : float
        z-score
    n : int or float
        Number of samples that were used to generate the z-score.

    Returns
    -------
    zprob : float
        Probability of getting a zscore as or more extreme than the passed z
        given the total numger of samples that generated it (n).

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    '''
    if n <= 3:  # sample size must be greater than 3 otherwise this transform
        # isn't supported.
        return np.nan
    return zprob(z * ((n - 3) ** .5))


def fisher_population_correlation(corrcoefs, sample_sizes):
    """Calculate population rho, homogeneity from corrcoefs using Z transform.

    Parameters
    ----------
    corrcoefs : array-like
        A list or array of floats.
    sample_sizes : array-like
        A list or array of ints.

    Returns
    -------
    rho : float
        Combined rho for the population.
    h_val : float
        probablity that the rhos of the different samples were homogenous.

    Notes
    -----
    This function calculates the correlation of a population that is relying
    on multiple different studies that have calculated correlation coefficients
    of their own. The procedure is detailed in [1]_ pgs
    576 - 578. Pvals that are nans will be excluded by this function.

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117

    Examples
    --------
    >>> from skbio.math.stats.test import fisher_population_correlation
    >>> cc = [.4, .6, .7, .9]
    >>> ss = [10, 25, 14, 50]
    >>> fisher_population_correlation(cc, ss)
    (0.80559851846605035, 0.0029974748507499201)
    """
    tmp_rs = np.array(corrcoefs)
    tmp_ns = np.array(sample_sizes)
    # make checks for nans and exclude them as they will cause things to break
    rs = tmp_rs[~np.isnan(tmp_rs)]
    ns = tmp_ns[~np.isnan(tmp_rs)]
    if not (ns > 3).all():
        # not all samples have size > 3 which causes 0 varaince estimation.
        # thus we must return nan for pval and h_val
        return np.nan, np.nan
    if not len(ns) > 1:
        # only one sample, because of reduced degrees of freedom must have at
        # least two samples to calculate the homogeneity.
        return np.nan, np.nan
    if (rs >= 1.0).any():
        # a failure will occur in chi_high calculation where an non-terminating
        # loop will be initiated.
        raise ValueError('A correlation coefficient >= 1 was passed. This is '
                         'a non real valured correlation coefficient and it\'s'
                         ' Fisher Z transform cannot be computed.')
    # calculate zs
    zs = np.array([fisher_z_transform(float(i)) for i in rs])
    # calculate variance weighted z average = z_bar
    z_bar = (zs * (ns - 3)).sum() / float((ns - 3).sum())
    rho = inverse_fisher_z_transform(z_bar)
    # calculate homogeneity
    x_2 = ((ns - 3) * (zs - z_bar) ** 2).sum()
    h_val = chi_high(x_2, len(ns) - 1)
    return rho, h_val


def cscore(v1, v2):
    '''Calculate C-score between v1 and v2 according to Stone and Roberts 1990.

    Parameters
    ----------
    v1 : array-like
        List or array of numeric values to be tested.
    v2 : array-like
        List or array of numeric values to be tested.

    Returns
    -------
    cscore : float
        C-score between v1 and v2

    Notes
    -----
    This function calculates the C-score between equal length vectors v1 and v2
    according to the formulation given in [1]_.

    References
    ----------
    .. [1] Stone and Roberts. 1990, Oecologia 85:74-79
    '''
    v1_b = v1.astype(bool)
    v2_b = v2.astype(bool)
    sij = (v1_b * v2_b).sum()
    return (v1_b.sum() - sij) * (v2_b.sum() - sij)