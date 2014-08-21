"""Empirical Power Estimation
==============================

The purpose of this module is to provide empirical, post-hoc power estimation
of microbiome data.

The power estimates here can use be used in conjunction with the statsmodel
power module to estimate effect size, based on the results.

The general format is to define a test for the data which takes a list of
values (i.e. sample ids) and returns a p-value.

"""

import random
from random import sample
from numpy import ndarray, array, zeros, ones, arange, round as nround
from scipy.stats import sem, t


def confidence_bound(vec, alpha=0.05, df=None):
    """Calculates a confidence bound assuming a normal distribution

    Parameters
    ----------
    vec : 1d array-like

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

    # Gets the df if not supplied
    if df is None:
        df = len(vec) - 1

    # Calculates the bound
    bound = sem(vec)*t.ppf((1-alpha/2), df)

    return bound


def calculate_power(p_values, alpha=0.05):
    """Calculates statical power empirically

    Parameters
    ----------
    p_values : 1d array-like
    alpha : float

    Returns
    -------
    power : float

    """
    if isinstance(p_values, list):
        p_values = array(p_values)

    w = (p_values < float(alpha)).sum()/float(len(p_values))

    return w


def compare_distributions(test, samples, num_samps=5, num_iter=1000):
    """Compares two distribution arrays iteratively

    Parameters
    ----------
    test : function
        the statistical test which accepts an array-like of sample ids
        (list of lists) and returns a p-value.

    samples : array-like
        samples can be a list of lists or an array where each sublist or row in
        the array corresponds to a sampled group.

    num_samps : {10, int, vector}
        the number of samples to draw from each distribution. If this is a
        vector, the length must correspond to the number of samples.

    num_iter : {1000, int}
        the number of p-values to generate for each point on the curve

    Returns
    -------
    p_values : array
        the bootstrapped p-values


    """
    # Determines the number of groups
    num_groups = len(samples)

    # Handles the number of samples for later instances
    if isinstance(num_samps, int):
        num_samps = array([num_samps]*num_groups)
    elif isinstance(num_samps, int) and len(num_samps) == 1:
        num_samps = array([num_samps*num_groups])
    elif isinstance(num_samps, ndarray) and len(num_samps) == 1:
        num_samps = num_samps*ones(num_groups)

    # Prealocates the pvalue matrix
    p_values = zeros((num_iter))

    for idx in range(num_iter):
        subs = [sample(pop, num_samps[i]) for i, pop in enumerate(samples)]
        p_values[idx] = test(subs)

    return p_values


def calculate_power_curve(test, samples, sample_counts, ratio=None,
                          num_iter=1000, alpha=0.05):
    """Generates an empirical power curve for the samples

    Parameters
    ----------
    test : function
        the statistical test which accepts an array-like of sample ids
        (list of lists) and returns a p-value.

    samples : array-like
        samples can be a list of lists or an array where each sublist or row in
        the array corresponds to a sampled group.

    sample_counts : {vector}
        A vector of the number of samples which should be sampled in each curve

    ratio : {None, vector}
        The fraction of the sample counts which should be assigned to each
        group. This must be a none-type object, or the same length as samples.

    num_iter : {1000, int}
        the number of p-values to generate for each point on the curve

    Returns
    -------
    p_values : array
        the bootstrapped p-values
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
                pwr[id2] = calculate_power(ps, a)
            else:
                pwr[id1, id2] = calculate_power(ps, a)

    return pwr


def bootstrap_power_curve(test, samples, sample_counts, ratio=None,
                          alpha=0.05, num_iter=500, num_runs=10):
    """
    Repeatedly calculates the power curve for a specified alpha level

    Parameters
    ----------
    test : function
        the statistical test which accepts an array-like of sample ids
        (list of lists) and returns a p-value.

    samples : array-like
        samples can be a list of lists or an array where each sublist or row in
        the array corresponds to a sampled group.

    sample_counts : {vector}
        A vector of the number of samples which should be sampled in each curve

    ratio : {None, vector}
        The fraction of the sample counts which should be assigned to each
        group. This must be a none-type object, or the same length as samples.

    num_iter : {1000, int}
        the number of p-values to generate for each point on the curve

    num_runs : {5, int}
        the number of times to calculate each curve

    Returns
    -------
    p_mean : vector
        the mean p-values from the iterations

    p_std : vector
        the variance in the p-values
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
    power_bound = confidence_bound(power, alpha=alpha[0])

    # Calculates summary statitics
    return power_mean, power_bound
