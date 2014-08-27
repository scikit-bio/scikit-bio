"""Empirical Power Estimation
==============================

The purpose of this module is to provide empirical, post-hoc power estimation
of microbiome data.

The power estimates here can use be used in conjunction with the statsmodel
power module to estimate effect size, based on the results.

The general format is to define a test for the data which takes a list of
values (i.e. sample ids) and returns a p-value.

"""

from numpy import ndarray, array, zeros, ones, round as nround
from numpy.random import choice
from scipy.stats import sem, t


def confidence_bound(vec, alpha=0.05, df=None):
    """
    Calculates a confidence bound assuming a normal distribution

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
    """
    Calculates statical power empirically

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
    """
    Compares two distribution arrays iteratively

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
        subs = [choice(array(pop), num_samps[i], replace=False)
                for i, pop in enumerate(samples)]
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
