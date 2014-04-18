#!/usr/bin/env python
from __future__ import division

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
from scipy.special import gammaln
from scipy.optimize import fmin_powell

from skbio.maths.subsample import subsample


def _validate(counts, suppress_cast=False):
    """Validate and convert input to counts vector.

    Note: may not always return a copy of `counts`!

    """
    counts = np.asarray(counts)

    if not suppress_cast:
        counts = counts.astype(int, casting='safe', copy=False)

    if counts.ndim != 1:
        raise ValueError("Only 1-D vectors are supported.")

    return counts


def _indices_to_counts(indices, result=None):
    """Converts vector of indices to counts of each index.

    WARNING: does not check that 'result' array is big enough to store new
    counts, suggest preallocating based on whole dataset if doing cumulative
    analysis.

    """
    indices = _validate(indices)
    if result is None:
        max_val = indices.max()
        result = np.zeros(max_val + 1)
    for i in indices:
        result[i] += 1
    return result


def berger_parker_d(counts):
    """Fraction of the sample that belongs to the most abundant species.

    References
    ----------
    .. [1] Berger & Parker 1970, by way of SDR-IV online help.

    """
    counts = _validate(counts)
    return counts.max() / counts.sum()


def brillouin_d(counts):
    """Calculate Brilloun index of alpha diversity.

    References
    ----------
    .. [1] Pielou 1975, by way of SDR-IV.

    """
    counts = _validate(counts)
    nz = counts[counts.nonzero()]
    n = nz.sum()
    return (gammaln(n + 1) - gammaln(nz + 1).sum()) / n


def dominance(counts):
    """Calculate probability that two species sampled are the same.

    Dominance = 1 - Simpson's index, sum of squares of probabilities.

    """
    counts = _validate(counts)
    freqs = counts / counts.sum()
    return (freqs * freqs).sum()


def doubles(counts):
    """Return count of double occurrences."""
    counts = _validate(counts)
    return (counts == 2).sum()


def enspie(counts):
    """Calculate ENS_pie alpha diversity measure.

    ENS_pie = 1 / sum(pi ^ 2) with the sum occurring over all ``S`` species in
    the pool. ``pi`` is the proportion of the entire community that species
    ``i`` represents.

    Notes
    -----
    For more information about ENS_pie, see [1]_.

    References
    ----------
    .. [1] "Scale-dependent effect sizes of ecological drivers on biodiversity:
       why standardised sampling is not enough". Chase and Knight. Ecology
       Letters, Volume 16, Issue Supplement s1, pgs 17-26 May 2013.

    """
    counts = _validate(counts)
    return 1 / dominance(counts)

# For backwards-compatibility with QIIME.
simpson_reciprocal = enspie


def equitability(counts, base=2):
    """Calculate Shannon index corrected for number of species, pure evenness.

    """
    counts = _validate(counts)
    numerator = shannon(counts, base)
    denominator = np.log(observed_species(counts)) / np.log(base)
    return numerator / denominator


def esty_ci(counts):
    """Esty's CI for (1-m).

    counts: Vector of counts (NOT the sample)

    Esty's CI is defined in
    Esty WW (1983) A Normal limit law for a nonparametric estimator of the
    coverage of a random sample. Ann Statist 11: 905-912.

    n1 / n  +/- z * square-root(W);

    where
    n1 = number of species observed once in n samples;
    n = sample size;
    z = a constant that depends on the targeted confidence and based on
        the Normal distribution. For a 95% CI, z=1.959963985;
    n2 = number of species observed twice in n samples;
    W = [ n1*(n - n1)  +  2*n*n2 ] / (n**3).

    Note: for other confidence levels we first need the appropriate z,
          Not yet hooked up to CLI.

    Returns: (upper bound, lower bound)
    """
    counts = _validate(counts)
    n1 = singles(counts)
    n2 = doubles(counts)
    n = counts.sum()
    z = 1.959963985
    W = (n1 * (n - n1) + 2 * n * n2) / (n ** 3)

    return n1 / n + z * np.sqrt(W), n1 / n - z * np.sqrt(W)


def fisher_alpha(counts, bounds=(1e-3, 1e12)):
    """Fisher's alpha: S = alpha ln(1+N/alpha) where S=species, N=individuals

    bounds are guess for Powell optimizer bracket search.

    WARNING: may need to adjust bounds for some datasets.
    """
    counts = _validate(counts)
    n = counts.sum()
    s = observed_species(counts)

    def f(alpha):
        if alpha >= bounds[0] and alpha <= bounds[1]:
            return (alpha * np.log(1 + (n / alpha)) - s) ** 2
        else:
            return np.inf

    alpha = fmin_powell(f, 1.0, disp=False)

    if f(alpha) > 1.0:
        raise RuntimeError("Optimizer failed to converge (error > 1.0), so "
                           "could not compute Fisher's alpha.")
    return alpha


def goods_coverage(counts):
    """Return Good's Coverage of counts.

    C = 1 - (n1/N)
    n1 = number of OTUs with abundance of 1
    N = number of individuals (sum of abundances for all OTUs)

    """
    counts = _validate(counts)
    n1 = singles(counts)
    N = counts.sum()
    return 1 - (n1 / N)


def heip_e(counts):
    """Calculate Heip's evenness measure.

    References
    ----------
    .. [1] Heip & Engels 1974.

    """
    counts = _validate(counts)
    return (np.exp(shannon(counts, base=np.e) - 1) /
            (observed_species(counts) - 1))


def kempton_taylor_q(counts, lower_quantile=.25, upper_quantile=.75):
    """Kempton-Taylor (1976) q index of alpha diversity, by way of SDR-IV.

    Estimates the slope of the cumulative abundance curve in the interquantile
    range. By default, uses lower and upper quartiles, rounding inwards.

    Note: this differs slightly from the results given in Magurran 1998.
    Specifically, we have 14 in the numerator rather than 15. Magurran
    recommends counting half of the species with the same # counts as the
    point where the UQ falls and the point where the LQ falls, but the
    justification for this is unclear (e.g. if there were a very large #
    species that just overlapped one of the quantiles, the results would
    be considerably off). Leaving the calculation as-is for now, but consider
    changing.
    """
    counts = _validate(counts)
    n = len(counts)
    lower = int(np.ceil(n * lower_quantile))
    upper = int(n * upper_quantile)
    sorted_counts = np.sort(counts)
    return (upper - lower) / np.log(sorted_counts[upper] /
                                    sorted_counts[lower])


def margalef(counts):
    """Margalef's index, assumes log accumulation.

    References
    ----------
    Magurran 2004, p 77.

    """
    counts = _validate(counts)
    return (observed_species(counts) - 1) / np.log(counts.sum())


def mcintosh_d(counts):
    """Calculate McIntosh index of alpha diversity.

    References
    ----------
    .. [1] McIntosh 1967, by way of SDR-IV.

    """
    counts = _validate(counts)
    u = np.sqrt((counts * counts).sum())
    n = counts.sum()
    return (n - u) / (n - np.sqrt(n))


def mcintosh_e(counts):
    """Calculate McIntosh's evenness measure.

    References
    ----------
    .. [1] Heip & Engels 1974 p 560 (wrong in SDR-IV).

    """
    counts = _validate(counts)
    numerator = np.sqrt((counts * counts).sum())
    n = counts.sum()
    s = observed_species(counts)
    denominator = np.sqrt((n - s + 1) ** 2 + s - 1)
    return numerator / denominator


def menhinick(counts):
    """Menhinick's index, assumes sqrt accumulation.

    References
    ----------
    .. [1] Magurran 2004, p 77.

    """
    counts = _validate(counts)
    return observed_species(counts) / np.sqrt(counts.sum())


def michaelis_menten_fit(counts, num_repeats=1, params_guess=None,
                         return_b=False):
    """Michaelis-Menten fit to rarefaction curve of observed species

    Note: there is some controversy about how to do the fitting. The ML model
    givem by Raaijmakers 1987 is based on the assumption that error is roughly
    proportional to magnitude of observation, reasonable for enzyme kinetics
    but not reasonable for rarefaction data. Here we just do a nonlinear
    curve fit for the parameters using least-squares.

    S = Smax*n/(B + n) . n: number of individuals, S: # of species
    returns Smax

    inputs:
    num_repeats: will perform rarefaction (subsampling without replacement)
    this many times at each value of n
    params_guess: intial guess of Smax, B (None => default)
    return_b: if True will return the estimate for Smax, B. Default: just Smax

    the fit is made to datapoints where n = 1,2,...counts.sum(),
    S = species represented in random sample of n individuals

    """
    counts = _validate(counts)

    if params_guess is None:
        params_guess = np.array([100, 500])

    # observed # of species vs # of individuals sampled, S vs n
    xvals = np.arange(1, counts.sum() + 1)
    ymtx = []
    for i in range(num_repeats):
        ymtx.append(np.array([observed_species(subsample(counts, n))
                              for n in xvals]))
    ymtx = np.asarray(ymtx)
    yvals = ymtx.mean(0)

    # fit to obs_sp = max_sp * num_idiv / (num_indiv + B)
    # return max_sp
    def fitfn(p, n):  # works with vectors of n, returns vector of S
        return p[0] * n / (p[1] + n)

    def errfn(p, n, y):  # vectors of actual vals y and number of individuals n
        return ((fitfn(p, n) - y) ** 2).sum()

    p1 = fmin_powell(errfn, params_guess, ftol=1e-5, args=(xvals, yvals),
                     disp=False)
    if return_b:
        return p1
    else:
        return p1[0]  # return only S_max, not the K_m (B) param


def observed_species(counts):
    """Calculate number of distinct species."""
    counts = _validate(counts)
    return (counts != 0).sum()


def osd(counts):
    """Calculate **o**bserved, **s**ingles and **d**oubles from counts."""
    counts = _validate(counts)
    return observed_species(counts), singles(counts), doubles(counts)


def robbins(counts):
    """Robbins 1968 estimator for Pr(unobserved) at n trials.

    probability_of_unobserved_colors = S/(n+1),

    Notes
    -----
    This is the estimate for ``(n-1)`` counts, i.e. x-axis is off by 1.

    References
    ----------
    .. [1] H. E. Robbins (1968, Ann. of Stats. Vol 36, pp. 256-257)
    (where s = singletons).

    """
    counts = _validate(counts)
    return singles(counts) / counts.sum()


def shannon(counts, base=2):
    """Calculate Shannon entropy of counts, default in bits."""
    counts = _validate(counts)
    freqs = counts / counts.sum()
    nonzero_freqs = freqs[freqs.nonzero()]
    return -(nonzero_freqs * np.log(nonzero_freqs)).sum() / np.log(base)


def simpson(counts):
    """Calculate Simpson's index.

    Simpson's index = 1 - dominance.

    """
    counts = _validate(counts)
    return 1 - dominance(counts)


def simpson_e(counts):
    """Calculate Simpson's evenness."""
    counts = _validate(counts)
    return enspie(counts) / observed_species(counts)


def singles(counts):
    """Return count of single occurrences."""
    counts = _validate(counts)
    return (counts == 1).sum()


def strong(counts):
    """Calculate Strong's 2002 dominance index, by way of SDR-IV."""
    counts = _validate(counts)
    n = counts.sum()
    s = observed_species(counts)
    i = np.arange(1, len(counts) + 1)
    sorted_sum = np.sort(counts)[::-1].cumsum()
    return (sorted_sum / n - (i / s)).max()
