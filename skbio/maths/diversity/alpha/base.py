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
    """Validate and convert input to an acceptable counts vector type.

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

    This function can also be used to convert a counts vector to a vector of
    count frequencies (e.g., a vector containing the number of singletons,
    doubletons, tripletons, etc.).

    Notes
    -----
    Does not check that `result` is big enough to store new counts.

    """
    indices = _validate(indices)
    if result is None:
        max_val = indices.max()
        result = np.zeros(max_val + 1)
    for i in indices:
        result[i] += 1
    return result


def berger_parker_d(counts):
    """Calculate Berger-Parker dominance.

    Berger-Parker dominance is defined as the fraction of the sample that
    belongs to the most abundant species.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Berger-Parker dominance.

    Notes
    -----
    Berger-Parker dominance is defined in [1]_. The implementation here is
    based on the description given in the SDR-IV online manual [2]_.

    References
    ----------
    .. [1] Berger & Parker (1970). SDR-IV online help.
    .. [2] http://www.pisces-conservation.com/sdrhelp/index.html

    """
    counts = _validate(counts)
    return counts.max() / counts.sum()


def brillouin_d(counts):
    """Calculate Brillouin index of alpha diversity.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Brillouin index.

    Notes
    -----
    The implementation here is based on the description given in the SDR-IV
    online manual [1]_.

    References
    ----------
    .. [1] http://www.pisces-conservation.com/sdrhelp/index.html

    """
    counts = _validate(counts)
    nz = counts[counts.nonzero()]
    n = nz.sum()
    return (gammaln(n + 1) - gammaln(nz + 1).sum()) / n


def dominance(counts):
    """Calculate dominance.

    Dominance is defined as 1 - Simpson's index and ranges between 0 and 1.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Dominance.

    See Also
    --------
    simpson

    Notes
    -----
    The implementation here is based on the description given in [1]_.

    References
    ----------
    .. [1] http://folk.uio.no/ohammer/past/diversity.html

    """
    counts = _validate(counts)
    freqs = counts / counts.sum()
    return (freqs * freqs).sum()


def doubles(counts):
    """Calculate number of double occurrences (doubletons).

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    int
        Doubleton count.

    """
    counts = _validate(counts)
    return (counts == 2).sum()


def enspie(counts):
    """Calculate ENS_pie alpha diversity measure.

    ENS_pie is defined as ``1 / sum(pi ^ 2)`` with the sum occurring over all
    ``S`` species in the pool. ``pi`` is the proportion of the entire community
    that species ``i`` represents.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        ENS_pie alpha diversity measure.

    See Also
    --------
    dominance

    Notes
    -----
    ENS_pie is defined in [1]_. ENS_pie is equivalent to ``1 / dominance``.

    References
    ----------
    .. [1] Chase and Knight (2013). "Scale-dependent effect sizes of ecological
       drivers on biodiversity: why standardised sampling is not enough".
       Ecology Letters, Volume 16, Issue Supplement s1, pgs 17-26.

    """
    counts = _validate(counts)
    return 1 / dominance(counts)


def equitability(counts, base=2):
    """Calculate equitability (Shannon index corrected for number of species).

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Measure of equitability.

    See Also
    --------
    shannon

    Notes
    -----
    The implementation here is based on the description given in the SDR-IV
    online manual [1]_.

    References
    ----------
    .. [1] http://www.pisces-conservation.com/sdrhelp/index.html

    """
    counts = _validate(counts)
    numerator = shannon(counts, base)
    denominator = np.log(observed_species(counts)) / np.log(base)
    return numerator / denominator


def esty_ci(counts):
    """Calculate Esty's CI.

    Esty's CI is defined as ``n1/n  +/- z*sqrt(W)`` where ``n1`` is the number
    of species observed once in ``n`` samples; ``n`` is the sample size; and
    ``z`` is a constant that depends on the targeted confidence and based on
    the Normal distribution. For a 95% CI, z=1.959963985; ``n2`` is the
    number of species observed twice in ``n`` samples; ``W`` is
    ``[n1*(n-n1) + 2*n*n2] / (n**3)``.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    tuple
        Esty's confidence interval as ``(upper_bound, lower_bound)``.

    Notes
    -----
    Esty's CI is defined in [1]_. ``z`` is hardcoded for a 95% confidence
    interval.

    References
    ----------
    .. [1] Esty, W. W. (1983). "A normal limit law for a nonparametric
       estimator of the coverage of a random sample". Ann Statist 11: 905-912.

    """
    counts = _validate(counts)

    n1 = singles(counts)
    n2 = doubles(counts)
    n = counts.sum()
    z = 1.959963985
    W = (n1 * (n - n1) + 2 * n * n2) / (n ** 3)

    return n1 / n + z * np.sqrt(W), n1 / n - z * np.sqrt(W)


def fisher_alpha(counts, bounds=(1e-3, 1e12)):
    """Calculate Fisher's alpha.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.
    bounds : tuple, optional
        Bounds for Powell optimizer (may need to adjust bounds for some
        datasets).

    Returns
    -------
    double
        Fisher's alpha.

    Notes
    -----
    The implementation here is based on the description given in the SDR-IV
    online manual [1]_.

    References
    ----------
    .. [1] http://www.pisces-conservation.com/sdrhelp/index.html

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
    """Calculate Good's coverage of counts.

    Good's coverage estimator is defined as ``1 - (n1/N)`` where ``n1`` is
    the number of singletons and ``N`` is the number of individuals (sum of
    abundances for all species).

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Good's coverage estimator.

    """
    counts = _validate(counts)
    n1 = singles(counts)
    N = counts.sum()
    return 1 - (n1 / N)


def heip_e(counts):
    """Calculate Heip's evenness measure.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Heip's evenness measure.

    Notes
    -----
    The implementation here is based on the description in [1]_.

    References
    ----------
    .. [1] Heip, C. 1974. A new index measuring evenness. J. Mar. Biol. Ass.
       UK., 54, 555-557.

    """
    counts = _validate(counts)
    return (np.exp(shannon(counts, base=np.e) - 1) /
            (observed_species(counts) - 1))


def kempton_taylor_q(counts, lower_quantile=.25, upper_quantile=.75):
    """Calculate Kempton-Taylor Q index of alpha diversity.

    Estimates the slope of the cumulative abundance curve in the interquantile
    range. By default, uses lower and upper quartiles, rounding inwards.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Kempton-Taylor Q index of alpha diversity.

    Notes
    -----
    The index is defined in [1]_. The implementation here is based on the
    description given in the SDR-IV online manual [2]_.

    The implementation provided here differs slightly from the results given in
    Magurran 1998. Specifically, we have 14 in the numerator rather than 15.
    Magurran recommends counting half of the species with the same # counts as
    the point where the UQ falls and the point where the LQ falls, but the
    justification for this is unclear (e.g. if there were a very large #
    species that just overlapped one of the quantiles, the results would be
    considerably off). Leaving the calculation as-is for now, but consider
    changing.

    References
    ----------
    .. [1] Kempton, R. A. and Taylor, L. R. (1976) Models and statistics for
       species diversity. Nature, 262, 818-820.
    .. [2] http://www.pisces-conservation.com/sdrhelp/index.html

    """
    counts = _validate(counts)
    n = len(counts)
    lower = int(np.ceil(n * lower_quantile))
    upper = int(n * upper_quantile)
    sorted_counts = np.sort(counts)
    return (upper - lower) / np.log(sorted_counts[upper] /
                                    sorted_counts[lower])


def margalef(counts):
    """Calculate Margalef's richness index.

    Assumes log accumulation.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Margalef's richness index.

    Notes
    -----
    Based on the description in [1]_.

    References
    ----------
    .. [1] Magurran, A E 2004. Measuring biological diversity. Blackwell. pp.
       76-77.

    """
    counts = _validate(counts)
    return (observed_species(counts) - 1) / np.log(counts.sum())


def mcintosh_d(counts):
    """Calculate McIntosh dominance index D.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        McIntosh dominance index D.

    See Also
    --------
    mcintosh_e

    Notes
    -----
    The index was proposed in [1]_. The implementation here is based on the
    description given in the SDR-IV online manual [2]_.

    References
    ----------
    .. [1] McIntosh, R. P. 1967 An index of diversity and the relation of
       certain concepts to diversity. Ecology 48, 1115-1126.
    .. [2] http://www.pisces-conservation.com/sdrhelp/index.html

    """
    counts = _validate(counts)
    u = np.sqrt((counts * counts).sum())
    n = counts.sum()
    return (n - u) / (n - np.sqrt(n))


def mcintosh_e(counts):
    """Calculate McIntosh's evenness measure E.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        McIntosh evenness measure E.

    See Also
    --------
    mcintosh_d

    Notes
    -----
    The implementation here is based on the description given in [1]_, *NOT*
    the one in the SDR-IV online manual, which is wrong.

    References
    ----------
    .. [1] Heip & Engels 1974 p 560.

    """
    counts = _validate(counts)
    numerator = np.sqrt((counts * counts).sum())
    n = counts.sum()
    s = observed_species(counts)
    denominator = np.sqrt((n - s + 1) ** 2 + s - 1)
    return numerator / denominator


def menhinick(counts):
    """Calculate Menhinick's richness index.

    Assumes square-root accumulation.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Menhinick's richness index.

    Notes
    -----
    Based on the description in [1]_.

    References
    ----------
    .. [1] Magurran, A E 2004. Measuring biological diversity. Blackwell. pp.
       76-77.

    """
    counts = _validate(counts)
    return observed_species(counts) / np.sqrt(counts.sum())


def michaelis_menten_fit(counts, num_repeats=1, params_guess=None,
                         return_b=False):
    """Calculate Michaelis-Menten fit to rarefaction curve of observed species.

    The Michaelis-Menten equation is defined as ``S = (S_max*n) / (B+n)``,
    where ``n`` is the number of individuals and ``S`` is the number of
    species. This function estimates the ``S_max`` parameter, and also provides
    an option to return ``B``.

    The fit is made to datapoints where n = 1, 2, ..., counts.sum(), and ``S``
    is the number of species represented in a random sample of ``n``
    individuals.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.
    num_repeats : int, optional
        The number of times to perform rarefaction (subsampling without
        replacement) at each value of ``n``.
    params_guess : tuple, optional
        Initial guess of ``S_max`` and ``B``. Default is ``(100, 500)``.
    return_b : bool, optional
        If ``True``, return the estimate for both ``S_max`` and ``B``. The
        default is to just return ``S_max``.

    Returns
    -------
    S_max : double
        Estimate of the ``S_max`` parameter in the Michaelis-Menten equation.
    B : double
        If ``return_b`` is ``True``, will also return the estimate for ``B`` in
        the Michaelis-Menten equation.

    See Also
    --------
    skbio.maths.subsample

    Notes
    -----
    There is some controversy about how to do the fitting. The ML model given
    in [1]_ is based on the assumption that error is roughly proportional to
    magnitude of observation, reasonable for enzyme kinetics but not reasonable
    for rarefaction data. Here we just do a nonlinear curve fit for the
    parameters using least-squares.

    References
    ----------
    .. [1] Raaijmakers, J. G. W. 1987 Statistical analysis of the
       Michaelis-Menten equation. Biometrics 43, 793-803.

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
    """Calculate the number of distinct species.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    int
        Distinct species count.

    """
    counts = _validate(counts)
    return (counts != 0).sum()


def osd(counts):
    """Calculate observed species, singles, and doubles.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    osd : tuple
        Observed species, singles, and doubles.

    See Also
    --------
    observed_species
    singles
    doubles

    Notes
    -----
    This is a convenience function used by many of the other measures that rely
    on these three measures.

    """
    counts = _validate(counts)
    return observed_species(counts), singles(counts), doubles(counts)


def robbins(counts):
    """Calculate Robbins' estimator for the probability of unobserved outcomes.

    Robbins' estimator is defined as ``S / (n+1)`` where ``S`` is the number of
    singletons.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Robbins' estimate.

    Notes
    -----
    Robbins' estimator is defined in [1]_. The estimate computed here is for
    ``n-1`` counts, i.e. the x-axis is off by 1.

    References
    ----------
    .. [1] Robbins, H. E (1968). Ann. of Stats. Vol 36, pp. 256-257.

    """
    counts = _validate(counts)
    return singles(counts) / counts.sum()


def shannon(counts, base=2):
    """Calculate Shannon entropy of counts (H), default in bits.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.
    base : scalar, optional
        Logarithm base to use in the calculations.

    Returns
    -------
    double
        Shannon diversity index H.

    Notes
    -----
    The implementation here is based on the description given in the SDR-IV
    online manual [1]_, except that the default logarithm base used here is 2
    instead of e.

    References
    ----------
    .. [1] http://www.pisces-conservation.com/sdrhelp/index.html

    """
    counts = _validate(counts)
    freqs = counts / counts.sum()
    nonzero_freqs = freqs[freqs.nonzero()]
    return -(nonzero_freqs * np.log(nonzero_freqs)).sum() / np.log(base)


def simpson(counts):
    """Calculate Simpson's index.

    Simpson's index = 1 - dominance.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Simpson's index.

    See Also
    --------
    dominance

    Notes
    -----
    The implementation here is ``1 - dominance`` as described in [1]_. Other
    references (such as [2]_) define Simpson's index as ``1 / dominance``.

    References
    ----------
    .. [1] http://folk.uio.no/ohammer/past/diversity.html
    .. [2] http://www.pisces-conservation.com/sdrhelp/index.html

    """
    counts = _validate(counts)
    return 1 - dominance(counts)


def simpson_e(counts):
    """Calculate Simpson's evenness measure E.

    Simpson's E = ``(1 / dominance) / observed_species``.

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Simpson's evenness measure E.

    See Also
    --------
    dominance
    enspie
    simpson

    Notes
    -----
    The implementation here is based on the description given in [1]_.

    References
    ----------
    .. [1] http://www.tiem.utk.edu/~gross/bioed/bealsmodules/simpsonDI.html

    """
    counts = _validate(counts)
    return enspie(counts) / observed_species(counts)


def singles(counts):
    """Calculate number of single occurrences (singletons).

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    int
        Singleton count.

    """
    counts = _validate(counts)
    return (counts == 1).sum()


def strong(counts):
    """Calculate Strong's dominance index (Dw).

    Parameters
    ----------
    counts : (N,) array_like, int
        Vector of counts.

    Returns
    -------
    double
        Strong's dominance index (Dw).

    Notes
    -----
    Strong's dominance index is defined in [1]_. The implementation here is
    based on the description given in the SDR-IV online manual [2]_.

    References
    ----------
    .. [1] Strong, W. L., 2002 Assessing species abundance uneveness within and
       between plant communities. Community Ecology, 3, 237-246.
    .. [2] http://www.pisces-conservation.com/sdrhelp/index.html

    """
    counts = _validate(counts)
    n = counts.sum()
    s = observed_species(counts)
    i = np.arange(1, len(counts) + 1)
    sorted_sum = np.sort(counts)[::-1].cumsum()
    return (sorted_sum / n - (i / s)).max()
