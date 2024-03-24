# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn

import numpy as np
from scipy.special import gammaln
from scipy.optimize import fmin_powell, minimize_scalar

from skbio.stats import subsample_counts
from skbio.diversity._util import _validate_counts_vector
from skbio.util._warning import _warn_deprecated


def berger_parker_d(counts):
    r"""Calculate Berger-Parker dominance index.

    Berger-Parker dominance index :math:`d` is defined as the fraction of the
    sample that belongs to the most abundant taxon:

    .. math::

       d = \frac{n_{max}}{N}

    where :math:`n_{max}` is the number of individuals in the most abundant
    taxon (or any of the most abundant taxa in the case of ties), and :math:`N`
    is the total number of individuals in the sample.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        Berger-Parker dominance index.

    Notes
    -----
    Berger-Parker dominance index was originally described in [1]_.

    References
    ----------
    .. [1] Berger, W. H., & Parker, F. L. (1970). Diversity of planktonic
       foraminifera in deep-sea sediments. Science, 168(3937), 1345-1347.

    """
    counts = _validate_counts_vector(counts)
    if (N := counts.sum()) == 0:
        return 0.0
    return counts.max() / N


def brillouin_d(counts):
    r"""Calculate Brillouin's diversity index.

    Brillouin's diversity index (:math:`H_B`) is defined as

    .. math::

       H_B = \frac{\ln N!-\sum^s_{i=1}{\ln n_i!}}{N}

    where :math:`N` is the total number of individuals in the sample, :math:`s`
    is the number of taxa, and :math:`n_i` is the number of individuals in the
    :math:`i^{\text{th}}` taxon.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        Brillouin's diversity index.

    Notes
    -----
    Brillouin's diversity index was originally described in [1]_.

    References
    ----------
    .. [1] Brillouin, L. (1956). Science and Information Theory. Academic
       Press. New York.

    """
    counts = _validate_counts_vector(counts)
    if (N := counts.sum()) == 0:
        return 0.0
    nz = counts[counts.nonzero()]
    return (gammaln(N + 1) - gammaln(nz + 1).sum()) / N


def dominance(counts):
    r"""Calculate Simpson's dominance index.

    Simpson's dominance index, a.k.a. Simpson's :math:`D`, measures the degree
    of concentration of taxon composition of a sample. It is defined as

    .. math::

       D = \sum{p_i^2}

    where :math:`p_i` is the proportion of the entire sample that taxon
    :math:`i` represents.

    Simpson's :math:`D` can be interpreted as the probability that two randomly
    selected individuals belong to the same taxon. It ranges between 0 and 1.

    Simpson's :math:`D` is sometimes referred to as "Simpson's index". It
    should be noted that :math:`D` is not a measure of community diversity. It
    is also important to distinguish :math:`D` from Simpson's diversity index
    (:math:`1 - D`) and Simpson's reciprocal index (:math:`1 / D`), both of
    which are measures of community diversity.

    Discrepancy exists among literature in using the term "Simpson index" and
    the denotion :math:`D`. It is therefore important to distinguish these
    metrics according to their mathematic definition.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        Simpson's dominance index.

    See Also
    --------
    simpson

    Notes
    -----
    Simpson's dominance index was originally described in [1]_.

    References
    ----------
    .. [1] Simpson, E. H. (1949). Measurement of diversity. Nature, 163(4148),
       688-688.

    """
    counts = _validate_counts_vector(counts)
    if (N := counts.sum()) == 0:
        return 0.0
    return ((counts / N) ** 2).sum()


def doubles(counts):
    """Calculate number of double-occurrence taxa (doubletons).

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    int
        Doubleton count.

    """
    counts = _validate_counts_vector(counts)
    return (counts == 2).sum()


def enspie(counts):
    r"""Calculate ENS_pie alpha diversity measure.

    ENS_pie is equivalent to ``1 / dominance``:

    .. math::

       ENS_{pie} = \frac{1}{\sum_{i=1}^s{p_i^2}}

    where :math:`s` is the number of taxa and :math:`p_i` is the proportion
    of the sample represented by taxon :math:`i`.

    Parameters
    ----------
    counts : 1-D array_like, int
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
    ENS_pie is defined in [1]_.

    References
    ----------
    .. [1] Chase and Knight (2013). "Scale-dependent effect sizes of ecological
       drivers on biodiversity: why standardised sampling is not enough".
       Ecology Letters, Volume 16, Issue Supplement s1, pgs 17-26.

    """
    counts = _validate_counts_vector(counts)
    return 1 / dominance(counts)


def esty_ci(counts):
    r"""Calculate Esty's confidence interval of Good's coverage estimator.

    Esty's confidence interval is defined as

    .. math::

       F_1/N \pm z\sqrt{W}

    where :math:`F_1` is the number of singleton taxa, :math:`N` is the
    total number of individuals, and :math:`z` is a constant that depends on
    the targeted confidence and based on the normal distribution.

    :math:`W` is defined as

    .. math::

       \frac{F_1(N-F_1)+2NF_2}{N^3}

    where :math:`F_2` is the number of doubleton taxa.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    tuple
        Esty's confidence interval as ``(lower_bound, upper_bound)``.

    Notes
    -----
    Esty's confidence interval was originally described in [1]_.

    :math:`z` is hardcoded for a 95% confidence interval.

    References
    ----------
    .. [1] Esty, W. W. (1983). "A normal limit law for a nonparametric
       estimator of the coverage of a random sample". Ann Statist 11: 905-912.

    """
    counts = _validate_counts_vector(counts)
    if (N := counts.sum()) == 0:
        return 0.0

    f1 = singles(counts)
    f2 = doubles(counts)
    z = 1.959963985
    W = (f1 * (N - f1) + 2 * N * f2) / (N**3)

    return f1 / N - z * np.sqrt(W), f1 / N + z * np.sqrt(W)


@np.errstate(invalid="ignore")
def fisher_alpha(counts):
    r"""Calculate Fisher's alpha, a metric of diversity.

    Fisher's alpha is estimated by solving the following equation for
    :math:`\alpha`:

    .. math::

       S=\alpha\ln(1+\frac{N}{\alpha})

    where :math:`S` is the number of taxa and :math:`N` is the total number
    of individuals in the sample.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        Fisher's alpha.

    Raises
    ------
    RuntimeError
        If the optimizer fails to solve for Fisher's alpha.

    Notes
    -----
    Fisher's alpha is defined in [1]_.

    There is no analytical solution to Fisher's alpha. However, one can use
    optimization techniques to obtain a numeric solution. This function calls
    SciPy's ``minimize_scalar`` to find alpha. It is deterministic. The result
    should be reasonably close to the true alpha.

    Alpha can become large when most taxa are singletons. Alpha = +inf when
    all taxa are singletons.

    When the sample is empty (i.e., all counts are zero), alpha = 0.

    References
    ----------
    .. [1] Fisher, R.A., Corbet, A.S. and Williams, C.B., 1943. The relation
       between the number of taxa and the number of individuals in a random
       sample of an animal population. The Journal of Animal Ecology, pp.42-58.

    """
    counts = _validate_counts_vector(counts)

    # alpha = 0 when sample has no individual
    if (N := counts.sum()) == 0:
        return 0.0

    # alpha = +inf when all taxa are singletons
    if N == (S := sobs(counts)):
        return np.inf

    # objective function to minimize:
    # S = alpha * ln (1 + N / alpha), where alpha > 0
    def f(x):
        return (x * np.log(1 + (N / x)) - S) ** 2 if x > 0 else np.inf

    # minimize the function using the default method (Brent's algorithm)
    res = minimize_scalar(f)

    # there is a chance optimization could fail
    if res.success is False:
        raise RuntimeError("Optimizer failed to solve for Fisher's alpha.")

    return res.x


def goods_coverage(counts):
    r"""Calculate Good's coverage estimator.

    Good's coverage estimator :math:`C` is an estimation of the proportion of
    the population represented in the sample. It is defined as

    .. math::

       C = 1 - \frac{F_1}{N}

    where :math:`F_1` is the number of taxa observed only once (i.e.,
    singletons) and :math:`N` is the total number of individuals.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        Good's coverage estimator.

    Notes
    -----
    Good's coverage estimator was originally described in [1]_.

    References
    ----------
    .. [1] Good, I. J. (1953). The population frequencies of species and the
       estimation of population parameters. Biometrika, 40(3-4), 237-264.

    """
    counts = _validate_counts_vector(counts)
    if (N := counts.sum()) == 0:
        return 0.0
    f1 = singles(counts)
    return 1 - (f1 / N)


def heip_e(counts):
    r"""Calculate Heip's evenness measure.

    Heip's evenness is defined as:

    .. math::

       \frac{(e^H-1)}{(S-1)}

    where :math:`H` is the Shannon-Wiener entropy of counts (using logarithm
    base :math:`e`) and :math:`S` is the number of taxa in the sample.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        Heip's evenness measure.

    See Also
    --------
    shannon
    pielou_e

    Notes
    -----
    Heip's evenness measure was originally described in [1]_.

    References
    ----------
    .. [1] Heip, C. 1974. A new index measuring evenness. J. Mar. Biol. Ass.
       UK., 54, 555-557.

    """
    counts = _validate_counts_vector(counts)
    return (np.exp(shannon(counts, base=np.e)) - 1) / (sobs(counts) - 1)


def kempton_taylor_q(counts, lower_quantile=0.25, upper_quantile=0.75):
    """Calculate Kempton-Taylor Q index of alpha diversity.

    Estimates the slope of the cumulative abundance curve in the interquantile
    range. By default, uses lower and upper quartiles, rounding inwards.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    lower_quantile : float, optional
        Lower bound of the interquantile range. Defaults to lower quartile.
    upper_quantile : float, optional
        Upper bound of the interquantile range. Defaults to upper quartile.

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
    Magurran recommends counting half of the taxa with the same # counts as the
    point where the UQ falls and the point where the LQ falls, but the
    justification for this is unclear (e.g. if there were a very large # taxa
    that just overlapped one of the quantiles, the results would be
    considerably off). Leaving the calculation as-is for now, but consider
    changing.

    References
    ----------
    .. [1] Kempton, R. A. and Taylor, L. R. (1976) Models and statistics for
       species diversity. Nature, 262, 818-820.
    .. [2] http://www.pisces-conservation.com/sdrhelp/index.html

    """
    counts = _validate_counts_vector(counts)
    n = len(counts)
    lower = int(np.ceil(n * lower_quantile))
    upper = int(n * upper_quantile)
    sorted_counts = np.sort(counts)
    return (upper - lower) / np.log(sorted_counts[upper] / sorted_counts[lower])


def margalef(counts):
    r"""Calculate Margalef's richness index.

    Margalef's richness index :math:`D` is defined as:

    .. math::

       D = \frac{(S - 1)}{\ln N}

    where :math:`S` is the number of taxa and :math:`N` is the total number
    of individuals in the sample.

    Assumes log accumulation.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        Margalef's richness index.

    Notes
    -----
    Margalef's richness index was originally described in [1]_.

    References
    ----------
    .. [1] Margalef, R. (1958) Information Theory in Ecology. General Systems,
       3, 36-71.

    """
    counts = _validate_counts_vector(counts)
    if (N := counts.sum()) == 0:
        return 0.0
    return (sobs(counts) - 1) / np.log(N)


def mcintosh_d(counts):
    r"""Calculate McIntosh dominance index.

    McIntosh dominance index :math:`D` is defined as:

    .. math::

       D = \frac{N - U}{N - \sqrt{N}}

    where :math:`N` is the total number of individuals in the sample and
    :math:`U` is defined as:

    .. math::

       U = \sqrt{\sum{{n_i}^2}}

    where :math:`n_i` is the number of individuals in the :math:`i^{\text{th}}`
    taxon.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        McIntosh dominance index.

    See Also
    --------
    mcintosh_e

    Notes
    -----
    McIntosh dominance index was originally described in [1]_.

    References
    ----------
    .. [1] McIntosh, R. P. 1967 An index of diversity and the relation of
       certain concepts to diversity. Ecology 48, 1115-1126.

    """
    counts = _validate_counts_vector(counts)
    if (N := counts.sum()) == 0:
        return 0.0
    u = np.sqrt((counts**2).sum())
    return (N - u) / (N - np.sqrt(N))


def mcintosh_e(counts):
    r"""Calculate McIntosh's evenness measure E.

    McIntosh evenness measure E is defined as:

    .. math::

       E = \frac{\sqrt{\sum{n_i^2}}}{\sqrt{((N-S+1)^2 + S -1}}

    where :math:`n_i` is the number of individuals in the :math:`i^{\text{th}}`
    taxon, :math:`N` is the total number of individuals, and :math:`S` is the
    number of taxa in the sample.

    Parameters
    ----------
    counts : 1-D array_like, int
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
    The implementation here is based on the description given in [1]_, **NOT**
    the one in the SDR-IV online manual, which is wrong.

    References
    ----------
    .. [1] Heip & Engels (1974) Comparing Species Diversity and Evenness
       Indices. p 560.

    """
    counts = _validate_counts_vector(counts)
    numerator = np.sqrt((counts * counts).sum())
    N = counts.sum()
    S = sobs(counts)
    denominator = np.sqrt((N - S + 1) ** 2 + S - 1)
    return numerator / denominator


def menhinick(counts):
    r"""Calculate Menhinick's richness index.

    Menhinick's richness index is defined as:

    .. math::

       D_{Mn} = \frac{S}{\sqrt{N}}

    where :math:`S` is the number of taxa and :math:`N` is the total number
    of individuals in the sample.

    Assumes square-root accumulation.

    Parameters
    ----------
    counts : 1-D array_like, int
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
    counts = _validate_counts_vector(counts)
    return sobs(counts) / np.sqrt(counts.sum())


def michaelis_menten_fit(counts, num_repeats=1, params_guess=None):
    r"""Calculate Michaelis-Menten fit to rarefaction curve of observed taxa.

    The Michaelis-Menten equation is defined as:

    .. math::

       S=\frac{nS_{max}}{n+B}

    where :math:`n` is the number of individuals and :math:`S` is the number of
    taxa. This function estimates the :math:`S_{max}` parameter.

    The fit is made to datapoints for :math:`n=1,2,...,N`, where :math:`N` is
    the total number of individuals (sum of abundances for all taxa).
    :math:`S` is the number of taxa represented in a random sample of
    :math:`n` individuals.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    num_repeats : int, optional
        The number of times to perform rarefaction (subsampling without
        replacement) at each value of :math:`n`.
    params_guess : tuple, optional
        Initial guess of :math:`S_{max}` and :math:`B`. If ``None``, default
        guess for :math:`S_{max}` is :math:`S` (as :math:`S_{max}` should
        be >= :math:`S`) and default guess for :math:`B` is ``round(N / 2)``.

    Returns
    -------
    S_max : double
        Estimate of the :math:`S_{max}` parameter in the Michaelis-Menten
        equation.

    See Also
    --------
    skbio.stats.subsample_counts

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
    counts = _validate_counts_vector(counts)

    n_indiv = counts.sum()
    if params_guess is None:
        S_max_guess = sobs(counts)
        B_guess = int(round(n_indiv / 2))
        params_guess = (S_max_guess, B_guess)

    # observed # of taxa vs # of individuals sampled, S vs n
    xvals = np.arange(1, n_indiv + 1)
    ymtx = np.empty((num_repeats, len(xvals)), dtype=int)
    for i in range(num_repeats):
        ymtx[i] = np.asarray(
            [sobs(subsample_counts(counts, n)) for n in xvals], dtype=int
        )
    yvals = ymtx.mean(0)

    # Vectors of actual vals y and number of individuals n.
    def errfn(p, n, y):
        return (((p[0] * n / (p[1] + n)) - y) ** 2).sum()

    # Return S_max.
    return fmin_powell(errfn, params_guess, ftol=1e-5, args=(xvals, yvals), disp=False)[
        0
    ]


def sobs(counts):
    """Calculate the observed species richness of a sample.

    Observed species richness, usually denoted as :math:`S_{obs}` or simply
    :math:`S`, is the number of distinct species (i.e., taxa), or any discrete
    groups of biological entities found in a sample.

    It should be noted that observed species richness is smaller than or equal
    to the true species richness of a population from which the sample is
    collected.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    int
        Observed species richness.

    See Also
    --------
    observed_features

    """
    counts = _validate_counts_vector(counts)
    return (counts != 0).sum()


def observed_features(counts):
    """Calculate the number of distinct features.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    int
        Distinct feature count.

    See Also
    --------
    sobs

    Notes
    -----
    `observed_features` is an alias for `sobs`.

    """
    return sobs(counts)


def observed_otus(counts):
    """Calculate the number of distinct OTUs.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    int
        Distinct OTU count.

    Warnings
    --------
    ``observed_otus`` is deprecated as of ``0.6.0`` due to its usage of the
    historical term "OTU".

    See Also
    --------
    sobs

    Notes
    -----
    `observed_otus` is an alias for `sobs`.

    """
    # @deprecated
    _warn_deprecated(observed_otus, "0.6.0")

    return sobs(counts)


def osd(counts):
    """Calculate observed taxa, singletons, and doubletons.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    osd : tuple
        Numbers of observed taxa, singletons, and doubletons.

    See Also
    --------
    sobs
    singles
    doubles

    Notes
    -----
    This is a convenience function used by many of the other measures that rely
    on these three measures.

    """
    counts = _validate_counts_vector(counts)
    return sobs(counts), singles(counts), doubles(counts)


def pielou_e(counts):
    r"""Calculate Pielou's evenness index.

    Pielou's evenness index (:math:`J'`), a.k.a., Shannon's equitability index
    (:math:`E_H`), is defined as

    .. math::

       J' = \frac{(H)}{\ln(S)}

    where :math:`H` is the Shannon index of the sample and :math:`S` is the
    number of taxa in the sample.

    That is, :math:`J'` is the ratio of the actual Shannon index of the sample
    versus the maximum-possible Shannon index when all taxa have the same
    number of individuals. :math:`J'` ranges between 0 and 1.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        Pielou's evenness index.

    See Also
    --------
    shannon
    heip_e

    Notes
    -----
    Pielou's evenness index was originally described in [1]_.

    References
    ----------
    .. [1] Pielou, E. C., 1966. The measurement of diversity in different types
       of biological collections. Journal of Theoretical Biology, 13, 131-44.

    """
    counts = _validate_counts_vector(counts)
    return 0.0 if (H := shannon(counts, base=np.e)) == 0.0 else H / np.log(sobs(counts))


def robbins(counts):
    r"""Calculate Robbins' estimator for probability of unobserved outcomes.

    Robbins' estimator is defined as:

    .. math::

       \frac{F_1}{n+1}

    where :math:`F_1` is the number of singleton taxa.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        Robbins' estimate.

    Notes
    -----
    Robbins' estimator is defined in [1]_. The estimate computed here is for
    :math:`n-1` counts, i.e. the x-axis is off by 1.

    References
    ----------
    .. [1] Robbins, H. E (1968). Ann. of Stats. Vol 36, pp. 256-257.

    """
    counts = _validate_counts_vector(counts)
    return singles(counts) / counts.sum()


def shannon(counts, base=2):
    r"""Calculate Shannon's diversity index, default in bits.

    Shannon's diversity index, :math:`H'`, a.k.a., Shannon index, or Shannon-
    Wiener index, is defined as

    .. math::

       H' = -\sum_{i=1}^s\left(p_i\log_2 p_i\right)

    where :math:`s` is the number of taxa and :math:`p_i` is the proportion
    of the sample represented by taxon :math:`i`.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    base : scalar, optional
        Logarithm base to use in the calculation.

    Returns
    -------
    double
        Shannon's diversity index.

    Notes
    -----
    Shannon's diversity index was originally proposed in [1]_ as a measure of
    entropy.

    The default logarithm base used here is 2 instead of :math:`e`.

    References
    ----------
    .. [1] Shannon, C. E. (1948). A mathematical theory of communication. The
       Bell system technical journal, 27(3), 379-423.

    """
    counts = _validate_counts_vector(counts)
    if (N := counts.sum()) == 0:
        return 0.0
    freqs = counts / N
    nonzero_freqs = freqs[freqs.nonzero()]
    if nonzero_freqs.size <= 1:
        return 0.0
    return -(nonzero_freqs * np.log(nonzero_freqs)).sum() / np.log(base)


def simpson(counts):
    r"""Calculate Simpson's diversity index.

    Simpson's diversity index, a.k.a., Gini-Simpson index, or Gini impurity,
    is defined as

    .. math::

       1 - \sum{p_i^2}

    where :math:`p_i` is the proportion of the sample represented by taxon
    :math:`i`.

    Therefore, Simpson's diversity index is also denoted as :math:`1 - D`, in
    which :math:`D` is the Simpson's dominance index.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        Simpson's diversity index.

    See Also
    --------
    dominance

    Notes
    -----
    Simpson's diversity index was originally described in [1]_.

    References
    ----------
    .. [1] Simpson, E. H. (1949). Measurement of diversity. Nature, 163(4148),
       688-688.

    """
    counts = _validate_counts_vector(counts)
    return 1 - dominance(counts)


def simpson_e(counts):
    r"""Calculate Simpson's evenness index.

    Simpson's evenness (a.k.a., equitability) index :math:`E_D` is defined as

    .. math::

       E_D = \frac{1}{D \times S}

    where :math:`D` is the Simpson's dominance index and :math:`S` is the
    number of taxa in the sample.

    That is, :math:`E_D` is the ratio of the minimum-possible Simpson's
    dominance index when all taxa have the same number of individuals:
    :math:`D_{min} = 1 / S`, versus the actual Simpson's dominance index of the
    sample.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        Simpson's evenness index.

    See Also
    --------
    dominance
    enspie
    simpson

    Notes
    -----
    The implementation here is based on the description given in [1]_ and [2]_.

    References
    ----------
    .. [1] Simpson, E. H. (1949). Measurement of diversity. nature, 163(4148), 688-688.
    .. [2] Pielou, E. C. (1966). The measurement of diversity in different types of
       biological collections. Journal of theoretical biology, 13, 131-144.

    """
    counts = _validate_counts_vector(counts)
    return enspie(counts) / sobs(counts)


def singles(counts):
    """Calculate number of single-occurrence taxa (singletons).

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    int
        Singleton count.

    """
    counts = _validate_counts_vector(counts)
    return (counts == 1).sum()


def strong(counts):
    r"""Calculate Strong's dominance index.

    Strong's dominance index (:math:`D_w`) is defined as

    .. math::

       D_w = max_i[(\frac{b_i}{N})-\frac{i}{S}]

    where :math:`b_i` is the sequential cumulative totaling of the
    :math:`i^{\text{th}}` taxon abundance values ranked from largest to
    smallest, :math:`N` is the total number of individuals in the sample, and
    :math:`S` is the number of taxa in the sample. The expression in
    brackets is computed for all taxa, and :math:`max_i` denotes the maximum
    value in brackets for any taxa.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    double
        Strong's dominance index.

    Notes
    -----
    Strong's dominance index is defined in [1]_.

    References
    ----------
    .. [1] Strong, W. L., 2002 Assessing species abundance unevenness within
       and between plant communities. Community Ecology, 3, 237-246.

    """
    counts = _validate_counts_vector(counts)
    if (N := counts.sum()) == 0:
        return 0.0
    S = sobs(counts)
    i = np.arange(1, len(counts) + 1)
    sorted_sum = np.sort(counts)[::-1].cumsum()
    return (sorted_sum / N - (i / S)).max()
