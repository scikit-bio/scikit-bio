# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from warnings import warn
import functools

import numpy as np
from scipy.special import gammaln
from scipy.optimize import fmin_powell, minimize_scalar

from skbio.stats import subsample_counts
from skbio.diversity._util import _validate_counts_vector
from skbio.util._warning import _warn_deprecated


def _validate_alpha(empty=None, cast_int=False):
    """Validate counts vector for an alpha diversity metric.

    Parameters
    ----------
    func : callable
        Function that calculates an alpha diversity metric.
    empty : any, optional
        Return this value if set instead of calling the function when an input
        community is empty (i.e., no taxon, or all taxa have zero counts).
    cast_int : bool, optional
        Cast values into integers, if not already. ``False`` by default.

    Returns
    -------
    callable
        Decorated function.

    Notes
    -----
    This function serves as a decorator for individual functions that calculate
    alpha diversity metrics. The first positional argument of a decorated
    function must be a 1-D vector of counts/abundances of taxa in a community.
    Additional arguments may follow.

    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(counts, *args, **kwargs):
            counts = _validate_counts_vector(counts, cast_int)

            # drop zero values, as these represent taxa that are absent from
            # the community
            if not (nonzero := counts != 0).all():
                counts = counts[nonzero]

            # return a value if community is empty (after dropping zeros)
            if empty is not None and counts.size == 0:
                return empty

            # call function to calculate alpha diversity metric
            return func(counts, *args, **kwargs)

        return wrapper

    return decorator


@_validate_alpha(empty=np.nan)
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
    float
        Berger-Parker dominance index.

    Notes
    -----
    Berger-Parker dominance index was originally described in [1]_.

    References
    ----------
    .. [1] Berger, W. H., & Parker, F. L. (1970). Diversity of planktonic
       foraminifera in deep-sea sediments. Science, 168(3937), 1345-1347.

    """
    return counts.max() / counts.sum()


@_validate_alpha(empty=np.nan)
def brillouin_d(counts):
    r"""Calculate Brillouin's diversity index.

    Brillouin's diversity index (:math:`H_B`) is defined as:

    .. math::

       H_B = \frac{\ln N!-\sum_{i=1}^S{\ln n_i!}}{N}

    where :math:`N` is the total number of individuals in the sample, :math:`S`
    is the number of taxa, and :math:`n_i` is the number of individuals in the
    :math:`i^{\text{th}}` taxon.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    float
        Brillouin's diversity index.

    Notes
    -----
    Brillouin's diversity index was originally described in [1]_.

    References
    ----------
    .. [1] Brillouin, L. (1956). Science and Information Theory. Academic
       Press. New York.

    """
    return (gammaln((N := counts.sum()) + 1) - gammaln(counts + 1).sum()) / N


@_validate_alpha(empty=np.nan)
def dominance(counts, finite=False):
    r"""Calculate Simpson's dominance index.

    Simpson's dominance index, a.k.a. Simpson's :math:`D`, measures the degree
    of concentration of taxon composition of a sample. It is defined as:

    .. math::

       D = \sum_{i=1}^S{p_i^2}

    where :math:`S` is the number of taxa and :math:`p_i` is the proportion
    of the sample represented by taxon :math:`i`.

    Simpson's :math:`D` ranges from 0 (infinite diversity; no dominance) and 1
    (complete dominance, no diversity).

    Simpson's :math:`D` can be interpreted as the probability that two randomly
    selected individuals belong to the same taxon.

    Simpson's :math:`D` may be corrected for finite samples to account for the
    effect of sampling without replacement. This more accurately represents the
    above probability when the sample is small. It is calculated as:

    .. math::

       D = \frac{\sum_{i=1}^s{n_i(n_i - 1))}}{N(N - 1)}

    where :math:`n_i` is the number of individuals in the :math:`i^{\text{th}}`
    taxon and :math:`N` is the total number of individuals in the sample.

    Simpson's :math:`D` is sometimes referred to as "Simpson's index". It
    should be noted that :math:`D` is not a measure of community diversity. It
    is also important to distinguish :math:`D` from Simpson's diversity index
    (:math:`1 - D`) and inverse Simpson index (:math:`1 / D`), both of which
    are measures of community diversity.

    Discrepancy exists among literature in using the term "Simpson index" and
    the denotion :math:`D`. It is therefore important to distinguish these
    metrics according to their mathematic definition.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    finite : bool, optional
        If ``True``, correct for finite sampling.

    Returns
    -------
    float
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
    if finite:
        D = (counts * (counts - 1)).sum() / ((N := counts.sum()) * (N - 1))
    else:
        D = ((counts / counts.sum()) ** 2).sum()
    return D


@_validate_alpha()
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
    return (counts == 2).sum()


def enspie(counts, finite=False):
    r"""Calculate ENS_pie alpha diversity measure.

    The effective number of species (ENS) derived from Hurlbert's probability
    of interspecific encounter (PIE) ([1]_, [2]_) is defined as:

    .. math::

       ENS_{pie} = \frac{1}{\sum_{i=1}^S{p_i^2}}

    where :math:`S` is the number of taxa and :math:`p_i` is the proportion of
    the sample represented by taxon :math:`i`.

    Therefore, :math:`ENS_{pie}` is equivalent to the inverse Simpson index
    (``1 / D``).

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    finite : bool, optional
        If ``True``, correct for finite sampling.

    Returns
    -------
    float
        ENS_pie alpha diversity measure.

    See Also
    --------
    inv_simpson
    dominance

    Notes
    -----
    ``enspie`` is an alias for ``inv_simpson``.

    References
    ----------
    .. [1] Chase, J. M., & Knight, T. M. (2013). Scale-dependent effect sizes
       of ecological drivers on biodiversity: why standardised sampling is not
       enough. Ecology letters, 16, 17-26.

    .. [2] Hurlbert, S. H. (1971). The nonconcept of species diversity: a
       critique and alternative parameters. Ecology, 52(4), 577-586.

    """
    return inv_simpson(counts, finite=finite)


@_validate_alpha(empty=np.nan)
def esty_ci(counts):
    r"""Calculate Esty's confidence interval of Good's coverage estimator.

    Esty's confidence interval is defined as:

    .. math::

       F_1/N \pm z\sqrt{W}

    where :math:`F_1` is the number of singleton taxa, :math:`N` is the
    total number of individuals, and :math:`z` is a constant that depends on
    the targeted confidence and based on the normal distribution.

    :math:`W` is defined as:

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

    See Also
    --------
    goods_coverage

    Notes
    -----
    Esty's confidence interval was originally described in [1]_.

    :math:`z` is hardcoded for a 95% confidence interval.

    References
    ----------
    .. [1] Esty, W. W. (1983). "A normal limit law for a nonparametric
       estimator of the coverage of a random sample". Ann Statist 11: 905-912.

    """
    N = counts.sum()
    f1 = (counts == 1).sum()
    f2 = (counts == 2).sum()
    z = 1.959963985
    W = (f1 * (N - f1) + 2 * N * f2) / (N**3)
    return f1 / N - z * np.sqrt(W), f1 / N + z * np.sqrt(W)


@_validate_alpha(empty=np.nan)
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
    float
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
    # alpha = +inf when all taxa are singletons
    if (N := counts.sum()) == (S := counts.size):
        return np.inf

    # objective function to minimize:
    # S = alpha * ln (1 + N / alpha), where alpha > 0
    def f(x):
        return (x * np.log(1 + (N / x)) - S) ** 2 if x > 0 else np.inf

    # minimize the function using the default method (Brent's algorithm)
    with np.errstate(invalid="ignore"):
        res = minimize_scalar(f)

    # there is a chance optimization could fail
    if res.success is False:
        raise RuntimeError("Optimizer failed to solve for Fisher's alpha.")

    return res.x


@_validate_alpha(empty=np.nan)
def goods_coverage(counts):
    r"""Calculate Good's coverage estimator.

    Good's coverage estimator :math:`C`, a.k.a. Turing estimator or Good-
    Turing (GT) estimator, is an estimation of the proportion of the
    population represented in the sample. It is defined as:

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
    float
        Good's coverage estimator.

    See Also
    --------
    esty_ci
    robbins

    Notes
    -----
    Good's coverage estimator was originally described in [1]_.

    References
    ----------
    .. [1] Good, I. J. (1953). The population frequencies of species and the
       estimation of population parameters. Biometrika, 40(3-4), 237-264.

    """
    return 1 - ((counts == 1).sum() / counts.sum())


@_validate_alpha()
def heip_e(counts):
    r"""Calculate Heip's evenness measure.

    Heip's evenness is defined as:

    .. math::

       \frac{(e^H-1)}{(S-1)}

    where :math:`H` is Shannon's diversity index and :math:`S` is the number
    of taxa in the sample.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    float
        Heip's evenness measure.

    See Also
    --------
    shannon
    pielou_e

    Notes
    -----
    Heip's evenness measure was originally described in [1]_.

    When there is only one taxon, the return value is 1.0.

    References
    ----------
    .. [1] Heip, C. 1974. A new index measuring evenness. J. Mar. Biol. Ass.
       UK., 54, 555-557.

    """
    if (S := counts.size) == 0:
        return np.nan
    elif S == 1:
        return 1.0
    return (shannon(counts, exp=True) - 1) / (S - 1)


@_validate_alpha(empty=np.nan)
def hill(counts, order=2):
    r"""Calculate Hill number.

    Hill number (:math:`^qD`) is a generalized measure of the effective number
    of species. It is defined as:

    .. math::

       ^qD = (\sum_{i=1}^S p_i^q)^{\frac{1}{1-q}}

    where :math:`S` is the number of taxa and :math:`p_i` is the proportion
    of the sample represented by taxon :math:`i`.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    order : int or float, optional
        Order (:math:`q`). Ranges between 0 and infinity. Default is 2.

    Returns
    -------
    float
        Hill number.

    See Also
    --------
    inv_simpson
    renyi
    shannon
    sobs

    Notes
    -----
    Hill number was originally defined in [1]_. It is a measurement of "true
    diversity", or the effective number of species (ENS) ([2]_), which is
    defined as the number of equally abundant taxa that would make the same
    diversity measurement given the observed total abundance of the community.

    Hill number is a generalization of multiple diversity metrics. Depending on
    the order :math:`q`, it is equivalent to:

    - :math:`q=0`: Observed species richness (:math:`S_{obs}`).
    - :math:`q \to 1`: The exponential of Shannon index (:math:`\exp{H'}`),
      i.e., perplexity.
    - :math:`q=2`: Inverse Simpson index (:math:`1 / D`).
    - :math:`q \to \infty`: :math:`1/\max{p}`, i.e., the inverse of
      Berger-Parker dominance index.

    The order :math:`q` determines the influence of taxon abundance on the
    metric. A larger (or smaller) :math:`q` puts more weight on the abundant
    (or rare) taxa.

    Hill number is equivalent to the exponential of Renyi entropy.

    References
    ----------
    .. [1] Hill, M. O. (1973). Diversity and evenness: a unifying notation and
       its consequences. Ecology, 54(2), 427-432.

    .. [2] Jost, L. (2006). Entropy and diversity. Oikos, 113(2), 363-375.

    """
    probs = counts / counts.sum()
    if order == 1:
        return _perplexity(probs)
    elif np.isposinf(order):
        return 1 / probs.max()
    else:
        return (probs**order).sum() ** (1 / (1 - order))


@_validate_alpha(empty=np.nan)
def kempton_taylor_q(counts, lower_quantile=0.25, upper_quantile=0.75):
    r"""Calculate Kempton-Taylor Q index of alpha diversity.

    Kempton-Taylor Q index measures diversity based on the middle-ranking taxa
    in the abundance distribution. Specifically, it estimates the slope of the
    cumulative abundance curve in the interquantile range. It is defined as:

    .. math::

       Q = \frac{S_{lower..upper}}{\ln n_{lower} - \ln n_{upper}}

    where "lower" and "upper" are the taxa at the lower and upper quantiles of
    the abundance distribution, :math:`S` is the number of taxa, and :math:`n`
    is the number of individuals.

    By default, the lower and upper quartiles are used. Therefore:

    .. math::

       Q = \frac{S}{2(\ln n_{0.25} - \ln n_{0.75})}

    The quantiles are rounded inwards in this implementation.

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
    float
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
    S = counts.size
    lower = int(np.ceil(S * lower_quantile))
    upper = int(S * upper_quantile)
    sorted_counts = np.sort(counts)
    return (upper - lower) / np.log(sorted_counts[upper] / sorted_counts[lower])


def inv_simpson(counts, finite=False):
    r"""Calculate inverse Simpson index.

    The inverse Simpson index (:math:`1 / D`), a.k.a., Simpson's reciprocal
    index, is defined as:

    .. math::

       1 / D = \frac{1}{\sum_{i=1}^S{p_i^2}}

    where :math:`S` is the number of taxa and :math:`p_i` is the proportion
    of the sample represented by taxon :math:`i`.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    finite : bool, optional
        If ``True``, correct for finite sampling when calculating :math:`D`.

    Returns
    -------
    float
        Inverse Simpson index.

    See Also
    --------
    dominance

    Notes
    -----
    :math:`1 / D` is a measurement of the effective number of species (ENS).
    It is equivalent to Hill number with order 2 (:math:`^2D`).

    Inverse Simpson index was originally described in [1]_.

    References
    ----------
    .. [1] Simpson, E. H. (1949). Measurement of diversity. Nature, 163(4148),
       688-688.

    """
    return 1 / dominance(counts, finite=finite)


@_validate_alpha(empty=np.nan)
def margalef(counts):
    r"""Calculate Margalef's richness index.

    Margalef's richness index :math:`D` is defined as:

    .. math::

       D = \frac{(S - 1)}{\ln N}

    where :math:`S` is the number of taxa and :math:`N` is the total number
    of individuals in the sample.

    Margalef's richness index assumes log accumulation.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    float
        Margalef's richness index.

    See Also
    --------
    menhinick

    Notes
    -----
    Margalef's richness index was originally described in [1]_.

    References
    ----------
    .. [1] Margalef, R. (1958) Information Theory in Ecology. General Systems,
       3, 36-71.

    """
    if (N := counts.sum()) == 1:
        return np.nan
    return (counts.size - 1) / np.log(N)


@_validate_alpha(empty=np.nan)
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
    float
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
    if (N := counts.sum()) == 1:
        return np.nan
    u = np.sqrt((counts**2).sum())
    return (N - u) / (N - np.sqrt(N))


@_validate_alpha(empty=np.nan)
def mcintosh_e(counts):
    r"""Calculate McIntosh's evenness measure.

    McIntosh's evenness measure :math:`E` is defined as:

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
    float
        McIntosh evenness measure.

    See Also
    --------
    mcintosh_d

    Notes
    -----
    McIntosh's evenness measure was originally described in [1]_.

    References
    ----------
    .. [1] Heip & Engels (1974) Comparing Species Diversity and Evenness
       Indices. p 560.

    """
    S = counts.size
    N = counts.sum()
    numerator = np.sqrt((counts * counts).sum())
    denominator = np.sqrt((N - S + 1) ** 2 + S - 1)
    return numerator / denominator


@_validate_alpha(empty=np.nan)
def menhinick(counts):
    r"""Calculate Menhinick's richness index.

    Menhinick's richness index is defined as:

    .. math::

       D_{Mn} = \frac{S}{\sqrt{N}}

    where :math:`S` is the number of taxa and :math:`N` is the total number
    of individuals in the sample.

    Menhinick's richness index assumes square-root accumulation.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    float
        Menhinick's richness index.

    See Also
    --------
    margalef

    Notes
    -----
    Based on the description in [1]_.

    References
    ----------
    .. [1] Magurran, A E 2004. Measuring biological diversity. Blackwell. pp.
       76-77.

    """
    return counts.size / np.sqrt(counts.sum())


@_validate_alpha(empty=np.nan)
def michaelis_menten_fit(counts, num_repeats=1, params_guess=None):
    r"""Calculate Michaelis-Menten fit to rarefaction curve of observed taxa.

    The Michaelis-Menten equation estimates the asymptote of the rarefaction
    curve. It is an estimator of the true richness of a community given the
    observation. It is defined as:

    .. math::

       S = \frac{nS_{max}}{n+B}

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
    float
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
    ``observed_features`` is an alias for ``sobs``.

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
    ``observed_otus`` is an alias for ``sobs``.

    """
    # @deprecated
    _warn_deprecated(observed_otus, "0.6.0")

    return sobs(counts)


@_validate_alpha()
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
    return counts.size, (counts == 1).sum(), (counts == 2).sum()


@_validate_alpha()
def pielou_e(counts, base=None):
    r"""Calculate Pielou's evenness index.

    Pielou's evenness index (:math:`J'`), a.k.a., Shannon's equitability index
    (:math:`E_H`), is defined as:

    .. math::

       J' = \frac{(H)}{\log(S)}

    where :math:`H` is the Shannon index of the sample and :math:`S` is the
    number of taxa in the sample.

    That is, :math:`J'` is the ratio of the actual Shannon index of the sample
    versus the maximum-possible Shannon index when all taxa have the same
    number of individuals. :math:`J'` ranges between 0 and 1.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    base : int or float, optional
        Logarithm base to use in the calculation. Default is ``e``.

    Returns
    -------
    float
        Pielou's evenness index.

    See Also
    --------
    shannon
    heip_e

    Notes
    -----
    Pielou's evenness index was originally described in [1]_.

    When there is only one taxon, the return value is 1.0.

    References
    ----------
    .. [1] Pielou, E. C., 1966. The measurement of diversity in different types
       of biological collections. Journal of Theoretical Biology, 13, 131-44.

    """
    if (S := counts.size) == 0:
        return np.nan
    elif S == 1:
        return 1.0
    H = shannon(counts, base=base)
    H_max = np.log(S)
    if base is not None:
        H_max /= np.log(base)
    return H / H_max


@_validate_alpha()
def renyi(counts, order=2, base=None):
    r"""Calculate Renyi entropy.

    Renyi entropy (:math:`^qH`) is a generalization of Shannon index, with an
    exponent (order) :math:`q` instead of 1. It is defined as:

    .. math::

       ^qH = \frac{1}{1-q}\log_b{(\sum_{i=1}^S p_i^q)}

    where :math:`S` is the number of taxa and :math:`p_i` is the proportion
    of the sample represented by taxon :math:`i`.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    order : int or float, optional
        Order (:math:`q`). Ranges between 0 and infinity. Default is 2.
    base : int or float, optional
        Logarithm base to use in the calculation. Default is ``e``.

    Returns
    -------
    float
        Renyi entropy.

    See Also
    --------
    hill
    inv_simpson
    shannon
    tsallis

    Notes
    -----
    Renyi entropy was originally defined in [1]_. It is a generalization of
    multiple entropy notions, as determined by the order (:math:`q`). Special
    cases of Renyi entropy include:

    - :math:`q=0`: Max-entropy (:math:`\log{S}`).
    - :math:`q \to 1`: Shannon entropy (index).
    - :math:`q=2`: Collision entropy, a.k.a, Renyi's quadratic entropy, or
      "Renyi entropy". Equivalent to the logarithm of inverse Simpson index.
    - :math:`q \to \infty`: Min-entropy (:math:`-\log{\max{p}}`).

    Renyi entropy is equivalent to the logarithm of Hill number.

    References
    ----------
    .. [1] Rényi, A. (1961, January). On measures of entropy and information.
       In Proceedings of the fourth Berkeley symposium on mathematical
       statistics and probability, volume 1: contributions to the theory of
       statistics (Vol. 4, pp. 547-562). University of California Press.

    """
    if (S := counts.size) == 0:
        return np.nan
    elif S == 1:
        return 0.0

    probs = counts / counts.sum()

    # max-entropy
    if order == 0:
        qH = np.log(S)
    # Shannon entropy
    elif order == 1:
        qH = _entropy(probs)
    # min-entropy
    elif np.isposinf(order):
        qH = -np.log(probs.max())
    else:
        qH = np.log((probs**order).sum()) / (1 - order)

    if base is not None:
        qH /= np.log(base)
    return qH


@_validate_alpha(empty=np.nan)
def robbins(counts):
    r"""Calculate Robbins' estimator for probability of unobserved outcomes.

    Robbins' estimator is defined as:

    .. math::

       \frac{F_1}{N}

    where :math:`F_1` is the number of singleton taxa and and :math:`N` is the
    total number of individuals in the sample.

    The result can be interpreted as the probability of discovering a new taxon
    at the :math:`N`-th individual given the current :math:`N - 1` individuals.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.

    Returns
    -------
    float
        Robbins' estimator.

    See Also
    --------
    goods_coverage

    Notes
    -----
    Robbins' estimator is defined in [1]_.

    References
    ----------
    .. [1] Robbins, H. E. (1968). Estimating the total probability of the
       unobserved outcomes of an experiment. Ann. Math. Statist., 39(6),
       256-257.

    """
    return (counts == 1).sum() / counts.sum()


def _entropy(probs):
    """Calculate entropy."""
    return (-probs * np.log(probs)).sum()


def _perplexity(probs):
    """Calculate perplexity."""
    return (probs**-probs).prod()


@_validate_alpha(empty=np.nan)
def shannon(counts, base=None, exp=False):
    r"""Calculate Shannon's diversity index.

    Shannon's diversity index, :math:`H'`, a.k.a., Shannon index, or Shannon-
    Wiener index, is equivalent to entropy in information theory. It is defined
    as:

    .. math::

       H' = -\sum_{i=1}^S\left(p_i\log_b(p_i)\right)

    where :math:`S` is the number of taxa and :math:`p_i` is the proportion
    of the sample represented by taxon :math:`i`.

    The logarithm base :math:`b` defaults to ``e``, but may be 2, 10 or other
    custom values.

    The exponential of Shannon index, :math:`exp(H')`, measures the effective
    number of species (a.k.a., true diversity). It is equivalent to perplexity
    in information theory, or Hill number with order 1 (:math:`^1D`). The value
    is independent from the base:

    .. math::

       exp(H') = b ^ {-\sum_{i=1}^S\left(p_i\log_b(p_i)\right)} = \prod_{i=1}
       ^{S}p_i^{-p_i}

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    base : int or float, optional
        Logarithm base to use in the calculation. Default is ``e``.

        .. versionchanged:: 0.6.1

            The default logarithm base was changed from 2 to :math:`e` for
            consistency with the majority of literature.

    exp : bool, optional
        If ``True``, return the exponential of Shannon index.

    Returns
    -------
    float
        Shannon's diversity index.

    Notes
    -----
    Shannon index (i.e., entropy) was originally proposed in [1]_. The
    exponential of Shannon index (i.e., perplexity) was discussed in [2]_ in
    the context of community diversity.

    References
    ----------
    .. [1] Shannon, C. E. (1948). A mathematical theory of communication. The
       Bell system technical journal, 27(3), 379-423.

    .. [2] Jost, L. (2006). Entropy and diversity. Oikos, 113(2), 363-375.

    """
    probs = counts / counts.sum()

    # perplexity
    if exp is True:
        return _perplexity(probs)

    # entropy
    else:
        H = _entropy(probs)
        if base is not None:
            H /= np.log(base)
        return H


def simpson(counts, finite=False):
    r"""Calculate Simpson's diversity index.

    Simpson's diversity index, a.k.a., Gini-Simpson index, or Gini impurity,
    is defined as:

    .. math::

       1 - \sum_{i=1}^S{p_i^2}

    where :math:`S` is the number of taxa and :math:`p_i` is the proportion
    of the sample represented by taxon :math:`i`.

    Therefore, Simpson's diversity index is also denoted as :math:`1 - D`, in
    which :math:`D` is the Simpson's dominance index.

    Simpson's diversity index can be interpreted as the probability that two
    randomly selected individuals belong to different taxa. It is also known
    as Hurlbert's probability of interspecific encounter (PIE).

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    finite : bool, optional
        If ``True``, correct for finite sampling when calculating :math:`D`.

    Returns
    -------
    float
        Simpson's diversity index.

    See Also
    --------
    dominance

    Notes
    -----
    Simpson's diversity index was originally described in [1]_.

    Hurlbert's probability of interspecific encounter was described in [2]_.

    References
    ----------
    .. [1] Simpson, E. H. (1949). Measurement of diversity. Nature, 163(4148),
       688-688.

    .. [2] Hurlbert, S. H. (1971). The nonconcept of species diversity: a
       critique and alternative parameters. Ecology, 52(4), 577-586.

    """
    return 1 - dominance(counts, finite=finite)


def simpson_d(counts, finite=False):
    """Calculate Simpson's dominance index, a.k.a. Simpson's D.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    finite : bool, optional
        If ``True``, correct for finite sampling.

    Returns
    -------
    int
        Simpson's dominance index.

    See Also
    --------
    dominance
    simpson
    simpson_e

    Notes
    -----
    ``simpson_d`` is an alias for ``dominance``.

    """
    return dominance(counts, finite=finite)


@_validate_alpha(empty=np.nan)
def simpson_e(counts):
    r"""Calculate Simpson's evenness index.

    Simpson's evenness (a.k.a., equitability) index :math:`E_D` is defined as:

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
    float
        Simpson's evenness index.

    See Also
    --------
    dominance
    simpson

    Notes
    -----
    The implementation here is based on the description given in [1]_ and [2]_.

    References
    ----------
    .. [1] Simpson, E. H. (1949). Measurement of diversity. nature, 163(4148),
       688-688.

    .. [2] Pielou, E. C. (1966). The measurement of diversity in different
       types of biological collections. Journal of theoretical biology, 13,
       131-144.

    """
    # Note: the finite version of simpson_e might be: 1 / (D(S + 1)), because
    # S + 1 is the maximum possible finite D given S. Otherwise, the result can
    # be greater than 1 for small samples. However, I didn't find literature
    # stating this. Therefore, the `finite` parameter is not used here.
    return 1 / (counts.size * dominance(counts))


@_validate_alpha()
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
    return (counts == 1).sum()


@_validate_alpha()
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
    return counts.size


@_validate_alpha(empty=np.nan)
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
    float
        Strong's dominance index.

    Notes
    -----
    Strong's dominance index is defined in [1]_.

    References
    ----------
    .. [1] Strong, W. L., 2002 Assessing species abundance unevenness within
       and between plant communities. Community Ecology, 3, 237-246.

    """
    S = counts.size
    sorted_sum = np.sort(counts)[::-1].cumsum()
    i = np.arange(1, S + 1)
    return (sorted_sum / counts.sum() - (i / S)).max()


@_validate_alpha()
def tsallis(counts, order=2):
    r"""Calculate Tsallis entropy.

    Tsallis entropy (:math:`^qH`), a.k.a. HCDT entropy, is a generalization of
    Boltzmann-Gibbs entropy with an exponent (order) :math:`q`. It is defined
    as:

    .. math::

       ^qH = \frac{1}{q - 1}(1 - \sum_{i=1}^S p_i^q)

    where :math:`S` is the number of taxa and :math:`p_i` is the proportion
    of the sample represented by taxon :math:`i`.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    order : int or float, optional
        Order (:math:`q`). Ranges between 0 and infinity. Default is 2.

    Returns
    -------
    float
        Tsallis entropy.

    See Also
    --------
    renyi
    shannon
    simpson
    sobs

    Notes
    -----
    Tsallis entropy was originally defined in [1]_. Special cases of Tsallis
    entropy given order :math:`q` include:

    - :math:`q=0`: Observed species richness (:math:`S_{obs}`) minus 1.
    - :math:`q \to 1`: Shannon index :math:`H'`.
    - :math:`q=2`: Simpson diversity index (:math:`1 - D`).
    - :math:`q \to \infty`: 0.

    References
    ----------
    .. [1] Tsallis, C. (1988). Possible generalization of Boltzmann-Gibbs
       statistics. Journal of statistical physics, 52, 479-487.

    """
    if (S := counts.size) == 0:
        return np.nan
    elif S == 1:
        return 0.0
    probs = counts / counts.sum()
    if order == 1:
        return _entropy(probs)
    elif np.isposinf(order):
        return 0.0
    else:
        return (1 - (probs**order).sum()) / (order - 1)
