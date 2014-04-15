#!/usr/bin/env python
r"""
Alpha diversity measures (:mod:`skbio.maths.diversity.alpha`)
=============================================================

.. currentmodule:: skbio.maths.diversity.alpha

This module provides implementations of various alpha diversity measures.

Functions
---------

.. autosummary::
   :toctree: generated/

   observed_species

"""
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


def berger_parker_d(counts):
    """Fraction of the sample that belongs to the most abundant species.

    References
    ----------
    .. [1] Berger & Parker 1970, by way of SDR-IV online help.

    """
    return counts.max() / counts.sum()


def brillouin_d(counts):
    """Calculate Brilloun index of alpha diversity.

    References
    ----------
    .. [1] Pielou 1975, by way of SDR-IV.

    """
    nz = counts[counts.nonzero()]
    n = nz.sum()
    return (gammaln(n + 1) - gammaln(nz + 1).sum()) / n


def dominance(counts):
    """Calculate probability that two species sampled are the same.

    Dominance = 1 - Simpson's index, sum of squares of probabilities.

    """
    freqs = counts / counts.sum()
    return (freqs * freqs).sum()


def doubles(counts):
    """Return count of double occurrences."""
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
    return 1 / dominance(counts)

# For backwards-compatibility with QIIME.
simpson_reciprocal = enspie


def equitability(counts, base=2):
    """Calculate Shannon index corrected for number of species, pure evenness.

    """
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
    n1 = singles(counts)
    n2 = doubles(counts)
    n = counts.sum()
    z = 1.959963985
    W = (n1 * (n - n1) + 2 * n * n2) / (n ** 3)

    return n1 / n + z * np.sqrt(W), n1 / n - z * np.sqrt(W)


def heip_e(counts):
    """Calculate Heip's evenness measure.

    References
    ----------
    .. [1] Heip & Engels 1974.

    """
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
    return (observed_species(counts) - 1) / np.log(counts.sum())


def mcintosh_d(counts):
    """Calculate McIntosh index of alpha diversity.

    References
    ----------
    .. [1] McIntosh 1967, by way of SDR-IV.

    """
    u = np.sqrt((counts * counts).sum())
    n = counts.sum()
    return (n - u) / (n - np.sqrt(n))


def mcintosh_e(counts):
    """Calculate McIntosh's evenness measure.

    References
    ----------
    .. [1] Heip & Engels 1974 p 560 (wrong in SDR-IV).

    """
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
    return observed_species(counts) / np.sqrt(counts.sum())


def observed_species(counts):
    """Calculate number of distinct species."""
    return (counts != 0).sum()


def osd(counts):
    """Calculate **o**bserved, **s**ingles and **d**oubles from counts."""
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
    return singles(counts) / counts.sum()


def shannon(counts, base=2):
    """Calculate Shannon entropy of counts, default in bits."""
    freqs = counts / counts.sum()
    nonzero_freqs = freqs[freqs.nonzero()]
    return -(nonzero_freqs * np.log(nonzero_freqs)).sum() / np.log(base)


def simpson(counts):
    """Calculate Simpson's index.

    Simpson's index = 1 - dominance.

    """
    return 1 - dominance(counts)


def simpson_e(counts):
    """Calculate Simpson's evenness."""
    return enspie(counts) / observed_species(counts)


def singles(counts):
    """Return count of single occurrences."""
    return (counts == 1).sum()


def strong(counts):
    """Calculate Strong's 2002 dominance index, by way of SDR-IV."""
    n = counts.sum()
    s = observed_species(counts)
    i = np.arange(1, len(counts) + 1)
    sorted_sum = np.sort(counts)[::-1].cumsum()
    return (sorted_sum / n - (i / s)).max()


def _indices_to_counts(indices, result=None):
    """Converts vector of indices to counts of each index.

    WARNING: does not check that 'result' array is big enough to store new
    counts, suggest preallocating based on whole dataset if doing cumulative
    analysis.

    """
    if result is None:
        max_val = indices.max()
        result = np.zeros(max_val + 1)
    for i in indices:
        result[i] += 1
    return result
