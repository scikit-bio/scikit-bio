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

from .base import _indices_to_counts, _validate, osd


def chao1(counts, bias_corrected=True):
    """Calculates chao1 according to table in EstimateS manual.

    Specifically, uses bias-corrected version unless bias_corrected is set
    to False _and_ there are both singletons and doubletons.

    Uncorrected:

    Calculates chao1 given counts. Eq. 1 in EstimateS manual.

    Formula: chao1 = S_obs + N_1^2/(2*N_2) where N_1 and N_2 are
    count of singletons and doubletons respectively.

    Note: this is the original formula from Chao 1984, not bias-corrected,
    and is Equation 1 in the EstimateS manual.

    Corrected:

    Calculates bias-corrected chao1 given counts: Eq. 2 in EstimateS manual.

    Formula: chao1 = S_obs + N_1(N_1-1)/(2*(N_2+1)) where N_1 and N_2 are
    count of singletons and doubletons respectively.

    Note: this is the bias-corrected formulat from Chao 1987, Eq. 2 in the
    EstimateS manual.

    """
    o, s, d = osd(counts)

    if not bias_corrected and s and d:
        return o + s ** 2 / (d * 2)
    else:
        return o + s * (s - 1) / (2 * (d + 1))


def chao1_confidence(counts, bias_corrected=True, zscore=1.96):
    """Returns chao1 confidence (lower, upper) from counts."""
    o, s, d = osd(counts)
    if s:
        chao = chao1(counts, bias_corrected)
        chaovar = _chao1_var(counts, bias_corrected)
        return _chao_confidence_with_singletons(chao, o, chaovar, zscore)
    else:
        n = counts.sum()
        return _chao_confidence_no_singletons(n, o, zscore)


def _chao1_var(counts, bias_corrected=True):
    """Calculates chao1 variance using decision rules in EstimateS."""
    o, s, d = osd(counts)
    if not d:
        c = chao1(counts, bias_corrected)
        return _chao1_var_no_doubletons(s, c)
    if not s:
        n = counts.sum()
        return _chao1_var_no_singletons(n, o)
    if bias_corrected:
        return _chao1_var_bias_corrected(s, d)
    else:
        return _chao1_var_uncorrected(s, d)


def _chao1_var_uncorrected(singles, doubles):
    """Calculates chao1, uncorrected.

    From EstimateS manual, equation 5.
    """
    r = float(singles) / doubles
    return doubles * (.5 * r ** 2 + r ** 3 + .24 * r ** 4)


def _chao1_var_bias_corrected(singles, doubles):
    """Calculates chao1 variance, bias-corrected.

    From EstimateS manual, equation 6.
    """
    s, d = float(singles), float(doubles)
    return s * (s - 1) / (2 * (d + 1)) + (s * (2 * s - 1) ** 2) / (4 * (d + 1) ** 2) + \
        (s ** 2 * d * (s - 1) ** 2) / (4 * (d + 1) ** 4)


def _chao1_var_no_doubletons(singles, chao1):
    """Calculates chao1 variance in absence of doubletons.

    From EstimateS manual, equation 7.

    chao1 is the estimate of the mean of Chao1 from the same dataset.
    """
    s = float(singles)
    return s * (s - 1) / 2 + s * (2 * s - 1) ** 2 / 4 - s ** 4 / (4 * chao1)


def _chao1_var_no_singletons(n, observed):
    """Calculates chao1 variance in absence of singletons. n = # individuals.

    From EstimateS manual, equation 8.
    """
    o = float(observed)
    return o * np.exp(-n / o) * (1 - np.exp(-n / o))


def _chao_confidence_with_singletons(chao, observed, var_chao, zscore=1.96):
    """Calculates confidence bounds for chao1 or chao2.

    Uses Eq. 13 of EstimateS manual.

    zscore = score to use for confidence, default = 1.96 for 95% confidence.
    """
    T = float(chao - observed)
    # if no diff betweeh chao and observed, CI is just point estimate of
    # observed
    if T == 0:
        return observed, observed
    K = np.exp(abs(zscore) * np.sqrt(np.log(1 + (var_chao / T ** 2))))
    return observed + T / K, observed + T * K


def _chao_confidence_no_singletons(n, observed, zscore=1.96):
    """Calculates confidence bounds for chao1/chao2 in absence of singletons.

    Uses Eq. 14 of EstimateS manual.

    n = number of individuals, observed = number of species.
    """
    s = float(observed)
    P = np.exp(-n / s)
    return max(s, s / (1 - P) - zscore * np.sqrt((s * P / (1 - P)))), \
        s / (1 - P) + zscore * np.sqrt(s * P / (1 - P))
