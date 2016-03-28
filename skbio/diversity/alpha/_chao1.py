# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from ._base import osd
from skbio.diversity._util import _validate_counts_vector
from skbio.util._decorator import experimental


@experimental(as_of="0.4.0")
def chao1(counts, bias_corrected=True):
    r"""Calculate chao1 richness estimator.

    Uses the bias-corrected version unless `bias_corrected` is ``False`` *and*
    there are both singletons and doubletons.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    bias_corrected : bool, optional
        Indicates whether or not to use the bias-corrected version of the
        equation. If ``False`` *and* there are both singletons and doubletons,
        the uncorrected version will be used. The biased-corrected version will
        be used otherwise.

    Returns
    -------
    double
        Computed chao1 richness estimator.

    See Also
    --------
    chao1_ci

    Notes
    -----
    The uncorrected version is based on Equation 6 in [1]_:

    .. math::

       chao1=S_{obs}+\frac{F_1^2}{2F_2}

    where :math:`F_1` and :math:`F_2` are the count of singletons and
    doubletons, respectively.

    The bias-corrected version is defined as

    .. math::

       chao1=S_{obs}+\frac{F_1(F_1-1)}{2(F_2+1)}

    References
    ----------
    .. [1] Chao, A. 1984. Non-parametric estimation of the number of classes in
       a population. Scandinavian Journal of Statistics 11, 265-270.

    """
    counts = _validate_counts_vector(counts)
    o, s, d = osd(counts)

    if not bias_corrected and s and d:
        return o + s ** 2 / (d * 2)
    else:
        return o + s * (s - 1) / (2 * (d + 1))


@experimental(as_of="0.4.0")
def chao1_ci(counts, bias_corrected=True, zscore=1.96):
    """Calculate chao1 confidence interval.

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    bias_corrected : bool, optional
        Indicates whether or not to use the bias-corrected version of the
        equation. If ``False`` *and* there are both singletons and doubletons,
        the uncorrected version will be used. The biased-corrected version will
        be used otherwise.
    zscore : scalar, optional
        Score to use for confidence. Default of 1.96 is for a 95% confidence
        interval.

    Returns
    -------
    tuple
        chao1 confidence interval as ``(lower_bound, upper_bound)``.

    See Also
    --------
    chao1

    Notes
    -----
    The implementation here is based on the equations in the EstimateS manual
    [1]_. Different equations are employed to calculate the chao1 variance and
    confidence interval depending on `bias_corrected` and the presence/absence
    of singletons and/or doubletons.

    Specifically, the following EstimateS equations are used:

    1. No singletons, Equation 14.
    2. Singletons but no doubletons, Equations 7, 13.
    3. Singletons and doubletons, ``bias_corrected=True``, Equations 6, 13.
    4. Singletons and doubletons, ``bias_corrected=False``, Equations 5, 13.

    References
    ----------
    .. [1] http://viceroy.eeb.uconn.edu/estimates/

    """
    counts = _validate_counts_vector(counts)
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
    r = singles / doubles
    return doubles * (.5 * r ** 2 + r ** 3 + .24 * r ** 4)


def _chao1_var_bias_corrected(s, d):
    """Calculates chao1 variance, bias-corrected.

    `s` is the number of singletons and `d` is the number of doubletons.

    From EstimateS manual, equation 6.

    """
    return (s * (s - 1) / (2 * (d + 1)) + (s * (2 * s - 1) ** 2) /
            (4 * (d + 1) ** 2) + (s ** 2 * d * (s - 1) ** 2) /
            (4 * (d + 1) ** 4))


def _chao1_var_no_doubletons(s, chao1):
    """Calculates chao1 variance in absence of doubletons.

    From EstimateS manual, equation 7.

    `s` is the number of singletons, and `chao1` is the estimate of the mean of
    Chao1 from the same dataset.

    """
    return s * (s - 1) / 2 + s * (2 * s - 1) ** 2 / 4 - s ** 4 / (4 * chao1)


def _chao1_var_no_singletons(n, o):
    """Calculates chao1 variance in absence of singletons.

    `n` is the number of individuals and `o` is the number of observed OTUs.

    From EstimateS manual, equation 8.

    """
    return o * np.exp(-n / o) * (1 - np.exp(-n / o))


def _chao_confidence_with_singletons(chao, observed, var_chao, zscore=1.96):
    """Calculates confidence bounds for chao1 or chao2.

    Uses Eq. 13 of EstimateS manual.

    `zscore` is the score to use for confidence. The default of 1.96 is for 95%
    confidence.

    """
    T = chao - observed
    # if no diff betweeh chao and observed, CI is just point estimate of
    # observed
    if T == 0:
        return observed, observed
    K = np.exp(abs(zscore) * np.sqrt(np.log(1 + (var_chao / T ** 2))))
    return observed + T / K, observed + T * K


def _chao_confidence_no_singletons(n, s, zscore=1.96):
    """Calculates confidence bounds for chao1/chao2 in absence of singletons.

    Uses Eq. 14 of EstimateS manual.

    `n` is the number of individuals and `s` is the number of OTUs.

    """
    P = np.exp(-n / s)
    return (max(s, s / (1 - P) - zscore * np.sqrt((s * P / (1 - P)))),
            s / (1 - P) + zscore * np.sqrt(s * P / (1 - P)))
