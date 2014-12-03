#!/usr/bin/env python
from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np

from ._base import _validate


def ace(counts, rare_threshold=10):
    """Calculate the ACE metric (Abundance-based Coverage Estimator).

    Parameters
    ----------
    counts : 1-D array_like, int
        Vector of counts.
    rare_threshold : int, optional
        Threshold at which an OTU containing as many or fewer individuals will
        be considered rare.

    Returns
    -------
    double
        Computed ACE metric.

    Raises
    ------
    ValueError
        If every rare OTU is a singleton.

    Notes
    -----
    ACE was first introduced in [1]_ and [2]_. The implementation here is based
    on the description given in the EstimateS manual [3]_.

    If no rare OTUs exist, returns the number of abundant OTUs. The default
    value of 10 for `rare_threshold` is based on [4]_.

    If `counts` contains zeros, indicating OTUs which are known to exist in the
    environment but did not appear in the sample, they will be ignored for the
    purpose of calculating the number of rare OTUs.

    References
    ----------
    .. [1] Chao, A. & S.-M Lee. 1992 Estimating the number of classes via
       sample coverage. Journal of the American Statistical Association 87,
       210-217.
    .. [2] Chao, A., M.-C. Ma, & M. C. K. Yang. 1993. Stopping rules and
       estimation for recapture debugging with unequal failure rates.
       Biometrika 80, 193-201.
    .. [3] http://viceroy.eeb.uconn.edu/estimates/
    .. [4] Chao, A., W.-H. Hwang, Y.-C. Chen, and C.-Y. Kuo. 2000. Estimating
       the number of shared species in two communities. Statistica Sinica
       10:227-246.

    """
    counts = _validate(counts)
    freq_counts = np.bincount(counts)
    s_rare = _otus_rare(freq_counts, rare_threshold)
    singles = freq_counts[1]

    if singles > 0 and singles == s_rare:
        raise ValueError("The only rare OTUs are singletons, so the ACE "
                         "metric is undefined. EstimateS suggests using "
                         "bias-corrected Chao1 instead.")

    s_abun = _otus_abundant(freq_counts, rare_threshold)
    if s_rare == 0:
        return s_abun

    n_rare = _number_rare(freq_counts, rare_threshold)
    c_ace = 1 - singles / n_rare

    top = s_rare * _number_rare(freq_counts, rare_threshold, gamma=True)
    bottom = c_ace * n_rare * (n_rare - 1)
    gamma_ace = (top / bottom) - 1

    if gamma_ace < 0:
        gamma_ace = 0

    return s_abun + (s_rare / c_ace) + ((singles / c_ace) * gamma_ace)


def _otus_rare(freq_counts, rare_threshold):
    """Count number of rare OTUs."""
    return freq_counts[1:rare_threshold + 1].sum()


def _otus_abundant(freq_counts, rare_threshold):
    """Count number of abundant OTUs."""
    return freq_counts[rare_threshold + 1:].sum()


def _number_rare(freq_counts, rare_threshold, gamma=False):
    """Return number of individuals in rare OTUs.

    ``gamma=True`` generates the ``n_rare`` used for the variation coefficient.

    """
    n_rare = 0

    if gamma:
        for i, j in enumerate(freq_counts[:rare_threshold + 1]):
            n_rare = n_rare + (i * j) * (i - 1)
    else:
        for i, j in enumerate(freq_counts[:rare_threshold + 1]):
            n_rare = n_rare + (i * j)

    return n_rare
