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

from .base import _indices_to_counts, _validate


def ace(counts, rare_threshold=10):
    """Implements the ACE metric from EstimateS. Based on the equations
    given under ACE:Abundance-based Coverage Estimator.

    counts = an OTU by sample vector
    rare_threshold = threshold at which a species containing as many or
    fewer individuals will be considered rare.

    IMPORTANT NOTES:

    Raises a value error if every rare species is a singleton.

    if no rare species exist, just returns the number of abundant species

    rare_threshold default value is 10. Based on Chao 2000 in Statistica
    Sinica pg. 229 citing empirical observations by Chao, Ma, Yang 1993.

    If counts contains 0's, indicating species which are known
    to exist in the environment but did not appear in the sample, they
    will be ignored for the purpose of calculating s_rare.

    """
    counts = _validate(counts)
    freq_counts = _indices_to_counts(counts)

    if freq_counts[1:rare_threshold].sum() == 0:
        return _species_abundant(freq_counts, rare_threshold)

    if freq_counts[1] == freq_counts[1:rare_threshold].sum():
        raise ValueError("The only rare species are singletons, so the ACE "
                         "metric is undefined. EstimateS suggests using "
                         "bias-corrected Chao1 instead.")

    s_abun = _species_abundant(freq_counts, rare_threshold)
    s_rare = freq_counts[1:rare_threshold + 1].sum()
    n_rare = _number_rare(freq_counts, rare_threshold)
    c_ace = 1 - (freq_counts[1]).sum() / float(n_rare)

    top = s_rare * _number_rare(freq_counts, rare_threshold, gamma=True)
    bottom = c_ace * n_rare * (n_rare - 1.0)
    gamma_ace = (top / bottom) - 1.0

    if gamma_ace < 0:
        gamma_ace = 0

    return s_abun + (s_rare / c_ace) + ((freq_counts[1] / c_ace) * gamma_ace)


def _species_abundant(freq_counts, rare_threshold):
    """Count number of abundant species."""
    return freq_counts[rare_threshold + 1:].sum()


def _number_rare(freq_counts, rare_threshold, gamma=False):
    """Return number of individuals in rare species.

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
