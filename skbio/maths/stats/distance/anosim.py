#! /usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import absolute_import, division

import numpy as np
from scipy.stats import rankdata

from .base import CategoricalStats


class ANOSIM(CategoricalStats):
    """ANOSIM statistical method executor.

    Analysis of Similarities (ANOSIM) is a non-parametric method that tests
    whether two or more groups of objects are significantly different based on
    a categorical factor. The ranks of the distances in the distance matrix are
    used to calculate an R statistic, which ranges between -1 (anti-grouping)
    to +1 (strong grouping), with an R value of 0 indicating random grouping.

    Notes
    -----
    See [1]_ for the original ANOSIM reference. The general algorithm and
    interface are similar to ``vegan::anosim``, available in R's vegan package
    [2]_.

    References
    ----------
    .. [1] Clarke, KR. "Non-parametric multivariate analyses of changes in
       community structure." Australian journal of ecology 18.1 (1993):
       117-143.

    .. [2] http://cran.r-project.org/web/packages/vegan/index.html

    """

    short_method_name = 'ANOSIM'
    long_method_name = 'Analysis of Similarities'
    test_statistic_name = 'R statistic'

    def __init__(self, distance_matrix, grouping):
        super(ANOSIM, self).__init__(distance_matrix, grouping)

        self._divisor = self._dm.shape[0] * ((self._dm.shape[0] - 1) / 4)
        self._ranked_dists = rankdata(self._dm.condensed_form(),
                                      method='average')

    def _run(self, grouping):
        """Compute ANOSIM R statistic (between -1 and +1)."""
        # Create a matrix where True means that the two objects are in the same
        # group. This ufunc requires that grouping is a numeric vector (e.g.,
        # it won't work with a grouping vector of strings).
        grouping_matrix = np.equal.outer(grouping, grouping)

        # Extract upper triangle from the grouping matrix. It is important to
        # extract the values in the same order that the distances are extracted
        # from the distance matrix (see self._ranked_dists). Extracting the
        # upper triangle (excluding the diagonal) preserves this order.
        grouping_tri = grouping_matrix[self._tri_idxs]

        return self._compute_r_stat(grouping_tri)

    def _compute_r_stat(self, grouping_tri):
        # within
        r_W = np.mean(self._ranked_dists[grouping_tri])
        # between
        r_B = np.mean(self._ranked_dists[np.invert(grouping_tri)])
        return (r_B - r_W) / self._divisor
