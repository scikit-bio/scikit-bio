#! /usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The bipy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import absolute_import, division

import numpy as np
from scipy.stats import rankdata

from bipy.core.distance import SymmetricDistanceMatrix
from .base import CategoricalStatsResults


class ANOSIM(object):
    short_method_name = 'ANOSIM'
    long_method_name = 'Analysis of Similarities'

    def __init__(self, distance_matrix, grouping):
        if not isinstance(distance_matrix, SymmetricDistanceMatrix):
            raise TypeError("Input must be a SymmetricDistanceMatrix.")
        if len(grouping) != distance_matrix.num_samples:
            raise ValueError("Grouping vector size must match the number of "
                             "sample IDs in the distance matrix.")

        # Find the group labels and convert grouping to an integer vector
        # (factor).
        groups, grouping = np.unique(grouping, return_inverse=True)

        if len(groups) == len(grouping):
            raise ValueError("All values in the grouping vector are unique. "
                             "ANOSIM cannot operate on a grouping vector with "
                             "only unique values (e.g., there are no 'within' "
                             "distances because each group of samples "
                             "contains only a single sample).")
        if len(groups) == 1:
            raise ValueError("All values in the grouping vector are the same. "
                             "ANOSIM cannot operate on a grouping vector with "
                             "only a single group of samples (e.g., there are "
                             "no 'between' distances because there is only a "
                             "single group).")

        self._dm = distance_matrix
        self._divisor = self._dm.num_samples * ((self._dm.num_samples - 1) / 4)
        self._grouping = grouping
        self._groups = groups
        self._ranked_dists = rankdata(self._dm.condensed_form(),
                                      method='average')
        self._tri_idxs = np.triu_indices(self._dm.num_samples, k=1)

    def __call__(self, permutations=999):
        if permutations < 0:
            raise ValueError("Number of permutations must be greater than or "
                             "equal to zero.")

        r_stat = self._anosim(self._grouping)

        p_value = None
        if permutations > 0:
            perm_stats = np.empty(permutations, dtype=np.float64)

            for i in range(permutations):
                perm_grouping = np.random.permutation(self._grouping)
                perm_stats[i] = self._anosim(perm_grouping)

            p_value = ((perm_stats >= r_stat).sum() + 1) / (permutations + 1)

        return CategoricalStatsResults(self.short_method_name,
                                       self.long_method_name,
                                       self._dm.num_samples, self._groups,
                                       r_stat, p_value, permutations)

    def _anosim(self, grouping):
        # Create a matrix where True means that the two samples are in the same
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
        """Return ANOSIM R statistic (between -1 and +1)."""
        # within
        r_W = np.mean(self._ranked_dists[grouping_tri])
        # between
        r_B = np.mean(self._ranked_dists[np.invert(grouping_tri)])
        return (r_B - r_W) / self._divisor
