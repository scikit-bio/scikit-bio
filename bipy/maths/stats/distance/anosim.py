#! /usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The bipy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import absolute_import, division
from collections import namedtuple
from itertools import combinations

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

        grouping = np.asarray(grouping)
        groups = np.unique(grouping)

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
        # Create grouping matrix, where True means that the two samples are in
        # the same group.
        grouping_matrix = np.zeros(self._dm.shape, dtype=bool)
        for group in self._groups:
            within_indices = self._index_combinations(
                np.where(grouping == group)[0])
            grouping_matrix[within_indices] = True

        # Extract upper triangle from the grouping matrix. It is important to
        # extract the values in the same order that the distances are extracted
        # from the distance matrix (see self._ranked_dists). Extracting the
        # upper triangle (excluding the diagonal) preserves this order.
        grouping_tri = grouping_matrix[self._tri_idxs]

        return self._compute_r_stat(grouping_tri)

    def _index_combinations(self, indices):
        # Modified from http://stackoverflow.com/a/11144716
        return np.tile(indices, len(indices)), np.repeat(indices, len(indices))

    def _compute_r_stat(self, grouping_tri):
        """Return ANOSIM R statistic (between -1 and +1)."""
        # within
        r_W = np.mean(self._ranked_dists[grouping_tri])
        # between
        r_B = np.mean(self._ranked_dists[np.invert(grouping_tri)])
        return (r_B - r_W) / self._divisor
