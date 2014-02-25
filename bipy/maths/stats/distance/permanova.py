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


class PERMANOVA(object):
    short_method_name = 'PERMANOVA'
    long_method_name = 'Permutational Multivariate Analysis of Variance'

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
                             "PERMANOVA cannot operate on a grouping vector with "
                             "only unique values (e.g., there are no 'within' "
                             "distances because each group of samples "
                             "contains only a single sample).")
        if len(groups) == 1:
            raise ValueError("All values in the grouping vector are the same. "
                             "PERMANOVA cannot operate on a grouping vector with "
                             "only a single group of samples (e.g., there are "
                             "no 'between' distances because there is only a "
                             "single group).")

        self._dm = distance_matrix
        self._grouping = grouping
        self._groups = groups
        self._num_groups = len(self._groups)
        self._distances = self._dm.condensed_form()
        self._s_T = (self._distances ** 2).sum() / self._dm.num_samples
        self._tri_idxs = np.triu_indices(self._dm.num_samples, k=1)

    def __call__(self, permutations=999):
        if permutations < 0:
            raise ValueError("Number of permutations must be greater than or "
                             "equal to zero.")

        f_stat = self._permanova(self._grouping)

        p_value = None
        if permutations > 0:
            perm_stats = np.empty(permutations, dtype=np.float64)

            for i in range(permutations):
                perm_grouping = np.random.permutation(self._grouping)
                perm_stats[i] = self._permanova(perm_grouping)

            p_value = ((perm_stats >= f_stat).sum() + 1) / (permutations + 1)

        return CategoricalStatsResults(self.short_method_name,
                                       self.long_method_name,
                                       self._dm.num_samples, self._groups,
                                       f_stat, p_value, permutations)

    def _permanova(self, grouping):
        """Compute PERMANOVA pseudo-F statistic."""
        # Create grouping matrix.
        grouping_matrix = -1 * np.ones(self._dm.shape, dtype=int)
        for group_idx, group in enumerate(self._groups):
            within_indices = self._index_combinations(
                np.where(grouping == group)[0])
            grouping_matrix[within_indices] = group_idx

        # Extract upper triangle.
        grouping_tri = grouping_matrix[self._tri_idxs]

        # Calculate number of samples in each group.
        unique_n = np.bincount(np.unique(grouping, return_inverse=True)[1])

        return self._compute_f_stat(grouping_tri, unique_n)

    def _index_combinations(self, indices):
        # Modified from http://stackoverflow.com/a/11144716
        return np.tile(indices, len(indices)), np.repeat(indices, len(indices))

    def _compute_f_stat(self, grouping_tri, unique_n):
        a = self._num_groups
        N = self._dm.num_samples

        # Calculate s_W for each group, accounting for different group sizes.
        s_W = 0
        for i in range(a):
            s_W += (self._distances[grouping_tri == i] ** 2).sum() / unique_n[i]

        s_A = self._s_T - s_W
        return (s_A / (a - 1)) / (s_W / (N - a))
