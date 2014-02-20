#! /usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from collections import namedtuple

import numpy as np
from scipy.stats import rankdata

from bipy.core.distance import SymmetricDistanceMatrix


# TODO: store sample size and number of groups
ANOSIMResults = namedtuple('ANOSIMResults', ('short_method_name',
                                             'long_method_name',
                                             'r_statistic', 'p_value',
                                             'permutations'))


class ANOSIM(object):
    short_method_name = 'ANOSIM'
    long_method_name = 'Analysis of Similarities'

    def __init__(self, distance_matrix, grouping):
        if not isinstance(distance_matrix, SymmetricDistanceMatrix):
            raise TypeError("Input must be a SymmetricDistanceMatrix.")
        if len(grouping) != distance_matrix.num_samples:
            raise ValueError("Grouping vector size must match the number of "
                             "sample IDs in the distance matrix.")

        # TODO: test for uniqueness and single-value-only in grouping vector
        self._dm = distance_matrix
        self._grouping = grouping

    def __call__(self, permutations=999):
        # TODO: test for invalid number of permutations
        r_stat = self._anosim(self._grouping)

        p_value = None
        if permutations > 0:
            perm_stats = np.empty(permutations, dtype=np.float64)

            for i in range(permutations):
                perm_grouping = np.random.permutation(self._grouping)
                perm_stats[i] = self._anosim(perm_grouping)

            p_value = ((perm_stats >= r_stat).sum() + 1) / (permutations + 1)

        return ANOSIMResults(self.short_method_name, self.long_method_name,
                             r_stat, p_value, permutations)

    def _anosim(self, grouping):
        # Create grouping matrix, where a one means that the two samples are in
        # the same group (e.g. control) and a zero means that they aren't.
        within_between = np.zeros(self._dm.shape, dtype=np.int32)
        for i, i_value in enumerate(grouping):
            for j, j_value in enumerate(grouping):
                if i_value == j_value:
                    within_between[i][j] = 1

        # Extract upper triangle from the distance and grouping matrices.
        num_samples = self._dm.num_samples
        ranked_distances = rankdata(self._dm.condensed_form(),
                                    method='average')
        grouping_upper = within_between[np.triu_indices(num_samples, k=1)]

        return self._compute_r_stat(ranked_distances, grouping_upper,
                                    num_samples)

    def _compute_r_stat(self, ranked_distances, grouping_upper, num_samples):
        """Code that performs the actual math involved in solving ANOSIM.

        Returns the ANOSIM R value (between -1 and 1).

        Arguments:
            adjusted_ranks - list of the ranks, adjusted for ties
            sorted_groups - list associating distances to groups
            num_samps: how many total samples
        """
        r_W = np.mean(ranked_distances[grouping_upper == 1])
        r_B = np.mean(ranked_distances[grouping_upper == 0])
        divisor = num_samples * ((num_samples - 1) / 4)
        return (r_B - r_W) / divisor
