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
        dm_size = self._dm.num_samples
        distances = self._dm.data[np.triu_indices(dm_size, k=1)]
        grouping = within_between[np.triu_indices(dm_size, k=1)]
        #distances = self._dm.data[np.tri(dm_size) == 0]
        #grouping = within_between[np.tri(dm_size) == 0]

        # Sort extracted data.
        sorted_distances = []
        sorted_grouping = []
        for idx in np.argsort(distances):
            sorted_distances.append(distances[idx])
            sorted_grouping.append(grouping[idx])

        # Account for rank ties, then compute R statistic.
        rank_list = range(1, len(sorted_distances) + 1)
        adjusted_rank_list = self._remove_ties(sorted_distances, rank_list)
        return self._compute_r_stat(adjusted_rank_list, sorted_grouping,
                                    dm_size)

    def _remove_ties(self, sorted_dists, ranks):
        """Replaces repeat values with the average of them.

        Returns a list containing the adjusted ranks.

        Arguments:
            sorted_dists: list of the sorted distances
            ranks: list containing the ranks of each of the distances
        """
        result = []
        ties = []
        tie_count = 0
        tie_flag = 0

        for i in range(len(sorted_dists) - 1):
            # Store state information.
            curr_dist = sorted_dists[i]
            next_dist = sorted_dists[i + 1]
            rank_val = ranks[i]

            # A tie has not occured yet.
            if tie_flag == 0:
                if curr_dist == next_dist:
                    # We have a tie, so add the current rank to the tie list.
                    tie_count = tie_count + 1
                    ties.append(rank_val)
                    first_tie_index = i
                    tie_flag = 1
                else:
                    # If no tie, fill in the list with the current rank.
                    result.append(rank_val)
            else:
                # A tie has already occured.
                if curr_dist == next_dist:
                    # If another tie occurs, add the current rank to the tie
                    # list.
                    tie_count = tie_count + 1
                    ties.append(rank_val)
                else:
                    # No more ties, average their values and attach to adjusted
                    # rank list.
                    ties.append(rank_val)
                    last_tie_index = i
                    result.extend(self._get_adjusted_vals(ties,
                                                          first_tie_index, last_tie_index))
                    tie_flag = 0
                    tie_count = 0
                    ties = []
        # If there is a tie that extends to the final position, we must process
        # it here to avoid out of list bounds errors.
        if tie_flag == 1:
            ties.append(ranks[i + 1])
            last_tie_index = i + 1
            result.extend(self._get_adjusted_vals(ties, first_tie_index,
                                                  last_tie_index))
        else:
            result.append(ranks[i + 1])
        return result

    def _get_adjusted_vals(self, ties, first_tie_idx, last_tie_idx):
        """Helper function to _remove_ties. Consolidates repeated code."""
        adjusted_val = sum(ties) / len(ties)
        return [adjusted_val] * ((last_tie_idx - first_tie_idx) + 1)

    def _compute_r_stat(self, adjusted_ranks, sorted_groups, num_samps):
        """Code that performs the actual math involved in solving ANOSIM.

        Returns the ANOSIM R value (between -1 and 1).

        Arguments:
            adjusted_ranks - list of the ranks, adjusted for ties
            sorted_groups - list associating distances to groups
            num_samps: how many total samples
        """
        adjusted_ranks = np.array(adjusted_ranks)
        sorted_groups = np.array(sorted_groups)

        # Compute r_W and r_B.
        r_W = np.mean(adjusted_ranks[sorted_groups == 1])
        r_B = np.mean(adjusted_ranks[sorted_groups == 0])
        divisor = num_samps * ((num_samps - 1) / 4)
        return (r_B - r_W) / divisor
