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
from itertools import combinations

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
        self._divisor = self._dm.num_samples * ((self._dm.num_samples - 1) / 4)
        self._grouping = np.asarray(grouping)
        self._groups = np.unique(self._grouping)
        self._ranked_dists = rankdata(self._dm.condensed_form(),
                                      method='average')
        self._tri_idxs = np.triu_indices(self._dm.num_samples, k=1)

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
        within_between = np.zeros(self._dm.shape, dtype=bool)
        for group in self._groups:
            combs = self._cartesian(np.where(grouping == group)[0])
            within_between[combs] = True

        # Extract triangle from the distance and grouping matrices. TODO: add
        # note about importance of extraction order (must match condensed dm
        # form).
        grouping_tri = within_between[self._tri_idxs]

        return self._compute_r_stat(grouping_tri)

    def _cartesian(self, a):
        # Modified from http://stackoverflow.com/a/11144716
        return np.tile(a, len(a)), np.repeat(a, len(a))

    def _compute_r_stat(self, grouping_upper):
        """Code that performs the actual math involved in solving ANOSIM.

        Returns the ANOSIM R value (between -1 and 1).

        Arguments:
            adjusted_ranks - list of the ranks, adjusted for ties
            sorted_groups - list associating distances to groups
            num_samps: how many total samples
        """
        r_W = np.mean(self._ranked_dists[grouping_upper])
        r_B = np.mean(self._ranked_dists[np.invert(grouping_upper)])
        return (r_B - r_W) / self._divisor
