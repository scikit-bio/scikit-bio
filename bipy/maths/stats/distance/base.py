#! /usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The bipy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import absolute_import, division
from StringIO import StringIO
import csv

import numpy as np

from bipy.core.distance import SymmetricDistanceMatrix


class CategoricalStats(object):
    short_method_name = ''
    long_method_name = ''
    test_statistic_name = ''

    def __init__(self, distance_matrix, grouping):
        if not isinstance(distance_matrix, SymmetricDistanceMatrix):
            raise TypeError("Input must be a SymmetricDistanceMatrix.")
        if len(grouping) != distance_matrix.num_samples:
            raise ValueError("Grouping vector size must match the number of "
                             "IDs in the distance matrix.")

        # Find the group labels and convert grouping to an integer vector
        # (factor).
        groups, grouping = np.unique(grouping, return_inverse=True)

        if len(groups) == len(grouping):
            raise ValueError("All values in the grouping vector are unique. "
                             "This method cannot operate on a grouping vector "
                             "with only unique values (e.g., there are no "
                             "'within' distances because each group of "
                             "samples contains only a single sample).")
        if len(groups) == 1:
            raise ValueError("All values in the grouping vector are the same. "
                             "This method cannot operate on a grouping vector "
                             "with only a single group of samples (e.g., "
                             "there are no 'between' distances because there "
                             "is only a single group).")

        self._dm = distance_matrix
        self._grouping = grouping
        self._groups = groups
        self._tri_idxs = np.triu_indices(self._dm.num_samples, k=1)

    def __call__(self, permutations=999):
        if permutations < 0:
            raise ValueError("Number of permutations must be greater than or "
                             "equal to zero.")

        stat = self._run(self._grouping)

        p_value = None
        if permutations > 0:
            perm_stats = np.empty(permutations, dtype=np.float64)

            for i in range(permutations):
                perm_grouping = np.random.permutation(self._grouping)
                perm_stats[i] = self._run(perm_grouping)

            p_value = ((perm_stats >= stat).sum() + 1) / (permutations + 1)

        return CategoricalStatsResults(self.short_method_name,
                                       self.long_method_name,
                                       self.test_statistic_name,
                                       self._dm.num_samples, self._groups,
                                       stat, p_value, permutations)


class CategoricalStatsResults(object):

    def __init__(self, short_method_name, long_method_name,
                 test_statistic_name, sample_size, groups, statistic, p_value,
                 permutations):
        self.short_method_name = short_method_name
        self.long_method_name = long_method_name
        self.test_statistic_name = test_statistic_name
        self.sample_size = sample_size
        self.groups = groups
        self.statistic = statistic
        self.p_value = p_value
        self.permutations = permutations

    def summary(delimiter='\t'):
        """Return a formatted summary of results as a string."""
        p_value_str = self._format_p_value()

        summary = StringIO()
        csv_writer = csv.writer(summary, delimiter=delimiter)
        csv_writer.writerow(('Method name', 'Sample size',
                             'Number of groups', self.test_statistic_name,
                             'p-value', 'Number of permutations'))
        csv_writer.writerow((self.short_method_name,
                             '%d' % self.sample_size,
                             '%d' % len(self.groups), self.statistic,
                             p_value_str, '%d' % self.permutations))

        return summary.getvalue()

    def _format_p_value(self):
        """Format p-value as a string with the correct number of decimals.

        Number of decimals is determined by the number of permutations.
        """
        if self.permutations < 10:
            # This can be the last step of a long process, so we don't want to
            # fail.
            result = ('Too few permutations to compute p-value (permutations '
                      '= %d)' % self.permutations)
        elif self.p_value is None:
            result = 'N/A'
        else:
            decimal_places = int(np.log10(self.permutations + 1))
            result = ('%1.' + '%df' % decimal_places) % self.p_value

        return result
