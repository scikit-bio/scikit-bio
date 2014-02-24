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


class CategoricalStatsResults(object):

    def __init__(self, short_method_name, long_method_name, sample_size,
                 groups, statistic, p_value, permutations):
        self.short_method_name = short_method_name
        self.long_method_name = long_method_name
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
                             'Number of groups', 'Test statistic',
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
