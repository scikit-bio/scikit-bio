#! /usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import absolute_import, division

import numpy as np

from .base import CategoricalStats


class PERMANOVA(CategoricalStats):
    """PERMANOVA statistical method executor.

    Permutational Multivariate Analysis of Variance (PERMANOVA) is a
    non-parametric method that tests whether two or more groups of objects are
    significantly different based on a categorical factor. It is conceptually
    similar to ANOVA except that it operates on a distance matrix, which allows
    for multivariate analysis. PERMANOVA computes a pseudo-F statistic and
    tests the significance through permutations.

    Notes
    -----
    See [1]_ for the original PERMANOVA reference, as well as
    ``vegan::adonis``, available in R's vegan package [2]_.

    References
    ----------
    .. [1] Anderson, Marti J. "A new method for non-parametric multivariate
       analysis of variance." Austral Ecology 26.1 (2001): 32-46.

    .. [2] http://cran.r-project.org/web/packages/vegan/index.html

    """

    short_method_name = 'PERMANOVA'
    long_method_name = 'Permutational Multivariate Analysis of Variance'
    test_statistic_name = 'pseudo-F statistic'

    def __init__(self, distance_matrix, grouping):
        super(PERMANOVA, self).__init__(distance_matrix, grouping)

        # Calculate number of objects in each group.
        self._group_sizes = np.bincount(self._grouping)
        self._num_groups = len(self._groups)
        self._distances = self._dm.condensed_form()
        self._s_T = (self._distances ** 2).sum() / self._dm.shape[0]

    def _run(self, grouping):
        """Compute PERMANOVA pseudo-F statistic."""
        # Create a matrix where objects in the same group are marked with the
        # group index (e.g. 0, 1, 2, etc.). objects that are not in the same
        # group are marked with -1.
        grouping_matrix = -1 * np.ones(self._dm.shape, dtype=int)
        for group_idx in range(len(self._groups)):
            within_indices = self._index_combinations(
                np.where(grouping == group_idx)[0])
            grouping_matrix[within_indices] = group_idx

        # Extract upper triangle (in same order as distances were extracted
        # from full distance matrix).
        grouping_tri = grouping_matrix[self._tri_idxs]

        return self._compute_f_stat(grouping_tri)

    def _index_combinations(self, indices):
        # Modified from http://stackoverflow.com/a/11144716
        return np.tile(indices, len(indices)), np.repeat(indices, len(indices))

    def _compute_f_stat(self, grouping_tri):
        a = self._num_groups
        N = self._dm.shape[0]

        # Calculate s_W for each group, accounting for different group sizes.
        s_W = 0
        for i in range(a):
            s_W += ((self._distances[grouping_tri == i] ** 2).sum() /
                    self._group_sizes[i])

        s_A = self._s_T - s_W
        return (s_A / (a - 1)) / (s_W / (N - a))
