#! /usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division
from future.utils.six import StringIO
from unittest import TestCase, main

import pandas as pd

from skbio.core.distance import DissimilarityMatrix, DistanceMatrix
from skbio.maths.stats.distance.base import (CategoricalStats,
                                             CategoricalStatsResults)


class CategoricalStatsTests(TestCase):
    def setUp(self):
        self.dm = DistanceMatrix([[0.0, 1.0, 2.0], [1.0, 0.0, 3.0],
                                  [2.0, 3.0, 0.0]], ['a', 'b', 'c'])
        self.grouping = [1, 2, 1]
        # Ordering of IDs shouldn't matter, nor should extra IDs.
        self.df = pd.read_csv(
            StringIO('ID,Group\nb,Group1\na,Group2\nc,Group1\nd,Group3'),
            index_col=0)
        self.df_missing_id = pd.read_csv(
            StringIO('ID,Group\nb,Group1\nc,Group1'), index_col=0)
        self.categorical_stats = CategoricalStats(self.dm, self.grouping)
        self.categorical_stats_from_df = CategoricalStats(self.dm, self.df,
                                                          column='Group')

    def test_init_invalid_input(self):
        # Requires a DistanceMatrix.
        with self.assertRaises(TypeError):
            _ = CategoricalStats(DissimilarityMatrix([[0, 2], [3, 0]],
                                                     ['a', 'b']), [1, 2])

        # Requires column if DataFrame.
        with self.assertRaises(ValueError):
            _ = CategoricalStats(self.dm, self.df)

        # Cannot provide column if not data frame.
        with self.assertRaises(ValueError):
            _ = CategoricalStats(self.dm, self.grouping, column='Group')

        # Column must exist in data frame.
        with self.assertRaises(ValueError):
            _ = CategoricalStats(self.dm, self.df, column='foo')

        # All distance matrix IDs must be in data frame.
        with self.assertRaises(ValueError):
            _ = CategoricalStats(self.dm, self.df_missing_id, column='Group')

        # Grouping vector length must match number of objects in dm.
        with self.assertRaises(ValueError):
            _ = CategoricalStats(self.dm, [1, 2])

        # Grouping vector cannot have only unique values.
        with self.assertRaises(ValueError):
            _ = CategoricalStats(self.dm, [1, 2, 3])

        # Grouping vector cannot have only a single group.
        with self.assertRaises(ValueError):
            _ = CategoricalStats(self.dm, [1, 1, 1])

    def test_call(self):
        with self.assertRaises(NotImplementedError):
            _ = self.categorical_stats()

    def test_call_invalid_permutations(self):
        with self.assertRaises(ValueError):
            _ = self.categorical_stats(-1)


class CategoricalStatsResultsTests(TestCase):

    def setUp(self):
        self.results = CategoricalStatsResults('foo', 'Foo', 'my stat', 42,
                                               ['a', 'b', 'c', 'd'],
                                               0.01234567890, 0.1151111, 99)
        self.p_value = 0.119123123123

    def test_str(self):
        exp = ('Method name  Sample size  Number of groups       my stat  '
               'p-value  Number of permutations\n        foo           42'
               '                 4  0.0123456789     0.12'
               '                      99\n')
        obs = str(self.results)
        self.assertEqual(obs, exp)

    def test_summary(self):
        exp = ('Method name\tSample size\tNumber of groups\tmy stat\tp-value\t'
               'Number of permutations\nfoo\t42\t4\t0.0123456789\t0.12\t99\n')
        obs = self.results.summary()
        self.assertEqual(obs, exp)

    def test_format_p_value(self):
        obs = self.results._format_p_value(self.p_value, 100)
        self.assertEqual(obs, '0.12')

        obs = self.results._format_p_value(self.p_value, 250)
        self.assertEqual(obs, '0.12')

        obs = self.results._format_p_value(self.p_value, 1000)
        self.assertEqual(obs, '0.119')

    def test_format_p_value_few_perms(self):
        obs = self.results._format_p_value(self.p_value, 9)
        self.assertEqual(obs, 'Too few permutations to compute p-value '
                              '(permutations = 9)')

        obs = self.results._format_p_value(self.p_value, 1)
        self.assertEqual(obs, 'Too few permutations to compute p-value '
                              '(permutations = 1)')

        obs = self.results._format_p_value(self.p_value, 0)
        self.assertEqual(obs, 'Too few permutations to compute p-value '
                              '(permutations = 0)')

    def test_format_p_value_none(self):
        obs = self.results._format_p_value(None, 0)
        self.assertEqual(obs, 'N/A')


if __name__ == '__main__':
    main()
