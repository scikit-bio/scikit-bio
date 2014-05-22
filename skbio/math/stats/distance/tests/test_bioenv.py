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
from pandas.util.testing import assert_frame_equal

from skbio.core.distance import DistanceMatrix
from skbio.math.stats.distance.bioenv import bioenv
from skbio.util.testing import get_data_path


class BIOENVTests(TestCase):
    """All results were verified with R (vegan::bioenv)."""

    def setUp(self):
        self.dm_88_soils = DistanceMatrix.from_file(
            get_data_path('88_soils_dm.txt'))
        self.df_88_soils = pd.read_csv(get_data_path('88_soils_df.txt'),
                                       sep='\t', index_col=0)
        self.cols = ['TOT_ORG_CARB', 'SILT_CLAY', 'ELEVATION',
                     'SOIL_MOISTURE_DEFICIT', 'CARB_NITRO_RATIO',
                     'ANNUAL_SEASON_TEMP', 'ANNUAL_SEASON_PRECPT', 'PH',
                     'CMIN_RATE', 'LONGITUDE', 'LATITUDE']
        self.exp_88_soils = pd.read_csv(
            get_data_path('88_soils_exp_results.txt'), sep='\t', index_col=0)
        self.exp_88_soils_single_column = pd.read_csv(
            get_data_path('88_soils_single_column_exp_results.txt'), sep='\t',
            index_col=0)

    def test_bioenv_all_columns_implicit(self):
        # Test with all columns (implicitly, not specified).
        obs = bioenv(self.dm_88_soils, self.df_88_soils)
        #obs.to_csv('tests/data/88_soils_exp_results.txt', sep='\t')
        assert_frame_equal(obs, self.exp_88_soils)

    def test_bioenv_all_columns_explicit(self):
        # Test with all columns being specified.
        obs = bioenv(self.dm_88_soils, self.df_88_soils, columns=self.cols)
        assert_frame_equal(obs, self.exp_88_soils)

    def test_bioenv_single_column(self):
        obs = bioenv(self.dm_88_soils, self.df_88_soils, columns=['PH'])
        assert_frame_equal(obs, self.exp_88_soils_single_column)

    def test_bioenv_no_distance_matrix(self):
        with self.assertRaises(TypeError):
            bioenv('breh', self.df_88_soils)

    def test_bioenv_no_data_frame(self):
        with self.assertRaises(TypeError):
            bioenv(self.dm_88_soils, None)

    def test_bioenv_duplicate_columns(self):
        with self.assertRaises(ValueError):
            bioenv(self.dm_88_soils, self.df_88_soils,
                   columns=self.cols + ['PH'])

    def test_bioenv_no_columns(self):
        with self.assertRaises(ValueError):
            bioenv(self.dm_88_soils, self.df_88_soils, columns=[])


if __name__ == '__main__':
    main()
