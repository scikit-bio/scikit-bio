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

import numpy as np
import numpy.testing as npt
import pandas as pd

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

    def test_bioenv(self):
        obs = bioenv(self.dm_88_soils, self.df_88_soils, self.cols)

        exp = {'method_name': 'BEST',
               'rho_vals': [
                   (0.75, '8'),
                   (0.5, '1,11'),
                   (0.5107142857142857, '1,8,11'),
                   (0.47857142857142854, '1,6,8,11'),
                   (0.46071428571428574, '1,5,6,8,9'),
                   (0.4357142857142857, '1,5,6,8,9,11'),
                   (0.38928571428571423, '1,2,4,5,6,7,8'),
                   (0.38928571428571423, '1,2,4,5,6,7,8,9'),
                   (0.38928571428571423, '1,2,4,5,6,7,8,9,10'),
                   (0.38928571428571423, '1,2,4,5,6,7,8,9,10,11'),
                   (0.16785714285714282, '1,2,3,4,5,6,7,8,9,10,11')],
               'num_vars': 11,
               'vars': ['TOT_ORG_CARB = 1', 'SILT_CLAY = 2',
                        'ELEVATION = 3', 'SOIL_MOISTURE_DEFICIT = 4',
                        'CARB_NITRO_RATIO = 5', 'ANNUAL_SEASON_TEMP = 6',
                        'ANNUAL_SEASON_PRECPT = 7', 'PH = 8',
                        'CMIN_RATE = 9', 'LONGITUDE = 10',
                        'LATITUDE = 11']}

        # the difference caused by the new spearmans_rho function is in the
        # final decimal place of the rho_vals. its not a significant difference,
        # but the test is failing because we are doing a dict comparison for
        # float values.
        self.assertEqual(exp['method_name'], obs['method_name'])
        self.assertEqual(exp['num_vars'], obs['num_vars'])
        self.assertEqual(exp['vars'], obs['vars'])
        for i, j in zip(exp['rho_vals'], obs['rho_vals']):
            self.assertAlmostEqual(i[0], j[0])
            self.assertEqual(i[1], j[1])
        # check that the keys are the same since we have checked all the values
        # we expect to be there are the same
        self.assertTrue(set(exp.keys()) == set(obs.keys()))


if __name__ == '__main__':
    main()
