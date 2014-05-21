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


class BIOENVTests(TestCase):
    """All results were verified with R (vegan::bioenv)."""

    def setUp(self):
        self.bv_dm_88soils_str = ["\tMT2.141698\tCA1.141704\tBB2.141659\t"
                                  "CO2.141657\tTL3.141709\tSN3.141650", "MT2.141698\t0.0\t"
                                  "0.623818643706\t0.750015427505\t0.585201193913\t0.729023583672\t"
                                  "0.622135587669", "CA1.141704\t0.623818643706\t0.0\t0.774881224555"
                                  "\t0.649822398416\t0.777203137034\t0.629507320436", "BB2.141659\t"
                                  "0.750015427505\t0.774881224555\t0.0\t0.688845424001\t0.567470311282"
                                  "\t0.721707516043", "CO2.141657\t0.585201193913\t0.649822398416\t"
                                  "0.688845424001\t0.0\t0.658853575764\t0.661223617505", "TL3.141709\t"
                                  "0.729023583672\t0.777203137034\t0.567470311282\t0.658853575764\t0.0\t"
                                  "0.711173405838", "SN3.141650\t0.622135587669\t0.629507320436\t"
                                  "0.721707516043\t0.661223617505\t0.711173405838\t0.0"]
        self.bv_dm_88soils = DistanceMatrix.from_file(self.bv_dm_88soils_str)

        self.bv_map_88soils_str = StringIO('\n'.join(["#SampleId\tTOT_ORG_CARB\tSILT_CLAY\t"
                                   "ELEVATION\tSOIL_MOISTURE_DEFICIT\tCARB_NITRO_RATIO\t"
                                   "ANNUAL_SEASON_TEMP\tANNUAL_SEASON_PRECPT\tPH\tCMIN_RATE\tLONGITUDE\t"
                                   "LATITUDE", "MT2.141698\t39.1\t35\t1000\t70\t23.087\t7\t450\t6.66\t"
                                   "19.7\t-114\t46.8", "CA1.141704\t16.7\t73\t2003\t198\t13\t10.3\t400\t"
                                   "7.27\t2.276\t-111.7666667\t36.05", "BB2.141659\t52.2\t44\t400\t-680\t"
                                   "21.4\t6.1\t1200\t4.6\t2.223\t-68.1\t44.86666667", "CO2.141657\t18.1\t"
                                   "24\t2400\t104\t31.8\t6.1\t350\t5.68\t9.223\t-105.3333333\t"
                                   "40.58333333", "TL3.141709\t53.9\t52\t894\t-212\t24.6\t-9.3\t400\t"
                                   "4.23\t16.456\t-149.5833333\t68.63333333", "SN3.141650\t16.6\t20\t"
                                   "3000\t-252\t13.9\t3.6\t600\t5.74\t6.289\t-118.1666667\t36.45"]))
        self.bv_map_88soils = pd.read_csv(self.bv_map_88soils_str, sep='\t', index_col=0)

        self.cats = ['TOT_ORG_CARB', 'SILT_CLAY', 'ELEVATION',
                     'SOIL_MOISTURE_DEFICIT', 'CARB_NITRO_RATIO',
                     'ANNUAL_SEASON_TEMP', 'ANNUAL_SEASON_PRECPT', 'PH',
                     'CMIN_RATE', 'LONGITUDE', 'LATITUDE']

    def test_bioenv(self):
        obs = bioenv(self.bv_dm_88soils, self.bv_map_88soils, self.cats)

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
