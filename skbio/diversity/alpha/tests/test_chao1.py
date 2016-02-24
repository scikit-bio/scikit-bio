# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio.diversity.alpha import chao1, chao1_ci
from skbio.diversity.alpha._chao1 import _chao1_var


class Chao1Tests(TestCase):
    def setUp(self):
        self.counts = np.array([0, 1, 1, 4, 2, 5, 2, 4, 1, 2])
        self.no_singles = np.array([0, 2, 2, 4, 5, 0, 0, 0, 0, 0])
        self.no_doubles = np.array([0, 1, 1, 4, 5, 0, 0, 0, 0, 0])

    def test_chao1(self):
        self.assertEqual(chao1(self.counts), 9.75)
        self.assertEqual(chao1(self.counts, bias_corrected=False), 10.5)

        self.assertEqual(chao1(self.no_singles), 4)
        self.assertEqual(chao1(self.no_singles, bias_corrected=False), 4)

        self.assertEqual(chao1(self.no_doubles), 5)
        self.assertEqual(chao1(self.no_doubles, bias_corrected=False), 5)

    def test_chao1_ci(self):
        # Should match observed results from EstimateS. NOTE: EstimateS rounds
        # to 2 dp.
        obs = chao1_ci(self.counts)
        npt.assert_allclose(obs, (9.07, 17.45), rtol=0.01)

        obs = chao1_ci(self.counts, bias_corrected=False)
        npt.assert_allclose(obs, (9.17, 21.89), rtol=0.01)

        obs = chao1_ci(self.no_singles)
        npt.assert_array_almost_equal(obs, (4, 4.95), decimal=2)

        obs = chao1_ci(self.no_singles, bias_corrected=False)
        npt.assert_array_almost_equal(obs, (4, 4.95), decimal=2)

        obs = chao1_ci(self.no_doubles)
        npt.assert_array_almost_equal(obs, (4.08, 17.27), decimal=2)

        obs = chao1_ci(self.no_doubles, bias_corrected=False)
        npt.assert_array_almost_equal(obs, (4.08, 17.27), decimal=2)

    def test_chao1_var(self):
        # Should match observed results from EstimateS.NOTE: EstimateS reports
        # sd, not var, and rounds to 2 dp.
        obs = _chao1_var(self.counts)
        npt.assert_allclose(obs, 1.42 ** 2, rtol=0.01)

        obs = _chao1_var(self.counts, bias_corrected=False)
        npt.assert_allclose(obs, 2.29 ** 2, rtol=0.01)

        obs = _chao1_var(self.no_singles)
        self.assertAlmostEqual(obs, 0.39 ** 2, delta=0.01)

        obs = _chao1_var(self.no_singles, bias_corrected=False)
        self.assertAlmostEqual(obs, 0.39 ** 2, delta=0.01)

        obs = _chao1_var(self.no_doubles)
        self.assertAlmostEqual(obs, 2.17 ** 2, delta=0.01)

        obs = _chao1_var(self.no_doubles, bias_corrected=False)
        self.assertAlmostEqual(obs, 2.17 ** 2, delta=0.01)


if __name__ == '__main__':
    main()
