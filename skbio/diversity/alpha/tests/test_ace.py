# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np
import numpy.testing as npt

from skbio.diversity import alpha_diversity
from skbio.diversity.alpha import ace


class AceTests(unittest.TestCase):
    def test_ace(self):
        npt.assert_almost_equal(ace(np.array([2, 0])), 1.0)
        npt.assert_almost_equal(ace(np.array([12, 0, 9])), 2.0)
        npt.assert_almost_equal(ace(np.array([12, 2, 8])), 3.0)
        npt.assert_almost_equal(ace(np.array([12, 2, 1])), 4.0)
        npt.assert_almost_equal(ace(np.array([12, 1, 2, 1])), 7.0)
        npt.assert_almost_equal(ace(np.array([12, 3, 2, 1])), 4.6)
        npt.assert_almost_equal(ace(np.array([12, 3, 6, 1, 10])), 5.62749672)

        # Just return the number of taxa when all are abundant.
        npt.assert_almost_equal(ace(np.array([12, 12, 13, 14])), 4.0)

        # Border case: only singletons and 10-tons, no abundant taxa.
        npt.assert_almost_equal(ace([0, 1, 1, 0, 0, 10, 10, 1, 0, 0]),
                                9.35681818182)

    def test_ace_empty_sample(self):
        # A sample with no observed taxon has nothing to estimate.
        npt.assert_almost_equal(ace(np.array([0])), 0.0)
        npt.assert_almost_equal(ace(np.array([0, 0, 0])), 0.0)

        # An empty sample must not abort the calculation for the whole table.
        counts = np.array([[12, 3, 2, 1], [0, 0, 0, 0], [12, 12, 13, 14]])
        obs = alpha_diversity("ace", counts)
        npt.assert_almost_equal(obs.values, [4.6, 0.0, 4.0])

    def test_ace_only_rare_singletons(self):
        with self.assertRaises(ValueError):
            ace([0, 0, 43, 0, 1, 0, 1, 42, 1, 43])


if __name__ == '__main__':
    unittest.main()
