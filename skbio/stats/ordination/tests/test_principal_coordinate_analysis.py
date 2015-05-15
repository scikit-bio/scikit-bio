# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest
import warnings

import numpy as np
import numpy.testing as npt

from skbio import DistanceMatrix
from skbio.stats.distance import DissimilarityMatrixError
from skbio.stats.ordination import pcoa
from skbio.util import get_data_path


class TestPCoAResults(unittest.TestCase):
    def setUp(self):
        """Sample data set from page 111 of W.J Krzanowski. Principles
        of multivariate analysis, 2000, Oxford University Press."""
        matrix = np.loadtxt(get_data_path('PCoA_sample_data'))
        dist_matrix = DistanceMatrix(matrix, map(str, range(matrix.shape[0])))
        self.dist_matrix = dist_matrix

    def test_negative_eigenvalue_warning(self):
        """This data has some small negative eigenvalues."""
        npt.assert_warns(RuntimeWarning, pcoa, self.dist_matrix)

    def test_values(self):
        """Adapted from cogent's `test_principal_coordinate_analysis`:
        "I took the example in the book (see intro info), and did the
        principal coordinates analysis, plotted the data and it looked
        right"."""
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            scores = pcoa(self.dist_matrix)

        exp_eigvals = np.array([0.73599103, 0.26260032, 0.14926222, 0.06990457,
                                0.02956972, 0.01931184, 0., 0., 0., 0., 0., 0.,
                                0., 0.])
        exp_site = np.loadtxt(get_data_path('exp_PCoAzeros_site'))
        exp_prop_expl = np.array([0.58105792, 0.20732046, 0.1178411,
                                  0.05518899, 0.02334502, 0.01524651, 0., 0.,
                                  0., 0., 0., 0., 0., 0.])
        exp_site_ids = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                        '10', '11', '12', '13']
        # Note the absolute value because column can have signs swapped
        npt.assert_almost_equal(scores.eigvals, exp_eigvals)
        npt.assert_almost_equal(np.abs(scores.site), exp_site)
        npt.assert_almost_equal(scores.proportion_explained, exp_prop_expl)
        npt.assert_equal(scores.site_ids, exp_site_ids)


class TestPCoAEigenResults(unittest.TestCase):
    def setUp(self):
        dist_matrix = DistanceMatrix.read(get_data_path('PCoA_sample_data_3'))
        self.scores = pcoa(dist_matrix)

        self.ids = ['PC.636', 'PC.635', 'PC.356', 'PC.481', 'PC.354', 'PC.593',
                    'PC.355', 'PC.607', 'PC.634']

    def test_values(self):
        exp = OrdinationResults(
            short_method_name='PCoA',
            long_method_name='Principal Coordinate Analysis',
            eigvals=pd.Series([0.51236726, 0.30071909, 0.26791207, 0.20898868,
                               0.19169895, 0.16054235,  0.15017696,  0.12245775,
                               0.0])
        )



        results = self.scores

        npt.assert_almost_equal(len(results.eigvals), len(results.samples[0]))

        expected = np.loadtxt(get_data_path('exp_PCoAEigenResults_site'))
        npt.assert_almost_equal(*normalize_signs(expected, results.site))

        expected = np.array([0.51236726, 0.30071909, 0.26791207, 0.20898868,
                             0.19169895, 0.16054235,  0.15017696,  0.12245775,
                             0.0])
        npt.assert_almost_equal(results.eigvals, expected)

        expected = np.array([0.2675738328, 0.157044696, 0.1399118638,
                             0.1091402725, 0.1001110485, 0.0838401162,
                             0.0784269939, 0.0639511764, 0.0])
        npt.assert_almost_equal(results.proportion_explained, expected)

        npt.assert_equal(results.site_ids, self.ids)


class TestPCoAResultsExtensive(unittest.TestCase):
    def setUp(self):
        matrix = np.loadtxt(get_data_path('PCoA_sample_data_2'))
        self.ids = [str(i) for i in range(matrix.shape[0])]
        dist_matrix = DistanceMatrix(matrix, self.ids)
        self.scores = pcoa(dist_matrix)

    def test_values(self):
        results = self.scores
        npt.assert_equal(len(results.eigvals), len(results.site[0]))

        expected = np.array([[-0.028597, 0.22903853, 0.07055272,
                              0.26163576, 0.28398669, 0.0],
                             [0.37494056, 0.22334055, -0.20892914,
                              0.05057395, -0.18710366, 0.0],
                             [-0.33517593, -0.23855979, -0.3099887,
                              0.11521787, -0.05021553, 0.0],
                             [0.25412394, -0.4123464, 0.23343642,
                              0.06403168, -0.00482608, 0.0],
                             [-0.28256844, 0.18606911, 0.28875631,
                              -0.06455635, -0.21141632, 0.0],
                             [0.01727687, 0.012458, -0.07382761,
                              -0.42690292, 0.1695749, 0.0]])
        npt.assert_almost_equal(*normalize_signs(expected, results.site))

        expected = np.array([0.3984635, 0.36405689, 0.28804535, 0.27479983,
                            0.19165361, 0.0])
        npt.assert_almost_equal(results.eigvals, expected)

        expected = np.array([0.2626621381, 0.2399817314, 0.1898758748,
                             0.1811445992, 0.1263356565, 0.0])
        npt.assert_almost_equal(results.proportion_explained, expected)

        npt.assert_equal(results.site_ids, self.ids)


class TestPCoAErrors(unittest.TestCase):
    def test_input(self):
        with npt.assert_raises(DissimilarityMatrixError):
            pcoa([[1, 2], [3, 4]])


if __name__ == "__main__":
    unittest.main()
