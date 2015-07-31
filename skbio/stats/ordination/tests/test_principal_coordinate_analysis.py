# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import pandas as pd
import numpy as np
import numpy.testing as npt
from unittest import TestCase, main

from skbio import DistanceMatrix, OrdinationResults
from skbio.stats.distance import DissimilarityMatrixError
from skbio.stats.ordination import pcoa
from skbio.util import get_data_path, assert_ordination_results_equal


class TestPCoA(TestCase):
    def setUp(self):
        # Sample data set from page 111 of W.J Krzanowski. Principles
        # of multivariate analysis, 2000, Oxford University Press.
        self.dm = DistanceMatrix(
            np.loadtxt(get_data_path('PCoA_sample_data')))

    def test_simple(self):
        eigvals = [0.51236726, 0.30071909, 0.26791207, 0.20898868,
                   0.19169895, 0.16054235,  0.15017696,  0.12245775,
                   0.0]
        proportion_explained = [0.2675738328, 0.157044696, 0.1399118638,
                                0.1091402725, 0.1001110485,
                                0.0838401162, 0.0784269939,
                                0.0639511764, 0.0]
        sample_ids = ['PC.636', 'PC.635', 'PC.356', 'PC.481', 'PC.354',
                      'PC.593', 'PC.355', 'PC.607', 'PC.634']
        axis_labels = ['PC%d' % i for i in range(1, 10)]

        expected_results = OrdinationResults(
            short_method_name='PCoA',
            long_method_name='Principal Coordinate Analysis',
            eigvals=pd.Series(eigvals, index=axis_labels),
            samples=pd.DataFrame(
                np.loadtxt(get_data_path('exp_PCoAEigenResults_site')),
                index=sample_ids, columns=axis_labels),
            proportion_explained=pd.Series(proportion_explained,
                                           index=axis_labels))

        dm = DistanceMatrix.read(get_data_path('PCoA_sample_data_3'))
        results = pcoa(dm)

        assert_ordination_results_equal(results, expected_results,
                                        ignore_directionality=True)

    def test_extensive(self):
        eigvals = [0.3984635, 0.36405689, 0.28804535, 0.27479983,
                   0.19165361, 0.0]
        proportion_explained = [0.2626621381, 0.2399817314,
                                0.1898758748, 0.1811445992,
                                0.1263356565, 0.0]
        sample_ids = [str(i) for i in range(6)]
        axis_labels = ['PC%d' % i for i in range(1, 7)]
        samples = [[-0.028597, 0.22903853, 0.07055272, 0.26163576,
                    0.28398669, 0.0],
                   [0.37494056, 0.22334055, -0.20892914, 0.05057395,
                    -0.18710366, 0.0],
                   [-0.33517593, -0.23855979, -0.3099887, 0.11521787,
                    -0.05021553, 0.0],
                   [0.25412394, -0.4123464, 0.23343642, 0.06403168,
                    -0.00482608, 0.0],
                   [-0.28256844, 0.18606911, 0.28875631, -0.06455635,
                    -0.21141632, 0.0],
                   [0.01727687, 0.012458, -0.07382761, -0.42690292,
                    0.1695749, 0.0]]

        expected_results = OrdinationResults(
            short_method_name='PCoA',
            long_method_name='Principal Coordinate Analysis',
            eigvals=pd.Series(eigvals, index=axis_labels),
            samples=pd.DataFrame(samples, index=sample_ids,
                                 columns=axis_labels),
            proportion_explained=pd.Series(proportion_explained,
                                           index=axis_labels))

        data = np.loadtxt(get_data_path('PCoA_sample_data_2'))
        # test passing a numpy.ndarray and a DistanceMatrix to pcoa
        # gives same results
        for dm in (data, DistanceMatrix(data)):
            results = pcoa(dm)
            assert_ordination_results_equal(results, expected_results,
                                            ignore_directionality=True)

    def test_book_example_dataset(self):
        # Adapted from PyCogent's `test_principal_coordinate_analysis`:
        #   "I took the example in the book (see intro info), and did
        #   the principal coordinates analysis, plotted the data and it
        #   looked right".
        eigvals = [0.73599103, 0.26260032, 0.14926222, 0.06990457,
                   0.02956972, 0.01931184, 0., 0., 0., 0., 0., 0., 0.,
                   0.]
        proportion_explained = [0.58105792, 0.20732046, 0.1178411,
                                0.05518899, 0.02334502, 0.01524651, 0.,
                                0., 0., 0., 0., 0., 0., 0.]
        sample_ids = [str(i) for i in range(14)]
        axis_labels = ['PC%d' % i for i in range(1, 15)]

        expected_results = OrdinationResults(
            short_method_name='PCoA',
            long_method_name='Principal Coordinate Analysis',
            eigvals=pd.Series(eigvals, index=axis_labels),
            samples=pd.DataFrame(
                np.loadtxt(get_data_path('exp_PCoAzeros_site')),
                index=sample_ids, columns=axis_labels),
            proportion_explained=pd.Series(proportion_explained,
                                           index=axis_labels))

        results = npt.assert_warns(RuntimeWarning, pcoa, self.dm)

        # Note the absolute value because column can have signs swapped
        results.samples = np.abs(results.samples)
        assert_ordination_results_equal(results, expected_results,
                                        ignore_directionality=True)

    def test_invalid_input(self):
        with npt.assert_raises(DissimilarityMatrixError):
            pcoa([[1, 2], [3, 4]])


if __name__ == "__main__":
    main()
