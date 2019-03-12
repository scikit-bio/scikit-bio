# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import numpy.testing as npt
import pandas as pd
from copy import deepcopy
from unittest import TestCase, main

from skbio import DistanceMatrix, OrdinationResults
from skbio.stats.distance import DissimilarityMatrixError
from skbio.stats.ordination import pcoa, pcoa_biplot
from skbio.util import get_data_path, assert_ordination_results_equal


class TestPCoA(TestCase):
    def setUp(self):
        # Sample data set from page 111 of W.J Krzanowski. Principles
        # of multivariate analysis, 2000, Oxford University Press.
        self.dm = DistanceMatrix(
            np.loadtxt(get_data_path('PCoA_sample_data')))

    def test_simple(self):
        eigvals = [0.51236726, 0.30071909, 0.26791207, 0.20898868,
                   0.19169895, 0.16054235, 0.15017696, 0.12245775,
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

    def test_fsvd_inplace(self):
        dm1 = DistanceMatrix.read(get_data_path('PCoA_sample_data_3'))
        dm2 = DistanceMatrix.read(get_data_path('PCoA_sample_data_3'))

        expected_results = pcoa(dm1, method="eigh", number_of_dimensions=3,
                                inplace=True)

        results = pcoa(dm2, method="fsvd", number_of_dimensions=3,
                       inplace=True)

        assert_ordination_results_equal(results, expected_results,
                                        ignore_directionality=True,
                                        ignore_method_names=True)

    def test_fsvd(self):
        dm1 = DistanceMatrix.read(get_data_path('PCoA_sample_data_3'))
        dm2 = DistanceMatrix.read(get_data_path('PCoA_sample_data_3'))
        dm3 = DistanceMatrix.read(get_data_path('PCoA_sample_data_3'))

        # Test eigh vs. fsvd pcoa and inplace parameter
        expected_results = pcoa(dm1, method="eigh", number_of_dimensions=3,
                                inplace=False)

        results = pcoa(dm2, method="fsvd", number_of_dimensions=3,
                       inplace=False)

        results_inplace = pcoa(dm2, method="fsvd", number_of_dimensions=3,
                               inplace=True)

        assert_ordination_results_equal(results, expected_results,
                                        ignore_directionality=True,
                                        ignore_method_names=True)

        assert_ordination_results_equal(results, results_inplace,
                                        ignore_directionality=True,
                                        ignore_method_names=True)

        # Test number_of_dimensions edge cases
        results2 = pcoa(dm3, method="fsvd", number_of_dimensions=0,
                        inplace=False)
        expected_results2 = pcoa(dm3, method="fsvd",
                                 number_of_dimensions=dm3.data.shape[0],
                                 inplace=False)

        assert_ordination_results_equal(results2, expected_results2,
                                        ignore_directionality=True,
                                        ignore_method_names=True)

        with self.assertRaises(ValueError):
            dim_too_large = dm1.data.shape[0] + 10
            pcoa(dm2, method="fsvd", number_of_dimensions=dim_too_large)

        with self.assertRaises(ValueError):
            pcoa(dm2, method="fsvd", number_of_dimensions=-1)

        with self.assertRaises(ValueError):
            dim_too_large = dm1.data.shape[0] + 10
            pcoa(dm2, method="eigh", number_of_dimensions=dim_too_large)

        with self.assertRaises(ValueError):
            pcoa(dm2, method="eigh", number_of_dimensions=-1)

        dm_big = DistanceMatrix.read(get_data_path('PCoA_sample_data_12dim'))
        with self.assertWarnsRegex(RuntimeWarning,
                                   r"no value for number_of_dimensions"):
            pcoa(dm_big, method="fsvd", number_of_dimensions=0)

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


class TestPCoABiplot(TestCase):
    def setUp(self):
        # Crawford dataset for unweighted UniFrac
        fp = get_data_path('PCoA_sample_data_3')
        self.ordination = pcoa(DistanceMatrix.read(fp))

        fp = get_data_path('PCoA_biplot_descriptors')
        self.descriptors = pd.read_table(fp, index_col='Taxon').T

    def test_pcoa_biplot_from_ape(self):
        """Test against a reference implementation from R's ape package

        The test data was generated with the R script below and using a
        modified version of pcoa.biplot that returns the U matrix.

        library(ape)
        # files can be found in the test data folder of the ordination module
        y = t(read.table('PCoA_biplot_descriptors', row.names = 1, header = 1))
        dm = read.table('PCoA_sample_data_3', row.names = 1, header = 1)

        h = pcoa(dm)

        # biplot.pcoa will only calculate the biplot for two axes at a time
        acc = NULL
        for (axes in c(1, 3, 5, 7)) {
            new = biplot.pcoa(h, y, plot.axes=c(axes, axes+1),
                              rn = rep('.', length(colnames(dm))) )

            if(is.null(acc)) {
                acc = new
            }
            else {
                b = acc
                acc <- cbind(acc, new)
            }
        }
        write.csv(acc, file='PCoA_biplot_projected_descriptors')
        """
        obs = pcoa_biplot(self.ordination, self.descriptors)

        # we'll build a dummy ordination results object based on the expected
        # the main thing we'll compare and modify is the features dataframe
        exp = deepcopy(obs)

        fp = get_data_path('PCoA_biplot_projected_descriptors')
        # R won't calculate the last dimension, so pad with zeros to make the
        # arrays comparable
        exp.features = pd.read_table(fp, sep=',', index_col=0)
        exp.features['Axis.9'] = np.zeros_like(exp.features['Axis.8'])

        # make the order comparable
        exp.features = exp.features.reindex(obs.features.index)

        assert_ordination_results_equal(obs, exp, ignore_directionality=True,
                                        ignore_axis_labels=True)

    def test_pcoa_biplot_subset_input(self):
        # create a 2D copy of the full ordination
        two_dims = deepcopy(self.ordination)
        two_dims.eigvals = two_dims.eigvals[:2]
        two_dims.samples = two_dims.samples.iloc[:, :2]
        two_dims.proportion_explained = two_dims.proportion_explained[:2]

        # only look at the features
        subset = pcoa_biplot(two_dims, self.descriptors).features
        full = pcoa_biplot(self.ordination, self.descriptors).features

        # the biplot should be identical regardless of the number of axes used
        pd.util.testing.assert_almost_equal(subset, full.iloc[:, :2])

    def test_mismatching_samples(self):
        new_index = self.descriptors.index.tolist()
        new_index[3] = 'Not.an.id'
        self.descriptors.index = pd.Index(new_index)

        with self.assertRaisesRegex(ValueError, r'The eigenvectors and the '
                                                'descriptors must describe '
                                                'the same '
                                                'samples.'):
            pcoa_biplot(self.ordination, self.descriptors)

    def test_not_a_pcoa(self):
        self.ordination.short_method_name = 'RDA'
        self.ordination.long_method_name = 'Redundancy Analysis'
        with self.assertRaisesRegex(ValueError, r'This biplot computation can'
                                                ' only be performed in a '
                                                'PCoA matrix.'):
            pcoa_biplot(self.ordination, self.descriptors)

    def test_from_seralized_results(self):
        # the current implementation of ordination results loses some
        # information, test that pcoa_biplot works fine regardless
        results = OrdinationResults.read(get_data_path('PCoA_skbio'))

        serialized = pcoa_biplot(results, self.descriptors)
        in_memory = pcoa_biplot(self.ordination, self.descriptors)

        assert_ordination_results_equal(serialized, in_memory,
                                        ignore_directionality=True,
                                        ignore_axis_labels=True,
                                        ignore_method_names=True)


if __name__ == "__main__":
    main()
