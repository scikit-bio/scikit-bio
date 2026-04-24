# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import numpy.testing as npt
import pandas as pd
from copy import deepcopy
from warnings import catch_warnings
from unittest import TestCase, main

from skbio import DistanceMatrix, OrdinationResults
from skbio.stats.distance import PairwiseMatrixError
from skbio.stats.ordination import pcoa, pcoa_biplot
from skbio.util import (get_data_path, assert_ordination_results_equal,
                        assert_data_frame_almost_equal, ArrayAPITestMixin,
                        array_backends)
from skbio.stats.ordination._principal_coordinate_analysis import (
    e_matrix,
    f_matrix,
    center_distance_matrix,
    _host_partial_eigh,
    _fsvd,
)


class TestPCoA(TestCase):
    def setUp(self):
        # Sample data set from page 111 of W.J Krzanowski. Principles
        # of multivariate analysis, 2000, Oxford University Press.
        self.dm = DistanceMatrix(np.loadtxt(get_data_path('PCoA_sample_data')))
        self.dm3 = DistanceMatrix.read(get_data_path('PCoA_sample_data_3'))
        self.dm_invalid = np.array([[1, 2], [3, 4], [5, 6]])

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

        results = pcoa(self.dm3)

        assert_ordination_results_equal(results, expected_results,
                                        ignore_directionality=True)

    def test_fsvd_inplace(self):
        expected_results = pcoa(
            self.dm3.copy(), method="eigh", dimensions=3, inplace=True)

        results = pcoa(self.dm3.copy(), method="fsvd", dimensions=3, inplace=True)

        assert_ordination_results_equal(results, expected_results,
                                        ignore_directionality=True,
                                        ignore_method_names=True)

    def test_fsvd(self):
        # Test eigh vs. fsvd pcoa and inplace parameter
        expected_results = pcoa(self.dm3, method="eigh", dimensions=3, inplace=False)

        results = pcoa(self.dm3, method="fsvd", dimensions=3, inplace=False)

        results_inplace = pcoa(self.dm3.copy(), method="fsvd", dimensions=3,
                               inplace=True)

        assert_ordination_results_equal(results, expected_results,
                                        ignore_directionality=True,
                                        ignore_method_names=True)

        assert_ordination_results_equal(results, results_inplace,
                                        ignore_directionality=True,
                                        ignore_method_names=True)

        # Test dimensions edge cases
        results2 = pcoa(self.dm3, method="fsvd", dimensions=0, inplace=False)
        expected_results2 = pcoa(self.dm3, method="fsvd", dimensions=self.dm3.shape[0],
                                 inplace=False)

        assert_ordination_results_equal(results2, expected_results2,
                                        ignore_directionality=True,
                                        ignore_method_names=True)

        with self.assertRaises(ValueError):
            dim_too_large = self.dm3.shape[0] + 10
            pcoa(self.dm3, method="fsvd", dimensions=dim_too_large)

        with self.assertRaises(ValueError):
            pcoa(self.dm3, method="fsvd", dimensions=-1)

        with self.assertRaises(ValueError):
            dim_too_large = self.dm3.shape[0] + 10
            pcoa(self.dm3, method="eigh", dimensions=dim_too_large)

        with self.assertRaises(ValueError):
            pcoa(self.dm3, method="eigh", dimensions=-1)

        dm_big = DistanceMatrix.read(get_data_path('PCoA_sample_data_12dim'))
        with self.assertWarnsRegex(RuntimeWarning,
                                   r"no value for dimensions"):
            pcoa(dm_big, method="fsvd", dimensions=0)

    def test_permutted(self):
        # this should not throw
        pcoa(self.dm3, method="fsvd", dimensions=3, inplace=False)

        # some operations, like permute, will change memory structure
        # we want to test that this does not break pcoa
        permutted = self.dm3.permute()
        # we just want to assure it does not throw
        pcoa(permutted, method="fsvd", dimensions=3, inplace=False)

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

        with self.assertWarns(RuntimeWarning):
            results = pcoa(self.dm)

        # Note the absolute value because column can have signs swapped
        results.samples = np.abs(results.samples)
        assert_ordination_results_equal(results, expected_results,
                                        ignore_directionality=True)

    def test_invalid_input(self):
        with self.assertRaises(PairwiseMatrixError):
            pcoa([[1, 2], [3, 4]])

    def test_warn_neg_eigval(self):
        """Test warnings of negative eigenvalues."""
        # In this example, negative-most: -0.109, positive-most: 0.736, ratio: 0.148,
        # which is above the threshold, therefore a warning is raised by default.
        # Note: Explicitly pass the dimensions, to avoid spurious "large DM" warning
        with self.assertWarns(RuntimeWarning):
            results = pcoa(self.dm, dimensions=self.dm.data.shape[0])

        # warn regardless of magnitude
        with self.assertWarns(RuntimeWarning):
            results = pcoa(self.dm, dimensions=self.dm.data.shape[0], warn_neg_eigval=True)

        # disable warning
        with catch_warnings(record=True) as obs:
            results = pcoa(self.dm, dimensions=self.dm.data.shape[0], warn_neg_eigval=False)
        self.assertEqual(obs, [])

        # larger (more stringent) threshold
        with catch_warnings(record=True) as obs:
            results = pcoa(self.dm, dimensions=self.dm.data.shape[0], warn_neg_eigval=0.2)
        self.assertEqual(obs, [])

        # In this example, all eigenvalues are zero or positive, therefore no warning
        # will be raised.
        with catch_warnings(record=True) as obs:
            results = pcoa(self.dm3)
        self.assertEqual(obs, [])

        with catch_warnings(record=True) as obs:
            results = pcoa(self.dm3, warn_neg_eigval=True)
        self.assertEqual(obs, [])

        # invalid parameter
        msg = ("warn_neg_eigval must be Boolean or a floating-point number between 0 "
               "and 1.")
        with self.assertRaisesRegex(ValueError, msg):
            pcoa(self.dm3, warn_neg_eigval=5.0)
        with self.assertRaisesRegex(ValueError, msg):
            pcoa(self.dm3, warn_neg_eigval=-2.5)

    def test_integer_dimensions(self):
        """Test with an integer dimensions."""
        results = pcoa(self.dm3, dimensions=3)
        self.assertEqual(results.samples.shape[1], 3)

    def test_large_float(self):
        """Test with a float dimensions > 1."""
        msg = "A floating-point number greater than 1 cannot be"
        with self.assertRaisesRegex(ValueError, msg):
            pcoa(self.dm3, dimensions=2.5)

    def test_eigh_method_with_float(self):
        """Test with a float dimensions with eigh method to retain 80%
        variance."""
        results = pcoa(self.dm3, method="eigh", dimensions=0.8,
                       inplace=False, seed=None)
        cumulative_variance = np.cumsum(results.proportion_explained.values)
        self.assertGreaterEqual(cumulative_variance[-1], 0.8)

    def test_edge_case_for_all_variance(self):
        """Test with a dimensions close to 1 to retain nearly all variance.
        """
        results = pcoa(self.dm3, dimensions=0.9999)
        cumulative_variance = np.cumsum(results.proportion_explained.values)
        self.assertGreaterEqual(cumulative_variance[-1], 0.9999)

    def test_fsvd_method_with_float(self):
        """Test FSVD with float dimensions for variance threshold."""
        results = pcoa(self.dm3, method="fsvd", dimensions=0.7,
                       inplace=False, seed=None)
        cumulative_variance = np.cumsum(results.proportion_explained.values)
        self.assertGreaterEqual(cumulative_variance[-1], 0.7)

    def test_invalid_method(self):
        """Test that correct error is raised when method is invalid."""
        with self.assertRaisesRegex(
            ValueError,
            "PCoA eigendecomposition method asdf not supported."
            ):
            results = pcoa(self.dm3, method="asdf")

    def test_fsvd_non_square_input(self):
        with self.assertRaisesRegex(ValueError, "FSVD expects square distance matrix"):
            results = _fsvd(self.dm_invalid)

    def test_fsvd_invalid_dimensions(self):
        with self.assertRaisesRegex(
            ValueError,
            ("Invalid operation: cannot reduce distance matrix "
            "to negative dimensions using PCoA. Did you intend " 
            'to specify the default value "0", which sets ' 
            "the dimensions equal to the " 
            "dimensionality of the given distance matrix?")
            ):
            results = _fsvd(self.dm_invalid[1:4], dimensions=-1)


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
        assert_data_frame_almost_equal(subset, full.iloc[:, :2])

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


# A small symmetric distance matrix used by most helpers.
_DM = [
    [0.0, 0.7, 0.6, 0.5],
    [0.7, 0.0, 0.4, 0.3],
    [0.6, 0.4, 0.0, 0.2],
    [0.5, 0.3, 0.2, 0.0],
]

class TestEMatrixBackends(TestCase, ArrayAPITestMixin):
    @array_backends("numpy", "jax", "torch", "cupy")
    def test_e_matrix_backends(self, xp, device):
        data = self.make_array(xp, device, _DM)
        result = e_matrix(data)
 
        # Reference computed in numpy.
        data_np = np.asarray(_DM, dtype=np.float64)
        expected = -(data_np * data_np) / 2
 
        self.assert_type_preserved(result, xp, device)
        self.assert_close(result, expected)
 
 
class TestFMatrixBackends(TestCase, ArrayAPITestMixin):
    @array_backends("numpy", "jax", "torch", "cupy")
    def test_f_matrix_backends(self, xp, device):
        data_np = np.asarray(_DM, dtype=np.float64)
        e_np = -(data_np * data_np) / 2
 
        e = self.make_array(xp, device, e_np.tolist())
        result = f_matrix(e)
 
        # Reference: manual double-centering in numpy.
        row = e_np.mean(axis=1, keepdims=True)
        col = e_np.mean(axis=0, keepdims=True)
        expected = e_np - row - col + e_np.mean()
 
        self.assert_type_preserved(result, xp, device)
        self.assert_close(result, expected)
 
    @array_backends("numpy", "jax", "torch", "cupy")
    def test_f_matrix_is_symmetric(self, xp, device):
        data = self.make_array(xp, device, _DM)
        result = f_matrix(e_matrix(data))
        # Double-centered matrix of a symmetric input should be symmetric.
        # Using .T instead of xp.matrix_transpose (not available in torch).
        self.assert_close(result, result.T)
 
 
class TestCenterDistanceMatrixBackends(TestCase, ArrayAPITestMixin):
    @array_backends("numpy", "jax", "torch", "cupy")
    def test_center_distance_matrix_backends(self, xp, device):
        data = self.make_array(xp, device, _DM)
        result = center_distance_matrix(data)
 
        # Reference: run the numpy path.
        data_np = np.asarray(_DM, dtype=np.float64)
        expected = center_distance_matrix(data_np)
 
        self.assert_type_preserved(result, xp, device)
        self.assert_close(result, expected)
 
    @array_backends("numpy", "jax", "torch", "cupy")
    def test_center_distance_matrix_inplace_ignored_for_non_numpy(self, xp, device):
        # `inplace=True` is documented to be ignored for non-numpy arrays;
        # the call should still succeed and produce the centered matrix.
        data = self.make_array(xp, device, _DM)
        result = center_distance_matrix(data, inplace=True)
 
        data_np = np.asarray(_DM, dtype=np.float64)
        expected = center_distance_matrix(data_np)
 
        self.assert_close(result, expected)
 
 
class TestHostPartialEighBackends(TestCase, ArrayAPITestMixin):
    @array_backends("numpy", "jax", "torch", "cupy")
    def test_host_partial_eigh_returns_to_device(self, xp, device):
        # Build a symmetric positive-ish matrix on the backend.
        centered_np = center_distance_matrix(np.asarray(_DM, dtype=np.float64))
        centered = self.make_array(xp, device, centered_np.tolist())
 
        n = centered_np.shape[0]
        subidx = [n - 2, n - 1]  # top two eigenpairs
 
        eigvals, eigvecs = _host_partial_eigh(centered, subidx)
 
        # Results should live on the original backend/device.
        self.assert_type_preserved(eigvals, xp, device)
        self.assert_type_preserved(eigvecs, xp, device)
 
        # And match a direct scipy call on the numpy reference.
        from scipy.linalg import eigh as scipy_eigh
 
        exp_vals, exp_vecs = scipy_eigh(centered_np, subset_by_index=subidx)
        self.assert_close(eigvals, exp_vals)
        # Eigenvectors can differ by sign; compare |vecs|.
        self.assert_close(xp.abs(eigvecs), np.abs(exp_vecs))
 
 
class TestFSVDBackends(TestCase, ArrayAPITestMixin):
    @array_backends("numpy", "jax", "torch", "cupy")
    def test_fsvd_shape_and_backend(self, xp, device):
        centered_np = center_distance_matrix(np.asarray(_DM, dtype=np.float64))
        centered = self.make_array(xp, device, centered_np.tolist())
 
        dims = 2
        eigvals, eigvecs = _fsvd(centered, dimensions=dims, seed=42)
 
        self.assert_type_preserved(eigvals, xp, device)
        self.assert_type_preserved(eigvecs, xp, device)
        self.assertEqual(tuple(eigvals.shape), (dims,))
        self.assertEqual(tuple(eigvecs.shape), (centered_np.shape[0], dims))
 
    @array_backends("numpy", "jax", "torch", "cupy")
    def test_fsvd_matches_numpy(self, xp, device):
        # Same seed on both paths → same (up to sign) subspace.
        centered_np = center_distance_matrix(np.asarray(_DM, dtype=np.float64))
        dims = 2
 
        ref_vals, ref_vecs = _fsvd(centered_np, dimensions=dims, seed=0)
 
        centered = self.make_array(xp, device, centered_np.tolist())
        res_vals, res_vecs = _fsvd(centered, dimensions=dims, seed=0)
 
        # Eigenvalues are unique; eigenvectors can flip sign per column.
        self.assert_close(res_vals, ref_vals, rtol=1e-5, atol=1e-5)
        self.assert_close(xp.abs(res_vecs), np.abs(ref_vecs), rtol=1e-5, atol=1e-5)


if __name__ == "__main__":
    main()
