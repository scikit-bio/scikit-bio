# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from six import binary_type, text_type, StringIO

import warnings
import unittest

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import numpy.testing as npt
import pandas as pd
from IPython.core.display import Image, SVG
from nose.tools import assert_is_instance, assert_raises_regexp, assert_true
from scipy.spatial.distance import pdist

from skbio import DistanceMatrix
from skbio.stats.ordination import (
    CA, RDA, CCA, PCoA, OrdinationResults, corr, mean_and_std,
    assert_ordination_results_equal)
from skbio.util import get_data_path


def normalize_signs(arr1, arr2):
    """Change column signs so that "column" and "-column" compare equal.

    This is needed because results of eigenproblmes can have signs
    flipped, but they're still right.

    Notes
    =====

    This function tries hard to make sure that, if you find "column"
    and "-column" almost equal, calling a function like np.allclose to
    compare them after calling `normalize_signs` succeeds.

    To do so, it distinguishes two cases for every column:

    - It can be all almost equal to 0 (this includes a column of
      zeros).
    - Otherwise, it has a value that isn't close to 0.

    In the first case, no sign needs to be flipped. I.e., for
    |epsilon| small, np.allclose(-epsilon, 0) is true if and only if
    np.allclose(epsilon, 0) is.

    In the second case, the function finds the number in the column
    whose absolute value is largest. Then, it compares its sign with
    the number found in the same index, but in the other array, and
    flips the sign of the column as needed.
    """
    # Let's convert everyting to floating point numbers (it's
    # reasonable to assume that eigenvectors will already be floating
    # point numbers). This is necessary because np.array(1) /
    # np.array(0) != np.array(1.) / np.array(0.)
    arr1 = np.asarray(arr1, dtype=np.float64)
    arr2 = np.asarray(arr2, dtype=np.float64)

    if arr1.shape != arr2.shape:
        raise ValueError(
            "Arrays must have the same shape ({0} vs {1}).".format(arr1.shape,
                                                                   arr2.shape)
            )

    # To avoid issues around zero, we'll compare signs of the values
    # with highest absolute value
    max_idx = np.abs(arr1).argmax(axis=0)
    max_arr1 = arr1[max_idx, range(arr1.shape[1])]
    max_arr2 = arr2[max_idx, range(arr2.shape[1])]

    sign_arr1 = np.sign(max_arr1)
    sign_arr2 = np.sign(max_arr2)

    # Store current warnings, and ignore division by zero (like 1. /
    # 0.) and invalid operations (like 0. / 0.)
    wrn = np.seterr(invalid='ignore', divide='ignore')
    differences = sign_arr1 / sign_arr2
    # The values in `differences` can be:
    #    1 -> equal signs
    #   -1 -> diff signs
    #   Or nan (0/0), inf (nonzero/0), 0 (0/nonzero)
    np.seterr(**wrn)

    # Now let's deal with cases where `differences != \pm 1`
    special_cases = (~np.isfinite(differences)) | (differences == 0)
    # In any of these cases, the sign of the column doesn't matter, so
    # let's just keep it
    differences[special_cases] = 1

    return arr1 * differences, arr2


def chi_square_distance(data_table, between_rows=True):
    """Computes the chi-square distance between two rows or columns of input.

    It is a measure that has no upper limit, and it excludes double-zeros.

    Parameters
    ----------
    data_table : 2D array_like
        An array_like object of shape (n, p). The input must be a
        frequency table (so that the sum of all cells equals 1, and
        all values are non-negative).
    between_rows : bool (defaults to True)
        Indicates whether distance is computed between rows (default)
        or columns.

    Returns
    -------
    Y : ndarray
        Returns a condensed distance matrix. For each i and j (where
        i<j<n), the chi square distance between u=X[i] and v=X[j] is
        computed and stored in `Y[(n choose 2) - (n - i choose 2) + (j
        - i - 1)]`.

    See Also
    --------
    scipy.spatial.distance.squareform

    References
    ----------
    This coefficient appears in Legendre and Legendre (1998) as
    formula 7.54 (as D_{16}). Another source is
    http://www.springerreference.com/docs/html/chapterdbid/60817.html
    """
    data_table = np.asarray(data_table, dtype=np.float64)
    if not np.allclose(data_table.sum(), 1):
        raise ValueError("Input is not a frequency table: if it is an"
                         " abundance table you could scale it as"
                         " `data_table / data_table.sum()`.")
    if np.any(data_table < 0):
        raise ValueError("A frequency table can't have negative values.")

    # The distances are always computed between the rows of F
    F = data_table if between_rows else data_table.T

    row_sums = F.sum(axis=1, keepdims=True)
    column_sums = F.sum(axis=0)
    scaled_F = F / (row_sums * np.sqrt(column_sums))

    return pdist(scaled_F, 'euclidean')


class TestNormalizeSigns(object):
    def test_shapes_and_nonarray_input(self):
        with npt.assert_raises(ValueError):
            normalize_signs([[1, 2], [3, 5]], [[1, 2]])

    def test_works_when_different(self):
        """Taking abs value of everything would lead to false
        positives."""
        a = np.array([[1, -1],
                      [2, 2]])
        b = np.array([[-1, -1],
                      [2, 2]])
        with npt.assert_raises(AssertionError):
            npt.assert_equal(*normalize_signs(a, b))

    def test_easy_different(self):
        a = np.array([[1, 2],
                      [3, -1]])
        b = np.array([[-1, 2],
                      [-3, -1]])
        npt.assert_equal(*normalize_signs(a, b))

    def test_easy_already_equal(self):
        a = np.array([[1, -2],
                      [3, 1]])
        b = a.copy()
        npt.assert_equal(*normalize_signs(a, b))

    def test_zeros(self):
        a = np.array([[0, 3],
                      [0, -1]])
        b = np.array([[0, -3],
                      [0, 1]])
        npt.assert_equal(*normalize_signs(a, b))

    def test_hard(self):
        a = np.array([[0, 1],
                      [1, 2]])
        b = np.array([[0, 1],
                      [-1, 2]])
        npt.assert_equal(*normalize_signs(a, b))

    def test_harder(self):
        """We don't want a value that might be negative due to
        floating point inaccuracies to make a call to allclose in the
        result to be off."""
        a = np.array([[-1e-15, 1],
                      [5, 2]])
        b = np.array([[1e-15, 1],
                      [5, 2]])
        # Clearly a and b would refer to the same "column
        # eigenvectors" but a slopppy implementation of
        # normalize_signs could change the sign of column 0 and make a
        # comparison fail
        npt.assert_almost_equal(*normalize_signs(a, b))

    def test_column_zeros(self):
        a = np.array([[0, 1],
                      [0, 2]])
        b = np.array([[0, -1],
                      [0, -2]])
        npt.assert_equal(*normalize_signs(a, b))

    def test_column_almost_zero(self):
        a = np.array([[1e-15, 3],
                      [-2e-14, -6]])
        b = np.array([[0, 3],
                      [-1e-15, -6]])
        npt.assert_almost_equal(*normalize_signs(a, b))


class TestChiSquareDistance(object):
    def test_errors(self):
        a = np.array([[-0.5, 0],
                      [1, 0.5]])
        with npt.assert_raises(ValueError):
            chi_square_distance(a)
        b = np.array([[0.5, 0],
                      [0.5, 0.1]])
        with npt.assert_raises(ValueError):
            chi_square_distance(b)

    def test_results(self):
        """Some random numbers."""
        a = np.array([[0.02808988764,  0.056179775281,  0.084269662921,
                       0.140449438202],
                      [0.01404494382,  0.196629213483,  0.109550561798,
                       0.033707865169],
                      [0.02808988764,  0.112359550562,  0.056179775281,
                       0.140449438202]])
        dist = chi_square_distance(a)
        expected = [0.91413919964333856,
                    0.33651110106124049,
                    0.75656884966269089]
        npt.assert_almost_equal(dist, expected)

    def test_results2(self):
        """A tiny example from Legendre & Legendre 1998, p. 285."""
        a = np.array([[0, 1, 1],
                      [1, 0, 0],
                      [0, 4, 4]])
        dist = chi_square_distance(a / a.sum())
        # Note L&L used a terrible calculator because they got a wrong
        # number (says it's 3.477) :(
        expected = [3.4785054261852175, 0, 3.4785054261852175]
        npt.assert_almost_equal(dist, expected)


class TestUtils(object):
    def setup(self):
        self.x = np.array([[1, 2, 3], [4, 5, 6]])
        self.y = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])

    def test_mean_and_std_no_mean_no_std(self):
        with npt.assert_raises(ValueError):
            mean_and_std(self.x, with_mean=False, with_std=False)

    def test_corr_shape_mismatch(self):
        with npt.assert_raises(ValueError):
            corr(self.x, self.y)

    def test_assert_ordination_results_equal(self):
        minimal1 = OrdinationResults([1, 2])

        # a minimal set of results should be equal to itself
        assert_ordination_results_equal(minimal1, minimal1)

        # type mismatch
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal(minimal1, 'foo')

        # numeric values should be checked that they're almost equal
        almost_minimal1 = OrdinationResults([1.0000001, 1.9999999])
        assert_ordination_results_equal(minimal1, almost_minimal1)

        # species_ids missing in one, present in the other
        almost_minimal1.species_ids = ['abc', 'def']
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal(minimal1, almost_minimal1)
        almost_minimal1.species_ids = None

        # site_ids missing in one, present in the other
        almost_minimal1.site_ids = ['abc', 'def']
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal(minimal1, almost_minimal1)
        almost_minimal1.site_ids = None

        # test each of the optional numeric attributes
        for attr in ('species', 'site', 'biplot', 'site_constraints',
                     'proportion_explained'):
            # missing optional numeric attribute in one, present in the other
            setattr(almost_minimal1, attr, [[1, 2], [3, 4]])
            with npt.assert_raises(AssertionError):
                assert_ordination_results_equal(minimal1, almost_minimal1)
            setattr(almost_minimal1, attr, None)

            # optional numeric attributes present in both, but not almost equal
            setattr(minimal1, attr, [[1, 2], [3, 4]])
            setattr(almost_minimal1, attr, [[1, 2], [3.00002, 4]])
            with npt.assert_raises(AssertionError):
                assert_ordination_results_equal(minimal1, almost_minimal1)
            setattr(minimal1, attr, None)
            setattr(almost_minimal1, attr, None)

            # optional numeric attributes present in both, and almost equal
            setattr(minimal1, attr, [[1, 2], [3, 4]])
            setattr(almost_minimal1, attr, [[1, 2], [3.00000002, 4]])
            assert_ordination_results_equal(minimal1, almost_minimal1)
            setattr(minimal1, attr, None)
            setattr(almost_minimal1, attr, None)


class TestCAResults(object):
    def setup(self):
        """Data from table 9.11 in Legendre & Legendre 1998."""
        self.X = np.loadtxt(get_data_path('L&L_CA_data'))
        self.ordination = CA(self.X, ['Site1', 'Site2', 'Site3'],
                             ['Species1', 'Species2', 'Species3'])

    def test_scaling2(self):
        scores = self.ordination.scores(scaling=2)
        # p. 460 L&L 1998
        F_hat = np.array([[0.40887, -0.06955],
                          [-0.11539,  0.29977],
                          [-0.30997, -0.18739]])
        npt.assert_almost_equal(*normalize_signs(F_hat, scores.species),
                                decimal=5)
        V_hat = np.array([[-0.84896, -0.88276],
                          [-0.22046,  1.34482],
                          [1.66697, -0.47032]])
        npt.assert_almost_equal(*normalize_signs(V_hat, scores.site),
                                decimal=5)

    def test_scaling1(self):
        scores = self.ordination.scores(scaling=1)
        # p. 458
        V = np.array([[1.31871, -0.34374],
                      [-0.37215,  1.48150],
                      [-0.99972, -0.92612]])
        npt.assert_almost_equal(*normalize_signs(V, scores.species), decimal=5)
        F = np.array([[-0.26322, -0.17862],
                      [-0.06835,  0.27211],
                      [0.51685, -0.09517]])
        npt.assert_almost_equal(*normalize_signs(F, scores.site), decimal=5)

    def test_maintain_chi_square_distance_scaling1(self):
        """In scaling 1, chi^2 distance among rows (sites) is equal to
        euclidean distance between them in transformed space."""
        frequencies = self.X / self.X.sum()
        chi2_distances = chi_square_distance(frequencies)
        transformed_sites = self.ordination.scores(1).site
        euclidean_distances = pdist(transformed_sites, 'euclidean')
        npt.assert_almost_equal(chi2_distances, euclidean_distances)

    def test_maintain_chi_square_distance_scaling2(self):
        """In scaling 2, chi^2 distance among columns (species) is
        equal to euclidean distance between them in transformed space."""
        frequencies = self.X / self.X.sum()
        chi2_distances = chi_square_distance(frequencies, between_rows=False)
        transformed_species = self.ordination.scores(2).species
        euclidean_distances = pdist(transformed_species, 'euclidean')
        npt.assert_almost_equal(chi2_distances, euclidean_distances)


class TestCAErrors(object):
    def test_negative(self):
        X = np.array([[1, 2], [-0.1, -2]])
        with npt.assert_raises(ValueError):
            CA(X, None, None)


class TestRDAErrors(object):
    def test_shape(self):
        for n, p, n_, m in [(3, 4, 2, 1), (3, 4, 3, 10)]:
            Y = np.random.randn(n, p)
            X = np.random.randn(n_, m)
            yield npt.assert_raises, ValueError, RDA, Y, X, None, None


class TestRDAResults(object):
    # STATUS: L&L only shows results with scaling 1, and they agree
    # with vegan's (module multiplying by a constant). I can also
    # compute scaling 2, agreeing with vegan, but there are no written
    # results in L&L.
    def setup(self):
        """Data from table 11.3 in Legendre & Legendre 1998."""
        Y = np.loadtxt(get_data_path('example2_Y'))
        X = np.loadtxt(get_data_path('example2_X'))
        self.ordination = RDA(Y, X,
                              ['Site0', 'Site1', 'Site2', 'Site3', 'Site4',
                               'Site5', 'Site6', 'Site7', 'Site8', 'Site9'],
                              ['Species0', 'Species1', 'Species2', 'Species3',
                               'Species4', 'Species5'])

    def test_scaling1(self):
        scores = self.ordination.scores(1)

        # Load data as computed with vegan 2.0-8
        vegan_species = np.loadtxt(get_data_path(
            'example2_species_scaling1_from_vegan'))
        npt.assert_almost_equal(scores.species, vegan_species, decimal=6)

        vegan_site = np.loadtxt(get_data_path(
            'example2_site_scaling1_from_vegan'))
        npt.assert_almost_equal(scores.site, vegan_site, decimal=6)

    def test_scaling2(self):
        scores = self.ordination.scores(2)

        # Load data as computed with vegan 2.0-8
        vegan_species = np.loadtxt(get_data_path(
            'example2_species_scaling2_from_vegan'))
        npt.assert_almost_equal(scores.species, vegan_species, decimal=6)

        vegan_site = np.loadtxt(get_data_path(
            'example2_site_scaling2_from_vegan'))
        npt.assert_almost_equal(scores.site, vegan_site, decimal=6)


class TestCCAErrors(object):
    def setup(self):
        """Data from table 11.3 in Legendre & Legendre 1998."""
        self.Y = np.loadtxt(get_data_path('example3_Y'))
        self.X = np.loadtxt(get_data_path('example3_X'))

    def test_shape(self):
        X, Y = self.X, self.Y
        with npt.assert_raises(ValueError):
            CCA(Y, X[:-1], None, None)

    def test_Y_values(self):
        X, Y = self.X, self.Y
        Y[0, 0] = -1
        with npt.assert_raises(ValueError):
            CCA(Y, X, None, None)
        Y[0] = 0
        with npt.assert_raises(ValueError):
            CCA(Y, X, None, None)


class TestCCAResults(object):
    def setup(self):
        """Data from table 11.3 in Legendre & Legendre 1998
        (p. 590). Loaded results as computed with vegan 2.0-8 and
        compared with table 11.5 if also there."""
        Y = np.loadtxt(get_data_path('example3_Y'))
        X = np.loadtxt(get_data_path('example3_X'))
        self.ordination = CCA(Y, X[:, :-1],
                              ['Site0', 'Site1', 'Site2', 'Site3', 'Site4',
                               'Site5', 'Site6', 'Site7', 'Site8', 'Site9'],
                              ['Species0', 'Species1', 'Species2', 'Species3',
                               'Species4', 'Species5', 'Species6', 'Species7',
                               'Species8'])

    def test_scaling1_species(self):
        scores = self.ordination.scores(1)

        vegan_species = np.loadtxt(get_data_path(
            'example3_species_scaling1_from_vegan'))
        npt.assert_almost_equal(scores.species, vegan_species, decimal=6)

    def test_scaling1_site(self):
        scores = self.ordination.scores(1)

        vegan_site = np.loadtxt(get_data_path(
            'example3_site_scaling1_from_vegan'))
        npt.assert_almost_equal(scores.site, vegan_site, decimal=4)

    def test_scaling2_species(self):
        scores = self.ordination.scores(2)

        vegan_species = np.loadtxt(get_data_path(
            'example3_species_scaling2_from_vegan'))
        npt.assert_almost_equal(scores.species, vegan_species, decimal=5)

    def test_scaling2_site(self):
        scores = self.ordination.scores(2)

        vegan_site = np.loadtxt(get_data_path(
            'example3_site_scaling2_from_vegan'))
        npt.assert_almost_equal(scores.site, vegan_site, decimal=4)


class TestPCoAResults(object):
    def setup(self):
        """Sample data set from page 111 of W.J Krzanowski. Principles
        of multivariate analysis, 2000, Oxford University Press."""
        matrix = np.loadtxt(get_data_path('PCoA_sample_data'))
        dist_matrix = DistanceMatrix(matrix, map(str, range(matrix.shape[0])))
        self.dist_matrix = dist_matrix

    def test_negative_eigenvalue_warning(self):
        """This data has some small negative eigenvalues."""
        npt.assert_warns(RuntimeWarning, PCoA, self.dist_matrix)

    def test_values(self):
        """Adapted from cogent's `test_principal_coordinate_analysis`:
        "I took the example in the book (see intro info), and did the
        principal coordinates analysis, plotted the data and it looked
        right"."""
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', category=RuntimeWarning)
            ordination = PCoA(self.dist_matrix)
        scores = ordination.scores()

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


class TestPCoAResultsExtensive(object):
    def setup(self):
        matrix = np.loadtxt(get_data_path('PCoA_sample_data_2'))
        self.ids = [str(i) for i in range(matrix.shape[0])]
        dist_matrix = DistanceMatrix(matrix, self.ids)
        self.ordination = PCoA(dist_matrix)

    def test_values(self):
        results = self.ordination.scores()

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


class TestPCoAEigenResults(object):
    def setup(self):
        dist_matrix = DistanceMatrix.read(get_data_path('PCoA_sample_data_3'))
        self.ordination = PCoA(dist_matrix)

        self.ids = ['PC.636', 'PC.635', 'PC.356', 'PC.481', 'PC.354', 'PC.593',
                    'PC.355', 'PC.607', 'PC.634']

    def test_values(self):
        results = self.ordination.scores()

        npt.assert_almost_equal(len(results.eigvals), len(results.site[0]))

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


class TestPCoAPrivateMethods(object):
    def setup(self):
        self.matrix = np.arange(1, 7).reshape(2, 3)
        self.matrix2 = np.arange(1, 10).reshape(3, 3)

    def test_E_matrix(self):
        E = PCoA._E_matrix(self.matrix)
        expected_E = np.array([[-0.5,  -2.,  -4.5],
                               [-8., -12.5, -18.]])
        npt.assert_almost_equal(E, expected_E)

    def test_F_matrix(self):
        F = PCoA._F_matrix(self.matrix2)
        expected_F = np.zeros((3, 3))
        # Note that `test_make_F_matrix` in cogent is wrong
        npt.assert_almost_equal(F, expected_F)


class TestPCoAErrors(object):
    def test_input(self):
        with npt.assert_raises(TypeError):
            PCoA([[1, 2], [3, 4]])


class TestOrdinationResults(unittest.TestCase):
    def setUp(self):
        # Define in-memory CA results to serialize and deserialize.
        eigvals = np.array([0.0961330159181, 0.0409418140138])
        species = np.array([[0.408869425742, 0.0695518116298],
                            [-0.1153860437, -0.299767683538],
                            [-0.309967102571, 0.187391917117]])
        site = np.array([[-0.848956053187, 0.882764759014],
                         [-0.220458650578, -1.34482000302],
                         [1.66697179591, 0.470324389808]])
        biplot = None
        site_constraints = None
        prop_explained = None
        species_ids = ['Species1', 'Species2', 'Species3']
        site_ids = ['Site1', 'Site2', 'Site3']

        self.ordination_results = OrdinationResults(
            eigvals=eigvals, species=species, site=site, biplot=biplot,
            site_constraints=site_constraints,
            proportion_explained=prop_explained, species_ids=species_ids,
            site_ids=site_ids)

        # DataFrame for testing plot method. Has a categorical column with a
        # mix of numbers and strings. Has a numeric column with a mix of ints,
        # floats, and strings that can be converted to floats. Has a numeric
        # column with missing data (np.nan).
        self.df = pd.DataFrame([['foo', '42', 10],
                                [22, 0, 8],
                                [22, -4.2, np.nan],
                                ['foo', '42.19', 11]],
                               index=['A', 'B', 'C', 'D'],
                               columns=['categorical', 'numeric', 'nancolumn'])

        # Minimal ordination results for easier testing of plotting method.
        # Paired with df above.
        eigvals = np.array([0.50, 0.25, 0.25])
        site = np.array([[0.1, 0.2, 0.3],
                         [0.2, 0.3, 0.4],
                         [0.3, 0.4, 0.5],
                         [0.4, 0.5, 0.6]])
        self.min_ord_results = OrdinationResults(eigvals=eigvals, site=site,
                                                 site_ids=['A', 'B', 'C', 'D'])

    def test_deprecated_io(self):
        fh = StringIO()
        npt.assert_warns(UserWarning, self.ordination_results.to_file, fh)
        fh.seek(0)
        deserialized = npt.assert_warns(UserWarning,
                                        OrdinationResults.from_file, fh)
        assert_ordination_results_equal(deserialized, self.ordination_results)
        self.assertTrue(type(deserialized) == OrdinationResults)

    def test_str(self):
        exp = ("Ordination results:\n"
               "\tEigvals: 2\n"
               "\tProportion explained: N/A\n"
               "\tSpecies: 3x2\n"
               "\tSite: 3x2\n"
               "\tBiplot: N/A\n"
               "\tSite constraints: N/A\n"
               "\tSpecies IDs: 'Species1', 'Species2', 'Species3'\n"
               "\tSite IDs: 'Site1', 'Site2', 'Site3'")
        obs = str(self.ordination_results)
        self.assertEqual(obs, exp)

        # all optional attributes missing
        exp = ("Ordination results:\n"
               "\tEigvals: 1\n"
               "\tProportion explained: N/A\n"
               "\tSpecies: N/A\n"
               "\tSite: N/A\n"
               "\tBiplot: N/A\n"
               "\tSite constraints: N/A\n"
               "\tSpecies IDs: N/A\n"
               "\tSite IDs: N/A")
        obs = str(OrdinationResults(np.array([4.2])))
        self.assertEqual(obs, exp)

    def check_basic_figure_sanity(self, fig, exp_num_subplots, exp_title,
                                  exp_legend_exists, exp_xlabel, exp_ylabel,
                                  exp_zlabel):
        # check type
        assert_is_instance(fig, mpl.figure.Figure)

        # check number of subplots
        axes = fig.get_axes()
        npt.assert_equal(len(axes), exp_num_subplots)

        # check title
        ax = axes[0]
        npt.assert_equal(ax.get_title(), exp_title)

        # shouldn't have tick labels
        for tick_label in (ax.get_xticklabels() + ax.get_yticklabels() +
                           ax.get_zticklabels()):
            npt.assert_equal(tick_label.get_text(), '')

        # check if legend is present
        legend = ax.get_legend()
        if exp_legend_exists:
            assert_true(legend is not None)
        else:
            assert_true(legend is None)

        # check axis labels
        npt.assert_equal(ax.get_xlabel(), exp_xlabel)
        npt.assert_equal(ax.get_ylabel(), exp_ylabel)
        npt.assert_equal(ax.get_zlabel(), exp_zlabel)

    def test_plot_no_metadata(self):
        fig = self.min_ord_results.plot()
        self.check_basic_figure_sanity(fig, 1, '', False, '0', '1', '2')

    def test_plot_with_numeric_metadata_and_plot_options(self):
        fig = self.min_ord_results.plot(
            self.df, 'numeric', axes=(1, 0, 2),
            axis_labels=['PC 2', 'PC 1', 'PC 3'], title='a title', cmap='Reds')
        self.check_basic_figure_sanity(
            fig, 2, 'a title', False, 'PC 2', 'PC 1', 'PC 3')

    def test_plot_with_categorical_metadata_and_plot_options(self):
        fig = self.min_ord_results.plot(
            self.df, 'categorical', axes=[2, 0, 1], title='a title',
            cmap='Accent')
        self.check_basic_figure_sanity(fig, 1, 'a title', True, '2', '0', '1')

    def test_plot_with_invalid_axis_labels(self):
        with assert_raises_regexp(ValueError, 'axis_labels.*4'):
            self.min_ord_results.plot(axes=[2, 0, 1],
                                      axis_labels=('a', 'b', 'c', 'd'))

    def test_validate_plot_axes_valid_input(self):
        # shouldn't raise an error on valid input. nothing is returned, so
        # nothing to check here
        self.min_ord_results._validate_plot_axes(self.min_ord_results.site.T,
                                                 (1, 2, 0))

    def test_validate_plot_axes_invalid_input(self):
        # not enough dimensions
        with assert_raises_regexp(ValueError, '2 dimension\(s\)'):
            self.min_ord_results._validate_plot_axes(
                np.asarray([[0.1, 0.2, 0.3], [0.2, 0.3, 0.4]]), (0, 1, 2))

        coord_matrix = self.min_ord_results.site.T

        # wrong number of axes
        with assert_raises_regexp(ValueError, 'exactly three.*found 0'):
            self.min_ord_results._validate_plot_axes(coord_matrix, [])
        with assert_raises_regexp(ValueError, 'exactly three.*found 4'):
            self.min_ord_results._validate_plot_axes(coord_matrix,
                                                     (0, 1, 2, 3))

        # duplicate axes
        with assert_raises_regexp(ValueError, 'must be unique'):
            self.min_ord_results._validate_plot_axes(coord_matrix, (0, 1, 0))

        # out of range axes
        with assert_raises_regexp(ValueError, 'axes\[1\].*3'):
            self.min_ord_results._validate_plot_axes(coord_matrix, (0, -1, 2))
        with assert_raises_regexp(ValueError, 'axes\[2\].*3'):
            self.min_ord_results._validate_plot_axes(coord_matrix, (0, 2, 3))

    def test_get_plot_point_colors_invalid_input(self):
        # column provided without df
        with npt.assert_raises(ValueError):
            self.min_ord_results._get_plot_point_colors(None, 'numeric',
                                                        ['B', 'C'], 'jet')

        # df provided without column
        with npt.assert_raises(ValueError):
            self.min_ord_results._get_plot_point_colors(self.df, None,
                                                        ['B', 'C'], 'jet')

        # column not in df
        with assert_raises_regexp(ValueError, 'missingcol'):
            self.min_ord_results._get_plot_point_colors(self.df, 'missingcol',
                                                        ['B', 'C'], 'jet')

        # id not in df
        with assert_raises_regexp(ValueError, 'numeric'):
            self.min_ord_results._get_plot_point_colors(
                self.df, 'numeric', ['B', 'C', 'missingid', 'A'], 'jet')

        # missing data in df
        with assert_raises_regexp(ValueError, 'nancolumn'):
            self.min_ord_results._get_plot_point_colors(self.df, 'nancolumn',
                                                        ['B', 'C', 'A'], 'jet')

    def test_get_plot_point_colors_no_df_or_column(self):
        obs = self.min_ord_results._get_plot_point_colors(None, None,
                                                          ['B', 'C'], 'jet')
        npt.assert_equal(obs, (None, None))

    def test_get_plot_point_colors_numeric_column(self):
        # subset of the ids in df
        exp = [0.0, -4.2, 42.0]
        obs = self.min_ord_results._get_plot_point_colors(
            self.df, 'numeric', ['B', 'C', 'A'], 'jet')
        npt.assert_almost_equal(obs[0], exp)
        assert_true(obs[1] is None)

        # all ids in df
        exp = [0.0, 42.0, 42.19, -4.2]
        obs = self.min_ord_results._get_plot_point_colors(
            self.df, 'numeric', ['B', 'A', 'D', 'C'], 'jet')
        npt.assert_almost_equal(obs[0], exp)
        assert_true(obs[1] is None)

    def test_get_plot_point_colors_categorical_column(self):
        # subset of the ids in df
        exp_colors = [[0., 0., 0.5, 1.], [0., 0., 0.5, 1.], [0.5, 0., 0., 1.]]
        exp_color_dict = {
            'foo': [0.5, 0., 0., 1.],
            22: [0., 0., 0.5, 1.]
        }
        obs = self.min_ord_results._get_plot_point_colors(
            self.df, 'categorical', ['B', 'C', 'A'], 'jet')
        npt.assert_almost_equal(obs[0], exp_colors)
        npt.assert_equal(obs[1], exp_color_dict)

        # all ids in df
        exp_colors = [[0., 0., 0.5, 1.], [0.5, 0., 0., 1.], [0.5, 0., 0., 1.],
                      [0., 0., 0.5, 1.]]
        obs = self.min_ord_results._get_plot_point_colors(
            self.df, 'categorical', ['B', 'A', 'D', 'C'], 'jet')
        npt.assert_almost_equal(obs[0], exp_colors)
        # should get same color dict as before
        npt.assert_equal(obs[1], exp_color_dict)

    def test_plot_categorical_legend(self):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # we shouldn't have a legend yet
        assert_true(ax.get_legend() is None)

        self.min_ord_results._plot_categorical_legend(
            ax, {'foo': 'red', 'bar': 'green'})

        # make sure we have a legend now
        legend = ax.get_legend()
        assert_true(legend is not None)

        # do some light sanity checking to make sure our input labels and
        # colors are present. we're not using nose.tools.assert_items_equal
        # because it isn't available in Python 3.
        labels = [t.get_text() for t in legend.get_texts()]
        npt.assert_equal(sorted(labels), ['bar', 'foo'])

        colors = [l.get_color() for l in legend.get_lines()]
        npt.assert_equal(sorted(colors), ['green', 'red'])

    def test_repr_png(self):
        obs = self.min_ord_results._repr_png_()
        assert_is_instance(obs, binary_type)
        assert_true(len(obs) > 0)

    def test_repr_svg(self):
        obs = self.min_ord_results._repr_svg_()
        assert_is_instance(obs, text_type)
        assert_true(len(obs) > 0)

    def test_png(self):
        assert_is_instance(self.min_ord_results.png, Image)

    def test_svg(self):
        assert_is_instance(self.min_ord_results.svg, SVG)


if __name__ == '__main__':
    import nose
    nose.runmodule()
