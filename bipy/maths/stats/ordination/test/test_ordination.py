#! /usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013--, bipy development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import division

import os
import warnings

import numpy as np
import numpy.testing as npt
from scipy.spatial.distance import pdist

from bipy.maths.stats.ordination import CA, RDA, CCA, PCoA
from bipy.core.distance import SymmetricDistanceMatrix


def get_data_path(fn):
    """Return path to filename `fn` in the data folder."""
    path = os.path.dirname(os.path.abspath(__file__))
    data_path = os.path.join(path, 'data', fn)
    return data_path


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
    """Computes the chi-square distance between two rows or columns of
    input.

    It is a measure that has no upper limit, and it excludes
    double-zeros.

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
    `scipy.spatial.distance.squareform` for a routine that transforms
    condensed distance matrices to square-form distance matrices (and
    vice-versa).

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


class TestCAResults(object):
    def setup(self):
        """Data from table 9.11 in Legendre & Legendre 1998."""
        self.X = np.loadtxt(get_data_path('L&L_CA_data'))
        self.ordination = CA(self.X)

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
            CA(X)


class TestRDAErrors(object):
    def test_shape(self):
        for n, p, n_, m in [(3, 4, 2, 1), (3, 4, 3, 10)]:
            Y = np.random.randn(n, p)
            X = np.random.randn(n_, m)
            yield npt.assert_raises, ValueError, RDA, Y, X


class TestRDAResults(object):
    # STATUS: L&L only shows results with scaling 1, and they agree
    # with vegan's (module multiplying by a constant). I can also
    # compute scaling 2, agreeing with vegan, but there're no written
    # results in L&L.
    def setup(self):
        """Data from table 9.11 in Legendre & Legendre 1998."""
        Y = np.loadtxt(get_data_path('example2_Y'))
        X = np.loadtxt(get_data_path('example2_X'))
        self.ordination = RDA(Y, X)

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
        """Data from table 9.11 in Legendre & Legendre 1998."""
        self.Y = np.loadtxt(get_data_path('example3_Y'))
        self.X = np.loadtxt(get_data_path('example3_X'))

    def test_shape(self):
        X, Y = self.X, self.Y
        with npt.assert_raises(ValueError):
            CCA(Y, X[:-1])

    def test_Y_values(self):
        X, Y = self.X, self.Y
        Y[0, 0] = -1
        with npt.assert_raises(ValueError):
            CCA(Y, X)
        Y[0] = 0
        with npt.assert_raises(ValueError):
            CCA(Y, X)


class TestCCAResults(object):
    def setup(self):
        """Data from table 9.11 in Legendre & Legendre 1998. Loaded
        results as computed with vegan 2.0-8"""
        Y = np.loadtxt(get_data_path('example3_Y'))
        X = np.loadtxt(get_data_path('example3_X'))
        self.ordination = CCA(Y, X)

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
        dist_matrix = SymmetricDistanceMatrix(matrix,
                                              map(str, range(matrix.shape[0])))
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

        # Note the absolute value because column can have signs swapped
        npt.assert_almost_equal(np.abs(scores.species[0, 0]),
                                0.24078813304509292)

        # cogent returned the scores transposed
        npt.assert_almost_equal(np.abs(scores.species[0, 1]),
                                0.23367716219400031)


class TestPCoAResultsExtensive(object):
    def setup(self):
        matrix = np.loadtxt(get_data_path('PCoA_sample_data_2'))
        dist_matrix = SymmetricDistanceMatrix(matrix,
                                              map(str, range(matrix.shape[0])))
        self.ordination = PCoA(dist_matrix)

    def test_values(self):
        expected = np.array([[-0.028597, 0.22903853, 0.07055272,
                              0.26163576, 0.28398669],
                             [0.37494056, 0.22334055, -0.20892914,
                              0.05057395, -0.18710366],
                             [-0.33517593, -0.23855979, -0.3099887,
                              0.11521787, -0.05021553],
                             [0.25412394, -0.4123464, 0.23343642,
                              0.06403168, -0.00482608],
                             [-0.28256844, 0.18606911, 0.28875631,
                              -0.06455635, -0.21141632],
                             [0.01727687, 0.012458, -0.07382761,
                              -0.42690292, 0.1695749]])
        actual = self.ordination.scores().species
        npt.assert_almost_equal(*normalize_signs(expected, actual))


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
