#! /usr/bin/env python

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import division

import warnings
from StringIO import StringIO
from itertools import izip

import numpy as np
import numpy.testing as npt
from scipy.spatial.distance import pdist

from skbio.maths.stats.ordination import CA, RDA, CCA, PCoA, OrdinationResults
from skbio.core.distance import DistanceMatrix
from skbio.util.testing import get_data_path
from skbio.core.exception import FileFormatError


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

        # Note the absolute value because column can have signs swapped
        npt.assert_almost_equal(np.abs(scores.site[0, 0]),
                                0.24078813304509292)

        # cogent returned the scores transposed
        npt.assert_almost_equal(np.abs(scores.site[0, 1]),
                                0.23367716219400031)


class TestPCoAResultsExtensive(object):
    def setup(self):
        matrix = np.loadtxt(get_data_path('PCoA_sample_data_2'))
        self.ids = map(str, range(matrix.shape[0]))
        dist_matrix = DistanceMatrix(matrix, self.ids)
        self.ordination = PCoA(dist_matrix)

    def test_values(self):
        results = self.ordination.scores()

        npt.assert_almost_equal(len(results.eigvals), len(results.site[0]))

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
        with open(get_data_path('PCoA_sample_data_3'), 'U') as lines:
            dist_matrix = DistanceMatrix.from_file(lines)
        self.ordination = PCoA(dist_matrix)
        self.ids = ['PC.636', 'PC.635', 'PC.356', 'PC.481', 'PC.354', 'PC.593',
                    'PC.355', 'PC.607', 'PC.634']

    def test_values(self):
        results = self.ordination.scores()

        npt.assert_almost_equal(len(results.eigvals), len(results.site[0]))

        expected = np.array([[-0.25846546, 0.17399955, 0.03828758, -0.19447751,
                              0.0831176, 0.26243033, -0.02316364, -0.0184794,
                              0.0],
                             [-0.27100114, -0.01859513, -0.08648419,
                              0.11806425, -0.19880836, -0.02117236,
                              -0.19102403, 0.15564659, 0.0],
                             [0.2350779, 0.09625193, -0.34579273, -0.00320863,
                              -0.09637777, 0.04570254, 0.18547281, 0.0404094,
                              0.0],
                             [0.02614077, -0.01114597, 0.1476606, 0.29087661,
                              0.20394547, 0.06197124, 0.10164133, 0.105691,
                              0.0],
                             [0.28500755, -0.01925499, 0.06232634, 0.1381268,
                              -0.1047986, 0.09517207, -0.1296361, -0.22068717,
                              0.0],
                             [0.20463633, -0.13936115, 0.29151382, -0.18156679,
                              -0.15958013, -0.02464121, 0.08662524,
                              0.09962215, 0.0],
                             [0.2334824, 0.22525797, -0.01886231, -0.10772998,
                              0.177109, -0.19290584, -0.14981947, 0.0383549,
                              0.0],
                             [-0.09496319, -0.4209748, -0.15486945,
                              -0.08984275, 0.15261819, -0.03342327,
                              -0.02512248, -0.05089885, 0.0],
                             [-0.35991516, 0.1138226, 0.06622034, 0.029758,
                              -0.05722541, -0.19313351, 0.14502633,
                              -0.14965861, 0.0]])
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


class TestOrdinationResults(object):
    def setup(self):
        # CA results
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
        ca_scores = OrdinationResults(eigvals=eigvals, species=species,
                                      site=site, biplot=biplot,
                                      site_constraints=site_constraints,
                                      proportion_explained=prop_explained,
                                      species_ids=species_ids,
                                      site_ids=site_ids)
        # CCA results
        eigvals = np.array([0.366135830393, 0.186887643052, 0.0788466514249,
                            0.082287840501, 0.0351348475787, 0.0233265839374,
                            0.0099048981912, 0.00122461669234,
                            0.000417454724117])
        species = np.array([[0.110350890177, 0.282399990052, -0.203028976154,
                             -0.00192462284409, -0.082232863384,
                             0.0857314258364, -0.0122038907184,
                             -0.0425198793666, 0.00466719926338],
                            [0.141359038961, 0.303495645402, 0.395441211576,
                             -0.14126625534, -0.0268859204718, 0.143253061936,
                             0.0430260301697, 0.0476377655759,
                             -0.00228172378295],
                            [-1.01552204222, 0.0958317865043, -0.198262718034,
                             -0.104801030067, 0.130025239749, 0.0244045261332,
                             0.0464701211285, 0.0269279200532,
                             0.0350103177576],
                            [-1.03620650502, 0.109624974112, 0.220984718362,
                             0.223640072997, -0.243745876054, -0.0259064859794,
                             -0.0534088909011, -0.0315611195993,
                             0.0256448427552],
                            [1.05371722248, 0.537178749104, -0.438075060322,
                             0.223480553581, -0.323948461806, 0.124644870822,
                             -0.119275907223, 0.0416254660329,
                             -0.0381955235096],
                            [0.998558655, 0.573960582723, 0.679918103399,
                             -0.389963380717, 0.299077945999, 0.328451006171,
                             0.21215881857, -0.0829871883001,
                             -0.0439653996462],
                            [0.255245719879, -0.178168259149, -0.204127155429,
                             0.433397565801, 0.0707099230629, -0.18817306522,
                             0.126908756045, 0.0044937289123,
                             -0.0122511718244],
                            [0.146555872394, -0.857362497037, -0.0152499051659,
                             0.0527604990862, 0.354475793915, -0.0416813697787,
                             -0.199011239586, -0.00213723187073,
                             -0.00782946141667],
                            [0.413705102117, -0.707948964322, 0.21569736034,
                             -0.690314241725, -0.148431001217, -0.334251868558,
                             -0.00628707445028, -0.00364123416731,
                             -0.0122722164511]])
        site = np.array([[0.710587311248, -3.08166800613, 0.219651379947,
                          -1.24528801163, -1.07293546227, -0.506241907472,
                          0.244126652455, -3.63164833508, 1.16311896657],
                         [0.584771352278, -3.00669301091, -0.947448656768,
                          2.69965142856, 2.13682885838, 0.813520011254,
                          0.471530333182, 0.908423015086, -1.34724387844],
                         [0.762734278287, -3.15258603503, 2.13924426714,
                          -3.1162748358, -2.30660936925, -0.698929858809,
                          -1.39062619586, 4.84117591747, 0.562102984837],
                         [1.11230735331, 1.07150585141, -1.87527740873,
                          0.666370241998, -1.10153224699, 1.43517552491,
                          -1.10619960297, 0.0137029328454, -0.0371803939101],
                         [-0.979116769996, -0.0603144289026, -0.696277367656,
                          -0.612646703308, 0.983006619615, 0.315662442163,
                          0.574110232297, 0.328630035672, 0.868027697443],
                         [1.04322560423, 0.459426970165, -0.639802790578,
                          0.287156643872, -0.573935423877, -1.44980634943,
                          1.70166994063, 0.306164261447, -0.442115969758],
                         [-0.954490118162, -0.0847021660539, 0.132509124411,
                          -0.42143341064, -0.111552348931, -0.394242454835,
                          -0.673963982894, -0.379018566362, -1.7472502885],
                         [0.947268764751, -0.108370567311, 0.526107182587,
                          -0.00565282365567, 1.26272400228, -1.06565692165,
                          -1.46326596729, -0.154459216567, 0.778139732463],
                         [-1.14808173207, 0.490449274267, 0.478353666755,
                          1.17015870919, -1.00599224074, 0.0735071441404,
                          0.0860462673715, 0.0417647558417, 0.935819560428],
                         [1.03291557934, 1.0350490304, 2.74691777314,
                          -1.28083971649, 0.363002636972, 1.98647950015,
                          1.05356145232, -0.24813142226, -0.463165215106]])
        biplot = np.array([[-0.169746767979, 0.63069090084, 0.760769036049],
                           [-0.994016563505, 0.0609533148724,
                            -0.0449369418179],
                           [0.184352565909, -0.974867543612, 0.0309865007541]])
        site_constraints = np.array([[0.692138797603, -3.08053663489,
                                      -0.328747278055, -1.24528801163,
                                      -1.07293546227, -0.506241907472,
                                      0.244126652455, -3.63164833508,
                                      1.16311896657],
                                     [0.664559513865, -3.06214571808,
                                      0.230249303805, 2.69965142856,
                                      2.13682885838, 0.813520011254,
                                      0.471530333182, 0.908423015086,
                                      -1.34724387844],
                                     [0.636980230127, -3.04375480127,
                                      0.789245885666, -3.1162748358,
                                      -2.30660936925, -0.698929858809,
                                      -1.39062619586, 4.84117591747,
                                      0.562102984837],
                                     [1.10887578995, 0.500396915484,
                                      -1.55606822404, 0.666370241998,
                                      -1.10153224699, 1.43517552491,
                                      -1.10619960297, 0.0137029328454,
                                      -0.0371803939101],
                                     [-0.970016224052, 0.0654867737684,
                                      -1.1206070781, -0.612646703308,
                                      0.983006619615, 0.315662442163,
                                      0.574110232297, 0.328630035672,
                                      0.868027697443],
                                     [1.05371722248, 0.537178749104,
                                      -0.438075060322, 0.287156643872,
                                      -0.573935423877, -1.44980634943,
                                      1.70166994063, 0.306164261447,
                                      -0.442115969758],
                                     [-1.02517479153, 0.102268607388,
                                      -0.00261391438256, -0.42143341064,
                                      -0.111552348931, -0.394242454835,
                                      -0.673963982894, -0.379018566362,
                                      -1.7472502885],
                                     [0.998558655, 0.573960582723,
                                      0.679918103399, -0.00565282365567,
                                      1.26272400228, -1.06565692165,
                                      -1.46326596729, -0.154459216567,
                                      0.778139732463],
                                     [-1.080333359, 0.139050441007,
                                      1.11537924934, 1.17015870919,
                                      -1.00599224074, 0.0735071441404,
                                      0.0860462673715, 0.0417647558417,
                                      0.935819560428],
                                     [0.943400087524, 0.610742416342,
                                      1.79791126712, -1.28083971649,
                                      0.363002636972, 1.98647950015,
                                      1.05356145232, -0.24813142226,
                                      -0.463165215106]])
        prop_explained = None
        species_ids = ['Species0', 'Species1', 'Species2', 'Species3',
                       'Species4', 'Species5', 'Species6', 'Species7',
                       'Species8']
        site_ids = ['Site0', 'Site1', 'Site2', 'Site3', 'Site4', 'Site5',
                    'Site6', 'Site7', 'Site8', 'Site9']
        cca_scores = OrdinationResults(eigvals=eigvals, species=species,
                                       site=site, biplot=biplot,
                                       site_constraints=site_constraints,
                                       proportion_explained=prop_explained,
                                       species_ids=species_ids,
                                       site_ids=site_ids)
        # PCoA results
        eigvals = np.array([0.512367260461, 0.300719094427, 0.267912066004,
                            0.208988681078, 0.19169895326, 0.16054234528,
                            0.15017695712, 0.122457748167, 0.0])
        species = None
        site = np.array([[-0.258465461183, 0.173999546883, 0.0382875792552,
                          -0.19447750562, 0.0831176020844, 0.262430333201,
                          -0.0231636392235, -0.0184794039581, 0.0],
                         [-0.271001135391, -0.0185951319063, -0.0864841926349,
                          0.118064245315, -0.198808358437, -0.0211723599535,
                          -0.191024027565, 0.155646592377, 0.0],
                         [0.235077898175, 0.0962519254489, -0.345792726714,
                          -0.00320862577619, -0.0963777675519, 0.0457025386953,
                          0.185472813286, 0.0404093971793, 0.0],
                         [0.0261407664325, -0.0111459676533, 0.147660603015,
                          0.29087660853, 0.203945472801, 0.0619712384758,
                          0.101641328709, 0.105690998719, 0.0],
                         [0.285007552283, -0.0192549888483, 0.0623263375385,
                          0.138126799852, -0.104798602423, 0.0951720730628,
                          -0.129636097542, -0.220687170372, 0.0],
                         [0.204636326241, -0.139361150932, 0.291513819623,
                          -0.181566786821, -0.159580132715, -0.0246412130162,
                          0.0866252404441, 0.0996221476871, 0.0],
                         [0.233482403212, 0.225257974068, -0.0188623096268,
                          -0.107729981831, 0.177108999572, -0.192905835151,
                          -0.149819471408, 0.0383549037465, 0.0],
                         [-0.0949631911323, -0.420974802495, -0.154869454869,
                          -0.0898427509281, 0.152618194488, -0.0334232691501,
                          -0.0251224777303, -0.0508988536409, 0.0],
                         [-0.359915158638, 0.113822595435, 0.0662203444138,
                          0.0297579972788, -0.0572254078183, -0.193133506163,
                          0.145026331031, -0.149658611738, 0.0]])
        biplot = None
        site_constraints = None
        prop_explained = np.array([0.267573832777, 0.15704469605,
                                   0.139911863774, 0.109140272454,
                                   0.100111048503, 0.0838401161912,
                                   0.0784269939011, 0.0639511763509, 0.0])
        species_ids = None
        site_ids = ['PC.636', 'PC.635', 'PC.356', 'PC.481', 'PC.354', 'PC.593',
                    'PC.355', 'PC.607', 'PC.634']
        pcoa_scores = OrdinationResults(eigvals=eigvals, species=species,
                                        site=site, biplot=biplot,
                                        site_constraints=site_constraints,
                                        proportion_explained=prop_explained,
                                        species_ids=species_ids,
                                        site_ids=site_ids)
        # RDA results
        eigvals = np.array([25.8979540892, 14.9825779819, 8.93784077262,
                            6.13995623072, 1.68070536498, 0.57735026919,
                            0.275983624351])
        species = np.array([[1.38198713901, -1.71496426179, 0.632272455288,
                             0.00712898231575, 0.120512431133,
                             -0.0723104306179, -0.00815886062344],
                            [0.919178380672, -1.25430767906, -1.1787426896,
                             -0.00712898231576, -0.120512431133,
                             -0.0723104306179, 0.00815886062344],
                            [3.39897234869, 0.446168315515, 0.406691610423,
                             0.749336668014, 0.0793892812781,
                             7.37971401683e-17, 0.0329418170936],
                            [2.52353261895, 0.446932822723, -0.413412046583,
                             -0.639449029945, -0.0640330006084,
                             3.40602185392e-17, 0.0335491330226],
                            [-0.53155341411, -1.34263985744, 0.464155649196,
                             -0.412041388665, 0.198336195195,
                             7.37971401683e-17, 0.00604836743485],
                            [-0.288618167117, -0.571491852197, -0.406527290424,
                             0.206020694333, -0.0991680975973,
                             -1.13534061797e-17, -0.00302418371743]])
        site = np.array([[-1.48848983495, 2.12675623514, 0.727805340002,
                          -0.227234564008, -3.8413042049, -2.30487725273,
                          0.260061682644],
                         [-1.5541678384, 2.37027298265, 0.475523558326,
                          1.58712629997e-16, 1.39853499536e-15, 4.60975450547,
                          -1.41948353841e-14],
                         [-1.51048450796, 2.19216727329, 0.00519576944216,
                          0.227234564008, 3.8413042049, -2.30487725273,
                          -0.260061682644],
                         [-0.872786591764, -2.6271708553, 2.68871897067,
                          -1.97005774092, 0.948287641474, -2.0356145959e-14,
                          0.0289185344306],
                         [2.97228673755, 0.322310666722, 2.50294580667,
                          3.50264153009, 0.489477682536, -1.25529566747e-14,
                          2.11938273809],
                         [-0.879968888341, -2.19620098193, 0.710888524695,
                          -0.656685913639, 0.316095880491, -4.47835211098e-15,
                          0.00963951147681],
                         [2.64194948913, 0.390104638861, -0.086230363198,
                          -0.211189359785, -0.298609965083, -3.88762243221e-15,
                          -4.5952222736],
                         [-0.887151184918, -1.76523110855, -1.26694192128,
                          0.656685913639, -0.316095880491, 1.21458337555e-14,
                          -0.00963951147698],
                         [2.47314610115, 0.521252384288, -2.51313331808,
                          -3.29145217031, -0.190867717454, 1.65563320466e-14,
                          2.4758395355],
                         [-0.894333481495, -1.33426123517, -3.24477236725,
                          1.97005774092, -0.948287641474, 3.0262803659e-14,
                          -0.0289185344308]])
        biplot = np.array([[0.422650019179, -0.559142585857, -0.713250678211],
                           [0.988495963777, 0.150787422017, -0.0117848614073],
                           [-0.556516618887, 0.817599992718, 0.147714267459],
                           [-0.404079676685, -0.9058434809, -0.127150316558]])
        site_constraints = np.array([[-1.48131076339, 2.07063239013,
                                      1.42061063192, -0.227234564008,
                                      -3.8413042049, -2.30487725273,
                                      0.260061682644],
                                     [-1.51771406044, 2.22973216369,
                                      0.402841555923, 1.58712629997e-16,
                                      1.39853499536e-15, 4.60975450547,
                                      -1.41948353841e-14],
                                     [-1.55411735749, 2.38883193726,
                                      -0.61492752007, 0.227234564008,
                                      3.8413042049, -2.30487725273,
                                      -0.260061682644],
                                     [-0.774350145471, -2.45801536594,
                                      2.77528052969, -1.97005774092,
                                      0.948287641474, -2.0356145959e-14,
                                      0.0289185344306],
                                     [2.76860070338, 0.0930230161545,
                                      2.00339886045, 3.50264153009,
                                      0.489477682536, -1.25529566747e-14,
                                      2.11938273809],
                                     [-0.847156739577, -2.13981581881,
                                      0.739742377702, -0.656685913639,
                                      0.316095880491, -4.47835211098e-15,
                                      0.00963951147681],
                                     [2.69579410928, 0.41122256329,
                                      -0.0321392915344, -0.211189359785,
                                      -0.298609965083, -3.88762243221e-15,
                                      -4.5952222736],
                                     [-0.919963333683, -1.82161627167,
                                      -1.29579577429, 0.656685913639,
                                      -0.316095880491, 1.21458337555e-14,
                                      -0.00963951147698],
                                     [2.62298751517, 0.729422110426,
                                      -2.06767744352, -3.29145217031,
                                      -0.190867717454, 1.65563320466e-14,
                                      2.4758395355],
                                     [-0.992769927788, -1.50341672453,
                                      -3.33133392627, 1.97005774092,
                                      -0.948287641474, 3.0262803659e-14,
                                      -0.0289185344308]])
        prop_explained = None
        species_ids = ['Species0', 'Species1', 'Species2', 'Species3',
                       'Species4', 'Species5']
        site_ids = ['Site0', 'Site1', 'Site2', 'Site3', 'Site4', 'Site5',
                    'Site6', 'Site7', 'Site8', 'Site9']
        rda_scores = OrdinationResults(eigvals=eigvals, species=species,
                                       site=site, biplot=biplot,
                                       site_constraints=site_constraints,
                                       proportion_explained=prop_explained,
                                       species_ids=species_ids,
                                       site_ids=site_ids)

        self.scores = [ca_scores, cca_scores, pcoa_scores, rda_scores]
        self.test_paths = ['L&L_CA_data_scores', 'example3_scores',
                           'PCoA_sample_data_3_scores', 'example2_scores']

        self.fferror_test_paths = ['error1', 'error2', 'error3', 'error4',
                                   'error5', 'error6']
        self.verror_test_paths = ['v_error1', 'v_error2', 'v_error3',
                                  'v_error4', 'v_error5', 'v_error6',
                                  'v_error7', 'v_error8', 'v_error9',
                                  'v_error10']

    def test_to_file(self):
        for scores, test_path in izip(self.scores, self.test_paths):
            obs_f = StringIO()
            scores.to_file(obs_f)
            obs = obs_f.getvalue()
            obs_f.close()

            with open(get_data_path(test_path), 'U') as f:
                exp = f.read()

            npt.assert_equal(obs, exp)

    def test_from_file(self):
        for scores, test_path in izip(self.scores, self.test_paths):
            with open(get_data_path(test_path), 'U') as f:
                obs = OrdinationResults.from_file(f)

            npt.assert_almost_equal(obs.eigvals, scores.eigvals)
            if scores.species is not None:
                npt.assert_almost_equal(obs.species, scores.species)
            else:
                npt.assert_equal(obs.species, scores.species)
            npt.assert_equal(obs.species_ids, scores.species_ids)

            if scores.site is not None:
                npt.assert_almost_equal(obs.site, scores.site)
            else:
                npt.assert_equal(obs.site, scores.site)
            npt.assert_equal(obs.site_ids, scores.site_ids)

            if scores.biplot is not None:
                npt.assert_almost_equal(obs.biplot, scores.biplot)
            else:
                npt.assert_equal(obs.biplot, scores.biplot)

            if scores.site_constraints is not None:
                npt.assert_almost_equal(obs.site_constraints,
                                        scores.site_constraints)
            else:
                npt.assert_equal(obs.site_constraints, scores.site_constraints)

            if scores.proportion_explained is not None:
                npt.assert_almost_equal(obs.proportion_explained,
                                        scores.proportion_explained)
            else:
                npt.assert_equal(obs.proportion_explained,
                                 scores.proportion_explained)

    def test_from_file_error(self):
        for test_path in self.fferror_test_paths:
            with open(get_data_path(test_path), 'U') as f:
                with npt.assert_raises(FileFormatError):
                    OrdinationResults.from_file(f)

        for test_path in self.verror_test_paths:
            with open(get_data_path(test_path), 'U') as f:
                with npt.assert_raises(ValueError):
                    OrdinationResults.from_file(f)
