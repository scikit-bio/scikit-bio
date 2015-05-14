# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import warnings

import numpy as np
import numpy.testing as npt
from scipy.spatial.distance import pdist

from skbio import DistanceMatrix, OrdinationResults
from skbio.stats.ordination import CA, rda, CCA, PCoA, corr, mean_and_std
from skbio.util import get_data_path, assert_ordination_results_equal
import pandas as pd

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
            yield npt.assert_raises, ValueError, rda, Y, X, None, None


class TestRDAResults(object):
    # STATUS: L&L only shows results with scaling 1, and they agree
    # with vegan's (module multiplying by a constant). I can also
    # compute scaling 2, agreeing with vegan, but there are no written
    # results in L&L.
    def setup(self):
        """Data from table 11.3 in Legendre & Legendre 1998."""
        sample_ids = ['Site0', 'Site1', 'Site2', 'Site3', 'Site4',
                      'Site5', 'Site6', 'Site7', 'Site8', 'Site9']
        features_ids = ['Species0', 'Species1', 'Species2', 'Species3',
                       'Species4', 'Species5']

        self.Y = pd.DataFrame(np.loadtxt(get_data_path('example2_Y')),
                              index=sample_ids,columns=features_ids)
        self.X = pd.DataFrame(np.loadtxt(get_data_path('example2_X')),
                              index=sample_ids)
    def test_scaling1(self):

        scores = rda(self.Y, self.X, scaling=1)
        pc_ids = ['RDA1', 'RDA2', 'RDA3', 'RDA4', 'RDA5', 'RDA6', 'RDA7']
        biplot_scores = np.array(
            [[ 0.422650019, -0.559142586, -0.713250678, 1.165734176e-16, 1.471045508e-16  , 1.831867991e-16],
             [ 0.988495964, 0.150787422, -0.011784861, 6.106226635e-17, 6.661338148e-17   , 8.326672685e-17],
             [ -0.556516619, 0.817599993, 0.147714267, -4.996003611e-17, 4.440892099e-17  , -7.216449660e-17],
             [ -0.404079677, -0.905843481, -0.127150317, 2.775557562e-18, -2.220446049e-17, 0.000000000e+00]])

        sample_constraints = np.array(
            [[-1.203551785, 0.973290974, 0.398346330, -4.377163939e-02, -2.025458896e-01, -4.174844658e-02, 2.251711925e-03],
             [-1.233129139, 1.048075071, 0.112958788, 1.946349502e-16, -3.553871934e-16, 8.349689316e-02, -1.554395167e-16],
             [-1.262706493, 1.122859169, -0.172428754, 4.377163939e-02,  2.025458896e-01, -4.174844658e-02, -2.251711925e-03],
             [-0.629152587, -1.155378512, 0.778202548, -3.794874137e-01,  5.000170610e-02, 3.937851438e-16, 2.503875549e-04],
             [ 2.249463380, 0.043725029, 0.561763065, 6.747052880e-01,  2.580938331e-02, 6.726138671e-16, 1.835041340e-02],
             [-0.688307296, -1.005810318, 0.207427464, -1.264958046e-01,  1.666723537e-02, -6.333664505e-17, 8.346251830e-05],
             [ 2.190308672, 0.193293223, -0.009012018, -4.068089086e-02, -1.574523073e-02, -6.651371118e-18, -3.978716391e-02],
             [-0.747462004, -0.856242124, -0.363347619, 1.264958046e-01, -1.666723537e-02, -4.098446294e-16, -8.346251830e-05],
             [ 2.131153964, 0.342861418, -0.579787102, -6.340243972e-01, -1.006415258e-02, -4.849800803e-16, 2.143675051e-02],
             [-0.806616713, -0.706673930, -0.934122703, 3.794874137e-01, -5.000170610e-02, -7.280846416e-16, -2.503875549e-04]])

         # Load data as computed with vegan 2.0-8
        vegan_features = np.loadtxt(get_data_path(
            'example2_species_scaling1_from_vegan'))
        #npt.assert_almost_equal(scores.features, vegan_features, decimal=6)

        vegan_samples = np.loadtxt(get_data_path(
            'example2_site_scaling1_from_vegan'))
        #npt.assert_almost_equal(scores.samples, vegan_samples, decimal=6)


        exp = OrdinationResults(
            'RDA','Redundancy Analysis',
            samples=pd.DataFrame(vegan_samples,
                                 index=['Site0', 'Site1', 'Site2',
                                        'Site3', 'Site4',
                                        'Site5', 'Site6', 'Site7',
                                        'Site8', 'Site9'],
                                 columns=pc_ids),
            features=pd.DataFrame(vegan_features,
                                  index=['Species0', 'Species1',
                                         'Species2', 'Species3',
                                         'Species4', 'Species5'],
                                  columns=pc_ids),
            sample_constraints=pd.DataFrame(sample_constraints,
                                            index=['Site0', 'Site1', 'Site2',
                                                   'Site3', 'Site4',
                                                   'Site5', 'Site6', 'Site7',
                                                   'Site8', 'Site9'],
                                            columns=pc_ids),
            biplot_scores = pd.DataFrame(biplot_scores),
            proportion_explained = pd.Series([0.44275783, 0.25614586,
                                             0.15280354, 0.10497021,
                                             0.02873375, 0.00987052,
                                             0.00471828], index=pc_ids),
            eigvals = pd.Series([25.897954, 14.982578, 8.937841, 6.139956,
                                 1.680705, 0.577350, 0.275984], index=pc_ids)
            )

        npt.assert_almost_equal(scores.features, vegan_features, decimal=6)
        npt.assert_almost_equal(scores.samples, vegan_samples, decimal=6)
        assert_ordination_results_equal(scores, exp, ignore_biplot_scores_labels=True)

    def test_scaling2(self):
        scores = rda(self.Y, self.X, scaling=2)

        # Load data as computed with vegan 2.0-8
        vegan_features = np.loadtxt(get_data_path(
            'example2_species_scaling2_from_vegan'))
        vegan_samples = np.loadtxt(get_data_path(
            'example2_site_scaling2_from_vegan'))

        sample_constraints = np.array(
                [[-1.48131076e+00, 2.07063239e+00, 1.42061063e+00,
                  -2.27234564e-01, -3.84130420e+00, -2.30487725e+00,
                  2.60061683e-01],
                 [-1.51771406e+00, 2.22973216e+00, 4.02841556e-01,
                  1.01042110e-15, -6.73995569e-15, 4.60975451e+00,
                  -1.79525017e-14],
                 [-1.55411736e+00, 2.38883194e+00, -6.14927520e-01,
                   2.27234564e-01, 3.84130420e+00, -2.30487725e+00,
                  -2.60061683e-01],
                 [-7.74350145e-01, -2.45801537e+00, 2.77528053e+00,
                  -1.97005774e+00, 9.48287641e-01, 2.17403639e-14,
                  2.89185344e-02],
                 [2.76860070e+00, 9.30230162e-02, 2.00339886e+00,
                  3.50264153e+00, 4.89477683e-01, 3.71341338e-14,
                  2.11938274e+00],
                 [-8.47156740e-01, -2.13981582e+00, 7.39742378e-01,
                  -6.56685914e-01, 3.16095880e-01, -3.49673352e-15,
                  9.63951148e-03],
                 [2.69579411e+00, 4.11222563e-01, -3.21392915e-02,
                  -2.11189360e-01, -2.98609965e-01, -3.67213519e-16,
                  -4.59522227e+00],
                 [-9.19963334e-01, -1.82161627e+00, -1.29579577e+00,
                  6.56685914e-01, -3.16095880e-01, -2.26269871e-14,
                  -9.63951148e-03],
                 [2.62298752e+00, 7.29422110e-01, -2.06767744e+00,
                  -3.29145217e+00, -1.90867717e-01, -2.67751173e-14,
                  2.47583954e+00],
                 [-9.92769928e-01, -1.50341672e+00, -3.33133393e+00,
                  1.97005774e+00, -9.48287641e-01, -4.01966029e-14,
                  -2.89185344e-02]])
        biplot_scores = pd.DataFrame(
            [[4.22650019e-01, -5.59142586e-01, -7.13250678e-01,
              1.16573418e-16, 1.47104551e-16, 1.83186799e-16],
             [9.88495964e-01, 1.50787422e-01, -1.17848614e-02,
              6.10622664e-17, 6.66133815e-17, 8.32667268e-17],
             [-5.56516619e-01, 8.17599993e-01, 1.47714267e-01,
              -4.99600361e-17, 4.44089210e-17, -7.21644966e-17],
             [-4.04079677e-01, -9.05843481e-01, -1.27150317e-01,
              2.77555756e-18, -2.22044605e-17, 0.00000000e+00]])
        pc_ids = ['RDA1', 'RDA2', 'RDA3', 'RDA4', 'RDA5', 'RDA6', 'RDA7']
        exp = OrdinationResults(
            'RDA','Redundancy Analysis',
            samples=pd.DataFrame(vegan_samples,
                                 index=['Site0', 'Site1', 'Site2',
                                        'Site3', 'Site4',
                                        'Site5', 'Site6', 'Site7',
                                        'Site8', 'Site9'],
                                 columns=pc_ids),
            features=pd.DataFrame(vegan_features,
                                  index=['Species0', 'Species1',
                                         'Species2', 'Species3',
                                         'Species4', 'Species5'],
                                  columns=pc_ids),
            sample_constraints=pd.DataFrame(sample_constraints,
                                            index=['Site0', 'Site1', 'Site2',
                                                   'Site3', 'Site4',
                                                   'Site5', 'Site6', 'Site7',
                                                   'Site8', 'Site9'],
                                            columns=pc_ids),
            biplot_scores = pd.DataFrame(biplot_scores),
            proportion_explained = pd.Series([0.44275783, 0.25614586,
                                             0.15280354, 0.10497021,
                                             0.02873375, 0.00987052,
                                             0.00471828], index=pc_ids),
            eigvals = pd.Series([25.89795409, 14.98257798, 8.93784077, 6.13995623,
                                 1.68070536, 0.57735027, 0.27598362],
                                 index=pc_ids)

            )
        npt.assert_almost_equal(scores.features, vegan_features, decimal=6)
        npt.assert_almost_equal(scores.samples, vegan_samples, decimal=6)
        assert_ordination_results_equal(scores, exp, ignore_biplot_scores_labels=True)

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


if __name__ == '__main__':
    import nose
    nose.runmodule()
