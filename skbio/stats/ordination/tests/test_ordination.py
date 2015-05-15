# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import warning

import numpy as np
import numpy.testing as npt
import pandas as pd
from scipy.spatial.distance import pdist
import pandas as pd

from skbio import OrdinationResults
from skbio.stats.ordination import (
    ca, rda, cca, corr, mean_and_std, e_matrix, f_matrix)
from skbio.util import get_data_path, assert_ordination_results_equal


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

        self.contingency = pd.DataFrame(self.X, ['Site1', 'Site2', 'Site3'],
                                        ['Species1', 'Species2', 'Species3'])

    def test_scaling2(self):

        eigvals = pd.Series(np.array([0.09613302, 0.04094181]), ['CA1', 'CA2'])
        # p. 460 L&L 1998
        features = pd.DataFrame(np.array([[0.40887, -0.06955],  # F_hat
                                          [-0.11539, 0.29977],
                                          [-0.30997, -0.18739]]),
                                ['Species1', 'Species2', 'Species3'],
                                ['CA1', 'CA2'])
        samples = pd.DataFrame(np.array([[-0.84896, -0.88276],  # V_hat
                                         [-0.22046, 1.34482],
                                         [1.66697, -0.47032]]),
                               ['Site1', 'Site2', 'Site3'],
                               ['CA1', 'CA2'])
        exp = OrdinationResults('CA', 'Correspondance Analysis',
                                eigvals=eigvals, features=features,
                                samples=samples)

        scores = ca(self.contingency, 2)

        assert_ordination_results_equal(exp, scores, decimal=5,
                                        ignore_directionality=True)

    def test_scaling1(self):
        eigvals = pd.Series(np.array([0.09613302, 0.04094181]), ['CA1', 'CA2'])
        # p. 458
        features = pd.DataFrame(np.array([[1.31871, -0.34374],  # V
                                          [-0.37215, 1.48150],
                                          [-0.99972, -0.92612]]),
                                ['Species1', 'Species2', 'Species3'],
                                ['CA1', 'CA2'])
        samples = pd.DataFrame(np.array([[-0.26322, -0.17862],  # F
                                         [-0.06835, 0.27211],
                                         [0.51685, -0.09517]]),
                               ['Site1', 'Site2', 'Site3'],
                               ['CA1', 'CA2'])
        exp = OrdinationResults('CA', 'Correspondance Analysis',
                                eigvals=eigvals, features=features,
                                samples=samples)
        scores = ca(self.contingency, 1)

        assert_ordination_results_equal(exp, scores, decimal=5,
                                        ignore_directionality=True)

    def test_maintain_chi_square_distance_scaling1(self):
        """In scaling 1, chi^2 distance among rows (sites) is equal to
        euclidean distance between them in transformed space."""
        frequencies = self.X / self.X.sum()
        chi2_distances = chi_square_distance(frequencies)
        transformed_sites = ca(self.contingency, 1).samples.values
        euclidean_distances = pdist(transformed_sites, 'euclidean')
        npt.assert_almost_equal(chi2_distances, euclidean_distances)

    def test_maintain_chi_square_distance_scaling2(self):
        """In scaling 2, chi^2 distance among columns (species) is
        equal to euclidean distance between them in transformed space."""
        frequencies = self.X / self.X.sum()
        chi2_distances = chi_square_distance(frequencies, between_rows=False)
        transformed_species = ca(self.contingency, 2).features.values
        euclidean_distances = pdist(transformed_species, 'euclidean')
        npt.assert_almost_equal(chi2_distances, euclidean_distances)


class TestCAErrors(object):
    def test_negative(self):
        X = np.array([[1, 2], [-0.1, -2]])
        with npt.assert_raises(ValueError):
            ca(pd.DataFrame(X))


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
            [[0.422650019, -0.559142586, -0.713250678,
              1.165734176e-16, 1.471045508e-16, 1.831867991e-16],
             [0.988495964, 0.150787422, -0.011784861,
              6.106226635e-17, 6.661338148e-17, 8.326672685e-17],
             [-0.556516619, 0.817599993, 0.147714267,
              -4.996003611e-17, 4.440892099e-17, -7.216449660e-17],
             [-0.404079677, -0.905843481, -0.127150317,
              2.775557562e-18, -2.220446049e-17, 0.000000000e+00]])

        sample_constraints = np.array(
            [[-1.203551785, 0.973290974, 0.398346330, -4.377163939e-02,
              -2.025458896e-01, -4.174844658e-02, 2.251711925e-03],
             [-1.233129139, 1.048075071, 0.112958788, 1.946349502e-16,
              -3.553871934e-16, 8.349689316e-02, -1.554395167e-16],
             [-1.262706493, 1.122859169, -0.172428754, 4.377163939e-02,
              2.025458896e-01, -4.174844658e-02, -2.251711925e-03],
             [-0.629152587, -1.155378512, 0.778202548, -3.794874137e-01,
              5.000170610e-02, 3.937851438e-16, 2.503875549e-04],
             [2.249463380, 0.043725029, 0.561763065, 6.747052880e-01,
              2.580938331e-02, 6.726138671e-16, 1.835041340e-02],
             [-0.688307296, -1.005810318, 0.207427464, -1.264958046e-01,
              1.666723537e-02, -6.333664505e-17, 8.346251830e-05],
             [2.190308672, 0.193293223, -0.009012018, -4.068089086e-02,
              -1.574523073e-02, -6.651371118e-18, -3.978716391e-02],
             [-0.747462004, -0.856242124, -0.363347619, 1.264958046e-01,
              -1.666723537e-02, -4.098446294e-16, -8.346251830e-05],
             [2.131153964, 0.342861418, -0.579787102, -6.340243972e-01,
              -1.006415258e-02, -4.849800803e-16, 2.143675051e-02],
             [-0.806616713, -0.706673930, -0.934122703, 3.794874137e-01,
              -5.000170610e-02, -7.280846416e-16, -2.503875549e-04]])

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
        """Data from R's vegan package"""
        self.Y = pd.read_csv(get_data_path('varespec'), sep='\t', index_col=0)
        self.X = pd.read_csv(get_data_path('varechem'), sep='\t', index_col=0)

    def test_shape(self):
        X, Y = self.X, self.Y
        with npt.assert_raises(ValueError):
            cca(Y, X[:-1])

    def test_Y_values(self):
        X, Y = self.X, self.Y
        Y[0, 0] = -1
        with npt.assert_raises(ValueError):
            cca(Y, X)
        Y[0] = 0
        with npt.assert_raises(ValueError):
            cca(Y, X)

    def test_scaling(self):
        X, Y = self.X, self.Y
        with npt.assert_raises(NotImplementedError):
            cca(Y, X, 3)


class TestCCAResults(object):
    def setup(self):
        """Data from table 11.3 in Legendre & Legendre 1998
        (p. 590). Loaded results as computed with vegan 2.0-8 and
        compared with table 11.5 if also there."""
        self.Y = pd.DataFrame(
            np.loadtxt(get_data_path('example3_Y')),
            columns=['Feature0', 'Feature1', 'Feature2', 'Feature3',
                     'Feature4', 'Feature5', 'Feature6', 'Feature7',
                     'Feature8'],
            index=['Sample0', 'Sample1', 'Sample2', 'Sample3', 'Sample4',
                   'Sample5', 'Sample6', 'Sample7', 'Sample8', 'Sample9'])
        self.X = pd.DataFrame(
            np.loadtxt(get_data_path('example3_X')),
            columns=['Constraint0', 'Constraint1',
                     'Constraint2', 'Constraint3'],
            index=['Sample0', 'Sample1', 'Sample2', 'Sample3', 'Sample4',
                   'Sample5', 'Sample6', 'Sample7', 'Sample8', 'Sample9']
            ).ix[:, :-1]

    def test_scaling1_species(self):
        scores = cca(self.Y, self.X, 1)

        vegan_species = np.loadtxt(get_data_path(
            'example3_species_scaling1_from_vegan'))
        npt.assert_almost_equal(scores.species, vegan_species, decimal=6)

    def test_scaling1_site(self):
        scores = cca(self.Y, self.X, 1)

        vegan_site = np.loadtxt(get_data_path(
            'example3_site_scaling1_from_vegan'))
        npt.assert_almost_equal(scores.site, vegan_site, decimal=4)

    def test_scaling2_species(self):
        scores = cca(self.Y, self.X, 2)

        vegan_species = np.loadtxt(get_data_path(
            'example3_species_scaling2_from_vegan'))
        npt.assert_almost_equal(scores.species, vegan_species, decimal=5)

    def test_scaling2_site(self):
        scores = cca(self.Y, self.X, 2)

        vegan_site = np.loadtxt(get_data_path(
            'example3_site_scaling2_from_vegan'))
        npt.assert_almost_equal(scores.site, vegan_site, decimal=4)


class TestPCoAPrivateMethods(object):
    def setup(self):
        self.matrix = np.arange(1, 7).reshape(2, 3)
        self.matrix2 = np.arange(1, 10).reshape(3, 3)

    def test_e_matrix(self):
        E = e_matrix(self.matrix)
        expected_E = np.array([[-0.5,  -2.,  -4.5],
                               [-8., -12.5, -18.]])
        npt.assert_almost_equal(E, expected_E)

    def test_f_matrix(self):
        F = f_matrix(self.matrix2)
        expected_F = np.zeros((3, 3))
        # Note that `test_make_F_matrix` in cogent is wrong
        npt.assert_almost_equal(F, expected_F)


if __name__ == '__main__':
    import nose
    nose.runmodule()
