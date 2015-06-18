# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import numpy as np
import numpy.testing as npt
import pandas as pd
from scipy.spatial.distance import pdist

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
        self.sample_ids = ['Site1', 'Site2', 'Site3']
        self.feature_ids = ['Species1', 'Species2', 'Species3']
        self.pc_ids = ['CA1', 'CA2']
        self.contingency = pd.DataFrame(self.X, self.sample_ids,
                                        self.feature_ids)

    def test_scaling2(self):

        eigvals = pd.Series(np.array([0.09613302, 0.04094181]), self.pc_ids)
        # p. 460 L&L 1998
        features = pd.DataFrame(np.array([[0.40887, -0.06955],  # F_hat
                                          [-0.11539, 0.29977],
                                          [-0.30997, -0.18739]]),
                                self.feature_ids,
                                self.pc_ids)
        samples = pd.DataFrame(np.array([[-0.84896, -0.88276],  # V_hat
                                         [-0.22046, 1.34482],
                                         [1.66697, -0.47032]]),
                               self.sample_ids,
                               self.pc_ids)
        exp = OrdinationResults('CA', 'Correspondance Analysis',
                                eigvals=eigvals, features=features,
                                samples=samples)

        scores = ca(self.contingency, 2)

        assert_ordination_results_equal(exp, scores, decimal=5,
                                        ignore_directionality=True)

    def test_scaling1(self):
        eigvals = pd.Series(np.array([0.09613302, 0.04094181]), self.pc_ids)
        # p. 458
        features = pd.DataFrame(np.array([[1.31871, -0.34374],  # V
                                          [-0.37215, 1.48150],
                                          [-0.99972, -0.92612]]),
                                self.feature_ids,
                                self.pc_ids)
        samples = pd.DataFrame(np.array([[-0.26322, -0.17862],  # F
                                         [-0.06835, 0.27211],
                                         [0.51685, -0.09517]]),
                               self.sample_ids,
                               self.pc_ids)
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
            Y = pd.DataFrame(np.random.randn(n, p))
            X = pd.DataFrame(np.random.randn(n_, m))
            yield npt.assert_raises, ValueError, rda, Y, X, None, None


class TestRDAResults(object):
    # STATUS: L&L only shows results with scaling 1, and they agree
    # with vegan's (module multiplying by a constant). I can also
    # compute scaling 2, agreeing with vegan, but there are no written
    # results in L&L.
    def setup(self):
        """Data from table 11.3 in Legendre & Legendre 1998."""
        self.sample_ids = ['Site0', 'Site1', 'Site2', 'Site3', 'Site4',
                           'Site5', 'Site6', 'Site7', 'Site8', 'Site9']
        self.feature_ids = ['Species0', 'Species1', 'Species2', 'Species3',
                            'Species4', 'Species5']
        self.env_ids = map(str, range(4))
        self.pc_ids = ['RDA1', 'RDA2', 'RDA3', 'RDA4', 'RDA5', 'RDA6', 'RDA7']

        self.Y = pd.DataFrame(
            np.loadtxt(get_data_path('example2_Y')),
            index=self.sample_ids, columns=self.feature_ids)

        self.X = pd.DataFrame(
            np.loadtxt(get_data_path('example2_X')),
            index=self.sample_ids, columns=self.env_ids)

    def test_scaling1(self):

        scores = rda(self.Y, self.X, scaling=1)

        biplot_scores = pd.DataFrame(np.loadtxt(
            get_data_path('example2_biplot_scaling1')))

        sample_constraints = pd.DataFrame(np.loadtxt(
            get_data_path('example2_sample_constraints_scaling1')))

        # Load data as computed with vegan 2.0-8
        vegan_features = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_species_scaling1_from_vegan')),
            index=self.feature_ids,
            columns=self.pc_ids)

        vegan_samples = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_site_scaling1_from_vegan')),
            index=self.sample_ids,
            columns=self.pc_ids)

        sample_constraints = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_sample_constraints_scaling1')),
            index=self.sample_ids,
            columns=self.pc_ids)

        biplot_scores = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_biplot_scaling1')))

        # These are wrong. @JorgeC don't forget to fix these
        proportion_explained = pd.Series([0.44275783, 0.25614586,
                                          0.15280354, 0.10497021,
                                          0.02873375, 0.00987052,
                                          0.00471828],
                                         index=self.pc_ids)

        eigvals = pd.Series([25.897954, 14.982578, 8.937841, 6.139956,
                             1.680705, 0.577350, 0.275984],
                            index=self.pc_ids)

        exp = OrdinationResults(
            'RDA', 'Redundancy Analysis',
            samples=vegan_samples,
            features=vegan_features,
            sample_constraints=sample_constraints,
            biplot_scores=biplot_scores,
            proportion_explained=proportion_explained,
            eigvals=eigvals)

        assert_ordination_results_equal(scores, exp,
                                        ignore_biplot_scores_labels=True,
                                        decimal=6)

    def test_scaling2(self):

        scores = rda(self.Y, self.X, scaling=2)

        biplot_scores = pd.DataFrame(np.loadtxt(
            get_data_path('example2_biplot_scaling2')))

        sample_constraints = pd.DataFrame(np.loadtxt(
            get_data_path('example2_sample_constraints_scaling2')))

        # Load data as computed with vegan 2.0-8
        vegan_features = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_species_scaling2_from_vegan')),
            index=self.feature_ids,
            columns=self.pc_ids)

        vegan_samples = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_site_scaling2_from_vegan')),
            index=self.sample_ids,
            columns=self.pc_ids)

        sample_constraints = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_sample_constraints_scaling2')),
            index=self.sample_ids,
            columns=self.pc_ids)

        biplot_scores = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example2_biplot_scaling2')))

        # These are wrong. @JorgeC don't forget to fix these
        proportion_explained = pd.Series([0.44275783, 0.25614586,
                                          0.15280354, 0.10497021,
                                          0.02873375, 0.00987052,
                                          0.00471828],
                                         index=self.pc_ids)

        eigvals = pd.Series([25.897954, 14.982578, 8.937841, 6.139956,
                             1.680705, 0.577350, 0.275984],
                            index=self.pc_ids)

        exp = OrdinationResults(
            'RDA', 'Redundancy Analysis',
            samples=vegan_samples,
            features=vegan_features,
            sample_constraints=sample_constraints,
            biplot_scores=biplot_scores,
            proportion_explained=proportion_explained,
            eigvals=eigvals)

        assert_ordination_results_equal(scores, exp,
                                        ignore_biplot_scores_labels=True,
                                        decimal=6)


class TestCCAErrors(object):
    def setup(self):
        """Data from table 11.3 in Legendre & Legendre 1998."""
        self.Y = pd.DataFrame(np.loadtxt(get_data_path('example3_Y')))
        self.X = pd.DataFrame(np.loadtxt(get_data_path('example3_X')))

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


class TestCCAResults1(object):
    def setup(self):
        """Data from table 11.3 in Legendre & Legendre 1998
        (p. 590). Loaded results as computed with vegan 2.0-8 and
        compared with table 11.5 if also there."""
        self.feature_ids = ['Feature0', 'Feature1', 'Feature2', 'Feature3',
                            'Feature4', 'Feature5', 'Feature6', 'Feature7',
                            'Feature8']
        self.sample_ids = ['Sample0', 'Sample1', 'Sample2', 'Sample3',
                           'Sample4', 'Sample5', 'Sample6', 'Sample7',
                           'Sample8', 'Sample9']
        self.env_ids = ['Constraint0', 'Constraint1',
                        'Constraint2', 'Constraint3']
        self.pc_ids = ['CCA1', 'CCA2', 'CCA3', 'CCA4', 'CCA5', 'CCA6', 'CCA7',
                       'CCA8', 'CCA9']
        self.Y = pd.DataFrame(
            np.loadtxt(get_data_path('example3_Y')),
            columns=self.feature_ids,
            index=self.sample_ids)
        self.X = pd.DataFrame(
            np.loadtxt(get_data_path('example3_X')),
            columns=self.env_ids,
            index=self.sample_ids
            ).ix[:, :-1]

    def test_scaling1(self):
        scores = cca(self.Y, self.X, scaling=1)

        # Load data as computed with vegan 2.0-8
        vegan_features = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_species_scaling1_from_vegan')),
            index=self.feature_ids,
            columns=self.pc_ids)

        vegan_samples = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_site_scaling1_from_vegan')),
            index=self.sample_ids,
            columns=self.pc_ids)

        sample_constraints = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_sample_constraints_scaling1')),
            index=self.sample_ids,
            columns=self.pc_ids)

        biplot_scores = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_biplot_scaling1')))

        proportion_explained = pd.Series([0.466911, 0.238327, 0.100548,
                                          0.104937, 0.044805, 0.029747,
                                          0.012631, 0.001562, 0.000532],
                                         index=self.pc_ids)
        eigvals = pd.Series([0.366136, 0.186888, 0.078847, 0.082288,
                             0.035135, 0.023327, 0.009905, 0.001225,
                             0.000417], index=self.pc_ids)

        exp = OrdinationResults(
            'CCA', 'Canonical Correspondence Analysis',
            samples=vegan_samples,
            features=vegan_features,
            sample_constraints=sample_constraints,
            biplot_scores=biplot_scores,
            proportion_explained=proportion_explained,
            eigvals=eigvals)

        assert_ordination_results_equal(scores, exp,
                                        ignore_biplot_scores_labels=True,
                                        decimal=6)

    def test_scaling2(self):
        scores = cca(self.Y, self.X, scaling=2)

        # Load data as computed with vegan 2.0-8
        vegan_features = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_species_scaling2_from_vegan')),
            index=self.feature_ids,
            columns=self.pc_ids)

        vegan_samples = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_site_scaling2_from_vegan')),
            index=self.sample_ids,
            columns=self.pc_ids)

        sample_constraints = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_sample_constraints_scaling2')),
            index=self.sample_ids,
            columns=self.pc_ids)

        biplot_scores = pd.DataFrame(
            np.loadtxt(get_data_path(
                'example3_biplot_scaling2')))

        proportion_explained = pd.Series([0.466911, 0.238327, 0.100548,
                                          0.104937, 0.044805, 0.029747,
                                          0.012631, 0.001562, 0.000532],
                                         index=self.pc_ids)
        eigvals = pd.Series([0.366136, 0.186888, 0.078847, 0.082288,
                             0.035135, 0.023327, 0.009905, 0.001225,
                             0.000417], index=self.pc_ids)

        exp = OrdinationResults(
            'CCA', 'Canonical Correspondence Analysis',
            samples=vegan_samples,
            features=vegan_features,
            sample_constraints=sample_constraints,
            biplot_scores=biplot_scores,
            proportion_explained=proportion_explained,
            eigvals=eigvals)

        assert_ordination_results_equal(scores, exp,
                                        ignore_biplot_scores_labels=True,
                                        decimal=6)


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
