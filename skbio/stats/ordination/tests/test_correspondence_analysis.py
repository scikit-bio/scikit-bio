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
from scipy.spatial.distance import pdist
from unittest import TestCase, main

from skbio import OrdinationResults
from skbio.stats.ordination import ca
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


class TestChiSquareDistance(TestCase):
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


class TestCAResults(TestCase):
    def setUp(self):
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

        proportion_explained = pd.Series(np.array([0.701318, 0.298682]),
                                         self.pc_ids)

        exp = OrdinationResults('CA', 'Correspondance Analysis',
                                eigvals=eigvals, features=features,
                                samples=samples,
                                proportion_explained=proportion_explained)

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
        proportion_explained = pd.Series(np.array([0.701318, 0.298682]),
                                         self.pc_ids)

        exp = OrdinationResults('CA', 'Correspondance Analysis',
                                eigvals=eigvals, features=features,
                                samples=samples,
                                proportion_explained=proportion_explained)
        scores = ca(self.contingency, 1)

        assert_ordination_results_equal(exp, scores, decimal=5,
                                        ignore_directionality=True)

    def test_maintain_chi_square_distance_scaling1(self):
        """In scaling 1, chi^2 distance among rows (samples) is equal to
        euclidean distance between them in transformed space."""
        frequencies = self.X / self.X.sum()
        chi2_distances = chi_square_distance(frequencies)
        transformed_sites = ca(self.contingency, 1).samples.values
        euclidean_distances = pdist(transformed_sites, 'euclidean')
        npt.assert_almost_equal(chi2_distances, euclidean_distances)

    def test_maintain_chi_square_distance_scaling2(self):
        """In scaling 2, chi^2 distance among columns (features) is
        equal to euclidean distance between them in transformed space."""
        frequencies = self.X / self.X.sum()
        chi2_distances = chi_square_distance(frequencies, between_rows=False)
        transformed_species = ca(self.contingency, 2).features.values
        euclidean_distances = pdist(transformed_species, 'euclidean')
        npt.assert_almost_equal(chi2_distances, euclidean_distances)


class TestCAErrors(TestCase):
    def setUp(self):
        pass

    def test_negative(self):
        X = np.array([[1, 2], [-0.1, -2]])
        with npt.assert_raises(ValueError):
            ca(pd.DataFrame(X))


if __name__ == '__main__':
    main()
