#! /usr/bin/env python

#-----------------------------------------------------------------------------
# Copyright (c) 2013, The BiPy Developers.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

from __future__ import print_function
import os
import numpy as np
import numpy.testing as npt

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
    flipped, but they're still right."""
    sign_arr1 = np.sign(arr1[0])
    sign_arr2 = np.sign(arr2[0])
    differences = sign_arr1 / sign_arr2  # 1 -> equal signs, -1 -> diff signs
    return arr1 * differences, arr2


class TestCAResults(object):
    def setup(self):
        """Data from table 9.11 in Legendre & Legendre 1998."""
        X = np.loadtxt(get_data_path('L&L_CA_data'))
        self.ordination = CA(X)

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
    # compute scaling 2, agreeing with vegan, but there no written
    # results in L&L.
    def setup(self):
        """Data from table 9.11 in Legendre & Legendre 1998."""
        Y = np.loadtxt(get_data_path('example2_Y'))
        X = np.loadtxt(get_data_path('example2_X')).reshape(-1, 4, order='F')
        self.ordination = RDA(Y, X)


class TestCCAErrors(object):
    def setup(self):
        """Data from table 9.11 in Legendre & Legendre 1998."""
        self.Y = np.loadtxt(get_data_path('example3_Y'))
        self.X = np.loadtxt(get_data_path('example3_X')).reshape(-1, 4,
                                                                 order='F')

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
    # TODO: Either hardcode some results or call vegan? Hardcoding
    # them sounds better than requiring R and vegan to run tests.
    def setup(self):
        pass


class TestPCoAResults(object):
    def setup(self):
        """Sample data set from page 111 of W.J Krzanowski. Principles
        of multivariate analysis, 2000, Oxford University Press."""
        matrix = np.loadtxt(get_data_path('PCoA_sample_data'))
        dist_matrix = SymmetricDistanceMatrix(matrix,
                                              map(str, range(matrix.shape[0])))
        self.ordination = PCoA(dist_matrix)

    def test_values(self):
        """Adapted from cogent's `test_principal_coordinate_analysis`:
        "I took the example in the book (see intro info), and did the
        principal coordinates analysis, plotted the data and it looked
        right"."""
        scores = self.ordination.scores()

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
                              0.26163576, 0.28398669, 0],
                             [0.37494056, 0.22334055, -0.20892914,
                              0.05057395, -0.18710366, 0],
                             [-0.33517593, -0.23855979, -0.3099887,
                              0.11521787, -0.05021553, 0],
                             [0.25412394, -0.4123464, 0.23343642,
                              0.06403168, -0.00482608, 0],
                             [-0.28256844, 0.18606911, 0.28875631,
                              -0.06455635, -0.21141632, 0],
                             [0.01727687, 0.012458, -0.07382761,
                              -0.42690292, 0.1695749, 0]])
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
