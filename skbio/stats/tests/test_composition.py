from __future__ import absolute_import, division, print_function

# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
from numpy.testing import assert_array_almost_equal
from skbio.stats.composition import (_closure, multiplicative_replacement,
                                     perturb, perturb_inv, power,
                                     clr, centralize)


class CompositionTests(TestCase):

    def setUp(self):
        self.data1 = np.vstack((
            np.array([2, 2, 6]),
            np.array([4, 4, 2])))
        self.data2 = np.array([2, 2, 6])

        self.data3 = np.vstack((
            np.array([1, 2, 3, 0, 5]),
            np.array([1, 0, 0, 4, 5]),
            np.array(range(1, 6))))
        self.data4 = np.array([1, 2, 3, 0, 5])
        self.data5 = np.array([[[1, 2, 3, 0, 5]]])

    def test_closure(self):
        assert_array_almost_equal(_closure(self.data1),
                                  np.vstack((np.array([.2, .2, .6]),
                                             np.array([.4, .4, .2]))))
        assert_array_almost_equal(_closure(self.data2),
                                  np.array([.2, .2, .6]))
        with self.assertRaises(ValueError):
            _closure(self.data5)

    def test_perturb(self):
        pmat = perturb(_closure(self.data1), np.array([1, 1, 1]))
        assert_array_almost_equal(pmat,
                                  np.vstack((np.array([.2, .2, .6]),
                                             np.array([.4, .4, .2]))))

        pmat = perturb(_closure(self.data1), np.array([10, 10, 20]))
        assert_array_almost_equal(pmat,
                                  np.vstack((
                                      np.array([.125, .125, .75]),
                                      np.array([1./3, 1./3, 1./3]))))

        pmat = perturb(_closure(self.data1), np.array([10, 10, 20]))
        assert_array_almost_equal(pmat,
                                  np.vstack((
                                      np.array([.125, .125, .75]),
                                      np.array([1./3, 1./3, 1./3]))))

        pmat = perturb(_closure(self.data2), np.array([1, 2, 1]))
        assert_array_almost_equal(pmat, np.array([1./6, 2./6, 3./6]))

    def test_power(self):
        pmat = power(_closure(self.data1), 2)
        assert_array_almost_equal(pmat,
                                  np.vstack((
                                      np.array([.04, .04, .36])/.44,
                                      np.array([.16, .16, .04])/.36)))

        pmat = power(_closure(self.data2), 2)
        assert_array_almost_equal(pmat, np.array([.04, .04, .36])/.44)

    def test_perturb_inv(self):
        pmat = perturb_inv(_closure(self.data1), np.array([.1, .1, .1]))
        imat = perturb(_closure(self.data1), np.array([10, 10, 10]))
        assert_array_almost_equal(pmat, imat)
        pmat = perturb_inv(_closure(self.data1), np.array([1, 1, 1]))
        assert_array_almost_equal(pmat,
                                  np.vstack((
                                      np.array([.2, .2, .6]),
                                      np.array([.4, .4, .2]))))

    def test_multiplicative_replacement(self):
        amat = multiplicative_replacement(self.data3)
        assert_array_almost_equal(amat,
                                  np.array([[0.09056604, 0.18113208,
                                             0.27169811, 0.00377358,
                                             0.45283019],
                                            [0.09913793, 0.00431034,
                                             0.00431034, 0.39655172,
                                             0.49568966],
                                            [0.06666667, 0.13333333,
                                             0.2, 0.26666667, 0.33333333]]))
        amat = multiplicative_replacement(self.data4)
        assert_array_almost_equal(amat,np.array([0.09056604, 0.18113208,
                                                 0.27169811, 0.00377358,
                                                 0.45283019]))
        with self.assertRaises(ValueError):
            multiplicative_replacement(self.data5)

    def test_clr(self):
        cmat = clr(_closure(self.data1))
        A = np.array([.2, .2, .6])
        B = np.array([.4, .4, .2])

        assert_array_almost_equal(cmat,
                                  np.vstack((
                                      np.log(A / np.exp(np.log(A).mean())),
                                      np.log(B / np.exp(np.log(B).mean())))))
        cmat = clr(_closure(self.data2))
        A = np.array([.2, .2, .6])
        assert_array_almost_equal(cmat,
                                  np.log(A / np.exp(np.log(A).mean())))
        with self.assertRaises(ValueError):
            clr(self.data5)

    def test_centralize(self):
        cmat = centralize(_closure(self.data1))
        assert_array_almost_equal(cmat,
                                  np.array([[0.22474487, 0.22474487,
                                             0.55051026],
                                            [0.41523958, 0.41523958,
                                             0.16952085]]))
        with self.assertRaises(ValueError):
            centralize(self.data2)
        with self.assertRaises(ValueError):
            centralize(self.data5)

if __name__ == "__main__":
    main()
