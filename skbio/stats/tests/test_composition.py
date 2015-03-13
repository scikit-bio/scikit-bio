# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main
import numpy as np
import numpy.testing as npt
from skbio.stats.composition import (_closure, multiplicative_replacement,
                                     perturb, perturb_inv, power,
                                     clr, centralize)


class CompositionTests(TestCase):

    def setUp(self):
        self.data1 = np.array([[2, 2, 6],
                               [4, 4, 2]])
        self.data2 = np.array([2, 2, 6])

        self.data3 = np.array([[1, 2, 3, 0, 5],
                               [1, 0, 0, 4, 5],
                               [1, 2, 3, 4, 5]])
        self.data4 = np.array([1, 2, 3, 0, 5])
        self.data5 = [[2, 2, 6], [4, 4, 2]]
        self.data6 = [[1, 2, 3, 0, 5],
                      [1, 0, 0, 4, 5],
                      [1, 2, 3, 4, 5]]
        # Bad datasets
        self.bad1 = np.array([1, 2, -1])
        self.bad2 = np.array([[[1, 2, 3, 0, 5]]])

    def test_closure(self):
        npt.assert_allclose(_closure(self.data1),
                            np.array([[.2, .2, .6],
                                      [.4, .4, .2]]))
        npt.assert_allclose(_closure(self.data2),
                            np.array([.2, .2, .6]))
        npt.assert_allclose(_closure(self.data5),
                            np.array([[.2, .2, .6],
                                      [.4, .4, .2]]))
        with self.assertRaises(ValueError):
            _closure(self.bad1)

        with self.assertRaises(ValueError):
            _closure(self.bad2)

    def test_perturb(self):
        pmat = perturb(_closure(self.data1),
                       _closure(np.array([1, 1, 1])))
        npt.assert_allclose(pmat,
                            np.array([[.2, .2, .6],
                                      [.4, .4, .2]]))

        pmat = perturb(_closure(self.data1),
                       _closure(np.array([10, 10, 20])))
        npt.assert_allclose(pmat,
                            np.array([[.125, .125, .75],
                                      [1./3, 1./3, 1./3]]))

        pmat = perturb(_closure(self.data1),
                       _closure(np.array([10, 10, 20])))
        npt.assert_allclose(pmat,
                            np.array([[.125, .125, .75],
                                      [1./3, 1./3, 1./3]]))

        pmat = perturb(_closure(self.data2),
                       _closure([1, 2, 1]))
        npt.assert_allclose(pmat, np.array([1./6, 2./6, 3./6]))

        pmat = perturb(_closure(self.data5),
                       _closure(np.array([1, 1, 1])))
        npt.assert_allclose(pmat,
                            np.array([[.2, .2, .6],
                                      [.4, .4, .2]]))

        with self.assertRaises(ValueError):
            perturb(_closure(self.data5), self.bad1)

    def test_power(self):
        pmat = power(_closure(self.data1), 2)
        npt.assert_allclose(pmat,
                            np.array([[.04/.44, .04/.44, .36/.44],
                                      [.16/.36, .16/.36, .04/.36]]))

        pmat = power(_closure(self.data2), 2)
        npt.assert_allclose(pmat, np.array([.04, .04, .36])/.44)

        pmat = power(_closure(self.data5), 2)
        npt.assert_allclose(pmat,
                            np.array([[.04/.44, .04/.44, .36/.44],
                                      [.16/.36, .16/.36, .04/.36]]))

        with self.assertRaises(ValueError):
            power(self.data1, 2)
        with self.assertRaises(ValueError):
            power(self.bad1, 2)

    def test_perturb_inv(self):
        pmat = perturb_inv(_closure(self.data1),
                           _closure([.1, .1, .1]))
        imat = perturb(_closure(self.data1),
                       _closure([10, 10, 10]))
        npt.assert_allclose(pmat, imat)
        pmat = perturb_inv(_closure(self.data1),
                           _closure([1, 1, 1]))
        npt.assert_allclose(pmat,
                            _closure([[.2, .2, .6],
                                      [.4, .4, .2]]))
        pmat = perturb_inv(_closure(self.data5),
                           _closure([.1, .1, .1]))
        imat = perturb(_closure(self.data1), _closure([10, 10, 10]))
        npt.assert_allclose(pmat, imat)

        with self.assertRaises(ValueError):
            perturb_inv(_closure(self.data1), self.bad1)

    def test_multiplicative_replacement(self):
        np.set_printoptions(precision=10)
        amat = multiplicative_replacement(_closure(self.data3))
        npt.assert_allclose(amat,
                            np.array([[0.087273, 0.174545, 0.261818,
                                       0.04, 0.436364],
                                      [0.092, 0.04, 0.04, 0.368, 0.46],
                                      [0.066667, 0.133333, 0.2,
                                       0.266667, 0.333333]]),
                            rtol=1e-5, atol=1e-5)

        amat = multiplicative_replacement(_closure(self.data4))
        npt.assert_allclose(amat,
                            np.array([0.087273, 0.174545, 0.261818,
                                      0.04, 0.436364]),
                            rtol=1e-5, atol=1e-5)

        amat = multiplicative_replacement(_closure(self.data6))
        npt.assert_allclose(amat,
                            np.array([[0.087273, 0.174545, 0.261818,
                                       0.04, 0.436364],
                                      [0.092, 0.04, 0.04, 0.368, 0.46],
                                      [0.066667, 0.133333, 0.2,
                                       0.266667, 0.333333]]),
                            rtol=1e-5, atol=1e-5)

        with self.assertRaises(ValueError):
            multiplicative_replacement(self.bad1)
        with self.assertRaises(ValueError):
            multiplicative_replacement(self.bad2)
        with self.assertRaises(ValueError):
            multiplicative_replacement(self.data1)

    def test_clr(self):
        cmat = clr(_closure(self.data1))
        A = np.array([.2, .2, .6])
        B = np.array([.4, .4, .2])

        npt.assert_allclose(cmat,
                            [np.log(A / np.exp(np.log(A).mean())),
                             np.log(B / np.exp(np.log(B).mean()))])
        cmat = clr(_closure(self.data2))
        A = np.array([.2, .2, .6])
        npt.assert_allclose(cmat,
                            np.log(A / np.exp(np.log(A).mean())))

        cmat = clr(_closure(self.data5))
        A = np.array([.2, .2, .6])
        B = np.array([.4, .4, .2])

        npt.assert_allclose(cmat,
                            [np.log(A / np.exp(np.log(A).mean())),
                             np.log(B / np.exp(np.log(B).mean()))])

        with self.assertRaises(ValueError):
            clr(self.data1)
        with self.assertRaises(ValueError):
            clr(self.bad1)
        with self.assertRaises(ValueError):
            clr(self.bad2)

    def test_centralize(self):
        cmat = centralize(_closure(self.data1))
        npt.assert_allclose(cmat,
                            np.array([[0.22474487, 0.22474487, 0.55051026],
                                      [0.41523958, 0.41523958, 0.16952085]]))
        cmat = centralize(_closure(self.data5))
        npt.assert_allclose(cmat,
                            np.array([[0.22474487, 0.22474487, 0.55051026],
                                      [0.41523958, 0.41523958, 0.16952085]]))

        with self.assertRaises(ValueError):
            centralize(self.data1)
        with self.assertRaises(ValueError):
            centralize(self.bad1)
        with self.assertRaises(ValueError):
            centralize(self.bad2)

if __name__ == "__main__":
    main()
