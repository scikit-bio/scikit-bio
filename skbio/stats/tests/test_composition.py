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
from skbio.stats.composition import (closure, multiplicative_replacement,
                                     perturb, perturb_inv, power, inner,
                                     clr, clr_inv, ilr, ilr_inv,
                                     centralize)


class CompositionTests(TestCase):

    def setUp(self):
        # Compositional data
        self.cdata1 = np.array([[2, 2, 6],
                                [4, 4, 2]])
        self.cdata2 = np.array([2, 2, 6])

        self.cdata3 = np.array([[1, 2, 3, 0, 5],
                                [1, 0, 0, 4, 5],
                                [1, 2, 3, 4, 5]])
        self.cdata4 = np.array([1, 2, 3, 0, 5])
        self.cdata5 = [[2, 2, 6], [4, 4, 2]]
        self.cdata6 = [[1, 2, 3, 0, 5],
                       [1, 0, 0, 4, 5],
                       [1, 2, 3, 4, 5]]
        self.cdata7 = [np.exp(1), 1, 1]
        self.cdata8 = [np.exp(1), 1, 1, 1]

        # Simplicial orthonormal basis obtained from Gram-Schmidt
        self.ortho1 = [[0.44858053, 0.10905743, 0.22118102, 0.22118102],
                       [0.3379924, 0.3379924, 0.0993132, 0.22470201],
                       [0.3016453, 0.3016453, 0.3016453, 0.09506409]]

        # Real data
        self.rdata1 = [[0.70710678, -0.70710678, 0., 0.],
                       [0.40824829, 0.40824829, -0.81649658, 0.],
                       [0.28867513, 0.28867513, 0.28867513, -0.8660254]]

        # Bad datasets
        self.bad1 = np.array([1, 2, -1])
        self.bad2 = np.array([[[1, 2, 3, 0, 5]]])

    def test_closure(self):

        npt.assert_allclose(closure(self.cdata1),
                            np.array([[.2, .2, .6],
                                      [.4, .4, .2]]))
        npt.assert_allclose(closure(self.cdata2),
                            np.array([.2, .2, .6]))
        npt.assert_allclose(closure(self.cdata5),
                            np.array([[.2, .2, .6],
                                      [.4, .4, .2]]))
        with self.assertRaises(ValueError):
            closure(self.bad1)

        with self.assertRaises(ValueError):
            closure(self.bad2)

        # make sure that inplace modification is not occurring
        closure(self.cdata2)
        npt.assert_allclose(self.cdata2, np.array([2, 2, 6]))

    def test_perturb(self):
        pmat = perturb(closure(self.cdata1),
                       closure(np.array([1, 1, 1])))
        npt.assert_allclose(pmat,
                            np.array([[.2, .2, .6],
                                      [.4, .4, .2]]))

        pmat = perturb(closure(self.cdata1),
                       closure(np.array([10, 10, 20])))
        npt.assert_allclose(pmat,
                            np.array([[.125, .125, .75],
                                      [1./3, 1./3, 1./3]]))

        pmat = perturb(closure(self.cdata1),
                       closure(np.array([10, 10, 20])))
        npt.assert_allclose(pmat,
                            np.array([[.125, .125, .75],
                                      [1./3, 1./3, 1./3]]))

        pmat = perturb(closure(self.cdata2),
                       closure([1, 2, 1]))
        npt.assert_allclose(pmat, np.array([1./6, 2./6, 3./6]))

        pmat = perturb(closure(self.cdata5),
                       closure(np.array([1, 1, 1])))
        npt.assert_allclose(pmat,
                            np.array([[.2, .2, .6],
                                      [.4, .4, .2]]))

        with self.assertRaises(ValueError):
            perturb(closure(self.cdata5), self.bad1)

        # make sure that inplace modification is not occurring
        perturb(self.cdata2, [1, 2, 3])
        npt.assert_allclose(self.cdata2, np.array([2, 2, 6]))

    def test_power(self):
        pmat = power(closure(self.cdata1), 2)
        npt.assert_allclose(pmat,
                            np.array([[.04/.44, .04/.44, .36/.44],
                                      [.16/.36, .16/.36, .04/.36]]))

        pmat = power(closure(self.cdata2), 2)
        npt.assert_allclose(pmat, np.array([.04, .04, .36])/.44)

        pmat = power(closure(self.cdata5), 2)
        npt.assert_allclose(pmat,
                            np.array([[.04/.44, .04/.44, .36/.44],
                                      [.16/.36, .16/.36, .04/.36]]))

        with self.assertRaises(ValueError):
            power(self.bad1, 2)

        # make sure that inplace modification is not occurring
        power(self.cdata2, 4)
        npt.assert_allclose(self.cdata2, np.array([2, 2, 6]))

    def test_perturb_inv(self):
        pmat = perturb_inv(closure(self.cdata1),
                           closure([.1, .1, .1]))
        imat = perturb(closure(self.cdata1),
                       closure([10, 10, 10]))
        npt.assert_allclose(pmat, imat)
        pmat = perturb_inv(closure(self.cdata1),
                           closure([1, 1, 1]))
        npt.assert_allclose(pmat,
                            closure([[.2, .2, .6],
                                     [.4, .4, .2]]))
        pmat = perturb_inv(closure(self.cdata5),
                           closure([.1, .1, .1]))
        imat = perturb(closure(self.cdata1), closure([10, 10, 10]))
        npt.assert_allclose(pmat, imat)

        with self.assertRaises(ValueError):
            perturb_inv(closure(self.cdata1), self.bad1)

        # make sure that inplace modification is not occurring
        perturb_inv(self.cdata2, [1, 2, 3])
        npt.assert_allclose(self.cdata2, np.array([2, 2, 6]))

    def test_inner(self):
        a = inner(self.cdata5, self.cdata5)
        npt.assert_allclose(a, np.array([[0.80463264, -0.50766667],
                                         [-0.50766667, 0.32030201]]))

        b = inner(self.cdata7, self.cdata7)
        npt.assert_allclose(b, 0.66666666666666663)

        # Make sure that orthogonality holds
        npt.assert_allclose(inner(self.ortho1, self.ortho1), np.identity(3),
                            rtol=1e-04, atol=1e-06)

        with self.assertRaises(ValueError):
            inner(self.cdata1, self.cdata8)

        # make sure that inplace modification is not occurring
        inner(self.cdata1, self.cdata1)
        npt.assert_allclose(self.cdata1,
                            np.array([[2, 2, 6],
                                      [4, 4, 2]]))

    def test_multiplicative_replacement(self):
        amat = multiplicative_replacement(closure(self.cdata3))
        npt.assert_allclose(amat,
                            np.array([[0.087273, 0.174545, 0.261818,
                                       0.04, 0.436364],
                                      [0.092, 0.04, 0.04, 0.368, 0.46],
                                      [0.066667, 0.133333, 0.2,
                                       0.266667, 0.333333]]),
                            rtol=1e-5, atol=1e-5)

        amat = multiplicative_replacement(closure(self.cdata4))
        npt.assert_allclose(amat,
                            np.array([0.087273, 0.174545, 0.261818,
                                      0.04, 0.436364]),
                            rtol=1e-5, atol=1e-5)

        amat = multiplicative_replacement(closure(self.cdata6))
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

        # make sure that inplace modification is not occurring
        multiplicative_replacement(self.cdata4)
        npt.assert_allclose(self.cdata4, np.array([1, 2, 3, 0, 5]))

    def test_clr(self):
        cmat = clr(closure(self.cdata1))
        A = np.array([.2, .2, .6])
        B = np.array([.4, .4, .2])

        npt.assert_allclose(cmat,
                            [np.log(A / np.exp(np.log(A).mean())),
                             np.log(B / np.exp(np.log(B).mean()))])
        cmat = clr(closure(self.cdata2))
        A = np.array([.2, .2, .6])
        npt.assert_allclose(cmat,
                            np.log(A / np.exp(np.log(A).mean())))

        cmat = clr(closure(self.cdata5))
        A = np.array([.2, .2, .6])
        B = np.array([.4, .4, .2])

        npt.assert_allclose(cmat,
                            [np.log(A / np.exp(np.log(A).mean())),
                             np.log(B / np.exp(np.log(B).mean()))])
        with self.assertRaises(ValueError):
            clr(self.bad1)
        with self.assertRaises(ValueError):
            clr(self.bad2)

        # make sure that inplace modification is not occurring
        clr(self.cdata2)
        npt.assert_allclose(self.cdata2, np.array([2, 2, 6]))

    def test_clr_inv(self):
        npt.assert_allclose(clr_inv(self.rdata1), self.ortho1)

        npt.assert_allclose(clr(clr_inv(self.rdata1)), self.rdata1)

        # make sure that inplace modification is not occurring
        clr_inv(self.rdata1)
        npt.assert_allclose(self.rdata1,
                            np.array([[0.70710678, -0.70710678, 0., 0.],
                                      [0.40824829, 0.40824829,
                                       -0.81649658, 0.],
                                      [0.28867513, 0.28867513,
                                       0.28867513, -0.8660254]]))

    def test_centralize(self):
        cmat = centralize(closure(self.cdata1))
        npt.assert_allclose(cmat,
                            np.array([[0.22474487, 0.22474487, 0.55051026],
                                      [0.41523958, 0.41523958, 0.16952085]]))
        cmat = centralize(closure(self.cdata5))
        npt.assert_allclose(cmat,
                            np.array([[0.22474487, 0.22474487, 0.55051026],
                                      [0.41523958, 0.41523958, 0.16952085]]))

        with self.assertRaises(ValueError):
            centralize(self.bad1)
        with self.assertRaises(ValueError):
            centralize(self.bad2)

        # make sure that inplace modification is not occurring
        centralize(self.cdata1)
        npt.assert_allclose(self.cdata1,
                            np.array([[2, 2, 6],
                                      [4, 4, 2]]))

    def test_ilr(self):
        mat = closure(self.cdata7)
        npt.assert_array_almost_equal(ilr(mat),
                                      np.array([0.70710678, 0.40824829]))

        # Should give same result as inner
        npt.assert_allclose(ilr(self.ortho1), np.identity(3),
                            rtol=1e-04, atol=1e-06)

        with self.assertRaises(ValueError):
            ilr(self.cdata1, basis=self.cdata1)

        # make sure that inplace modification is not occurring
        ilr(self.cdata1)
        npt.assert_allclose(self.cdata1,
                            np.array([[2, 2, 6],
                                      [4, 4, 2]]))

    def test_ilr_inv(self):
        mat = closure(self.cdata7)
        npt.assert_array_almost_equal(ilr_inv(ilr(mat)), mat)

        npt.assert_allclose(ilr_inv(np.identity(3)), self.ortho1,
                            rtol=1e-04, atol=1e-06)

        with self.assertRaises(ValueError):
            ilr_inv(self.cdata1, basis=self.cdata1)

        # make sure that inplace modification is not occurring
        ilr_inv(self.cdata1)
        npt.assert_allclose(self.cdata1,
                            np.array([[2, 2, 6],
                                      [4, 4, 2]]))

if __name__ == "__main__":
    main()
