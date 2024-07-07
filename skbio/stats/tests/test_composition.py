# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
import copy

import numpy as np
import numpy.testing as npt
from numpy.random import normal
import pandas as pd
import pandas.testing as pdt
from scipy.sparse import coo_matrix
from scipy.stats import f_oneway

from skbio import TreeNode
from skbio.util import assert_data_frame_almost_equal
from skbio.stats.distance import DistanceMatrixError
from skbio.stats.composition import (
    closure, multi_replace, multiplicative_replacement, perturb, perturb_inv, power,
    inner, clr, clr_inv, ilr, ilr_inv, alr, alr_inv, sbp_basis, _gram_schmidt_basis,
    centralize, _calc_p_adjust, ancom, vlr, pairwise_vlr, tree_basis, dirmult_ttest)


def assert_coo_allclose(res, exp, rtol=1e-7, atol=1e-7):
    res_data = np.vstack((res.row, res.col, res.data)).T
    exp_data = np.vstack((exp.row, exp.col, exp.data)).T

    # sort by row and col
    res_data = res_data[res_data[:, 1].argsort()]
    res_data = res_data[res_data[:, 0].argsort()]
    exp_data = exp_data[exp_data[:, 1].argsort()]
    exp_data = exp_data[exp_data[:, 0].argsort()]
    npt.assert_allclose(res_data, exp_data, rtol=rtol, atol=atol)


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
        # negative count
        self.bad1 = np.array([1, 2, -1])
        # zero count
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

    def test_closure_warning(self):
        with self.assertRaises(ValueError):
            closure([0., 0., 0.])

        with self.assertRaises(ValueError):
            closure([[0., 0., 0.],
                     [0., 5., 5.]])

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

    def test_multi_replace(self):
        amat = multi_replace(closure(self.cdata3))
        npt.assert_allclose(amat,
                            np.array([[0.087273, 0.174545, 0.261818,
                                       0.04, 0.436364],
                                      [0.092, 0.04, 0.04, 0.368, 0.46],
                                      [0.066667, 0.133333, 0.2,
                                       0.266667, 0.333333]]),
                            rtol=1e-5, atol=1e-5)

        amat = multi_replace(closure(self.cdata4))
        npt.assert_allclose(amat,
                            np.array([0.087273, 0.174545, 0.261818,
                                      0.04, 0.436364]),
                            rtol=1e-5, atol=1e-5)

        amat = multi_replace(closure(self.cdata6))
        npt.assert_allclose(amat,
                            np.array([[0.087273, 0.174545, 0.261818,
                                       0.04, 0.436364],
                                      [0.092, 0.04, 0.04, 0.368, 0.46],
                                      [0.066667, 0.133333, 0.2,
                                       0.266667, 0.333333]]),
                            rtol=1e-5, atol=1e-5)

        with self.assertRaises(ValueError):
            multi_replace(self.bad1)
        with self.assertRaises(ValueError):
            multi_replace(self.bad2)

        # make sure that inplace modification is not occurring
        multi_replace(self.cdata4)
        npt.assert_allclose(self.cdata4, np.array([1, 2, 3, 0, 5]))

    def multi_replace_warning(self):
        with self.assertRaises(ValueError):
            multi_replace([0, 1, 2], delta=1)

    def test_multiplicative_replacement(self):
        mat = closure(self.cdata3)
        npt.assert_allclose(multiplicative_replacement(mat), multi_replace(mat))

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
        npt.assert_allclose(clr(clr_inv(self.rdata1)), self.rdata1,
                            rtol=1e-4, atol=1e-5)

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

        # no check
        npt.assert_array_almost_equal(ilr(mat, check=False),
                                      np.array([0.70710678, 0.40824829]))

        with self.assertRaises(ValueError):
            ilr(self.cdata1, basis=self.cdata1)

        # make sure that inplace modification is not occurring
        ilr(self.cdata1)
        npt.assert_allclose(self.cdata1,
                            np.array([[2, 2, 6],
                                      [4, 4, 2]]))

    def test_ilr_basis(self):
        table = np.array([[1., 10.],
                          [1.14141414, 9.90909091],
                          [1.28282828, 9.81818182],
                          [1.42424242, 9.72727273],
                          [1.56565657, 9.63636364]])
        basis = np.atleast_2d(clr([[0.80442968, 0.19557032]]))
        res = ilr(table, basis=basis)
        exp = np.array([[np.log(1/10)*np.sqrt(1/2)],
                        [np.log(1.14141414 / 9.90909091)*np.sqrt(1/2)],
                        [np.log(1.28282828 / 9.81818182)*np.sqrt(1/2)],
                        [np.log(1.42424242 / 9.72727273)*np.sqrt(1/2)],
                        [np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]])

        npt.assert_allclose(res, exp)

    def test_ilr_basis_one_dimension_error(self):
        table = np.array([[1., 10.],
                          [1.14141414, 9.90909091],
                          [1.28282828, 9.81818182],
                          [1.42424242, 9.72727273],
                          [1.56565657, 9.63636364]])
        basis = np.array([0.80442968, 0.19557032])
        with self.assertRaises(ValueError):
            ilr(table, basis=basis)

    def test_ilr_inv(self):
        mat = closure(self.cdata7)
        npt.assert_array_almost_equal(ilr_inv(ilr(mat)), mat)

        npt.assert_allclose(ilr_inv(np.identity(3)), self.ortho1,
                            rtol=1e-04, atol=1e-06)

        # no check
        npt.assert_array_almost_equal(ilr_inv(ilr(mat), check=False), mat)

        with self.assertRaises(ValueError):
            ilr_inv(self.cdata1, basis=self.cdata1)

        # make sure that inplace modification is not occurring
        ilr_inv(self.cdata1)
        npt.assert_allclose(self.cdata1,
                            np.array([[2, 2, 6],
                                      [4, 4, 2]]))

    def test_ilr_basis_isomorphism(self):
        # tests to make sure that the isomorphism holds
        # with the introduction of the basis.
        basis = np.atleast_2d(clr([[0.80442968, 0.19557032]]))
        table = np.array([[np.log(1/10)*np.sqrt(1/2),
                           np.log(1.14141414 / 9.90909091)*np.sqrt(1/2),
                           np.log(1.28282828 / 9.81818182)*np.sqrt(1/2),
                           np.log(1.42424242 / 9.72727273)*np.sqrt(1/2),
                           np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]]).T
        lr = ilr_inv(table, basis=basis)
        res = ilr(lr, basis=basis)
        npt.assert_allclose(res, table)

        table = np.array([[1., 10.],
                          [1.14141414, 9.90909091],
                          [1.28282828, 9.81818182],
                          [1.42424242, 9.72727273],
                          [1.56565657, 9.63636364]])

        res = ilr_inv(ilr(table, basis=basis), basis=basis)
        npt.assert_allclose(res, closure(table.squeeze()))

    def test_ilr_inv_basis(self):
        exp = closure(np.array([[1., 10.],
                                [1.14141414, 9.90909091],
                                [1.28282828, 9.81818182],
                                [1.42424242, 9.72727273],
                                [1.56565657, 9.63636364]]))
        basis = np.atleast_2d(clr([[0.80442968, 0.19557032]]))
        table = np.array([[np.log(1/10)*np.sqrt(1/2),
                           np.log(1.14141414 / 9.90909091)*np.sqrt(1/2),
                           np.log(1.28282828 / 9.81818182)*np.sqrt(1/2),
                           np.log(1.42424242 / 9.72727273)*np.sqrt(1/2),
                           np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]]).T

        res = ilr_inv(table, basis=basis)
        npt.assert_allclose(res, exp)

    def test_ilr_inv_basis_one_dimension_error(self):
        basis = clr([0.80442968, 0.19557032])
        table = np.array([[np.log(1/10)*np.sqrt(1/2),
                           np.log(1.14141414 / 9.90909091)*np.sqrt(1/2),
                           np.log(1.28282828 / 9.81818182)*np.sqrt(1/2),
                           np.log(1.42424242 / 9.72727273)*np.sqrt(1/2),
                           np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]]).T
        with self.assertRaises(ValueError):
            ilr_inv(table, basis=basis)

    def test_alr(self):
        # 2d-composition
        comp1 = closure(self.cdata1)
        alr2d_byhand = np.array([np.log(comp1[:, 0]/comp1[:, 1]),
                                 np.log(comp1[:, 2]/comp1[:, 1])]).T
        alr2d_method = alr(comp1, denominator_idx=1)
        npt.assert_allclose(alr2d_byhand, alr2d_method)

        # 1d-composition
        comp2 = closure(self.cdata2)
        alr1d_byhand = np.array([np.log(comp2[0]/comp2[1]),
                                 np.log(comp2[2]/comp2[1])]).T
        alr1d_method = alr(comp2, denominator_idx=1)
        npt.assert_allclose(alr1d_byhand, alr1d_method)

        with self.assertRaises(ValueError):
            alr(self.bad1)
        with self.assertRaises(ValueError):
            alr(self.bad2)

        # make sure that inplace modification is not occurring
        alr(self.cdata2)
        npt.assert_allclose(self.cdata2, np.array([2, 2, 6]))

        # matrix must be 1d or 2d
        with self.assertRaises(ValueError):
            alr(np.atleast_3d(self.cdata2))

    def test_alr_inv(self):
        # 2d-composition
        comp1 = closure(self.cdata1)
        alr2d_byhand = np.array([np.log(comp1[:, 0]/comp1[:, 1]),
                                 np.log(comp1[:, 2]/comp1[:, 1])]).T
        alr2d_method = alr(comp1, denominator_idx=1)
        B = 1/(1 + np.exp(alr2d_byhand[:, 0]) + np.exp(alr2d_byhand[:, 1]))
        A = B * np.exp(alr2d_byhand[:, 0])
        C = B * np.exp(alr2d_byhand[:, 1])
        alrinv2d_byhand = np.column_stack((A, B, C))
        alrinv2d_method = alr_inv(alr2d_method, denominator_idx=1)
        npt.assert_allclose(alrinv2d_byhand, alrinv2d_method)

        # 1d-composition
        comp2 = closure(self.cdata2)
        alr1d_byhand = np.array([np.log(comp2[0]/comp2[1]),
                                 np.log(comp2[2]/comp2[1])]).T
        alr1d_method = alr(comp2, denominator_idx=1)
        B = 1/(1 + np.exp(alr1d_byhand[0]) + np.exp(alr1d_byhand[1]))
        A = B * np.exp(alr1d_byhand[0])
        C = B * np.exp(alr1d_byhand[1])
        alrinv1d_byhand = np.column_stack((A, B, C))[0, :]
        alrinv1d_method = alr_inv(alr1d_method, denominator_idx=1)
        npt.assert_allclose(alrinv1d_byhand, alrinv1d_method)

        # make sure that inplace modification is not occurring
        alr_inv(self.rdata1)
        npt.assert_allclose(self.rdata1,
                            np.array([[0.70710678, -0.70710678, 0., 0.],
                                      [0.40824829, 0.40824829,
                                       -0.81649658, 0.],
                                      [0.28867513, 0.28867513,
                                       0.28867513, -0.8660254]]))

        with self.assertRaises(ValueError):
            alr_inv(self.bad2)

    def test_sbp_basis_gram_schmidt(self):
        gsbasis = _gram_schmidt_basis(5)
        sbp = np.array([[1, -1, 0, 0, 0],
                        [1, 1, -1, 0, 0],
                        [1, 1, 1, -1, 0],
                        [1, 1, 1, 1, -1]])
        sbpbasis = sbp_basis(sbp)
        npt.assert_allclose(gsbasis, sbpbasis)

    def test_sbp_basis_elementwise(self):
        sbp = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1],
                        [1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, 0],
                        [1, 1, 1, 1, 1, 1, -1, 0, 0, 0, 0, 0],
                        [1, 1, -1, -1, -1, 1, 0, 0, 0, 0, 0, 0],
                        [1, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                        [1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
                        [0, 0, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, 0],
                        [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0]])
        sbpbasis = sbp_basis(sbp)
        # by hand, element-wise
        r = np.apply_along_axis(func1d=lambda x: np.sum(x > 0),
                                axis=1, arr=sbp)
        s = np.apply_along_axis(func1d=lambda x: np.sum(x < 0),
                                axis=1, arr=sbp)
        psi = np.zeros(sbp.shape)
        for i in range(0, sbp.shape[0]):
            for j in range(0, sbp.shape[1]):
                if sbp[i, j] == 1:
                    psi[i, j] = np.sqrt(s[i]/(r[i]*(r[i]+s[i])))
                elif sbp[i, j] == -1:
                    psi[i, j] = -np.sqrt(r[i]/(s[i]*(r[i]+s[i])))
        npt.assert_allclose(psi, sbpbasis)


class TestTreeBasis(TestCase):

    def test_tree_basis_base_case(self):
        tree = u"(a,b);"
        t = TreeNode.read([tree])

        exp_basis = coo_matrix(
            np.array([[-np.sqrt(1. / 2),
                       np.sqrt(1. / 2)]]))
        exp_keys = [t.name]
        res_basis, res_keys = tree_basis(t)

        assert_coo_allclose(exp_basis, res_basis)
        self.assertListEqual(exp_keys, res_keys)

    def test_tree_basis_invalid(self):
        with self.assertRaises(ValueError):
            tree = u"(a,b,c);"
            t = TreeNode.read([tree])
            tree_basis(t)

    def test_tree_basis_unbalanced(self):
        tree = u"((a,b)c, d);"
        t = TreeNode.read([tree])
        exp_basis = coo_matrix(np.array(
            [[-np.sqrt(1. / 6), -np.sqrt(1. / 6), np.sqrt(2. / 3)],
             [-np.sqrt(1. / 2), np.sqrt(1. / 2), 0]]
        ))
        exp_keys = [t.name, t[0].name]
        res_basis, res_keys = tree_basis(t)

        assert_coo_allclose(exp_basis, res_basis)
        self.assertListEqual(exp_keys, res_keys)

    def test_tree_basis_unbalanced2(self):
        tree = u"(d, (a,b)c);"

        t = TreeNode.read([tree])

        exp_basis = coo_matrix(np.array(
            [
                [-np.sqrt(2. / 3), np.sqrt(1. / 6), np.sqrt(1. / 6)],
                [0, -np.sqrt(1. / 2), np.sqrt(1. / 2)]
            ]
        ))

        exp_keys = [t.name, t[1].name]
        res_basis, res_keys = tree_basis(t)
        assert_coo_allclose(exp_basis, res_basis, atol=1e-7, rtol=1e-7)
        self.assertListEqual(exp_keys, res_keys)


class AncomTests(TestCase):
    def setUp(self):
        # Basic count data with 2 groupings
        self.table1 = pd.DataFrame([
            [10, 10, 10, 20, 20, 20],
            [11, 12, 11, 21, 21, 21],
            [10, 11, 10, 10, 11, 10],
            [10, 11, 10, 10, 10, 9],
            [10, 11, 10, 10, 10, 10],
            [10, 11, 10, 10, 10, 11],
            [10, 13, 10, 10, 10, 12]]).T
        self.cats1 = pd.Series([0, 0, 0, 1, 1, 1])
        # Real valued data with 2 groupings
        D, L = 40, 80
        np.random.seed(0)
        self.table2 = np.vstack((np.concatenate((normal(10, 1, D),
                                                 normal(200, 1, D))),
                                 np.concatenate((normal(20, 1, D),
                                                 normal(100000, 1, D))),
                                 normal(10, 1, L),
                                 normal(10, 1, L),
                                 np.concatenate((normal(20, 1, D),
                                                 normal(100000, 1, D))),
                                 normal(10, 1, L),
                                 normal(10, 1, L),
                                 normal(10, 1, L),
                                 normal(10, 1, L)))
        self.table2 = np.absolute(self.table2)
        self.table2 = pd.DataFrame(self.table2.astype(int).T)
        self.cats2 = pd.Series([0]*D + [1]*D)

        # Real valued data with 2 groupings and no significant difference
        self.table3 = pd.DataFrame([
            [10, 10.5, 10, 10, 10.5, 10.3],
            [11, 11.5, 11, 11, 11.5, 11.3],
            [10, 10.5, 10, 10, 10.5, 10.2],
            [10, 10.5, 10, 10, 10.5, 10.3],
            [10, 10.5, 10, 10, 10.5, 10.1],
            [10, 10.5, 10, 10, 10.5, 10.6],
            [10, 10.5, 10, 10, 10.5, 10.4]]).T
        self.cats3 = pd.Series([0, 0, 0, 1, 1, 1])

        # Real valued data with 3 groupings
        D, L = 40, 120
        np.random.seed(0)
        self.table4 = np.vstack((np.concatenate((normal(10, 1, D),
                                                 normal(200, 1, D),
                                                 normal(400, 1, D))),
                                 np.concatenate((normal(20, 1, D),
                                                 normal(100000, 1, D),
                                                 normal(2000, 1, D))),
                                 normal(10, 1, L),
                                 normal(10, 1, L),
                                 np.concatenate((normal(20, 1, D),
                                                 normal(100000, 1, D),
                                                 normal(2000, 1, D))),
                                 normal(10, 1, L),
                                 normal(10, 1, L),
                                 normal(10, 1, L),
                                 normal(10, 1, L)))
        self.table4 = np.absolute(self.table4)
        self.table4 = pd.DataFrame(self.table4.astype(int).T)
        self.cats4 = pd.Series([0]*D + [1]*D + [2]*D)

        # Noncontiguous case
        self.table5 = pd.DataFrame([
            [11, 12, 21, 11, 21, 21],
            [10, 11, 10, 10, 11, 10],
            [10, 11, 10, 10, 10, 9],
            [10, 11, 10, 10, 10, 10],
            [10, 11, 10, 10, 10, 11],
            [10, 10, 20, 9,  20, 20],
            [10, 13, 10, 10, 10, 12]]).T
        self.cats5 = pd.Series([0, 0, 1, 0, 1, 1])

        # Different number of classes case
        self.table6 = pd.DataFrame([
            [11, 12, 9, 11, 21, 21],
            [10, 11, 10, 10, 11, 10],
            [10, 11, 10, 10, 10, 9],
            [10, 11, 10, 10, 10, 10],
            [10, 11, 10, 10, 10, 11],
            [10, 10, 10, 9,  20, 20],
            [10, 13, 10, 10, 10, 12]]).T
        self.cats6 = pd.Series([0, 0, 0, 0, 1, 1])

        # Categories are letters
        self.table7 = pd.DataFrame([
            [11, 12, 9, 11, 21, 21],
            [10, 11, 10, 10, 11, 10],
            [10, 11, 10, 10, 10, 9],
            [10, 11, 10, 10, 10, 10],
            [10, 11, 10, 10, 10, 11],
            [10, 10, 10, 9,  20, 20],
            [10, 13, 10, 10, 10, 12]]).T
        self.cats7 = pd.Series(['a', 'a', 'a', 'a', 'b', 'b'])

        # Swap samples
        self.table8 = pd.DataFrame([
            [10, 10, 10, 20, 20, 20],
            [11, 12, 11, 21, 21, 21],
            [10, 11, 10, 10, 11, 10],
            [10, 11, 10, 10, 10, 9],
            [10, 11, 10, 10, 10, 10],
            [10, 11, 10, 10, 10, 11],
            [10, 13, 10, 10, 10, 12]]).T
        self.table8.index = ['a', 'b', 'c',
                             'd', 'e', 'f']
        self.cats8 = pd.Series([0, 0, 1, 0, 1, 1],
                               index=['a', 'b', 'd',
                                      'c', 'e', 'f'])

        # Real valued data with 3 groupings
        D, L = 40, 120
        np.random.seed(0)
        self.table9 = np.vstack((np.concatenate((normal(10, 1, D),
                                                 normal(200, 1, D),
                                                 normal(400, 1, D))),
                                 np.concatenate((normal(200000, 1, D),
                                                 normal(10, 1, D),
                                                 normal(2000, 1, D))),
                                 normal(10, 10, L),
                                 normal(10, 10, L),
                                 np.concatenate((normal(2000, 1, D),
                                                 normal(100000, 1, D),
                                                 normal(2000, 1, D))),
                                 normal(10000, 1000, L),
                                 normal(10, 10, L),
                                 normal(10, 10, L),
                                 normal(10, 10, L),
                                 normal(10000, 1000, L),
                                 normal(10, 10, L),
                                 normal(10, 10, L),
                                 normal(10, 10, L),
                                 np.concatenate((normal(2000, 1, D),
                                                 normal(100000, 1, D),
                                                 normal(2000, 1, D))),
                                 normal(10000, 1000, L),
                                 normal(10, 10, L),
                                 normal(10, 10, L),
                                 normal(10, 10, L)))
        self.table9 = np.absolute(self.table9)+1
        self.table9 = pd.DataFrame(self.table9.astype(int).T)
        self.cats9 = pd.Series([0]*D + [1]*D + [2]*D)

        # Real valued data with 2 groupings
        D, L = 40, 80
        np.random.seed(0)
        self.table10 = np.vstack((np.concatenate((normal(10, 1, D),
                                                  normal(200, 1, D))),
                                  np.concatenate((normal(10, 1, D),
                                                  normal(200, 1, D))),
                                  np.concatenate((normal(20, 10, D),
                                                  normal(100, 10, D))),
                                  normal(10, 1, L),
                                  np.concatenate((normal(200, 100, D),
                                                  normal(100000, 100, D))),
                                  np.concatenate((normal(200000, 100, D),
                                                  normal(300, 100, D))),
                                  np.concatenate((normal(200000, 100, D),
                                                  normal(300, 100, D))),
                                  np.concatenate((normal(20, 20, D),
                                                  normal(40, 10, D))),
                                  np.concatenate((normal(20, 20, D),
                                                  normal(40, 10, D))),
                                  np.concatenate((normal(20, 20, D),
                                                  normal(40, 10, D))),
                                  normal(100, 10, L),
                                  normal(100, 10, L),
                                  normal(1000, 10, L),
                                  normal(1000, 10, L),
                                  normal(10, 10, L),
                                  normal(10, 10, L),
                                  normal(10, 10, L),
                                  normal(10, 10, L)))
        self.table10 = np.absolute(self.table10) + 1
        self.table10 = pd.DataFrame(self.table10.astype(int).T)
        self.cats10 = pd.Series([0]*D + [1]*D)

        # zero count
        self.bad1 = pd.DataFrame(np.array([
            [10, 10, 10, 20, 20, 0],
            [11, 11, 11, 21, 21, 21],
            [10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 10]]).T)
        # negative count
        self.bad2 = pd.DataFrame(np.array([
            [10, 10, 10, 20, 20, 1],
            [11, 11, 11, 21, 21, 21],
            [10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, -1],
            [10, 10, 10, 10, 10, 10]]).T)

        # missing count
        self.bad3 = pd.DataFrame(np.array([
            [10, 10, 10, 20, 20, 1],
            [11, 11, 11, 21, 21, 21],
            [10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, 10],
            [10, 10, 10, 10, 10, np.nan],
            [10, 10, 10, 10, 10, 10]]).T)
        self.badcats1 = pd.Series([0, 0, 0, 1, np.nan, 1])
        self.badcats2 = pd.Series([0, 0, 0, 0, 0, 0])
        self.badcats3 = pd.Series([0, 0, 1, 1])
        self.badcats4 = pd.Series(range(len(self.table1)))
        self.badcats5 = pd.Series([1]*len(self.table1))

    def test_ancom_basic_counts(self):
        test_table = pd.DataFrame(self.table1)
        original_table = copy.deepcopy(test_table)
        test_cats = pd.Series(self.cats1)
        original_cats = copy.deepcopy(test_cats)
        result = ancom(test_table, test_cats, p_adjust=None)
        # Test to make sure that the input table hasn't be altered
        assert_data_frame_almost_equal(original_table, test_table)
        # Test to make sure that the input table hasn't be altered
        pdt.assert_series_equal(original_cats, test_cats)
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Reject null hypothesis': np.array([True, True, False, False,
                                                 False, False, False],
                                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_percentiles(self):
        table = pd.DataFrame([[12, 11],
                              [9, 11],
                              [1, 11],
                              [22, 100],
                              [20, 53],
                              [23, 1]],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'],
                             columns=['b1', 'b2'])
        grouping = pd.Series(['a', 'a', 'a', 'b', 'b', 'b'],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'])

        percentiles = [0.0, 25.0, 50.0, 75.0, 100.0]
        groups = ['a', 'b']
        tuples = [(p, g) for g in groups for p in percentiles]
        exp_mi = pd.MultiIndex.from_tuples(tuples,
                                           names=['Percentile', 'Group'])
        exp_data = np.array(
            [[1.0, 11.0], [5.0, 11.0], [9.0, 11.0], [10.5, 11.0], [12.0, 11.0],
             [20.0, 1.0], [21.0, 27.0], [22.0, 53.0], [22.5, 76.5],
             [23.0, 100.0]])
        exp = pd.DataFrame(exp_data.T, columns=exp_mi, index=['b1', 'b2'])

        result = ancom(table, grouping)[1]
        assert_data_frame_almost_equal(result, exp)

    def test_ancom_percentiles_alt_categories(self):
        table = pd.DataFrame([[12],
                              [9],
                              [1],
                              [22],
                              [20],
                              [23]],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'],
                             columns=['b1'])
        grouping = pd.Series(['a', 'a', 'c', 'b', 'b', 'c'],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'])

        percentiles = [0.0, 25.0, 50.0, 75.0, 100.0]
        groups = ['a', 'b', 'c']
        tuples = [(p, g) for g in groups for p in percentiles]
        exp_mi = pd.MultiIndex.from_tuples(tuples,
                                           names=['Percentile', 'Group'])
        exp_data = np.array([[9.0], [9.75], [10.5], [11.25], [12.0],  # a
                             [20.0], [20.5], [21.0], [21.5], [22.0],  # b
                             [1.0], [6.5], [12.0], [17.5], [23.0]])   # c
        exp = pd.DataFrame(exp_data.T, columns=exp_mi, index=['b1'])

        result = ancom(table, grouping, percentiles=percentiles)[1]
        assert_data_frame_almost_equal(result, exp)

    def test_ancom_alt_percentiles(self):
        table = pd.DataFrame([[12],
                              [9],
                              [1],
                              [22],
                              [20],
                              [23]],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'],
                             columns=['b1'])
        grouping = pd.Series(['a', 'a', 'a', 'b', 'b', 'b'],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'])

        percentiles = [42.0, 50.0]
        groups = ['a', 'b']
        tuples = [(p, g) for g in groups for p in percentiles]
        exp_mi = pd.MultiIndex.from_tuples(tuples,
                                           names=['Percentile', 'Group'])
        exp_data = np.array([[7.71999999], [9.0],  # a
                             [21.68], [22.0]])     # b
        exp = pd.DataFrame(exp_data.T, columns=exp_mi, index=['b1'])

        result = ancom(table, grouping, percentiles=percentiles)[1]
        assert_data_frame_almost_equal(result, exp)

    def test_ancom_percentiles_swapped(self):
        table = pd.DataFrame([[12],
                              [9],
                              [1],
                              [22],
                              [20],
                              [23]],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'],
                             columns=['b1'])
        grouping = pd.Series(['a', 'a', 'b', 'a', 'b', 'b'],
                             index=['s1', 's2', 's4', 's3', 's5', 's6'])

        percentiles = [42.0, 50.0]
        groups = ['a', 'b']
        tuples = [(p, g) for g in groups for p in percentiles]
        exp_mi = pd.MultiIndex.from_tuples(tuples,
                                           names=['Percentile', 'Group'])
        exp_data = np.array([[7.71999999], [9.0],  # a
                             [21.68], [22.0]])     # b
        exp = pd.DataFrame(exp_data.T, columns=exp_mi, index=['b1'])

        result = ancom(table, grouping, percentiles=percentiles)[1]
        assert_data_frame_almost_equal(result, exp)

    def test_ancom_percentile_order_unimportant(self):
        table = pd.DataFrame([[12],
                              [9],
                              [1],
                              [22],
                              [20],
                              [23]],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'],
                             columns=['b1'])
        grouping = pd.Series(['a', 'a', 'a', 'b', 'b', 'b'],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'])
        # order of percentiles in unimportant after sorting
        result1 = ancom(table, grouping, percentiles=[50.0, 42.0])[1]
        result2 = ancom(table, grouping, percentiles=[42.0, 50.0])[1]
        assert_data_frame_almost_equal(
            result1.sort_index(axis=1), result2.sort_index(axis=1))

    def test_ancom_percentiles_iterator(self):
        table = pd.DataFrame([[12],
                              [9],
                              [1],
                              [22],
                              [20],
                              [23]],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'],
                             columns=['b1'])
        grouping = pd.Series(['a', 'a', 'a', 'b', 'b', 'b'],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'])

        percentiles = [42.0, 50.0]
        groups = ['a', 'b']
        tuples = [(p, g) for g in groups for p in percentiles]
        exp_mi = pd.MultiIndex.from_tuples(tuples,
                                           names=['Percentile', 'Group'])
        exp_data = np.array([[7.71999999], [9.0],  # a
                             [21.68], [22.0]])     # b
        exp = pd.DataFrame(exp_data.T, columns=exp_mi, index=['b1'])

        result = ancom(table, grouping, percentiles=iter(percentiles))[1]
        assert_data_frame_almost_equal(result, exp)

    def test_ancom_no_percentiles(self):
        table = pd.DataFrame([[12],
                              [9],
                              [1],
                              [22],
                              [20],
                              [23]],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'],
                             columns=['b1'])
        grouping = pd.Series(['a', 'a', 'a', 'b', 'b', 'b'],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'])
        result = ancom(table, grouping, percentiles=[])[1]
        assert_data_frame_almost_equal(result, pd.DataFrame())

    def test_ancom_percentile_out_of_range(self):
        table = pd.DataFrame([[12],
                              [9],
                              [1],
                              [22],
                              [20],
                              [23]],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'],
                             columns=['b1'])
        grouping = pd.Series(['a', 'a', 'a', 'b', 'b', 'b'],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'])
        with self.assertRaises(ValueError):
            ancom(table, grouping, percentiles=[-1.0])
        with self.assertRaises(ValueError):
            ancom(table, grouping, percentiles=[100.1])
        with self.assertRaises(ValueError):
            ancom(table, grouping, percentiles=[10.0, 3.0, 101.0, 100])

    def test_ancom_duplicate_percentiles(self):
        table = pd.DataFrame([[12],
                              [9],
                              [1],
                              [22],
                              [20],
                              [23]],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'],
                             columns=['b1'])
        grouping = pd.Series(['a', 'a', 'a', 'b', 'b', 'b'],
                             index=['s1', 's2', 's3', 's4', 's5', 's6'])
        with self.assertRaises(ValueError):
            ancom(table, grouping, percentiles=[10.0, 10.0])

    def test_ancom_basic_proportions(self):
        # Converts from counts to proportions
        test_table = pd.DataFrame(closure(self.table1))
        original_table = copy.deepcopy(test_table)
        test_cats = pd.Series(self.cats1)
        original_cats = copy.deepcopy(test_cats)
        result = ancom(test_table, test_cats, p_adjust=None)
        # Test to make sure that the input table hasn't be altered
        assert_data_frame_almost_equal(original_table, test_table)
        # Test to make sure that the input table hasn't be altered
        pdt.assert_series_equal(original_cats, test_cats)
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Reject null hypothesis': np.array([True, True, False, False,
                                                 False, False, False],
                                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_multiple_groups(self):
        test_table = pd.DataFrame(self.table4)
        original_table = copy.deepcopy(test_table)
        test_cats = pd.Series(self.cats4)
        original_cats = copy.deepcopy(test_cats)
        result = ancom(test_table, test_cats)
        # Test to make sure that the input table hasn't be altered
        assert_data_frame_almost_equal(original_table, test_table)
        # Test to make sure that the input table hasn't be altered
        pdt.assert_series_equal(original_cats, test_cats)
        exp = pd.DataFrame(
            {'W': np.array([8, 7, 3, 3, 7, 3, 3, 3, 3]),
             'Reject null hypothesis': np.array([True, True, False, False,
                                                 True, False, False, False,
                                                 False], dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_noncontiguous(self):
        result = ancom(self.table5, self.cats5, p_adjust=None)
        exp = pd.DataFrame(
            {'W': np.array([6, 2, 2, 2, 2, 6, 2]),
             'Reject null hypothesis': np.array([True, False, False, False,
                                                 False, True, False],
                                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_unbalanced(self):
        result = ancom(self.table6, self.cats6, p_adjust=None)
        exp = pd.DataFrame(
            {'W': np.array([5, 3, 3, 2, 2, 5, 2]),
             'Reject null hypothesis': np.array([True, False, False, False,
                                                 False, True, False],
                                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_letter_categories(self):
        result = ancom(self.table7, self.cats7, p_adjust=None)
        exp = pd.DataFrame(
            {'W': np.array([5, 3, 3, 2, 2, 5, 2]),
             'Reject null hypothesis': np.array([True, False, False, False,
                                                 False, True, False],
                                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_significance_test_none(self):
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Reject null hypothesis': np.array([True, True, False, False,
                                                 False, False, False],
                                                dtype=bool)})
        result = ancom(self.table1, self.cats1, significance_test=None)
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_significance_test_callable(self):
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Reject null hypothesis': np.array([True, True, False, False,
                                                 False, False, False],
                                                dtype=bool)})
        result = ancom(self.table1, self.cats1, significance_test=f_oneway)
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_multiple_comparisons(self):
        exp = pd.DataFrame(
            {'W': np.array([0] * 7),
             'Reject null hypothesis': np.array([False] * 7, dtype=bool)})
        for method in 'holm', 'bh':
            result = ancom(self.table1, self.cats1, p_adjust=method,
                           significance_test='mannwhitneyu')
            assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_multiple_comparisons_deprecated(self):
        # @deprecated
        exp = pd.DataFrame(
            {'W': np.array([0] * 7),
             'Reject null hypothesis': np.array([False] * 7, dtype=bool)})
        result = ancom(self.table1, self.cats1,
                       significance_test='mannwhitneyu',
                       multiple_comparisons_correction=None)
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_alternative_test(self):
        result = ancom(self.table1, self.cats1, p_adjust=None,
                       significance_test="ttest_ind")
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Reject null hypothesis': np.array([True,  True, False, False,
                                                 False, False, False],
                                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_incorrect_test(self):
        with self.assertRaises(ValueError) as cm:
            ancom(self.table1, self.cats1, significance_test="not_a_test")
        msg = 'Function "not_a_test" does not exist under scipy.stats.'
        self.assertEqual(str(cm.exception), msg)

    def test_ancom_normal_data(self):
        result = ancom(self.table2, self.cats2, p_adjust=None,
                       significance_test="ttest_ind")
        exp = pd.DataFrame(
            {'W': np.array([8, 8, 3, 3, 8, 3, 3, 3, 3]),
             'Reject null hypothesis': np.array([True, True, False, False,
                                                 True, False, False,
                                                 False, False],
                                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_basic_counts_swapped(self):
        result = ancom(self.table8, self.cats8)
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Reject null hypothesis': np.array([True, True, False, False,
                                                 False, False, False],
                                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_no_signal(self):
        result = ancom(self.table3, self.cats3, p_adjust=None)
        exp = pd.DataFrame(
            {'W': np.array([0]*7),
             'Reject null hypothesis': np.array([False]*7, dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_tau(self):
        exp1 = pd.DataFrame(
            {'W': np.array([8, 7, 3, 3, 7, 3, 3, 3, 3]),
             'Reject null hypothesis': np.array([True, False, False, False,
                                                 False, False, False, False,
                                                 False], dtype=bool)})
        exp2 = pd.DataFrame(
            {'W': np.array([17, 17, 5, 6, 16, 5, 7, 5,
                            4, 5, 8, 4, 5, 16, 5, 11, 4, 6]),
             'Reject null hypothesis': np.array([True, True, False, False,
                                                 True, False, False, False,
                                                 False, False, False, False,
                                                 False, True, False, False,
                                                 False, False],  dtype=bool)})
        exp3 = pd.DataFrame(
            {'W': np.array([16, 16, 17, 10, 17, 16, 16,
                            15, 15, 15, 13, 10, 10, 10,
                            9, 9, 9, 9]),
             'Reject null hypothesis': np.array([True, True, True, False,
                                                 True, True, True, True,
                                                 True, True, True, False,
                                                 False, False, False, False,
                                                 False, False], dtype=bool)})

        result1 = ancom(self.table4, self.cats4, p_adjust=None, tau=0.25)
        result2 = ancom(self.table9, self.cats9, p_adjust=None, tau=0.02)
        result3 = ancom(self.table10, self.cats10, p_adjust=None, tau=0.02)

        assert_data_frame_almost_equal(result1[0], exp1)
        assert_data_frame_almost_equal(result2[0], exp2)
        assert_data_frame_almost_equal(result3[0], exp3)

    def test_ancom_theta(self):
        result = ancom(self.table1, self.cats1, theta=0.3)
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Reject null hypothesis': np.array([True, True, False, False,
                                                 False, False, False],
                                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_alpha(self):
        result = ancom(self.table1, self.cats1, p_adjust=None, alpha=0.5)
        exp = pd.DataFrame(
            {'W': np.array([6, 6, 4, 5, 5, 4, 2]),
             'Reject null hypothesis': np.array([True, True, False, True,
                                                 True, False, False],
                                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_fail_type(self):
        with self.assertRaises(TypeError):
            ancom(self.table1.values, self.cats1)
        with self.assertRaises(TypeError):
            ancom(self.table1, self.cats1.values)

    def test_ancom_fail_zeros(self):
        with self.assertRaises(ValueError):
            ancom(self.bad1, self.cats2, p_adjust=None)

    def test_ancom_fail_negative(self):
        with self.assertRaises(ValueError):
            ancom(self.bad2, self.cats2, p_adjust=None)

    def test_ancom_fail_not_implemented_p_adjust(self):
        with self.assertRaises(ValueError):
            ancom(self.table2, self.cats2, p_adjust='fdr')

    def test_ancom_fail_missing(self):
        with self.assertRaises(ValueError):
            ancom(self.bad3, self.cats1)

        with self.assertRaises(ValueError):
            ancom(self.table1, self.badcats1)

    def test_ancom_fail_groups(self):
        with self.assertRaises(ValueError):
            ancom(self.table1, self.badcats2)

    def test_ancom_fail_size_mismatch(self):
        with self.assertRaises(ValueError):
            ancom(self.table1, self.badcats3)

    def test_ancom_fail_group_unique(self):
        with self.assertRaises(ValueError):
            ancom(self.table1, self.badcats4)

    def test_ancom_fail_1_group(self):
        with self.assertRaises(ValueError):
            ancom(self.table1, self.badcats5)

    def test_ancom_fail_tau(self):
        with self.assertRaises(ValueError):
            ancom(self.table1, self.cats1, tau=-1)
        with self.assertRaises(ValueError):
            ancom(self.table1, self.cats1, tau=1.1)

    def test_ancom_fail_theta(self):
        with self.assertRaises(ValueError):
            ancom(self.table1, self.cats1, theta=-1)
        with self.assertRaises(ValueError):
            ancom(self.table1, self.cats1, theta=1.1)

    def test_ancom_fail_alpha(self):
        with self.assertRaises(ValueError):
            ancom(self.table1, self.cats1, alpha=-1)
        with self.assertRaises(ValueError):
            ancom(self.table1, self.cats1, alpha=1.1)

    def test_ancom_fail_multiple_groups(self):
        # np.exceptions was introduced in NumPy 1.25, before which errors were
        # members of np. The following code is for backward compatibility.
        npe = getattr(np, 'exceptions', np)
        with self.assertRaises((TypeError, npe.AxisError)):
            ancom(self.table4, self.cats4,
                  significance_test="ttest_ind")


class FDRTests(TestCase):
    def test_holm_bonferroni(self):
        p = [0.005, 0.011, 0.02, 0.04, 0.13]
        obs = _calc_p_adjust("holm-bonferroni", p)
        exp = p * np.arange(1, 6)[::-1]
        for a, b in zip(obs, exp):
            self.assertAlmostEqual(a, b)

    def test_benjamini_hochberg(self):
        p = [0.005, 0.011, 0.02, 0.04, 0.13]
        obs = _calc_p_adjust("benjamini-hochberg", p)
        exp = [0.025, 0.0275, 0.03333333, 0.05, 0.13]
        for a, b in zip(obs, exp):
            self.assertAlmostEqual(a, b)


class VLRTests(TestCase):
    def setUp(self):
        self.mat = np.array([[1, 1, 2], [3, 5, 8], [13, 21, 55]])
        self.mat_neg = np.array([[-1, 1, 2], [3, -5, 8], [13, 21, -55]])
        self.mat_with_zero = np.array([[0, 1, 2], [3, 5, 8], [13, 21, 55]])

    def test_vlr(self):
        # No zeros
        output = vlr(
            x=self.mat[0],
            y=self.mat[1],
            ddof=1,
            robust=False,
        )
        self.assertAlmostEqual(output, 0.0655828061998637)

        # With zeros
        output = vlr(
            x=self.mat_with_zero[0],
            y=self.mat_with_zero[1],
            ddof=1,
            robust=False,
        )
        assert np.isnan(output)

        # assert raises error
        with self.assertRaises(ValueError):
            vlr(
                x=self.mat_neg[0],
                y=self.mat_neg[1],
                ddof=1,
                robust=False,
            )

    def test_robust_vlr(self):
        # No zeros
        output = vlr(
            x=self.mat[0],
            y=self.mat[1],
            ddof=1,
            robust=True,
        )
        self.assertAlmostEqual(output, 0.0655828061998637)

        # With zeros
        output = vlr(
            x=self.mat_with_zero[0],
            y=self.mat_with_zero[1],
            ddof=1,
            robust=True,
        )
        self.assertAlmostEqual(output, 0.024896522246558722)

    def test_pairwise_vlr(self):

        # No zeros
        dism = pairwise_vlr(self.mat, ids=None, ddof=1, robust=False)
        output = dism.condensed_form().sum()
        self.assertAlmostEqual(output, 0.2857382286903922)

        # With zeros
        with self.assertRaises(DistanceMatrixError):
            pairwise_vlr(self.mat_with_zero, ids=None, ddof=1, robust=False)

        # no validation
        dism = pairwise_vlr(self.mat, ids=None, ddof=1, robust=False,
                            validate=False)
        output = dism.data.sum() / 2
        self.assertAlmostEqual(output, 0.2857382286903922)


class DirMultTTestTests(TestCase):
    def setUp(self):
        np.random.seed(0)
        # Create sample data for testing
        self.data = {
            'feature1': [5, 8, 12, 15, 20],
            'feature2': [3, 6, 9, 12, 15],
            'feature3': [10, 15, 20, 25, 30],
        }
        self.table = pd.DataFrame(self.data)
        self.grouping = pd.Series(['Group1', 'Group1', 'Group2', 'Group2', 'Group2'])
        self.treatment = 'Group2'
        self.reference = 'Group1'

        d = 50
        n = 200
        self.depth = depth = 1000
        p1 = np.random.lognormal(0, 1, size=d) * 10
        p2 = np.random.lognormal(0.01, 1, size=d) * 10
        self.p1, self.p2 = p1 / p1.sum(), p2 / p2.sum()
        self.data2 = np.vstack(
            (
                [np.random.multinomial(depth, self.p1) for _ in range(n)],
                [np.random.multinomial(depth, self.p2) for _ in range(n)]
            )
        )
        self.table2 = pd.DataFrame(self.data2)
        self.grouping2 = pd.Series(['Group1'] * n + ['Group2'] * n)

    def test_dirmult_ttest_toy(self):
        p1 = np.array([5, 6, 7])
        p2 = np.array([4, 7, 7])
        p1, p2 = p1 / p1.sum(), p2 / p2.sum()
        depth = 1000
        n = 100
        data = np.vstack(
            (
                [np.random.multinomial(depth, p1) for _ in range(n)],
                [np.random.multinomial(depth, p2) for _ in range(n)]
            )
        )
        table = pd.DataFrame(data)
        grouping = pd.Series(['Group1'] * n + ['Group2'] * n)

        exp_lfc = np.log2([4/5, 7/6, 7/7])
        exp_lfc = (exp_lfc - exp_lfc.mean())  # convert to CLR coordinates

        res = dirmult_ttest(table, grouping, self.treatment, self.reference)

        npt.assert_array_less(exp_lfc, res['CI(97.5)'])
        npt.assert_array_less(res['CI(2.5)'], exp_lfc)

    def test_dirmult_ttest_toy_depth(self):
        p1 = np.array([5, 6, 7, 8, 9, 4])
        p2 = np.array([4, 7, 7, 6, 5, 7])
        p1, p2 = p1 / p1.sum(), p2 / p2.sum()
        depth = 100
        n = 100
        data = np.vstack(
            (
                [np.random.multinomial(depth, p1) for _ in range(n)],
                [np.random.multinomial(depth, p2) for _ in range(n)]
            )
        )
        table = pd.DataFrame(data)
        grouping = pd.Series(['Group1'] * n + ['Group2'] * n)
        exp_lfc = np.log2([4/5, 7/6, 7/7, 6/8, 5/9, 7/4])
        exp_lfc = (exp_lfc - exp_lfc.mean())  # convert to CLR coordinates
        res_100 = dirmult_ttest(table, grouping, self.treatment, self.reference)

        # increase sequencing depth by 100 fold
        depth = 10000
        data = np.vstack(
            (
                [np.random.multinomial(depth, p1) for _ in range(n)],
                [np.random.multinomial(depth, p2) for _ in range(n)]
            )
        )
        table = pd.DataFrame(data)
        res_10000 = dirmult_ttest(table, grouping, self.treatment, self.reference)

        # when the sequencing depth increases, the confidence intervals
        # should also shrink

        npt.assert_array_less(res_100['CI(2.5)'], res_10000['CI(2.5)'])
        npt.assert_array_less(res_10000['CI(97.5)'], res_100['CI(97.5)'])

    def test_dirmult_ttest_output(self):
        exp_lfc = np.log2(self.p2 / self.p1)
        exp_lfc = exp_lfc - exp_lfc.mean()
        res = dirmult_ttest(self.table2, self.grouping2,
                            self.treatment, self.reference)

        npt.assert_array_less(res['Log2(FC)'], res['CI(97.5)'])
        npt.assert_array_less(res['CI(2.5)'], res['Log2(FC)'])

        # a couple of things that complicate the tests
        # first, there is going to be a little bit of a fudge factor due
        # to the pseudocount, so we will define it via log2(0.5)
        eps = np.abs(np.log2(0.5))

        # second, the confidence interval is expected to be inaccurate
        # for (1/20) of the tests. So we should double check to
        # see if the confidence intervals were able to capture
        # 95% of the log-fold changes correctly
        self.assertGreater(np.mean(res['CI(2.5)'] - eps < exp_lfc), 0.95)
        self.assertGreater(np.mean(res['CI(97.5)'] + eps > exp_lfc), 0.95)

    def test_dirmult_ttest_valid_input(self):
        result = dirmult_ttest(self.table, self.grouping, self.treatment, self.reference)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape[1], 8)  # Expected number of columns
        pdt.assert_index_equal(result.index,
                               pd.Index(['feature1', 'feature2', 'feature3']))

    def test_dirmult_ttest_no_p_adjust(self):
        result = dirmult_ttest(self.table, self.grouping, self.treatment, self.reference,
                               p_adjust=None)
        pdt.assert_series_equal(result['pvalue'], result['qvalue'], check_names=False)

    def test_dirmult_ttest_invalid_table_type(self):
        with self.assertRaises(TypeError):
            dirmult_ttest("invalid_table", self.grouping, self.treatment, self.reference)

    def test_dirmult_ttest_invalid_grouping_type(self):
        with self.assertRaises(TypeError):
            dirmult_ttest(self.table, "invalid_grouping", self.treatment, self.reference)

    def test_dirmult_ttest_negative_values_in_table(self):
        self.table.iloc[0, 0] = -5  # Modify a value to be negative
        with self.assertRaises(ValueError):
            dirmult_ttest(self.table, self.grouping, self.treatment, self.reference)

    def test_dirmult_ttest_missing_values_in_grouping(self):
        self.grouping[1] = np.nan  # Introduce a missing value in grouping
        with self.assertRaises(ValueError):
            dirmult_ttest(self.table, self.grouping, self.treatment, self.reference)

    def test_dirmult_ttest_missing_values_in_table(self):
        self.table.iloc[2, 1] = np.nan  # Introduce a missing value in the table
        with self.assertRaises(ValueError):
            dirmult_ttest(self.table, self.grouping, self.treatment, self.reference)

    def test_dirmult_ttest_inconsistent_indexes(self):
        self.table.index = ['a', 'b', 'c', 'd', 'e']  # Change table index
        with self.assertRaises(ValueError):
            dirmult_ttest(self.table, self.grouping, self.treatment, self.reference)


if __name__ == "__main__":
    main()
