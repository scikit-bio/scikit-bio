# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main, skipIf

import numpy as np
import numpy.testing as npt

from skbio.stats.composition import (
    closure, clr, clr_inv, ilr, ilr_inv, alr, alr_inv, sbp_basis, _gram_schmidt_basis
)

import array_api_compat as aac
import warnings

try:
    import torch
    no_torch = False
except ImportError:
    no_torch = True

import subprocess
def no_gpu_available():
    try:
        result = subprocess.run(
            ["nvidia-smi"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return not result.returncode == 0
    except FileNotFoundError:
        return True

# class name_space_config(TestCase)
Array = object


def assert_allclose(x:Array, y:Array, rtol=1e-07, atol=1e-7):
    # check namespace eq
    if not isinstance(y,type(x)):
        raise TypeError("x and y different types")
    if x.dtype is not y.dtype:
        warnings.warn(
                "different dtypes",
                UserWarning,
            )
    if not x.shape==y.shape: # shapes are tuples not array
        raise ValueError(f"x and y different shapes, {x.shape} at\
            the first array but {y.shape} at the second ")
    xp = aac.array_namespace(x)
    if not xp.all(xp.abs(x-y)<=(atol+rtol*xp.abs(y))):
        raise ValueError(f'Not equal to tolerance rtol={rtol}, atol={atol}')


@skipIf(no_torch, "Skipping all tests: no torch dependency")
class torch_cpu(TestCase):
    def setUp(self):
        # args
        self.namespace = torch
        self.namespace_asarray = torch.tensor
        self.device = "cpu"

        # data
        # 2d ndarray, shape (2,3)
        self.cdata1 = self.namespace_asarray([[2, 2, 6],
                                        [4, 4, 2]], device = self.device)
        # 1d ndarray, shape (3,)
        self.cdata2 = self.namespace_asarray([2, 2, 6], device = self.device)
        self.cdata3 = self.namespace_asarray([[1, 2, 3, 0, 5],
                                         [1, 0, 0, 4, 5],
                                         [1, 2, 3, 4, 5]],
                                        device = self.device)

        self.cdata4 = self.namespace_asarray([1, 2, 3, 0, 5],
                                        device = self.device)

        # 2d nested list, shape (2,3), containing zeros
        self.cdata5 = [[2, 2, 6], [4, 4, 2]]
        # 2d nested list, shape (3,5), containing zeros
        self.cdata6 = [[1, 2, 3, 0, 5],
                       [1, 0, 0, 4, 5],
                       [1, 2, 3, 4, 5]]

        self.cdata7 = [np.exp(1), 1, 1]
        self.cdata8 = [np.exp(1), 1, 1, 1]

        # Simplicial orthonormal basis obtained from Gram-Schmidt
        self.ortho1 = self.namespace_asarray(
                    [[0.44858053, 0.10905743, 0.22118102, 0.22118102],
                     [0.3379924, 0.3379924, 0.0993132, 0.22470201],
                     [0.3016453, 0.3016453, 0.3016453, 0.09506409]],
                    device=self.device
                    )

        # Real data
        self.rdata1 = self.namespace_asarray([
                        [0.70710678, -0.70710678, 0., 0.],
                        [0.40824829, 0.40824829, -0.81649658, 0.],
                        [0.28867513, 0.28867513, 0.28867513, -0.8660254]],
                                             device=self.device)

        # Bad datasets
        # negative count
        self.bat1 = self.namespace_asarray([1, 2, -1],
                                           device=self.device)
        # zero count
        self.bad2 = self.namespace_asarray([[[1, 2, 3, 0, 5]]],
                                           device=self.device)
        # singleton array
        self.bad3 = self.namespace_asarray([[1], [2], [3], [4], [5]],
                                      device=self.device)

    def test_closure(self):
        self.cdata1 = self.namespace_asarray(
                                        [[2, 2, 6],
                                         [4, 4, 2]],
                                        device = self.device)
        rst_1 = closure(self.cdata1)
        rst_1_ = self.namespace_asarray([[.2, .2, .6],
                                      [.4, .4, .2]], device = self.device)
        assert_allclose(rst_1, rst_1_)

        self.cdata2 = self.namespace_asarray([2, 2, 6], device = self.device)
        rst_2 = closure(self.cdata2)
        rst_2_ = self.namespace_asarray([.2, .2, .6], device = self.device)
        assert_allclose(rst_2, rst_2_)

        self.cdata5 = [[2, 2, 6], [4, 4, 2]]
        rst_3 = closure(self.cdata5)
        rst_3_ = np.array([[.2, .2, .6],
                        [.4, .4, .2]])
        assert_allclose(rst_3, rst_3_)

        with self.assertRaises(ValueError):
            self.bad1 = self.namespace_asarray([1, 2, -1],
                                               device=self.device)
            closure(self.bad1)

        # make sure that inplace modification is not occurring
        closure(self.cdata2)
        rst_4_ = self.namespace_asarray([2, 2, 6],
                                        device=self.device)
        assert_allclose(self.cdata2, rst_4_)

    def test_closure_warning(self):
        with self.assertRaises(ValueError):
            closure([0., 0., 0.])

        with self.assertRaises(ValueError):
            closure([[0., 0., 0.],
                     [0., 5., 5.]])

    def test_clr(self):
        #
        self.cdata1 = self.namespace_asarray([[2, 2, 6],
                                              [4, 4, 2]],
                                             device = self.device)
        cmat = clr(closure(self.cdata1))
        A = np.array([.2, .2, .6])
        B = np.array([.4, .4, .2])

        cmat_ = self.namespace_asarray([np.log(A / np.exp(np.log(A).mean())),
                                            np.log(B / np.exp(np.log(B).mean()))],
                                        device = self.device)
        assert_allclose(cmat, cmat_)

        #
        self.cdata2 = self.namespace_asarray([2, 2, 6], device = self.device)
        cmat = clr(closure(self.cdata2))
        A = np.array([.2, .2, .6])
        cmat_ = self.namespace_asarray(np.log(A / np.exp(np.log(A).mean())),
                                       device=self.device)
        assert_allclose(cmat, cmat_)

        #
        self.cdata5 = [[2, 2, 6], [4, 4, 2]]
        cmat = clr(closure(self.cdata5))
        A = np.array([.2, .2, .6])
        B = np.array([.4, .4, .2])
        npt.assert_allclose(cmat,
                            [np.log(A / np.exp(np.log(A).mean())),
                             np.log(B / np.exp(np.log(B).mean()))])

        #
        self.bad1 = self.namespace_asarray([1, 2, -1],
                                           device=self.device)
        with self.assertRaises(ValueError):
            clr(self.bad1)
        self.bad2 = self.namespace_asarray([[[1, 2, 3, 0, 5]]],
                                           device=self.device)
        with self.assertRaises(ValueError):
            clr(self.bad2)

        # make sure that inplace modification is not occurring
        clr(self.cdata2)
        rst = self.namespace_asarray([2, 2, 6], device=self.device)
        npt.assert_allclose(self.cdata2, rst)

    def test_clr_inv(self):
        self.rdata1 = self.namespace_asarray([
                        [0.70710678, -0.70710678, 0., 0.],
                        [0.40824829, 0.40824829, -0.81649658, 0.],
                        [0.28867513, 0.28867513, 0.28867513, -0.8660254]],
                                             device=self.device)
        rst_1 = clr_inv(self.rdata1)
        self.ortho1 = self.namespace_asarray(
                    [[0.44858053, 0.10905743, 0.22118102, 0.22118102],
                     [0.3379924, 0.3379924, 0.0993132, 0.22470201],
                     [0.3016453, 0.3016453, 0.3016453, 0.09506409]],
                    device=self.device
                    )
        assert_allclose(rst_1, self.ortho1)
        rst_2 = clr(clr_inv(self.rdata1))
        assert_allclose(rst_2, self.rdata1,
                        rtol=1e-4, atol=1e-5)

        # make sure that inplace modification is not occurring
        clr_inv(self.rdata1)
        rst_3_ = self.namespace_asarray([[0.70710678, -0.70710678, 0., 0.],
                                      [0.40824829, 0.40824829, -0.81649658, 0.],
                                      [0.28867513, 0.28867513,0.28867513, -0.8660254]],
                                        device= self.device)
        assert_allclose(self.rdata1, rst_3_)

    def test_ilr(self):
        self.cdata7 = self.namespace_asarray([np.exp(1), 1, 1],
                                             device=self.device)
        mat = closure(self.cdata7)
        mat_ = self.namespace_asarray([0.70710678, 0.40824829],
                                      device=self.device)
        assert_allclose(ilr(mat), mat_)

        # Should give same result as inner
        self.ortho1 = self.namespace_asarray(
                    [[0.44858053, 0.10905743, 0.22118102, 0.22118102],
                     [0.3379924, 0.3379924, 0.0993132, 0.22470201],
                     [0.3016453, 0.3016453, 0.3016453, 0.09506409]],
                    device=self.device
                    )
        eyes = self.namespace_asarray(np.identity(3),
                                        device=self.device)
        assert_allclose(ilr(self.ortho1), eyes)

        # no check
        self.cdata1 = self.namespace_asarray(
                                        [[2, 2, 6],
                                         [4, 4, 2]],
                                        device = self.device)
        assert_allclose(ilr(mat, validate=False),
                        self.namespace_asarray([0.70710678, 0.40824829],
                                               device=self.device)
                        )

        with self.assertRaises(ValueError):
            ilr(self.cdata1, basis=self.cdata1)

        # make sure that inplace modification is not occurring
        ilr(self.cdata1)
        rst_1_ = self.namespace_asarray([[2, 2, 6],
                                         [4, 4, 2]],
                                        device=self.device)
        assert_allclose(self.cdata1, rst_1_)

    def test_ilr_basis(self):
        table = self.namespace_asarray([[1., 10.],
                                        [1.14141414, 9.90909091],
                                        [1.28282828, 9.81818182],
                                        [1.42424242, 9.72727273],
                                        [1.56565657, 9.63636364]],
                                       device=self.device)
        basis = self.namespace_asarray([[0.80442968, 0.19557032]],
                                     device=self.device)
        basis = clr(basis).reshape(1,-1) # NOTE: _check_org will do the reshape
                                        # however ilr will evaluate
        res = ilr(table, basis=basis)
        exp = self.namespace_asarray([[np.log(1/10)*np.sqrt(1/2)],
                        [np.log(1.14141414 / 9.90909091)*np.sqrt(1/2)],
                        [np.log(1.28282828 / 9.81818182)*np.sqrt(1/2)],
                        [np.log(1.42424242 / 9.72727273)*np.sqrt(1/2)],
                        [np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]],
                                     device=self.device)
        assert_allclose(res, exp)

    def test_ilr_basis_one_dimension_error(self):
        table = self.namespace_asarray([[1., 10.],
                          [1.14141414, 9.90909091],
                          [1.28282828, 9.81818182],
                          [1.42424242, 9.72727273],
                          [1.56565657, 9.63636364]],
                                       device=self.device)
        basis = self.namespace_asarray([0.80442968, 0.19557032],
                                       device=self.device)
        basis = clr(basis) # clr: has squeeze
        with self.assertRaises(ValueError):
            ilr(table, basis=basis)

    def test_ilr_inv(self):
        self.cdata7 = self.namespace_asarray([np.exp(1), 1, 1],
                                             device=self.device)
        mat = closure(self.cdata7)
        assert_allclose(ilr_inv(ilr(mat)), mat)
        eye = self.namespace_asarray(np.identity(3),
                                     device=self.device)
        self.ortho1 = self.namespace_asarray(
                    [[0.44858053, 0.10905743, 0.22118102, 0.22118102],
                     [0.3379924, 0.3379924, 0.0993132, 0.22470201],
                     [0.3016453, 0.3016453, 0.3016453, 0.09506409]],
                    device=self.device
                    )
        assert_allclose(ilr_inv(eye), self.ortho1)

        # no check
        assert_allclose(ilr_inv(ilr(mat), validate=False), mat)

        #
        self.cdata1 = self.namespace_asarray([[2, 2, 6],
                                              [4, 4, 2]],
                                            device = self.device)
        with self.assertRaises(ValueError):
            ilr_inv(self.cdata1,
                    basis=self.cdata1)

        # make sure that inplace modification is not occurring
        ilr_inv(self.cdata1)
        rst_ =  self.namespace_asarray([[2, 2, 6],
                                        [4, 4, 2]],
                                       device=self.device)
        assert_allclose(self.cdata1, rst_)

    def test_ilr_basis_isomorphism(self):
        # tests to make sure that the isomorphism holds
        # with the introduction of the basis.
        basis = self.namespace_asarray([[0.80442968, 0.19557032]],
                                       device=self.device)
        basis = clr(basis).reshape(1,-1)
        table = self.namespace_asarray([[np.log(1/10)*np.sqrt(1/2),
                           np.log(1.14141414 / 9.90909091)*np.sqrt(1/2),
                           np.log(1.28282828 / 9.81818182)*np.sqrt(1/2),
                           np.log(1.42424242 / 9.72727273)*np.sqrt(1/2),
                           np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]],
                                       device=self.device
                                       ).T
        lr = ilr_inv(table, basis=basis)
        res = ilr(lr, basis=basis)
        assert_allclose(res, table)

        #
        table = self.namespace_asarray([[1., 10.],
                          [1.14141414, 9.90909091],
                          [1.28282828, 9.81818182],
                          [1.42424242, 9.72727273],
                          [1.56565657, 9.63636364]],
                                       device=self.device)

        res = ilr_inv(ilr(table, basis=basis), basis=basis)
        assert_allclose(res, closure(table))

    def test_ilr_inv_basis(self):
        exp = self.namespace_asarray([[1., 10.],
                                [1.14141414, 9.90909091],
                                [1.28282828, 9.81818182],
                                [1.42424242, 9.72727273],
                                [1.56565657, 9.63636364]],
                                     device=self.device)
        exp = closure(exp)
        basis = self.namespace_asarray([[0.80442968, 0.19557032]],
                                       device=self.device)
        basis = clr(basis).reshape(1,-1)
        table = self.namespace_asarray([[np.log(1/10)*np.sqrt(1/2),
                           np.log(1.14141414 / 9.90909091)*np.sqrt(1/2),
                           np.log(1.28282828 / 9.81818182)*np.sqrt(1/2),
                           np.log(1.42424242 / 9.72727273)*np.sqrt(1/2),
                           np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]],
                                       device=self.device).T

        res = ilr_inv(table, basis=basis)
        assert_allclose(res, exp)

    def test_ilr_inv_basis_one_dimension_error(self):
        basis = self.namespace_asarray([0.80442968, 0.19557032],
                                       device=self.device)
        basis = clr([0.80442968, 0.19557032])
        table = self.namespace_asarray([[np.log(1/10)*np.sqrt(1/2),
                           np.log(1.14141414 / 9.90909091)*np.sqrt(1/2),
                           np.log(1.28282828 / 9.81818182)*np.sqrt(1/2),
                           np.log(1.42424242 / 9.72727273)*np.sqrt(1/2),
                           np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]],
                                       device=self.device).T
        with self.assertRaises(ValueError):
            ilr_inv(table, basis=basis)

    def test_alr(self):
        # 2d-composition
        self.cdata1 = self.namespace_asarray([[2, 2, 6],
                                              [4, 4, 2]],
                                            device = self.device)
        comp1 = closure(self.cdata1)
        alr2d_method = alr(comp1, denominator_idx=1)

        comp_1_byhand = np.array([[2, 2, 6],\
                                  [4, 4, 2]])
        alr2d_byhand = self.namespace_asarray(
                        [np.log(comp_1_byhand[:, 0]/comp_1_byhand[:, 1]),
                            np.log(comp_1_byhand[:, 2]/comp_1_byhand[:, 1])],
                        device=self.device
                    ).T

        assert_allclose(alr2d_byhand, alr2d_method)

        # 1d-composition
        self.cdata2 = self.namespace_asarray([2, 2, 6],
                                             device = self.device)
        comp2 = closure(self.cdata2)
        alr1d_byhand = self.namespace_asarray([np.log(comp2[0]/comp2[1]),
                                                np.log(comp2[2]/comp2[1])],
                                              device=self.device
                                              ).T
        alr1d_method = alr(comp2, denominator_idx=1)
        assert_allclose(alr1d_byhand, alr1d_method)

        #
        self.bad1 = self.namespace_asarray([1, 2, -1],
                                           device=self.device)
        with self.assertRaises(ValueError):
            alr(self.bad1)
        self.bad2 = self.namespace_asarray([[[1, 2, 3, 0, 5]]],
                                           device=self.device)
        with self.assertRaises(ValueError):
            alr(self.bad2)

        # make sure that inplace modification is not occurring
        self.cdata2 = self.namespace_asarray([2, 2, 6],
                                             device = self.device)
        alr(self.cdata2)
        rst_ = self.namespace_asarray([2, 2, 6],
                                      device=self.device)
        assert_allclose(self.cdata2, rst_)

    def test_alr_inv(self):
        # 2d-composition
        self.cdata1 = self.namespace_asarray(
                                        [[2, 2, 6],
                                         [4, 4, 2]],
                                        device = self.device)
        comp1 = closure(self.cdata1)
        alr2d_method = alr(comp1, denominator_idx=1)
        alrinv2d_method = alr_inv(alr2d_method, denominator_idx=1)

        comp1_byhand = closure(np.array([[2, 2, 6],
                                 [4, 4, 2]]))
        alr2d_byhand = np.array([np.log(comp1_byhand[:, 0]/comp1_byhand[:, 1]),
                                 np.log(comp1_byhand[:, 2]/comp1_byhand[:, 1])],
                                ).T
        B = 1/(1 + np.exp(alr2d_byhand[:, 0]) + np.exp(alr2d_byhand[:, 1]))
        A = B * np.exp(alr2d_byhand[:, 0])
        C = B * np.exp(alr2d_byhand[:, 1])
        alrinv2d_byhand = np.column_stack((A, B, C))
        alrinv2d_byhand = self.namespace_asarray(alrinv2d_byhand,
                                                 device=self.device)
        assert_allclose(alrinv2d_byhand, alrinv2d_method)

        # 1d-composition
        self.cdata2 = self.namespace_asarray([2, 2, 6],
                                             device = self.device)
        comp2 = closure(self.cdata2)
        alr1d_method = alr(comp2, denominator_idx=1)
        alrinv1d_method = alr_inv(alr1d_method, denominator_idx=1)

        alr1d_byhand = np.array([np.log(comp2[0]/comp2[1]),
                                 np.log(comp2[2]/comp2[1])]).T
        B = 1/(1 + np.exp(alr1d_byhand[0]) + np.exp(alr1d_byhand[1]))
        A = B * np.exp(alr1d_byhand[0])
        C = B * np.exp(alr1d_byhand[1])
        alrinv1d_byhand = np.column_stack((A, B, C))[0, :]
        alrinv1d_byhand = self.namespace_asarray(alrinv1d_byhand,
                                                 device=self.device)
        assert_allclose(alrinv1d_byhand, alrinv1d_method)

        # make sure that inplace modification is not occurring
        self.rdata1 = self.namespace_asarray([
                        [0.70710678, -0.70710678, 0., 0.],
                        [0.40824829, 0.40824829, -0.81649658, 0.],
                        [0.28867513, 0.28867513, 0.28867513, -0.8660254]],
                                             device=self.device)
        alr_inv(self.rdata1)
        rst = self.namespace_asarray([[0.70710678, -0.70710678,\
                                            0., 0.],
                                      [0.40824829, 0.40824829,\
                                            -0.81649658, 0.],
                                      [0.28867513, 0.28867513,\
                                            0.28867513, -0.8660254]],
                                     device=self.device)
        assert_allclose(self.rdata1, rst)

    def test_sbp_basis_gram_schmidt(self):
        gsbasis = _gram_schmidt_basis(5)
        gsbasis = self.namespace_asarray(gsbasis,
                                         device=self.device)
        sbp = np.array([[1, -1, 0, 0, 0],
                        [1, 1, -1, 0, 0],
                        [1, 1, 1, -1, 0],
                        [1, 1, 1, 1, -1]])
        sbpbasis = sbp_basis(sbp)
        sbpbasis = self.namespace_asarray(sbpbasis,
                                         device=self.device)
        assert_allclose(gsbasis, sbpbasis)

    # def test_sbp_basis_elementwise(self):
    #     sbp = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1],
    #                     [1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, 0],
    #                     [1, 1, 1, 1, 1, 1, -1, 0, 0, 0, 0, 0],
    #                     [1, 1, -1, -1, -1, 1, 0, 0, 0, 0, 0, 0],
    #                     [1, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    #                     [1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0],
    #                     [0, 0, 1, -1, -1, 0, 0, 0, 0, 0, 0, 0],
    #                     [0, 0, 0, 1, -1, 0, 0, 0, 0, 0, 0, 0],
    #                     [0, 0, 0, 0, 0, 0, 0, 1, -1, -1, 1, 0],
    #                     [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, -1, 0],
    #                     [0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0]])
    #     sbpbasis = sbp_basis(sbp)
    #     # by hand, element-wise
    #     r = np.apply_along_axis(func1d=lambda x: np.sum(x > 0),
    #                             axis=1, arr=sbp)
    #     s = np.apply_along_axis(func1d=lambda x: np.sum(x < 0),
    #                             axis=1, arr=sbp)
    #     psi = np.zeros(sbp.shape)
    #     for i in range(0, sbp.shape[0]):
    #         for j in range(0, sbp.shape[1]):
    #             if sbp[i, j] == 1:
    #                 psi[i, j] = np.sqrt(s[i]/(r[i]*(r[i]+s[i])))
    #             elif sbp[i, j] == -1:
    #                 psi[i, j] = -np.sqrt(r[i]/(s[i]*(r[i]+s[i])))
    #     npt.assert_allclose(psi, sbpbasis)


@skipIf(
    no_gpu_available() or no_torch,
    "Skipping all tests: no GPU available or no PyTorch dependency",
)
class torch_cuda(TestCase):
    def setUp(self):
        # args
        self.namespace = torch
        self.namespace_asarray = torch.tensor
        self.device = "cuda"

        # data
        # 2d ndarray, shape (2,3)
        self.cdata1 = self.namespace_asarray([[2, 2, 6],
                                        [4, 4, 2]], device = self.device)
        # 1d ndarray, shape (3,)
        self.cdata2 = self.namespace_asarray([2, 2, 6], device = self.device)
        self.cdata3 = self.namespace_asarray([[1, 2, 3, 0, 5],
                                         [1, 0, 0, 4, 5],
                                         [1, 2, 3, 4, 5]],
                                        device = self.device)

        self.cdata4 = self.namespace_asarray([1, 2, 3, 0, 5],
                                        device = self.device)

        # 2d nested list, shape (2,3), containing zeros
        self.cdata5 = [[2, 2, 6], [4, 4, 2]]
        # 2d nested list, shape (3,5), containing zeros
        self.cdata6 = [[1, 2, 3, 0, 5],
                       [1, 0, 0, 4, 5],
                       [1, 2, 3, 4, 5]]

        self.cdata7 = [np.exp(1), 1, 1]
        self.cdata8 = [np.exp(1), 1, 1, 1]

        # Simplicial orthonormal basis obtained from Gram-Schmidt
        self.ortho1 = self.namespace_asarray(
                    [[0.44858053, 0.10905743, 0.22118102, 0.22118102],
                     [0.3379924, 0.3379924, 0.0993132, 0.22470201],
                     [0.3016453, 0.3016453, 0.3016453, 0.09506409]],
                    device=self.device
                    )

        # Real data
        self.rdata1 = [[0.70710678, -0.70710678, 0., 0.],
                       [0.40824829, 0.40824829, -0.81649658, 0.],
                       [0.28867513, 0.28867513, 0.28867513, -0.8660254]]

        # Bad datasets
        # negative count
        self.bat1 = self.namespace_asarray([1, 2, -1],
                                           device=self.device)
        # zero count
        self.bad2_np = np.array([[[1, 2, 3, 0, 5]]])
        self.bad2 = self.namespace_asarray([[[1, 2, 3, 0, 5]]],
                                           device=self.device)
        # singleton array
        self.bad3 = self.namespace_asarray([[1], [2], [3], [4], [5]],
                                      device=self.device)

    def test_closure(self):
        self.cdata1 = self.namespace_asarray(
                                        [[2, 2, 6],
                                         [4, 4, 2]],
                                        device = self.device)
        rst_1 = closure(self.cdata1)
        rst_1_ = self.namespace_asarray([[.2, .2, .6],
                                      [.4, .4, .2]], device = self.device)
        assert_allclose(rst_1, rst_1_)

        self.cdata2 = self.namespace_asarray([2, 2, 6], device = self.device)
        rst_2 = closure(self.cdata2)
        rst_2_ = self.namespace_asarray([.2, .2, .6], device = self.device)
        assert_allclose(rst_2, rst_2_)

        self.cdata5 = [[2, 2, 6], [4, 4, 2]]
        rst_3 = closure(self.cdata5)
        rst_3_ = np.array([[.2, .2, .6],
                        [.4, .4, .2]])
        assert_allclose(rst_3, rst_3_)

        with self.assertRaises(ValueError):
            self.bad1 = self.namespace_asarray([1, 2, -1],
                                               device=self.device)
            closure(self.bad1)

        with self.assertRaises(ValueError):
            self.bad2 = self.namespace_asarray([[[1, 2, 3, 0, 5]]],
                                               device=self.device)
            closure(self.bad2)

        # make sure that inplace modification is not occurring
        closure(self.cdata2)
        rst_4_ = self.namespace_asarray([2, 2, 6],
                                        device=self.device)
        assert_allclose(self.cdata2, rst_4_)

    def test_closure_warning(self):
        with self.assertRaises(ValueError):
            closure([0., 0., 0.])

        with self.assertRaises(ValueError):
            closure([[0., 0., 0.],
                     [0., 5., 5.]])

    def test_clr(self):
        #
        self.cdata1 = self.namespace_asarray([[2, 2, 6],
                                        [4, 4, 2]], device = self.device)
        cmat = clr(closure(self.cdata1))
        A = np.array([.2, .2, .6])
        B = np.array([.4, .4, .2])

        cmat_ = self.namespace_asarray([np.log(A / np.exp(np.log(A).mean())),
                                            np.log(B / np.exp(np.log(B).mean()))],
                                        device = self.device)
        assert_allclose(cmat, cmat_)

        #
        self.cdata2 = self.namespace_asarray([2, 2, 6], device = self.device)
        cmat = clr(closure(self.cdata2))
        A = np.array([.2, .2, .6])
        cmat_ = self.namespace_asarray(np.log(A / np.exp(np.log(A).mean())),
                                       device=self.device)
        assert_allclose(cmat, cmat_)

        #
        self.cdata5 = [[2, 2, 6], [4, 4, 2]]
        cmat = clr(closure(self.cdata5))
        A = np.array([.2, .2, .6])
        B = np.array([.4, .4, .2])
        npt.assert_allclose(cmat,
                            [np.log(A / np.exp(np.log(A).mean())),
                             np.log(B / np.exp(np.log(B).mean()))])

        #
        self.bad1 = self.namespace_asarray([1, 2, -1],
                                           device=self.device)
        with self.assertRaises(ValueError):
            clr(self.bad1)
        self.bad2 = self.namespace_asarray([[[1, 2, 3, 0, 5]]],
                                           device=self.device)
        with self.assertRaises(ValueError):
            clr(self.bad2)

        # make sure that inplace modification is not occurring
        clr(self.cdata2)
        rst = self.namespace_asarray([2, 2, 6], device=self.device)
        assert_allclose(self.cdata2, rst)

    def test_clr_inv(self):
        self.rdata1 = self.namespace_asarray([
                        [0.70710678, -0.70710678, 0., 0.],
                        [0.40824829, 0.40824829, -0.81649658, 0.],
                        [0.28867513, 0.28867513, 0.28867513, -0.8660254]],
                                             device=self.device)
        rst_1 = clr_inv(self.rdata1)
        self.ortho1 = self.namespace_asarray(
                    [[0.44858053, 0.10905743, 0.22118102, 0.22118102],
                     [0.3379924, 0.3379924, 0.0993132, 0.22470201],
                     [0.3016453, 0.3016453, 0.3016453, 0.09506409]],
                    device=self.device

        )
        assert_allclose(rst_1, self.ortho1)
        rst_2 = clr(clr_inv(self.rdata1))
        assert_allclose(rst_2, self.rdata1,
                        rtol=1e-4, atol=1e-5)

        # make sure that inplace modification is not occurring
        clr_inv(self.rdata1)
        rst_3_ = self.namespace_asarray([[0.70710678, -0.70710678, 0., 0.],
                                      [0.40824829, 0.40824829, -0.81649658, 0.],
                                      [0.28867513, 0.28867513,0.28867513, -0.8660254]],
                                        device= self.device)
        assert_allclose(self.rdata1, rst_3_)

    def test_ilr(self):
        self.cdata7 = self.namespace_asarray([np.exp(1), 1, 1],
                                             device=self.device)
        mat = closure(self.cdata7)
        mat_ = self.namespace_asarray([0.70710678, 0.40824829],
                                      device=self.device)
        assert_allclose(ilr(mat), mat_)

        # Should give same result as inner
        self.ortho1 = self.namespace_asarray(
                    [[0.44858053, 0.10905743, 0.22118102, 0.22118102],
                     [0.3379924, 0.3379924, 0.0993132, 0.22470201],
                     [0.3016453, 0.3016453, 0.3016453, 0.09506409]],
                    device=self.device
                    )
        eyes = self.namespace_asarray(np.identity(3),
                                        device=self.device)
        assert_allclose(ilr(self.ortho1), eyes)

        # no check
        self.cdata1 = self.namespace_asarray([[2, 2, 6],
                                              [4, 4, 2]],
                                        device = self.device,
                                        dtype=torch.float32)
        assert_allclose(ilr(mat, validate=False),
                        self.namespace_asarray([0.70710678, 0.40824829],
                                               device=self.device)
                        )

        with self.assertRaises(ValueError):
            ilr(self.cdata1, basis=self.cdata1)

        # make sure that inplace modification is not occurring
        ilr(self.cdata1)
        rst_1_ = self.namespace_asarray([[2, 2, 6],
                                         [4, 4, 2]],
                                        device=self.device)
        assert_allclose(self.cdata1, rst_1_)

    def test_ilr_basis(self):
        table = self.namespace_asarray([[1., 10.],
                                        [1.14141414, 9.90909091],
                                        [1.28282828, 9.81818182],
                                        [1.42424242, 9.72727273],
                                        [1.56565657, 9.63636364]],
                                       device=self.device)
        basis = self.namespace_asarray([[0.80442968, 0.19557032]],
                                     device=self.device)
        basis = clr(basis).reshape(1,-1)
        res = ilr(table, basis=basis)
        exp = self.namespace_asarray([[np.log(1/10)*np.sqrt(1/2)],
                        [np.log(1.14141414 / 9.90909091)*np.sqrt(1/2)],
                        [np.log(1.28282828 / 9.81818182)*np.sqrt(1/2)],
                        [np.log(1.42424242 / 9.72727273)*np.sqrt(1/2)],
                        [np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]],
                                     device=self.device)
        assert_allclose(res, exp)

    def test_ilr_basis_one_dimension_error(self):
        table = self.namespace_asarray([[1., 10.],
                          [1.14141414, 9.90909091],
                          [1.28282828, 9.81818182],
                          [1.42424242, 9.72727273],
                          [1.56565657, 9.63636364]],
                                       device=self.device)
        basis = self.namespace_asarray([[0.80442968, 0.19557032]],
                                       device=self.device)
        with self.assertRaises(ValueError):
            ilr(table, basis=basis)

    def test_ilr_inv(self):
        self.cdata7 = self.namespace_asarray([np.exp(1), 1, 1],
                                             device=self.device)
        mat = closure(self.cdata7)
        assert_allclose(ilr_inv(ilr(mat)), mat)
        eye = self.namespace_asarray(np.identity(3),
                                     device=self.device)
        self.ortho1 = self.namespace_asarray(
                    [[0.44858053, 0.10905743, 0.22118102, 0.22118102],
                     [0.3379924, 0.3379924, 0.0993132, 0.22470201],
                     [0.3016453, 0.3016453, 0.3016453, 0.09506409]],
                    device=self.device
                    )
        assert_allclose(ilr_inv(eye), self.ortho1)

        # no check
        assert_allclose(ilr_inv(ilr(mat), validate=False), mat)

        #
        self.cdata1 = self.namespace_asarray(
                                            [[2, 2, 6],
                                             [4, 4, 2]],
                                            device = self.device,
                                            dtype=torch.float32)
        with self.assertRaises(ValueError):
            ilr_inv(self.cdata1, basis=self.cdata1)

        # make sure that inplace modification is not occurring
        ilr_inv(self.cdata1)
        rst_ =  self.namespace_asarray([[2, 2, 6],
                                        [4, 4, 2]],
                                       device=self.device)
        assert_allclose(self.cdata1, rst_)

    def test_ilr_basis_isomorphism(self):
        # tests to make sure that the isomorphism holds
        # with the introduction of the basis.
        basis = self.namespace_asarray([[0.80442968, 0.19557032]],
                                       device=self.device)
        basis = clr(basis).reshape(1,-1)
        table = self.namespace_asarray([[np.log(1/10)*np.sqrt(1/2),
                           np.log(1.14141414 / 9.90909091)*np.sqrt(1/2),
                           np.log(1.28282828 / 9.81818182)*np.sqrt(1/2),
                           np.log(1.42424242 / 9.72727273)*np.sqrt(1/2),
                           np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]],
                                       device=self.device
                                       ).T
        lr = ilr_inv(table, basis=basis)
        res = ilr(lr, basis=basis)
        assert_allclose(res, table)

        #
        table = self.namespace_asarray([[1., 10.],
                          [1.14141414, 9.90909091],
                          [1.28282828, 9.81818182],
                          [1.42424242, 9.72727273],
                          [1.56565657, 9.63636364]],
                                       device=self.device)

        res = ilr_inv(ilr(table, basis=basis), basis=basis)
        assert_allclose(res, closure(table))

    def test_ilr_inv_basis(self):
        exp = self.namespace_asarray([[1., 10.],
                                [1.14141414, 9.90909091],
                                [1.28282828, 9.81818182],
                                [1.42424242, 9.72727273],
                                [1.56565657, 9.63636364]],
                                     device=self.device)
        exp = closure(exp)
        basis = self.namespace_asarray([[0.80442968, 0.19557032]],
                                       device=self.device)
        basis = clr(basis).reshape(1,-1)
        table = self.namespace_asarray([[np.log(1/10)*np.sqrt(1/2),
                           np.log(1.14141414 / 9.90909091)*np.sqrt(1/2),
                           np.log(1.28282828 / 9.81818182)*np.sqrt(1/2),
                           np.log(1.42424242 / 9.72727273)*np.sqrt(1/2),
                           np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]],
                                       device=self.device).T

        res = ilr_inv(table, basis=basis)
        assert_allclose(res, exp)

    def test_ilr_inv_basis_one_dimension_error(self):
        basis = self.namespace_asarray([0.80442968, 0.19557032],
                                       device=self.device)
        basis = clr([0.80442968, 0.19557032])
        table = self.namespace_asarray([[np.log(1/10)*np.sqrt(1/2),
                           np.log(1.14141414 / 9.90909091)*np.sqrt(1/2),
                           np.log(1.28282828 / 9.81818182)*np.sqrt(1/2),
                           np.log(1.42424242 / 9.72727273)*np.sqrt(1/2),
                           np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]],
                                       device=self.device).T
        with self.assertRaises(ValueError):
            ilr_inv(table, basis=basis)

    def test_alr(self):
        # 2d-composition
        self.cdata1 = self.namespace_asarray(
                                            [[2, 2, 6],
                                             [4, 4, 2]],
                                            device = self.device)
        comp1 = closure(self.cdata1)
        alr2d_method = alr(comp1, denominator_idx=1)

        comp1_byhand = closure(np.array(
                                [[2, 2, 6],
                                [4, 4, 2]])) # the support for closure
                                             # is tested in another file
        alr2d_byhand = self.namespace_asarray(
                        [np.log(comp1_byhand[:, 0]/comp1_byhand[:, 1]),
                         np.log(comp1_byhand[:, 2]/comp1_byhand[:, 1])],
                        device=self.device
                    ).T

        assert_allclose(alr2d_byhand, alr2d_method)

        # 1d-composition
        self.cdata2 = self.namespace_asarray([2, 2, 6],
                                             device = self.device)
        comp2 = closure(self.cdata2)
        alr1d_method = alr(comp2, denominator_idx=1)
        comp2_byhand = closure(np.array([2, 2, 6]))
        alr1d_byhand = self.namespace_asarray([
                        np.log(comp2_byhand[0]/comp2_byhand[1]),
                        np.log(comp2_byhand[2]/comp2_byhand[1])],
                                              device=self.device
                                              ).T

        assert_allclose(alr1d_byhand, alr1d_method)

        #
        self.bad1 = self.namespace_asarray([1, 2, -1],
                                           device=self.device)
        with self.assertRaises(ValueError):
            alr(self.bad1)
        self.bad2 = self.namespace_asarray([[[1, 2, 3, 0, 5]]],
                                           device=self.device)
        with self.assertRaises(ValueError):
            alr(self.bad2)

        # make sure that inplace modification is not occurring
        self.cdata2 = self.namespace_asarray([2, 2, 6],
                                             device = self.device)
        alr(self.cdata2)
        rst_ = self.namespace_asarray([2, 2, 6],
                                      device=self.device)
        assert_allclose(self.cdata2, rst_)

    def test_alr_inv(self):
        # 2d-composition
        self.cdata1 = self.namespace_asarray(
                                        [[2, 2, 6],
                                         [4, 4, 2]],
                                        device = self.device)
        comp1 = closure(self.cdata1)
        alr2d_method = alr(comp1, denominator_idx=1)
        alrinv2d_method = alr_inv(alr2d_method, denominator_idx=1)

        comp1_byhand = closure(np.array([[2, 2, 6],
                                  [4, 4, 2]]))
        alr2d_byhand = np.array([np.log(comp1_byhand[:, 0]/comp1_byhand[:, 1]),
                                 np.log(comp1_byhand[:, 2]/comp1_byhand[:, 1])],
                                ).T
        B = 1/(1 + np.exp(alr2d_byhand[:, 0]) + np.exp(alr2d_byhand[:, 1]))
        A = B * np.exp(alr2d_byhand[:, 0])
        C = B * np.exp(alr2d_byhand[:, 1])
        alrinv2d_byhand = np.column_stack((A, B, C))
        alrinv2d_byhand = self.namespace_asarray(alrinv2d_byhand,
                                                 device=self.device)
        assert_allclose(alrinv2d_byhand, alrinv2d_method)

        # 1d-composition
        self.cdata2 = self.namespace_asarray([2, 2, 6],
                                             device = self.device)
        comp2 = closure(self.cdata2)
        alr1d_method = alr(comp2, denominator_idx=1)
        alrinv1d_method = alr_inv(alr1d_method, denominator_idx=1)

        comp2_byhand = closure(np.array([2, 2, 6]))
        alr1d_byhand = np.array([
                        np.log(comp2_byhand[0]/comp2_byhand[1]),
                        np.log(comp2_byhand[2]/comp2_byhand[1])]
                                ).T
        B = 1/(1 + np.exp(alr1d_byhand[0]) + np.exp(alr1d_byhand[1]))
        A = B * np.exp(alr1d_byhand[0])
        C = B * np.exp(alr1d_byhand[1])
        alrinv1d_byhand = np.column_stack((A, B, C))[0, :]
        alrinv1d_byhand = self.namespace_asarray(alrinv1d_byhand,
                                                 device=self.device)
        assert_allclose(alrinv1d_byhand, alrinv1d_method)

        # make sure that inplace modification is not occurring
        self.rdata1 = self.namespace_asarray([
                        [0.70710678, -0.70710678, 0., 0.],
                        [0.40824829, 0.40824829, -0.81649658, 0.],
                        [0.28867513, 0.28867513, 0.28867513, -0.8660254]],
                                             device=self.device)
        alr_inv(self.rdata1)
        rst = self.namespace_asarray([[0.70710678, -0.70710678,\
                                            0., 0.],
                                      [0.40824829, 0.40824829,\
                                            -0.81649658, 0.],
                                      [0.28867513, 0.28867513,\
                                            0.28867513, -0.8660254]],
                                     device=self.device)
        assert_allclose(self.rdata1, rst)

    def test_sbp_basis_gram_schmidt(self):
        gsbasis = _gram_schmidt_basis(5)
        gsbasis = self.namespace_asarray(gsbasis,
                                         device=self.device)
        sbp = np.array([[1, -1, 0, 0, 0],
                        [1, 1, -1, 0, 0],
                        [1, 1, 1, -1, 0],
                        [1, 1, 1, 1, -1]])
        sbpbasis = sbp_basis(sbp)
        sbpbasis = self.namespace_asarray(sbpbasis,
                                         device=self.device)
        assert_allclose(gsbasis, sbpbasis)


if __name__ == "__main__":
    main()
