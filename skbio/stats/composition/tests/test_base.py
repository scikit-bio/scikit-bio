# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
from numpy.exceptions import AxisError
import pandas as pd
from scipy.sparse import coo_matrix

from skbio import TreeNode
from skbio.stats.distance import DistanceMatrixError
from skbio.stats.composition import (
    closure, multi_replace, perturb, perturb_inv, power, inner, clr, clr_inv, rclr,
    ilr, ilr_inv, alr, alr_inv, sbp_basis, centralize, vlr, pairwise_vlr, tree_basis)
from skbio.stats.composition._base import (
    _check_composition, _check_basis, _gram_schmidt_basis)


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

        # 3-D array (tensor) of 2 x 3 x 4
        self.cdata9 = np.array([[[1, 2, 6, 1],
                                 [1, 5, 3, 1],
                                 [5, 1, 2, 2]],
                                [[2, 4, 1, 3],
                                 [3, 1, 3, 3],
                                 [4, 1, 1, 4]]])

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
        # all-zero rows
        self.bad3 = np.array([[0, 1, 2], [0, 0, 0], [3, 0, 4]])


    def test_check_composition(self):
        self.assertIsNone(_check_composition(np, self.cdata1))
        self.assertIsNone(_check_composition(np, self.cdata2))
        self.assertIsNone(_check_composition(np, self.cdata3))

        msg = "Input matrix must have a numeric data type."
        with self.assertRaises(TypeError) as cm:
            _check_composition(np, np.array(['a', 'b', 'c']))
        self.assertEqual(str(cm.exception), msg)

        msg = "Input matrix cannot have infinite or NaN values."
        with self.assertRaises(ValueError) as cm:
            _check_composition(np, np.array([1., np.nan, 2.]))
        self.assertEqual(str(cm.exception), msg)
        with self.assertRaises(ValueError) as cm:
            _check_composition(np, np.array([1., np.inf, 2.]))
        self.assertEqual(str(cm.exception), msg)

        msg = "Input matrix cannot have negative components."
        with self.assertRaises(ValueError) as cm:
            _check_composition(np, self.bad1)
        self.assertEqual(str(cm.exception), msg)

        msg = "Input matrix cannot have compositions with all zeros."
        with self.assertRaises(ValueError) as cm:
            _check_composition(np, self.bad3)
        self.assertEqual(str(cm.exception), msg)

        # all-zero composition in column not in row
        mat = np.array([[1, 5, 0, 3], [2, 0, 0, 4], [3, 8, 0, 0]])
        self.assertIsNone(_check_composition(np, mat))
        self.assertIsNone(_check_composition(np, mat, axis=1))
        with self.assertRaises(ValueError) as cm:
            _check_composition(np, mat, axis=0)
        self.assertEqual(str(cm.exception), msg)

        # single vector with a zero value
        self.assertIsNone(_check_composition(np, self.cdata4))
        self.assertIsNone(_check_composition(np, np.atleast_2d(self.cdata4)))
        with self.assertRaises(ValueError) as cm:
            _check_composition(np, self.cdata4.reshape(-1, 1))
        self.assertEqual(str(cm.exception), msg)

        # edge case: single scalar
        self.assertIsNone(_check_composition(np, np.array(5)))
        with self.assertRaises(ValueError) as cm:
            _check_composition(np, np.array(0))
        self.assertEqual(str(cm.exception), msg)

    def test_check_basis(self):
        # as range
        basis_non_orthongnal = np.array([[2, -2, 0], [2, 2, 4], [2, 2, -1]])
        basis_non_orthonormal = np.array([[2, -2, 0], [2, 2, 4], [2, 2, -1]])

        basis_unmatch_subspace_dim = np.array([[1, 0, 0]])

        # generic basis, not necessarily the basis of unit ball subspace S^2
        basis_int = np.array([[1, 0, 0], [0, 1, 0]])
        basis_real= np.array([[1.0, 0, 0], [0, 1.0, 0]])

        # action + assert
        with self.assertRaises(ValueError) as cm:
            _check_basis(np, basis_non_orthongnal, orthonormal=True)
        self.assertEqual(str(cm.exception), "Basis is not orthonormal.")

        with self.assertRaises(ValueError) as cm:
            _check_basis(np, basis_non_orthonormal, orthonormal=True)
        self.assertEqual(str(cm.exception), "Basis is not orthonormal.")

        msg = "Number of basis 1 not match to the subspace dim 2."
        with self.assertRaises(ValueError) as cm:
            _check_basis(np, basis_unmatch_subspace_dim, orthonormal=True,
                         subspace_dim=2)
        self.assertEqual(str(cm.exception), msg)

        self.assertIsNone(_check_basis(np, basis_int, orthonormal=True))
        self.assertIsNone(_check_basis(np, basis_real, orthonormal=True))

        # old test
        basis = np.array([[0.80442968, 0.19557032]])
        with self.assertRaises(ValueError) as cm:
            _check_basis(np, basis, orthonormal=True)
        self.assertEqual(str(cm.exception), "Basis is not orthonormal.")

        basis = clr(basis)
        self.assertIsNone(_check_basis(np, basis, orthonormal=True))

    def test_closure(self):
        # 2-D matrix
        mat = self.cdata1
        obs = closure(mat)
        exp = np.array([[.2, .2, .6],
                        [.4, .4, .2]])
        npt.assert_allclose(obs, exp)

        # confirm that compositions sum to 1
        npt.assert_allclose(obs.sum(axis=-1), 1.)

        # custom axis
        obs = closure(mat, axis=1)
        npt.assert_allclose(obs, exp)

        obs = closure(mat, axis=0)
        exp = np.array([[0.333, 0.333, 0.75 ],
                        [0.667, 0.667, 0.25 ]])
        npt.assert_array_equal(obs.round(3), exp)
        npt.assert_allclose(obs.sum(axis=0), 1.)

        obs = closure(mat, axis=-2)
        npt.assert_array_equal(obs.round(3), exp)

        # invalid axis
        self.assertRaises(AxisError, closure, mat, axis=3)

        # 1-D vector
        vec = self.cdata2
        obs = closure(vec)
        exp = np.array([.2, .2, .6])
        npt.assert_allclose(obs, exp)

        # make sure that inplace modification is not occurring
        self.assertIsNot(obs, vec)
        npt.assert_array_equal(vec, np.array([2, 2, 6]))

        # input is a list
        lst = self.cdata1.tolist()
        obs = closure(lst)
        exp = np.array([[.2, .2, .6],
                        [.4, .4, .2]])
        npt.assert_allclose(obs, exp)

        # input is a dataframe
        df = pd.DataFrame(self.cdata1)
        obs = closure(df)
        npt.assert_allclose(obs, exp)

        # negative value is prohibited
        msg = "Input matrix cannot have negative components."
        with self.assertRaises(ValueError) as cm:
            closure(self.bad1)
        self.assertEqual(str(cm.exception), msg)

        # zero value is allowed
        obs = closure(self.bad2)
        exp = np.array([[[0.091, 0.182, 0.273, 0.   , 0.455]]])
        npt.assert_array_equal(obs.round(3), exp)

        # all-zero composition
        msg = "Input matrix cannot have compositions with all zeros."
        with self.assertRaises(ValueError) as cm:
            closure(self.bad3)
        self.assertEqual(str(cm.exception), msg)

        # not all-zero in another axis
        obs = closure(self.bad3, axis=0)
        exp = np.array([[0.   , 1.   , 0.333],
                        [0.   , 0.   , 0.   ],
                        [1.   , 0.   , 0.667]])
        npt.assert_array_equal(obs.round(3), exp)

        # 3-D tensor
        ten = self.cdata9
        obs = closure(ten)
        exp = np.array([[[.1, .2, .6, .1],
                         [.1, .5, .3, .1],
                         [.5, .1, .2, .2]],
                        [[.2, .4, .1, .3],
                         [.3, .1, .3, .3],
                         [.4, .1, .1, .4]]])
        npt.assert_allclose(obs, exp)
        npt.assert_allclose(obs.sum(axis=-1), 1.)

        # middle axis
        obs = closure(ten, axis=1)
        exp = np.array([[[0.143, 0.25 , 0.545, 0.25 ],
                         [0.143, 0.625, 0.273, 0.25 ],
                         [0.714, 0.125, 0.182, 0.5  ]],
                        [[0.222, 0.667, 0.2  , 0.3  ],
                         [0.333, 0.167, 0.6  , 0.3  ],
                         [0.444, 0.167, 0.2  , 0.4  ]]])
        npt.assert_array_equal(obs.round(3), exp)
        npt.assert_allclose(obs.sum(axis=1), 1.)

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
        self.assertIsNotNone(perturb(
            closure(self.cdata5), self.bad1, validate=False))

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
        self.assertIsNotNone(power(self.bad1, 2, validate=False))

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
        self.assertIsNotNone(perturb_inv(
            closure(self.cdata1), self.bad1, validate=False))

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

        # just for test
        self.assertIsNotNone(inner(self.cdata7, self.cdata7, validate=False))

        # dimension not match
        with self.assertRaises(ValueError):
            inner(self.cdata1, self.cdata8)

        # invalid compositions
        with self.assertRaises(ValueError):
            inner(self.bad1, self.bad1)

        # make sure that inplace modification is not occurring
        inner(self.cdata1, self.cdata1)
        npt.assert_allclose(self.cdata1,
                            np.array([[2, 2, 6],
                                      [4, 4, 2]]))

    def test_multi_replace(self):
        obs = multi_replace(closure(self.cdata3))
        exp = np.array([[0.087273, 0.174545, 0.261818, 0.04, 0.436364],
                        [0.092, 0.04, 0.04, 0.368, 0.46],
                        [0.066667, 0.133333, 0.2, 0.266667, 0.333333]])
        npt.assert_allclose(obs, exp, rtol=1e-5, atol=1e-5)

        obs = multi_replace(closure(self.cdata6))
        npt.assert_allclose(obs, exp, rtol=1e-5, atol=1e-5)

        obs = multi_replace(closure(self.cdata4))
        exp = np.array([0.087273, 0.174545, 0.261818, 0.04, 0.436364])
        npt.assert_allclose(obs, exp, rtol=1e-5, atol=1e-5)

        # manually specify auto-calculated delta
        obs = multi_replace(closure(self.cdata4), delta=0.04)
        npt.assert_allclose(obs, exp, rtol=1e-5, atol=1e-5)

        # non-default delta
        obs = multi_replace(closure(self.cdata4), delta=0.05)
        exp = np.array([0.086364, 0.172727, 0.259091, 0.05, 0.431818])
        npt.assert_allclose(obs, exp, rtol=1e-5, atol=1e-5)

        msg = "Consider using a smaller `delta`."
        with self.assertRaisesRegex(ValueError, msg):
            obs = multi_replace(closure(self.cdata4), delta=2.0)

        self.assertRaises(ValueError, multi_replace, self.bad1)

        # make sure that inplace modification is not occurring
        multi_replace(self.cdata4)
        npt.assert_allclose(self.cdata4, np.array([1, 2, 3, 0, 5]))

    def test_clr(self):
        # 2-D matrix
        mat = self.cdata1
        cmat = clr(mat)

        # calculation by hand
        A = np.array([.2, .2, .6])
        B = np.array([.4, .4, .2])
        exp = [np.log(A / np.exp(np.log(A).mean())),
               np.log(B / np.exp(np.log(B).mean()))]
        npt.assert_allclose(cmat, exp)

        # results are 0-centered
        npt.assert_allclose(cmat.sum(axis=1), 0, atol=1e-8)

        # closure has no effect on result
        cmat = clr(closure(mat))
        npt.assert_allclose(cmat, exp)

        # CLR is not sensitive to scale
        cmat = clr(mat * 100)
        npt.assert_allclose(cmat, exp)

        # custom axis
        cmat = clr(mat, axis=0)
        exp = np.vstack([clr(x) for x in mat.T]).T
        npt.assert_allclose(cmat, exp)

        # 1-D vector
        cmat = clr(closure(self.cdata2))
        A = np.array([.2, .2, .6])
        exp = np.log(A / np.exp(np.log(A).mean()))
        npt.assert_allclose(cmat, exp)

        # invalid input matrix
        msg = "Input matrix cannot have negative or zero components."

        # negative value
        with self.assertRaises(ValueError) as cm:
            clr(self.bad1)
        self.assertEqual(str(cm.exception), msg)

        # zero value
        with self.assertRaises(ValueError) as cm:
            clr(self.bad2)
        self.assertEqual(str(cm.exception), msg)

        # all-zero composition
        with self.assertRaises(ValueError) as cm:
            clr(self.bad3)
        self.assertEqual(str(cm.exception), msg)

        # make sure that inplace modification is not occurring
        clr(self.cdata2)
        npt.assert_allclose(self.cdata2, np.array([2, 2, 6]))

        # 3-D tensor as input
        ten = self.cdata9
        obs = clr(ten)
        exp = np.array([[[-0.62123,  0.07192,  1.17053, -0.62123],
                         [-0.67701,  0.93243,  0.4216 , -0.67701],
                         [ 0.8605 , -0.74893, -0.05579, -0.05579]],
                        [[-0.10137,  0.59178, -0.79451,  0.3041 ],
                         [ 0.27465, -0.82396,  0.27465,  0.27465],
                         [ 0.69315, -0.69315, -0.69315,  0.69315]]])
        npt.assert_array_equal(obs.round(5), exp)

        # The result should be identical to applying clr to each matrix separately,
        for obs2d, mat in zip(obs, ten):
            npt.assert_allclose(obs2d, clr(mat))

        # ...and identical to applying clr to each row separately.
        for obs2d, mat in zip(obs, ten):
            for obs1d, vec in zip(obs2d, mat):
                npt.assert_allclose(obs1d, clr(vec))

        # middle axis
        obs = clr(ten, axis=1)
        exp = np.vstack([np.expand_dims(clr(mat, axis=0), axis=0) for mat in ten])
        npt.assert_allclose(obs, exp)

    def test_clr_inv(self):
        mat = self.rdata1.copy()
        obs = clr_inv(mat)
        npt.assert_allclose(obs, self.ortho1)

        # check that clr_inv is the inverse of clr (if mat is already closure)
        npt.assert_allclose(clr(obs), self.rdata1, rtol=1e-4, atol=1e-5)

        # make sure that inplace modification is not occurring
        self.assertIsNot(obs, mat)
        npt.assert_allclose(mat, self.rdata1)

        # custom axis
        mat = self.cdata1
        obs = clr_inv(clr(mat, axis=0), axis=0)
        exp = closure(mat, axis=0)
        npt.assert_allclose(obs, exp)

        # 3-D tensor as input (see `test_clr` above)
        arr3d = self.cdata9
        obs = clr_inv(clr(arr3d))
        exp = closure(arr3d)
        npt.assert_allclose(obs, exp)

        # input not centered
        with self.assertWarns(UserWarning):
            clr_inv(self.cdata1)

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

        # make sure that inplace modification is not occurring
        centralize(self.cdata1)
        npt.assert_allclose(self.cdata1,
                            np.array([[2, 2, 6],
                                      [4, 4, 2]]))

    def test_ilr(self):
        mat = closure(self.cdata7)
        exp = np.array([0.70710678, 0.40824829])
        npt.assert_allclose(ilr(mat), exp)

        # Should give same result as inner
        npt.assert_allclose(ilr(self.ortho1), np.identity(3),
                            rtol=1e-04, atol=1e-06)

        # no check
        npt.assert_allclose(ilr(mat, validate=False), exp)

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
        obs = ilr(table, basis=basis)
        exp = np.array([[np.log(1/10)*np.sqrt(1/2)],
                        [np.log(1.14141414 / 9.90909091)*np.sqrt(1/2)],
                        [np.log(1.28282828 / 9.81818182)*np.sqrt(1/2)],
                        [np.log(1.42424242 / 9.72727273)*np.sqrt(1/2)],
                        [np.log(1.56565657 / 9.63636364)*np.sqrt(1/2)]])

        npt.assert_allclose(obs, exp)

        obs = ilr(table, basis=basis, validate=False)
        npt.assert_allclose(obs, exp)

    def test_ilr_errors(self):
        msg = "Input matrix cannot have negative or zero components."
        with self.assertRaises(ValueError) as cm:
            ilr(self.bad1)
        self.assertEqual(str(cm.exception), msg)

        # msg = "Input matrix can only have two dimensions or less."
        # with self.assertRaises(ValueError) as cm:
        #     ilr(np.array([[[1, 2, 3]]]))
        # self.assertEqual(str(cm.exception), msg)

        basis = np.array([[0.80442968, 0.19557032]])
        msg = "Number of basis 1 not match to the subspace dim 2."
        with self.assertRaises(ValueError) as cm:
            ilr(self.cdata1, basis=basis)
        self.assertEqual(str(cm.exception), msg)

        basis = np.squeeze(clr(basis))
        msg = "Basis needs to be a 2-D matrix, not a 1-D matrix."
        with self.assertRaises(ValueError) as cm:
            ilr(self.cdata1, basis=basis)
        self.assertEqual(str(cm.exception), msg)

    def test_ilr_inv(self):
        mat = closure(self.cdata7)
        npt.assert_allclose(ilr_inv(ilr(mat)), mat)

        npt.assert_allclose(ilr_inv(np.identity(3)), self.ortho1,
                            rtol=1e-04, atol=1e-06)

        # no check
        npt.assert_allclose(ilr_inv(ilr(mat), validate=False), mat)

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
        alr2d_method = alr(comp1, ref_idx=1)
        npt.assert_allclose(alr2d_byhand, alr2d_method)

        # 1d-composition
        comp2 = closure(self.cdata2)
        alr1d_byhand = np.array([np.log(comp2[0]/comp2[1]),
                                 np.log(comp2[2]/comp2[1])]).T
        alr1d_method = alr(comp2, ref_idx=1)
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
        alr2d_method = alr(comp1, ref_idx=1)
        B = 1/(1 + np.exp(alr2d_byhand[:, 0]) + np.exp(alr2d_byhand[:, 1]))
        A = B * np.exp(alr2d_byhand[:, 0])
        C = B * np.exp(alr2d_byhand[:, 1])
        alrinv2d_byhand = np.column_stack((A, B, C))
        alrinv2d_method = alr_inv(alr2d_method, ref_idx=1)
        npt.assert_allclose(alrinv2d_byhand, alrinv2d_method)

        # 1d-composition
        comp2 = closure(self.cdata2)
        alr1d_byhand = np.array([np.log(comp2[0]/comp2[1]),
                                 np.log(comp2[2]/comp2[1])]).T
        alr1d_method = alr(comp2, ref_idx=1)
        B = 1/(1 + np.exp(alr1d_byhand[0]) + np.exp(alr1d_byhand[1]))
        A = B * np.exp(alr1d_byhand[0])
        C = B * np.exp(alr1d_byhand[1])
        alrinv1d_byhand = np.column_stack((A, B, C))[0, :]
        alrinv1d_method = alr_inv(alr1d_method, ref_idx=1)
        npt.assert_allclose(alrinv1d_byhand, alrinv1d_method)

        # make sure that inplace modification is not occurring
        alr_inv(self.rdata1)
        npt.assert_allclose(self.rdata1,
                            np.array([[0.70710678, -0.70710678, 0., 0.],
                                      [0.40824829, 0.40824829,
                                       -0.81649658, 0.],
                                      [0.28867513, 0.28867513,
                                       0.28867513, -0.8660254]]))

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
        dism = pairwise_vlr(self.mat, ids=None, ddof=1, robust=False, validate=False)
        output = dism.data.sum() / 2
        self.assertAlmostEqual(output, 0.2857382286903922)


class TestRclr(TestCase):
    """Tests for the robust centered log-ratio transformation."""

    def test_basic_rclr(self):
        """Test basic rclr transformation."""
        # Simple matrix with no zeros
        mat = np.array([[1, 2, 3], [4, 5, 6]])
        result = rclr(mat)

        # For each row, the mean of transformed values should be ~0
        # (only over observed values)
        for i in range(mat.shape[0]):
            observed = ~np.isnan(result[i])
            npt.assert_almost_equal(result[i, observed].mean(), 0.0)

    def test_rclr_with_zeros(self):
        """Test rclr handles zeros by producing NaN."""
        mat = np.array([[1, 0, 3], [4, 5, 0]])
        result = rclr(mat)

        # Zeros should become NaN
        self.assertTrue(np.isnan(result[0, 1]))
        self.assertTrue(np.isnan(result[1, 2]))

        # Non-zero positions should not be NaN
        self.assertFalse(np.isnan(result[0, 0]))
        self.assertFalse(np.isnan(result[0, 2]))

    def test_rclr_centering(self):
        """Test that rclr centers each row correctly."""
        mat = np.array([[1, 2, 0, 4], [0, 3, 3, 0], [2, 2, 2, 2]])
        result = rclr(mat)

        # For rows with observed values, mean should be 0
        for i in range(mat.shape[0]):
            observed = ~np.isnan(result[i])
            if np.any(observed):
                npt.assert_almost_equal(result[i, observed].mean(), 0.0)

    def test_rclr_negative_values_error(self):
        """Test that negative values raise an error."""
        mat = np.array([[1, -2, 3], [4, 5, 6]])

        with self.assertRaises(ValueError) as context:
            rclr(mat)

        self.assertIn("negative", str(context.exception))

    def test_rclr_inf_error(self):
        """Test that infinite values raise an error."""
        mat = np.array([[1, np.inf, 3], [4, 5, 6]])

        with self.assertRaises(ValueError) as context:
            rclr(mat)

        self.assertIn("infinite", str(context.exception))

    def test_rclr_nan_preserved(self):
        """Test that NaN values (missing entries) are preserved."""
        mat = np.array([[1, np.nan, 3], [4, 5, 6]])
        result = rclr(mat)

        # NaN should remain NaN
        self.assertTrue(np.isnan(result[0, 1]))

        # Non-NaN values should be transformed
        self.assertFalse(np.isnan(result[0, 0]))
        self.assertFalse(np.isnan(result[0, 2]))

    def test_rclr_1d_input(self):
        """Test that 1D input works correctly."""
        vec = np.array([1, 2, 3])
        result = rclr(vec)

        # Should handle 1D gracefully
        self.assertEqual(result.shape, vec.shape)
        npt.assert_almost_equal(np.nanmean(result), 0.0)

    def test_rclr_uniform_row(self):
        """Test rclr on uniform row (all same values)."""
        mat = np.array([[2, 2, 2, 2]])
        result = rclr(mat)

        # Uniform row should have all zeros (log-ratio of equal values)
        npt.assert_almost_equal(result[0], [0.0, 0.0, 0.0, 0.0])

    def test_rclr_preserves_ratios(self):
        """Test that rclr preserves log-ratios between features."""
        mat = np.array([[1, 2, 4]])
        result = rclr(mat)

        # log(2) - log(1) = log(2)
        expected_ratio = np.log(2)
        observed_ratio = result[0, 1] - result[0, 0]
        npt.assert_almost_equal(observed_ratio, expected_ratio)

    def test_rclr_3d_tensor(self):
        """Test rclr on 3D tensor."""
        tensor = np.array([
            [[1, 2, 3], [4, 5, 6]],
            [[7, 8, 9], [10, 11, 12]]
        ])
        result = rclr(tensor)

        # Shape should be preserved
        self.assertEqual(result.shape, tensor.shape)

        # Each sample's mean should be approximately 0
        result_2d = result.reshape(-1, tensor.shape[-1])
        for i in range(result_2d.shape[0]):
            observed = ~np.isnan(result_2d[i])
            if np.any(observed):
                npt.assert_almost_equal(result_2d[i, observed].mean(), 0.0,
                                        decimal=5)

    def test_rclr_3d_with_zeros(self):
        """Test rclr handles zeros in 3D tensor."""
        tensor = np.array([
            [[1, 0, 3], [0, 5, 6]],
            [[7, 8, 0], [10, 0, 12]]
        ])
        result = rclr(tensor)

        # Zeros should become NaN
        self.assertTrue(np.isnan(result[0, 0, 1]))
        self.assertTrue(np.isnan(result[0, 1, 0]))
        self.assertTrue(np.isnan(result[1, 0, 2]))
        self.assertTrue(np.isnan(result[1, 1, 1]))

    def test_rclr_various_shapes(self):
        """Test that various tensor shapes work correctly."""
        shapes = [(2, 3, 4), (5, 2, 6), (3, 3, 3)]

        for shape in shapes:
            tensor = np.random.rand(*shape) + 0.1  # Avoid zeros
            result = rclr(tensor)
            self.assertEqual(result.shape, shape)

    def test_rclr_dense_equals_clr(self):
        """Test that rclr equals clr on dense data without zeros."""
        # Dense count data with no zeros
        count_data = np.array([[2, 2, 6], [4, 4, 2]])

        # Apply rclr
        rclr_result = rclr(count_data)

        # Apply closure then clr
        closed_data = closure(count_data)
        clr_result = clr(closed_data)

        # Results should be identical for dense data
        npt.assert_allclose(rclr_result, clr_result, rtol=1e-10)

    def test_rclr_sparse_expected_values(self):
        """Test rclr on sparse data with known expected values."""
        # Sparse count data with zeros
        count_data = np.array([[3, 3, 0], [0, 4, 2]])

        # Expected values (zeros become NaN, observed values centered)
        expected = np.array([[0.0, 0.0, np.nan],
                            [np.nan, 0.34657359, -0.34657359]])

        result = rclr(count_data)

        # Check non-NaN values match expected
        npt.assert_allclose(result[~np.isnan(result)],
                           expected[~np.isnan(expected)],
                           rtol=1e-5)

        # Check NaN positions match
        npt.assert_array_equal(np.isnan(result), np.isnan(expected))

    def test_rclr_axis_parameter(self):
        """Test rclr with different axis values."""
        mat = np.array([[1, 2, 3], [4, 5, 6]])

        # Default axis=-1 (features on last axis)
        result_default = rclr(mat)
        result_last = rclr(mat, axis=-1)
        npt.assert_array_equal(result_default, result_last)

        # axis=0 (features on first axis)
        result_axis0 = rclr(mat, axis=0)
        self.assertEqual(result_axis0.shape, mat.shape)

        # Each column should be centered
        for j in range(mat.shape[1]):
            observed = ~np.isnan(result_axis0[:, j])
            if np.any(observed):
                npt.assert_almost_equal(result_axis0[observed, j].mean(), 0.0)

    def test_rclr_log_ratio_preservation(self):
        """Test that log-ratios are preserved exactly."""
        # Data where we can compute expected values analytically
        data = np.array([[1, 2, 4, 8]])

        result = rclr(data)

        # Log-ratio between adjacent pairs should be log(2)
        expected_ratio = np.log(2)
        for i in range(3):
            observed_ratio = result[0, i + 1] - result[0, i]
            npt.assert_almost_equal(observed_ratio, expected_ratio)

    def test_rclr_geometric_mean_centering(self):
        """Test that each row is centered at geometric mean."""
        np.random.seed(42)
        data = np.random.rand(5, 10) * 100 + 1  # No zeros

        result = rclr(data)

        # Each row should have mean of 0 (centered at geometric mean)
        for i in range(data.shape[0]):
            row_mean = np.nanmean(result[i])
            npt.assert_almost_equal(row_mean, 0.0, decimal=10)

    def test_rclr_validate_false(self):
        """Test that validation can be disabled."""
        # This should work with validate=False even with potential edge cases
        mat = np.array([[1, 2, 3], [4, 5, 6]])
        result = rclr(mat, validate=False)
        self.assertEqual(result.shape, mat.shape)


if __name__ == "__main__":
    main()
