# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from copy import deepcopy

import numpy as np
import numpy.testing as npt
from numpy.exceptions import AxisError
from numpy.random import normal
import pandas as pd
import pandas.testing as pdt
from scipy.sparse import coo_matrix
from scipy.stats import f_oneway, ConstantInputWarning

from skbio import TreeNode
from skbio.util import assert_data_frame_almost_equal
from skbio.stats.distance import DistanceMatrixError
from skbio.stats.composition import (
    _check_composition, _check_basis, _check_grouping, _check_trt_ref_groups,
    _check_metadata, _type_cast_to_float,
    closure, multi_replace, perturb, perturb_inv, power, inner, clr, clr_inv, ilr,
    ilr_inv, alr, alr_inv, sbp_basis, _gram_schmidt_basis, centralize, _check_sig_test,
    _check_p_adjust, ancom, vlr, pairwise_vlr, tree_basis, dirmult_ttest, dirmult_lme)


def assert_coo_allclose(res, exp, rtol=1e-7, atol=1e-7):
    res_data = np.vstack((res.row, res.col, res.data)).T
    exp_data = np.vstack((exp.row, exp.col, exp.data)).T

    # sort by row and col
    res_data = res_data[res_data[:, 1].argsort()]
    res_data = res_data[res_data[:, 0].argsort()]
    exp_data = exp_data[exp_data[:, 1].argsort()]
    exp_data = exp_data[exp_data[:, 0].argsort()]
    npt.assert_allclose(res_data, exp_data, rtol=rtol, atol=atol)


class MiscTests(TestCase):
    def test_check_grouping(self):
        matrix = np.array([[1, 2], [3, 4], [5, 6]])
        grouping = [0, 0, 1]
        obs = _check_grouping(grouping, matrix)
        npt.assert_array_equal(obs[0], [0, 1])
        npt.assert_array_equal(obs[1], grouping)

        grouping = [5, 2, 5]
        obs = _check_grouping(grouping, matrix)
        npt.assert_array_equal(obs[0], [2, 5])
        npt.assert_array_equal(obs[1], [1, 0, 1])

        grouping = ["b", "b", "a"]
        obs = _check_grouping(grouping, matrix)
        npt.assert_array_equal(obs[0], ["a", "b"])
        npt.assert_array_equal(obs[1], [1, 1, 0])

        grouping = pd.Series(grouping)
        obs = _check_grouping(grouping, matrix)
        npt.assert_array_equal(obs[0], ["a", "b"])
        npt.assert_array_equal(obs[1], [1, 1, 0])

        msg = "`table` contains sample IDs that are absent in `grouping`."
        samples = ["x", "y", "z"]
        with self.assertRaises(ValueError) as cm:
            _check_grouping(grouping, matrix, samples=samples)
        self.assertEqual(str(cm.exception), msg)

        grouping.index = ["x", "y", "z"]
        obs = _check_grouping(grouping, matrix, samples=samples)
        npt.assert_array_equal(obs[0], ["a", "b"])
        npt.assert_array_equal(obs[1], [1, 1, 0])

        msg = "Sample counts in `table` and `grouping` are not consistent."
        grouping = ["b", "c", "a", "b"]
        with self.assertRaises(ValueError) as cm:
            _check_grouping(grouping, matrix)
        self.assertEqual(str(cm.exception), msg)

        grouping = pd.Series(grouping)
        with self.assertRaises(ValueError) as cm:
            _check_grouping(grouping, matrix)
        self.assertEqual(str(cm.exception), msg)

        grouping.index = ["y", "z", "x", "w"]
        samples = ["x", "y", "z"]
        obs = _check_grouping(grouping, matrix, samples=samples)
        npt.assert_array_equal(obs[0], ["a", "b", "c"])
        npt.assert_array_equal(obs[1], [0, 1, 2])

        msg = "Cannot handle missing values in `grouping`."
        grouping = np.array([1., np.nan, 3.])
        with self.assertRaises(ValueError) as cm:
            _check_grouping(grouping, matrix)
        self.assertEqual(str(cm.exception), msg)

        grouping = [1, None, 3]
        with self.assertRaises(ValueError) as cm:
            _check_grouping(grouping, matrix)
        self.assertEqual(str(cm.exception), msg)

        msg = "`grouping` must be convertible to a 1-D vector."
        grouping = np.array([["a", "b"], ["c", "g"], ["e", "d"]])
        with self.assertRaises(ValueError) as cm:
            _check_grouping(grouping, matrix)
        self.assertEqual(str(cm.exception), msg)

        grouping = 123
        with self.assertRaises(ValueError) as cm:
            _check_grouping(grouping, matrix)
        self.assertEqual(str(cm.exception), msg)

    def test_check_trt_ref_groups(self):
        # two groups
        grouping = ["B", "A", "B", "B", "A", "A", "A", "B"]
        groups, labels = _check_grouping(grouping, np.empty((8, 1)))

        obs = _check_trt_ref_groups("A", "B", groups, labels)
        npt.assert_array_equal(obs[0], [1, 4, 5, 6])
        npt.assert_array_equal(obs[1], [0, 2, 3, 7])

        obs = _check_trt_ref_groups("B", "A", groups, labels)
        npt.assert_array_equal(obs[0], [0, 2, 3, 7])
        npt.assert_array_equal(obs[1], [1, 4, 5, 6])

        # default groups
        obs = _check_trt_ref_groups(None, None, groups, labels)
        npt.assert_array_equal(obs[0], [1, 4, 5, 6])
        npt.assert_array_equal(obs[1], [0, 2, 3, 7])

        obs = _check_trt_ref_groups("A", None, groups, labels)
        npt.assert_array_equal(obs[0], [1, 4, 5, 6])
        npt.assert_array_equal(obs[1], [0, 2, 3, 7])

        obs = _check_trt_ref_groups("B", None, groups, labels)
        npt.assert_array_equal(obs[0], [0, 2, 3, 7])
        npt.assert_array_equal(obs[1], [1, 4, 5, 6])

        msg = "Treatment group C is not found in grouping."
        with self.assertRaises(ValueError) as cm:
            _check_trt_ref_groups("C", None, groups, labels)
        self.assertEqual(str(cm.exception), msg)

        msg = "Reference group D is not found in grouping."
        with self.assertRaises(ValueError) as cm:
            _check_trt_ref_groups("A", "D", groups, labels)
        self.assertEqual(str(cm.exception), msg)

        msg = "Treatment and reference groups must not be identical."
        with self.assertRaises(ValueError) as cm:
            _check_trt_ref_groups("A", "A", groups, labels)
        self.assertEqual(str(cm.exception), msg)
        with self.assertRaises(ValueError) as cm:
            _check_trt_ref_groups(None, "A", groups, labels)
        self.assertEqual(str(cm.exception), msg)

        # one group
        grouping = ["A", "A", "A", "A"]
        groups, labels = _check_grouping(grouping, np.empty((4, 1)))

        msg = "There must be at least two groups in grouping."
        with self.assertRaises(ValueError) as cm:
            _check_trt_ref_groups(None, None, groups, labels)
        self.assertEqual(str(cm.exception), msg)

        # three groups
        grouping = ["A", "C", "B", "B", "C", "A", "A", "C"]
        groups, labels = _check_grouping(grouping, np.empty((8, 1)))

        obs = _check_trt_ref_groups("A", "B", groups, labels)
        npt.assert_array_equal(obs[0], [0, 5, 6])
        npt.assert_array_equal(obs[1], [2, 3])

        obs = _check_trt_ref_groups("C", "A", groups, labels)
        npt.assert_array_equal(obs[0], [1, 4, 7])
        npt.assert_array_equal(obs[1], [0, 5, 6])

        obs = _check_trt_ref_groups("B", None, groups, labels)
        npt.assert_array_equal(obs[0], [2, 3])
        npt.assert_array_equal(obs[1], [0, 1, 4, 5, 6, 7])

        obs = _check_trt_ref_groups(None, "C", groups, labels)
        npt.assert_array_equal(obs[0], [0, 5, 6])
        npt.assert_array_equal(obs[1], [1, 4, 7])

        obs = _check_trt_ref_groups(None, None, groups, labels)
        npt.assert_array_equal(obs[0], [0, 5, 6])
        npt.assert_array_equal(obs[1], [1, 2, 3, 4, 7])

    def test_check_metadata(self):
        mat = np.empty(12).reshape(3, 4)
        df = pd.DataFrame([("Alice", 20, 28.0),
                           ("Bob",   32, 33.0),
                           ("Carol", 25, 26.5)],
                          columns=["name", "age", "bmi"])
        obs = _check_metadata(df, mat)
        self.assertIs(obs, df)

        lst = [("Alice", 20, 28.0),
               ("Bob",   32, 33.0),
               ("Carol", 25, 26.5)]
        obs = _check_metadata(lst, mat)
        self.assertIsInstance(obs, pd.DataFrame)

        dic = {"name": ["Alice", "Bob", "Carol"],
               "age":  [20, 32, 25],
               "bmi":  [28.0, 33.0, 26.5]}
        obs = _check_metadata(dic, mat)
        self.assertIsInstance(obs, pd.DataFrame)

        arr = np.array([("Alice", 20, 28.0),
                        ("Bob",   32, 33.0),
                        ("Carol", 25, 26.5)],
                       dtype=[("name", "U10"),
                              ("age", "i4"),
                              ("bmi", "f4")])
        obs = _check_metadata(arr, mat)
        self.assertIsInstance(obs, pd.DataFrame)

        msg = "Metadata must be a pandas DataFrame"
        with self.assertRaisesRegex(TypeError, msg):
            _check_metadata(42, mat)
        with self.assertRaisesRegex(TypeError, msg):
            _check_metadata("hello", mat)

        msg = "Sample counts in table and metadata are not consistent."
        with self.assertRaises(ValueError) as cm:
            _check_metadata(df, mat.reshape(4, 3))
        self.assertEqual(str(cm.exception), msg)

        df.index = ["a", "b", "c"]
        obs = _check_metadata(df, mat)
        self.assertIs(obs, df)

        # check sample IDs
        samples = ["a", "b", "c"]
        obs = _check_metadata(df, mat, samples=samples)
        self.assertIs(obs, df)

        # reorder samples
        samples = ["b", "c", "a"]
        obs = _check_metadata(df, mat, samples=samples)
        self.assertIsNot(obs, df)
        pdt.assert_index_equal(obs.index, pd.Index(samples))

        # filter and reorder samples
        samples = ["c", "b"]
        obs = _check_metadata(df, mat, samples=samples)
        self.assertIsNot(obs, df)
        pdt.assert_index_equal(obs.index, pd.Index(samples))

        msg = "Metadata contains sample IDs that are absent in the table."
        samples = ["a", "b", "x"]
        with self.assertRaises(ValueError) as cm:
            _check_metadata(df, mat, samples=samples)
        self.assertEqual(str(cm.exception), msg)

        msg = "Cannot handle missing values in metadata."
        df = pd.DataFrame(np.array([1.0, np.nan, 2.0]).reshape(3, -1))
        with self.assertRaises(ValueError) as cm:
            _check_metadata(df, mat)
        self.assertEqual(str(cm.exception), msg)
        df = pd.DataFrame(np.array([1.0, None, 2.0]).reshape(3, -1))
        with self.assertRaises(ValueError) as cm:
            _check_metadata(df, mat)
        self.assertEqual(str(cm.exception), msg)

    def test_type_cast_to_float(self):
        df = pd.DataFrame([("Alice", 20, 28.0),
                           ("Bob",   32, 33.0),
                           ("Carol", 25, 26.5)],
                          columns=["name", "age", "bmi"])
        obs = _type_cast_to_float(df)
        self.assertIsInstance(obs, pd.DataFrame)
        self.assertIsNot(obs, df)
        self.assertEqual(obs["name"].dtype, np.object_)
        self.assertEqual(obs["age"].dtype, np.float64)
        self.assertEqual(obs["bmi"].dtype, np.float64)

    def test_check_sig_test(self):
        from scipy.stats import ttest_ind, mannwhitneyu, f_oneway, kruskal

        obs = _check_sig_test(ttest_ind)
        self.assertIs(obs, ttest_ind)

        obs = _check_sig_test("ttest_ind")
        self.assertIs(obs, ttest_ind)

        obs = _check_sig_test(f_oneway)
        self.assertIs(obs, f_oneway)

        obs = _check_sig_test("f_oneway")
        self.assertIs(obs, f_oneway)

        msg = 'Function "not_a_test" does not exist under scipy.stats.'
        with self.assertRaises(ValueError) as cm:
            _check_sig_test("not_a_test")
        self.assertEqual(str(cm.exception), msg)

        msg = "`sig_test` must be a function or a string."
        with self.assertRaises(TypeError) as cm:
            _check_sig_test(123)
        self.assertEqual(str(cm.exception), msg)

        msg = ('"mannwhitneyu" is a two-way statistical test whereas 3 sample '
               "groups were provided.")
        with self.assertRaises(ValueError) as cm:
            _check_sig_test(mannwhitneyu, n_groups=3)
        self.assertEqual(str(cm.exception), msg)

        obs = _check_sig_test(mannwhitneyu, n_groups=2)
        obs = _check_sig_test(kruskal, n_groups=5)

    def test_check_p_adjust(self):
        p = [0.005, 0.011, 0.02, 0.04, 0.13]
        obs = _check_p_adjust("holm-bonferroni")(p)
        exp = p * np.arange(1, 6)[::-1]
        for a, b in zip(obs, exp):
            self.assertAlmostEqual(a, b)

        p = [0.005, 0.011, 0.02, 0.04, 0.13]
        obs = _check_p_adjust("benjamini-hochberg")(p)
        exp = [0.025, 0.0275, 0.03333333, 0.05, 0.13]
        for a, b in zip(obs, exp):
            self.assertAlmostEqual(a, b)


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
        #asrange
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
        original_table = deepcopy(test_table)
        test_cats = pd.Series(self.cats1)
        original_cats = deepcopy(test_cats)
        result = ancom(test_table, test_cats, p_adjust=None)
        # Test to make sure that the input table hasn't be altered
        assert_data_frame_almost_equal(original_table, test_table)
        # Test to make sure that the input table hasn't be altered
        pdt.assert_series_equal(original_cats, test_cats)
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Signif': np.array([True, True, False, False, False, False, False],
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

        percentiles = np.array([0.0, 25.0, 50.0, 75.0, 100.0])
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
        original_table = deepcopy(test_table)
        test_cats = pd.Series(self.cats1)
        original_cats = deepcopy(test_cats)
        result = ancom(test_table, test_cats, p_adjust=None)
        # Test to make sure that the input table hasn't be altered
        assert_data_frame_almost_equal(original_table, test_table)
        # Test to make sure that the input table hasn't be altered
        pdt.assert_series_equal(original_cats, test_cats)
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Signif': np.array([True, True, False, False, False, False, False],
                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_multiple_groups(self):
        test_table = pd.DataFrame(self.table4)
        original_table = deepcopy(test_table)
        test_cats = pd.Series(self.cats4)
        original_cats = deepcopy(test_cats)
        result = ancom(test_table, test_cats)
        # Test to make sure that the input table hasn't be altered
        assert_data_frame_almost_equal(original_table, test_table)
        # Test to make sure that the input table hasn't be altered
        pdt.assert_series_equal(original_cats, test_cats)
        exp = pd.DataFrame(
            {'W': np.array([8, 7, 3, 3, 7, 3, 3, 3, 3]),
             'Signif': np.array([True, True, False, False, True, False, False, False,
                                 False], dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_noncontiguous(self):
        result = ancom(self.table5, self.cats5, p_adjust=None)
        exp = pd.DataFrame(
            {'W': np.array([6, 2, 2, 2, 2, 6, 2]),
             'Signif': np.array([True, False, False, False, False, True, False],
                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_unbalanced(self):
        result = ancom(self.table6, self.cats6, p_adjust=None)
        exp = pd.DataFrame(
            {'W': np.array([5, 3, 3, 2, 2, 5, 2]),
             'Signif': np.array([True, False, False, False, False, True, False],
                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_letter_categories(self):
        result = ancom(self.table7, self.cats7, p_adjust=None)
        exp = pd.DataFrame(
            {'W': np.array([5, 3, 3, 2, 2, 5, 2]),
             'Signif': np.array([True, False, False, False, False, True, False],
                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_sig_test_none(self):
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Signif': np.array([True, True, False, False, False, False, False],
                                dtype=bool)})
        result = ancom(self.table1, self.cats1, sig_test=None)
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_sig_test_callable(self):
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Signif': np.array([True, True, False, False, False, False, False],
                                dtype=bool)})
        result = ancom(self.table1, self.cats1, sig_test=f_oneway)
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_multiple_comparisons(self):
        exp = pd.DataFrame(
            {'W': np.array([0] * 7),
             'Signif': np.array([False] * 7, dtype=bool)})
        for method in 'holm', 'bh':
            result = ancom(self.table1, self.cats1, p_adjust=method,
                           sig_test='mannwhitneyu')
            assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_alternative_test(self):
        result = ancom(self.table1, self.cats1, p_adjust=None,
                       sig_test="ttest_ind")
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Signif': np.array([True,  True, False, False, False, False, False],
                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_incorrect_test(self):
        with self.assertRaises(ValueError) as cm:
            ancom(self.table1, self.cats1, sig_test="not_a_test")
        msg = 'Function "not_a_test" does not exist under scipy.stats.'
        self.assertEqual(str(cm.exception), msg)

    def test_ancom_normal_data(self):
        result = ancom(self.table2, self.cats2, p_adjust=None,
                       sig_test="ttest_ind")
        exp = pd.DataFrame(
            {'W': np.array([8, 8, 3, 3, 8, 3, 3, 3, 3]),
             'Signif': np.array([True, True, False, False, True, False, False,
                                 False, False], dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_basic_counts_swapped(self):
        result = ancom(self.table8, self.cats8)
        exp = pd.DataFrame(
            {'W': np.array([5, 5, 2, 2, 2, 2, 2]),
             'Signif': np.array([True, True, False, False, False, False, False],
                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_no_signal(self):
        with self.assertWarns(ConstantInputWarning):
            result = ancom(self.table3, self.cats3, p_adjust=None)
        exp = pd.DataFrame(
            {'W': np.array([0]*7),
             'Signif': np.array([False]*7, dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_tau(self):
        exp1 = pd.DataFrame(
            {'W': np.array([8, 7, 3, 3, 7, 3, 3, 3, 3]),
             'Signif': np.array([True, False, False, False, False, False, False,
                                 False, False], dtype=bool)})
        exp2 = pd.DataFrame(
            {'W': np.array([17, 17, 5, 6, 16, 5, 7, 5,
                            4, 5, 8, 4, 5, 16, 5, 11, 4, 6]),
             'Signif': np.array([True, True, False, False,
                                                 True, False, False, False,
                                                 False, False, False, False,
                                                 False, True, False, False,
                                                 False, False],  dtype=bool)})
        exp3 = pd.DataFrame(
            {'W': np.array([16, 16, 17, 10, 17, 16, 16,
                            15, 15, 15, 13, 10, 10, 10,
                            9, 9, 9, 9]),
             'Signif': np.array([True, True, True, False,
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
             'Signif': np.array([True, True, False, False, False, False, False],
                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

    def test_ancom_alpha(self):
        result = ancom(self.table1, self.cats1, p_adjust=None, alpha=0.5)
        exp = pd.DataFrame(
            {'W': np.array([6, 6, 4, 5, 5, 4, 2]),
             'Signif': np.array([True, True, False, True, True, False, False],
                                dtype=bool)})
        assert_data_frame_almost_equal(result[0], exp)

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
        msg = ('"ttest_ind" is a two-way statistical test whereas 3 sample '
               "groups were provided.")
        with self.assertRaises(ValueError) as cm:
            ancom(self.table4, self.cats4, sig_test="ttest_ind")
        self.assertEqual(str(cm.exception), msg)


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
        self.data2 = np.vstack((
            np.random.multinomial(depth, self.p1, size=n),
            np.random.multinomial(depth, self.p2, size=n)))
        self.table2 = pd.DataFrame(self.data2)
        self.grouping2 = pd.Series(['Group1'] * n + ['Group2'] * n)

    def test_dirmult_ttest_demo(self):
        # The same example as in doctest
        data = np.array([
            [ 20, 110, 100, 101, 100, 103, 104],
            [ 33, 110, 120, 100, 101, 100, 102],
            [ 12, 110, 100, 110, 100,  50,  90],
            [202, 201,   9,  10,  10,  11,  11],
            [200, 202,  10,  10,  13,  10,  10],
            [203, 201,  14,  10,  10,  13,  12],
        ])
        samples = ["s1", "s2", "s3", "s4", "s5", "s6"]
        features = ["b1", "b2", "b3", "b4", "b5", "b6", "b7"]
        table = pd.DataFrame(data, index=samples, columns=features)
        labels = ["treatment", "treatment", "treatment",
                  "placebo", "placebo", "placebo"]
        grouping = pd.Series(labels, index=samples)

        obs = dirmult_ttest(table, grouping, 'treatment', 'placebo', seed=0)
        self.assertTupleEqual(obs.shape, (7, 7))
        self.assertListEqual(obs.index.to_list(), features)
        exp = {
            "T-statistic": [-17.179, -16.873,  6.943,  6.523,  6.654,  3.84,   7.601],
            "Log2(FC)":    [ -4.992,  -2.534,  1.628,  1.707,  1.528,  1.182,  1.48 ],
            "CI(2.5)":     [ -7.884,  -3.595, -1.048, -0.467, -1.037, -0.703, -0.601],
            "CI(97.5)":    [ -2.293,  -1.462,  4.751,  4.165,  3.978,  3.556,  4.044],
            "pvalue":      [  0.003,   0.001,  0.021,  0.013,  0.019,  0.045,  0.017],
            "qvalue":      [  0.02 ,   0.007,  0.068,  0.066,  0.068,  0.068,  0.068],
        }
        for key, value in exp.items():
            npt.assert_array_equal(obs[key].to_numpy().round(3), np.array(value))
        exp = np.array([True, True, False, False, False, False, False])
        npt.assert_array_equal(obs["Signif"].to_numpy(), exp)

    def test_dirmult_ttest_toy(self):
        p1 = np.array([5, 6, 7])
        p2 = np.array([4, 7, 7])
        p1, p2 = p1 / p1.sum(), p2 / p2.sum()
        depth = 1000
        n = 100
        data = np.vstack((
            np.random.multinomial(depth, p1, size=n),
            np.random.multinomial(depth, p2, size=n)))
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
        data = np.vstack((
            np.random.multinomial(depth, p1, size=n),
            np.random.multinomial(depth, p2, size=n)))
        table = pd.DataFrame(data)
        grouping = pd.Series(['Group1'] * n + ['Group2'] * n)
        exp_lfc = np.log2([4/5, 7/6, 7/7, 6/8, 5/9, 7/4])
        exp_lfc = (exp_lfc - exp_lfc.mean())  # convert to CLR coordinates
        res_100 = dirmult_ttest(table, grouping, self.treatment, self.reference)

        # increase sequencing depth by 100 fold
        depth = 10000
        data = np.vstack((
            np.random.multinomial(depth, p1, size=n),
            np.random.multinomial(depth, p2, size=n)))
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
        result = dirmult_ttest(self.table, self.grouping, self.treatment,
                               self.reference)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result.shape[1], 7)  # Expected number of columns
        pdt.assert_index_equal(result.index,
                               pd.Index(['feature1', 'feature2', 'feature3']))

    def test_dirmult_ttest_array_input(self):
        result = dirmult_ttest(self.table.to_numpy(), self.grouping, self.treatment,
                               self.reference)
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIsInstance(result.index, pd.RangeIndex)

    def test_dirmult_ttest_no_group(self):
        result = dirmult_ttest(self.table, self.grouping)
        self.assertIsInstance(result, pd.DataFrame)
        result = dirmult_ttest(self.table, self.grouping, treatment=self.treatment)
        self.assertIsInstance(result, pd.DataFrame)
        result = dirmult_ttest(self.table, self.grouping, reference=self.treatment)
        self.assertIsInstance(result, pd.DataFrame)

    def test_dirmult_ttest_no_pseudocount(self):
        result = dirmult_ttest(self.table, self.grouping, self.treatment,
                               self.reference, pseudocount=None)
        self.assertIsInstance(result, pd.DataFrame)

    def test_dirmult_ttest_no_p_adjust(self):
        result = dirmult_ttest(self.table, self.grouping, self.treatment,
                               self.reference, p_adjust=None)
        pdt.assert_series_equal(result['pvalue'], result['qvalue'], check_names=False)

    def test_dirmult_ttest_invalid_table_type(self):
        with self.assertRaises(TypeError):
            dirmult_ttest("invalid_table", self.grouping, self.treatment,
                          self.reference)

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

class DirMultLMETests(TestCase):
    # NOTE: `dirmult_lme` performs numerical optimization, which might (though rarely)
    # generate slightly different results on different platforms. The following tests
    # have specific numbers commented out, just to be safe in the CI workflow. But one
    # may check the accuracy of results locally by restoring the commented code.
    def setUp(self):
        np.random.seed(0)
        index = ["subject1", "subject2", "subject3", "subject4", "subject5", "subject6"]
        columns = ["feature1", "feature2", "feature3", "feature4"]

        # create sample data for testing
        self.table = pd.DataFrame(
            [[20, 110, 100, 101],
             [33, 110, 120, 100],
             [12, 110, 100, 110],
             [202, 201, 9, 10],
             [200, 202, 10, 10],
             [203, 201, 14, 10]],
            index=index,
            columns=columns)

        self.metadata = pd.DataFrame(
            {"Covar1": [1,1,2,2,3,3],
             "Covar2": [1,1,1,1,2,2],
             "Covar3": [1,2,1,2,1,2]},
            index=index)

    def test_dirmult_lme_demo(self):
        # a regular analysis
        res = dirmult_lme(
            table=self.table, metadata=self.metadata, formula="Covar2 + Covar3",
            grouping="Covar1", draws=1, seed=0, p_adjust="sidak")
        exp = """
  FeatureID Covariate  Reps  Log2(FC)   CI(2.5)  CI(97.5)    pvalue    qvalue  Signif
0  feature1    Covar2     1  2.708376 -0.818699  6.235451  0.132319  0.247129   False
1  feature1    Covar3     1  1.696770 -1.224053  4.617594  0.254876  0.444790   False
2  feature2    Covar2     1  0.956017 -0.736692  2.648726  0.268312  0.464632   False
3  feature2    Covar3     1  0.325451 -0.657935  1.308836  0.516565  0.766291   False
4  feature3    Covar2     1 -1.990268 -4.362246  0.381711  0.100061  0.190110   False
5  feature3    Covar3     1 -0.812892 -2.910615  1.284830  0.447548  0.694797   False
6  feature4    Covar2     1 -1.674125 -3.885791  0.537540  0.137915  0.256810   False
7  feature4    Covar3     1 -1.209329 -3.204871  0.786213  0.234925  0.414660   False
""".strip("\n")
        # self.assertEqual(str(res), exp)
        self.assertIsInstance(res, pd.DataFrame)
        self.assertTupleEqual(res.shape, (8, 9))
        pdt.assert_index_equal(res.columns, pd.Index([
            "FeatureID", "Covariate", "Reps", "Log2(FC)", "CI(2.5)", "CI(97.5)",
            "pvalue", "qvalue", "Signif"]))
        self.assertListEqual(res["FeatureID"].tolist(), [
            "feature1", "feature1", "feature2", "feature2", "feature3", "feature3",
            "feature4", "feature4"])
        self.assertListEqual(res["Covariate"].tolist(), [
            "Covar2", "Covar3", "Covar2", "Covar3", "Covar2", "Covar3", "Covar2",
            "Covar3"])
        self.assertTrue((res["Reps"] == 1).all())
        # npt.assert_array_equal(res["Log2(FC)"].round(5), np.array([
        #     2.70838, 1.69677, 0.95602, 0.32545, -1.99027, -0.81289, -1.67413,
        #     -1.20933]))
        # npt.assert_array_equal(res["CI(2.5)"].round(5), np.array([
        #     -0.8187, -1.22405, -0.73669, -0.65793, -4.36225, -2.91061, -3.88579,
        #     -3.20487]))
        # npt.assert_array_equal(res["CI(97.5)"].round(5), np.array([
        #     6.23545, 4.61759, 2.64873, 1.30884, 0.38171, 1.28483, 0.53754, 0.78621]))
        # npt.assert_array_equal(res["pvalue"].round(5), np.array([
        #     0.13232, 0.25488, 0.26831, 0.51657, 0.10006, 0.44755, 0.13792, 0.23492]))
        # npt.assert_array_equal(res["qvalue"].round(5), np.array([
        #     0.24713, 0.44479, 0.46463, 0.76629, 0.19011, 0.6948 , 0.25681, 0.41466]))
        self.assertTrue((res["Signif"] == False).all())

        # confirm that 2.5% < fold-change < 97.5%
        npt.assert_array_less(res["Log2(FC)"], res["CI(97.5)"])
        npt.assert_array_less(res["CI(2.5)"], res["Log2(FC)"])

    def test_dirmult_lme_array_input(self):
        res = dirmult_lme(
            table=self.table.to_numpy(), metadata=self.metadata,
            formula="Covar2 + Covar3", grouping="Covar1", draws=1, seed=0,
            p_adjust="sidak")
        self.assertIsInstance(res, pd.DataFrame)
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.247, 0.445, 0.465, 0.766, 0.190, 0.695, 0.257, 0.415]))

    def test_dirmult_lme_alt_grouping(self):
        res = dirmult_lme(
            table=self.table, metadata=self.metadata, formula="Covar2 + Covar3",
            grouping=self.metadata["Covar1"].to_numpy(), draws=1, seed=0,
            p_adjust="sidak")
        self.assertIsInstance(res, pd.DataFrame)
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.247, 0.445, 0.465, 0.766, 0.190, 0.695, 0.257, 0.415]))

    def test_dirmult_lme_no_pseudocount(self):
        res = dirmult_lme(
            table=self.table + 0.5, metadata=self.metadata, formula="Covar2 + Covar3",
            grouping="Covar1", draws=1, seed=0, p_adjust="sidak", pseudocount=None)
        self.assertIsInstance(res, pd.DataFrame)
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.247, 0.445, 0.465, 0.766, 0.190, 0.695, 0.257, 0.415]))

    def test_dirmult_lme_no_intercept(self):
        # "-1" at the end of formula suppresses intercept
        res = dirmult_lme(
            table=self.table, metadata=self.metadata, formula="Covar2 + Covar3 - 1",
            grouping="Covar1", draws=1, seed=0, p_adjust="sidak")
        self.assertIsInstance(res, pd.DataFrame)
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.988, 0.998, 0.129, 0.761, 0.766, 0.988, 0.971, 0.965]))

    def test_dirmult_lme_re_formula(self):
        # "1": random effect only in intercept (the default scenario)
        res = dirmult_lme(
            table=self.table, metadata=self.metadata, formula="Covar2 + Covar3",
            grouping="Covar1", draws=1, seed=0, p_adjust="sidak", re_formula="1")
        self.assertIsInstance(res, pd.DataFrame)
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.247, 0.445, 0.465, 0.766, 0.190, 0.695, 0.257, 0.415]))

        # add a random slope for Covar2
        res = dirmult_lme(
            table=self.table, metadata=self.metadata, formula="Covar2 + Covar3",
            grouping="Covar1", draws=1, seed=0, p_adjust="sidak",
            re_formula="1 + Covar2", fit_method="bfgs")
        self.assertIsInstance(res, pd.DataFrame)
        # NOTE: This test function is prone to generate different results on different
        # platforms. If the following test is to be executed, be prepared that lower
        # precision is needed to get the test pass.
        npt.assert_array_equal(res["qvalue"].round(3), np.array([
            0.274, 0.463, 0.645, 0.766, 0.251, 0.703, 0.400, 0.415]))

    def test_dirmult_lme_vc_formula(self):
        res = dirmult_lme(
            table=self.table, metadata=self.metadata, formula="Covar2 + Covar3",
            grouping="Covar1", draws=1, seed=0, p_adjust="sidak",
            vc_formula={"Covar2": "0 + C(Covar2)"})
        self.assertIsInstance(res, pd.DataFrame)
        # npt.assert_array_equal(res["qvalue"].round(3), np.array([
        #     0.247, 0.445, 0.465, 0.766, 0.190, 0.695, 0.257, 0.415]))

    def test_dirmult_lme_no_p_adjust_and_reml(self):
        result = dirmult_lme(
            self.table, self.metadata, "Covar2 + Covar3", "Covar1",
            draws=1, seed=0, p_adjust=None)
        pdt.assert_series_equal(result['pvalue'], result['qvalue'], check_names=False)

        fit_kwargs = {"reml": False}
        res_ml = dirmult_lme(
            self.table, self.metadata, "Covar2 + Covar3", "Covar1",
            draws=1, seed=0, p_adjust=None, fit_kwargs=fit_kwargs)

        with self.assertRaises(AssertionError):
            npt.assert_allclose(result['CI(97.5)'], res_ml['CI(97.5)'])
            npt.assert_allclose(result['CI(2.5)'], res_ml['CI(2.5)'])
            npt.assert_allclose(result['pvalue'], res_ml['pvalue'])

    def test_dirmult_lme_fit_warnings(self):
        # Convergence warnings are frequently raised during model fitting.
        from statsmodels.tools.sm_exceptions import ConvergenceWarning

        with self.assertWarns(ConvergenceWarning):
            dirmult_lme(table=self.table, metadata=self.metadata,
                        formula="Covar2 + Covar3", grouping="Covar1",
                        draws=1, seed=0, p_adjust="sidak",
                        fit_warnings=True)

    def test_dirmult_lme_fail_all(self):
        # Supply a non-existent optimization method to make it fail.
        msg = "LME fit failed for all features in all replicates."
        with self.assertRaises(ValueError) as cm:
            dirmult_lme(table=self.table, metadata=self.metadata,
                        formula="Covar2 + Covar3", grouping="Covar1",
                        draws=1, seed=0, p_adjust="sidak",
                        fit_method="not_a_method")
        self.assertEqual(str(cm.exception), msg)

    # def test_dirmult_lme_fail_some(self):
    #     # With the BFGS method, LME model fitting will not converge on two features.
    #     # Output will be NaN.
    #     msg = "LME fit failed for 2 features in all replicates, reporting NaNs."
    #     with self.assertWarns(UserWarning) as cm:
    #         res = dirmult_lme(table=self.table, metadata=self.metadata,
    #                           formula="Covar2 + Covar3", grouping="Covar1",
    #                           draws=1, seed=0, p_adjust="sidak",
    #                           fit_method="bfgs", fit_converge=True)
    #     self.assertEqual(str(cm.warning), msg)
    #     self.assertTrue(res.query(
    #         "FeatureID == ['feature1', 'feature3']")["Log2(FC)"].isnull().all())
    #     self.assertTrue(res.query(
    #         "FeatureID == ['feature2', 'feature4']")["Log2(FC)"].notnull().all())

    def test_dirmult_lme_invalid_table_type(self):
        with self.assertRaises(TypeError):
            dirmult_lme("not a table", self.metadata, "Covar2 + Covar3", "Covar1")

    def test_dirmult_lme_invalid_metadata_type(self):
        with self.assertRaises(TypeError):
            dirmult_lme(self.table, "not metadata", "Covar2 + Covar3", "Covar1")

    def test_dirmult_lme_invalid_grouping(self):
        with self.assertRaises(ValueError):
            dirmult_lme(self.table, self.metadata, "Covar2 + Covar3", "hello")

    def test_dirmult_lme_inconsistent_indexes(self):
        # change table index
        self.table.index = ["a", "b", "c", "d", "e", "f"]
        with self.assertRaises(ValueError):
            dirmult_lme(self.table, self.metadata, "Covar2 + Covar3", "Covar1")

    def test_dirmult_lme_single_feature(self):
        # only keep the first column of the table
        self.table = self.table.iloc[:, :1]
        with self.assertRaises(ValueError):
            dirmult_lme(self.table, self.metadata, "Covar1", "Covar1")

    def test_dirmult_lme_toy_data(self):
        # simulate a dataset of 20 samples by 3 features, with 10 repeated measurements
        # (covar1) and 2 sample groups (covar2)
        p1 = np.array([5, 6, 7])
        p2 = np.array([4, 7, 7])
        p1, p2 = p1 / p1.sum(), p2 / p2.sum()
        depth = 1000
        n = 10
        index_range = range(1, n * 2 + 1)
        data = np.vstack((
            np.random.multinomial(depth, p1, size=n),
            np.random.multinomial(depth, p2, size=n)))
        table = pd.DataFrame(data)
        table.columns = ["feature1", "feature2", "feature3"]
        table.index = [f"subject{i}" for i in index_range]

        exp_lfc = np.log2([4/5, 7/6, 7/7])
        exp_lfc = (exp_lfc - exp_lfc.mean())  # convert to CLR coordinates

        metadata = pd.DataFrame({
            "covar1": [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10],
            "covar2": [1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2]},
            index=[f"subject{i}" for i in index_range])

        res = dirmult_lme(
            table=table, metadata=metadata, formula="covar2", grouping="covar1",
            draws=8, seed=0, p_adjust="sidak")

        npt.assert_array_equal(res["Log2(FC)"].round(5), [-0.28051, 0.32118, -0.04067])
        npt.assert_array_equal(res["CI(2.5)"].round(5), [-0.41639, 0.18745, -0.16542])
        npt.assert_array_equal(res["CI(97.5)"].round(5), [-0.14359, 0.46473, 0.07706])

        # confirm expected fold change is within confidence interval
        npt.assert_array_less(exp_lfc, res['CI(97.5)'])
        npt.assert_array_less(res['CI(2.5)'], exp_lfc)

    def test_dirmult_lme_toy_data_depth(self):
        p1 = np.array([5, 6, 7, 8, 9, 4])
        p2 = np.array([4, 7, 7, 6, 5, 7])
        p1, p2 = p1 / p1.sum(), p2 / p2.sum()
        depth = 100
        n = 10
        index_range = range(1, n * 2 + 1)
        data = np.vstack((
            np.random.multinomial(depth, p1, size=n),
            np.random.multinomial(depth, p2, size=n)))
        table = pd.DataFrame(data)
        table.columns = [
            "feature1", "feature2", "feature3", "feature4", "feature5", "feature6"]
        table.index = [f"subject{i}" for i in index_range]

        metadata = pd.DataFrame({
            "covar1": [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10],
            "covar2": [1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2]},
            index=[f"subject{i}" for i in index_range])

        exp_lfc = np.log2([4/5, 7/6, 7/7, 6/8, 5/9, 7/4])
        exp_lfc = (exp_lfc - exp_lfc.mean())  # convert to CLR coordinates

        res_100 = dirmult_lme(
            table=table, metadata=metadata, formula="covar2", grouping="covar1",
            draws=8, seed=0, p_adjust="sidak")

        # increase sequencing depth by 100 fold
        depth = 10000
        data = np.vstack((
            np.random.multinomial(depth, p1, size=n),
            np.random.multinomial(depth, p2, size=n)))
        table = pd.DataFrame(data)
        table.columns = [
            "feature1", "feature2", "feature3", "feature4", "feature5", "feature6"]
        table.index = [f"subject{i}" for i in index_range]
        metadata.index = [f"subject{i}" for i in index_range]
        res_10000 = dirmult_lme(
            table=table, metadata=metadata, formula="covar2", grouping="covar1",
            draws=8, seed=0, p_adjust="sidak")

        # when the sequencing depth increases, the confidence intervals
        # should also shrink
        npt.assert_array_less(res_100['CI(2.5)'], res_10000['CI(2.5)'])
        npt.assert_array_less(res_10000['CI(97.5)'], res_100['CI(97.5)'])


if __name__ == "__main__":
    main()
