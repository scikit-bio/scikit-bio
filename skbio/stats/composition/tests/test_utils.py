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
import pandas as pd
import pandas.testing as pdt

from skbio.stats.composition._utils import (
    _check_grouping,
    _check_trt_ref_groups,
    _check_metadata,
    _build_dmatrix,
    _check_sig_test,
    _adjust_pvalues,
    _sm_p_adjust,
)


class UtilsTests(TestCase):

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

    def test_build_dmatrix(self):
        meta = pd.DataFrame({
            "num": [1, 2, 3],  # as continuous
            "str": ["orange", "apple", "pear"],  # sorted alphabetically
            "num_str": ["1", "2", "3"],  # categorical
            "bool": [True, False, True],  # categorical
            "cat": pd.Categorical(["red", "blue", "blue"], categories=[
                "red", "blue"], ordered=True)})  # no sorting
        obs = _build_dmatrix(
            "num + str + num_str + bool + cat", meta)
        exp_covars = [
            "Intercept",
            "str[T.orange]",
            "str[T.pear]",
            "num_str[T.2]",
            "num_str[T.3]",
            "bool[T.True]",
            "cat[T.blue]",
            "num"]  # continuous variable moves after categorical
        exp_mat = np.array([
            [1, 1, 0, 0, 0, 1, 0, 1],
            [1, 0, 0, 1, 0, 0, 1, 2],
            [1, 0, 1, 0, 1, 1, 1, 3]], dtype=float)
        self.assertListEqual(obs.design_info.column_names, exp_covars)
        npt.assert_array_equal(obs, exp_mat)

        # specify data type
        obs = _build_dmatrix("num", meta, dtype=np.float32)
        self.assertEqual(obs.dtype, np.float32)
        npt.assert_array_equal(obs, [[1, 1], [1, 2], [1, 3]])

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


class AdjustPvaluesTests(TestCase):

    def test_methods(self):
        pval = np.array([0.04, 0.01, 0.03, 0.2])
        original = pval.copy()
        for method in ("holm", "holm-bonferroni", "HOLM"):
            obs = _adjust_pvalues(pval, method)
            npt.assert_allclose(obs, [0.09, 0.04, 0.09, 0.2])
        for method in ("bh", "benjamini-hochberg", "BH"):
            obs = _adjust_pvalues(pval, method)
            npt.assert_allclose(obs, [0.16 / 3, 0.04, 0.16 / 3, 0.2])
        npt.assert_array_equal(_adjust_pvalues(pval), obs)
        npt.assert_array_equal(pval, original)
        self.assertFalse(np.shares_memory(obs, pval))

    def test_bonferroni_and_by(self):
        pval = np.array([0.04, 0.01, 0.03, 0.2])
        for method in ("bonferroni", "bonf", "BONFERRONI"):
            npt.assert_allclose(_adjust_pvalues(pval, method),
                                [0.16, 0.04, 0.12, 0.8])
        # H_4 = 25/12; scale the BH values by this factor.
        for method in ("by", "benjamini-yekutieli", "BY"):
            npt.assert_allclose(_adjust_pvalues(pval, method),
                                [1 / 9, 1 / 12, 1 / 9, 5 / 12])

    def test_ties_and_bounds(self):
        for method in ("holm", "bh", "bonf", "by"):
            obs = _adjust_pvalues([1, 0.5, 0, 0.5], method)
            exp = ([1, 2 / 3, 0, 2 / 3] if method == "bh" else [1, 1, 0, 1])
            npt.assert_allclose(obs, exp)
            npt.assert_array_equal(_adjust_pvalues(np.zeros(5), method), 0)
            npt.assert_array_equal(_adjust_pvalues(np.ones(5), method), 1)

    def test_nan(self):
        # Expected results match R p.adjust with its default comparison count.
        pval = np.array([[0.01, np.nan, 0.4],
                         [np.nan, np.nan, 0.2],
                         [0.04, np.nan, np.nan],
                         [np.nan, np.nan, 0.1]])
        original = pval.copy()
        for method in ("holm", "bh"):
            exp = np.array([[0.02, np.nan, 0.4],
                            [np.nan, np.nan, 0.4],
                            [0.04, np.nan, np.nan],
                            [np.nan, np.nan, 0.3]])
            if method == "bh":
                exp[1, 2] = 0.3
            npt.assert_allclose(_adjust_pvalues(pval, method), exp)
            npt.assert_array_equal(pval, original)
            obs = _adjust_pvalues([np.nan, 0.04, np.nan], method)
            npt.assert_array_equal(obs, [np.nan, 0.04, np.nan])

    def test_axes(self):
        from statsmodels.stats.multitest import multipletests

        rng = np.random.default_rng(42)
        pval = rng.random((3, 4, 5)) ** 4
        pval[0, 1, 2] = np.nan
        original = pval.copy()
        for method, sm_method in (("holm", "holm"), ("bh", "fdr_bh"),
                                  ("bonferroni", "bonferroni"), ("by", "fdr_by")):
            for axis in range(3):
                exp = np.empty_like(pval)
                data = np.moveaxis(pval, axis, -1)
                result = np.moveaxis(exp, axis, -1)
                for idx in np.ndindex(data.shape[:-1]):
                    col = data[idx]
                    valid = ~np.isnan(col)
                    result[idx] = np.nan
                    result[idx][valid] = multipletests(
                        col[valid], method=sm_method)[1]
                for ax in (axis, axis - 3, np.int64(axis)):
                    obs = _adjust_pvalues(pval, method, axis=ax)
                    npt.assert_allclose(obs, exp, rtol=1e-14, atol=0)
                view = pval.transpose(2, 1, 0)[::-1]
                out = np.empty_like(pval).transpose(2, 1, 0)[::-1]
                obs = _adjust_pvalues(view, method, axis=2 - axis, out=out)
                self.assertIs(obs, out)
                npt.assert_allclose(obs, exp.transpose(2, 1, 0)[::-1],
                                    rtol=1e-14, atol=0)
        npt.assert_array_equal(pval, original)

    def test_axis_none(self):
        pval = np.array([[0.04, 0.01], [0.03, 0.2]])
        exp = [[0.09, 0.04], [0.09, 0.2]]
        npt.assert_allclose(_adjust_pvalues(pval, "holm", axis=None), exp)
        for value in (0.04, np.nan):
            obs = _adjust_pvalues(value, axis=None)
            self.assertEqual(obs.shape, ())
            npt.assert_array_equal(obs, value)
        out = np.array(0.04)
        self.assertIs(_adjust_pvalues(out, axis=None, out=out), out)
        npt.assert_array_equal(out, 0.04)
        pval[0, 0] = np.nan
        exp = [[np.nan, 0.03], [0.06, 0.2]]
        npt.assert_allclose(_adjust_pvalues(pval, "holm", axis=None), exp)

    def test_empty_and_singleton(self):
        for method in ("holm", "bh", "bonferroni", "by"):
            for shape in ((0,), (0, 3), (3, 0), (2, 0, 4)):
                pval = np.empty(shape)
                for axis in (*range(len(shape)), None):
                    obs = _adjust_pvalues(pval, method, axis=axis)
                    npt.assert_array_equal(obs, pval)
            pval = np.array([[0.01, np.nan, 0.8]])
            npt.assert_array_equal(_adjust_pvalues(pval, method), pval)
            out = np.empty_like(pval)
            self.assertIs(_adjust_pvalues(pval, method, out=out), out)
            npt.assert_array_equal(out, pval)

    def test_dtypes(self):
        for dtype in (np.float16, np.float32, np.float64, np.longdouble):
            pval = np.array([0.04, 0.01, 0.03, 0.2], dtype=dtype)
            for method in ("holm", "bh", "bonferroni", "by"):
                obs = _adjust_pvalues(pval, method)
                self.assertEqual(obs.dtype, dtype)
                exp = _adjust_pvalues(pval.astype(np.float64), method)
                # The reference is limited to float64 even for longdouble.
                tol = 2 * max(np.finfo(dtype).eps, np.finfo(np.float64).eps)
                npt.assert_allclose(obs, exp, rtol=tol)

        # Non-floating point data types are prohibited.
        msg = "`pval` must have a floating-point data type."
        for dtype in (np.int8, np.uint64):
            with self.assertRaises(TypeError) as cm:
                _adjust_pvalues(np.array([0, 1, 1], dtype=dtype))
            self.assertEqual(str(cm.exception), msg)

        # Low-precision output must not overflow during intermediate scaling.
        for method in ("holm", "bh", "bonferroni", "by"):
            with np.errstate(over="raise", invalid="raise"):
                obs = _adjust_pvalues(np.ones(70000, dtype=np.float16), method)
            npt.assert_array_equal(obs, 1)

    def test_layouts_and_out(self):
        pval = np.array([[0.04, np.nan, 0.1],
                         [0.01, 0.04, 0.4],
                         [0.03, 0.01, 0.2],
                         [0.2, 0.2, np.nan]])
        for method in ("holm", "bh", "bonferroni", "by"):
            for p in (pval.copy(), np.asfortranarray(pval), pval.T,
                      pval[::-1, ::-1], pval[::2, ::2]):
                original = p.copy()
                for axis in (0, 1, None):
                    exp = _adjust_pvalues(p.copy(), method, axis=axis)
                    npt.assert_allclose(_adjust_pvalues(p, method, axis=axis), exp)
                    # Use a strided buffer so reshaping it would silently copy.
                    backing = np.full((p.shape[0] * 2, p.shape[1] * 2), -1.)
                    out = backing[::2, ::-2]
                    obs = _adjust_pvalues(p, method, axis=axis, out=out)
                    self.assertIs(obs, out)
                    npt.assert_allclose(out, exp)
                    npt.assert_array_equal(backing[1::2], -1)
                    npt.assert_array_equal(p, original)
                    work = p.copy(order="K")
                    self.assertIs(_adjust_pvalues(
                        work, method, axis=axis, out=work), work)
                    npt.assert_allclose(work, exp)
        # Exact aliasing also works on a reversed, noncontiguous input.
        work = pval.copy()[::-1, ::-1]
        exp = _adjust_pvalues(work)
        _adjust_pvalues(work, out=work)
        npt.assert_allclose(work, exp)
        pval.flags.writeable = False
        npt.assert_allclose(_adjust_pvalues(pval), _adjust_pvalues(pval.copy()))
        pval = np.array([0.01, 0.04], dtype=np.float32)
        out = np.empty_like(pval)
        self.assertIs(_adjust_pvalues(pval, out=out), out)
        npt.assert_allclose(out, [0.02, 0.04])

    def test_varying_family_sizes(self):
        from statsmodels.stats.multitest import multipletests

        pval = np.tile([0.04, 0.01, 0.03, 0.2, 0.6], (6, 1)).T
        pval[:2, 1] = np.nan
        pval[:4, 2] = np.nan
        pval[:, 3] = np.nan
        pval[::2, 4] = np.nan
        for method, sm_method in (("holm", "holm"), ("bh", "fdr_bh"),
                                  ("bonferroni", "bonferroni"), ("by", "fdr_by")):
            exp = np.full_like(pval, np.nan)
            for col in range(pval.shape[1]):
                valid = ~np.isnan(pval[:, col])
                if valid.any():
                    values = pval[valid, col]
                    exp[valid, col] = multipletests(
                        values, method=sm_method)[1]
            npt.assert_allclose(_adjust_pvalues(pval, method), exp)
            work = pval.copy()
            _adjust_pvalues(work, method, out=work)
            npt.assert_allclose(work, exp)

    def test_n_tests(self):
        from statsmodels.stats.multitest import multipletests

        pval = np.array([[0.01, 0.4, np.nan], [0.04, np.nan, np.nan],
                         [0.2, 0.01, np.nan]])
        original = pval.copy()
        for method, sm_method in (("bonf", "bonferroni"), ("holm", "holm"),
                                  ("bh", "fdr_bh"), ("by", "fdr_by"),
                                  ("sidak", "sidak"), ("hommel", "hommel")):
            for axis in (0, 1, None):
                for n_tests in (5, 10, 10.0):
                    exp = np.full_like(pval, np.nan)
                    if axis is None:
                        families = [(pval.ravel(), exp.reshape(-1))]
                    else:
                        families = zip(np.moveaxis(pval, axis, -1),
                                       np.moveaxis(exp, axis, -1))
                    for col, dest in families:
                        valid = ~np.isnan(col)
                        n = valid.sum()
                        if n:
                            padded = np.pad(col[valid], (0, int(n_tests) - n),
                                            constant_values=1.)
                            dest[valid] = multipletests(padded, method=sm_method)[1][:n]
                    obs = _adjust_pvalues(pval, method, axis=axis, n_tests=n_tests)
                    npt.assert_allclose(obs, exp, rtol=1e-14, atol=0)
                    work = pval.copy()
                    self.assertIs(_adjust_pvalues(
                        work, method, axis=axis, n_tests=n_tests, out=work), work)
                    npt.assert_allclose(work, exp, rtol=1e-14, atol=0)
        npt.assert_array_equal(pval, original)

    def test_n_tests_fractional(self):
        # R: p.adjust(c(.01, .02, .03), method, n = 7.5).
        # BY uses sum(1 / (1:7.5)) = H_7, not a continuous harmonic number.
        expected = {
            "bonf": [0.075, 0.15, 0.225],
            "holm": [0.075, 0.13, 0.165],
            "bh": [0.075] * 3,
            "by": [0.1944642857142857] * 3,
        }
        aliases = {"bonf": "bonferroni", "holm": "holm-bonferroni",
                   "bh": "benjamini-hochberg", "by": "benjamini-yekutieli"}
        for method, exp in expected.items():
            for name in (method, aliases[method]):
                pval = np.array([[0.03, np.nan, 0.01, 0.02]])
                target = np.array([[exp[2], np.nan, exp[0], exp[1]]])
                for axis in (1, None):
                    work = pval.copy()
                    obs = _adjust_pvalues(work, name, axis=axis, n_tests=7.5,
                                         out=work)
                    self.assertIs(obs, work)
                    npt.assert_allclose(obs, target, rtol=1e-14)
        for method in ("sidak", "hommel", "b", "h", "fdr_bh", "fdr_by"):
            with self.assertRaisesRegex(ValueError, "Fractional `n_tests`"):
                _adjust_pvalues(np.array([0.01]), method, n_tests=7.5)

    def test_n_tests_counts_and_clipping(self):
        pval = np.array([np.nan, 0.01, np.nan, 0.04, np.nan])
        for method in ("bonf", "holm", "bh", "by", "sidak"):
            # The specified count need only cover nonmissing entries.
            obs = _adjust_pvalues(pval, method, n_tests=2)
            npt.assert_allclose(obs, _adjust_pvalues(pval, method), rtol=1e-14)
            obs = _adjust_pvalues(np.array([0.6, 0.8]), method, n_tests=5)
            if method != "sidak":
                npt.assert_array_equal(obs, [1., 1.])
            obs = _adjust_pvalues(np.array([np.nan, np.nan]), method, n_tests=0)
            npt.assert_array_equal(obs, [np.nan, np.nan])
            npt.assert_array_equal(_adjust_pvalues(np.array([]), method, n_tests=5), [])
        npt.assert_array_equal(_adjust_pvalues(pval, None, n_tests=10), pval)
        npt.assert_allclose(_adjust_pvalues(np.array([0.04]), "by", n_tests=1), [0.04])

    def test_n_tests_large(self):
        # Native methods must not allocate arrays proportional to n_tests.
        n_tests = 10**12
        pval = np.array([0.01, 0.02, 0.03]) / n_tests
        harmonic = np.log(n_tests) + np.euler_gamma + 0.5 / n_tests
        expected = {
            "bonf": [0.01, 0.02, 0.03],
            "holm": pval * (n_tests - np.arange(3)),
            "bh": np.full(3, 0.01),
            "by": np.full(3, 0.01 * harmonic),
        }
        for method, exp in expected.items():
            npt.assert_allclose(_adjust_pvalues(pval, method, n_tests=n_tests),
                                exp, rtol=1e-14, atol=0)

    def test_fallback(self):
        from statsmodels.stats.multitest import multipletests

        pval = np.array([[0.01, np.nan, 0.03],
                         [np.nan, np.nan, 0.2],
                         [0.04, np.nan, 0.5]])
        original = pval.copy()
        for method in ("sidak", "holm-sidak"):
            exp = np.full_like(pval, np.nan)
            for col in (0, 2):
                valid = ~np.isnan(pval[:, col])
                exp[valid, col] = multipletests(
                    pval[valid, col], method=method)[1]
            npt.assert_allclose(_adjust_pvalues(pval, method), exp)
            work = pval.copy()
            self.assertIs(_adjust_pvalues(work, method, out=work), work)
            npt.assert_allclose(work, exp)
            npt.assert_array_equal(pval, original)
        npt.assert_array_equal(_adjust_pvalues(np.empty((0, 2)), "sidak"),
                               np.empty((0, 2)))

    def test_no_adjustment(self):
        pval = np.array([[0.01, np.nan], [0.04, 0.2]])
        obs = _adjust_pvalues(pval, None)
        npt.assert_array_equal(obs, pval)
        self.assertFalse(np.shares_memory(obs, pval))
        out = np.empty_like(pval)
        self.assertIs(_adjust_pvalues(pval, None, out=out), out)
        npt.assert_array_equal(out, pval)
        self.assertIs(_adjust_pvalues(pval, None, out=pval), pval)

    def test_sm_p_adjust(self):
        self.assertIsNone(_sm_p_adjust(None))

        pval = [0.005, 0.011, 0.02, 0.04, 0.13]
        obs = _sm_p_adjust("holm")(pval)
        exp = pval * np.arange(1, 6)[::-1]
        for a, b in zip(obs, exp):
            self.assertAlmostEqual(a, b)

        pval = [0.005, 0.011, 0.02, 0.04, 0.13]
        obs = _sm_p_adjust("fdr_bh")(pval)
        exp = [0.025, 0.0275, 0.03333333, 0.05, 0.13]
        for a, b in zip(obs, exp):
            self.assertAlmostEqual(a, b)

        msg = "'hello' is not an available multiple testing correction method."
        with self.assertRaisesRegex(ValueError, msg):
            _sm_p_adjust("hello")(pval)


if __name__ == "__main__":
    main()
