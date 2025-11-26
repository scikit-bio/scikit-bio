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
    _type_cast_to_float,
    _check_sig_test,
    _check_p_adjust,
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


if __name__ == "__main__":
    main()
