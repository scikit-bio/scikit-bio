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
from numpy.random import normal
import pandas as pd
import pandas.testing as pdt
from scipy.stats import f_oneway, ConstantInputWarning

from skbio.util import assert_data_frame_almost_equal
from skbio.stats.composition import closure
from skbio.stats.composition._ancom import ancom


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


if __name__ == "__main__":
    main()
