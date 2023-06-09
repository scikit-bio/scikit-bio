# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import os
import itertools
import unittest

import pandas as pd
import numpy as np
import numpy.testing as npt

from skbio import OrdinationResults
from skbio.util import (get_data_path, assert_ordination_results_equal,
                        assert_data_frame_almost_equal)
from skbio.util._testing import _normalize_signs, assert_series_almost_equal


class TestGetDataPath(unittest.TestCase):
    def test_get_data_path(self):
        fn = 'parrot'
        path = os.path.dirname(os.path.abspath(__file__))
        data_path = os.path.join(path, 'data', fn)
        data_path_2 = get_data_path(fn)
        self.assertEqual(data_path_2, data_path)


class TestAssertOrdinationResultsEqual(unittest.TestCase):
    def test_assert_ordination_results_equal(self):
        minimal1 = OrdinationResults('foo', 'bar', pd.Series([1.0, 2.0]),
                                     pd.DataFrame([[1, 2, 3], [4, 5, 6]]))

        # a minimal set of results should be equal to itself
        assert_ordination_results_equal(minimal1, minimal1)

        # type mismatch
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal(minimal1, 'foo')

        # numeric values should be checked that they're almost equal
        almost_minimal1 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.0000001, 1.9999999]),
            pd.DataFrame([[1, 2, 3], [4, 5, 6]]))
        assert_ordination_results_equal(minimal1, almost_minimal1)

        # test each of the optional numeric attributes
        for attr in ('features', 'samples', 'biplot_scores',
                     'sample_constraints'):
            # missing optional numeric attribute in one, present in the other
            setattr(almost_minimal1, attr, pd.DataFrame([[1, 2], [3, 4]]))
            with npt.assert_raises(AssertionError):
                assert_ordination_results_equal(minimal1, almost_minimal1)
            setattr(almost_minimal1, attr, None)

            # optional numeric attributes present in both, but not almost equal
            setattr(minimal1, attr, pd.DataFrame([[1, 2], [3, 4]]))
            setattr(almost_minimal1, attr, pd.DataFrame([[1, 2],
                                                         [3.00002, 4]]))
            with npt.assert_raises(AssertionError):
                assert_ordination_results_equal(minimal1, almost_minimal1)
            setattr(minimal1, attr, None)
            setattr(almost_minimal1, attr, None)

            # optional numeric attributes present in both, and almost equal
            setattr(minimal1, attr, pd.DataFrame([[1.0, 2.0], [3.0, 4.0]]))
            setattr(almost_minimal1, attr,
                    pd.DataFrame([[1.0, 2.0], [3.00000002, 4]]))
            assert_ordination_results_equal(minimal1, almost_minimal1)
            setattr(minimal1, attr, None)
            setattr(almost_minimal1, attr, None)

        # missing optional numeric attribute in one, present in the other
        almost_minimal1.proportion_explained = pd.Series([1, 2, 3])
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal(minimal1, almost_minimal1)
        almost_minimal1.proportion_explained = None

        # optional numeric attributes present in both, but not almost equal
        minimal1.proportion_explained = pd.Series([1, 2, 3])
        almost_minimal1.proportion_explained = pd.Series([1, 2, 3.00002])
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal(minimal1, almost_minimal1)
        almost_minimal1.proportion_explained = None
        almost_minimal1.proportion_explained = None

        # optional numeric attributes present in both, and almost equal
        minimal1.proportion_explained = pd.Series([1, 2, 3])
        almost_minimal1.proportion_explained = pd.Series([1, 2, 3.00000002])
        assert_ordination_results_equal(minimal1, almost_minimal1)
        almost_minimal1.proportion_explained = None
        almost_minimal1.proportion_explained = None


class TestNormalizeSigns(unittest.TestCase):
    def test_shapes_and_nonarray_input(self):
        with self.assertRaises(ValueError):
            _normalize_signs([[1, 2], [3, 5]], [[1, 2]])

    def test_works_when_different(self):
        """Taking abs value of everything would lead to false
        positives."""
        a = np.array([[1, -1],
                      [2, 2]])
        b = np.array([[-1, -1],
                      [2, 2]])
        with self.assertRaises(AssertionError):
            npt.assert_equal(*_normalize_signs(a, b))

    def test_easy_different(self):
        a = np.array([[1, 2],
                      [3, -1]])
        b = np.array([[-1, 2],
                      [-3, -1]])
        npt.assert_equal(*_normalize_signs(a, b))

    def test_easy_already_equal(self):
        a = np.array([[1, -2],
                      [3, 1]])
        b = a.copy()
        npt.assert_equal(*_normalize_signs(a, b))

    def test_zeros(self):
        a = np.array([[0, 3],
                      [0, -1]])
        b = np.array([[0, -3],
                      [0, 1]])
        npt.assert_equal(*_normalize_signs(a, b))

    def test_hard(self):
        a = np.array([[0, 1],
                      [1, 2]])
        b = np.array([[0, 1],
                      [-1, 2]])
        npt.assert_equal(*_normalize_signs(a, b))

    def test_harder(self):
        """We don't want a value that might be negative due to
        floating point inaccuracies to make a call to allclose in the
        result to be off."""
        a = np.array([[-1e-15, 1],
                      [5, 2]])
        b = np.array([[1e-15, 1],
                      [5, 2]])
        # Clearly a and b would refer to the same "column
        # eigenvectors" but a slopppy implementation of
        # _normalize_signs could change the sign of column 0 and make a
        # comparison fail
        npt.assert_almost_equal(*_normalize_signs(a, b))

    def test_column_zeros(self):
        a = np.array([[0, 1],
                      [0, 2]])
        b = np.array([[0, -1],
                      [0, -2]])
        npt.assert_equal(*_normalize_signs(a, b))

    def test_column_almost_zero(self):
        a = np.array([[1e-15, 3],
                      [-2e-14, -6]])
        b = np.array([[0, 3],
                      [-1e-15, -6]])
        npt.assert_almost_equal(*_normalize_signs(a, b))


class TestAssertDataFrameAlmostEqual(unittest.TestCase):
    def setUp(self):
        self.df = pd.DataFrame({'bar': ['a', 'b', 'cd', 'e'],
                                'foo': [42, 42.0, np.nan, 0]})

    def test_not_equal(self):
        unequal_dfs = [
            self.df,
            # floating point error too large to be "almost equal"
            pd.DataFrame({'bar': ['a', 'b', 'cd', 'e'],
                          'foo': [42, 42.001, np.nan, 0]}),
            # extra NaN
            pd.DataFrame({'bar': ['a', 'b', 'cd', 'e'],
                          'foo': [42, np.nan, np.nan, 0]}),
            # different column order
            pd.DataFrame(self.df, columns=['foo', 'bar']),
            # different index order
            pd.DataFrame(self.df, index=np.arange(4)[::-1]),
            # different index type
            pd.DataFrame(self.df, index=np.arange(4).astype(float)),
            # various forms of "empty" DataFrames that are not equivalent
            pd.DataFrame(),
            pd.DataFrame(index=np.arange(10)),
            pd.DataFrame(columns=np.arange(10)),
            pd.DataFrame(index=np.arange(10), columns=np.arange(10)),
            pd.DataFrame(index=np.arange(9)),
            pd.DataFrame(columns=np.arange(9)),
            pd.DataFrame(index=np.arange(9), columns=np.arange(9))
        ]

        # each df should compare equal to itself and a copy of itself
        for df in unequal_dfs:
            assert_data_frame_almost_equal(df, df)
            assert_data_frame_almost_equal(df, pd.DataFrame(df, copy=True))

        # every pair of dfs should not compare equal. use permutations instead
        # of combinations to test that comparing df1 to df2 and df2 to df1 are
        # both not equal
        for df1, df2 in itertools.permutations(unequal_dfs, 2):
            with self.assertRaises(AssertionError):
                assert_data_frame_almost_equal(df1, df2)

    def test_equal(self):
        equal_dfs = [
            self.df,
            # floating point error small enough to be "almost equal"
            pd.DataFrame({'bar': ['a', 'b', 'cd', 'e'],
                          'foo': [42, 42.00001, np.nan, 0]})
        ]

        for df in equal_dfs:
            assert_data_frame_almost_equal(df, df)

        for df1, df2 in itertools.permutations(equal_dfs, 2):
            assert_data_frame_almost_equal(df1, df2)


class TestAssertSeriesAlmostEqual(unittest.TestCase):

    def setUp(self):
        self.series = [
            pd.Series(dtype='float64'),
            pd.Series(dtype=object),
            pd.Series(dtype='int64'),
            pd.Series([1, 2, 3]),
            pd.Series([3, 2, 1]),
            pd.Series([1, 2, 3, 4]),
            pd.Series([1., 2., 3.]),
            pd.Series([1, 2, 3], [1.0, 2.0, 3.0]),
            pd.Series([1, 2, 3], [1, 2, 3]),
            pd.Series([1, 2, 3], ['c', 'b', 'a']),
            pd.Series([3, 2, 1], ['c', 'b', 'a']),
        ]

    def test_not_equal(self):
        # no pair of series should compare equal
        for s1, s2 in itertools.permutations(self.series, 2):
            with self.assertRaises(AssertionError):
                assert_series_almost_equal(s1, s2)

    def test_equal(self):
        s1 = pd.Series([1., 2., 3.])
        s2 = pd.Series([1.000001, 2., 3.])
        assert_series_almost_equal(s1, s2)

        # all series should be equal to themselves and copies of themselves
        for s in self.series:
            assert_series_almost_equal(s, s)
            assert_series_almost_equal(s, pd.Series(s, copy=True))


if __name__ == '__main__':
    unittest.main()
