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
from unittest.mock import patch

import pandas as pd
import numpy as np
import numpy.testing as npt

from skbio import OrdinationResults
from skbio.util import (get_data_path, assert_ordination_results_equal,
                        assert_data_frame_almost_equal)
from skbio.util._testing import (
    _normalize_signs,
    assert_series_almost_equal,
    assert_index_equal,
    assert_ordination_results_equal_np,
    _assert_series_equal,
    _assert_frame_dists_equal,
    _data_frame_to_default_int_type,
    ReallyEqualMixin,
    _read_env,
    _get_available_backends,
    _should_run,
    xp_assert_close,
    xp_assert_equal,
    backends,
    ArrayAPITestMixin,
)


class TestGetDataPath(unittest.TestCase):
    def test_get_data_path(self):
        fn = 'parrot'
        path = os.path.dirname(os.path.abspath(__file__))
        data_path = os.path.join(path, 'data', fn)
        data_path_2 = get_data_path(fn)
        self.assertEqual(data_path_2, data_path)

    def test_get_data_path_subfolder(self):
        fn = 'parrot'
        path = os.path.dirname(os.path.abspath(__file__))
        data_path = os.path.join(path, 'custom_subfolder', fn)
        data_path_2 = get_data_path(fn, subfolder='custom_subfolder')
        self.assertEqual(data_path_2, data_path)


class TestReallyEqualMixin(unittest.TestCase, ReallyEqualMixin):
    def test_assert_really_equal(self):
        self.assertReallyEqual(1, 1)
        self.assertReallyEqual('a', 'a')
        self.assertReallyEqual([1, 2], [1, 2])

    def test_assert_really_not_equal(self):
        self.assertReallyNotEqual(1, 2)
        self.assertReallyNotEqual('a', 'b')
        self.assertReallyNotEqual(1, '1')

    def test_assert_really_equal_fails_on_unequal(self):
        with self.assertRaises(AssertionError):
            self.assertReallyEqual(1, 2)

    def test_assert_really_not_equal_fails_on_equal(self):
        with self.assertRaises(AssertionError):
            self.assertReallyNotEqual(1, 1)


class TestAssertOrdinationResultsEqual(unittest.TestCase):
    def setUp(self):
        self.minimal1 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.0, 2.0]),
            pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        self.minimal2 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.0000001, 1.9999999]),
            pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))

    def test_equal_to_itself(self):
        assert_ordination_results_equal(self.minimal1, self.minimal1)

    def test_type_mismatch(self):
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal(self.minimal1, 'foo')

    def test_almost_equal_numeric(self):
        assert_ordination_results_equal(self.minimal1, self.minimal2)

    def test_optional_numeric_attributes(self):
        for attr in ('features', 'samples', 'biplot_scores',
                     'sample_constraints'):
            m1 = OrdinationResults(
                'foo', 'bar',
                pd.Series([1.0, 2.0]),
                pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
            m2 = OrdinationResults(
                'foo', 'bar',
                pd.Series([1.0, 2.0]),
                pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))

            # missing in one, present in other
            setattr(m2, attr, pd.DataFrame([[1.0, 2.0], [3.0, 4.0]]))
            with npt.assert_raises(AssertionError):
                assert_ordination_results_equal(m1, m2)
            setattr(m2, attr, None)

            # both present, not almost equal (different pairwise distances)
            setattr(m1, attr, pd.DataFrame([[1.0, 2.0], [3.0, 4.0]]))
            setattr(m2, attr, pd.DataFrame([[1.0, 2.0], [3.5, 4.0]]))
            with npt.assert_raises(AssertionError):
                assert_ordination_results_equal(m1, m2)
            setattr(m1, attr, None)
            setattr(m2, attr, None)

            # both present, almost equal
            setattr(m1, attr, pd.DataFrame([[1.0, 2.0], [3.0, 4.0]]))
            setattr(m2, attr, pd.DataFrame([[1.0, 2.0], [3.00000002, 4.0]]))
            assert_ordination_results_equal(m1, m2)
            setattr(m1, attr, None)
            setattr(m2, attr, None)

    def test_proportion_explained_mismatch(self):
        m2 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.0, 2.0]),
            pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        m2.proportion_explained = pd.Series([1.0, 2.0, 3.0])
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal(self.minimal1, m2)

    def test_proportion_explained_almost_equal(self):
        m1 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.0, 2.0]),
            pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        m2 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.0, 2.0]),
            pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        m1.proportion_explained = pd.Series([1.0, 2.0, 3.0])
        m2.proportion_explained = pd.Series([1.0, 2.0, 3.00000002])
        assert_ordination_results_equal(m1, m2)

    def test_ignore_method_names(self):
        m1 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.0, 2.0]),
            pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        m2 = OrdinationResults(
            'baz', 'qux',
            pd.Series([1.0, 2.0]),
            pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        # Should fail without ignore
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal(m1, m2)
        # Should pass with ignore
        assert_ordination_results_equal(m1, m2, ignore_method_names=True)

    def test_ignore_axis_labels(self):
        m1 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.0, 2.0], index=['a', 'b']),
            pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                         columns=['x', 'y', 'z']))
        m2 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.0, 2.0], index=['c', 'd']),
            pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                         columns=['p', 'q', 'r']))
        assert_ordination_results_equal(m1, m2, ignore_axis_labels=True)

    def test_ignore_directionality(self):
        m1 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.0, 2.0]),
            pd.DataFrame([[1.0, 2.0], [3.0, 4.0]]))
        m2 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.0, 2.0]),
            pd.DataFrame([[-1.0, 2.0], [-3.0, 4.0]]))
        # Pairwise distances should be the same even with sign flip
        assert_ordination_results_equal(m1, m2, ignore_directionality=True)

    def test_decimal_parameter(self):
        m1 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.0, 2.0]),
            pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        m2 = OrdinationResults(
            'foo', 'bar',
            pd.Series([1.001, 2.0]),
            pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]))
        # Should fail with high precision
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal(m1, m2, decimal=7)
        # Should pass with low precision
        assert_ordination_results_equal(m1, m2, decimal=2)


class TestAssertOrdinationResultsEqualNp(unittest.TestCase):
    """Tests for the NumPy version of ordination results comparison.

    ``assert_ordination_results_equal_np`` uses ``npt.assert_almost_equal``
    directly on ``eigvals`` and ``proportion_explained``, which requires
    NumPy arrays (not pd.Series) to avoid ambiguous truth-value errors in
    newer NumPy/Pandas. We therefore convert these attributes to raw arrays
    after constructing OrdinationResults.
    """

    @staticmethod
    def _make_ord(short='foo', long='bar',
                  eigvals=None, samples=None, **overrides):
        """Create an OrdinationResults with NumPy-array eigvals/prop_expl."""
        if eigvals is None:
            eigvals = pd.Series([1.0, 2.0])
        if samples is None:
            samples = pd.DataFrame([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        obj = OrdinationResults(short, long, eigvals, samples)
        # Convert to numpy arrays so npt.assert_almost_equal works
        obj.eigvals = obj.eigvals.values
        obj.samples = obj.samples.values
        for attr, val in overrides.items():
            if isinstance(val, pd.Series):
                setattr(obj, attr, val.values)
            elif isinstance(val, pd.DataFrame):
                setattr(obj, attr, val.values)
            else:
                setattr(obj, attr, val)
        return obj

    def test_basic_equal(self):
        m1 = self._make_ord()
        assert_ordination_results_equal_np(m1, m1)

    def test_ignore_method_names(self):
        m1 = self._make_ord('foo', 'bar')
        m2 = self._make_ord('baz', 'qux')
        with npt.assert_raises(AssertionError):
            assert_ordination_results_equal_np(m1, m2)
        assert_ordination_results_equal_np(m1, m2, ignore_method_names=True)

    def test_features_none(self):
        m1 = self._make_ord(features=None)
        m2 = self._make_ord(features=None)
        assert_ordination_results_equal_np(m1, m2)

    def test_features_present(self):
        feat = np.array([[1.0, 2.0], [3.0, 4.0]])
        m1 = self._make_ord(features=feat)
        m2 = self._make_ord(features=feat.copy())
        assert_ordination_results_equal_np(m1, m2)

    def test_biplot_scores_none(self):
        m1 = self._make_ord(biplot_scores=None)
        m2 = self._make_ord(biplot_scores=None)
        assert_ordination_results_equal_np(m1, m2)

    def test_biplot_scores_present(self):
        bp = np.array([[0.1, 0.2], [0.3, 0.4]])
        m1 = self._make_ord(biplot_scores=bp)
        m2 = self._make_ord(biplot_scores=bp.copy())
        assert_ordination_results_equal_np(m1, m2)

    def test_sample_constraints_none(self):
        m1 = self._make_ord(sample_constraints=None)
        m2 = self._make_ord(sample_constraints=None)
        assert_ordination_results_equal_np(m1, m2)

    def test_sample_constraints_present(self):
        sc = np.array([[1.0, 2.0], [3.0, 4.0]])
        m1 = self._make_ord(sample_constraints=sc)
        m2 = self._make_ord(sample_constraints=sc.copy())
        assert_ordination_results_equal_np(m1, m2)

    def test_proportion_explained_none(self):
        m1 = self._make_ord(proportion_explained=None)
        m2 = self._make_ord(proportion_explained=None)
        assert_ordination_results_equal_np(m1, m2)

    def test_proportion_explained_present(self):
        pe = np.array([0.5, 0.3])
        m1 = self._make_ord(proportion_explained=pe)
        m2 = self._make_ord(proportion_explained=pe.copy())
        assert_ordination_results_equal_np(m1, m2)


class TestAssertSeriesEqual(unittest.TestCase):
    """Tests for the private _assert_series_equal helper."""

    def test_both_none(self):
        _assert_series_equal(None, None, ignore_index=False)

    def test_one_none(self):
        with self.assertRaises(AssertionError):
            _assert_series_equal(pd.Series([1.0]), None, ignore_index=False)
        with self.assertRaises(AssertionError):
            _assert_series_equal(None, pd.Series([1.0]), ignore_index=False)

    def test_almost_equal(self):
        s1 = pd.Series([1.0, 2.0, 3.0])
        s2 = pd.Series([1.0, 2.00000001, 3.0])
        _assert_series_equal(s1, s2, ignore_index=False)

    def test_not_almost_equal(self):
        s1 = pd.Series([1.0, 2.0, 3.0])
        s2 = pd.Series([1.0, 2.1, 3.0])
        with self.assertRaises(AssertionError):
            _assert_series_equal(s1, s2, ignore_index=False)

    def test_ignore_index(self):
        s1 = pd.Series([1.0, 2.0], index=['a', 'b'])
        s2 = pd.Series([1.0, 2.0], index=['x', 'y'])
        _assert_series_equal(s1, s2, ignore_index=True)

    def test_different_index_fails(self):
        s1 = pd.Series([1.0, 2.0], index=['a', 'b'])
        s2 = pd.Series([1.0, 2.0], index=['x', 'y'])
        with self.assertRaises(AssertionError):
            _assert_series_equal(s1, s2, ignore_index=False)


class TestAssertFrameDistsEqual(unittest.TestCase):
    """Tests for the private _assert_frame_dists_equal helper."""

    def test_both_none(self):
        _assert_frame_dists_equal(None, None)

    def test_one_none(self):
        df = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]])
        with self.assertRaises(AssertionError):
            _assert_frame_dists_equal(df, None)
        with self.assertRaises(AssertionError):
            _assert_frame_dists_equal(None, df)

    def test_equal_dists(self):
        df = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]])
        _assert_frame_dists_equal(df, df)

    def test_different_dists(self):
        df1 = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]])
        df2 = pd.DataFrame([[1.0, 2.0], [10.0, 4.0]])
        with self.assertRaises(AssertionError):
            _assert_frame_dists_equal(df1, df2)

    def test_ignore_index(self):
        df1 = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]], index=['a', 'b'])
        df2 = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]], index=['x', 'y'])
        _assert_frame_dists_equal(df1, df2, ignore_index=True)
        with self.assertRaises(AssertionError):
            _assert_frame_dists_equal(df1, df2, ignore_index=False)

    def test_ignore_columns(self):
        df1 = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]], columns=['a', 'b'])
        df2 = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]], columns=['x', 'y'])
        _assert_frame_dists_equal(df1, df2, ignore_columns=True)
        with self.assertRaises(AssertionError):
            _assert_frame_dists_equal(df1, df2, ignore_columns=False)

    def test_ignore_directionality(self):
        df1 = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]])
        df2 = pd.DataFrame([[-1.0, 2.0], [-3.0, 4.0]])
        # pairwise dists are the same regardless of column sign flips
        _assert_frame_dists_equal(df1, df2, ignore_directionality=True)


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
            assert_data_frame_almost_equal(df, pd.DataFrame(df).copy())

        # every pair of dfs should not compare equal
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

    def test_rtol_parameter(self):
        df1 = pd.DataFrame({'a': [1.0, 2.0]})
        df2 = pd.DataFrame({'a': [1.01, 2.0]})
        with self.assertRaises(AssertionError):
            assert_data_frame_almost_equal(df1, df2, rtol=1e-5)
        assert_data_frame_almost_equal(df1, df2, rtol=0.1)


class TestDataFrameToDefaultIntType(unittest.TestCase):
    def test_converts_int_columns(self):
        df = pd.DataFrame({'a': np.array([1, 2, 3], dtype=np.int64),
                           'b': [1.0, 2.0, 3.0],
                           'c': ['x', 'y', 'z']})
        _data_frame_to_default_int_type(df)
        # After conversion, int columns should have platform default int type
        self.assertEqual(df['a'].dtype, np.dtype(int))
        # Non-int columns should be unchanged
        self.assertEqual(df['b'].dtype, np.float64)
        self.assertEqual(df['c'].dtype, object)

    def test_no_int_columns(self):
        df = pd.DataFrame({'a': [1.0, 2.0], 'b': ['x', 'y']})
        _data_frame_to_default_int_type(df)
        self.assertEqual(df['a'].dtype, np.float64)


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
            assert_series_almost_equal(s, pd.Series(s).copy())


class TestAssertIndexEqual(unittest.TestCase):
    def test_equal_indices(self):
        idx1 = pd.Index([1, 2, 3])
        idx2 = pd.Index([1, 2, 3])
        assert_index_equal(idx1, idx2)

    def test_different_values(self):
        idx1 = pd.Index([1, 2, 3])
        idx2 = pd.Index([1, 2, 4])
        with self.assertRaises(AssertionError):
            assert_index_equal(idx1, idx2)

    def test_different_types(self):
        idx1 = pd.Index([1, 2, 3])
        idx2 = pd.Index([1.0, 2.0, 3.0])
        with self.assertRaises(AssertionError):
            assert_index_equal(idx1, idx2)

    def test_different_names(self):
        idx1 = pd.Index([1, 2, 3], name='a')
        idx2 = pd.Index([1, 2, 3], name='b')
        with self.assertRaises(AssertionError):
            assert_index_equal(idx1, idx2)


class TestReadEnv(unittest.TestCase):
    def test_default_empty(self):
        with patch.dict(os.environ, {}, clear=True):
            os.environ.pop('SKBIO_ARRAY_API_BACKEND', None)
            os.environ.pop('SKBIO_DEVICE', None)
            backend, device = _read_env()
            self.assertEqual(backend, '')
            self.assertIsNone(device)

    def test_custom_values(self):
        with patch.dict(os.environ, {'SKBIO_ARRAY_API_BACKEND': ' jax ',
                                     'SKBIO_DEVICE': ' cuda '}):
            backend, device = _read_env()
            self.assertEqual(backend, 'jax')
            self.assertEqual(device, 'cuda')

    def test_backend_only(self):
        with patch.dict(os.environ, {'SKBIO_ARRAY_API_BACKEND': 'torch'},
                        clear=False):
            os.environ.pop('SKBIO_DEVICE', None)
            backend, device = _read_env()
            self.assertEqual(backend, 'torch')
            self.assertIsNone(device)


class TestGetAvailableBackends(unittest.TestCase):
    def test_numpy_always_present(self):
        backends_dict = _get_available_backends()
        self.assertIn('numpy', backends_dict)
        xp, devices = backends_dict['numpy']
        self.assertIs(xp, np)
        self.assertIn('cpu', devices)


class TestShouldRun(unittest.TestCase):
    def test_default_numpy_cpu(self):
        """With no env vars, only numpy/cpu should run."""
        with patch('skbio.util._testing._read_env', return_value=('', None)):
            self.assertTrue(_should_run('numpy', 'cpu'))
            self.assertFalse(_should_run('jax', 'cpu'))
            self.assertFalse(_should_run('numpy', 'gpu'))

    def test_all_backends(self):
        with patch('skbio.util._testing._read_env', return_value=('all', None)):
            self.assertTrue(_should_run('numpy', 'cpu'))
            self.assertTrue(_should_run('jax', 'cpu'))
            self.assertTrue(_should_run('torch', 'cuda'))

    def test_all_with_device_filter(self):
        with patch('skbio.util._testing._read_env', return_value=('all', 'cpu')):
            self.assertTrue(_should_run('numpy', 'cpu'))
            self.assertFalse(_should_run('torch', 'cuda'))

    def test_specific_backend(self):
        with patch('skbio.util._testing._read_env', return_value=('torch', None)):
            self.assertFalse(_should_run('numpy', 'cpu'))
            self.assertTrue(_should_run('torch', 'cpu'))
            self.assertTrue(_should_run('torch', 'cuda'))

    def test_specific_backend_and_device(self):
        with patch('skbio.util._testing._read_env', return_value=('torch', 'cuda')):
            self.assertFalse(_should_run('torch', 'cpu'))
            self.assertTrue(_should_run('torch', 'cuda'))


class TestXpAssertClose(unittest.TestCase):
    def test_close_arrays(self):
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([1.0, 2.0000001, 3.0])
        xp_assert_close(a, b)

    def test_not_close_arrays(self):
        a = np.array([1.0, 2.0, 3.0])
        b = np.array([1.0, 3.0, 3.0])
        with self.assertRaises(AssertionError):
            xp_assert_close(a, b)

    def test_custom_tolerances(self):
        a = np.array([1.0])
        b = np.array([1.1])
        with self.assertRaises(AssertionError):
            xp_assert_close(a, b, rtol=1e-7, atol=0)
        xp_assert_close(a, b, rtol=0.2, atol=0)


class TestXpAssertEqual(unittest.TestCase):
    def test_equal_arrays(self):
        a = np.array([1, 2, 3])
        b = np.array([1, 2, 3])
        xp_assert_equal(a, b)

    def test_not_equal_arrays(self):
        a = np.array([1, 2, 3])
        b = np.array([1, 2, 4])
        with self.assertRaises(AssertionError):
            xp_assert_equal(a, b)


class TestBackendsDecorator(unittest.TestCase):
    def test_runs_numpy_by_default(self):
        """With default env, the decorator should run the test on numpy/cpu."""
        ran = []

        class FakeTest(unittest.TestCase):
            @backends('numpy')
            def test_it(self, xp, device):
                ran.append((xp, device))

        with patch('skbio.util._testing._read_env', return_value=('', None)):
            t = FakeTest()
            t.test_it()
        self.assertEqual(len(ran), 1)
        self.assertIs(ran[0][0], np)
        self.assertEqual(ran[0][1], 'cpu')

    def test_missing_backend_raises(self):
        """Requesting a non-installed backend should raise RuntimeError."""
        class FakeTest(unittest.TestCase):
            @backends('numpy')
            def test_it(self, xp, device):
                pass

        with patch('skbio.util._testing._read_env',
                   return_value=('nonexistent', None)):
            with patch('skbio.util._testing._BACKENDS', {'numpy': (np, ['cpu'])}):
                t = FakeTest()
                with self.assertRaises(RuntimeError):
                    t.test_it()

    def test_no_matching_backends_skips(self):
        """If no backends match, should raise SkipTest."""
        class FakeTest(unittest.TestCase):
            @backends('jax')  # jax likely not installed in test env
            def test_it(self, xp, device):
                pass

        with patch('skbio.util._testing._read_env', return_value=('', None)):
            with patch('skbio.util._testing._BACKENDS', {'numpy': (np, ['cpu'])}):
                t = FakeTest()
                with self.assertRaises(unittest.SkipTest):
                    t.test_it()

    def test_cpu_only(self):
        """cpu_only=True should skip non-cpu devices."""
        ran = []

        class FakeTest(unittest.TestCase):
            @backends('numpy', cpu_only=True)
            def test_it(self, xp, device):
                ran.append(device)

        with patch('skbio.util._testing._read_env', return_value=('all', None)):
            with patch('skbio.util._testing._BACKENDS',
                       {'numpy': (np, ['cpu', 'cuda'])}):
                t = FakeTest()
                t.test_it()
        self.assertEqual(ran, ['cpu'])

    def test_default_all_backends(self):
        """No backend_names argument should default to all four."""
        ran = []

        class FakeTest(unittest.TestCase):
            @backends()
            def test_it(self, xp, device):
                ran.append('ok')

        with patch('skbio.util._testing._read_env', return_value=('all', None)):
            with patch('skbio.util._testing._BACKENDS',
                       {'numpy': (np, ['cpu'])}):
                with patch('skbio.util._testing._should_run', return_value=True):
                    t = FakeTest()
                    t.test_it()
        self.assertTrue(len(ran) >= 1)

    def test_no_matching_device_raises(self):
        """If SKBIO_DEVICE is set but no backend has it, RuntimeError."""
        class FakeTest(unittest.TestCase):
            @backends('numpy')
            def test_it(self, xp, device):
                pass

        with patch('skbio.util._testing._read_env',
                   return_value=('all', 'nonexistent')):
            with patch('skbio.util._testing._BACKENDS',
                       {'numpy': (np, ['cpu'])}):
                t = FakeTest()
                with self.assertRaises(RuntimeError):
                    t.test_it()


class TestArrayAPITestMixin(unittest.TestCase, ArrayAPITestMixin):
    def test_make_array(self):
        arr = self.make_array(np, 'cpu', [[1.0, 2.0], [3.0, 4.0]])
        npt.assert_array_equal(arr, np.array([[1.0, 2.0], [3.0, 4.0]]))

    def test_make_array_with_dtype(self):
        arr = self.make_array(np, 'cpu', [1, 2, 3], dtype=np.float32)
        self.assertEqual(arr.dtype, np.float32)

    def test_assert_close(self):
        a = np.array([1.0, 2.0])
        b = np.array([1.0, 2.0000001])
        self.assert_close(a, b)

    def test_assert_close_fails(self):
        a = np.array([1.0, 2.0])
        b = np.array([1.0, 3.0])
        with self.assertRaises(AssertionError):
            self.assert_close(a, b)

    def test_normalize_device(self):
        self.assertEqual(self._normalize_device('cuda:0'), 'gpu')
        self.assertEqual(self._normalize_device('gpu'), 'gpu')
        self.assertEqual(self._normalize_device('CUDA'), 'gpu')
        self.assertEqual(self._normalize_device('cpu'), 'cpu')


if __name__ == '__main__':
    unittest.main()