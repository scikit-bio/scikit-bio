# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.utils import PY3

import copy
import os
import inspect

import six
import pandas as pd
from nose import core
from nose.tools import nottest

import numpy as np
import numpy.testing as npt
import pandas.util.testing as pdt

from ._decorator import experimental


class ReallyEqualMixin(object):
    """Use this for testing __eq__/__ne__.

    Taken and modified from the following public domain code:
      https://ludios.org/testing-your-eq-ne-cmp/

    """

    def assertReallyEqual(self, a, b):
        # assertEqual first, because it will have a good message if the
        # assertion fails.
        self.assertEqual(a, b)
        self.assertEqual(b, a)
        self.assertTrue(a == b)
        self.assertTrue(b == a)
        self.assertFalse(a != b)
        self.assertFalse(b != a)

        # We do not support cmp/__cmp__ because they do not exist in Python 3.
        # However, we still test this to catch potential bugs where the
        # object's parent class defines a __cmp__.
        if not PY3:
            self.assertEqual(0, cmp(a, b))  # noqa
            self.assertEqual(0, cmp(b, a))  # noqa

    def assertReallyNotEqual(self, a, b):
        # assertNotEqual first, because it will have a good message if the
        # assertion fails.
        self.assertNotEqual(a, b)
        self.assertNotEqual(b, a)
        self.assertFalse(a == b)
        self.assertFalse(b == a)
        self.assertTrue(a != b)
        self.assertTrue(b != a)

        # We do not support cmp/__cmp__ because they do not exist in Python 3.
        # However, we still test this to catch potential bugs where the
        # object's parent class defines a __cmp__.
        if not PY3:
            self.assertNotEqual(0, cmp(a, b))  # noqa
            self.assertNotEqual(0, cmp(b, a))  # noqa


class MetadataMixinTests(object):
    def test_constructor_invalid_type(self):
        for md in (0, 'a', ('f', 'o', 'o'), np.array([]), pd.DataFrame()):
            with six.assertRaisesRegex(self, TypeError,
                                       'metadata must be a dict'):
                self._metadata_constructor_(metadata=md)

    def test_constructor_no_metadata(self):
        for md in None, {}:
            obj = self._metadata_constructor_(metadata=md)

            self.assertFalse(obj.has_metadata())
            self.assertEqual(obj.metadata, {})

    def test_constructor_with_metadata(self):
        obj = self._metadata_constructor_(metadata={'foo': 'bar'})
        self.assertEqual(obj.metadata, {'foo': 'bar'})

        obj = self._metadata_constructor_(
                metadata={'': '', 123: {'a': 'b', 'c': 'd'}})
        self.assertEqual(obj.metadata, {'': '', 123: {'a': 'b', 'c': 'd'}})

    def test_constructor_handles_missing_metadata_efficiently(self):
        self.assertIsNone(self._metadata_constructor_()._metadata)
        self.assertIsNone(self._metadata_constructor_(metadata=None)._metadata)

    def test_constructor_makes_shallow_copy_of_metadata(self):
        md = {'foo': 'bar', 42: []}
        obj = self._metadata_constructor_(metadata=md)

        self.assertEqual(obj.metadata, md)
        self.assertIsNot(obj.metadata, md)

        md['foo'] = 'baz'
        self.assertEqual(obj.metadata, {'foo': 'bar', 42: []})

        md[42].append(True)
        self.assertEqual(obj.metadata, {'foo': 'bar', 42: [True]})

    def test_eq(self):
        self.assertReallyEqual(
                self._metadata_constructor_(metadata={'foo': 42}),
                self._metadata_constructor_(metadata={'foo': 42}))

        self.assertReallyEqual(
                self._metadata_constructor_(metadata={'foo': 42, 123: {}}),
                self._metadata_constructor_(metadata={'foo': 42, 123: {}}))

    def test_eq_missing_metadata(self):
        self.assertReallyEqual(self._metadata_constructor_(),
                               self._metadata_constructor_())
        self.assertReallyEqual(self._metadata_constructor_(),
                               self._metadata_constructor_(metadata={}))
        self.assertReallyEqual(self._metadata_constructor_(metadata={}),
                               self._metadata_constructor_(metadata={}))

    def test_eq_handles_missing_metadata_efficiently(self):
        obj1 = self._metadata_constructor_()
        obj2 = self._metadata_constructor_()
        self.assertReallyEqual(obj1, obj2)

        self.assertIsNone(obj1._metadata)
        self.assertIsNone(obj2._metadata)

    def test_ne(self):
        # Both have metadata.
        obj1 = self._metadata_constructor_(metadata={'id': 'foo'})
        obj2 = self._metadata_constructor_(metadata={'id': 'bar'})
        self.assertReallyNotEqual(obj1, obj2)

        # One has metadata.
        obj1 = self._metadata_constructor_(metadata={'id': 'foo'})
        obj2 = self._metadata_constructor_()
        self.assertReallyNotEqual(obj1, obj2)

    def test_copy_metadata_none(self):
        obj = self._metadata_constructor_()
        obj_copy = copy.copy(obj)

        self.assertEqual(obj, obj_copy)
        self.assertIsNot(obj, obj_copy)

        self.assertIsNone(obj._metadata)
        self.assertIsNone(obj_copy._metadata)

    def test_copy_metadata_empty(self):
        obj = self._metadata_constructor_(metadata={})
        obj_copy = copy.copy(obj)

        self.assertEqual(obj, obj_copy)
        self.assertIsNot(obj, obj_copy)

        self.assertEqual(obj._metadata, {})
        self.assertIsNone(obj_copy._metadata)

    def test_copy_with_metadata(self):
        obj = self._metadata_constructor_(metadata={'foo': [1]})
        obj_copy = copy.copy(obj)

        self.assertEqual(obj, obj_copy)
        self.assertIsNot(obj, obj_copy)

        self.assertIsNot(obj._metadata, obj_copy._metadata)
        self.assertIs(obj._metadata['foo'], obj_copy._metadata['foo'])

        obj_copy.metadata['foo'].append(2)
        obj_copy.metadata['foo2'] = 42

        self.assertEqual(obj_copy.metadata, {'foo': [1, 2], 'foo2': 42})
        self.assertEqual(obj.metadata, {'foo': [1, 2]})

    def test_deepcopy_metadata_none(self):
        obj = self._metadata_constructor_()
        obj_copy = copy.deepcopy(obj)

        self.assertEqual(obj, obj_copy)
        self.assertIsNot(obj, obj_copy)

        self.assertIsNone(obj._metadata)
        self.assertIsNone(obj_copy._metadata)

    def test_deepcopy_metadata_empty(self):
        obj = self._metadata_constructor_(metadata={})
        obj_copy = copy.deepcopy(obj)

        self.assertEqual(obj, obj_copy)
        self.assertIsNot(obj, obj_copy)

        self.assertEqual(obj._metadata, {})
        self.assertIsNone(obj_copy._metadata)

    def test_deepcopy_with_metadata(self):
        obj = self._metadata_constructor_(metadata={'foo': [1]})
        obj_copy = copy.deepcopy(obj)

        self.assertEqual(obj, obj_copy)
        self.assertIsNot(obj, obj_copy)

        self.assertIsNot(obj._metadata, obj_copy._metadata)
        self.assertIsNot(obj._metadata['foo'], obj_copy._metadata['foo'])

        obj_copy.metadata['foo'].append(2)
        obj_copy.metadata['foo2'] = 42

        self.assertEqual(obj_copy.metadata, {'foo': [1, 2], 'foo2': 42})
        self.assertEqual(obj.metadata, {'foo': [1]})

    def test_deepcopy_memo_is_respected(self):
        # Basic test to ensure deepcopy's memo is passed through to recursive
        # deepcopy calls.
        obj = self._metadata_constructor_(metadata={'foo': 'bar'})
        memo = {}
        copy.deepcopy(obj, memo)
        self.assertGreater(len(memo), 2)

    def test_metadata_getter(self):
        obj = self._metadata_constructor_(
                metadata={42: 'foo', ('hello', 'world'): 43})

        self.assertIsInstance(obj.metadata, dict)
        self.assertEqual(obj.metadata, {42: 'foo', ('hello', 'world'): 43})

        obj.metadata[42] = 'bar'
        self.assertEqual(obj.metadata, {42: 'bar', ('hello', 'world'): 43})

    def test_metadata_getter_no_metadata(self):
        obj = self._metadata_constructor_()

        self.assertIsNone(obj._metadata)
        self.assertIsInstance(obj.metadata, dict)
        self.assertEqual(obj.metadata, {})
        self.assertIsNotNone(obj._metadata)

    def test_metadata_setter(self):
        obj = self._metadata_constructor_()

        self.assertFalse(obj.has_metadata())

        obj.metadata = {'hello': 'world'}
        self.assertTrue(obj.has_metadata())
        self.assertEqual(obj.metadata, {'hello': 'world'})

        obj.metadata = {}
        self.assertFalse(obj.has_metadata())
        self.assertEqual(obj.metadata, {})

    def test_metadata_setter_makes_shallow_copy(self):
        obj = self._metadata_constructor_()

        md = {'foo': 'bar', 42: []}
        obj.metadata = md

        self.assertEqual(obj.metadata, md)
        self.assertIsNot(obj.metadata, md)

        md['foo'] = 'baz'
        self.assertEqual(obj.metadata, {'foo': 'bar', 42: []})

        md[42].append(True)
        self.assertEqual(obj.metadata, {'foo': 'bar', 42: [True]})

    def test_metadata_setter_invalid_type(self):
        obj = self._metadata_constructor_(metadata={123: 456})

        for md in (None, 0, 'a', ('f', 'o', 'o'), np.array([]),
                   pd.DataFrame()):
            with six.assertRaisesRegex(self, TypeError,
                                       'metadata must be a dict'):
                obj.metadata = md
            self.assertEqual(obj.metadata, {123: 456})

    def test_metadata_deleter(self):
        obj = self._metadata_constructor_(metadata={'foo': 'bar'})

        self.assertEqual(obj.metadata, {'foo': 'bar'})

        del obj.metadata
        self.assertIsNone(obj._metadata)
        self.assertFalse(obj.has_metadata())

        # Delete again.
        del obj.metadata
        self.assertIsNone(obj._metadata)
        self.assertFalse(obj.has_metadata())

        obj = self._metadata_constructor_()

        self.assertIsNone(obj._metadata)
        self.assertFalse(obj.has_metadata())
        del obj.metadata
        self.assertIsNone(obj._metadata)
        self.assertFalse(obj.has_metadata())

    def test_has_metadata(self):
        obj = self._metadata_constructor_()

        self.assertFalse(obj.has_metadata())
        # Handles metadata efficiently.
        self.assertIsNone(obj._metadata)

        self.assertFalse(
                self._metadata_constructor_(metadata={}).has_metadata())

        self.assertTrue(
                self._metadata_constructor_(metadata={'': ''}).has_metadata())
        self.assertTrue(
                self._metadata_constructor_(
                        metadata={'foo': 42}).has_metadata())


class PositionalMetadataMixinTests(object):
    def test_constructor_invalid_positional_metadata_type(self):
        with six.assertRaisesRegex(self, TypeError,
                                   'Invalid positional metadata. Must be '
                                   'consumable by `pd.DataFrame` constructor. '
                                   'Original pandas error message: '):
            self._positional_metadata_constructor_(0, positional_metadata=2)

    def test_constructor_positional_metadata_len_mismatch(self):
        # Zero elements.
        with six.assertRaisesRegex(self, ValueError, '\(0\).*\(4\)'):
            self._positional_metadata_constructor_(4, positional_metadata=[])

        # Not enough elements.
        with six.assertRaisesRegex(self, ValueError, '\(3\).*\(4\)'):
            self._positional_metadata_constructor_(
                4, positional_metadata=[2, 3, 4])

        # Too many elements.
        with six.assertRaisesRegex(self, ValueError, '\(5\).*\(4\)'):
            self._positional_metadata_constructor_(
                4, positional_metadata=[2, 3, 4, 5, 6])

        # Series not enough rows.
        with six.assertRaisesRegex(self, ValueError, '\(3\).*\(4\)'):
            self._positional_metadata_constructor_(
                4, positional_metadata=pd.Series(range(3)))

        # Series too many rows.
        with six.assertRaisesRegex(self, ValueError, '\(5\).*\(4\)'):
            self._positional_metadata_constructor_(
                4, positional_metadata=pd.Series(range(5)))

        # DataFrame not enough rows.
        with six.assertRaisesRegex(self, ValueError, '\(3\).*\(4\)'):
            self._positional_metadata_constructor_(
                4, positional_metadata=pd.DataFrame({'quality': range(3)}))

        # DataFrame too many rows.
        with six.assertRaisesRegex(self, ValueError, '\(5\).*\(4\)'):
            self._positional_metadata_constructor_(
                4, positional_metadata=pd.DataFrame({'quality': range(5)}))

    def test_constructor_no_positional_metadata(self):
        # Length zero with missing/empty positional metadata.
        for empty in None, {}, pd.DataFrame():
            obj = self._positional_metadata_constructor_(
                0, positional_metadata=empty)

            self.assertFalse(obj.has_positional_metadata())
            assert_data_frame_almost_equal(obj.positional_metadata,
                                           pd.DataFrame(index=np.arange(0)))

        # Nonzero length with missing positional metadata.
        obj = self._positional_metadata_constructor_(
            3, positional_metadata=None)

        self.assertFalse(obj.has_positional_metadata())
        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame(index=np.arange(3)))

    def test_constructor_with_positional_metadata_len_zero(self):
        for data in [], (), np.array([]):
            obj = self._positional_metadata_constructor_(
                0, positional_metadata={'foo': data})

            self.assertTrue(obj.has_positional_metadata())
            assert_data_frame_almost_equal(
                obj.positional_metadata,
                pd.DataFrame({'foo': data}, index=np.arange(0)))

    def test_constructor_with_positional_metadata_len_one(self):
        for data in [2], (2, ), np.array([2]):
            obj = self._positional_metadata_constructor_(
                1, positional_metadata={'foo': data})

            self.assertTrue(obj.has_positional_metadata())
            assert_data_frame_almost_equal(
                obj.positional_metadata,
                pd.DataFrame({'foo': data}, index=np.arange(1)))

    def test_constructor_with_positional_metadata_len_greater_than_one(self):
        for data in ([0, 42, 42, 1, 0, 8, 100, 0, 0],
                     (0, 42, 42, 1, 0, 8, 100, 0, 0),
                     np.array([0, 42, 42, 1, 0, 8, 100, 0, 0])):
            obj = self._positional_metadata_constructor_(
                9, positional_metadata={'foo': data})

            self.assertTrue(obj.has_positional_metadata())
            assert_data_frame_almost_equal(
                obj.positional_metadata,
                pd.DataFrame({'foo': data}, index=np.arange(9)))

    def test_constructor_with_positional_metadata_multiple_columns(self):
        obj = self._positional_metadata_constructor_(
            5, positional_metadata={'foo': np.arange(5),
                                    'bar': np.arange(5)[::-1]})

        self.assertTrue(obj.has_positional_metadata())
        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'foo': np.arange(5),
                          'bar': np.arange(5)[::-1]}, index=np.arange(5)))

    def test_constructor_with_positional_metadata_custom_index(self):
        df = pd.DataFrame({'foo': np.arange(5), 'bar': np.arange(5)[::-1]},
                          index=['a', 'b', 'c', 'd', 'e'])
        obj = self._positional_metadata_constructor_(
            5, positional_metadata=df)

        self.assertTrue(obj.has_positional_metadata())
        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'foo': np.arange(5),
                          'bar': np.arange(5)[::-1]}, index=np.arange(5)))

    def test_constructor_handles_missing_positional_metadata_efficiently(self):
        obj = self._positional_metadata_constructor_(4)
        self.assertIsNone(obj._positional_metadata)

        obj = self._positional_metadata_constructor_(
            4, positional_metadata=None)
        self.assertIsNone(obj._positional_metadata)

    def test_constructor_makes_shallow_copy_of_positional_metadata(self):
        df = pd.DataFrame({'foo': [22, 22, 0], 'bar': [[], [], []]},
                          index=['a', 'b', 'c'])
        obj = self._positional_metadata_constructor_(
            3, positional_metadata=df)

        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'foo': [22, 22, 0], 'bar': [[], [], []]},
                         index=np.arange(3)))
        self.assertIsNot(obj.positional_metadata, df)

        # Original df is not mutated.
        orig_df = pd.DataFrame({'foo': [22, 22, 0], 'bar': [[], [], []]},
                               index=['a', 'b', 'c'])
        assert_data_frame_almost_equal(df, orig_df)

        # Change values of column (using same dtype).
        df['foo'] = [42, 42, 42]
        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'foo': [22, 22, 0], 'bar': [[], [], []]},
                         index=np.arange(3)))

        # Change single value of underlying data.
        df.values[0][0] = 10
        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'foo': [22, 22, 0], 'bar': [[], [], []]},
                         index=np.arange(3)))

        # Mutate list (not a deep copy).
        df['bar'][0].append(42)
        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'foo': [22, 22, 0], 'bar': [[42], [], []]},
                         index=np.arange(3)))

    def test_eq_basic(self):
        obj1 = self._positional_metadata_constructor_(
            3, positional_metadata={'foo': [1, 2, 3]})
        obj2 = self._positional_metadata_constructor_(
            3, positional_metadata={'foo': [1, 2, 3]})
        self.assertReallyEqual(obj1, obj2)

    def test_eq_from_different_source(self):
        obj1 = self._positional_metadata_constructor_(
            3, positional_metadata={'foo': np.array([1, 2, 3])})
        obj2 = self._positional_metadata_constructor_(
            3, positional_metadata=pd.DataFrame({'foo': [1, 2, 3]},
                                                index=['foo', 'bar', 'baz']))
        self.assertReallyEqual(obj1, obj2)

    def test_eq_missing_positional_metadata(self):
        for empty in None, {}, pd.DataFrame(), pd.DataFrame(index=[]):
            obj = self._positional_metadata_constructor_(
                0, positional_metadata=empty)

            self.assertReallyEqual(
                obj,
                self._positional_metadata_constructor_(0))
            self.assertReallyEqual(
                obj,
                self._positional_metadata_constructor_(
                    0, positional_metadata=empty))

        for empty in None, pd.DataFrame(index=['a', 'b']):
            obj = self._positional_metadata_constructor_(
                2, positional_metadata=empty)

            self.assertReallyEqual(
                obj,
                self._positional_metadata_constructor_(2))
            self.assertReallyEqual(
                obj,
                self._positional_metadata_constructor_(
                    2, positional_metadata=empty))

    def test_eq_handles_missing_positional_metadata_efficiently(self):
        obj1 = self._positional_metadata_constructor_(1)
        obj2 = self._positional_metadata_constructor_(1)
        self.assertReallyEqual(obj1, obj2)

        self.assertIsNone(obj1._positional_metadata)
        self.assertIsNone(obj2._positional_metadata)

    def test_ne_len_zero(self):
        # Both have positional metadata.
        obj1 = self._positional_metadata_constructor_(
            0, positional_metadata={'foo': []})
        obj2 = self._positional_metadata_constructor_(
            0, positional_metadata={'foo': [], 'bar': []})
        self.assertReallyNotEqual(obj1, obj2)

        # One has positional metadata.
        obj1 = self._positional_metadata_constructor_(
            0, positional_metadata={'foo': []})
        obj2 = self._positional_metadata_constructor_(0)
        self.assertReallyNotEqual(obj1, obj2)

    def test_ne_len_greater_than_zero(self):
        # Both have positional metadata.
        obj1 = self._positional_metadata_constructor_(
            3, positional_metadata={'foo': [1, 2, 3]})
        obj2 = self._positional_metadata_constructor_(
            3, positional_metadata={'foo': [1, 2, 2]})
        self.assertReallyNotEqual(obj1, obj2)

        # One has positional metadata.
        obj1 = self._positional_metadata_constructor_(
            3, positional_metadata={'foo': [1, 2, 3]})
        obj2 = self._positional_metadata_constructor_(3)
        self.assertReallyNotEqual(obj1, obj2)

    def test_copy_positional_metadata_none(self):
        obj = self._positional_metadata_constructor_(3)
        obj_copy = copy.copy(obj)

        self.assertEqual(obj, obj_copy)
        self.assertIsNot(obj, obj_copy)

        self.assertIsNone(obj._positional_metadata)
        self.assertIsNone(obj_copy._positional_metadata)

    def test_copy_positional_metadata_empty(self):
        obj = self._positional_metadata_constructor_(
            3, positional_metadata=pd.DataFrame(index=range(3)))
        obj_copy = copy.copy(obj)

        self.assertEqual(obj, obj_copy)
        self.assertIsNot(obj, obj_copy)

        assert_data_frame_almost_equal(obj._positional_metadata,
                                       pd.DataFrame(index=range(3)))
        self.assertIsNone(obj_copy._positional_metadata)

    def test_copy_with_positional_metadata(self):
        obj = self._positional_metadata_constructor_(
            4, positional_metadata={'bar': [[], [], [], []],
                                    'baz': [42, 42, 42, 42]})
        obj_copy = copy.copy(obj)

        self.assertEqual(obj, obj_copy)
        self.assertIsNot(obj, obj_copy)

        self.assertIsNot(obj._positional_metadata,
                         obj_copy._positional_metadata)
        self.assertIsNot(obj._positional_metadata.values,
                         obj_copy._positional_metadata.values)
        self.assertIs(obj._positional_metadata.loc[0, 'bar'],
                      obj_copy._positional_metadata.loc[0, 'bar'])

        obj_copy.positional_metadata.loc[0, 'bar'].append(1)
        obj_copy.positional_metadata.loc[0, 'baz'] = 43

        assert_data_frame_almost_equal(
            obj_copy.positional_metadata,
            pd.DataFrame({'bar': [[1], [], [], []],
                          'baz': [43, 42, 42, 42]}))
        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'bar': [[1], [], [], []],
                          'baz': [42, 42, 42, 42]}))

    def test_deepcopy_positional_metadata_none(self):
        obj = self._positional_metadata_constructor_(3)
        obj_copy = copy.deepcopy(obj)

        self.assertEqual(obj, obj_copy)
        self.assertIsNot(obj, obj_copy)

        self.assertIsNone(obj._positional_metadata)
        self.assertIsNone(obj_copy._positional_metadata)

    def test_deepcopy_positional_metadata_empty(self):
        obj = self._positional_metadata_constructor_(
            3, positional_metadata=pd.DataFrame(index=range(3)))
        obj_copy = copy.deepcopy(obj)

        self.assertEqual(obj, obj_copy)
        self.assertIsNot(obj, obj_copy)

        assert_data_frame_almost_equal(obj._positional_metadata,
                                       pd.DataFrame(index=range(3)))
        self.assertIsNone(obj_copy._positional_metadata)

    def test_deepcopy_with_positional_metadata(self):
        obj = self._positional_metadata_constructor_(
            4, positional_metadata={'bar': [[], [], [], []],
                                    'baz': [42, 42, 42, 42]})
        obj_copy = copy.deepcopy(obj)

        self.assertEqual(obj, obj_copy)
        self.assertIsNot(obj, obj_copy)

        self.assertIsNot(obj._positional_metadata,
                         obj_copy._positional_metadata)
        self.assertIsNot(obj._positional_metadata.values,
                         obj_copy._positional_metadata.values)
        self.assertIsNot(obj._positional_metadata.loc[0, 'bar'],
                         obj_copy._positional_metadata.loc[0, 'bar'])

        obj_copy.positional_metadata.loc[0, 'bar'].append(1)
        obj_copy.positional_metadata.loc[0, 'baz'] = 43

        assert_data_frame_almost_equal(
            obj_copy.positional_metadata,
            pd.DataFrame({'bar': [[1], [], [], []],
                          'baz': [43, 42, 42, 42]}))
        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'bar': [[], [], [], []],
                          'baz': [42, 42, 42, 42]}))

    def test_deepcopy_memo_is_respected(self):
        # Basic test to ensure deepcopy's memo is passed through to recursive
        # deepcopy calls.
        obj = self._positional_metadata_constructor_(
            3, positional_metadata={'foo': [1, 2, 3]})
        memo = {}
        copy.deepcopy(obj, memo)
        self.assertGreater(len(memo), 2)

    def test_positional_metadata_getter(self):
        obj = self._positional_metadata_constructor_(
            3, positional_metadata={'foo': [22, 22, 0]})

        self.assertIsInstance(obj.positional_metadata, pd.DataFrame)
        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame({'foo': [22, 22, 0]}))

        # Update existing column.
        obj.positional_metadata['foo'] = [42, 42, 43]
        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame({'foo': [42, 42, 43]}))

        # Add new column.
        obj.positional_metadata['foo2'] = [True, False, True]
        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'foo': [42, 42, 43],
                          'foo2': [True, False, True]}))

    def test_positional_metadata_getter_no_positional_metadata(self):
        obj = self._positional_metadata_constructor_(4)

        self.assertIsNone(obj._positional_metadata)
        self.assertIsInstance(obj.positional_metadata, pd.DataFrame)
        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame(index=np.arange(4)))
        self.assertIsNotNone(obj._positional_metadata)

    def test_positional_metadata_getter_set_column_series(self):
        length = 8
        obj = self._positional_metadata_constructor_(
            length, positional_metadata={'foo': range(length)})

        obj.positional_metadata['bar'] = pd.Series(range(length-3))
        # pandas.Series will be padded with NaN if too short.
        npt.assert_equal(obj.positional_metadata['bar'],
                         np.array(list(range(length-3)) + [np.nan]*3))

        obj.positional_metadata['baz'] = pd.Series(range(length+3))
        # pandas.Series will be truncated if too long.
        npt.assert_equal(obj.positional_metadata['baz'],
                         np.array(range(length)))

    def test_positional_metadata_getter_set_column_array(self):
        length = 8
        obj = self._positional_metadata_constructor_(
            length, positional_metadata={'foo': range(length)})

        # array-like objects will fail if wrong size.
        for array_like in (np.array(range(length-1)), range(length-1),
                           np.array(range(length+1)), range(length+1)):
            with six.assertRaisesRegex(self, ValueError,
                                       "Length of values does not match "
                                       "length of index"):
                obj.positional_metadata['bar'] = array_like

    def test_positional_metadata_setter_pandas_consumable(self):
        obj = self._positional_metadata_constructor_(3)

        self.assertFalse(obj.has_positional_metadata())

        obj.positional_metadata = {'foo': [3, 2, 1]}
        self.assertTrue(obj.has_positional_metadata())
        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame({'foo': [3, 2, 1]}))

        obj.positional_metadata = pd.DataFrame(index=np.arange(3))
        self.assertFalse(obj.has_positional_metadata())
        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame(index=np.arange(3)))

    def test_positional_metadata_setter_data_frame(self):
        obj = self._positional_metadata_constructor_(3)

        self.assertFalse(obj.has_positional_metadata())

        obj.positional_metadata = pd.DataFrame({'foo': [3, 2, 1]},
                                               index=['a', 'b', 'c'])
        self.assertTrue(obj.has_positional_metadata())
        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame({'foo': [3, 2, 1]}))

        obj.positional_metadata = pd.DataFrame(index=np.arange(3))
        self.assertFalse(obj.has_positional_metadata())
        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame(index=np.arange(3)))

    def test_positional_metadata_setter_none(self):
        obj = self._positional_metadata_constructor_(
            0, positional_metadata={'foo': []})

        self.assertTrue(obj.has_positional_metadata())
        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame({'foo': []}))

        # `None` behavior differs from constructor.
        obj.positional_metadata = None

        self.assertFalse(obj.has_positional_metadata())
        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame(index=np.arange(0)))

    def test_positional_metadata_setter_makes_shallow_copy(self):
        obj = self._positional_metadata_constructor_(3)

        df = pd.DataFrame({'foo': [22, 22, 0], 'bar': [[], [], []]},
                          index=['a', 'b', 'c'])
        obj.positional_metadata = df

        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'foo': [22, 22, 0], 'bar': [[], [], []]},
                         index=np.arange(3)))
        self.assertIsNot(obj.positional_metadata, df)

        # Original df is not mutated.
        orig_df = pd.DataFrame({'foo': [22, 22, 0], 'bar': [[], [], []]},
                               index=['a', 'b', 'c'])
        assert_data_frame_almost_equal(df, orig_df)

        # Change values of column (using same dtype).
        df['foo'] = [42, 42, 42]
        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'foo': [22, 22, 0], 'bar': [[], [], []]},
                         index=np.arange(3)))

        # Change single value of underlying data.
        df.values[0][0] = 10
        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'foo': [22, 22, 0], 'bar': [[], [], []]},
                         index=np.arange(3)))

        # Mutate list (not a deep copy).
        df['bar'][0].append(42)
        assert_data_frame_almost_equal(
            obj.positional_metadata,
            pd.DataFrame({'foo': [22, 22, 0], 'bar': [[42], [], []]},
                         index=np.arange(3)))

    def test_positional_metadata_setter_invalid_type(self):
        obj = self._positional_metadata_constructor_(
            3, positional_metadata={'foo': [1, 2, 42]})

        with six.assertRaisesRegex(self, TypeError,
                                   'Invalid positional metadata. Must be '
                                   'consumable by `pd.DataFrame` constructor. '
                                   'Original pandas error message: '):
            obj.positional_metadata = 2

        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame({'foo': [1, 2, 42]}))

    def test_positional_metadata_setter_len_mismatch(self):
        obj = self._positional_metadata_constructor_(
            3, positional_metadata={'foo': [1, 2, 42]})

        # `None` behavior differs from constructor.
        with six.assertRaisesRegex(self, ValueError, '\(0\).*\(3\)'):
            obj.positional_metadata = None

        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame({'foo': [1, 2, 42]}))

        with six.assertRaisesRegex(self, ValueError, '\(4\).*\(3\)'):
            obj.positional_metadata = [1, 2, 3, 4]

        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame({'foo': [1, 2, 42]}))

    def test_positional_metadata_deleter(self):
        obj = self._positional_metadata_constructor_(
            3, positional_metadata={'foo': [1, 2, 3]})

        assert_data_frame_almost_equal(obj.positional_metadata,
                                       pd.DataFrame({'foo': [1, 2, 3]}))

        del obj.positional_metadata
        self.assertIsNone(obj._positional_metadata)
        self.assertFalse(obj.has_positional_metadata())

        # Delete again.
        del obj.positional_metadata
        self.assertIsNone(obj._positional_metadata)
        self.assertFalse(obj.has_positional_metadata())

        obj = self._positional_metadata_constructor_(3)

        self.assertIsNone(obj._positional_metadata)
        self.assertFalse(obj.has_positional_metadata())
        del obj.positional_metadata
        self.assertIsNone(obj._positional_metadata)
        self.assertFalse(obj.has_positional_metadata())

    def test_has_positional_metadata(self):
        obj = self._positional_metadata_constructor_(4)
        self.assertFalse(obj.has_positional_metadata())
        self.assertIsNone(obj._positional_metadata)

        obj = self._positional_metadata_constructor_(0, positional_metadata={})
        self.assertFalse(obj.has_positional_metadata())

        obj = self._positional_metadata_constructor_(
            4, positional_metadata=pd.DataFrame(index=np.arange(4)))
        self.assertFalse(obj.has_positional_metadata())

        obj = self._positional_metadata_constructor_(
            4, positional_metadata=pd.DataFrame(index=['a', 'b', 'c', 'd']))
        self.assertFalse(obj.has_positional_metadata())

        obj = self._positional_metadata_constructor_(
            0, positional_metadata={'foo': []})
        self.assertTrue(obj.has_positional_metadata())

        obj = self._positional_metadata_constructor_(
            4, positional_metadata={'foo': [1, 2, 3, 4]})
        self.assertTrue(obj.has_positional_metadata())

        obj = self._positional_metadata_constructor_(
            2, positional_metadata={'foo': [1, 2], 'bar': ['abc', 'def']})
        self.assertTrue(obj.has_positional_metadata())


@nottest
class TestRunner(object):
    """Simple wrapper class around nosetests functionality.

    Parameters
    ----------
    filename : str
        __file__ attribute passed in from the caller. This tells the
        tester where to start looking for tests.

    Notes
    -----
    The primary purpose of this class is to create an interface which users
    of scikit-bio can use to run all of the built in tests. Normally this
    would be done by invoking nosetests directly from the command line, but
    scikit-bio needs several additional options which make the command long
    and ugly. This class invokes nose with the required options.

    """
    @experimental(as_of="0.4.0")
    def __init__(self, filename):
        self._filename = filename
        self._test_dir = os.path.dirname(filename)

    @experimental(as_of="0.4.0")
    def test(self, verbose=False):
        """Performs the actual running of the tests.

        Parameters
        ----------
        verbose : bool
            flag for running in verbose mode.

        Returns
        -------
        bool
            test run success status
        """
        # NOTE: it doesn't seem to matter what the first element of the argv
        # list is, there just needs to be something there.
        argv = [self._filename, '-I DO_NOT_IGNORE_ANYTHING']
        if PY3:
            argv.extend(['--with-doctest', '--doctest-tests'])
        if verbose:
            argv.append('-v')
        return core.run(argv=argv, defaultTest=self._test_dir)


@experimental(as_of="0.4.0")
def get_data_path(fn, subfolder='data'):
    """Return path to filename ``fn`` in the data folder.

    During testing it is often necessary to load data files. This
    function returns the full path to files in the ``data`` subfolder
    by default.

    Parameters
    ----------
    fn : str
        File name.

    subfolder : str, defaults to ``data``
        Name of the subfolder that contains the data.


    Returns
    -------
    str
        Inferred absolute path to the test data for the module where
        ``get_data_path(fn)`` is called.

    Notes
    -----
    The requested path may not point to an existing file, as its
    existence is not checked.

    """
    # getouterframes returns a list of tuples: the second tuple
    # contains info about the caller, and the second element is its
    # filename
    callers_filename = inspect.getouterframes(inspect.currentframe())[1][1]
    path = os.path.dirname(os.path.abspath(callers_filename))
    data_path = os.path.join(path, subfolder, fn)
    return data_path


@experimental(as_of="0.4.0")
def assert_ordination_results_equal(left, right, ignore_method_names=False,
                                    ignore_axis_labels=False,
                                    ignore_biplot_scores_labels=False,
                                    ignore_directionality=False,
                                    decimal=7):
    """Assert that ordination results objects are equal.

    This is a helper function intended to be used in unit tests that need to
    compare ``OrdinationResults`` objects.

    Parameters
    ----------
    left, right : OrdinationResults
        Ordination results to be compared for equality.
    ignore_method_names : bool, optional
        Ignore differences in `short_method_name` and `long_method_name`.
    ignore_axis_labels : bool, optional
        Ignore differences in axis labels (i.e., column labels).
    ignore_biplot_scores_labels : bool, optional
        Ignore differences in `biplot_scores` row and column labels.
    ignore_directionality : bool, optional
        Ignore differences in directionality (i.e., differences in signs) for
        attributes `samples` and `features`.

    Raises
    ------
    AssertionError
        If the two objects are not equal.

    """
    npt.assert_equal(type(left) is type(right), True)

    if not ignore_method_names:
        npt.assert_equal(left.short_method_name, right.short_method_name)
        npt.assert_equal(left.long_method_name, right.long_method_name)

    _assert_frame_equal(left.samples, right.samples,
                        ignore_columns=ignore_axis_labels,
                        ignore_directionality=ignore_directionality,
                        decimal=decimal)

    _assert_frame_equal(left.features, right.features,
                        ignore_columns=ignore_axis_labels,
                        ignore_directionality=ignore_directionality,
                        decimal=decimal)

    _assert_frame_equal(left.biplot_scores, right.biplot_scores,
                        ignore_biplot_scores_labels,
                        ignore_biplot_scores_labels,
                        decimal=decimal)

    _assert_frame_equal(left.sample_constraints, right.sample_constraints,
                        ignore_columns=ignore_axis_labels,
                        decimal=decimal)

    _assert_series_equal(left.eigvals, right.eigvals, ignore_axis_labels,
                         decimal=decimal)

    _assert_series_equal(left.proportion_explained, right.proportion_explained,
                         ignore_axis_labels,
                         decimal=decimal)


def _assert_series_equal(left_s, right_s, ignore_index=False, decimal=7):
    # assert_series_equal doesn't like None...
    if left_s is None or right_s is None:
        assert left_s is None and right_s is None
    else:
        npt.assert_almost_equal(left_s.values, right_s.values,
                                decimal=decimal)
        if not ignore_index:
            pdt.assert_index_equal(left_s.index, right_s.index)


def _assert_frame_equal(left_df, right_df, ignore_index=False,
                        ignore_columns=False, ignore_directionality=False,
                        decimal=7):
    # assert_frame_equal doesn't like None...
    if left_df is None or right_df is None:
        assert left_df is None and right_df is None
    else:
        left_values = left_df.values
        right_values = right_df.values

        if ignore_directionality:
            left_values, right_values = _normalize_signs(left_values,
                                                         right_values)
        npt.assert_almost_equal(left_values, right_values, decimal=decimal)

        if not ignore_index:
            pdt.assert_index_equal(left_df.index, right_df.index)
        if not ignore_columns:
            pdt.assert_index_equal(left_df.columns, right_df.columns)


def _normalize_signs(arr1, arr2):
    """Change column signs so that "column" and "-column" compare equal.

    This is needed because results of eigenproblmes can have signs
    flipped, but they're still right.

    Notes
    =====

    This function tries hard to make sure that, if you find "column"
    and "-column" almost equal, calling a function like np.allclose to
    compare them after calling `normalize_signs` succeeds.

    To do so, it distinguishes two cases for every column:

    - It can be all almost equal to 0 (this includes a column of
      zeros).
    - Otherwise, it has a value that isn't close to 0.

    In the first case, no sign needs to be flipped. I.e., for
    |epsilon| small, np.allclose(-epsilon, 0) is true if and only if
    np.allclose(epsilon, 0) is.

    In the second case, the function finds the number in the column
    whose absolute value is largest. Then, it compares its sign with
    the number found in the same index, but in the other array, and
    flips the sign of the column as needed.
    """
    # Let's convert everyting to floating point numbers (it's
    # reasonable to assume that eigenvectors will already be floating
    # point numbers). This is necessary because np.array(1) /
    # np.array(0) != np.array(1.) / np.array(0.)
    arr1 = np.asarray(arr1, dtype=np.float64)
    arr2 = np.asarray(arr2, dtype=np.float64)

    if arr1.shape != arr2.shape:
        raise ValueError(
            "Arrays must have the same shape ({0} vs {1}).".format(arr1.shape,
                                                                   arr2.shape)
            )

    # To avoid issues around zero, we'll compare signs of the values
    # with highest absolute value
    max_idx = np.abs(arr1).argmax(axis=0)
    max_arr1 = arr1[max_idx, range(arr1.shape[1])]
    max_arr2 = arr2[max_idx, range(arr2.shape[1])]

    sign_arr1 = np.sign(max_arr1)
    sign_arr2 = np.sign(max_arr2)

    # Store current warnings, and ignore division by zero (like 1. /
    # 0.) and invalid operations (like 0. / 0.)
    wrn = np.seterr(invalid='ignore', divide='ignore')
    differences = sign_arr1 / sign_arr2
    # The values in `differences` can be:
    #    1 -> equal signs
    #   -1 -> diff signs
    #   Or nan (0/0), inf (nonzero/0), 0 (0/nonzero)
    np.seterr(**wrn)

    # Now let's deal with cases where `differences != \pm 1`
    special_cases = (~np.isfinite(differences)) | (differences == 0)
    # In any of these cases, the sign of the column doesn't matter, so
    # let's just keep it
    differences[special_cases] = 1

    return arr1 * differences, arr2


@experimental(as_of="0.4.0")
def assert_data_frame_almost_equal(left, right):
    """Raise AssertionError if ``pd.DataFrame`` objects are not "almost equal".

    Wrapper of ``pd.util.testing.assert_frame_equal``. Floating point values
    are considered "almost equal" if they are within a threshold defined by
    ``assert_frame_equal``. This wrapper uses a number of
    checks that are turned off by default in ``assert_frame_equal`` in order to
    perform stricter comparisons (for example, ensuring the index and column
    types are the same). It also does not consider empty ``pd.DataFrame``
    objects equal if they have a different index.

    Other notes:

    * Index (row) and column ordering must be the same for objects to be equal.
    * NaNs (``np.nan``) in the same locations are considered equal.

    This is a helper function intended to be used in unit tests that need to
    compare ``pd.DataFrame`` objects.

    Parameters
    ----------
    left, right : pd.DataFrame
        ``pd.DataFrame`` objects to compare.

    Raises
    ------
    AssertionError
        If `left` and `right` are not "almost equal".

    See Also
    --------
    pandas.util.testing.assert_frame_equal

    """
    # pass all kwargs to ensure this function has consistent behavior even if
    # `assert_frame_equal`'s defaults change
    pdt.assert_frame_equal(left, right,
                           check_dtype=True,
                           check_index_type=True,
                           check_column_type=True,
                           check_frame_type=True,
                           check_less_precise=False,
                           check_names=True,
                           by_blocks=False,
                           check_exact=False)
    # this check ensures that empty DataFrames with different indices do not
    # compare equal. exact=True specifies that the type of the indices must be
    # exactly the same
    assert_index_equal(left.index, right.index)


def assert_series_almost_equal(left, right):
    # pass all kwargs to ensure this function has consistent behavior even if
    # `assert_series_equal`'s defaults change
    pdt.assert_series_equal(left, right,
                            check_dtype=True,
                            check_index_type=True,
                            check_series_type=True,
                            check_less_precise=False,
                            check_names=True,
                            check_exact=False,
                            check_datetimelike_compat=False,
                            obj='Series')
    # this check ensures that empty Series with different indices do not
    # compare equal.
    assert_index_equal(left.index, right.index)


def assert_index_equal(a, b):
    pdt.assert_index_equal(a, b,
                           exact=True,
                           check_names=True,
                           check_exact=True)
