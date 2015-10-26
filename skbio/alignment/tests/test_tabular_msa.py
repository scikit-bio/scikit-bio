# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest
import itertools

import six
import numpy as np
import pandas as pd
import numpy.testing as npt

from skbio import Sequence, DNA, RNA, Protein, TabularMSA
from skbio.util import OperationError, UniqueError
from skbio.util._testing import ReallyEqualMixin


class TabularMSASubclass(TabularMSA):
    """Used for testing purposes."""
    pass


class Unorderable(object):
    """For testing unorderable objects in Python 2 and 3."""
    def __lt__(self, other):
        raise TypeError()
    __cmp__ = __lt__


class TestTabularMSA(unittest.TestCase, ReallyEqualMixin):
    def test_from_dict_empty(self):
        self.assertEqual(TabularMSA.from_dict({}), TabularMSA([], keys=[]))

    def test_from_dict_single_sequence(self):
        self.assertEqual(TabularMSA.from_dict({'foo': DNA('ACGT')}),
                         TabularMSA([DNA('ACGT')], keys=['foo']))

    def test_from_dict_multiple_sequences(self):
        msa = TabularMSA.from_dict(
            {1: DNA('ACG'), 2: DNA('GGG'), 3: DNA('TAG')})
        # Sort because order is arbitrary.
        msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACG'), DNA('GGG'), DNA('TAG')], keys=[1, 2, 3]))

    def test_from_dict_invalid_input(self):
        # Basic test to make sure error-checking in the TabularMSA constructor
        # is being invoked.
        with six.assertRaisesRegex(self, ValueError, 'same length'):
            TabularMSA.from_dict({'a': DNA('ACG'), 'b': DNA('ACGT')})

    def test_constructor_invalid_dtype(self):
        with six.assertRaisesRegex(self, TypeError,
                                   'sequence.*alphabet.*Sequence'):
            TabularMSA([Sequence('')])

        with six.assertRaisesRegex(self, TypeError, 'sequence.*alphabet.*int'):
            TabularMSA([42, DNA('')])

    def test_constructor_not_monomorphic(self):
        with six.assertRaisesRegex(self, TypeError, 'mixed types.*RNA.*DNA'):
            TabularMSA([DNA(''), RNA('')])

        with six.assertRaisesRegex(self, TypeError,
                                   'mixed types.*float.*Protein'):
            TabularMSA([Protein(''), Protein(''), 42.0, Protein('')])

    def test_constructor_unequal_length(self):
        with six.assertRaisesRegex(self, ValueError, 'same length.*1 != 0'):
            TabularMSA([Protein(''), Protein('P')])

        with six.assertRaisesRegex(self, ValueError, 'same length.*1 != 3'):
            TabularMSA([Protein('PAW'), Protein('ABC'), Protein('A')])

    def test_constructor_non_iterable(self):
        with self.assertRaises(TypeError):
            TabularMSA(42)

    def test_constructor_non_unique_keys(self):
        with six.assertRaisesRegex(self, UniqueError, 'Duplicate keys:.*42'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], key=lambda x: 42)

        with six.assertRaisesRegex(self, UniqueError, "Duplicate keys:.*'a'"):
            TabularMSA([DNA('', metadata={'id': 'a'}),
                        DNA('', metadata={'id': 'b'}),
                        DNA('', metadata={'id': 'a'})],
                       key='id')

        with six.assertRaisesRegex(self, UniqueError, 'Duplicate keys:.*42'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], keys=iter([42, 42]))

    def test_constructor_non_hashable_keys(self):
        with self.assertRaises(TypeError):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], key=lambda x: [42])

        with self.assertRaises(TypeError):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], keys=iter([[42], [42]]))

    def test_constructor_key_and_keys_both_provided(self):
        with six.assertRaisesRegex(self, ValueError, 'both.*key.*keys'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], key=str, keys=['a', 'b'])

    def test_constructor_keys_length_mismatch(self):
        with six.assertRaisesRegex(self, ValueError,
                                   'Number.*keys.*number.*sequences: 0 != 2'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], keys=iter([]))

    def test_constructor_invalid_metadata(self):
        with self.assertRaises(TypeError):
            TabularMSA([], metadata=42)

    def test_constructor_empty_no_keys(self):
        # sequence empty
        msa = TabularMSA([])
        self.assertIsNone(msa.dtype)
        self.assertEqual(msa.shape, (0, 0))
        with self.assertRaises(OperationError):
            msa.keys
        with self.assertRaises(StopIteration):
            next(iter(msa))

        # position empty
        seqs = [DNA(''), DNA('')]
        msa = TabularMSA(seqs)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (2, 0))
        with self.assertRaises(OperationError):
            msa.keys
        self.assertEqual(list(msa), seqs)

    def test_constructor_empty_with_keys(self):
        # sequence empty
        msa = TabularMSA([], key=lambda x: x)
        npt.assert_array_equal(msa.keys, np.array([]))

        msa = TabularMSA([], keys=iter([]))
        npt.assert_array_equal(msa.keys, np.array([]))

        # position empty
        msa = TabularMSA([DNA('', metadata={'id': 42}),
                          DNA('', metadata={'id': 43})], key='id')
        npt.assert_array_equal(msa.keys, np.array([42, 43]))

        msa = TabularMSA([DNA(''), DNA('')], keys=iter([42, 43]))
        npt.assert_array_equal(msa.keys, np.array([42, 43]))

    def test_constructor_non_empty_no_keys(self):
        # 1x3
        seqs = [DNA('ACG')]
        msa = TabularMSA(seqs)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (1, 3))
        with self.assertRaises(OperationError):
            msa.keys
        self.assertEqual(list(msa), seqs)

        # 3x1
        seqs = [DNA('A'), DNA('C'), DNA('G')]
        msa = TabularMSA(seqs)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (3, 1))
        with self.assertRaises(OperationError):
            msa.keys
        self.assertEqual(list(msa), seqs)

    def test_constructor_non_empty_with_keys(self):
        seqs = [DNA('ACG'), DNA('CGA'), DNA('GTT')]
        msa = TabularMSA(seqs, key=str)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (3, 3))
        npt.assert_array_equal(msa.keys, np.array(['ACG', 'CGA', 'GTT']))
        self.assertEqual(list(msa), seqs)

        msa = TabularMSA(seqs, keys=iter([42, 43, 44]))
        npt.assert_array_equal(msa.keys, np.array([42, 43, 44]))

    def test_constructor_works_with_iterator(self):
        seqs = [DNA('ACG'), DNA('CGA'), DNA('GTT')]
        msa = TabularMSA(iter(seqs), key=str)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (3, 3))
        npt.assert_array_equal(msa.keys, np.array(['ACG', 'CGA', 'GTT']))
        self.assertEqual(list(msa), seqs)

    def test_constructor_handles_missing_metadata_efficiently(self):
        self.assertIsNone(TabularMSA([])._metadata)

    def test_constructor_no_metadata(self):
        self.assertFalse(TabularMSA([]).has_metadata())
        self.assertFalse(
            TabularMSA([DNA('', metadata={'id': 42})]).has_metadata())
        self.assertFalse(
            TabularMSA([DNA('AGC', metadata={'id': 42}),
                        DNA('---', metadata={'id': 43})]).has_metadata())

    def test_constructor_with_metadata(self):
        msa = TabularMSA([], metadata={'foo': 'bar'})
        self.assertEqual(msa.metadata, {'foo': 'bar'})

        msa = TabularMSA([DNA('', metadata={'id': 42})],
                         metadata={'foo': 'bar'})
        self.assertEqual(msa.metadata, {'foo': 'bar'})

        msa = TabularMSA([DNA('AGC'), DNA('---')], metadata={'foo': 'bar'})
        self.assertEqual(msa.metadata, {'foo': 'bar'})

    def test_constructor_makes_shallow_copy_of_metadata(self):
        md = {'foo': 'bar', 42: []}
        msa = TabularMSA([RNA('-.-'), RNA('.-.')], metadata=md)

        self.assertEqual(msa.metadata, md)
        self.assertIsNot(msa.metadata, md)

        md['foo'] = 'baz'
        self.assertEqual(msa.metadata, {'foo': 'bar', 42: []})

        md[42].append(True)
        self.assertEqual(msa.metadata, {'foo': 'bar', 42: [True]})

    def test_dtype(self):
        self.assertIsNone(TabularMSA([]).dtype)
        self.assertIs(TabularMSA([Protein('')]).dtype, Protein)

        with self.assertRaises(AttributeError):
            TabularMSA([]).dtype = DNA

        with self.assertRaises(AttributeError):
            del TabularMSA([]).dtype

    def test_shape(self):
        shape = TabularMSA([DNA('ACG'), DNA('GCA')]).shape
        self.assertEqual(shape, (2, 3))
        self.assertEqual(shape.sequence, shape[0])
        self.assertEqual(shape.position, shape[1])
        with self.assertRaises(TypeError):
            shape[0] = 3

        with self.assertRaises(AttributeError):
            TabularMSA([]).shape = (3, 3)

        with self.assertRaises(AttributeError):
            del TabularMSA([]).shape

    def test_keys_getter(self):
        with six.assertRaisesRegex(self, OperationError,
                                   'Keys do not exist.*reindex'):
            TabularMSA([]).keys

        keys = TabularMSA([DNA('AC'), DNA('AG'), DNA('AT')], key=str).keys
        self.assertIsInstance(keys, np.ndarray)
        self.assertEqual(keys.dtype, object)
        npt.assert_array_equal(keys, np.array(['AC', 'AG', 'AT']))

        # immutable
        with self.assertRaises(ValueError):
            keys[1] = 'AA'
        # original state is maintained
        npt.assert_array_equal(keys, np.array(['AC', 'AG', 'AT']))

    def test_keys_mixed_type(self):
        msa = TabularMSA([DNA('AC'), DNA('CA'), DNA('AA')],
                         keys=['abc', 'd', 42])
        npt.assert_array_equal(msa.keys,
                               np.array(['abc', 'd', 42], dtype=object))

    def test_keys_update_subset_of_keys(self):
        # keys can be copied, modified, then re-set
        msa = TabularMSA([DNA('AC'), DNA('AG'), DNA('AT')], key=str)
        npt.assert_array_equal(msa.keys, np.array(['AC', 'AG', 'AT']))

        new_keys = msa.keys.copy()
        new_keys[1] = 42
        msa.keys = new_keys
        npt.assert_array_equal(msa.keys,
                               np.array(['AC', 42, 'AT'], dtype=object))

        self.assertFalse(msa.keys.flags.writeable)
        self.assertTrue(new_keys.flags.writeable)
        new_keys[1] = 'GG'
        npt.assert_array_equal(msa.keys,
                               np.array(['AC', 42, 'AT'], dtype=object))

    def test_keys_setter_empty(self):
        msa = TabularMSA([])
        self.assertFalse(msa.has_keys())
        msa.keys = iter([])
        npt.assert_array_equal(msa.keys, np.array([]))

    def test_keys_setter_non_empty(self):
        msa = TabularMSA([DNA('AC'), DNA('AG'), DNA('AT')])
        self.assertFalse(msa.has_keys())
        msa.keys = range(3)
        npt.assert_array_equal(msa.keys, np.array([0, 1, 2]))
        msa.keys = range(3, 6)
        npt.assert_array_equal(msa.keys, np.array([3, 4, 5]))

    def test_keys_setter_length_mismatch(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], key=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with six.assertRaisesRegex(self, ValueError,
                                   'Number.*keys.*number.*sequences: 3 != 2'):
            msa.keys = iter(['ab', 'cd', 'ef'])

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

    def test_keys_setter_non_unique_keys(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], key=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with six.assertRaisesRegex(self, UniqueError, 'Duplicate keys:.*42'):
            msa.keys = [42, 42]

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

    def test_keys_setter_non_hashable_keys(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], key=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with self.assertRaises(TypeError):
            msa.keys = [[42], [42]]

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

    def test_keys_deleter_no_keys(self):
        msa = TabularMSA([])
        self.assertFalse(msa.has_keys())
        del msa.keys
        self.assertFalse(msa.has_keys())

    def test_keys_deleter_with_keys(self):
        msa = TabularMSA([RNA('UUU'), RNA('AAA')], key=str)
        self.assertTrue(msa.has_keys())
        del msa.keys
        self.assertFalse(msa.has_keys())

    def test_metadata_getter(self):
        msa = TabularMSA([])
        self.assertIsNone(msa._metadata)
        self.assertEqual(msa.metadata, {})
        self.assertIsNotNone(msa._metadata)
        self.assertIsInstance(msa.metadata, dict)

        msa = TabularMSA([], metadata={42: 'foo', ('hello', 'world'): 43})
        self.assertEqual(msa.metadata, {42: 'foo', ('hello', 'world'): 43})
        self.assertIsInstance(msa.metadata, dict)

        msa.metadata[42] = 'bar'
        self.assertEqual(msa.metadata, {42: 'bar', ('hello', 'world'): 43})

    def test_metadata_setter(self):
        msa = TabularMSA([DNA('A-A'), DNA('A-G')])
        self.assertFalse(msa.has_metadata())

        msa.metadata = {'hello': 'world'}
        self.assertTrue(msa.has_metadata())
        self.assertEqual(msa.metadata, {'hello': 'world'})

        msa.metadata = {}
        self.assertFalse(msa.has_metadata())

    def test_metadata_setter_makes_shallow_copy(self):
        msa = TabularMSA([RNA('-.-'), RNA('.-.')])
        md = {'foo': 'bar', 42: []}
        msa.metadata = md

        self.assertEqual(msa.metadata, md)
        self.assertIsNot(msa.metadata, md)

        md['foo'] = 'baz'
        self.assertEqual(msa.metadata, {'foo': 'bar', 42: []})

        md[42].append(True)
        self.assertEqual(msa.metadata, {'foo': 'bar', 42: [True]})

    def test_metadata_setter_invalid_type(self):
        msa = TabularMSA([Protein('PAW')], metadata={123: 456})

        for md in (None, 0, 'a', ('f', 'o', 'o'), np.array([]),
                   pd.DataFrame()):
            with six.assertRaisesRegex(self, TypeError,
                                       'metadata must be a dict'):
                msa.metadata = md
            self.assertEqual(msa.metadata, {123: 456})

    def test_metadata_deleter(self):
        msa = TabularMSA([Protein('PAW')], metadata={'foo': 'bar'})
        self.assertEqual(msa.metadata, {'foo': 'bar'})

        del msa.metadata
        self.assertIsNone(msa._metadata)
        self.assertFalse(msa.has_metadata())

        # delete again
        del msa.metadata
        self.assertIsNone(msa._metadata)
        self.assertFalse(msa.has_metadata())

        msa = TabularMSA([])
        self.assertIsNone(msa._metadata)
        self.assertFalse(msa.has_metadata())
        del msa.metadata
        self.assertIsNone(msa._metadata)
        self.assertFalse(msa.has_metadata())

    def test_bool(self):
        self.assertFalse(TabularMSA([]))
        self.assertFalse(TabularMSA([RNA('')]))
        self.assertFalse(
            TabularMSA([RNA('', metadata={'id': 1}),
                        RNA('', metadata={'id': 2})], key='id'))

        self.assertTrue(TabularMSA([RNA('U')]))
        self.assertTrue(TabularMSA([RNA('--'), RNA('..')]))
        self.assertTrue(TabularMSA([RNA('AUC'), RNA('GCA')]))

    def test_len(self):
        self.assertEqual(len(TabularMSA([])), 0)
        self.assertEqual(len(TabularMSA([DNA('')])), 1)
        self.assertEqual(len(TabularMSA([DNA('AT'), DNA('AG'), DNA('AT')])), 3)

    def test_iter(self):
        with self.assertRaises(StopIteration):
            next(iter(TabularMSA([])))

        seqs = [DNA(''), DNA('')]
        self.assertEqual(list(iter(TabularMSA(seqs))), seqs)

        seqs = [DNA('AAA'), DNA('GCT')]
        self.assertEqual(list(iter(TabularMSA(seqs))), seqs)

    def test_reversed(self):
        with self.assertRaises(StopIteration):
            next(reversed(TabularMSA([])))

        seqs = [DNA(''), DNA('', metadata={'id': 42})]
        self.assertEqual(list(reversed(TabularMSA(seqs))), seqs[::-1])

        seqs = [DNA('AAA'), DNA('GCT')]
        self.assertEqual(list(reversed(TabularMSA(seqs))), seqs[::-1])

    def test_eq_and_ne(self):
        # Each element contains the components necessary to construct a
        # TabularMSA object: seqs and kwargs. None of these objects (once
        # constructed) should compare equal to one another.
        components = [
            # empties
            ([], {}),
            ([], {'key': str}),
            ([RNA('')], {}),
            ([RNA('')], {'key': str}),

            # 1x1
            ([RNA('U')], {'key': str}),

            # 2x3
            ([RNA('AUG'), RNA('GUA')], {'key': str}),

            ([RNA('AG'), RNA('GG')], {}),
            # has keys
            ([RNA('AG'), RNA('GG')], {'key': str}),
            # different dtype
            ([DNA('AG'), DNA('GG')], {'key': str}),
            # different keys
            ([RNA('AG'), RNA('GG')], {'key': lambda x: str(x) + '42'}),
            # different sequence metadata
            ([RNA('AG', metadata={'id': 42}), RNA('GG')], {'key': str}),
            # different sequence data, same keys
            ([RNA('AG'), RNA('GA')],
             {'key': lambda x: 'AG' if 'AG' in x else 'GG'}),
            # different MSA metadata
            ([RNA('AG'), RNA('GG')], {'metadata': {'foo': 42}}),
            ([RNA('AG'), RNA('GG')], {'metadata': {'foo': 43}}),
            ([RNA('AG'), RNA('GG')], {'metadata': {'foo': 42, 'bar': 43}}),
        ]

        for seqs, kwargs in components:
            obj = TabularMSA(seqs, **kwargs)
            self.assertReallyEqual(obj, obj)
            self.assertReallyEqual(obj, TabularMSA(seqs, **kwargs))
            self.assertReallyEqual(obj, TabularMSASubclass(seqs, **kwargs))

        for (seqs1, kwargs1), (seqs2, kwargs2) in \
                itertools.combinations(components, 2):
            obj1 = TabularMSA(seqs1, **kwargs1)
            obj2 = TabularMSA(seqs2, **kwargs2)
            self.assertReallyNotEqual(obj1, obj2)
            self.assertReallyNotEqual(obj1,
                                      TabularMSASubclass(seqs2, **kwargs2))

        # completely different types
        msa = TabularMSA([])
        self.assertReallyNotEqual(msa, 42)
        self.assertReallyNotEqual(msa, [])
        self.assertReallyNotEqual(msa, {})
        self.assertReallyNotEqual(msa, '')

    def test_eq_constructed_from_different_iterables_compare_equal(self):
        msa1 = TabularMSA([DNA('ACGT')])
        msa2 = TabularMSA((DNA('ACGT'),))
        self.assertReallyEqual(msa1, msa2)

    def test_eq_missing_metadata(self):
        self.assertReallyEqual(TabularMSA([DNA('A')]),
                               TabularMSA([DNA('A')], metadata={}))

    def test_eq_handles_missing_metadata_efficiently(self):
        msa1 = TabularMSA([DNA('ACGT')])
        msa2 = TabularMSA([DNA('ACGT')])
        self.assertReallyEqual(msa1, msa2)

        self.assertIsNone(msa1._metadata)
        self.assertIsNone(msa2._metadata)

    def test_has_metadata(self):
        msa = TabularMSA([])
        self.assertFalse(msa.has_metadata())
        # Handles metadata efficiently.
        self.assertIsNone(msa._metadata)

        self.assertFalse(TabularMSA([], metadata={}).has_metadata())

        self.assertTrue(TabularMSA([], metadata={'': ''}).has_metadata())
        self.assertTrue(TabularMSA([], metadata={'foo': 42}).has_metadata())

    def test_has_keys(self):
        self.assertFalse(TabularMSA([]).has_keys())
        self.assertTrue(TabularMSA([], key=str).has_keys())

        self.assertFalse(TabularMSA([DNA('')]).has_keys())
        self.assertTrue(TabularMSA([DNA('')], key=str).has_keys())

        self.assertFalse(TabularMSA([DNA('ACG'), DNA('GCA')]).has_keys())
        self.assertTrue(
            TabularMSA([DNA('ACG', metadata={'id': 1}),
                        DNA('GCA', metadata={'id': 2})], key='id').has_keys())

        msa = TabularMSA([])
        self.assertFalse(msa.has_keys())
        msa.reindex(key=str)
        self.assertTrue(msa.has_keys())
        msa.reindex()
        self.assertFalse(msa.has_keys())

    def test_reindex_empty(self):
        # sequence empty
        msa = TabularMSA([])
        msa.reindex()
        self.assertEqual(msa, TabularMSA([]))
        self.assertFalse(msa.has_keys())

        msa.reindex(key=str)
        self.assertEqual(msa, TabularMSA([], key=str))
        npt.assert_array_equal(msa.keys, np.array([]))

        msa.reindex(keys=iter([]))
        self.assertEqual(msa, TabularMSA([], keys=iter([])))
        npt.assert_array_equal(msa.keys, np.array([]))

        # position empty
        msa = TabularMSA([DNA('')])
        msa.reindex()
        self.assertEqual(msa, TabularMSA([DNA('')]))
        self.assertFalse(msa.has_keys())

        msa.reindex(key=str)
        self.assertEqual(msa, TabularMSA([DNA('')], key=str))
        npt.assert_array_equal(msa.keys, np.array(['']))

        msa.reindex(keys=iter(['a']))
        self.assertEqual(msa, TabularMSA([DNA('')], keys=iter(['a'])))
        npt.assert_array_equal(msa.keys, np.array(['a']))

    def test_reindex_non_empty(self):
        msa = TabularMSA([DNA('ACG', metadata={'id': 1}),
                          DNA('AAA', metadata={'id': 2})], key=str)
        npt.assert_array_equal(msa.keys, np.array(['ACG', 'AAA']))

        msa.reindex(key='id')
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACG', metadata={'id': 1}),
                        DNA('AAA', metadata={'id': 2})], key='id'))
        npt.assert_array_equal(msa.keys, np.array([1, 2]))

        msa.reindex(keys=iter('ab'))
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACG', metadata={'id': 1}),
                        DNA('AAA', metadata={'id': 2})], keys=iter('ab')))
        npt.assert_array_equal(msa.keys, np.array(['a', 'b']))

        msa.reindex()
        self.assertFalse(msa.has_keys())

    def test_reindex_makes_copy_of_keys(self):
        msa = TabularMSA([DNA('AC'), DNA('AG'), DNA('AT')])
        keys = np.asarray([1, 2, 3])
        msa.reindex(keys=keys)
        npt.assert_array_equal(msa.keys, np.array([1, 2, 3]))

        self.assertFalse(msa.keys.flags.writeable)
        self.assertTrue(keys.flags.writeable)
        keys[1] = 42
        npt.assert_array_equal(msa.keys, np.array([1, 2, 3]))

    def test_reindex_key_and_keys_both_provided(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], key=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with six.assertRaisesRegex(self, ValueError, 'both.*key.*keys'):
            msa.reindex(key=str, keys=['a', 'b'])

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

    def test_reindex_keys_length_mismatch(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], key=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with six.assertRaisesRegex(self, ValueError,
                                   'Number.*keys.*number.*sequences: 0 != 2'):
            msa.reindex(keys=iter([]))

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

    def test_reindex_non_unique_keys(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], key=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with six.assertRaisesRegex(self, UniqueError, 'Duplicate keys:.*42'):
            msa.reindex(key=lambda x: 42)

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

        with six.assertRaisesRegex(self, UniqueError, 'Duplicate keys:.*42'):
            msa.reindex(keys=[42, 42])

        npt.assert_array_equal(msa.keys, keys)

    def test_reindex_non_hashable_keys(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], key=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with self.assertRaises(TypeError):
            msa.reindex(key=lambda x: [42])

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

        with self.assertRaises(TypeError):
            msa.reindex(keys=[[42], [42]])

        npt.assert_array_equal(msa.keys, keys)

    def test_sort_no_msa_keys_and_key_not_specified(self):
        msa = TabularMSA([])
        with self.assertRaises(OperationError):
            msa.sort()
        # original state is maintained
        self.assertEqual(msa, TabularMSA([]))

        msa = TabularMSA([DNA('TC'), DNA('AA')])
        with self.assertRaises(OperationError):
            msa.sort()
        self.assertEqual(msa, TabularMSA([DNA('TC'), DNA('AA')]))

    def test_sort_on_unorderable_msa_keys(self):
        unorderable = Unorderable()
        msa = TabularMSA([DNA('AAA'), DNA('ACG')], keys=[42, unorderable])
        with self.assertRaises(TypeError):
            msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([DNA('AAA'), DNA('ACG')], keys=[42, unorderable]))

    def test_sort_on_unorderable_key(self):
        unorderable = Unorderable()
        msa = TabularMSA([
            DNA('AAA', metadata={'id': 42}),
            DNA('ACG', metadata={'id': unorderable})], keys=[42, 43])
        with self.assertRaises(TypeError):
            msa.sort(key='id')
        self.assertEqual(
            msa,
            TabularMSA([
                DNA('AAA', metadata={'id': 42}),
                DNA('ACG', metadata={'id': unorderable})], keys=[42, 43]))

    def test_sort_on_invalid_key(self):
        msa = TabularMSA([DNA('AAA'), DNA('ACG')], keys=[42, 43])
        with self.assertRaises(KeyError):
            msa.sort(key='id')
        self.assertEqual(
            msa,
            TabularMSA([DNA('AAA'), DNA('ACG')], keys=[42, 43]))

    def test_sort_empty_on_msa_keys(self):
        msa = TabularMSA([], keys=[])
        msa.sort()
        self.assertEqual(msa, TabularMSA([], keys=[]))

        msa = TabularMSA([], keys=[])
        msa.sort(reverse=True)
        self.assertEqual(msa, TabularMSA([], keys=[]))

    def test_sort_single_sequence_on_msa_keys(self):
        msa = TabularMSA([DNA('ACGT')], keys=[42])
        msa.sort()
        self.assertEqual(msa, TabularMSA([DNA('ACGT')], keys=[42]))

        msa = TabularMSA([DNA('ACGT')], keys=[42])
        msa.sort(reverse=True)
        self.assertEqual(msa, TabularMSA([DNA('ACGT')], keys=[42]))

    def test_sort_multiple_sequences_on_msa_keys(self):
        msa = TabularMSA([
            DNA('TC'), DNA('GG'), DNA('CC')], keys=['z', 'a', 'b'])
        msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([
                DNA('GG'), DNA('CC'), DNA('TC')], keys=['a', 'b', 'z']))

        msa = TabularMSA([
            DNA('TC'), DNA('GG'), DNA('CC')], keys=['z', 'a', 'b'])
        msa.sort(reverse=True)
        self.assertEqual(
            msa,
            TabularMSA([
                DNA('TC'), DNA('CC'), DNA('GG')], keys=['z', 'b', 'a']))

    def test_sort_empty_no_msa_keys_on_metadata_key(self):
        msa = TabularMSA([])
        msa.sort(key='id')
        self.assertEqual(msa, TabularMSA([]))

        msa = TabularMSA([])
        msa.sort(key='id', reverse=True)
        self.assertEqual(msa, TabularMSA([]))

    def test_sort_empty_no_msa_keys_on_callable_key(self):
        msa = TabularMSA([])
        msa.sort(key=str)
        self.assertEqual(msa, TabularMSA([]))

        msa = TabularMSA([])
        msa.sort(key=str, reverse=True)
        self.assertEqual(msa, TabularMSA([]))

    def test_sort_empty_with_msa_keys_on_metadata_key(self):
        msa = TabularMSA([], keys=[])
        msa.sort(key='id')
        self.assertEqual(msa, TabularMSA([], keys=[]))

        msa = TabularMSA([], keys=[])
        msa.sort(key='id', reverse=True)
        self.assertEqual(msa, TabularMSA([], keys=[]))

    def test_sort_empty_with_msa_keys_on_callable_key(self):
        msa = TabularMSA([], keys=[])
        msa.sort(key=str)
        self.assertEqual(msa, TabularMSA([], keys=[]))

        msa = TabularMSA([], keys=[])
        msa.sort(key=str, reverse=True)
        self.assertEqual(msa, TabularMSA([], keys=[]))

    def test_sort_single_sequence_no_msa_keys_on_metadata_key(self):
        msa = TabularMSA([RNA('UCA', metadata={'id': 42})])
        msa.sort(key='id')
        self.assertEqual(msa, TabularMSA([RNA('UCA', metadata={'id': 42})]))

        msa = TabularMSA([RNA('UCA', metadata={'id': 42})])
        msa.sort(key='id', reverse=True)
        self.assertEqual(msa, TabularMSA([RNA('UCA', metadata={'id': 42})]))

    def test_sort_single_sequence_no_msa_keys_on_callable_key(self):
        msa = TabularMSA([RNA('UCA')])
        msa.sort(key=str)
        self.assertEqual(msa, TabularMSA([RNA('UCA')]))

        msa = TabularMSA([RNA('UCA')])
        msa.sort(key=str, reverse=True)
        self.assertEqual(msa, TabularMSA([RNA('UCA')]))

    def test_sort_single_sequence_with_msa_keys_on_metadata_key(self):
        msa = TabularMSA([RNA('UCA', metadata={'id': 42})], keys=['foo'])
        msa.sort(key='id')
        self.assertEqual(
            msa, TabularMSA([RNA('UCA', metadata={'id': 42})], keys=['foo']))

        msa = TabularMSA([RNA('UCA', metadata={'id': 42})], keys=['foo'])
        msa.sort(key='id', reverse=True)
        self.assertEqual(
            msa, TabularMSA([RNA('UCA', metadata={'id': 42})], keys=['foo']))

    def test_sort_single_sequence_with_msa_keys_on_callable_key(self):
        msa = TabularMSA([RNA('UCA')], keys=['foo'])
        msa.sort(key=str)
        self.assertEqual(msa, TabularMSA([RNA('UCA')], keys=['foo']))

        msa = TabularMSA([RNA('UCA')], keys=['foo'])
        msa.sort(key=str, reverse=True)
        self.assertEqual(msa, TabularMSA([RNA('UCA')], keys=['foo']))

    def test_sort_multiple_sequences_no_msa_keys_on_metadata_key(self):
        msa = TabularMSA([RNA('UCA', metadata={'id': 41}),
                          RNA('AAA', metadata={'id': 44}),
                          RNA('GAC', metadata={'id': -1}),
                          RNA('GAC', metadata={'id': 42})])
        msa.sort(key='id')
        self.assertEqual(msa, TabularMSA([RNA('GAC', metadata={'id': -1}),
                                          RNA('UCA', metadata={'id': 41}),
                                          RNA('GAC', metadata={'id': 42}),
                                          RNA('AAA', metadata={'id': 44})]))

        msa = TabularMSA([RNA('UCA', metadata={'id': 41}),
                          RNA('AAA', metadata={'id': 44}),
                          RNA('GAC', metadata={'id': -1}),
                          RNA('GAC', metadata={'id': 42})])
        msa.sort(key='id', reverse=True)
        self.assertEqual(msa, TabularMSA([RNA('AAA', metadata={'id': 44}),
                                          RNA('GAC', metadata={'id': 42}),
                                          RNA('UCA', metadata={'id': 41}),
                                          RNA('GAC', metadata={'id': -1})]))

    def test_sort_multiple_sequences_no_msa_keys_on_callable_key(self):
        msa = TabularMSA([RNA('UCC'),
                          RNA('UCG'),
                          RNA('UCA')])
        msa.sort(key=str)
        self.assertEqual(msa, TabularMSA([RNA('UCA'), RNA('UCC'), RNA('UCG')]))

        msa = TabularMSA([RNA('UCC'),
                          RNA('UCG'),
                          RNA('UCA')])
        msa.sort(key=str, reverse=True)
        self.assertEqual(msa, TabularMSA([RNA('UCG'), RNA('UCC'), RNA('UCA')]))

    def test_sort_multiple_sequences_with_msa_keys_on_metadata_key(self):
        msa = TabularMSA([DNA('TCA', metadata={'#': 41.2}),
                          DNA('AAA', metadata={'#': 44.5}),
                          DNA('GAC', metadata={'#': 42.999})],
                         keys=[None, ('hello', 'world'), True])
        msa.sort(key='#')
        self.assertEqual(
            msa,
            TabularMSA([DNA('TCA', metadata={'#': 41.2}),
                        DNA('GAC', metadata={'#': 42.999}),
                        DNA('AAA', metadata={'#': 44.5})],
                       keys=[None, True, ('hello', 'world')]))

        msa = TabularMSA([DNA('TCA', metadata={'#': 41.2}),
                          DNA('AAA', metadata={'#': 44.5}),
                          DNA('GAC', metadata={'#': 42.999})],
                         keys=[None, ('hello', 'world'), True])
        msa.sort(key='#', reverse=True)
        self.assertEqual(
            msa,
            TabularMSA([DNA('AAA', metadata={'#': 44.5}),
                        DNA('GAC', metadata={'#': 42.999}),
                        DNA('TCA', metadata={'#': 41.2})],
                       keys=[('hello', 'world'), True, None]))

    def test_sort_multiple_sequences_with_msa_keys_on_callable_key(self):
        msa = TabularMSA([RNA('UCC'),
                          RNA('UCG'),
                          RNA('UCA')], keys=[1, 'abc', None])
        msa.sort(key=str)
        self.assertEqual(msa, TabularMSA([RNA('UCA'), RNA('UCC'), RNA('UCG')],
                                         keys=[None, 1, 'abc']))

        msa = TabularMSA([RNA('UCC'),
                          RNA('UCG'),
                          RNA('UCA')], keys=[1, 'abc', None])
        msa.sort(key=str, reverse=True)
        self.assertEqual(msa, TabularMSA([RNA('UCG'), RNA('UCC'), RNA('UCA')],
                                         keys=['abc', 1, None]))

    def test_sort_on_key_with_some_repeats(self):
        msa = TabularMSA([
            DNA('TCCG', metadata={'id': 10}),
            DNA('TAGG', metadata={'id': 10}),
            DNA('GGGG', metadata={'id': 8}),
            DNA('ACGT', metadata={'id': 0}),
            DNA('TAGG', metadata={'id': 10})], keys=range(5))
        msa.sort(key='id')
        self.assertEqual(
            msa,
            TabularMSA([
                DNA('ACGT', metadata={'id': 0}),
                DNA('GGGG', metadata={'id': 8}),
                DNA('TCCG', metadata={'id': 10}),
                DNA('TAGG', metadata={'id': 10}),
                DNA('TAGG', metadata={'id': 10})], keys=[3, 2, 0, 1, 4]))

    def test_sort_on_key_with_all_repeats(self):
        msa = TabularMSA([
            DNA('TTT', metadata={'id': 'a'}),
            DNA('TTT', metadata={'id': 'b'}),
            DNA('TTT', metadata={'id': 'c'})], keys=range(3))
        msa.sort(key=str)
        self.assertEqual(
            msa,
            TabularMSA([
                DNA('TTT', metadata={'id': 'a'}),
                DNA('TTT', metadata={'id': 'b'}),
                DNA('TTT', metadata={'id': 'c'})], keys=range(3)))

    def test_sort_mixed_key_types(self):
        msa = TabularMSA([
            DNA('GCG', metadata={'id': 41}),
            DNA('CGC', metadata={'id': 42.2}),
            DNA('TCT', metadata={'id': 42})])
        msa.sort(key='id')
        self.assertEqual(
            msa,
            TabularMSA([
                DNA('GCG', metadata={'id': 41}),
                DNA('TCT', metadata={'id': 42}),
                DNA('CGC', metadata={'id': 42.2})]))

        msa = TabularMSA([
            DNA('GCG'),
            DNA('CGC'),
            DNA('TCT')], keys=[41, 42.2, 42])
        msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([
                DNA('GCG'),
                DNA('TCT'),
                DNA('CGC')], keys=[41, 42, 42.2]))

    def test_sort_already_sorted(self):
        msa = TabularMSA([DNA('TC'), DNA('GG'), DNA('CC')], keys=[1, 2, 3])
        msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([DNA('TC'), DNA('GG'), DNA('CC')], keys=[1, 2, 3]))

        msa = TabularMSA([DNA('TC'), DNA('GG'), DNA('CC')], keys=[3, 2, 1])
        msa.sort(reverse=True)
        self.assertEqual(
            msa,
            TabularMSA([DNA('TC'), DNA('GG'), DNA('CC')], keys=[3, 2, 1]))

    def test_sort_reverse_sorted(self):
        msa = TabularMSA([DNA('T'), DNA('G'), DNA('A')], keys=[3, 2, 1])
        msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([DNA('A'), DNA('G'), DNA('T')], keys=[1, 2, 3]))

        msa = TabularMSA([DNA('T'), DNA('G'), DNA('A')], keys=[1, 2, 3])
        msa.sort(reverse=True)
        self.assertEqual(
            msa,
            TabularMSA([DNA('A'), DNA('G'), DNA('T')], keys=[3, 2, 1]))

    def test_sort_identical_sequences(self):
        msa = TabularMSA([DNA(''), DNA(''), DNA('')], keys=['ab', 'aa', 'ac'])
        msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([DNA(''), DNA(''), DNA('')], keys=['aa', 'ab', 'ac']))

    def test_to_dict_no_keys(self):
        with self.assertRaises(OperationError):
            TabularMSA([]).to_dict()

        with self.assertRaises(OperationError):
            TabularMSA([DNA('AGCT'), DNA('TCGA')]).to_dict()

    def test_to_dict_empty(self):
        self.assertEqual(TabularMSA([], keys=[]).to_dict(), {})
        self.assertEqual(TabularMSA([RNA('')], keys=['foo']).to_dict(),
                         {'foo': RNA('')})

    def test_to_dict_non_empty(self):
        seqs = [Protein('PAW', metadata={'id': 42}),
                Protein('WAP', metadata={'id': -999})]
        msa = TabularMSA(seqs, key='id')
        self.assertEqual(msa.to_dict(), {42: seqs[0], -999: seqs[1]})

    def test_from_dict_to_dict_roundtrip(self):
        d = {}
        self.assertEqual(TabularMSA.from_dict(d).to_dict(), d)

        # can roundtrip even with mixed key types
        d1 = {'a': DNA('CAT'), 42: DNA('TAG')}
        d2 = TabularMSA.from_dict(d1).to_dict()
        self.assertEqual(d2, d1)
        self.assertIs(d1['a'], d2['a'])
        self.assertIs(d1[42], d2[42])


if __name__ == "__main__":
    unittest.main()
