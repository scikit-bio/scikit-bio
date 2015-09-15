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
        with six.assertRaisesRegex(self, ValueError, 'must match the length'):
            TabularMSA.from_dict({'a': DNA('ACG'), 'b': DNA('ACGT')})

    def test_constructor_invalid_dtype(self):
        with six.assertRaisesRegex(self, TypeError,
                                   'sequence.*alphabet.*Sequence'):
            TabularMSA([Sequence('')])

        with six.assertRaisesRegex(self, TypeError, 'sequence.*alphabet.*int'):
            TabularMSA([42, DNA('')])

    def test_constructor_not_monomorphic(self):
        with six.assertRaisesRegex(self, TypeError,
                                   'Must match the type.*RNA.*DNA'):
            TabularMSA([DNA(''), RNA('')])

        with six.assertRaisesRegex(self, TypeError,
                                   'Must match the type.*float.*Protein'):
            TabularMSA([Protein(''), Protein(''), 42.0, Protein('')])

    def test_constructor_unequal_length(self):
        with six.assertRaisesRegex(self, ValueError,
                                   'must match the length.*1 != 0'):
            TabularMSA([Protein(''), Protein('P')])

        with six.assertRaisesRegex(self, ValueError,
                                   'must match the length.*1 != 3'):
            TabularMSA([Protein('PAW'), Protein('ABC'), Protein('A')])

    def test_constructor_non_iterable(self):
        with self.assertRaises(TypeError):
            TabularMSA(42)

    def test_constructor_non_unique_keys(self):
        with six.assertRaisesRegex(self, UniqueError, 'Duplicate keys:.*42'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=lambda x: 42)

        with six.assertRaisesRegex(self, UniqueError, "Duplicate keys:.*'a'"):
            TabularMSA([DNA('', metadata={'id': 'a'}),
                        DNA('', metadata={'id': 'b'}),
                        DNA('', metadata={'id': 'a'})],
                       minter='id')

        with six.assertRaisesRegex(self, UniqueError, 'Duplicate keys:.*42'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], keys=iter([42, 42]))

    def test_constructor_non_hashable_keys(self):
        with self.assertRaises(TypeError):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=lambda x: [42])

        with self.assertRaises(TypeError):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], keys=iter([[42], [42]]))

    def test_constructor_minter_and_keys_both_provided(self):
        with six.assertRaisesRegex(self, ValueError, 'both.*minter.*keys'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str, keys=['a', 'b'])

    def test_constructor_keys_length_mismatch(self):
        with six.assertRaisesRegex(self, ValueError,
                                   'Number.*keys.*number.*sequences: 0 != 2'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], keys=iter([]))

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
        msa = TabularMSA([], minter=lambda x: x)
        npt.assert_array_equal(msa.keys, np.array([]))

        msa = TabularMSA([], keys=iter([]))
        npt.assert_array_equal(msa.keys, np.array([]))

        # position empty
        msa = TabularMSA([DNA('', metadata={'id': 42}),
                          DNA('', metadata={'id': 43})], minter='id')
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
        msa = TabularMSA(seqs, minter=str)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (3, 3))
        npt.assert_array_equal(msa.keys, np.array(['ACG', 'CGA', 'GTT']))
        self.assertEqual(list(msa), seqs)

        msa = TabularMSA(seqs, keys=iter([42, 43, 44]))
        npt.assert_array_equal(msa.keys, np.array([42, 43, 44]))

    def test_constructor_works_with_iterator(self):
        seqs = [DNA('ACG'), DNA('CGA'), DNA('GTT')]
        msa = TabularMSA(iter(seqs), minter=str)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (3, 3))
        npt.assert_array_equal(msa.keys, np.array(['ACG', 'CGA', 'GTT']))
        self.assertEqual(list(msa), seqs)

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

        keys = TabularMSA([DNA('AC'), DNA('AG'), DNA('AT')], minter=str).keys
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
        msa = TabularMSA([DNA('AC'), DNA('AG'), DNA('AT')], minter=str)
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
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with six.assertRaisesRegex(self, ValueError,
                                   'Number.*keys.*number.*sequences: 3 != 2'):
            msa.keys = iter(['ab', 'cd', 'ef'])

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

    def test_keys_setter_non_unique_keys(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with six.assertRaisesRegex(self, UniqueError, 'Duplicate keys:.*42'):
            msa.keys = [42, 42]

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

    def test_keys_setter_non_hashable_keys(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str)
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
        msa = TabularMSA([RNA('UUU'), RNA('AAA')], minter=str)
        self.assertTrue(msa.has_keys())
        del msa.keys
        self.assertFalse(msa.has_keys())

    def test_get_cached_minter_str(self):
        msa = TabularMSA([DNA('', metadata={'id': 42}),
                          DNA('', metadata={'id': 43})], minter='id')
        key = msa.minter
        self.assertEqual(key, 'id')

    def test_get_cached_minter_callable(self):
        def minter_func(x):
            return x.metadata['id']

        msa = TabularMSA([DNA('', metadata={'id': 42}),
                          DNA('', metadata={'id': 43})], minter=minter_func)
        key = msa.minter
        self.assertEqual(key, minter_func)

    def test_get_cached_minter_no_minter_exists(self):
        msa = TabularMSA([DNA(''), DNA('')])
        with six.assertRaisesRegex(
                self, OperationError,
                "MSA requires a key minter but none was provided, and no "
                "cached minter exists"):
            msa.minter

    def test_bool(self):
        self.assertFalse(TabularMSA([]))
        self.assertFalse(TabularMSA([RNA('')]))
        self.assertFalse(
            TabularMSA([RNA('', metadata={'id': 1}),
                        RNA('', metadata={'id': 2})], minter='id'))

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
            ([], {'minter': str}),
            ([RNA('')], {}),
            ([RNA('')], {'minter': str}),

            # 1x1
            ([RNA('U')], {'minter': str}),

            # 2x3
            ([RNA('AUG'), RNA('GUA')], {'minter': str}),

            ([RNA('AG'), RNA('GG')], {}),
            # has keys
            ([RNA('AG'), RNA('GG')], {'minter': str}),
            # different dtype
            ([DNA('AG'), DNA('GG')], {'minter': str}),
            # different keys
            ([RNA('AG'), RNA('GG')], {'minter': lambda x: str(x) + '42'}),
            # different sequence metadata
            ([RNA('AG', metadata={'id': 42}), RNA('GG')], {'minter': str}),
            # different sequence data, same keys
            ([RNA('AG'), RNA('GA')],
             {'minter': lambda x: 'AG' if 'AG' in x else 'GG'}),
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

    def test_has_keys(self):
        self.assertFalse(TabularMSA([]).has_keys())
        self.assertTrue(TabularMSA([], minter=str).has_keys())

        self.assertFalse(TabularMSA([DNA('')]).has_keys())
        self.assertTrue(TabularMSA([DNA('')], minter=str).has_keys())

        self.assertFalse(TabularMSA([DNA('ACG'), DNA('GCA')]).has_keys())
        self.assertTrue(
            TabularMSA([DNA('ACG', metadata={'id': 1}),
                        DNA('GCA', metadata={'id': 2})],
                       minter='id').has_keys())

        msa = TabularMSA([])
        self.assertFalse(msa.has_keys())
        msa.reindex(minter=str)
        self.assertTrue(msa.has_keys())
        msa.reindex()
        self.assertFalse(msa.has_keys())

    def test_reindex_empty(self):
        # sequence empty
        msa = TabularMSA([])
        msa.reindex()
        self.assertEqual(msa, TabularMSA([]))
        self.assertFalse(msa.has_keys())

        msa.reindex(minter=str)
        self.assertEqual(msa, TabularMSA([], minter=str))
        npt.assert_array_equal(msa.keys, np.array([]))

        msa.reindex(keys=iter([]))
        self.assertEqual(msa, TabularMSA([], keys=iter([])))
        npt.assert_array_equal(msa.keys, np.array([]))

        # position empty
        msa = TabularMSA([DNA('')])
        msa.reindex()
        self.assertEqual(msa, TabularMSA([DNA('')]))
        self.assertFalse(msa.has_keys())

        msa.reindex(minter=str)
        self.assertEqual(msa, TabularMSA([DNA('')], minter=str))
        npt.assert_array_equal(msa.keys, np.array(['']))

        msa.reindex(keys=iter(['a']))
        self.assertEqual(msa, TabularMSA([DNA('')], keys=iter(['a'])))
        npt.assert_array_equal(msa.keys, np.array(['a']))

    def test_reindex_non_empty(self):
        msa = TabularMSA([DNA('ACG', metadata={'id': 1}),
                          DNA('AAA', metadata={'id': 2})], minter=str)
        npt.assert_array_equal(msa.keys, np.array(['ACG', 'AAA']))

        msa.reindex(minter='id')
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACG', metadata={'id': 1}),
                        DNA('AAA', metadata={'id': 2})], minter='id'))
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

    def test_reindex_minter_and_keys_both_provided(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with six.assertRaisesRegex(self, ValueError, 'both.*minter.*keys'):
            msa.reindex(minter=str, keys=['a', 'b'])

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

    def test_reindex_keys_length_mismatch(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with six.assertRaisesRegex(self, ValueError,
                                   'Number.*keys.*number.*sequences: 0 != 2'):
            msa.reindex(keys=iter([]))

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

    def test_reindex_non_unique_keys(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with six.assertRaisesRegex(self, UniqueError, 'Duplicate keys:.*42'):
            msa.reindex(minter=lambda x: 42)

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

        with six.assertRaisesRegex(self, UniqueError, 'Duplicate keys:.*42'):
            msa.reindex(keys=[42, 42])

        npt.assert_array_equal(msa.keys, keys)

    def test_reindex_non_hashable_keys(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str)
        keys = np.array(['ACGT', 'TGCA'])
        npt.assert_array_equal(msa.keys, keys)

        with self.assertRaises(TypeError):
            msa.reindex(minter=lambda x: [42])

        # original state is maintained
        npt.assert_array_equal(msa.keys, keys)

        with self.assertRaises(TypeError):
            msa.reindex(keys=[[42], [42]])

        npt.assert_array_equal(msa.keys, keys)

    def test_sort_no_msa_keys_and_minter_not_specified(self):
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
        msa = TabularMSA(seqs, minter='id')
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


class TestAppend(unittest.TestCase):
    def setUp(self):
        self.msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

    def test_simple(self):
        expected = TabularMSA([DNA('ACGT'), DNA('TGCA'), DNA('AAAA')])
        self.msa.append(DNA('AAAA'))
        self.assertEqual(self.msa, expected)

    def test_to_empty_msa(self):
        expected = TabularMSA([DNA('ACGT')])
        msa = TabularMSA([])
        msa.append(DNA('ACGT'))
        self.assertEqual(msa, expected)

    def test_to_empty_msa_invalid_dtype(self):
        msa = TabularMSA([])
        with six.assertRaisesRegex(self, TypeError,
                                   'sequence.*alphabet.*Sequence'):
            msa.append(Sequence(''))

    def test_wrong_dtype_rna(self):
        with six.assertRaisesRegex(self, TypeError,
                                   'Must match the type.*RNA.*DNA'):
            self.msa.append(RNA('UUUU'))

    def test_wrong_dtype_float(self):
        with six.assertRaisesRegex(self, TypeError,
                                   'Must match the type.*float.*DNA'):
            self.msa.append(42.0)

    def test_wrong_length(self):
        with six.assertRaisesRegex(self, ValueError,
                                   'must match the length.*5 != 4'):
            self.msa.append(DNA('ACGTA'))

    def test_with_key(self):
        to_append = DNA('', metadata={'id': 'c'})
        expected = TabularMSA(
            [DNA('', metadata={'id': 'a'}),
             DNA('', metadata={'id': 'b'}),
             to_append],
            minter='id')
        msa = TabularMSA([DNA('', metadata={'id': 'a'}),
                          DNA('', metadata={'id': 'b'})],
                         minter='id')
        msa.append(to_append, minter='id')
        self.assertEqual(msa, expected)

    def test_with_minter_msa_has_no_keys(self):
        with six.assertRaisesRegex(self, OperationError,
                                   "key was provided but MSA does not have "
                                   "keys"):
            self.msa.append(DNA('AAAA'), 'id')

    def test_no_minter_msa_has_keys(self):
        to_append = DNA('', metadata={'id': 'c'})
        expected = TabularMSA(
            [DNA('', metadata={'id': 'a'}),
             DNA('', metadata={'id': 'b'}),
             to_append],
            minter='id')
        msa = TabularMSA([DNA('', metadata={'id': 'a'}),
                          DNA('', metadata={'id': 'b'})],
                         minter='id')
        msa.append(to_append)
        self.assertEqual(msa, expected)

    def test_no_minter_msa_has_keys_but_not_cached(self):
        msa = TabularMSA([DNA(''), DNA('')], keys=['a', 'b'])
        with six.assertRaisesRegex(self, OperationError,
                                   "MSA requires a key minter but none was "
                                   "provided, and no cached minter exists"):
            msa.append(DNA(''))


if __name__ == "__main__":
    unittest.main()
