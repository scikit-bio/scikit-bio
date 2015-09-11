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


class TestTabularMSA(unittest.TestCase, ReallyEqualMixin):
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
        npt.assert_array_equal(keys, np.array(['AC', 'AG', 'AT']))

        # immutable
        with self.assertRaises(ValueError):
            keys[1] = 'AA'
        # original state is maintained
        npt.assert_array_equal(keys, np.array(['AC', 'AG', 'AT']))

    def test_keys_update_subset_of_keys(self):
        # keys can be copied, modified, then re-set
        msa = TabularMSA([DNA('AC'), DNA('AG'), DNA('AT')], key=str)
        npt.assert_array_equal(msa.keys, np.array(['AC', 'AG', 'AT']))

        new_keys = msa.keys.copy()
        new_keys[1] = 'AA'
        msa.keys = new_keys
        npt.assert_array_equal(msa.keys, np.array(['AC', 'AA', 'AT']))

        self.assertFalse(msa.keys.flags.writeable)
        self.assertTrue(new_keys.flags.writeable)
        new_keys[1] = 'GG'
        npt.assert_array_equal(msa.keys, np.array(['AC', 'AA', 'AT']))

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

    def test_get_cached_key_str(self):
        msa = TabularMSA([DNA('', metadata={'id': 42}),
                          DNA('', metadata={'id': 43})], key='id')
        key = msa.get_cached_key()
        self.assertEqual(key, 'id')

    def test_get_cached_key_callable(self):
        def key_func(x):
            return x.metadata['id']

        msa = TabularMSA([DNA('', metadata={'id': 42}),
                          DNA('', metadata={'id': 43})], key=key_func)
        key = msa.get_cached_key()
        self.assertEqual(key, key_func)

    def test_get_cached_key_no_key_exists(self):
        msa = TabularMSA([DNA(''), DNA('')])
        with six.assertRaisesRegex(
                self, OperationError,
                "MSA requires a key but none was provided, and no "
                "cached key exists"):
            msa.get_cached_key()

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


class TestAppend(unittest.TestCase):
    def setUp(self):
        self.msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

    def test_simple(self):
        self.msa.append(DNA('AAAA'))
        self.assertEqual(self.msa.shape, (3, 4))
        # TODO: Hack to get last seq. once __getitem__ is implemented use
        # msa[-1]
        seq = None
        for seq in self.msa:
            pass
        self.assertEqual(seq, DNA('AAAA'))

    def test_to_empty_msa(self):
        msa = TabularMSA([])
        msa.append(DNA('ACGT'))
        self.assertEqual(len(msa), 1)

    def test_to_empty_msa_invalid_dtype(self):
        msa = TabularMSA([])
        with six.assertRaisesRegex(self, TypeError,
                                   'sequence.*alphabet.*Sequence'):
            msa.append(Sequence(''))

    def test_wrong_dtype_rna(self):
        with six.assertRaisesRegex(self, TypeError, 'mixed types.*RNA.*DNA'):
            self.msa.append(RNA('UUUU'))

    def test_wrong_dtype_float(self):
        with six.assertRaisesRegex(self, TypeError,
                                   'mixed types.*float.*DNA'):
            self.msa.append(42.0)

    def test_wrong_length(self):
        with six.assertRaisesRegex(self, ValueError, 'same length.*5 != 4'):
            self.msa.append(DNA('ACGTA'))

    def test_with_key(self):
        msa = TabularMSA([DNA('', metadata={'id': 'a'}),
                          DNA('', metadata={'id': 'b'})],
                         key='id')
        msa.append(DNA('', metadata={'id': 'c'}), key='id')
        npt.assert_array_equal(msa.keys, np.array(['a', 'b', 'c']))

    def test_with_key_msa_has_no_keys(self):
        with six.assertRaisesRegex(self, OperationError,
                                   "key was provided but MSA does not have "
                                   "keys"):
            self.msa.append(DNA('AAAA'), 'id')

    def test_no_key_msa_has_keys(self):
        msa = TabularMSA([DNA('', metadata={'id': 'a'}),
                          DNA('', metadata={'id': 'b'})],
                         key='id')
        msa.append(DNA('', metadata={'id': 'c'}))
        npt.assert_array_equal(msa.keys, np.array(['a', 'b', 'c']))

    def test_no_key_msa_has_keys_but_not_cached(self):
        msa = TabularMSA([DNA(''), DNA('')], keys=['a', 'b'])
        with six.assertRaisesRegex(self, OperationError,
                                   "MSA requires a key but none was "
                                   "provided, and no cached key exists"):
            msa.append(DNA(''))


if __name__ == "__main__":
    unittest.main()
