# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import collections.abc
import copy
import unittest
import functools
import itertools
import types

import numpy as np
import numpy.testing as npt
import pandas as pd
import scipy.stats

from skbio import Sequence, DNA, RNA, Protein, TabularMSA
from skbio.alignment import AlignPath
from skbio.sequence import GrammaredSequence
from skbio.util import classproperty
from skbio.util._decorator import overrides
from skbio.util._testing import ReallyEqualMixin
from skbio.metadata._testing import (MetadataMixinTests,
                                     PositionalMetadataMixinTests)
from skbio.util import assert_data_frame_almost_equal
from skbio.util._testing import assert_index_equal


class TabularMSASubclass(TabularMSA):
    """Used for testing purposes."""
    pass


class TestTabularMSAMetadata(unittest.TestCase, ReallyEqualMixin,
                             MetadataMixinTests):
    def setUp(self):
        self._metadata_constructor_ = functools.partial(TabularMSA, [])


class TestTabularMSAPositionalMetadata(unittest.TestCase, ReallyEqualMixin,
                                       PositionalMetadataMixinTests):
    def setUp(self):
        def factory(axis_len, positional_metadata=None):
            return TabularMSA([DNA('A' * axis_len)],
                              positional_metadata=positional_metadata)
        self._positional_metadata_constructor_ = factory


class TestTabularMSA(unittest.TestCase, ReallyEqualMixin):
    def test_from_dict_empty(self):
        self.assertEqual(TabularMSA.from_dict({}), TabularMSA([], index=[]))

    def test_from_dict_single_sequence(self):
        self.assertEqual(TabularMSA.from_dict({'foo': DNA('ACGT')}),
                         TabularMSA([DNA('ACGT')], index=['foo']))

    def test_from_dict_multiple_sequences(self):
        msa = TabularMSA.from_dict(
            {1: DNA('ACG'), 2: DNA('GGG'), 3: DNA('TAG')})
        # Sort because order is arbitrary.
        msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACG'), DNA('GGG'), DNA('TAG')], index=[1, 2, 3]))

    def test_from_dict_invalid_input(self):
        # Basic test to make sure error-checking in the TabularMSA constructor
        # is being invoked.
        with self.assertRaisesRegex(
                ValueError, r'must match the number of positions'):
            TabularMSA.from_dict({'a': DNA('ACG'), 'b': DNA('ACGT')})

    def test_constructor_invalid_dtype(self):
        with self.assertRaisesRegex(TypeError, r'GrammaredSequence.*Sequence'):
            TabularMSA([Sequence('')])

        with self.assertRaisesRegex(TypeError, r'GrammaredSequence.*int'):
            TabularMSA([42, DNA('')])

    def test_constructor_not_monomorphic(self):
        with self.assertRaisesRegex(TypeError,
                                    r'matching type.*RNA.*DNA'):
            TabularMSA([DNA(''), RNA('')])

        with self.assertRaisesRegex(TypeError,
                                    r'matching type.*float.*Protein'):
            TabularMSA([Protein(''), Protein(''), 42.0, Protein('')])

    def test_constructor_unequal_length(self):
        with self.assertRaisesRegex(
                ValueError,
                r'must match the number of positions.*1 != 0'):
            TabularMSA([Protein(''), Protein('P')])

        with self.assertRaisesRegex(
                ValueError,
                r'must match the number of positions.*1 != 3'):
            TabularMSA([Protein('PAW'), Protein('ABC'), Protein('A')])

    def test_constructor_non_iterable(self):
        with self.assertRaises(TypeError):
            TabularMSA(42)

    def test_constructor_minter_and_index_both_provided(self):
        with self.assertRaisesRegex(ValueError, r'both.*minter.*index'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str,
                       index=['a', 'b'])

    def test_constructor_invalid_minter_callable(self):
        with self.assertRaises(TypeError):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=float)

    def test_constructor_missing_minter_metadata_key(self):
        with self.assertRaises(KeyError):
            TabularMSA([DNA('ACGT', metadata={'foo': 'bar'}), DNA('TGCA')],
                       minter='foo')

    def test_constructor_unhashable_minter_metadata_key(self):
        with self.assertRaises(TypeError):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=[])

    def test_constructor_index_length_mismatch_iterable(self):
        with self.assertRaisesRegex(ValueError,
                                    r'sequences.*2.*index length.*0'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], index=iter([]))

    def test_constructor_index_length_mismatch_index_object(self):
        with self.assertRaisesRegex(ValueError,
                                    r'sequences.*2.*index length.*0'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], index=pd.Index([]))

    def test_constructor_invalid_index_scalar(self):
        with self.assertRaises(TypeError):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], index=42)

    def test_constructor_non_unique_labels(self):
        msa = TabularMSA([DNA('ACGT'), DNA('ACGT')], index=[1, 1])

        assert_index_equal(msa.index, pd.Index([1, 1], dtype=np.int64))

    def test_constructor_empty_no_index(self):
        # sequence empty
        msa = TabularMSA([])
        self.assertIsNone(msa.dtype)
        self.assertEqual(msa.shape, (0, 0))
        assert_index_equal(msa.index, pd.RangeIndex(0))
        with self.assertRaises(StopIteration):
            next(iter(msa))

        # position empty
        seqs = [DNA(''), DNA('')]
        msa = TabularMSA(seqs)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (2, 0))
        assert_index_equal(msa.index, pd.RangeIndex(2))
        self.assertEqual(list(msa), seqs)

    def test_constructor_empty_with_labels(self):
        # sequence empty
        msa = TabularMSA([], minter=lambda x: x)
        assert_index_equal(msa.index, pd.Index([]))

        msa = TabularMSA([], index=iter([]))
        assert_index_equal(msa.index, pd.Index([]))

        # position empty
        msa = TabularMSA([DNA('', metadata={'id': 42}),
                          DNA('', metadata={'id': 43})], minter='id')
        assert_index_equal(msa.index, pd.Index([42, 43]))

        msa = TabularMSA([DNA(''), DNA('')], index=iter([42, 43]))
        assert_index_equal(msa.index, pd.Index([42, 43]))

    def test_constructor_non_empty_no_labels_provided(self):
        # 1x3
        seqs = [DNA('ACG')]
        msa = TabularMSA(seqs)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (1, 3))
        assert_index_equal(msa.index, pd.RangeIndex(1))
        self.assertEqual(list(msa), seqs)

        # 3x1
        seqs = [DNA('A'), DNA('C'), DNA('G')]
        msa = TabularMSA(seqs)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (3, 1))
        assert_index_equal(msa.index, pd.RangeIndex(3))
        self.assertEqual(list(msa), seqs)

    def test_constructor_non_empty_with_labels_provided(self):
        seqs = [DNA('ACG'), DNA('CGA'), DNA('GTT')]
        msa = TabularMSA(seqs, minter=str)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (3, 3))
        assert_index_equal(msa.index, pd.Index(['ACG', 'CGA', 'GTT']))
        self.assertEqual(list(msa), seqs)

        msa = TabularMSA(seqs, index=iter([42, 43, 44]))
        assert_index_equal(msa.index, pd.Index([42, 43, 44]))

    def test_constructor_works_with_iterator(self):
        seqs = [DNA('ACG'), DNA('CGA'), DNA('GTT')]
        msa = TabularMSA(iter(seqs), minter=str)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (3, 3))
        assert_index_equal(msa.index, pd.Index(['ACG', 'CGA', 'GTT']))
        self.assertEqual(list(msa), seqs)

    def test_constructor_with_multiindex_index(self):
        msa = TabularMSA([DNA('AA'), DNA('GG')],
                         index=[('foo', 42), ('bar', 43)])

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(msa.index, pd.Index([('foo', 42), ('bar', 43)]))

    def test_constructor_with_multiindex_minter(self):
        def multiindex_minter(seq):
            if str(seq) == 'AC':
                return ('foo', 42)
            else:
                return ('bar', 43)

        msa = TabularMSA([DNA('AC'), DNA('GG')], minter=multiindex_minter)

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(msa.index, pd.Index([('foo', 42), ('bar', 43)]))

    def test_copy_constructor_respects_default_index(self):
        msa = TabularMSA([DNA('ACGT'), DNA('----'), DNA('AAAA')])

        copy = TabularMSA(msa)

        self.assertEqual(msa, copy)
        self.assertIsNot(msa, copy)
        assert_index_equal(msa.index, pd.RangeIndex(3))
        assert_index_equal(copy.index, pd.RangeIndex(3))

    def test_copy_constructor_without_metadata(self):
        msa = TabularMSA([DNA('ACGT'), DNA('----')])

        copy = TabularMSA(msa)

        self.assertEqual(msa, copy)
        self.assertIsNot(msa, copy)
        assert_index_equal(copy.index, pd.RangeIndex(2))

    def test_copy_constructor_with_metadata(self):
        msa = TabularMSA([DNA('ACGT'),
                          DNA('----')],
                         metadata={'foo': 42},
                         positional_metadata={'bar': range(4)},
                         index=['idx1', 'idx2'])

        copy = TabularMSA(msa)

        self.assertEqual(msa, copy)
        self.assertIsNot(msa, copy)
        self.assertIsNot(msa.metadata, copy.metadata)
        self.assertIsNot(msa.positional_metadata, copy.positional_metadata)
        # pd.Index is immutable, no copy necessary.
        self.assertIs(msa.index, copy.index)

    def test_copy_constructor_state_override_with_minter(self):
        msa = TabularMSA([DNA('ACGT'),
                          DNA('----')],
                         metadata={'foo': 42},
                         positional_metadata={'bar': range(4)},
                         index=['idx1', 'idx2'])

        copy = TabularMSA(msa, metadata={'foo': 43},
                          positional_metadata={'bar': range(4, 8)},
                          minter=str)

        self.assertNotEqual(msa, copy)

        self.assertEqual(
            copy,
            TabularMSA([DNA('ACGT'),
                        DNA('----')],
                       metadata={'foo': 43},
                       positional_metadata={'bar': range(4, 8)},
                       minter=str))

    def test_copy_constructor_state_override_with_index(self):
        msa = TabularMSA([DNA('ACGT'),
                          DNA('----')],
                         metadata={'foo': 42},
                         positional_metadata={'bar': range(4)},
                         index=['idx1', 'idx2'])

        copy = TabularMSA(msa, metadata={'foo': 43},
                          positional_metadata={'bar': range(4, 8)},
                          index=['a', 'b'])

        self.assertNotEqual(msa, copy)

        self.assertEqual(
            copy,
            TabularMSA([DNA('ACGT'),
                        DNA('----')],
                       metadata={'foo': 43},
                       positional_metadata={'bar': range(4, 8)},
                       index=['a', 'b']))

    def test_copy_constructor_with_minter_and_index(self):
        msa = TabularMSA([DNA('ACGT'), DNA('----')], index=['idx1', 'idx2'])

        with self.assertRaisesRegex(ValueError, r'both.*minter.*index'):
            TabularMSA(msa, index=['a', 'b'], minter=str)

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

    def test_index_getter_default_index(self):
        msa = TabularMSA([DNA('AC'), DNA('AG'), DNA('AT')])

        assert_index_equal(msa.index, pd.RangeIndex(3))

        # immutable
        with self.assertRaises(TypeError):
            msa.index[1] = 2
        # original state is maintained
        assert_index_equal(msa.index, pd.RangeIndex(3))

    def test_index_getter(self):
        index = TabularMSA([DNA('AC'), DNA('AG'), DNA('AT')], minter=str).index
        self.assertIsInstance(index, pd.Index)
        assert_index_equal(index, pd.Index(['AC', 'AG', 'AT']))

        # immutable
        with self.assertRaises(TypeError):
            index[1] = 'AA'
        # original state is maintained
        assert_index_equal(index, pd.Index(['AC', 'AG', 'AT']))

    def test_index_mixed_type(self):
        msa = TabularMSA([DNA('AC'), DNA('CA'), DNA('AA')],
                         index=['abc', 'd', 42])

        assert_index_equal(msa.index, pd.Index(['abc', 'd', 42]))

    def test_index_setter_empty(self):
        msa = TabularMSA([])
        msa.index = iter([])
        assert_index_equal(msa.index, pd.Index([]))

    def test_index_setter_non_empty(self):
        msa = TabularMSA([DNA('AC'), DNA('AG'), DNA('AT')])
        msa.index = range(3)
        assert_index_equal(msa.index, pd.RangeIndex(3))
        msa.index = range(3, 6)
        assert_index_equal(msa.index, pd.RangeIndex(3, 6))

    def test_index_setter_length_mismatch(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str)
        index = pd.Index(['ACGT', 'TGCA'])
        assert_index_equal(msa.index, index)

        with self.assertRaisesRegex(ValueError, r'Length mismatch.*2.*3'):
            msa.index = iter(['ab', 'cd', 'ef'])

        # original state is maintained
        assert_index_equal(msa.index, index)

    def test_index_setter_non_unique_index(self):
        msa = TabularMSA([RNA('UUU'), RNA('AAA')], minter=str)

        msa.index = ['1', '1']

        self.assertEqual(msa, TabularMSA([RNA('UUU'), RNA('AAA')],
                                         index=['1', '1']))

    def test_index_setter_tuples(self):
        msa = TabularMSA([RNA('UUU'), RNA('AAA')])

        msa.index = [('foo', 42), ('bar', 43)]

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(
            msa.index,
            pd.Index([('foo', 42), ('bar', 43)], tupleize_cols=True))

    def test_index_setter_preserves_range_index(self):
        msa = TabularMSA([RNA('UUU'), RNA('AAA')], minter=str)

        msa.index = pd.RangeIndex(2)

        self.assertEqual(msa, TabularMSA([RNA('UUU'), RNA('AAA')]))
        assert_index_equal(msa.index, pd.RangeIndex(2))

    def test_index_deleter(self):
        msa = TabularMSA([RNA('UUU'), RNA('AAA')], minter=str)
        assert_index_equal(msa.index, pd.Index(['UUU', 'AAA']))

        del msa.index
        assert_index_equal(msa.index, pd.RangeIndex(2))

        # Delete again.
        del msa.index
        assert_index_equal(msa.index, pd.RangeIndex(2))

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
            ([RNA('')], {}),
            ([RNA('')], {'minter': str}),

            # 1x1
            ([RNA('U')], {'minter': str}),

            # 2x3
            ([RNA('AUG'), RNA('GUA')], {'minter': str}),

            ([RNA('AG'), RNA('GG')], {}),
            # has labels
            ([RNA('AG'), RNA('GG')], {'minter': str}),
            # different dtype
            ([DNA('AG'), DNA('GG')], {'minter': str}),
            # different labels
            ([RNA('AG'), RNA('GG')], {'minter': lambda x: str(x) + '42'}),
            # different sequence metadata
            ([RNA('AG', metadata={'id': 42}), RNA('GG')], {'minter': str}),
            # different sequence data, same labels
            ([RNA('AG'), RNA('GA')],
             {'minter': lambda x: 'AG' if 'AG' in x else 'GG'}),
            # different MSA metadata
            ([RNA('AG'), RNA('GG')], {'metadata': {'foo': 42}}),
            ([RNA('AG'), RNA('GG')], {'metadata': {'foo': 43}}),
            ([RNA('AG'), RNA('GG')], {'metadata': {'foo': 42, 'bar': 43}}),
            # different MSA positional metadata
            ([RNA('AG'), RNA('GG')],
             {'positional_metadata': {'foo': [42, 43]}}),
            ([RNA('AG'), RNA('GG')],
             {'positional_metadata': {'foo': [43, 44]}}),
            ([RNA('AG'), RNA('GG')],
             {'positional_metadata': {'foo': [42, 43], 'bar': [43, 44]}}),
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

    def test_eq_ignores_minter_str_and_lambda(self):
        # as long as the labels generated by the minters are the same, it
        # doesn't matter whether the minters are equal.
        msa1 = TabularMSA([DNA('ACGT', metadata={'id': 'a'})], minter='id')
        msa2 = TabularMSA([DNA('ACGT', metadata={'id': 'a'})],
                          minter=lambda x: x.metadata['id'])
        self.assertReallyEqual(msa1, msa2)

    def test_eq_minter_and_index(self):
        # as long as the labels generated by the minters are the same, it
        # doesn't matter whether the minters are equal.
        msa1 = TabularMSA([DNA('ACGT', metadata={'id': 'a'})], index=['a'])
        msa2 = TabularMSA([DNA('ACGT', metadata={'id': 'a'})], minter='id')
        self.assertReallyEqual(msa1, msa2)

    def test_eq_default_index_and_equivalent_provided_index(self):
        msa1 = TabularMSA([DNA('ACGT'), DNA('----'), DNA('....')])
        msa2 = TabularMSA([DNA('ACGT'), DNA('----'), DNA('....')],
                          index=[0, 1, 2])

        self.assertReallyEqual(msa1, msa2)
        assert_index_equal(msa1.index, pd.RangeIndex(3))
        assert_index_equal(msa2.index, pd.Index([0, 1, 2], dtype=np.int64))

    def test_reassign_index_empty(self):
        # sequence empty
        msa = TabularMSA([])
        msa.reassign_index()
        self.assertEqual(msa, TabularMSA([]))
        assert_index_equal(msa.index, pd.RangeIndex(0))

        msa.reassign_index(minter=str)
        self.assertEqual(msa, TabularMSA([], minter=str))
        assert_index_equal(msa.index, pd.Index([]))

        # position empty
        msa = TabularMSA([DNA('')])
        msa.reassign_index()
        self.assertEqual(msa, TabularMSA([DNA('')]))
        assert_index_equal(msa.index, pd.RangeIndex(1))

        msa.reassign_index(minter=str)
        self.assertEqual(msa, TabularMSA([DNA('')], minter=str))
        assert_index_equal(msa.index, pd.Index(['']))

    def test_reassign_index_non_empty(self):
        msa = TabularMSA([DNA('ACG', metadata={'id': 1}),
                          DNA('AAA', metadata={'id': 2})], minter=str)
        assert_index_equal(msa.index, pd.Index(['ACG', 'AAA']))

        msa.reassign_index(minter='id')
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACG', metadata={'id': 1}),
                        DNA('AAA', metadata={'id': 2})], minter='id'))
        assert_index_equal(msa.index, pd.Index([1, 2]))

        msa.reassign_index(mapping={1: 5})
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACG', metadata={'id': 1}),
                        DNA('AAA', metadata={'id': 2})], index=[5, 2]))
        assert_index_equal(msa.index, pd.Index([5, 2]))

        msa.reassign_index()
        assert_index_equal(msa.index, pd.RangeIndex(2))

    def test_reassign_index_minter_and_mapping_both_provided(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str)

        with self.assertRaisesRegex(ValueError, r'both.*mapping.*minter.*'):
            msa.reassign_index(minter=str, mapping={"ACGT": "fleventy"})

        # original state is maintained
        assert_index_equal(msa.index, pd.Index(['ACGT', 'TGCA']))

    def test_reassign_index_mapping_invalid_type(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str)

        with self.assertRaisesRegex(TypeError,
                                    r'mapping.*dict.*callable.*list'):
            msa.reassign_index(mapping=['abc', 'def'])

        # original state is maintained
        assert_index_equal(msa.index, pd.Index(['ACGT', 'TGCA']))

    def test_reassign_index_with_mapping_dict_empty(self):
        seqs = [DNA("A"), DNA("C"), DNA("G")]
        msa = TabularMSA(seqs, index=[0.5, 1.5, 2.5])

        msa.reassign_index(mapping={})
        self.assertEqual(msa, TabularMSA(seqs, index=[0.5, 1.5, 2.5]))

    def test_reassign_index_with_mapping_dict_subset(self):
        seqs = [DNA("A"), DNA("C"), DNA("G")]
        mapping = {0.5: "a", 2.5: "c"}

        msa = TabularMSA(seqs, index=[0.5, 1.5, 2.5])
        msa.reassign_index(mapping=mapping)

        self.assertEqual(msa, TabularMSA(seqs, index=['a', 1.5, 'c']))

    def test_reassign_index_with_mapping_dict_superset(self):
        seqs = [DNA("A"), DNA("C"), DNA("G")]
        mapping = {0.5: "a", 1.5: "b", 2.5: "c", 3.5: "d"}

        msa = TabularMSA(seqs, index=[0.5, 1.5, 2.5])
        msa.reassign_index(mapping=mapping)

        self.assertEqual(msa, TabularMSA(seqs, index=['a', 'b', 'c']))

    def test_reassign_index_with_mapping_callable(self):
        seqs = [DNA("A"), DNA("C"), DNA("G")]

        msa = TabularMSA(seqs, index=[0, 1, 2])
        msa.reassign_index(mapping=str)

        self.assertEqual(msa, TabularMSA(seqs, index=['0', '1', '2']))

        msa.reassign_index(mapping=lambda e: int(e) + 42)

        self.assertEqual(msa, TabularMSA(seqs, index=[42, 43, 44]))

    def test_reassign_index_non_unique_existing_index(self):
        seqs = [DNA("A"), DNA("C"), DNA("G")]
        mapping = {0.5: "a", 1.5: "b", 2.5: "c", 3.5: "d"}

        msa = TabularMSA(seqs, index=[0.5, 0.5, 0.5])
        msa.reassign_index(mapping=mapping)

        self.assertEqual(msa, TabularMSA(seqs, index=['a', 'a', 'a']))

    def test_reassign_index_non_unique_new_index(self):
        seqs = [DNA("A"), DNA("C"), DNA("G")]
        mapping = {0.5: "a", 1.5: "a", 2.5: "a"}

        msa = TabularMSA(seqs, index=[0.5, 1.5, 2.5])
        msa.reassign_index(mapping=mapping)

        self.assertEqual(msa, TabularMSA(seqs, index=['a', 'a', 'a']))

    def test_reassign_index_to_multiindex_with_minter(self):
        msa = TabularMSA([DNA('AC'), DNA('.G')])

        def multiindex_minter(seq):
            if str(seq) == 'AC':
                return ('foo', 42)
            else:
                return ('bar', 43)

        msa.reassign_index(minter=multiindex_minter)

        self.assertIsInstance(msa.index, pd.MultiIndex)
        self.assertEqual(
            msa,
            TabularMSA([DNA('AC'), DNA('.G')],
                       index=[('foo', 42), ('bar', 43)]))

    def test_reassign_index_to_multiindex_with_mapping(self):
        msa = TabularMSA([DNA('AC'), DNA('.G')])
        mapping = {0: ('foo', 42), 1: ('bar', 43)}

        msa.reassign_index(mapping=mapping)

        self.assertIsInstance(msa.index, pd.MultiIndex)
        self.assertEqual(
            msa,
            TabularMSA([DNA('AC'), DNA('.G')],
                       index=[('foo', 42), ('bar', 43)]))

    def test_sort_on_unorderable_msa_index(self):
        msa = TabularMSA([DNA('AAA'), DNA('ACG'), DNA('---')],
                         index=[42, 41, 'foo'])
        with self.assertRaises(TypeError):
            msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([DNA('AAA'), DNA('ACG'), DNA('---')],
                       index=[42, 41, 'foo']))

    def test_sort_empty_on_msa_index(self):
        msa = TabularMSA([], index=[])
        msa.sort()
        self.assertEqual(msa, TabularMSA([], index=[]))

        msa = TabularMSA([], index=[])
        msa.sort(ascending=False)
        self.assertEqual(msa, TabularMSA([], index=[]))

    def test_sort_single_sequence_on_msa_index(self):
        msa = TabularMSA([DNA('ACGT')], index=[42])
        msa.sort()
        self.assertEqual(msa, TabularMSA([DNA('ACGT')], index=[42]))

        msa = TabularMSA([DNA('ACGT')], index=[42])
        msa.sort(ascending=False)
        self.assertEqual(msa, TabularMSA([DNA('ACGT')], index=[42]))

    def test_sort_multiple_sequences_on_msa_index(self):
        msa = TabularMSA([
            DNA('TC'), DNA('GG'), DNA('CC')], index=['z', 'a', 'b'])
        msa.sort(ascending=True)
        self.assertEqual(
            msa,
            TabularMSA([
                DNA('GG'), DNA('CC'), DNA('TC')], index=['a', 'b', 'z']))

        msa = TabularMSA([
            DNA('TC'), DNA('GG'), DNA('CC')], index=['z', 'a', 'b'])
        msa.sort(ascending=False)
        self.assertEqual(
            msa,
            TabularMSA([
                DNA('TC'), DNA('CC'), DNA('GG')], index=['z', 'b', 'a']))

    def test_sort_on_labels_with_some_repeats(self):
        msa = TabularMSA([
            DNA('TCCG', metadata={'id': 10}),
            DNA('TAGG', metadata={'id': 10}),
            DNA('GGGG', metadata={'id': 8}),
            DNA('TGGG', metadata={'id': 10}),
            DNA('ACGT', metadata={'id': 0}),
            DNA('TAGA', metadata={'id': 10})], minter='id')
        msa.sort()
        self.assertEqual(msa._seqs.index.to_list(), [0, 8, 10, 10, 10, 10])
        vals = list(msa._seqs.values)
        self.assertEqual(vals[0], DNA('ACGT', metadata={'id': 0}))
        self.assertEqual(vals[1], DNA('GGGG', metadata={'id': 8}))
        self.assertIn(DNA('TCCG', metadata={'id': 10}), vals)
        self.assertIn(DNA('TAGG', metadata={'id': 10}), vals[2:])
        self.assertIn(DNA('TGGG', metadata={'id': 10}), vals[2:])
        self.assertIn(DNA('TAGA', metadata={'id': 10}), vals[2:])

    def test_sort_on_key_with_all_repeats(self):
        msa = TabularMSA([
            DNA('TTT', metadata={'id': 'a'}),
            DNA('TTT', metadata={'id': 'b'}),
            DNA('TTT', metadata={'id': 'c'})], minter=str)
        msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([
                DNA('TTT', metadata={'id': 'a'}),
                DNA('TTT', metadata={'id': 'b'}),
                DNA('TTT', metadata={'id': 'c'})], minter=str))

    def test_sort_default_index(self):
        msa = TabularMSA([DNA('TC'), DNA('GG'), DNA('CC')])
        msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([DNA('TC'), DNA('GG'), DNA('CC')]))

    def test_sort_default_index_descending(self):
        msa = TabularMSA([DNA('TC'), DNA('GG'), DNA('CC')])
        msa.sort(ascending=False)
        self.assertEqual(
            msa,
            TabularMSA([DNA('CC'), DNA('GG'), DNA('TC')], index=[2, 1, 0]))

    def test_sort_already_sorted(self):
        msa = TabularMSA([DNA('TC'), DNA('GG'), DNA('CC')], index=[1, 2, 3])
        msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([DNA('TC'), DNA('GG'), DNA('CC')], index=[1, 2, 3]))

        msa = TabularMSA([DNA('TC'), DNA('GG'), DNA('CC')], index=[3, 2, 1])
        msa.sort(ascending=False)
        self.assertEqual(
            msa,
            TabularMSA([DNA('TC'), DNA('GG'), DNA('CC')], index=[3, 2, 1]))

    def test_sort_reverse_sorted(self):
        msa = TabularMSA([DNA('T'), DNA('G'), DNA('A')], index=[3, 2, 1])
        msa.sort()
        self.assertEqual(
            msa,
            TabularMSA([DNA('A'), DNA('G'), DNA('T')], index=[1, 2, 3]))

        msa = TabularMSA([DNA('T'), DNA('G'), DNA('A')], index=[1, 2, 3])
        msa.sort(ascending=False)
        self.assertEqual(
            msa,
            TabularMSA([DNA('A'), DNA('G'), DNA('T')], index=[3, 2, 1]))

    def test_sort_multiindex(self):
        multiindex = [(2, 'a'), (1, 'c'), (3, 'b')]
        sortedindex = [(1, 'c'), (2, 'a'), (3, 'b')]
        msa = TabularMSA([DNA('A'), DNA('C'), DNA('G')], index=multiindex)
        msa.sort()
        self.assertEqual(msa, TabularMSA([DNA('C'), DNA('A'), DNA('G')],
                                         index=sortedindex))

    def test_sort_multiindex_with_level(self):
        multiindex = [(2, 'a'), (1, 'c'), (3, 'b')]
        first_sorted = [(1, 'c'), (2, 'a'), (3, 'b')]
        second_sorted = [(2, 'a'), (3, 'b'), (1, 'c')]

        msa = TabularMSA([DNA('A'), DNA('C'), DNA('G')], index=multiindex)
        self.assertIsInstance(msa.index, pd.MultiIndex)

        msa.sort(level=0)
        self.assertEqual(msa, TabularMSA([DNA('C'), DNA('A'), DNA('G')],
                                         index=first_sorted))
        msa.sort(level=1)
        self.assertEqual(msa, TabularMSA([DNA('A'), DNA('G'), DNA('C')],
                                         index=second_sorted))

    def test_to_dict_falsey_msa(self):
        self.assertEqual(TabularMSA([]).to_dict(), {})
        self.assertEqual(TabularMSA([RNA('')], index=['foo']).to_dict(),
                         {'foo': RNA('')})

    def test_to_dict_non_empty(self):
        seqs = [Protein('PAW', metadata={'id': 42}),
                Protein('WAP', metadata={'id': -999})]
        msa = TabularMSA(seqs, minter='id')
        self.assertEqual(msa.to_dict(), {42: seqs[0], -999: seqs[1]})

    def test_to_dict_default_index(self):
        msa = TabularMSA([RNA('UUA'), RNA('-C-'), RNA('AAA')])

        d = msa.to_dict()

        self.assertEqual(d, {0: RNA('UUA'), 1: RNA('-C-'), 2: RNA('AAA')})

    def test_to_dict_duplicate_labels(self):
        msa = TabularMSA([DNA("A"), DNA("G")], index=[0, 0])

        with self.assertRaises(ValueError) as cm:
            msa.to_dict()

        self.assertIn("unique", str(cm.exception))

    def test_from_dict_to_dict_roundtrip(self):
        d = {}
        self.assertEqual(TabularMSA.from_dict(d).to_dict(), d)

        # can roundtrip even with mixed key types
        d1 = {'a': DNA('CAT'), 42: DNA('TAG')}
        d2 = TabularMSA.from_dict(d1).to_dict()
        self.assertEqual(d2, d1)
        self.assertIs(d1['a'], d2['a'])
        self.assertIs(d1[42], d2[42])
    
    def test_from_path_seqs(self):
        path = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                         states=[0, 2, 0, 6, 0, 1, 0],
                         starts=[0, 0, 0])
        seqs = [DNA("CGGTCGTAACGCGTACA"),
                DNA("CAGGTAAGCATACCTCA"),
                DNA("CGGTCGTCACTGTACACTA")]
        obj = TabularMSA.from_path_seqs(path=path, seqs=seqs)
        self.assertEqual(str(obj[0]), "CGGTCGTAACGCGTA---CA")
        self.assertEqual(str(obj[1]), "CAG--GTAAG-CATACCTCA")
        self.assertEqual(str(obj[2]), "CGGTCGTCAC-TGTACACTA")


class TestContains(unittest.TestCase):
    def test_no_sequences(self):
        msa = TabularMSA([], index=[])

        self.assertFalse('' in msa)
        self.assertFalse('foo' in msa)

    def test_with_str_labels(self):
        msa = TabularMSA([RNA('AU'), RNA('A.')], index=['foo', 'bar'])

        self.assertTrue('foo' in msa)
        self.assertTrue('bar' in msa)
        self.assertFalse('baz' in msa)
        self.assertFalse(0 in msa)

    def test_with_int_labels(self):
        msa = TabularMSA([RNA('AU'), RNA('A.')], index=[42, -1])

        self.assertTrue(42 in msa)
        self.assertTrue(-1 in msa)
        self.assertFalse(0 in msa)
        self.assertFalse('foo' in msa)


class TestCopy(unittest.TestCase):
    # Note: tests for metadata/positional_metadata are in mixin tests above.

    def test_no_sequences(self):
        msa = TabularMSA([])
        msa_copy = copy.copy(msa)

        self.assertEqual(msa, msa_copy)
        self.assertIsNot(msa, msa_copy)
        self.assertIsNot(msa._seqs, msa_copy._seqs)

    def test_with_sequences(self):
        msa = TabularMSA([DNA('ACGT', metadata={'foo': [1]}), DNA('TGCA')])
        msa_copy = copy.copy(msa)

        self.assertEqual(msa, msa_copy)
        self.assertIsNot(msa, msa_copy)
        self.assertIsNot(msa._seqs, msa_copy._seqs)
        self.assertIsNot(msa[0], msa_copy[0])
        self.assertIsNot(msa[1], msa_copy[1])

        msa_copy.append(DNA('AAAA'), reset_index=True)
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACGT', metadata={'foo': [1]}), DNA('TGCA')]))

        msa_copy._seqs[0].metadata['bar'] = 42
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACGT', metadata={'foo': [1]}), DNA('TGCA')]))

        msa_copy._seqs[0].metadata['foo'].append(2)
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACGT', metadata={'foo': [1, 2]}), DNA('TGCA')]))

    def test_with_index(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], index=['foo', 'bar'])
        msa_copy = copy.copy(msa)

        self.assertEqual(msa, msa_copy)
        self.assertIsNot(msa, msa_copy)
        # pd.Index is immutable, no copy necessary.
        self.assertIs(msa.index, msa_copy.index)

        msa_copy.index = [1, 2]
        assert_index_equal(msa_copy.index, pd.Index([1, 2]))
        assert_index_equal(msa.index, pd.Index(['foo', 'bar']))


class TestDeepCopy(unittest.TestCase):
    # Note: tests for metadata/positional_metadata are in mixin tests above.

    def test_no_sequences(self):
        msa = TabularMSA([])
        msa_copy = copy.deepcopy(msa)

        self.assertEqual(msa, msa_copy)
        self.assertIsNot(msa, msa_copy)
        self.assertIsNot(msa._seqs, msa_copy._seqs)

    def test_with_sequences(self):
        msa = TabularMSA([DNA('ACGT', metadata={'foo': [1]}), DNA('TGCA')])
        msa_copy = copy.deepcopy(msa)

        self.assertEqual(msa, msa_copy)
        self.assertIsNot(msa, msa_copy)
        self.assertIsNot(msa._seqs, msa_copy._seqs)
        self.assertIsNot(msa[0], msa_copy[0])
        self.assertIsNot(msa[1], msa_copy[1])

        msa_copy.append(DNA('AAAA'), reset_index=True)
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACGT', metadata={'foo': [1]}), DNA('TGCA')]))

        msa_copy._seqs[0].metadata['bar'] = 42
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACGT', metadata={'foo': [1]}), DNA('TGCA')]))

        msa_copy._seqs[0].metadata['foo'].append(2)
        self.assertEqual(
            msa,
            TabularMSA([DNA('ACGT', metadata={'foo': [1]}), DNA('TGCA')]))

    def test_with_index(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], index=['foo', 'bar'])
        msa_copy = copy.deepcopy(msa)

        self.assertEqual(msa, msa_copy)
        self.assertIsNot(msa, msa_copy)
        # pd.Index is immutable, no copy necessary.
        self.assertIs(msa.index, msa_copy.index)

        msa_copy.index = [1, 2]
        assert_index_equal(msa_copy.index, pd.Index([1, 2]))
        assert_index_equal(msa.index, pd.Index(['foo', 'bar']))


class SharedIndexTests:
    def get(self, obj, indexable):
        raise NotImplementedError()

    def test_tuple_too_big(self):
        with self.assertRaises(ValueError):
            self.get(TabularMSA([]), (None, None, None))

    def test_empty_msa_slice(self):
        msa = TabularMSA([])

        new = self.get(msa, slice(None, None))

        self.assertIsNot(msa, new)
        self.assertEqual(msa, new)

    def test_msa_slice_all_first_axis(self):
        msa = TabularMSA([RNA("AAA", metadata={1: 1}),
                          RNA("AAU", positional_metadata={0: [1, 2, 3]})],
                         metadata={0: 0}, positional_metadata={1: [3, 2, 1]})

        new_slice = self.get(msa, slice(None))
        new_ellipsis = self.get(msa, Ellipsis)

        self.assertIsNot(msa, new_slice)
        for s1, s2 in zip(msa, new_slice):
            self.assertIsNot(s1, s2)
        self.assertEqual(msa, new_slice)

        self.assertIsNot(msa, new_ellipsis)
        for s1, s2 in zip(msa, new_ellipsis):
            self.assertIsNot(s1, s2)
        self.assertEqual(msa, new_ellipsis)

    def test_msa_slice_all_both_axes(self):
        msa = TabularMSA([RNA("AAA", metadata={1: 1}),
                          RNA("AAU", positional_metadata={0: [1, 2, 3]})],
                         metadata={0: 0}, positional_metadata={1: [3, 2, 1]})

        new_slice = self.get(msa, (slice(None), slice(None)))
        new_ellipsis = self.get(msa, (Ellipsis, Ellipsis))

        self.assertIsNot(msa, new_slice)
        for s1, s2 in zip(msa, new_slice):
            self.assertIsNot(s1, s2)
        self.assertEqual(msa, new_slice)

        self.assertIsNot(msa, new_ellipsis)
        for s1, s2 in zip(msa, new_ellipsis):
            self.assertIsNot(s1, s2)
        self.assertEqual(msa, new_ellipsis)

    def test_bool_index_first_axis(self):
        a = DNA("AAA", metadata={1: 1})
        b = DNA("NNN", positional_metadata={1: ['x', 'y', 'z']})
        c = DNA("AAC")
        msa = TabularMSA([a, b, c], metadata={0: 'x'},
                         positional_metadata={0: [1, 2, 3]},
                         index=[True, False, True])

        new = self.get(msa, [True, True, False])

        self.assertEqual(new, TabularMSA([a, b], metadata={0: 'x'},
                                         positional_metadata={0: [1, 2, 3]},
                                         index=[True, False]))

    def test_bool_index_second_axis(self):
        a = DNA("AAA", metadata={1: 1})
        b = DNA("NNN", positional_metadata={1: ['x', 'y', 'z']})
        c = DNA("AAC")
        msa = TabularMSA([a, b, c], metadata={0: 'x'},
                         positional_metadata={0: [1, 2, 3]},
                         index=[True, False, True])

        new = self.get(msa, (Ellipsis, [True, True, False]))

        self.assertEqual(new, TabularMSA([a[0, 1], b[0, 1], c[0, 1]],
                                         metadata={0: 'x'},
                                         positional_metadata={0: [1, 2]},
                                         index=[True, False, True]))

    def test_bool_index_both_axes(self):
        a = DNA("AAA", metadata={1: 1})
        b = DNA("NNN", positional_metadata={1: ['x', 'y', 'z']})
        c = DNA("AAC")
        msa = TabularMSA([a, b, c], metadata={0: 'x'},
                         positional_metadata={0: [1, 2, 3]},
                         index=[True, False, True])

        new = self.get(msa, ([False, True, True], [True, True, False]))

        self.assertEqual(new, TabularMSA([b[0, 1], c[0, 1]],
                                         metadata={0: 'x'},
                                         positional_metadata={0: [1, 2]},
                                         index=[False, True]))

    def test_bool_index_too_big(self):
        msa = TabularMSA([DNA("ABCD"), DNA("GHKM"), DNA("NRST")],
                         index=[False, True, False])

        with self.assertRaises(IndexError):
            self.get(msa, [False, False, False, False])
        with self.assertRaises(IndexError):
            self.get(msa, [True, True, True, True])

        with self.assertRaises(IndexError):
            self.get(msa, (Ellipsis, [True, False, True, False, True]))

        with self.assertRaises(IndexError):
            self.get(msa, ([True, False, True, False],
                           [True, False, True, False, False]))

    def test_bool_index_too_small(self):
        msa = TabularMSA([DNA("ABCD"), DNA("GHKM"), DNA("NRST")],
                         index=[False, True, False])

        with self.assertRaises(IndexError):
            self.get(msa, [False])
        with self.assertRaises(IndexError):
            self.get(msa, [True])

        with self.assertRaises(IndexError):
            self.get(msa, (Ellipsis, [True]))

        with self.assertRaises(IndexError):
            self.get(msa, ([True, False], [True, False, True, False]))

    def test_bad_scalar(self):
        msa = TabularMSA([DNA("ABCD"), DNA("GHKM"), DNA("NRST")])

        with self.assertRaises((KeyError, TypeError)):
            self.get(msa, "foo")

        with self.assertRaises(IndexError):
            self.get(msa, (Ellipsis, "foo"))

    def test_bad_fancy_index(self):
        msa = TabularMSA([DNA("ABCD"), DNA("GHKM"), DNA("NRST")])

        with self.assertRaises((KeyError, TypeError, ValueError)):
            self.get(msa, [0, "foo"])

        with self.assertRaises(IndexError):
            self.get(msa, (Ellipsis, [0, "foo"]))

    def test_asburd_slice(self):
        msa = TabularMSA([DNA("ABCD"), DNA("GHKM"), DNA("NRST")])

        with self.assertRaises(TypeError):
            self.get(msa, {set(1): 0})


class SharedPropertyIndexTests(SharedIndexTests):
    def setUp(self):
        self.combo_msa = TabularMSA([
            DNA('ACGTA', metadata={0: 0},
                positional_metadata={0: [1, 2, 3, 4, 5]}),
            DNA('CGTAC', metadata={1: 1},
                positional_metadata={1: [1, 2, 3, 4, 5]}),
            DNA('GTACG', metadata={2: 2},
                positional_metadata={2: [1, 2, 3, 4, 5]}),
            DNA('TACGT', metadata={3: 3},
                positional_metadata={3: [1, 2, 3, 4, 5]}),
            DNA('ACGTT', metadata={4: 4},
                positional_metadata={4: [1, 2, 3, 4, 5]})
            ], index=list('ABCDE'), metadata={'x': 'x'},
            positional_metadata={'y': [5, 4, 3, 2, 1]})

        """First off, sorry to the next person who has to deal with this.

           The next few tests will try and slice by a bunch of stuff, with
           all combinations. Each element in the two lists is a tuple where
           the first element is the thing to slice with, and the second is
           the equivalent fancy index which describes the same range.

           This lets us describe the results a little more declaratively
           without setting up a thousand tests for each possible combination.
           This does mean the iloc via a fancy index and simple scalar must
           work correctly.
        """
        # This will be overriden for TestLoc because the first axis are labels
        self.combo_first_axis = [
            ([], []),
            (slice(0, 0), []),
            (Ellipsis, [0, 1, 2, 3, 4]),
            (slice(None), [0, 1, 2, 3, 4]),
            (slice(0, 10000), [0, 1, 2, 3, 4]),
            (3, 3),
            (-4, 1),
            ([0], [0]),
            ([2], [2]),
            (slice(1, 3), [1, 2]),
            (slice(3, 0, -1), [3, 2, 1]),
            ([-3, 2, 1], [2, 2, 1]),
            ([-4, -3, -2, -1], [1, 2, 3, 4]),
            (np.array([-3, 2, 1]), [2, 2, 1]),
            ([True, True, False, False, True], [0, 1, 4]),
            (np.array([True, True, False, True, False]), [0, 1, 3]),
            (range(3), [0, 1, 2]),
            ([slice(0, 2), slice(3, 4), 4], [0, 1, 3, 4])
        ]
        # Same in both TestLoc and TestILoc
        self.combo_second_axis = self.combo_first_axis

    def test_combo_single_axis_natural(self):
        for idx, exp in self.combo_first_axis:
            self.assertEqual(self.get(self.combo_msa, idx),
                             self.combo_msa.iloc[exp],
                             msg="%r did not match iloc[%r]" % (idx, exp))

    def test_combo_first_axis_only(self):
        for idx, exp in self.combo_first_axis:
            self.assertEqual(self.get(self.combo_msa, idx, axis=0),
                             self.combo_msa.iloc[exp, ...],
                             msg="%r did not match iloc[%r, ...]" % (idx, exp))

    def test_combo_second_axis_only(self):
        for idx, exp in self.combo_second_axis:
            self.assertEqual(self.get(self.combo_msa, idx, axis=1),
                             self.combo_msa.iloc[..., exp],
                             msg="%r did not match iloc[..., %r]" % (idx, exp))

    def test_combo_both_axes(self):
        for idx1, exp1 in self.combo_first_axis:
            for idx2, exp2 in self.combo_second_axis:
                self.assertEqual(self.get(self.combo_msa, (idx1, idx2)),
                                 self.combo_msa.iloc[exp1, exp2],
                                 msg=("%r did not match iloc[%r, %r]"
                                      % ((idx1, idx2), exp1, exp2)))


class TestLoc(SharedPropertyIndexTests, unittest.TestCase):
    def setUp(self):
        SharedPropertyIndexTests.setUp(self)
        self.combo_first_axis = [
            ([], []),
            (slice('X', "Z"), []),
            ('A', 0),
            ('E', 4),
            (['B'], [1]),
            (np.asarray(['B']), [1]),
            (slice('A', 'C', 2), [0, 2]),
            (slice('C', 'A', -2), [2, 0]),
            (slice('A', 'B'), [0, 1]),
            (slice(None), [0, 1, 2, 3, 4]),
            (slice('A', None), [0, 1, 2, 3, 4]),
            (slice(None, 'C'), [0, 1, 2]),
            (Ellipsis, [0, 1, 2, 3, 4]),
            (self.combo_msa.index, [0, 1, 2, 3, 4]),
            (['B', 'A', 'A', 'C'], [1, 0, 0, 2]),
            (np.asarray(['B', 'A', 'A', 'C']), [1, 0, 0, 2]),
            ([True, False, True, True, False], [0, 2, 3]),
            (np.asarray([True, False, True, True, False]), [0, 2, 3]),
        ]

    def test_forced_axis_returns_copy(self):
        msa = TabularMSA([Protein("EVANTHQMVS"), Protein("EVANTH*MVS")])

        self.assertIsNot(msa.loc(axis=1), msa.loc)

    def test_forced_axis_no_mutate(self):
        msa = TabularMSA([Protein("EVANTHQMVS"), Protein("EVANTH*MVS")])

        self.assertEqual(msa.loc(axis=1)[0], Sequence("EE"))
        self.assertEqual(msa.loc[0], Protein("EVANTHQMVS"))
        self.assertIsNone(msa.loc._axis)

    def get(self, obj, indexable, axis=None):
        if axis is None:
            return obj.loc[indexable]
        else:
            return obj.loc(axis=axis)[indexable]

    def test_complex_single_label(self):
        a = DNA("ACG")
        b = DNA("ACT")
        c = DNA("ACA")
        msa = TabularMSA([a, b, c], index=[('a', 0), ('a', 1), ('b', 0)])

        self.assertIs(a, self.get(msa, (('a', 0),)))
        self.assertIs(b, self.get(msa, (('a', 1),)))
        self.assertIs(c, self.get(msa, (('b', 0),)))

    def test_partial_label(self):
        a = DNA("ACG")
        b = DNA("ACT")
        c = DNA("ACA")
        msa = TabularMSA([a, b, c], index=[('a', 0), ('a', 1), ('b', 0)])
        exp_a = TabularMSA([a, b], index=[0, 1])
        exp_b = TabularMSA([c], index=[0])

        self.assertEqual(self.get(msa, 'a'), exp_a)
        self.assertEqual(self.get(msa, 'b'), exp_b)

    def test_label_not_exists(self):
        msa = TabularMSA([DNA("ACG")], index=['foo'])

        with self.assertRaises(KeyError):
            self.get(msa, 'bar')

    def test_duplicate_index_nonscalar_label(self):
        a = DNA("ACGA", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = DNA("A-GA", metadata={1: 1}, positional_metadata={1: [1, 2, 3, 4]})
        c = DNA("AAGA", metadata={2: 2}, positional_metadata={2: [1, 2, 3, 4]})
        d = DNA("ACCA", metadata={3: 3}, positional_metadata={3: [1, 2, 3, 4]})

        msa = TabularMSA([a, b, c, d], metadata={'x': 'y'},
                         positional_metadata={'z': [1, 2, 3, 4]},
                         index=[0, 0, 1, 2])

        self.assertEqual(self.get(msa, 0),
                         TabularMSA([a, b], metadata={'x': 'y'},
                                    positional_metadata={'z': [1, 2, 3, 4]},
                                    index=[0, 0]))

    def test_duplicate_index_scalar_label(self):
        a = DNA("ACGA", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = DNA("A-GA", metadata={1: 1}, positional_metadata={1: [1, 2, 3, 4]})
        c = DNA("AAGA", metadata={2: 2}, positional_metadata={2: [1, 2, 3, 4]})
        d = DNA("ACCA", metadata={3: 3}, positional_metadata={3: [1, 2, 3, 4]})

        msa = TabularMSA([a, b, c, d], metadata={'x': 'y'},
                         positional_metadata={'z': [1, 2, 3, 4]},
                         index=[0, 0, 1, 2])

        self.assertEqual(self.get(msa, 1), c)

    def test_multiindex_complex(self):
        a = DNA("ACG")
        b = DNA("ACT")
        c = DNA("ACA")
        msa = TabularMSA([a, b, c], index=[('a', 0), ('a', 1), ('b', 0)])
        exp = TabularMSA([a, c], index=[('a', 0), ('b', 0)])

        self.assertEqual(self.get(msa, [('a', 0), ('b', 0)]), exp)

    def test_fancy_index_missing_label(self):
        msa = TabularMSA([DNA("ACG")], index=['foo'])

        with self.assertRaises(KeyError):
            self.get(msa, ['foo', 'bar'])

        with self.assertRaises(KeyError):
            self.get(msa, ['bar'])

    def test_multiindex_fancy_indexing_incomplete_label(self):
        a = RNA("UUAG", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = RNA("UAAG", metadata={1: 0}, positional_metadata={1: [1, 2, 3, 4]})
        c = RNA("UAA-", metadata={2: 0}, positional_metadata={2: [1, 2, 3, 4]})
        d = RNA("UA-G", metadata={3: 0}, positional_metadata={3: [1, 2, 3, 4]})
        msa = TabularMSA([a, b, c, d], metadata={'x': 'y'},
                         positional_metadata={'c': ['a', 'b', 'c', 'd']},
                         index=[('a', 'x', 0), ('a', 'x', 1), ('a', 'y', 2),
                                ('b', 'x', 0)])

        self.assertEqual(self.get(msa, (('a', 'x'), Ellipsis)),
                         TabularMSA([a, b], metadata={'x': 'y'},
                                    positional_metadata={'c': ['a', 'b', 'c',
                                                               'd']},
                                    index=[0, 1]))

    def test_bool_index_scalar_bool_label(self):
        a = DNA("ACGA", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = DNA("A-GA", metadata={1: 1}, positional_metadata={1: [1, 2, 3, 4]})
        c = DNA("AAGA", metadata={2: 2}, positional_metadata={2: [1, 2, 3, 4]})
        d = DNA("ACCA", metadata={3: 3}, positional_metadata={3: [1, 2, 3, 4]})

        msa = TabularMSA([a, b, c, d], metadata={'x': 'y'},
                         positional_metadata={'z': [1, 2, 3, 4]},
                         index=[False, True, False, False])

        self.assertEqual(self.get(msa, True), b)

    def test_bool_index_nonscalar_bool_label(self):
        a = DNA("ACGA", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = DNA("A-GA", metadata={1: 1}, positional_metadata={1: [1, 2, 3, 4]})
        c = DNA("AAGA", metadata={2: 2}, positional_metadata={2: [1, 2, 3, 4]})
        d = DNA("ACCA", metadata={3: 3}, positional_metadata={3: [1, 2, 3, 4]})

        msa = TabularMSA([a, b, c, d], metadata={'x': 'y'},
                         positional_metadata={'z': [1, 2, 3, 4]},
                         index=[False, True, False, True])

        self.assertEqual(self.get(msa, True),
                         TabularMSA([b, d], metadata={'x': 'y'},
                                    positional_metadata={'z': [1, 2, 3, 4]},
                                    index=[True, True]))

    def test_categorical_index_scalar_label(self):
        msa = TabularMSA([RNA("ACUG"), RNA("ACUA"), RNA("AAUG"), RNA("AC-G")],
                         index=pd.CategoricalIndex(['a', 'b', 'b', 'c']))

        self.assertEqual(self.get(msa, 'a'), RNA("ACUG"))

    def test_categorical_index_nonscalar_label(self):
        msa = TabularMSA([RNA("ACUG"), RNA("ACUA"), RNA("AAUG"), RNA("AC-G")],
                         index=pd.CategoricalIndex(['a', 'b', 'b', 'c']))

        self.assertEqual(self.get(msa, 'b'),
                         TabularMSA([RNA("ACUA"), RNA("AAUG")],
                                    index=pd.CategoricalIndex(
                                        ['b', 'b'], categories=['a', 'b', 'c'])
                                    ))

    def test_float_index_out_of_order_slice(self):
        msa = TabularMSA([DNA("ACGG"), DNA("AAGC"), DNA("AAAA"), DNA("ACTC")],
                         index=[0.1, 2.4, 5.1, 2.6])

        with self.assertRaises(KeyError):
            self.get(msa, slice(0.1, 2.7))

        msa.sort()
        result = self.get(msa, slice(0.1, 2.7))

        self.assertEqual(result, TabularMSA([DNA("ACGG"), DNA("AAGC"),
                                             DNA("ACTC")],
                                            index=[0.1, 2.4, 2.6]))

    def test_nonscalar_fancy_index(self):
        msa = TabularMSA([DNA('ACGT'), DNA('ACGT'), DNA('ACGT')],
                         index=[('a', 0, 1), ('a', 1, 1), ('b', 0, 1)])

        with self.assertRaisesRegex(TypeError,
                                    r'tuple.*independent.*MultiIndex'):
            self.get(msa, ['a', 'b'])

    def test_missing_first_nonscalar_fancy_index(self):
        msa = TabularMSA([DNA('ACGT'), DNA('ACGT'), DNA('ACGT')],
                         index=[('a', 0, 1), ('a', 1, 1), ('b', 0, 1)])

        with self.assertRaises(KeyError):
            self.get(msa, ['x', 'a', 'b'])

    def test_tuple_fancy_index(self):
        msa = TabularMSA([DNA('ACGT'), DNA('ACGT'), DNA('ACGT')],
                         index=[('a', 0, 1), ('a', 1, 1), ('b', 0, 1)])

        with self.assertRaisesRegex(TypeError, r'tuple.*pd.MultiIndex.*label'):
            self.get(msa, ((('a', 0, 1), ('b', 0, 1)), Ellipsis))

    def test_non_multiindex_tuple(self):
        msa = TabularMSA([DNA('ACGT'), DNA('ACGT'), DNA('ACGT')])

        with self.assertRaisesRegex(TypeError, r'tuple.*first axis'):
            self.get(msa, ((0, 1), Ellipsis))

    def test_assertion_exists_for_future_failure_of_get_sequence_loc(self):
        # Ideally we wouldn't need this test or the branch, but the most common
        # failure for pandas would be returning a series instead of the value.
        # We should make sure that the user get's an error should this ever
        # happen again. Getting a series of DNA looks pretty weird...
        msa = TabularMSA([DNA('ACGT'), DNA('ACGT'), DNA('ACGT')])

        with self.assertRaises(AssertionError):
            msa._get_sequence_loc_([1, 2])


class TestILoc(SharedPropertyIndexTests, unittest.TestCase):
    def setUp(self):
        SharedPropertyIndexTests.setUp(self)
        self.combo_first_axis = self.combo_second_axis

    def test_forced_axis_returns_copy(self):
        msa = TabularMSA([Protein("EVANTHQMVS"), Protein("EVANTH*MVS")])

        self.assertIsNot(msa.iloc(axis=1), msa.iloc)

    def test_forced_axis_no_mutate(self):
        msa = TabularMSA([Protein("EVANTHQMVS"), Protein("EVANTH*MVS")])

        self.assertEqual(msa.iloc(axis=1)[0], Sequence("EE"))
        self.assertEqual(msa.iloc[0], Protein("EVANTHQMVS"))
        self.assertIsNone(msa.iloc._axis)

    def get(self, obj, indexable, axis=None):
        if axis is None:
            return obj.iloc[indexable]
        else:
            return obj.iloc(axis=axis)[indexable]

    def test_entire_fancy_first_axis(self):
        msa = TabularMSA([
            DNA("ACCA", metadata={'a': 'foo'},
                positional_metadata={'a': [7, 6, 5, 4]}),
            DNA("GGAA", metadata={'b': 'bar'},
                positional_metadata={'b': [3, 4, 5, 6]})
            ], metadata={'c': 'baz'},
            positional_metadata={'foo': [1, 2, 3, 4]})

        new_np_simple = self.get(msa, np.arange(2))
        new_list_simple = self.get(msa, [0, 1])
        new_list_backwards = self.get(msa, [-2, -1])

        self.assertIsNot(msa, new_np_simple)
        self.assertEqual(msa, new_np_simple)

        self.assertIsNot(msa, new_list_simple)
        self.assertEqual(msa, new_list_simple)

        self.assertIsNot(msa, new_list_backwards)
        self.assertEqual(msa, new_list_backwards)

    def test_fancy_entire_second_axis(self):
        msa = TabularMSA([
            DNA("ACCA", metadata={'a': 'foo'},
                positional_metadata={'a': [7, 6, 5, 4]}),
            DNA("GGAA", metadata={'b': 'bar'},
                positional_metadata={'b': [3, 4, 5, 6]})
            ], metadata={'c': 'baz'},
            positional_metadata={'foo': [1, 2, 3, 4]})

        new_np_simple = self.get(msa, (Ellipsis, np.arange(4)))
        new_list_simple = self.get(msa, (Ellipsis, [0, 1, 2, 3]))
        new_list_backwards = self.get(msa, (Ellipsis, [-4, -3, -2, -1]))

        self.assertIsNot(msa, new_np_simple)
        self.assertEqual(msa, new_np_simple)

        self.assertIsNot(msa, new_list_simple)
        self.assertEqual(msa, new_list_simple)

        self.assertIsNot(msa, new_list_backwards)
        self.assertEqual(msa, new_list_backwards)

    def test_fancy_entire_both_axes(self):
        msa = TabularMSA([
            DNA("ACCA", metadata={'a': 'foo'},
                positional_metadata={'a': [7, 6, 5, 4]}),
            DNA("GGAA", metadata={'b': 'bar'},
                positional_metadata={'b': [3, 4, 5, 6]})
            ], metadata={'c': 'baz'},
            positional_metadata={'foo': [1, 2, 3, 4]})

        new_np_simple = self.get(msa, (np.arange(2), np.arange(4)))
        new_list_simple = self.get(msa, ([0, 1], [0, 1, 2, 3]))
        new_list_backwards = self.get(msa, ([-2, -1], [-4, -3, -2, -1]))

        self.assertIsNot(msa, new_np_simple)
        self.assertEqual(msa, new_np_simple)

        self.assertIsNot(msa, new_list_simple)
        self.assertEqual(msa, new_list_simple)

        self.assertIsNot(msa, new_list_backwards)
        self.assertEqual(msa, new_list_backwards)

    def test_fancy_out_of_bound(self):
        with self.assertRaises(IndexError):
            self.get(TabularMSA([DNA('AC')]), [0, 1, 2])

        with self.assertRaises(IndexError):
            self.get(TabularMSA([DNA('AC')]), (Ellipsis, [0, 1, 2]))

    def test_fancy_empty_both_axis(self):
        msa = TabularMSA([DNA("ACGT", metadata={'x': 1}),
                          DNA("TGCA", metadata={'y': 2})], index=list("AB"))

        new_np_simple = self.get(msa, (np.arange(0), np.arange(0)))
        new_list_simple = self.get(msa, ([], []))

        self.assertEqual(TabularMSA([]), new_np_simple)
        self.assertEqual(TabularMSA([]), new_list_simple)

    def test_fancy_standard_first_axis(self):
        a = DNA("ACGT", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = DNA("ACGT", metadata={1: 1}, positional_metadata={1: [1, 2, 3, 4]})
        c = DNA("ACGT", metadata={2: 2}, positional_metadata={2: [1, 2, 3, 4]})

        msa = TabularMSA([a, b, c], metadata={3: 3},
                         positional_metadata={3: [1, 2, 3, 4]})

        self.assertEqual(self.get(msa, [0, 2]),
                         TabularMSA([a, c], metadata={3: 3},
                                    positional_metadata={3: [1, 2, 3, 4]},
                                    index=[0, 2]))

    def test_fancy_standard_second_axis(self):
        a = DNA("ACGT", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = DNA("ACGT", metadata={1: 1}, positional_metadata={1: [1, 2, 3, 4]})
        c = DNA("ACGT", metadata={2: 2}, positional_metadata={2: [1, 2, 3, 4]})

        msa = TabularMSA([a, b, c], metadata={3: 3},
                         positional_metadata={3: [1, 2, 3, 4]})

        self.assertEqual(self.get(msa, (Ellipsis, [0, 2])),
                         TabularMSA([a[0, 2], b[0, 2], c[0, 2]],
                                    metadata={3: 3},
                                    positional_metadata={3: [1, 3]},
                                    index=[0, 1, 2]))

    def test_fancy_standard_both_axes(self):
        a = DNA("ACGT", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = DNA("ACGT", metadata={1: 1}, positional_metadata={1: [1, 2, 3, 4]})
        c = DNA("ACGT", metadata={2: 2}, positional_metadata={2: [1, 2, 3, 4]})

        msa = TabularMSA([a, b, c], metadata={3: 3},
                         positional_metadata={3: [1, 2, 3, 4]})

        self.assertEqual(self.get(msa, ([0, 2], [0, 2])),
                         TabularMSA([a[0, 2], c[0, 2]],
                                    metadata={3: 3},
                                    positional_metadata={3: [1, 3]},
                                    index=[0, 2]))

    def test_fancy_empty_first_axis(self):
        a = DNA("ACGT", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = DNA("ACGT", metadata={1: 1}, positional_metadata={1: [1, 2, 3, 4]})
        c = DNA("ACGT", metadata={2: 2}, positional_metadata={2: [1, 2, 3, 4]})

        msa = TabularMSA([a, b, c], metadata={3: 3},
                         positional_metadata={3: [1, 2, 3, 4]})
        # TODO: Change for #1198
        self.assertEqual(self.get(msa, []),
                         TabularMSA([], metadata={3: 3}))

    def test_fancy_empty_second_axis(self):
        a = DNA("ACGT", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = DNA("ACGT", metadata={1: 1}, positional_metadata={1: [1, 2, 3, 4]})
        c = DNA("ACGT", metadata={2: 2}, positional_metadata={2: [1, 2, 3, 4]})

        msa = TabularMSA([a, b, c], metadata={3: 3},
                         positional_metadata={3: [1, 2, 3, 4]})

        self.assertEqual(self.get(msa, (Ellipsis, [])),
                         TabularMSA([a[0:0], b[0:0], c[0:0]],
                                    metadata={3: 3},
                                    positional_metadata={3: np.array(
                                        [], dtype=np.int64)}))

    def test_fancy_empty_both_axes(self):
        a = DNA("ACGT", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = DNA("ACGT", metadata={1: 1}, positional_metadata={1: [1, 2, 3, 4]})
        c = DNA("ACGT", metadata={2: 2}, positional_metadata={2: [1, 2, 3, 4]})

        msa = TabularMSA([a, b, c], metadata={3: 3},
                         positional_metadata={3: [1, 2, 3, 4]})
        # TODO: Change for #1198
        self.assertEqual(self.get(msa, ([], [])),
                         TabularMSA([], metadata={3: 3}))

    def test_fancy_out_of_bounds_first_axis(self):
        msa = TabularMSA([DNA("ACGT"), DNA("GCAT")])

        with self.assertRaises(IndexError):
            self.get(msa, [10])

        with self.assertRaises(IndexError):
            self.get(msa, [0, 1, 10])

    def test_fancy_out_of_bounds_second_axis(self):
        msa = TabularMSA([DNA("ACGT"), DNA("GCAT")])

        with self.assertRaises(IndexError):
            self.get(msa, (Ellipsis, [10]))

        with self.assertRaises(IndexError):
            self.get(msa, (Ellipsis, [1, 2, 4]))

    def test_get_scalar_first_axis(self):
        a = DNA("AA", metadata={'a': 'foo'}, positional_metadata={'x': [1, 2]})
        b = DNA("GG", metadata={'b': 'bar'}, positional_metadata={'y': [3, 4]})
        msa = TabularMSA([a, b])

        new0 = self.get(msa, 0)
        new1 = self.get(msa, 1)

        self.assertEqual(new0, a)
        self.assertEqual(new1, b)

    def test_get_scalar_second_axis(self):
        a = DNA("AA", metadata={'a': 'foo'}, positional_metadata={'x': [1, 2]})
        b = DNA("GC", metadata={'b': 'bar'}, positional_metadata={'y': [3, 4]})
        msa = TabularMSA([a, b], positional_metadata={'z': [5, 6]})

        new0 = self.get(msa, (Ellipsis, 0))
        new1 = self.get(msa, (Ellipsis, 1))

        self.assertEqual(new0,
                         Sequence("AG", metadata={'z': 5},
                                  positional_metadata={'x': [1, np.nan],
                                                       'y': [np.nan, 3]}))
        self.assertEqual(new1,
                         Sequence("AC", metadata={'z': 6},
                                  positional_metadata={'x': [2, np.nan],
                                                       'y': [np.nan, 4]}))

    def test_scalar_sliced_first_axis(self):
        a = DNA("ACGT", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = DNA("ACGT", metadata={1: 1}, positional_metadata={1: [1, 2, 3, 4]})
        c = DNA("ACGT", metadata={2: 2}, positional_metadata={2: [1, 2, 3, 4]})

        msa = TabularMSA([a, b, c], metadata={3: 3},
                         positional_metadata={3: [1, 2, 3, 4]})

        self.assertEqual(self.get(msa, (1, [1, 3])),
                         DNA("CT", metadata={1: 1},
                             positional_metadata={1: [2, 4]}))

    def test_scalar_sliced_second_axis(self):
        a = DNA("ACGT", metadata={0: 0}, positional_metadata={0: [1, 2, 3, 4]})
        b = DNA("ACGA", metadata={1: 1}, positional_metadata={1: [1, 2, 3, 4]})
        c = DNA("ACGT", metadata={2: 2}, positional_metadata={2: [1, 2, 3, 4]})

        msa = TabularMSA([a, b, c], metadata={3: 3},
                         positional_metadata={3: [1, 2, 3, 4]})

        self.assertEqual(self.get(msa, ([1, 2], 3)),
                         Sequence("AT", metadata={3: 4},
                                  positional_metadata={1: [4, np.nan],
                                                       2: [np.nan, 4]}))

    def test_get_scalar_out_of_bound_first_axis(self):
        a = DNA("AA", metadata={'a': 'foo'}, positional_metadata={'x': [1, 2]})
        b = DNA("GC", metadata={'b': 'bar'}, positional_metadata={'y': [3, 4]})
        msa = TabularMSA([a, b], positional_metadata={'z': [5, 6]})

        with self.assertRaises(IndexError):
            self.get(msa, 3)

    def test_get_scalar_out_of_bound_second_axis(self):
        a = DNA("AA", metadata={'a': 'foo'}, positional_metadata={'x': [1, 2]})
        b = DNA("GC", metadata={'b': 'bar'}, positional_metadata={'y': [3, 4]})
        msa = TabularMSA([a, b], positional_metadata={'z': [5, 6]})

        with self.assertRaises(IndexError):
            self.get(msa, (Ellipsis, 3))


class TestGetItem(SharedIndexTests, unittest.TestCase):
    def get(self, obj, indexable):
        return obj[indexable]

    def test_uses_iloc_not_loc(self):
        a = DNA("ACGA")
        b = DNA("ACGT")
        msa = TabularMSA([a, b], index=[1, 0])

        self.assertIs(msa[0], a)
        self.assertIs(msa[1], b)


class TestConstructor(unittest.TestCase):
    def setUp(self):
        self.seqs = [DNA("ACGT"), DNA("GCTA")]
        self.m = {'x': 'y', 0: 1}
        self.pm = pd.DataFrame({'foo': [1, 2, 3, 4]})
        self.index = pd.Index(['a', 'b'])
        self.msa = TabularMSA(self.seqs, metadata=self.m,
                              positional_metadata=self.pm, index=self.index)

    def test_no_override(self):
        result = self.msa._constructor_()

        self.assertEqual(self.msa, result)

        for seq1, seq2 in zip(result, self.msa):
            self.assertIsNot(seq1, seq2)

        self.assertIsNot(result.metadata, self.msa.metadata)
        self.assertIsNot(result.positional_metadata,
                         self.msa.positional_metadata)

    def test_sequence_override_same_seqs(self):
        result = self.msa._constructor_(sequences=self.seqs)

        self.assertEqual(self.msa, result)

        for seq1, seq2 in zip(result, self.msa):
            self.assertIsNot(seq1, seq2)

        self.assertIsNot(result.metadata, self.msa.metadata)
        self.assertIsNot(result.positional_metadata,
                         self.msa.positional_metadata)

    def test_sequence_override(self):
        seqs = [RNA("ACGU"), RNA("GCUA")]

        result = self.msa._constructor_(sequences=seqs)

        self.assertNotEqual(result, self.msa)
        self.assertEqual(list(result), seqs)
        assert_index_equal(result.index, self.index)
        self.assertEqual(result.metadata, self.m)
        assert_data_frame_almost_equal(result.positional_metadata, self.pm)

    def test_no_override_no_md(self):
        msa = TabularMSA(self.seqs, index=self.index)

        self.assertEqual(msa, msa._constructor_())

    def test_metadata_override(self):
        new_md = {'foo': {'x': 0}}

        result = self.msa._constructor_(metadata=new_md)

        self.assertNotEqual(result, self.msa)
        self.assertEqual(list(result), self.seqs)
        assert_index_equal(result.index, self.index)
        self.assertEqual(result.metadata, new_md)
        assert_data_frame_almost_equal(result.positional_metadata, self.pm)

    def test_positional_metadata_override(self):
        new_pm = pd.DataFrame({'x': [1, 2, 3, 4]})

        result = self.msa._constructor_(positional_metadata=new_pm)

        self.assertNotEqual(result, self.msa)
        self.assertEqual(list(result), self.seqs)
        assert_index_equal(result.index, self.index)
        self.assertEqual(result.metadata, self.m)
        assert_data_frame_almost_equal(result.positional_metadata, new_pm)

    def test_index_override(self):
        new_index = pd.Index([('a', 0), ('b', 1)])

        result = self.msa._constructor_(index=new_index)

        self.assertNotEqual(result, self.msa)
        self.assertEqual(list(result), self.seqs)
        assert_index_equal(result.index, new_index)
        self.assertEqual(result.metadata, self.m)
        assert_data_frame_almost_equal(result.positional_metadata, self.pm)


class TestAppend(unittest.TestCase):
    # Error cases
    def test_invalid_minter_index_reset_index_parameter_combos(self):
        msa = TabularMSA([])

        param_combos = (
            {},
            {'minter': str, 'index': 'foo', 'reset_index': True},
            {'minter': str, 'index': 'foo'},
            {'minter': str, 'reset_index': True},
            {'index': 'foo', 'reset_index': True}
        )

        for params in param_combos:
            with self.assertRaisesRegex(ValueError,
                                        r"one of.*minter.*index.*reset_index"):
                msa.append(DNA('ACGT'), **params)

            self.assertEqual(msa, TabularMSA([]))

    def test_invalid_dtype(self):
        msa = TabularMSA([])

        with self.assertRaisesRegex(TypeError, r'GrammaredSequence.*Sequence'):
            msa.append(Sequence(''), reset_index=True)

        self.assertEqual(msa, TabularMSA([]))

    def test_dtype_mismatch_rna(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        with self.assertRaisesRegex(TypeError, r'matching type.*RNA.*DNA'):
            msa.append(RNA('UUUU'), reset_index=True)

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_dtype_mismatch_float(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        with self.assertRaisesRegex(TypeError, r'matching type.*float.*DNA'):
            msa.append(42.0, reset_index=True)

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_length_mismatch(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        with self.assertRaisesRegex(
                ValueError, r'must match the number of positions.*5 != 4'):
            msa.append(DNA('ACGTA'), reset_index=True)

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_invalid_minter(self):
        msa = TabularMSA([DNA('ACGT')], index=['foo'])

        with self.assertRaises(KeyError):
            msa.append(DNA('AAAA'), minter='id')

        self.assertEqual(msa, TabularMSA([DNA('ACGT')], index=['foo']))

    # Valid cases: `minter`
    def test_minter_empty_msa(self):
        msa = TabularMSA([])

        msa.append(DNA('ACGT'), minter=str)

        self.assertEqual(msa, TabularMSA([DNA('ACGT')], minter=str))

    def test_minter_metadata_key(self):
        msa = TabularMSA([DNA('', metadata={'id': 'a'}),
                          DNA('', metadata={'id': 'b'})],
                         minter='id')

        msa.append(DNA('', metadata={'id': 'c'}), minter='id')

        self.assertEqual(
            msa,
            TabularMSA([
                DNA('', metadata={'id': 'a'}),
                DNA('', metadata={'id': 'b'}),
                DNA('', metadata={'id': 'c'})], minter='id'))

    def test_minter_callable(self):
        msa = TabularMSA([DNA('', metadata={'id': 'a'}),
                          DNA('', metadata={'id': 'b'})],
                         minter='id')

        msa.append(DNA(''), minter=str)

        self.assertEqual(
            msa,
            TabularMSA([
                DNA('', metadata={'id': 'a'}),
                DNA('', metadata={'id': 'b'}),
                DNA('')], index=['a', 'b', '']))

    def test_multiindex_minter_empty_msa(self):
        def multiindex_minter(seq):
            return ('foo', 42)

        msa = TabularMSA([])

        msa.append(DNA('AC'), minter=multiindex_minter)

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(msa.index, pd.Index([('foo', 42)]))

    def test_multiindex_minter_non_empty_msa(self):
        def multiindex_minter(seq):
            return ('baz', 44)

        msa = TabularMSA([RNA('UU'), RNA('CA')],
                         index=[('foo', 42), ('bar', 43)])

        msa.append(RNA('AC'), minter=multiindex_minter)

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(msa.index,
                           pd.Index([('foo', 42), ('bar', 43), ('baz', 44)]))

    # Valid cases: `index`
    def test_index_empty_msa(self):
        msa = TabularMSA([])

        msa.append(DNA('ACGT'), index='a')

        self.assertEqual(
            msa,
            TabularMSA([DNA('ACGT')], index=['a']))

    def test_index_non_empty_msa(self):
        msa = TabularMSA([DNA('AC'), DNA('GT')], index=['a', 'b'])

        msa.append(DNA('--'), index='foo')

        self.assertEqual(
            msa,
            TabularMSA([DNA('AC'), DNA('GT'), DNA('--')],
                       index=['a', 'b', 'foo']))

    def test_multiindex_index_empty_msa(self):
        msa = TabularMSA([])

        msa.append(DNA('AA'), index=('foo', 42))

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(msa.index, pd.Index([('foo', 42)]))

    def test_multiindex_index_non_empty_msa(self):
        msa = TabularMSA([RNA('A'), RNA('C')],
                         index=[('foo', 42), ('bar', 43)])

        msa.append(RNA('U'), index=('baz', 44))

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(msa.index,
                           pd.Index([('foo', 42), ('bar', 43), ('baz', 44)]))

    # Valid cases: `reset_index`
    def test_reset_index_empty_msa(self):
        msa = TabularMSA([])

        msa.append(DNA('ACGT'), reset_index=True)

        self.assertEqual(msa, TabularMSA([DNA('ACGT')]))
        assert_index_equal(msa.index, pd.RangeIndex(1))

    def test_reset_index_default_index(self):
        msa = TabularMSA([DNA('ACGT'), DNA('CCCC')])

        msa.append(DNA('ACGT'), reset_index=True)

        self.assertEqual(msa,
                         TabularMSA([DNA('ACGT'), DNA('CCCC'), DNA('ACGT')]))
        assert_index_equal(msa.index, pd.RangeIndex(3))

    def test_reset_index_non_default_index(self):
        msa = TabularMSA([DNA('ACGT'), DNA('CCCC')], index=['foo', 'bar'])

        msa.append(DNA('ACGT'), reset_index=True)

        self.assertEqual(msa,
                         TabularMSA([DNA('ACGT'), DNA('CCCC'), DNA('ACGT')]))
        assert_index_equal(msa.index, pd.RangeIndex(3))

    def test_reset_index_bool_cast(self):
        msa = TabularMSA([RNA('AC'), RNA('UU')], index=[42, 43])

        msa.append(RNA('..'), reset_index='abc')

        self.assertEqual(msa, TabularMSA([RNA('AC'), RNA('UU'), RNA('..')]))
        assert_index_equal(msa.index, pd.RangeIndex(3))

    # Valid cases (misc)
    def test_index_type_change(self):
        msa = TabularMSA([DNA('A'), DNA('.')])

        msa.append(DNA('C'), index='foo')

        self.assertEqual(
            msa,
            TabularMSA([DNA('A'), DNA('.'), DNA('C')], index=[0, 1, 'foo']))

    def test_duplicate_index(self):
        msa = TabularMSA([DNA('A'), DNA('.')], index=['foo', 'bar'])

        msa.append(DNA('C'), index='foo')

        self.assertEqual(
            msa,
            TabularMSA([DNA('A'), DNA('.'), DNA('C')],
                       index=['foo', 'bar', 'foo']))

    def test_empty_msa_with_positional_metadata_no_new_positions(self):
        msa = TabularMSA([], positional_metadata={'foo': []})

        msa.append(DNA(''), reset_index=True)

        self.assertEqual(
            msa,
            TabularMSA([DNA('')], positional_metadata={'foo': []}))

    def test_empty_msa_with_positional_metadata_add_new_positions(self):
        # bug in 0.4.2
        msa = TabularMSA([], positional_metadata={'foo': []})

        msa.append(DNA('AA'), reset_index=True)

        self.assertEqual(
            msa,
            TabularMSA([DNA('AA')]))


class TestExtend(unittest.TestCase):
    # Error cases
    #
    # Note: these tests check that the MSA isn't mutated when an error is
    # raised. Where applicable, the "invalid" sequence is preceded by valid
    # sequence(s) to test one possible (buggy) implementation of `extend`:
    # looping over `sequences` and calling `append`. These tests ensure that
    # valid sequences aren't appended to the MSA before the error is raised.
    def test_invalid_minter_index_reset_index_parameter_combos(self):
        msa = TabularMSA([])

        param_combos = (
            {},
            {'minter': str, 'index': 'foo', 'reset_index': True},
            {'minter': str, 'index': 'foo'},
            {'minter': str, 'reset_index': True},
            {'index': 'foo', 'reset_index': True}
        )

        for params in param_combos:
            with self.assertRaisesRegex(ValueError,
                                        r"one of.*minter.*index.*reset_index"):
                msa.extend([DNA('ACGT')], **params)

            self.assertEqual(msa, TabularMSA([]))

    def test_from_tabular_msa_index_param_still_required(self):
        msa = TabularMSA([DNA('AC'), DNA('TG')])

        with self.assertRaisesRegex(ValueError,
                                    r"one of.*minter.*index.*reset_index"):
            msa.extend(TabularMSA([DNA('GG'), DNA('CC')]))

        self.assertEqual(msa, TabularMSA([DNA('AC'), DNA('TG')]))

    def test_invalid_dtype(self):
        msa = TabularMSA([])

        with self.assertRaisesRegex(TypeError, r'GrammaredSequence.*Sequence'):
            msa.extend([Sequence('')], reset_index=True)

        self.assertEqual(msa, TabularMSA([]))

    def test_dtype_mismatch_rna(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        with self.assertRaisesRegex(TypeError, r'matching type.*RNA.*DNA'):
            msa.extend([DNA('----'), RNA('UUUU')], reset_index=True)

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_dtype_mismatch_float(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        with self.assertRaisesRegex(TypeError, r'matching type.*float.*DNA'):
            msa.extend([DNA('GGGG'), 42.0], reset_index=True)

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_length_mismatch(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        with self.assertRaisesRegex(
                ValueError, r'must match the number of positions.*5 != 4'):
            msa.extend([DNA('TTTT'), DNA('ACGTA')], reset_index=True)

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_invalid_minter(self):
        msa = TabularMSA([DNA('ACGT')], index=['foo'])

        with self.assertRaises(KeyError):
            msa.extend([DNA('AAAA', metadata={'id': 'foo'}),
                        DNA('----')], minter='id')

        self.assertEqual(msa, TabularMSA([DNA('ACGT')], index=['foo']))

    def test_invalid_index(self):
        msa = TabularMSA([DNA('ACGT')], index=['foo'])

        with self.assertRaises(TypeError):
            msa.extend([DNA('----')], index=42)

        self.assertEqual(msa, TabularMSA([DNA('ACGT')], index=['foo']))

    def test_sequences_index_length_mismatch(self):
        msa = TabularMSA([])

        with self.assertRaisesRegex(ValueError,
                                    r'sequences.*2.*index length.*3'):
            msa.extend([DNA('TTTT'), DNA('ACGT')], index=['a', 'b', 'c'])

        self.assertEqual(msa, TabularMSA([]))

    # Valid cases: `minter`
    def test_minter_empty_msa(self):
        msa = TabularMSA([])

        msa.extend([RNA('UU'), RNA('--')], minter=str)

        self.assertEqual(msa, TabularMSA([RNA('UU'), RNA('--')], minter=str))

    def test_minter_metadata_key(self):
        msa = TabularMSA([DNA('', metadata={'id': 'a'}),
                          DNA('', metadata={'id': 'b'})],
                         minter='id')

        msa.extend([DNA('', metadata={'id': 'c'}),
                    DNA('', metadata={'id': 'd'})], minter='id')

        self.assertEqual(
            msa,
            TabularMSA([
                DNA('', metadata={'id': 'a'}),
                DNA('', metadata={'id': 'b'}),
                DNA('', metadata={'id': 'c'}),
                DNA('', metadata={'id': 'd'})], minter='id'))

    def test_minter_callable(self):
        msa = TabularMSA([DNA('A', metadata={'id': 'a'}),
                          DNA('C', metadata={'id': 'b'})],
                         minter='id')

        msa.extend([DNA('G'), DNA('T')], minter=str)

        self.assertEqual(
            msa,
            TabularMSA([
                DNA('A', metadata={'id': 'a'}),
                DNA('C', metadata={'id': 'b'}),
                DNA('G'),
                DNA('T')], index=['a', 'b', 'G', 'T']))

    def test_multiindex_minter_empty_msa(self):
        def multiindex_minter(seq):
            if str(seq) == 'AC':
                return ('foo', 42)
            else:
                return ('bar', 43)

        msa = TabularMSA([])

        msa.extend([DNA('AC'), DNA('GG')], minter=multiindex_minter)

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(msa.index, pd.Index([('foo', 42), ('bar', 43)]))

    def test_multiindex_minter_non_empty_msa(self):
        def multiindex_minter(seq):
            if str(seq) == 'C':
                return ('baz', 44)
            else:
                return ('baz', 45)

        msa = TabularMSA([DNA('A'), DNA('G')],
                         index=[('foo', 42), ('bar', 43)])

        msa.extend([DNA('C'), DNA('T')], minter=multiindex_minter)

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(
            msa.index,
            pd.Index([('foo', 42), ('bar', 43), ('baz', 44), ('baz', 45)]))

    # Valid cases: `index`
    def test_index_empty_msa(self):
        msa = TabularMSA([])

        msa.extend([RNA('UAC'), RNA('AAU')], index=['foo', 'bar'])

        self.assertEqual(msa, TabularMSA([RNA('UAC'), RNA('AAU')],
                                         index=['foo', 'bar']))

    def test_index_non_empty_msa(self):
        msa = TabularMSA([DNA('AC'), DNA('GT')], index=['a', 'b'])

        msa.extend([DNA('--'), DNA('..')], index=['foo', 'bar'])

        self.assertEqual(
            msa,
            TabularMSA([DNA('AC'), DNA('GT'), DNA('--'), DNA('..')],
                       index=['a', 'b', 'foo', 'bar']))

    def test_multiindex_index_empty_msa(self):
        msa = TabularMSA([])

        msa.extend([DNA('AA'), DNA('GG')], index=[('foo', 42), ('bar', 43)])

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(msa.index, pd.Index([('foo', 42), ('bar', 43)]))

    def test_multiindex_index_non_empty_msa(self):
        msa = TabularMSA([DNA('.'), DNA('-')],
                         index=[('foo', 42), ('bar', 43)])

        msa.extend([DNA('A'), DNA('G')], index=[('baz', 44), ('baz', 45)])

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(
            msa.index,
            pd.Index([('foo', 42), ('bar', 43), ('baz', 44), ('baz', 45)]))

    def test_index_object_empty_msa(self):
        msa = TabularMSA([])

        msa.extend([DNA('AA'), DNA('GG')], index=pd.RangeIndex(2))

        self.assertEqual(msa, TabularMSA([DNA('AA'), DNA('GG')]))
        assert_index_equal(msa.index, pd.RangeIndex(2))

    def test_index_object_non_empty_msa(self):
        msa = TabularMSA([DNA('CT'), DNA('GG')])

        msa.extend([DNA('AA'), DNA('GG')], index=pd.RangeIndex(2))

        self.assertEqual(
            msa,
            TabularMSA([DNA('CT'), DNA('GG'), DNA('AA'), DNA('GG')],
                       index=[0, 1, 0, 1]))

    # Valid cases: `reset_index`
    def test_reset_index_empty_msa(self):
        msa = TabularMSA([])

        msa.extend([DNA('ACGT'), DNA('----')], reset_index=True)

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('----')]))
        assert_index_equal(msa.index, pd.RangeIndex(2))

    def test_reset_index_empty_msa_empty_iterable(self):
        msa = TabularMSA([])

        msa.extend([], reset_index=True)

        self.assertEqual(msa, TabularMSA([]))
        assert_index_equal(msa.index, pd.RangeIndex(0))

    def test_reset_index_non_empty_msa_empty_iterable(self):
        msa = TabularMSA([RNA('UU'), RNA('CC')], index=['a', 'b'])

        msa.extend([], reset_index=True)

        self.assertEqual(msa, TabularMSA([RNA('UU'), RNA('CC')]))
        assert_index_equal(msa.index, pd.RangeIndex(2))

    def test_reset_index_default_index(self):
        msa = TabularMSA([DNA('A'), DNA('G')])

        msa.extend([DNA('.'), DNA('-')], reset_index=True)

        self.assertEqual(msa,
                         TabularMSA([DNA('A'), DNA('G'), DNA('.'), DNA('-')]))
        assert_index_equal(msa.index, pd.RangeIndex(4))

    def test_reset_index_non_default_index(self):
        msa = TabularMSA([DNA('A'), DNA('G')], index=['a', 'b'])

        msa.extend([DNA('.'), DNA('-')], reset_index=True)

        self.assertEqual(msa,
                         TabularMSA([DNA('A'), DNA('G'), DNA('.'), DNA('-')]))
        assert_index_equal(msa.index, pd.RangeIndex(4))

    def test_reset_index_from_tabular_msa(self):
        msa = TabularMSA([DNA('AC'), DNA('TG')], index=[42, 43])

        msa.extend(TabularMSA([DNA('GG'), DNA('CC'), DNA('AA')],
                              index=['a', 'b', 'c']), reset_index=True)

        self.assertEqual(
            msa,
            TabularMSA([DNA('AC'), DNA('TG'), DNA('GG'), DNA('CC'),
                        DNA('AA')]))
        assert_index_equal(msa.index, pd.RangeIndex(5))

    def test_reset_index_bool_cast(self):
        msa = TabularMSA([RNA('AC'), RNA('UU')], index=[42, 43])

        msa.extend([RNA('..')], reset_index='abc')

        self.assertEqual(msa, TabularMSA([RNA('AC'), RNA('UU'), RNA('..')]))
        assert_index_equal(msa.index, pd.RangeIndex(3))

    # Valid cases (misc)
    def test_index_type_change(self):
        msa = TabularMSA([DNA('A'), DNA('.')])

        msa.extend([DNA('C')], index=['foo'])

        self.assertEqual(
            msa,
            TabularMSA([DNA('A'), DNA('.'), DNA('C')], index=[0, 1, 'foo']))

    def test_duplicate_index(self):
        msa = TabularMSA([DNA('A'), DNA('.')], index=['foo', 'bar'])

        msa.extend([DNA('C'), DNA('.')], index=['foo', 'baz'])

        self.assertEqual(
            msa,
            TabularMSA([DNA('A'), DNA('.'), DNA('C'), DNA('.')],
                       index=['foo', 'bar', 'foo', 'baz']))

    def test_empty_msa_with_positional_metadata_no_new_positions(self):
        msa = TabularMSA([], positional_metadata={'foo': []})

        msa.extend([DNA(''), DNA('')], reset_index=True)

        self.assertEqual(
            msa,
            TabularMSA([DNA(''), DNA('')], positional_metadata={'foo': []}))

    def test_empty_msa_with_positional_metadata_add_new_positions(self):
        # bug in 0.4.2
        msa = TabularMSA([], positional_metadata={'foo': []})

        msa.extend([DNA('AA'), DNA('GG')], reset_index=True)

        self.assertEqual(
            msa,
            TabularMSA([DNA('AA'),
                        DNA('GG')]))

    def test_empty_msa_empty_iterable(self):
        msa = TabularMSA([])

        msa.extend([], minter=str)

        self.assertEqual(msa, TabularMSA([]))

    def test_non_empty_msa_empty_iterable(self):
        msa = TabularMSA([DNA('AC')], index=['foo'])

        msa.extend([], index=[])

        self.assertEqual(msa, TabularMSA([DNA('AC')], index=['foo']))

    def test_single_sequence(self):
        msa = TabularMSA([DNA('AC')])

        msa.extend([DNA('-C')], minter=str)

        self.assertEqual(msa,
                         TabularMSA([DNA('AC'), DNA('-C')], index=[0, '-C']))

    def test_multiple_sequences(self):
        msa = TabularMSA([DNA('AC')])

        msa.extend([DNA('-C'), DNA('AG')], minter=str)

        self.assertEqual(msa,
                         TabularMSA([DNA('AC'), DNA('-C'), DNA('AG')],
                                    index=[0, '-C', 'AG']))

    def test_from_iterable(self):
        msa = TabularMSA([])

        msa.extend(iter([DNA('ACGT'), DNA('TGCA')]), reset_index=True)

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_from_tabular_msa_with_index(self):
        msa1 = TabularMSA([DNA('AC'), DNA('TG')])
        msa2 = TabularMSA([DNA('GG'), DNA('CC'), DNA('AA')])

        msa1.extend(msa2, index=msa2.index)

        self.assertEqual(
            msa1,
            TabularMSA([DNA('AC'), DNA('TG'), DNA('GG'), DNA('CC'),
                        DNA('AA')], index=[0, 1, 0, 1, 2]))


class TestJoin(unittest.TestCase):
    def assertEqualJoinedMSA(self, msa1, msa2):
        # `TabularMSA.join` doesn't guarantee index order in the joined MSA.
        # The order differs across pandas versions, so sort each MSA before
        # comparing for equality.

        # copy because `TabularMSA.sort` is in-place.
        msa1 = copy.copy(msa1)
        msa2 = copy.copy(msa2)
        msa1.sort()
        msa2.sort()

        self.assertEqual(msa1, msa2)

    def test_invalid_how(self):
        with self.assertRaisesRegex(ValueError, r'`how`'):
            TabularMSA([]).join(TabularMSA([]), how='really')

    def test_invalid_other_type(self):
        with self.assertRaisesRegex(TypeError, r'TabularMSA.*DNA'):
            TabularMSA([]).join(DNA('ACGT'))

    def test_dtype_mismatch(self):
        with self.assertRaisesRegex(TypeError, r'dtype.*RNA.*DNA'):
            TabularMSA([DNA('AC')]).join(TabularMSA([RNA('UG')]))

        with self.assertRaisesRegex(TypeError, r'dtype.*None.*DNA'):
            TabularMSA([DNA('AC')]).join(TabularMSA([]))

        with self.assertRaisesRegex(TypeError, r'dtype.*DNA.*None'):
            TabularMSA([]).join(TabularMSA([DNA('AC')]))

    def test_duplicate_index_labels(self):
        with self.assertRaisesRegex(ValueError,
                                    r"This MSA's index labels.*unique"):
            TabularMSA([DNA('AC'), DNA('--')], index=[0, 0]).join(
                TabularMSA([DNA('GT'), DNA('..')]))

        with self.assertRaisesRegex(ValueError,
                                    r"`other`'s index labels.*unique"):
            TabularMSA([DNA('AC'), DNA('--')]).join(
                TabularMSA([DNA('GT'), DNA('..')], index=[0, 0]))

    def test_no_metadata(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.')])
        msa2 = TabularMSA([DNA('-C'),
                           DNA('.G')])

        joined = msa1.join(msa2)

        self.assertEqualJoinedMSA(
            joined,
            TabularMSA([DNA('AC-C'),
                        DNA('G..G')]))

    def test_ignores_metadata(self):
        msa1 = TabularMSA([DNA('AC', metadata={'id': 'a'}),
                           DNA('G.', metadata={'id': 'b'}),
                           DNA('C-', metadata={'id': 'c'})],
                          metadata={'id': 'msa1'})
        msa2 = TabularMSA([DNA('-C', metadata={'id': 'd'}),
                           DNA('.G', metadata={'id': 'e'}),
                           DNA('CA', metadata={'id': 'f'})], index=[2, 1, 0],
                          metadata={'id': 'msa2'})

        joined = msa1.join(msa2)

        self.assertEqualJoinedMSA(
            joined,
            TabularMSA([DNA('ACCA'),
                        DNA('G..G'),
                        DNA('C--C')]))

    def test_outer_join_on_per_sequence_positional_metadata(self):
        msa1 = TabularMSA([
            DNA('AC', positional_metadata={'1': [1, 2], 'foo': ['a', 'b']}),
            DNA('GT', positional_metadata={'2': [3, 4], 'foo': ['c', 'd']})])
        msa2 = TabularMSA([
            DNA('CA', positional_metadata={'3': [5, 6], 'foo': ['e', 'f']}),
            DNA('TG', positional_metadata={'4': [7, 8], 'foo': ['g', 'h']})])

        joined = msa1.join(msa2)

        self.assertEqualJoinedMSA(
            joined,
            TabularMSA([
                DNA('ACCA',
                    positional_metadata={'1': [1, 2, np.nan, np.nan],
                                         '3': [np.nan, np.nan, 5, 6],
                                         'foo': ['a', 'b', 'e', 'f']}),
                DNA('GTTG',
                    positional_metadata={'2': [3, 4, np.nan, np.nan],
                                         '4': [np.nan, np.nan, 7, 8],
                                         'foo': ['c', 'd', 'g', 'h']})]))

    def test_no_sequences(self):
        msa1 = TabularMSA([], positional_metadata={'foo': []})
        msa2 = TabularMSA([], positional_metadata={'foo': []})

        joined = msa1.join(msa2)

        self.assertEqualJoinedMSA(joined, TabularMSA([]))

    def test_no_positions(self):
        msa1 = TabularMSA([DNA('', positional_metadata={'1': []}),
                           DNA('', positional_metadata={'2': []})],
                          positional_metadata={'foo': []})
        msa2 = TabularMSA([DNA('', positional_metadata={'3': []}),
                           DNA('', positional_metadata={'4': []})],
                          positional_metadata={'foo': []})

        joined = msa1.join(msa2)

        self.assertEqualJoinedMSA(
            joined,
            TabularMSA([DNA('', positional_metadata={'1': [], '3': []}),
                        DNA('', positional_metadata={'2': [], '4': []})],
                       positional_metadata={'foo': []}))

    def test_one_with_positions_one_without_positions(self):
        msa1 = TabularMSA([DNA('A', positional_metadata={'1': ['a']}),
                           DNA('C', positional_metadata={'2': ['b']})],
                          positional_metadata={'foo': ['bar']})
        msa2 = TabularMSA([DNA('', positional_metadata={'3': []}),
                           DNA('', positional_metadata={'4': []})],
                          positional_metadata={'foo': []})

        joined = msa1.join(msa2)

        self.assertEqualJoinedMSA(
            joined,
            TabularMSA([DNA('A', positional_metadata={'1': ['a'],
                                                      '3': [np.nan]}),
                        DNA('C', positional_metadata={'2': ['b'],
                                                      '4': [np.nan]})],
                       positional_metadata={'foo': ['bar']}))

    def test_how_strict(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.'),
                           DNA('C-')],
                          positional_metadata={'foo': [1, 2],
                                               'bar': ['a', 'b']})
        msa2 = TabularMSA([DNA('-C'),
                           DNA('.G'),
                           DNA('CA')], index=[2, 1, 0],
                          positional_metadata={'foo': [3, 4],
                                               'bar': ['c', 'd']})

        joined = msa1.join(msa2)

        self.assertEqualJoinedMSA(
            joined,
            TabularMSA([DNA('ACCA'),
                        DNA('G..G'),
                        DNA('C--C')],
                       positional_metadata={'bar': ['a', 'b', 'c', 'd'],
                                            'foo': [1, 2, 3, 4]}))

    def test_how_strict_failure_index_mismatch(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.'),
                           DNA('C-')])
        msa2 = TabularMSA([DNA('-C'),
                           DNA('.G'),
                           DNA('CA'),
                           DNA('--')])

        with self.assertRaisesRegex(ValueError, r'Index labels must all '
                                                'match'):
            msa1.join(msa2)

    def test_how_strict_failure_positional_metadata_mismatch(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.')],
                          positional_metadata={'foo': [1, 2],
                                               'bar': ['a', 'b']})
        msa2 = TabularMSA([DNA('-C'),
                           DNA('.G')],
                          positional_metadata={'foo': [3, 4]})

        with self.assertRaisesRegex(ValueError,
                                    r'Positional metadata columns.*match'):
            msa1.join(msa2)

    def test_how_inner(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.'),
                           DNA('C-'),
                           DNA('--')], index=[0, 1, 2, 3],
                          positional_metadata={'foo': [1, 2],
                                               'bar': ['a', 'b']})
        msa2 = TabularMSA([DNA('-C'),
                           DNA('.G'),
                           DNA('CA'),
                           DNA('..')], index=[2, 1, 0, -1],
                          positional_metadata={'foo': [3, 4],
                                               'baz': ['c', 'd']})

        joined = msa1.join(msa2, how='inner')

        self.assertEqualJoinedMSA(
            joined,
            TabularMSA([DNA('C--C'),
                        DNA('G..G'),
                        DNA('ACCA')], index=[2, 1, 0],
                       positional_metadata={'foo': [1, 2, 3, 4]}))

    def test_how_inner_no_positional_metadata_overlap(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.')], index=['b', 'a'],
                          positional_metadata={'foo': [1, 2]})
        msa2 = TabularMSA([DNA('-C'),
                           DNA('.G')], index=['a', 'b'],
                          positional_metadata={'bar': ['c', 'd']})

        joined = msa1.join(msa2, how='inner')

        self.assertEqualJoinedMSA(
            joined,
            TabularMSA([DNA('G.-C'),
                        DNA('AC.G')], index=['a', 'b']))

    def test_how_inner_no_index_overlap_with_positional_metadata_overlap(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.')],
                          positional_metadata={'foo': [1, 2]})
        msa2 = TabularMSA([DNA('-C'),
                           DNA('.G')], index=['a', 'b'],
                          positional_metadata={'foo': [3, 4]})

        joined = msa1.join(msa2, how='inner')

        self.assertEqualJoinedMSA(joined, TabularMSA([]))

    def test_how_outer(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.'),
                           DNA('C-'),
                           DNA('--')], index=[0, 1, 2, 3],
                          positional_metadata={'foo': [1, 2],
                                               'bar': ['a', 'b']})
        msa2 = TabularMSA([DNA('-CC'),
                           DNA('.GG'),
                           DNA('CAA'),
                           DNA('...')], index=[2, 1, 0, -1],
                          positional_metadata={'foo': [3, 4, 5],
                                               'baz': ['c', 'd', 'e']})

        joined = msa1.join(msa2, how='outer')

        self.assertEqualJoinedMSA(
            joined,
            TabularMSA([DNA('--...'),
                        DNA('ACCAA'),
                        DNA('G..GG'),
                        DNA('C--CC'),
                        DNA('-----')], index=range(-1, 4),
                       positional_metadata={
                           'bar': ['a', 'b', np.nan, np.nan, np.nan],
                           'baz': [np.nan, np.nan, 'c', 'd', 'e'],
                           'foo': [1, 2, 3, 4, 5]}))

    def test_how_left(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.'),
                           DNA('C-'),
                           DNA('--')], index=[0, 1, 2, 3],
                          positional_metadata={'foo': [1, 2],
                                               'bar': ['a', 'b']})
        msa2 = TabularMSA([DNA('-CC'),
                           DNA('.GG'),
                           DNA('CAA'),
                           DNA('...')], index=[2, 1, 0, -1],
                          positional_metadata={'foo': [3, 4, 5],
                                               'baz': ['c', 'd', 'e']})

        joined = msa1.join(msa2, how='left')

        self.assertEqualJoinedMSA(
            joined,
            TabularMSA([DNA('ACCAA'),
                        DNA('G..GG'),
                        DNA('C--CC'),
                        DNA('-----')],
                       positional_metadata={
                           'foo': [1, 2, 3, 4, 5],
                           'bar': ['a', 'b', np.nan, np.nan, np.nan]}))

    def test_how_right(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.'),
                           DNA('C-'),
                           DNA('--')], index=[0, 1, 2, 3],
                          positional_metadata={'foo': [1, 2],
                                               'bar': ['a', 'b']})
        msa2 = TabularMSA([DNA('-CC'),
                           DNA('.GG'),
                           DNA('CAA'),
                           DNA('...')], index=[2, 1, 0, -1],
                          positional_metadata={'foo': [3, 4, 5],
                                               'baz': ['c', 'd', 'e']})

        joined = msa1.join(msa2, how='right')

        self.assertEqualJoinedMSA(
            joined,
            TabularMSA([DNA('C--CC'),
                        DNA('G..GG'),
                        DNA('ACCAA'),
                        DNA('--...')], index=[2, 1, 0, -1],
                       positional_metadata={
                           'foo': [1, 2, 3, 4, 5],
                           'baz': [np.nan, np.nan, 'c', 'd', 'e']}))


class TestIterPositions(unittest.TestCase):
    def test_method_return_type(self):
        msa = TabularMSA([DNA('AC'),
                          DNA('GT')])

        obs = msa.iter_positions()

        self.assertIsInstance(obs, types.GeneratorType)

    def test_position_type(self):
        msa = TabularMSA([DNA('AC'),
                          DNA('GT')])

        first_position = next(msa.iter_positions())

        # Type should be *exactly* Sequence.
        self.assertIs(type(first_position), Sequence)

    def test_no_sequences(self):
        msa = TabularMSA([])

        obs = list(msa.iter_positions())

        self.assertEqual(obs, [])

    def test_no_sequences_ignore_metadata(self):
        msa = TabularMSA([])

        obs = list(msa.iter_positions(ignore_metadata=True))

        self.assertEqual(obs, [])

    def test_no_sequences_reverse(self):
        msa = TabularMSA([])

        obs = list(msa.iter_positions(reverse=True))

        self.assertEqual(obs, [])

    def test_no_sequences_reverse_ignore_metadata(self):
        msa = TabularMSA([])

        obs = list(msa.iter_positions(reverse=True, ignore_metadata=True))

        self.assertEqual(obs, [])

    def test_no_positions(self):
        msa = TabularMSA([DNA(''),
                          DNA('')])

        obs = list(msa.iter_positions())

        self.assertEqual(obs, [])

    def test_no_positions_ignore_metadata(self):
        msa = TabularMSA([DNA(''),
                          DNA('')])

        obs = list(msa.iter_positions(ignore_metadata=True))

        self.assertEqual(obs, [])

    def test_no_positions_reverse(self):
        msa = TabularMSA([DNA(''),
                          DNA('')])

        obs = list(msa.iter_positions(reverse=True))

        self.assertEqual(obs, [])

    def test_no_positions_reverse_ignore_metadata(self):
        msa = TabularMSA([DNA(''),
                          DNA('')])

        obs = list(msa.iter_positions(reverse=True, ignore_metadata=True))

        self.assertEqual(obs, [])

    def test_single_position(self):
        msa = TabularMSA([DNA('A')])

        obs = list(msa.iter_positions())

        self.assertEqual(obs, [Sequence('A')])

    def test_single_position_reverse(self):
        msa = TabularMSA([DNA('A'),
                          DNA('T')])

        obs = list(msa.iter_positions(reverse=True))

        self.assertEqual(obs, [Sequence('AT')])

    def test_multiple_positions(self):
        msa = TabularMSA([DNA('ACGT'),
                          DNA('A-G.'),
                          DNA('----')])

        obs = list(msa.iter_positions())

        self.assertEqual(obs,
                         [Sequence('AA-'), Sequence('C--'), Sequence('GG-'),
                          Sequence('T.-')])

    def test_multiple_positions_reverse(self):
        msa = TabularMSA([DNA('AC'),
                          DNA('A-'),
                          DNA('--')])

        obs = list(msa.iter_positions(reverse=True))

        self.assertEqual(obs,
                         [Sequence('C--'), Sequence('AA-')])

    def test_with_positional_metadata(self):
        # MSA *and* sequence positional metadata.
        msa_positional_metadata = {'pm1': [0.5, 1.5], 'foo': [9, 99]}
        seqs = [
            DNA('AC', positional_metadata={'foo': [42, 43]}),
            DNA('A-'),
            DNA('--', positional_metadata={'foo': [-1, -2],
                                           'bar': ['baz', 'bazz']})]
        msa = TabularMSA(seqs, positional_metadata=msa_positional_metadata)

        obs = list(msa.iter_positions())

        self.assertEqual(
            obs,
            [Sequence('AA-', metadata={'pm1': 0.5, 'foo': 9},
                      positional_metadata={'bar': [np.nan, np.nan, 'baz'],
                                           'foo': [42, np.nan, -1]}),
             Sequence('C--', metadata={'pm1': 1.5, 'foo': 99},
                      positional_metadata={'bar': [np.nan, np.nan, 'bazz'],
                                           'foo': [43, np.nan, -2]})])

    def test_with_positional_metadata_reverse(self):
        # MSA *and* sequence positional metadata.
        msa_positional_metadata = {'pm1': [0.5, 1.5], 'foo': [9, 99]}
        seqs = [
            DNA('AC', positional_metadata={'foo': [42, 43]}),
            DNA('A-'),
            DNA('--', positional_metadata={'foo': [-1, -2],
                                           'bar': ['baz', 'bazz']})]
        msa = TabularMSA(seqs, positional_metadata=msa_positional_metadata)

        obs = list(msa.iter_positions(reverse=True))

        self.assertEqual(
            obs,
            [Sequence('C--', metadata={'pm1': 1.5, 'foo': 99},
                      positional_metadata={'bar': [np.nan, np.nan, 'bazz'],
                                           'foo': [43, np.nan, -2]}),
             Sequence('AA-', metadata={'pm1': 0.5, 'foo': 9},
                      positional_metadata={'bar': [np.nan, np.nan, 'baz'],
                                           'foo': [42, np.nan, -1]})])

    def test_with_positional_metadata_ignore_metadata(self):
        # MSA *and* sequence positional metadata.
        msa_positional_metadata = {'pm1': [0.5, 1.5], 'foo': [9, 99]}
        seqs = [
            DNA('AC', positional_metadata={'foo': [42, 43]}),
            DNA('A-'),
            DNA('--', positional_metadata={'foo': [-1, -2],
                                           'bar': ['baz', 'bazz']})]
        msa = TabularMSA(seqs, positional_metadata=msa_positional_metadata)

        obs = list(msa.iter_positions(ignore_metadata=True))

        self.assertEqual(obs, [Sequence('AA-'), Sequence('C--')])


class TestConsensus(unittest.TestCase):
    def test_no_sequences(self):
        msa = TabularMSA([])

        cons = msa.consensus()

        self.assertEqual(cons, Sequence(''))

    def test_no_positions(self):
        msa = TabularMSA([DNA(''),
                          DNA('')])

        cons = msa.consensus()

        self.assertEqual(cons, DNA(''))

    def test_single_sequence(self):
        msa = TabularMSA([DNA('ACGT-.')])

        cons = msa.consensus()

        self.assertEqual(cons, DNA('ACGT--'))

    def test_multiple_sequences(self):
        msa = TabularMSA([DNA('ACGT'),
                          DNA('AG-.'),
                          DNA('AC-.')])

        cons = msa.consensus()

        self.assertEqual(cons, DNA('AC--'))

    def test_ties(self):
        msa = TabularMSA([DNA('A-'),
                          DNA('C-'),
                          DNA('G-')])

        cons = msa.consensus()

        self.assertTrue(cons in [DNA('A-'), DNA('C-'), DNA('G-')])

    def test_ties_with_gaps(self):
        msa = TabularMSA([DNA('-'),
                          DNA('.'),
                          DNA('T'),
                          DNA('T')])

        cons = msa.consensus()

        self.assertTrue(cons in [DNA('T'), DNA('-')])

    def test_default_gap_char(self):
        msa = TabularMSA([DNA('.'),
                          DNA('.'),
                          DNA('.')])

        cons = msa.consensus()

        self.assertEqual(cons, DNA('-'))

    def test_different_dtype(self):
        msa = TabularMSA([RNA('---'),
                          RNA('AG-'),
                          RNA('AGG')])

        cons = msa.consensus()

        self.assertEqual(cons, RNA('AG-'))

    def test_with_positional_metadata(self):
        # Defining *all* types of metadata to ensure correct metadata is
        # propagated to majority consensus sequence.
        seqs = [
            DNA('-.-', metadata={'id': 'seq1'},
                positional_metadata={'qual': range(0, 3)}),
            DNA('A.T', metadata={'id': 'seq2'},
                positional_metadata={'qual': range(3, 6)}),
            DNA('ACT', metadata={'id': 'seq3'},
                positional_metadata={'qual': range(6, 9)})
        ]
        msa = TabularMSA(seqs, metadata={'pubmed': 123456},
                         positional_metadata={'foo': [42, 43, 42],
                                              'bar': ['a', 'b', 'c']})

        cons = msa.consensus()

        self.assertEqual(
            cons,
            DNA('A-T', positional_metadata={'foo': [42, 43, 42],
                                            'bar': ['a', 'b', 'c']}))

    def test_mixed_gap_characters_as_majority(self):
        seqs = [
            DNA('A'),
            DNA('A'),
            DNA('A'),
            DNA('A'),
            DNA('.'),
            DNA('.'),
            DNA('.'),
            DNA('-'),
            DNA('-')
        ]
        msa = TabularMSA(seqs)

        cons = msa.consensus()

        self.assertEqual(cons, DNA('-'))


class TestConservation(unittest.TestCase):

    def test_no_sequences(self):
        msa = TabularMSA([])
        cons = msa.conservation()
        npt.assert_array_equal(cons, np.array([]))

    def test_shannon_entropy_dna(self):
        msa = TabularMSA([DNA('A'),
                          DNA('G')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.5], base=4)])
        npt.assert_array_equal(actual, expected)

        msa = TabularMSA([DNA('A'),
                          DNA('G'),
                          DNA('C'),
                          DNA('G')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.25, 0.25],
                                                      base=4)])
        npt.assert_array_equal(actual, expected)

        msa = TabularMSA([DNA('AAC'),
                          DNA('GAC')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.5], base=4),
                             1. - scipy.stats.entropy([1.0], base=4),
                             1. - scipy.stats.entropy([1.0], base=4)])
        npt.assert_array_equal(actual, expected)

        msa = TabularMSA([DNA('AACT'),
                          DNA('GACA')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.5], base=4),
                             1. - scipy.stats.entropy([1.0], base=4),
                             1. - scipy.stats.entropy([1.0], base=4),
                             1. - scipy.stats.entropy([0.5, 0.5], base=4)])
        npt.assert_array_equal(actual, expected)

    def test_shannon_entropy_rna(self):
        msa = TabularMSA([RNA('A'),
                          RNA('G')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.5], base=4)])
        npt.assert_array_equal(actual, expected)

        msa = TabularMSA([RNA('A'),
                          RNA('G'),
                          RNA('C'),
                          RNA('G')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.25, 0.25],
                                                      base=4)])
        npt.assert_array_equal(actual, expected)

        msa = TabularMSA([RNA('AAC'),
                          RNA('GAC')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.5], base=4),
                             1. - scipy.stats.entropy([1.0], base=4),
                             1. - scipy.stats.entropy([1.0], base=4)])
        npt.assert_array_equal(actual, expected)

        msa = TabularMSA([RNA('AACU'),
                          RNA('GACA')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.5], base=4),
                             1. - scipy.stats.entropy([1.0], base=4),
                             1. - scipy.stats.entropy([1.0], base=4),
                             1. - scipy.stats.entropy([0.5, 0.5], base=4)])
        npt.assert_array_equal(actual, expected)

    def test_shannon_entropy_protein(self):
        msa = TabularMSA([Protein('A'),
                          Protein('G')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.5], base=22)])
        npt.assert_array_equal(actual, expected)

        msa = TabularMSA([Protein('A'),
                          Protein('G'),
                          Protein('C'),
                          Protein('G')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.25, 0.25],
                                                      base=22)])
        npt.assert_array_equal(actual, expected)

        msa = TabularMSA([Protein('AAC'),
                          Protein('GAC')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.5], base=22),
                             1. - scipy.stats.entropy([1.0], base=22),
                             1. - scipy.stats.entropy([1.0], base=22)])
        npt.assert_array_equal(actual, expected)

        msa = TabularMSA([Protein('AACT'),
                          Protein('GACA')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.5], base=22),
                             1. - scipy.stats.entropy([1.0], base=22),
                             1. - scipy.stats.entropy([1.0], base=22),
                             1. - scipy.stats.entropy([0.5, 0.5], base=22)])
        npt.assert_array_equal(actual, expected)

    def test_degenerate_mode_nan(self):
        msa = TabularMSA([DNA('NAC'),
                          DNA('NNC')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty',
                                  degenerate_mode='nan')
        expected = np.array([np.nan,
                             np.nan,
                             1. - scipy.stats.entropy([1.0], base=4)])
        npt.assert_array_equal(actual, expected)

    def test_degenerate_mode_error(self):
        msa = TabularMSA([DNA('NACN'),
                          DNA('NNCA')])
        self.assertRaises(ValueError, msa.conservation,
                          metric='inverse_shannon_uncertainty',
                          degenerate_mode='error')

        msa = TabularMSA([DNA('AACA'),
                          DNA('ANCA')])
        self.assertRaises(ValueError, msa.conservation,
                          metric='inverse_shannon_uncertainty',
                          degenerate_mode='error')

    def test_error_on_degenerate_w_nan_on_gap(self):
        msa = TabularMSA([DNA('-ACA'),
                          DNA('-NCA')])
        self.assertRaises(ValueError, msa.conservation,
                          metric='inverse_shannon_uncertainty',
                          degenerate_mode='error',
                          gap_mode='nan')

    def test_column_with_degen_and_gap(self):
        msa = TabularMSA([DNA('N'),
                          DNA('-')])
        # test all eight combinations of gap_mode and degenerate_mode
        actual = msa.conservation(metric='inverse_shannon_uncertainty',
                                  degenerate_mode='nan',
                                  gap_mode='nan')
        npt.assert_array_equal(actual, np.array([np.nan]))

        actual = msa.conservation(metric='inverse_shannon_uncertainty',
                                  degenerate_mode='nan',
                                  gap_mode='ignore')
        npt.assert_array_equal(actual, np.array([np.nan]))

        actual = msa.conservation(metric='inverse_shannon_uncertainty',
                                  degenerate_mode='nan',
                                  gap_mode='include')
        npt.assert_array_equal(actual, np.array([np.nan]))

        self.assertRaises(ValueError, msa.conservation,
                          metric='inverse_shannon_uncertainty',
                          degenerate_mode='nan',
                          gap_mode='error')

        self.assertRaises(ValueError, msa.conservation,
                          metric='inverse_shannon_uncertainty',
                          degenerate_mode='error',
                          gap_mode='nan')

        self.assertRaises(ValueError, msa.conservation,
                          metric='inverse_shannon_uncertainty',
                          degenerate_mode='error',
                          gap_mode='error')

        self.assertRaises(ValueError, msa.conservation,
                          metric='inverse_shannon_uncertainty',
                          degenerate_mode='error',
                          gap_mode='include')

        self.assertRaises(ValueError, msa.conservation,
                          metric='inverse_shannon_uncertainty',
                          degenerate_mode='error',
                          gap_mode='ignore')

    def test_gap_mode_nan(self):
        msa = TabularMSA([DNA('-AC.'),
                          DNA('--CA')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty',
                                  gap_mode='nan')
        expected = np.array([np.nan,
                             np.nan,
                             1. - scipy.stats.entropy([1.0], base=4),
                             np.nan])
        npt.assert_array_equal(actual, expected)

    def test_gap_mode_include(self):
        msa = TabularMSA([DNA('AC'),
                          DNA('-G')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty',
                                  gap_mode='include')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.5], base=5),
                             1. - scipy.stats.entropy([0.5, 0.5], base=5)])
        npt.assert_array_equal(actual, expected)

        msa = TabularMSA([DNA('AC'),
                          DNA('.G')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty',
                                  gap_mode='include')
        expected = np.array([1. - scipy.stats.entropy([0.5, 0.5], base=5),
                             1. - scipy.stats.entropy([0.5, 0.5], base=5)])
        npt.assert_array_equal(actual, expected)

    def test_gap_mode_include_gaps_treated_as_single_char(self):
        msa = TabularMSA([DNA('.'),
                          DNA('-')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty',
                                  gap_mode='include')
        expected = np.array([1. - scipy.stats.entropy([1.0], base=5)])
        npt.assert_array_equal(actual, expected)

    def test_gap_mode_ignore(self):
        msa = TabularMSA([DNA('AC'),
                          DNA('-G')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty',
                                  gap_mode='ignore')
        expected = np.array([1. - scipy.stats.entropy([1.0], base=4),
                             1. - scipy.stats.entropy([0.5, 0.5], base=4)])
        npt.assert_array_equal(actual, expected)

        msa = TabularMSA([DNA('AC'),
                          DNA('.G')])
        actual = msa.conservation(metric='inverse_shannon_uncertainty',
                                  gap_mode='ignore')
        expected = np.array([1. - scipy.stats.entropy([1.0], base=4),
                             1. - scipy.stats.entropy([0.5, 0.5], base=4)])
        npt.assert_array_equal(actual, expected)

    def test_gap_mode_error(self):
        msa = TabularMSA([DNA('-AC-'),
                          DNA('--CA')])
        self.assertRaises(ValueError, msa.conservation,
                          metric='inverse_shannon_uncertainty',
                          gap_mode="error")

        msa = TabularMSA([DNA('AACA'),
                          DNA('A-CA')])
        self.assertRaises(ValueError, msa.conservation,
                          metric='inverse_shannon_uncertainty',
                          gap_mode="error")

        msa = TabularMSA([DNA('AACA'),
                          DNA('A.CA')])
        self.assertRaises(ValueError, msa.conservation,
                          metric='inverse_shannon_uncertainty',
                          gap_mode="error")

    def test_bad_metric(self):
        msa = TabularMSA([DNA('AA'),
                          DNA('A-')])
        with self.assertRaisesRegex(ValueError, r'xyz'):
            msa.conservation(metric='xyz')

        msa = TabularMSA([])
        with self.assertRaisesRegex(ValueError, r'xyz'):
            msa.conservation(metric='xyz')

    def test_bad_gap_mode(self):
        msa = TabularMSA([DNA('AA'),
                          DNA('A-')])
        with self.assertRaisesRegex(ValueError, r'xyz'):
            msa.conservation(gap_mode='xyz')

        msa = TabularMSA([])
        with self.assertRaisesRegex(ValueError, r'xyz'):
            msa.conservation(gap_mode='xyz')

    def test_bad_degenerate_mode(self):
        msa = TabularMSA([DNA('AA'),
                          DNA('A-')])
        with self.assertRaisesRegex(ValueError, r'xyz'):
            msa.conservation(degenerate_mode='xyz')

        msa = TabularMSA([])
        with self.assertRaisesRegex(ValueError, r'xyz'):
            msa.conservation(degenerate_mode='xyz')


class TestGapFrequencies(unittest.TestCase):
    def test_default_behavior(self):
        msa = TabularMSA([DNA('AA.'),
                          DNA('-A-')])

        freqs = msa.gap_frequencies()

        npt.assert_array_equal(np.array([1, 0, 2]), freqs)

    def test_invalid_axis_str(self):
        with self.assertRaisesRegex(ValueError, r"axis.*'foo'"):
            TabularMSA([]).gap_frequencies(axis='foo')

    def test_invalid_axis_int(self):
        with self.assertRaisesRegex(ValueError, r"axis.*2"):
            TabularMSA([]).gap_frequencies(axis=2)

    def test_position_axis_str_and_int_equivalent(self):
        msa = TabularMSA([DNA('ACGT'),
                          DNA('A.G-'),
                          DNA('----')])

        str_freqs = msa.gap_frequencies(axis='position')
        int_freqs = msa.gap_frequencies(axis=1)

        npt.assert_array_equal(str_freqs, int_freqs)
        npt.assert_array_equal(np.array([0, 2, 4]), str_freqs)

    def test_sequence_axis_str_and_int_equivalent(self):
        msa = TabularMSA([DNA('ACGT'),
                          DNA('A.G-'),
                          DNA('----')])

        str_freqs = msa.gap_frequencies(axis='sequence')
        int_freqs = msa.gap_frequencies(axis=0)

        npt.assert_array_equal(str_freqs, int_freqs)
        npt.assert_array_equal(np.array([1, 2, 1, 2]), str_freqs)

    def test_correct_dtype_absolute_empty(self):
        msa = TabularMSA([])

        freqs = msa.gap_frequencies(axis='position')

        npt.assert_array_equal(np.array([]), freqs)
        self.assertEqual(int, freqs.dtype)

    def test_correct_dtype_relative_empty(self):
        msa = TabularMSA([])

        freqs = msa.gap_frequencies(axis='position', relative=True)

        npt.assert_array_equal(np.array([]), freqs)
        self.assertEqual(float, freqs.dtype)

    def test_correct_dtype_absolute_non_empty(self):
        msa = TabularMSA([DNA('AC'),
                          DNA('-.')])

        freqs = msa.gap_frequencies(axis='position')

        npt.assert_array_equal(np.array([0, 2]), freqs)
        self.assertEqual(int, freqs.dtype)

    def test_correct_dtype_relative_non_empty(self):
        msa = TabularMSA([DNA('AC'),
                          DNA('-.')])

        freqs = msa.gap_frequencies(axis='position', relative=True)

        npt.assert_array_equal(np.array([0.0, 1.0]), freqs)
        self.assertEqual(float, freqs.dtype)

    def test_no_sequences_absolute(self):
        msa = TabularMSA([])

        seq_freqs = msa.gap_frequencies(axis='sequence')
        pos_freqs = msa.gap_frequencies(axis='position')

        npt.assert_array_equal(np.array([]), seq_freqs)
        npt.assert_array_equal(np.array([]), pos_freqs)

    def test_no_sequences_relative(self):
        msa = TabularMSA([])

        seq_freqs = msa.gap_frequencies(axis='sequence', relative=True)
        pos_freqs = msa.gap_frequencies(axis='position', relative=True)

        npt.assert_array_equal(np.array([]), seq_freqs)
        npt.assert_array_equal(np.array([]), pos_freqs)

    def test_no_positions_absolute(self):
        msa = TabularMSA([DNA('')])

        seq_freqs = msa.gap_frequencies(axis='sequence')
        pos_freqs = msa.gap_frequencies(axis='position')

        npt.assert_array_equal(np.array([]), seq_freqs)
        npt.assert_array_equal(np.array([0]), pos_freqs)

    def test_no_positions_relative(self):
        msa = TabularMSA([DNA('')])

        seq_freqs = msa.gap_frequencies(axis='sequence', relative=True)
        pos_freqs = msa.gap_frequencies(axis='position', relative=True)

        npt.assert_array_equal(np.array([]), seq_freqs)
        npt.assert_array_equal(np.array([np.nan]), pos_freqs)

    def test_single_sequence_absolute(self):
        msa = TabularMSA([DNA('.T')])

        seq_freqs = msa.gap_frequencies(axis='sequence')
        pos_freqs = msa.gap_frequencies(axis='position')

        npt.assert_array_equal(np.array([1, 0]), seq_freqs)
        npt.assert_array_equal(np.array([1]), pos_freqs)

    def test_single_sequence_relative(self):
        msa = TabularMSA([DNA('.T')])

        seq_freqs = msa.gap_frequencies(axis='sequence', relative=True)
        pos_freqs = msa.gap_frequencies(axis='position', relative=True)

        npt.assert_array_equal(np.array([1.0, 0.0]), seq_freqs)
        npt.assert_array_equal(np.array([0.5]), pos_freqs)

    def test_single_position_absolute(self):
        msa = TabularMSA([DNA('.'),
                          DNA('T')])

        seq_freqs = msa.gap_frequencies(axis='sequence')
        pos_freqs = msa.gap_frequencies(axis='position')

        npt.assert_array_equal(np.array([1]), seq_freqs)
        npt.assert_array_equal(np.array([1, 0]), pos_freqs)

    def test_single_position_relative(self):
        msa = TabularMSA([DNA('.'),
                          DNA('T')])

        seq_freqs = msa.gap_frequencies(axis='sequence', relative=True)
        pos_freqs = msa.gap_frequencies(axis='position', relative=True)

        npt.assert_array_equal(np.array([0.5]), seq_freqs)
        npt.assert_array_equal(np.array([1.0, 0.0]), pos_freqs)

    def test_position_axis_absolute(self):
        msa = TabularMSA([
                DNA('ACGT'),   # no gaps
                DNA('A.G-'),   # some gaps (mixed gap chars)
                DNA('----'),   # all gaps
                DNA('....')])  # all gaps

        freqs = msa.gap_frequencies(axis='position')

        npt.assert_array_equal(np.array([0, 2, 4, 4]), freqs)

    def test_position_axis_relative(self):
        msa = TabularMSA([DNA('ACGT'),
                          DNA('A.G-'),
                          DNA('CCC.'),
                          DNA('----'),
                          DNA('....')])

        freqs = msa.gap_frequencies(axis='position', relative=True)

        npt.assert_array_equal(np.array([0.0, 0.5, 0.25, 1.0, 1.0]), freqs)

    def test_sequence_axis_absolute(self):
        msa = TabularMSA([DNA('AC-.'),
                          DNA('A.-.'),
                          DNA('G--.')])

        freqs = msa.gap_frequencies(axis='sequence')

        npt.assert_array_equal(np.array([0, 2, 3, 3]), freqs)

    def test_sequence_axis_relative(self):
        msa = TabularMSA([DNA('AC--.'),
                          DNA('A.A-.'),
                          DNA('G-A-.')])

        freqs = msa.gap_frequencies(axis='sequence', relative=True)

        npt.assert_array_equal(np.array([0.0, 2/3, 1/3, 1.0, 1.0]), freqs)

    def test_relative_frequencies_precise(self):
        class CustomSequence(GrammaredSequence):
            @classproperty
            @overrides(GrammaredSequence)
            def gap_chars(cls):
                return set('0123456789')

            @classproperty
            @overrides(GrammaredSequence)
            def default_gap_char(cls):
                return '0'

            @classproperty
            @overrides(GrammaredSequence)
            def definite_chars(cls):
                return set('')

            @classproperty
            @overrides(GrammaredSequence)
            def degenerate_map(cls):
                return {}

        msa = TabularMSA([CustomSequence('0123456789')])

        freqs = msa.gap_frequencies(axis='position', relative=True)

        npt.assert_array_equal(np.array([1.0]), freqs)

    def test_custom_gap_characters(self):
        class CustomSequence(GrammaredSequence):
            @classproperty
            @overrides(GrammaredSequence)
            def gap_chars(cls):
                return set('#$*')

            @classproperty
            @overrides(GrammaredSequence)
            def default_gap_char(cls):
                return '#'

            @classproperty
            @overrides(GrammaredSequence)
            def definite_chars(cls):
                return set('ABC-.')

            @classproperty
            @overrides(GrammaredSequence)
            def degenerate_map(cls):
                return {'D': 'ABC-.'}

        msa = TabularMSA([CustomSequence('ABCD'),
                          CustomSequence('-.-.'),
                          CustomSequence('A#C*'),
                          CustomSequence('####'),
                          CustomSequence('$$$$')])

        freqs = msa.gap_frequencies(axis='position')

        npt.assert_array_equal(np.array([0, 0, 2, 4, 4]), freqs)


class TestGetPosition(unittest.TestCase):
    def test_without_positional_metadata(self):
        msa = TabularMSA([DNA('ACG'),
                          DNA('A-G')])

        position = msa._get_position_(1)

        self.assertEqual(position, Sequence('C-'))

    def test_with_positional_metadata(self):
        msa = TabularMSA([DNA('ACG'),
                          DNA('A-G')],
                         positional_metadata={'foo': [42, 43, 44],
                                              'bar': ['abc', 'def', 'ghi']})

        position = msa._get_position_(1)

        self.assertEqual(position,
                         Sequence('C-', metadata={'foo': 43, 'bar': 'def'}))


class TestIsSequenceAxis(unittest.TestCase):
    def setUp(self):
        self.msa = TabularMSA([])

    def test_invalid_str(self):
        with self.assertRaisesRegex(ValueError, r"axis.*'foo'"):
            self.msa._is_sequence_axis('foo')

    def test_invalid_int(self):
        with self.assertRaisesRegex(ValueError, r"axis.*2"):
            self.msa._is_sequence_axis(2)

    def test_positive_str(self):
        self.assertTrue(self.msa._is_sequence_axis('sequence'))

    def test_positive_int(self):
        self.assertTrue(self.msa._is_sequence_axis(0))

    def test_negative_str(self):
        self.assertFalse(self.msa._is_sequence_axis('position'))

    def test_negative_int(self):
        self.assertFalse(self.msa._is_sequence_axis(1))


class TestHashable(unittest.TestCase):
    def test_unhashable_type(self):
        self.assertNotIsInstance(TabularMSA([]), collections.abc.Hashable)

    def test_unhashable_object(self):
        with self.assertRaisesRegex(TypeError, r'unhashable'):
            hash(TabularMSA([]))


class TestRepr(unittest.TestCase):
    def test_repr(self):
        # basic sanity checks -- more extensive testing of formatting and
        # special cases is performed in TabularMSAReprDoctests below. here we
        # only test that pieces of the repr are present. these tests also
        # exercise coverage in case doctests stop counting towards coverage in
        # the future

        # str calls repr
        self.assertEqual(repr(TabularMSA([])), str(TabularMSA([])))
        self.assertEqual(repr(TabularMSA([DNA('')])),
                         str(TabularMSA([DNA('')])))
        self.assertEqual(repr(TabularMSA([DNA('ACGT')])),
                         str(TabularMSA([DNA('ACGT')])))
        self.assertEqual(repr(TabularMSA([DNA('ACGT'*25) for x in range(10)])),
                         str(TabularMSA([DNA('ACGT'*25) for x in range(10)])))

        # empty
        obs = repr(TabularMSA([]))
        self.assertEqual(obs.count('\n'), 5)
        self.assertTrue(obs.startswith('TabularMSA'))
        self.assertIn('sequence count: 0', obs)
        self.assertIn('position count: 0', obs)

        # minimal
        obs = repr(TabularMSA([DNA('')]))
        self.assertEqual(obs.count('\n'), 5)
        self.assertTrue(obs.startswith('TabularMSA'))
        self.assertIn('sequence count: 1', obs)
        self.assertIn('position count: 0', obs)
        self.assertIn('[DNA]', obs)

        # no metadata
        obs = repr(TabularMSA([DNA('ACGT')]))
        self.assertEqual(obs.count('\n'), 6)
        self.assertTrue(obs.startswith('TabularMSA'))
        self.assertIn('sequence count: 1', obs)
        self.assertIn('position count: 4', obs)
        self.assertIn('[DNA]', obs)
        self.assertTrue(obs.endswith('ACGT'))

        # sequence spanning > 5 lines
        obs = repr(TabularMSA([DNA('A' * 71) for x in range(6)]))
        self.assertEqual(obs.count('\n'), 10)
        self.assertTrue(obs.startswith('TabularMSA'))
        self.assertIn('sequence count: 6', obs)
        self.assertIn('position count: 71', obs)
        self.assertIn('\n...\n', obs)
        self.assertIn('[DNA]', obs)
        self.assertTrue(obs.endswith('AAAA'))

        # sequences overflowing
        obs = repr(TabularMSA([DNA('A' * 72)]))
        self.assertEqual(obs.count('\n'), 6)
        self.assertTrue(obs.startswith('TabularMSA'))
        self.assertIn('sequence count: 1', obs)
        self.assertIn('position count: 72', obs)
        self.assertIn('[DNA]', obs)
        self.assertTrue(obs.endswith(' ... ' + 'A'*33))


# NOTE: this must be a *separate* class for doctests only (no unit tests). nose
# will not run the unit tests otherwise
# TODO: check if this is still the case since nose is no longer used
#
# these doctests exercise the correct formatting of TabularMSA's repr in a
# variety of situations. they are more extensive than the unit tests above
# (TestRepr.test_repr) but cannot be relied upon for coverage (the unit tests
# take care of this)
class TabularMSAReprDoctests:
    r"""
    >>> from skbio import DNA, TabularMSA

    Empty (minimal) MSA:

    >>> TabularMSA([])
    TabularMSA
    ---------------------
    Stats:
        sequence count: 0
        position count: 0
    ---------------------

    MSA with single empty sequence:

    >>> TabularMSA([DNA('')])
    TabularMSA[DNA]
    ---------------------
    Stats:
        sequence count: 1
        position count: 0
    ---------------------

    MSA with single sequence with single character:

    >>> TabularMSA([DNA('G')])
    TabularMSA[DNA]
    ---------------------
    Stats:
        sequence count: 1
        position count: 1
    ---------------------
    G

    MSA with multicharacter sequence:

    >>> TabularMSA([DNA('ACGT')])
    TabularMSA[DNA]
    ---------------------
    Stats:
        sequence count: 1
        position count: 4
    ---------------------
    ACGT

    Full single line:

    >>> TabularMSA([DNA('A' * 71)])
    TabularMSA[DNA]
    -----------------------------------------------------------------------
    Stats:
        sequence count: 1
        position count: 71
    -----------------------------------------------------------------------
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    Full single line with 1 character overflow:

    >>> TabularMSA([DNA('A' * 72)])
    TabularMSA[DNA]
    -----------------------------------------------------------------------
    Stats:
        sequence count: 1
        position count: 72
    -----------------------------------------------------------------------
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA ... AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    Two sequences with full lines:

    >>> TabularMSA([DNA('T' * 71), DNA('T' * 71)])
    TabularMSA[DNA]
    -----------------------------------------------------------------------
    Stats:
        sequence count: 2
        position count: 71
    -----------------------------------------------------------------------
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

    Two sequences with full lines with 1 character overflow:

    >>> TabularMSA([DNA('T' * 72), DNA('T' * 72)])
    TabularMSA[DNA]
    -----------------------------------------------------------------------
    Stats:
        sequence count: 2
        position count: 72
    -----------------------------------------------------------------------
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT ... TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT ... TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

    Five full lines (maximum amount of information):

    >>> TabularMSA([DNA('A' * 71) for x in range(5)])
    TabularMSA[DNA]
    -----------------------------------------------------------------------
    Stats:
        sequence count: 5
        position count: 71
    -----------------------------------------------------------------------
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    Six lines starts "summarized" output:

    >>> TabularMSA([DNA('A' * 71) for x in range(6)])
    TabularMSA[DNA]
    -----------------------------------------------------------------------
    Stats:
        sequence count: 6
        position count: 71
    -----------------------------------------------------------------------
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    ...
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    Supply horrendous metadata and positional metadata to exercise a variety of
    metadata formatting cases and rules. Sorting should be by type, then by
    value within each type (Python 3 doesn't allow sorting of mixed types):

    >>> metadata = {
    ...     # str key, str value
    ...     'abc': 'some description',
    ...     # int value
    ...     'foo': 42,
    ...     # unsupported type (dict) value
    ...     'bar': {},
    ...     # int key, wrapped str (single line)
    ...     42: 'some words to test text wrapping and such... yada yada yada '
    ...         'yada yada yada yada yada.',
    ...     # bool key, wrapped str (multi-line)
    ...     True: 'abc ' * 34,
    ...     # float key, truncated str (too long)
    ...     42.5: 'abc ' * 200,
    ...     # unsupported type (tuple) key, unsupported type (list) value
    ...     ('foo', 'bar'): [1, 2, 3],
    ...     # bytes key, single long word that wraps
    ...     b'long word': 'abc' * 30,
    ...     # truncated key (too long), None value
    ...     'too long of a key name to display in repr': None,
    ...     # wrapped bytes value (has b'' prefix)
    ...     'bytes wrapped value': b'abcd' * 25,
    ...     # float value
    ...     0.1: 99.9999,
    ...     # bool value
    ...     43: False,
    ...     # None key, complex value
    ...     None: complex(-1.0, 0.0),
    ...     # nested quotes
    ...     10: '"\''
    ... }
    >>> positional_metadata = pd.DataFrame({
    ...     # str key, int list value
    ...     'foo': [1, 2, 3, 4],
    ...     # float key, float list value
    ...     42.5: [2.5, 3.0, 4.2, -0.00001],
    ...     # int key, object list value
    ...     42: [[], 4, 5, {}],
    ...     # truncated key (too long), bool list value
    ...     'abc' * 90: [True, False, False, True],
    ...     # None key
    ...     None: range(4)})
    >>> positional_metadata = positional_metadata.reindex(
    ...     columns=['foo', 42.5, 42, 'abc' * 90, None])
    >>> TabularMSA([DNA('ACGT')], metadata=metadata,
    ...            positional_metadata=positional_metadata)
    TabularMSA[DNA]
    -----------------------------------------------------------------------
    Metadata:
        None: (-1+0j)
        True: 'abc abc abc abc abc abc abc abc abc abc abc abc abc abc abc
               abc abc abc abc abc abc abc abc abc abc abc abc abc abc abc
               abc abc abc abc '
        b'long word': 'abcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabcabca
                       bcabcabcabcabcabcabcabcabcabcabcabcabc'
        0.1: 99.9999
        42.5: <class 'str'>
        10: '"\''
        42: 'some words to test text wrapping and such... yada yada yada
             yada yada yada yada yada.'
        43: False
        'abc': 'some description'
        'bar': <class 'dict'>
        'bytes wrapped value': b'abcdabcdabcdabcdabcdabcdabcdabcdabcdabcdab
                                 cdabcdabcdabcdabcdabcdabcdabcdabcdabcdabcd
                                 abcdabcdabcdabcd'
        'foo': 42
        <class 'str'>: None
        <class 'tuple'>: <class 'list'>
    Positional metadata:
        'foo': <dtype: int64>
        42.5: <dtype: float64>
        42: <dtype: object>
        <class 'str'>: <dtype: bool>
        None: <dtype: int64>
    Stats:
        sequence count: 1
        position count: 4
    -----------------------------------------------------------------------
    ACGT

    """
    pass


if __name__ == "__main__":
    unittest.main()
