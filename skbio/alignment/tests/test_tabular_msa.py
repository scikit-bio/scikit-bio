# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import copy
import unittest
import functools
import itertools
import types

import six
import numpy as np
import numpy.testing as npt
import pandas as pd

from skbio import Sequence, DNA, RNA, Protein, TabularMSA
from skbio.sequence._iupac_sequence import IUPACSequence
from skbio.util._decorator import classproperty, overrides
from skbio.util._testing import (ReallyEqualMixin, MetadataMixinTests,
                                 PositionalMetadataMixinTests,
                                 assert_index_equal)


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
        with six.assertRaisesRegex(
                self, ValueError, 'must match the number of positions'):
            TabularMSA.from_dict({'a': DNA('ACG'), 'b': DNA('ACGT')})

    def test_constructor_invalid_dtype(self):
        with six.assertRaisesRegex(self, TypeError,
                                   'sequence.*alphabet.*Sequence'):
            TabularMSA([Sequence('')])

        with six.assertRaisesRegex(self, TypeError, 'sequence.*alphabet.*int'):
            TabularMSA([42, DNA('')])

    def test_constructor_not_monomorphic(self):
        with six.assertRaisesRegex(self, TypeError,
                                   'matching type.*RNA.*DNA'):
            TabularMSA([DNA(''), RNA('')])

        with six.assertRaisesRegex(self, TypeError,
                                   'matching type.*float.*Protein'):
            TabularMSA([Protein(''), Protein(''), 42.0, Protein('')])

    def test_constructor_unequal_length(self):
        with six.assertRaisesRegex(
                self, ValueError,
                'must match the number of positions.*1 != 0'):
            TabularMSA([Protein(''), Protein('P')])

        with six.assertRaisesRegex(
                self, ValueError,
                'must match the number of positions.*1 != 3'):
            TabularMSA([Protein('PAW'), Protein('ABC'), Protein('A')])

    def test_constructor_non_iterable(self):
        with self.assertRaises(TypeError):
            TabularMSA(42)

    def test_constructor_non_unique_labels(self):
        msa = TabularMSA([DNA('ACGT'), DNA('ACGT')], index=[1, 1])

        assert_index_equal(msa.index, pd.Int64Index([1, 1]))

    def test_constructor_minter_and_index_both_provided(self):
        with six.assertRaisesRegex(self, ValueError, 'both.*minter.*index'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str,
                       index=['a', 'b'])

    def test_constructor_index_length_mismatch_iterable(self):
        with six.assertRaisesRegex(self, ValueError,
                                   'sequences.*2.*index length.*0'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], index=iter([]))

    def test_constructor_index_length_mismatch_index_object(self):
        with six.assertRaisesRegex(self, ValueError,
                                   'sequences.*2.*index length.*0'):
            TabularMSA([DNA('ACGT'), DNA('TGCA')], index=pd.Index([]))

    def test_constructor_empty_no_index(self):
        # sequence empty
        msa = TabularMSA([])
        self.assertIsNone(msa.dtype)
        self.assertEqual(msa.shape, (0, 0))
        assert_index_equal(msa.index, pd.Index([]))
        with self.assertRaises(StopIteration):
            next(iter(msa))

        # position empty
        seqs = [DNA(''), DNA('')]
        msa = TabularMSA(seqs)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (2, 0))
        assert_index_equal(msa.index, pd.Int64Index([0, 1]))
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
        assert_index_equal(msa.index, pd.Index([0]))
        self.assertEqual(list(msa), seqs)

        # 3x1
        seqs = [DNA('A'), DNA('C'), DNA('G')]
        msa = TabularMSA(seqs)
        self.assertIs(msa.dtype, DNA)
        self.assertEqual(msa.shape, (3, 1))
        assert_index_equal(msa.index, pd.Index([0, 1, 2]))
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

    def test_copy_constructor_handles_missing_metadata_efficiently(self):
        msa = TabularMSA([DNA('ACGT'), DNA('----')])

        copy = TabularMSA(msa)

        self.assertIsNone(msa._metadata)
        self.assertIsNone(msa._positional_metadata)
        self.assertIsNone(copy._metadata)
        self.assertIsNone(copy._positional_metadata)

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
        self.assertIsNot(msa.index, copy.index)

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
        assert_index_equal(msa.index, pd.Index([0, 1, 2]))
        msa.index = range(3, 6)
        assert_index_equal(msa.index, pd.Index([3, 4, 5]))

    def test_index_setter_length_mismatch(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str)
        index = pd.Index(['ACGT', 'TGCA'])
        assert_index_equal(msa.index, index)

        with six.assertRaisesRegex(self, ValueError,
                                   'Length mismatch.*2.*3'):
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

        self.assertIsInstance(msa.index, pd.Index)
        self.assertNotIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(
            msa.index,
            pd.Index([('foo', 42), ('bar', 43)], tupleize_cols=False))

    def test_index_deleter(self):
        msa = TabularMSA([RNA('UUU'), RNA('AAA')], minter=str)
        assert_index_equal(msa.index, pd.Index(['UUU', 'AAA']))
        del msa.index
        assert_index_equal(msa.index, pd.Index([0, 1]))

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

    def test_reassign_index_empty(self):
        # sequence empty
        msa = TabularMSA([])
        msa.reassign_index()
        self.assertEqual(msa, TabularMSA([]))
        assert_index_equal(msa.index, pd.Int64Index([]))

        msa.reassign_index(minter=str)
        self.assertEqual(msa, TabularMSA([], minter=str))
        assert_index_equal(msa.index, pd.Index([]))

        # position empty
        msa = TabularMSA([DNA('')])
        msa.reassign_index()
        self.assertEqual(msa, TabularMSA([DNA('')]))
        assert_index_equal(msa.index, pd.Index([0]))

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
        assert_index_equal(msa.index, pd.Index([0, 1]))

    def test_reassign_index_minter_and_mapping_both_provided(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')], minter=str)

        with six.assertRaisesRegex(self, ValueError,
                                   'both.*mapping.*minter.*'):
            msa.reassign_index(minter=str, mapping={"ACGT": "fleventy"})

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

    @unittest.skipIf(six.PY2, "Everything is orderable in Python 2.")
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
        self.assertEqual(
            msa,
            TabularMSA([
                DNA('ACGT', metadata={'id': 0}),
                DNA('GGGG', metadata={'id': 8}),
                DNA('TCCG', metadata={'id': 10}),
                DNA('TAGG', metadata={'id': 10}),
                DNA('TGGG', metadata={'id': 10}),
                DNA('TAGA', metadata={'id': 10})], minter='id'))

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
        # TODO: use __getitem__ when it exists.
        self.assertIsNot(msa._seqs[0], msa_copy._seqs[0])
        self.assertIsNot(msa._seqs[1], msa_copy._seqs[1])

        msa_copy.append(DNA('AAAA'))
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
        self.assertIsNot(msa.index, msa_copy.index)

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
        # TODO: use __getitem__ when it exists.
        self.assertIsNot(msa._seqs[0], msa_copy._seqs[0])
        self.assertIsNot(msa._seqs[1], msa_copy._seqs[1])

        msa_copy.append(DNA('AAAA'))
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
        self.assertIsNot(msa.index, msa_copy.index)

        msa_copy.index = [1, 2]
        assert_index_equal(msa_copy.index, pd.Index([1, 2]))
        assert_index_equal(msa.index, pd.Index(['foo', 'bar']))


class TestAppend(unittest.TestCase):
    def test_to_empty_msa(self):
        msa = TabularMSA([])

        msa.append(DNA('ACGT'))

        self.assertEqual(msa, TabularMSA([DNA('ACGT')]))

    def test_to_empty_with_minter(self):
        msa = TabularMSA([], minter=str)

        msa.append(DNA('ACGT'))

        self.assertEqual(msa, TabularMSA([DNA('ACGT')]))

    def test_to_empty_msa_with_index(self):
        msa = TabularMSA([])

        msa.append(DNA('ACGT'), index='a')

        self.assertEqual(
            msa,
            TabularMSA([DNA('ACGT')], index=['a']))

    def test_to_empty_msa_invalid_dtype(self):
        msa = TabularMSA([])

        with six.assertRaisesRegex(self, TypeError,
                                   'sequence.*alphabet.*Sequence'):
            msa.append(Sequence(''))

        self.assertEqual(msa, TabularMSA([]))

    def test_to_empty_msa_invalid_minter(self):
        msa = TabularMSA([])

        with self.assertRaises(KeyError):
            msa.append(DNA('ACGT'), minter='id')

        self.assertEqual(msa, TabularMSA([]))

    def test_to_non_empty_msa_invalid_minter(self):
        msa = TabularMSA([DNA('ACGT')], index=['foo'])

        with self.assertRaises(KeyError):
            msa.append(DNA('AAAA'), minter='id')

        self.assertEqual(msa, TabularMSA([DNA('ACGT')], index=['foo']))

    def test_wrong_dtype_rna(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        with six.assertRaisesRegex(self, TypeError,
                                   'matching type.*RNA.*DNA'):
            msa.append(RNA('UUUU'))

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_wrong_dtype_float(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        with six.assertRaisesRegex(self, TypeError,
                                   'matching type.*float.*DNA'):
            msa.append(42.0)

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_wrong_length(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        with six.assertRaisesRegex(
                self, ValueError,
                'must match the number of positions.*5 != 4'):
            msa.append(DNA('ACGTA'))

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_with_minter_metadata_key(self):
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

    def test_with_minter_callable(self):
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

    def test_with_index(self):
        msa = TabularMSA([DNA('AC'), DNA('GT')], index=['a', 'b'])

        msa.append(DNA('--'), index='foo')

        self.assertEqual(
            msa,
            TabularMSA([DNA('AC'), DNA('GT'), DNA('--')],
                       index=['a', 'b', 'foo']))

    def test_no_index_no_minter(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        msa.append(DNA('AAAA'))

        self.assertEqual(
            msa,
            TabularMSA([DNA('ACGT'), DNA('TGCA'), DNA('AAAA')]))

    def test_no_index_no_minter_msa_has_non_default_labels(self):
        msa = TabularMSA([DNA(''), DNA('')], index=['a', 'b'])

        with six.assertRaisesRegex(self, ValueError, "provide.*minter.*index"):
            msa.append(DNA(''))

        self.assertEqual(msa, TabularMSA([DNA(''), DNA('')], index=['a', 'b']))

    def test_with_index_type_change(self):
        msa = TabularMSA([DNA('A'), DNA('.')])

        msa.append(DNA('C'), index='foo')

        self.assertEqual(
            msa,
            TabularMSA([DNA('A'), DNA('.'), DNA('C')], index=[0, 1, 'foo']))

    def test_with_index_and_minter(self):
        msa = TabularMSA([])

        with six.assertRaisesRegex(self, ValueError, "both.*minter.*index"):
            msa.append(DNA(''), index='', minter=str)

        self.assertEqual(msa, TabularMSA([]))

    def test_multiple_appends_to_empty_msa_with_default_labels(self):
        msa = TabularMSA([])

        msa.append(RNA('U--'))
        msa.append(RNA('AA.'))

        self.assertEqual(msa, TabularMSA([RNA('U--'), RNA('AA.')]))

    def test_multiple_appends_to_non_empty_msa_with_default_labels(self):
        msa = TabularMSA([RNA('U--'), RNA('AA.')])

        msa.append(RNA('ACG'))
        msa.append(RNA('U-U'))

        self.assertEqual(
            msa,
            TabularMSA([RNA('U--'), RNA('AA.'), RNA('ACG'), RNA('U-U')]))

    def test_with_multiindex_index(self):
        msa = TabularMSA([])

        msa.append(DNA('AA'), index=('foo', 42))

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(msa.index, pd.Index([('foo', 42)]))

    def test_with_multiindex_minter(self):
        def multiindex_minter(seq):
            return ('foo', 42)

        msa = TabularMSA([])

        msa.append(DNA('AC'), minter=multiindex_minter)

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(msa.index, pd.Index([('foo', 42)]))


class TestExtend(unittest.TestCase):
    def test_empty_to_empty(self):
        msa = TabularMSA([])

        msa.extend([])

        self.assertEqual(msa, TabularMSA([]))

    def test_empty_to_non_empty(self):
        msa = TabularMSA([DNA('AC')])

        msa.extend([])

        self.assertEqual(msa, TabularMSA([DNA('AC')]))

    def test_single_sequence(self):
        msa = TabularMSA([DNA('AC')])

        msa.extend([DNA('-C')])

        self.assertEqual(msa, TabularMSA([DNA('AC'), DNA('-C')]))

    def test_multiple_sequences(self):
        msa = TabularMSA([DNA('AC')])

        msa.extend([DNA('-C'), DNA('AG')])

        self.assertEqual(msa, TabularMSA([DNA('AC'), DNA('-C'), DNA('AG')]))

    def test_from_iterable(self):
        msa = TabularMSA([])

        msa.extend(iter([DNA('ACGT'), DNA('TGCA')]))

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_from_tabular_msa_default_labels(self):
        msa = TabularMSA([DNA('AC'), DNA('TG')])

        msa.extend(TabularMSA([DNA('GG'), DNA('CC'), DNA('AA')],
                              index=['a', 'b', 'c']))

        self.assertEqual(
            msa,
            TabularMSA([DNA('AC'), DNA('TG'), DNA('GG'), DNA('CC'),
                        DNA('AA')]))

    def test_from_tabular_msa_non_default_labels(self):
        msa = TabularMSA([DNA('AC'), DNA('TG')], index=['a', 'b'])

        with six.assertRaisesRegex(self, ValueError, 'provide.*minter.*index'):
            msa.extend(TabularMSA([DNA('GG'), DNA('CC')]))

        self.assertEqual(
            msa,
            TabularMSA([DNA('AC'), DNA('TG')], index=['a', 'b']))

    def test_from_tabular_msa_with_index(self):
        msa1 = TabularMSA([DNA('AC'), DNA('TG')])
        msa2 = TabularMSA([DNA('GG'), DNA('CC'), DNA('AA')])

        msa1.extend(msa2, index=msa2.index)

        self.assertEqual(
            msa1,
            TabularMSA([DNA('AC'), DNA('TG'), DNA('GG'), DNA('CC'),
                        DNA('AA')], index=[0, 1, 0, 1, 2]))

    def test_minter_and_index(self):
        with six.assertRaisesRegex(self, ValueError, 'both.*minter.*index'):
            TabularMSA([]).extend([DNA('ACGT')], minter=str, index=['foo'])

    def test_no_minter_no_index_to_empty(self):
        msa = TabularMSA([])

        msa.extend([DNA('ACGT'), DNA('TGCA')])

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_no_minter_no_index_to_non_empty(self):
        msa = TabularMSA([DNA('ACGT')])

        msa.extend([DNA('TGCA'), DNA('--..')])

        self.assertEqual(msa,
                         TabularMSA([DNA('ACGT'), DNA('TGCA'), DNA('--..')]))

    def test_no_minter_no_index_msa_has_non_default_labels(self):
        msa = TabularMSA([DNA('ACGT')], index=[1])

        with six.assertRaisesRegex(self, ValueError, 'provide.*minter.*index'):
            msa.extend([DNA('TGCA')])

        self.assertEqual(msa, TabularMSA([DNA('ACGT')], index=[1]))

    def test_invalid_dtype(self):
        msa = TabularMSA([])

        with six.assertRaisesRegex(self, TypeError,
                                   'sequence.*alphabet.*Sequence'):
            msa.extend([Sequence('')])

        self.assertEqual(msa, TabularMSA([]))

    def test_invalid_minter(self):
        # This test (and the following error case tests) check that the MSA
        # isn't mutated when an error is raised. The "invalid" sequence is
        # preceded by valid sequence(s) to test one possible (buggy)
        # implementation of extend(): looping over sequences and calling
        # append(). These tests ensure that "valid" sequences aren't appended
        # to the MSA before the error is raised.
        msa = TabularMSA([DNA('ACGT')], index=['foo'])

        with self.assertRaises(KeyError):
            msa.extend([DNA('AAAA', metadata={'id': 'foo'}),
                        DNA('----')], minter='id')

        self.assertEqual(msa, TabularMSA([DNA('ACGT')], index=['foo']))

    def test_mismatched_dtype(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        with six.assertRaisesRegex(self, TypeError,
                                   'matching type.*RNA.*DNA'):
            msa.extend([DNA('----'), RNA('UUUU')])

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_wrong_dtype_float(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        with six.assertRaisesRegex(self, TypeError,
                                   'matching type.*float.*DNA'):
            msa.extend([DNA('GGGG'), 42.0])

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_wrong_length(self):
        msa = TabularMSA([DNA('ACGT'), DNA('TGCA')])

        with six.assertRaisesRegex(
                self, ValueError,
                'must match the number of positions.*5 != 4'):
            msa.extend([DNA('TTTT'), DNA('ACGTA')])

        self.assertEqual(msa, TabularMSA([DNA('ACGT'), DNA('TGCA')]))

    def test_sequences_index_length_mismatch(self):
        msa = TabularMSA([])

        with six.assertRaisesRegex(
                self, ValueError,
                'sequences.*2.*index length.*3'):
            msa.extend([DNA('TTTT'), DNA('ACGT')], index=['a', 'b', 'c'])

        self.assertEqual(msa, TabularMSA([]))

    def test_with_minter_metadata_key(self):
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

    def test_with_minter_callable(self):
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

    def test_with_index(self):
        msa = TabularMSA([DNA('AC'), DNA('GT')], index=['a', 'b'])

        msa.extend([DNA('--'), DNA('..')], index=['foo', 'bar'])

        self.assertEqual(
            msa,
            TabularMSA([DNA('AC'), DNA('GT'), DNA('--'), DNA('..')],
                       index=['a', 'b', 'foo', 'bar']))

    def test_with_index_type_change(self):
        msa = TabularMSA([DNA('A'), DNA('.')])

        msa.extend([DNA('C')], index=['foo'])

        self.assertEqual(
            msa,
            TabularMSA([DNA('A'), DNA('.'), DNA('C')], index=[0, 1, 'foo']))

    def test_multiple_extends_to_empty_msa_with_default_labels(self):
        msa = TabularMSA([])

        msa.extend([RNA('U-'), RNA('GG')])
        msa.extend([RNA('AA')])

        self.assertEqual(msa, TabularMSA([RNA('U-'), RNA('GG'), RNA('AA')]))

    def test_multiple_extends_to_non_empty_msa_with_default_labels(self):
        msa = TabularMSA([RNA('U--'), RNA('AA.')])

        msa.extend([RNA('ACG'), RNA('GCA')])
        msa.extend([RNA('U-U')])

        self.assertEqual(
            msa,
            TabularMSA([RNA('U--'),
                        RNA('AA.'),
                        RNA('ACG'),
                        RNA('GCA'),
                        RNA('U-U')]))

    def test_with_multiindex_index(self):
        msa = TabularMSA([])

        msa.extend([DNA('AA'), DNA('GG')], index=[('foo', 42), ('bar', 43)])

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(msa.index, pd.Index([('foo', 42), ('bar', 43)]))

    def test_with_multiindex_minter(self):
        def multiindex_minter(seq):
            if str(seq) == 'AC':
                return ('foo', 42)
            else:
                return ('bar', 43)

        msa = TabularMSA([])

        msa.extend([DNA('AC'), DNA('GG')], minter=multiindex_minter)

        self.assertIsInstance(msa.index, pd.MultiIndex)
        assert_index_equal(msa.index, pd.Index([('foo', 42), ('bar', 43)]))

    def test_with_index_object(self):
        msa = TabularMSA([])

        msa.extend([DNA('AA'), DNA('GG')],
                   index=pd.Index(['foo', 'bar']))

        self.assertEqual(
            msa,
            TabularMSA([DNA('AA'),
                        DNA('GG')], index=['foo', 'bar']))


class TestJoin(unittest.TestCase):
    def test_invalid_how(self):
        with six.assertRaisesRegex(self, ValueError, '`how`'):
            TabularMSA([]).join(TabularMSA([]), how='really')

    def test_invalid_other_type(self):
        with six.assertRaisesRegex(self, TypeError, 'TabularMSA.*DNA'):
            TabularMSA([]).join(DNA('ACGT'))

    def test_dtype_mismatch(self):
        with six.assertRaisesRegex(self, TypeError, 'dtype.*RNA.*DNA'):
            TabularMSA([DNA('AC')]).join(TabularMSA([RNA('UG')]))

        with six.assertRaisesRegex(self, TypeError, 'dtype.*None.*DNA'):
            TabularMSA([DNA('AC')]).join(TabularMSA([]))

        with six.assertRaisesRegex(self, TypeError, 'dtype.*DNA.*None'):
            TabularMSA([]).join(TabularMSA([DNA('AC')]))

    def test_duplicate_index_labels(self):
        with six.assertRaisesRegex(self, ValueError,
                                   "This MSA's index labels.*unique"):
            TabularMSA([DNA('AC'), DNA('--')], index=[0, 0]).join(
                TabularMSA([DNA('GT'), DNA('..')]))

        with six.assertRaisesRegex(self, ValueError,
                                   "`other`'s index labels.*unique"):
            TabularMSA([DNA('AC'), DNA('--')]).join(
                TabularMSA([DNA('GT'), DNA('..')], index=[0, 0]))

    def test_handles_missing_metadata_efficiently(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.')])
        msa2 = TabularMSA([DNA('-C'),
                           DNA('.G')])

        joined = msa1.join(msa2)

        self.assertEqual(
            joined,
            TabularMSA([DNA('AC-C'),
                        DNA('G..G')]))
        self.assertIsNone(msa1._metadata)
        self.assertIsNone(msa1._positional_metadata)
        self.assertIsNone(msa2._metadata)
        self.assertIsNone(msa2._positional_metadata)
        self.assertIsNone(joined._metadata)
        self.assertIsNone(joined._positional_metadata)

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

        self.assertEqual(
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

        self.assertEqual(
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

        self.assertEqual(joined, TabularMSA([]))

    def test_no_positions(self):
        msa1 = TabularMSA([DNA('', positional_metadata={'1': []}),
                           DNA('', positional_metadata={'2': []})],
                          positional_metadata={'foo': []})
        msa2 = TabularMSA([DNA('', positional_metadata={'3': []}),
                           DNA('', positional_metadata={'4': []})],
                          positional_metadata={'foo': []})

        joined = msa1.join(msa2)

        self.assertEqual(
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

        self.assertEqual(
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

        self.assertEqual(
            joined,
            TabularMSA([DNA('ACCA'),
                        DNA('G..G'),
                        DNA('C--C')],
                       positional_metadata={'foo': [1, 2, 3, 4],
                                            'bar': ['a', 'b', 'c', 'd']}))

    def test_how_strict_failure_index_mismatch(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.'),
                           DNA('C-')])
        msa2 = TabularMSA([DNA('-C'),
                           DNA('.G'),
                           DNA('CA'),
                           DNA('--')])

        with six.assertRaisesRegex(self, ValueError,
                                   'Index labels must all match'):
            msa1.join(msa2)

    def test_how_strict_failure_positional_metadata_mismatch(self):
        msa1 = TabularMSA([DNA('AC'),
                           DNA('G.')],
                          positional_metadata={'foo': [1, 2],
                                               'bar': ['a', 'b']})
        msa2 = TabularMSA([DNA('-C'),
                           DNA('.G')],
                          positional_metadata={'foo': [3, 4]})

        with six.assertRaisesRegex(self, ValueError,
                                   'Positional metadata columns.*match'):
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

        self.assertEqual(
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

        self.assertEqual(
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

        self.assertEqual(joined, TabularMSA([]))

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

        self.assertEqual(
            joined,
            TabularMSA([DNA('--...'),
                        DNA('ACCAA'),
                        DNA('G..GG'),
                        DNA('C--CC'),
                        DNA('-----')], index=range(-1, 4),
                       positional_metadata={
                           'foo': [1, 2, 3, 4, 5],
                           'bar': ['a', 'b', np.nan, np.nan, np.nan],
                           'baz': [np.nan, np.nan, 'c', 'd', 'e']}))

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

        self.assertEqual(
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

        self.assertEqual(
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

    def test_no_sequences_reverse(self):
        msa = TabularMSA([])

        obs = list(msa.iter_positions(reverse=True))

        self.assertEqual(obs, [])

    def test_no_positions(self):
        msa = TabularMSA([DNA(''),
                          DNA('')])

        obs = list(msa.iter_positions())

        self.assertEqual(obs, [])

    def test_no_positions_reverse(self):
        msa = TabularMSA([DNA(''),
                          DNA('')])

        obs = list(msa.iter_positions(reverse=True))

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
                      positional_metadata={'foo': [42, np.nan, -1],
                                           'bar': [np.nan, np.nan, 'baz']}),
             Sequence('C--', metadata={'pm1': 1.5, 'foo': 99},
                      positional_metadata={'foo': [43, np.nan, -2],
                                           'bar': [np.nan, np.nan, 'bazz']})])

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
                      positional_metadata={'foo': [43, np.nan, -2],
                                           'bar': [np.nan, np.nan, 'bazz']}),
             Sequence('AA-', metadata={'pm1': 0.5, 'foo': 9},
                      positional_metadata={'foo': [42, np.nan, -1],
                                           'bar': [np.nan, np.nan, 'baz']})])

    def test_handles_missing_positional_metadata_efficiently(self):
        msa = TabularMSA([DNA('AC'),
                          DNA('A-')])

        self.assertIsNone(msa._positional_metadata)

        list(msa.iter_positions())

        self.assertIsNone(msa._positional_metadata)


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

    def test_handles_missing_positional_metadata_efficiently(self):
        msa = TabularMSA([DNA('AC'),
                          DNA('AC')])

        self.assertIsNone(msa._positional_metadata)

        cons = msa.consensus()

        self.assertIsNone(msa._positional_metadata)
        self.assertIsNone(cons._positional_metadata)

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


class TestGapFrequencies(unittest.TestCase):
    def test_default_behavior(self):
        msa = TabularMSA([DNA('AA.'),
                          DNA('-A-')])

        freqs = msa.gap_frequencies()

        npt.assert_array_equal(np.array([1, 0, 2]), freqs)

    def test_invalid_axis_str(self):
        with six.assertRaisesRegex(self, ValueError, "axis.*'foo'"):
            TabularMSA([]).gap_frequencies(axis='foo')

    def test_invalid_axis_int(self):
        with six.assertRaisesRegex(self, ValueError, "axis.*2"):
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
        class CustomSequence(IUPACSequence):
            @classproperty
            @overrides(IUPACSequence)
            def gap_chars(cls):
                return set('0123456789')

            @classproperty
            @overrides(IUPACSequence)
            def nondegenerate_chars(cls):
                return set('')

            @classproperty
            @overrides(IUPACSequence)
            def degenerate_map(cls):
                return {}

        msa = TabularMSA([CustomSequence('0123456789')])

        freqs = msa.gap_frequencies(axis='position', relative=True)

        npt.assert_array_equal(np.array([1.0]), freqs)

    def test_custom_gap_characters(self):
        class CustomSequence(IUPACSequence):
            @classproperty
            @overrides(IUPACSequence)
            def gap_chars(cls):
                return set('#$*')

            @classproperty
            @overrides(IUPACSequence)
            def nondegenerate_chars(cls):
                return set('ABC-.')

            @classproperty
            @overrides(IUPACSequence)
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

        position = msa._get_position(1)

        self.assertEqual(position, Sequence('C-'))

    def test_with_positional_metadata(self):
        msa = TabularMSA([DNA('ACG'),
                          DNA('A-G')],
                         positional_metadata={'foo': [42, 43, 44],
                                              'bar': ['abc', 'def', 'ghi']})

        position = msa._get_position(1)

        self.assertEqual(position,
                         Sequence('C-', metadata={'foo': 43, 'bar': 'def'}))

    def test_handles_positional_metadata_efficiently(self):
        msa = TabularMSA([DNA('AA'),
                          DNA('--')])

        msa._get_position(1)

        self.assertIsNone(msa._positional_metadata)


class TestIsSequenceAxis(unittest.TestCase):
    def setUp(self):
        self.msa = TabularMSA([])

    def test_invalid_str(self):
        with six.assertRaisesRegex(self, ValueError, "axis.*'foo'"):
            self.msa._is_sequence_axis('foo')

    def test_invalid_int(self):
        with six.assertRaisesRegex(self, ValueError, "axis.*2"):
            self.msa._is_sequence_axis(2)

    def test_positive_str(self):
        self.assertTrue(self.msa._is_sequence_axis('sequence'))

    def test_positive_int(self):
        self.assertTrue(self.msa._is_sequence_axis(0))

    def test_negative_str(self):
        self.assertFalse(self.msa._is_sequence_axis('position'))

    def test_negative_int(self):
        self.assertFalse(self.msa._is_sequence_axis(1))


class TestRepr(unittest.TestCase):
    def test_repr(self):
        # basic sanity checks -- more extensive testing of formatting and
        # special cases is performed in TabularMSAReprDoctests below. here we
        # only test that pieces of the repr are present. these tests also
        # exercise coverage for py2/3 since the doctests in
        # TabularMSAReprDoctests only currently run in py3.

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
        self.assertIn('<DNA>', obs)

        # no metadata
        obs = repr(TabularMSA([DNA('ACGT')]))
        self.assertEqual(obs.count('\n'), 6)
        self.assertTrue(obs.startswith('TabularMSA'))
        self.assertIn('sequence count: 1', obs)
        self.assertIn('position count: 4', obs)
        self.assertIn('<DNA>', obs)
        self.assertTrue(obs.endswith('ACGT'))

        # sequence spanning > 5 lines
        obs = repr(TabularMSA([DNA('A' * 71) for x in range(6)]))
        self.assertEqual(obs.count('\n'), 10)
        self.assertTrue(obs.startswith('TabularMSA'))
        self.assertIn('sequence count: 6', obs)
        self.assertIn('position count: 71', obs)
        self.assertIn('\n...\n', obs)
        self.assertIn('<DNA>', obs)
        self.assertTrue(obs.endswith('AAAA'))

        # sequences overflowing
        obs = repr(TabularMSA([DNA('A' * 72)]))
        self.assertEqual(obs.count('\n'), 6)
        self.assertTrue(obs.startswith('TabularMSA'))
        self.assertIn('sequence count: 1', obs)
        self.assertIn('position count: 72', obs)
        self.assertIn('<DNA>', obs)
        self.assertTrue(obs.endswith(' ... ' + 'A'*33))


# NOTE: this must be a *separate* class for doctests only (no unit tests). nose
# will not run the unit tests otherwise
#
# these doctests exercise the correct formatting of TabularMSA's repr in a
# variety of situations. they are more extensive than the unit tests above
# (TestRepr.test_repr) but are only currently run in py3. thus, they cannot
# be relied upon for coverage (the unit tests take care of this)
class TabularMSAReprDoctests(object):
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
    TabularMSA<DNA>
    ---------------------
    Stats:
        sequence count: 1
        position count: 0
    ---------------------

    MSA with single sequence with single character:

    >>> TabularMSA([DNA('G')])
    TabularMSA<DNA>
    ---------------------
    Stats:
        sequence count: 1
        position count: 1
    ---------------------
    G

    MSA with multicharacter sequence:

    >>> TabularMSA([DNA('ACGT')])
    TabularMSA<DNA>
    ---------------------
    Stats:
        sequence count: 1
        position count: 4
    ---------------------
    ACGT

    Full single line:

    >>> TabularMSA([DNA('A' * 71)])
    TabularMSA<DNA>
    -----------------------------------------------------------------------
    Stats:
        sequence count: 1
        position count: 71
    -----------------------------------------------------------------------
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    Full single line with 1 character overflow:

    >>> TabularMSA([DNA('A' * 72)])
    TabularMSA<DNA>
    -----------------------------------------------------------------------
    Stats:
        sequence count: 1
        position count: 72
    -----------------------------------------------------------------------
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA ... AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    Two sequences with full lines:

    >>> TabularMSA([DNA('T' * 71), DNA('T' * 71)])
    TabularMSA<DNA>
    -----------------------------------------------------------------------
    Stats:
        sequence count: 2
        position count: 71
    -----------------------------------------------------------------------
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

    Two sequences with full lines with 1 character overflow:

    >>> TabularMSA([DNA('T' * 72), DNA('T' * 72)])
    TabularMSA<DNA>
    -----------------------------------------------------------------------
    Stats:
        sequence count: 2
        position count: 72
    -----------------------------------------------------------------------
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT ... TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
    TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT ... TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

    Five full lines (maximum amount of information):

    >>> TabularMSA([DNA('A' * 71) for x in range(5)])
    TabularMSA<DNA>
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
    TabularMSA<DNA>
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
    >>> positional_metadata = pd.DataFrame.from_items([
    ...     # str key, int list value
    ...     ('foo', [1, 2, 3, 4]),
    ...     # float key, float list value
    ...     (42.5, [2.5, 3.0, 4.2, -0.00001]),
    ...     # int key, object list value
    ...     (42, [[], 4, 5, {}]),
    ...     # truncated key (too long), bool list value
    ...     ('abc' * 90, [True, False, False, True]),
    ...     # None key
    ...     (None, range(4))])
    >>> TabularMSA([DNA('ACGT')], metadata=metadata,
    ...            positional_metadata=positional_metadata)
    TabularMSA<DNA>
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
