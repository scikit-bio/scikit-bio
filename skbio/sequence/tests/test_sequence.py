# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

import copy
import functools
import itertools
import re
from types import GeneratorType
from collections import Hashable
from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd
import scipy.spatial.distance

import skbio.sequence.distance
from skbio import Sequence, DNA
from skbio.util import assert_data_frame_almost_equal
from skbio.sequence._sequence import (_single_index_to_slice, _is_single_index,
                                      _as_slice_if_single_index)
from skbio.util._testing import ReallyEqualMixin
from skbio.metadata._testing import (MetadataMixinTests,
                                     IntervalMetadataMixinTests,
                                     PositionalMetadataMixinTests)
from skbio.metadata import IntervalMetadata


class SequenceSubclass(Sequence):
    """Used for testing purposes."""
    pass


class SequenceSubclassTwo(Sequence):
    """Used for testing purposes."""
    pass


class TestSequenceMetadata(TestCase, ReallyEqualMixin, MetadataMixinTests):
    def setUp(self):
        self._metadata_constructor_ = functools.partial(Sequence, '')


class TestSequencePositionalMetadata(TestCase, ReallyEqualMixin,
                                     PositionalMetadataMixinTests):
    def setUp(self):
        def factory(axis_len, positional_metadata=None):
            return Sequence('Z' * axis_len,
                            positional_metadata=positional_metadata)
        self._positional_metadata_constructor_ = factory


class TestSequenceIntervalMetadata(TestCase, ReallyEqualMixin,
                                   IntervalMetadataMixinTests):
    def setUp(self):
        super()._set_up()

        def factory(axis_len, interval_metadata=None):
            return Sequence('Z' * axis_len,
                            interval_metadata=interval_metadata)
        self._interval_metadata_constructor_ = factory


class TestSequenceBase(TestCase):
    def setUp(self):
        self.sequence_kinds = frozenset([
            str, Sequence,
            lambda s: np.frombuffer(s.encode('ascii'), dtype='|S1'),
            lambda s: np.frombuffer(s.encode('ascii'), dtype=np.uint8)])


class TestSequence(TestSequenceBase, ReallyEqualMixin):
    def setUp(self):
        super(TestSequence, self).setUp()

        self.lowercase_seq = Sequence('AAAAaaaa', lowercase='key')

        def empty_generator():
            yield from ()

        self.getitem_empty_indices = [
            [],
            (),
            {},
            empty_generator(),
            # ndarray of implicit float dtype
            np.array([]),
            np.array([], dtype=int)]

    def test_concat_bad_how(self):
        seq1 = seq2 = Sequence("123")
        with self.assertRaises(ValueError):
            Sequence.concat([seq1, seq2], how='foo')

    def test_concat_on_subclass(self):
        seq1 = SequenceSubclass("123")
        seq2 = Sequence("123")
        result = SequenceSubclass.concat([seq1, seq2])
        self.assertIs(type(result), SequenceSubclass)
        self.assertEqual(result, SequenceSubclass("123123"))

    def test_concat_on_empty_iterator(self):
        result = SequenceSubclass.concat((_ for _ in []))
        self.assertIs(type(result), SequenceSubclass)
        self.assertEqual(result, SequenceSubclass(""))

    def test_concat_on_bad_subclass(self):
        seq1 = Sequence("123")
        seq2 = SequenceSubclassTwo("123")
        with self.assertRaises(TypeError):
            SequenceSubclass.concat([seq1, seq2])

    def test_concat_interval_metadata(self):
        seq1 = Sequence("1234")
        seq1.interval_metadata.add(
            [(0, 2)], [(True, False)], {'gene': 'sagA'})
        seq2 = Sequence("5678")
        seq2.interval_metadata.add(
            [(1, 3)], [(False, True)], {'gene': 'sagB'})
        obs = Sequence.concat([seq1, seq2])
        exp = Sequence('12345678')
        exp.interval_metadata.add(
            [(0, 2)], [(True, False)], {'gene': 'sagA'})
        exp.interval_metadata.add(
            [(5, 7)], [(False, True)], {'gene': 'sagB'})
        self.assertEqual(exp, obs)

    def test_concat_one_seq_has_none_interval_metadata(self):
        seq1 = Sequence("1234")
        seq1.interval_metadata.add(
            [(0, 2)], [(True, False)], {'gene': 'sagA'})
        seq2 = Sequence("5678")
        seq3 = Sequence("910")
        seq3.interval_metadata.add(
            [(1, 3)], [(False, True)], {'gene': 'sagB'})
        obs = Sequence.concat([seq1, seq2, seq3])
        exp = Sequence('12345678910')
        exp.interval_metadata.add(
            [(0, 2)], [(True, False)], {'gene': 'sagA'})
        exp.interval_metadata.add(
            [(9, 11)], [(False, True)], {'gene': 'sagB'})
        self.assertEqual(exp, obs)

    def test_concat_default_how(self):
        seq1 = Sequence("1234", positional_metadata={'a': [1]*4})
        seq2 = Sequence("5678", positional_metadata={'a': [2]*4})
        seqbad = Sequence("9", positional_metadata={'b': [9]})
        result1 = Sequence.concat([seq1, seq2])
        result2 = Sequence.concat([seq1, seq2], how='strict')
        self.assertEqual(result1, result2)
        with self.assertRaisesRegex(ValueError,
                                    r'.*positional.*metadata.*inner.*outer.*'):
            Sequence.concat([seq1, seq2, seqbad])

    def test_concat_strict_simple(self):
        expected = Sequence(
            "12345678", positional_metadata={'a': [1, 1, 1, 1, 2, 2, 2, 2]})
        seq1 = Sequence("1234", positional_metadata={'a': [1]*4})
        seq2 = Sequence("5678", positional_metadata={'a': [2]*4})
        result = Sequence.concat([seq1, seq2], how='strict')
        self.assertEqual(result, expected)
        self.assertFalse(result.metadata)

    def test_concat_strict_many(self):
        odd_key = frozenset()
        expected = Sequence("13579",
                            positional_metadata={'a': list('skbio'),
                                                 odd_key: [1, 2, 3, 4, 5]})
        result = Sequence.concat([
                Sequence("1", positional_metadata={'a': ['s'], odd_key: [1]}),
                Sequence("3", positional_metadata={'a': ['k'], odd_key: [2]}),
                Sequence("5", positional_metadata={'a': ['b'], odd_key: [3]}),
                Sequence("7", positional_metadata={'a': ['i'], odd_key: [4]}),
                Sequence("9", positional_metadata={'a': ['o'], odd_key: [5]})
            ], how='strict')
        self.assertEqual(result, expected)
        self.assertFalse(result.metadata)

    def test_concat_strict_fail(self):
        seq1 = Sequence("1", positional_metadata={'a': [1]})
        seq2 = Sequence("2", positional_metadata={'b': [2]})
        with self.assertRaisesRegex(ValueError,
                                    r'.*positional.*metadata.*inner.*outer.*'):
            Sequence.concat([seq1, seq2], how='strict')

    def test_concat_outer_simple(self):
        seq1 = Sequence("1234")
        seq2 = Sequence("5678")
        result = Sequence.concat([seq1, seq2], how='outer')
        self.assertEqual(result, Sequence("12345678"))
        self.assertFalse(result.metadata)

    def test_concat_outer_missing(self):
        a = {}
        b = {}
        seq1 = Sequence("12", positional_metadata={'a': ['1', '2']})
        seq2 = Sequence("34", positional_metadata={'b': [3, 4], 'c': [a, b]})
        seq3 = Sequence("56")
        seq4 = Sequence("78", positional_metadata={'a': [7, 8]})
        seq5 = Sequence("90", positional_metadata={'b': [9, 0]})

        result = Sequence.concat([seq1, seq2, seq3, seq4, seq5], how='outer')
        expected = Sequence("1234567890", positional_metadata={
                                'a': ['1', '2', np.nan, np.nan, np.nan, np.nan,
                                      7, 8, np.nan, np.nan],
                                'b': [np.nan, np.nan, 3, 4, np.nan, np.nan,
                                      np.nan, np.nan, 9, 0],
                                'c': [np.nan, np.nan, a, b, np.nan, np.nan,
                                      np.nan, np.nan, np.nan, np.nan]
                            })
        self.assertEqual(result, expected)
        self.assertFalse(result.metadata)

    def test_concat_inner_simple(self):
        seq1 = Sequence("1234")
        seq2 = Sequence("5678", positional_metadata={'discarded': [1] * 4})
        result = Sequence.concat([seq1, seq2], how='inner')
        self.assertEqual(result, Sequence("12345678"))
        self.assertFalse(result.metadata)

    def test_concat_inner_missing(self):
        seq1 = Sequence("12", positional_metadata={'a': ['1', '2'],
                                                   'c': [{}, {}]})
        seq2 = Sequence("34", positional_metadata={'a': [3, 4], 'b': [3, 4]})
        seq3 = Sequence("56", positional_metadata={'a': [5, 6], 'b': [5, 6]})

        result = Sequence.concat([seq1, seq2, seq3], how='inner')
        expected = Sequence("123456", positional_metadata={'a': ['1', '2', 3,
                                                                 4, 5, 6]})
        self.assertEqual(result, expected)
        self.assertFalse(result.metadata)

    def test_init_default_parameters(self):
        seq = Sequence('.ABC123xyz-')

        npt.assert_equal(seq.values, np.array('.ABC123xyz-', dtype='c'))
        self.assertEqual('.ABC123xyz-', str(seq))
        self.assertFalse(seq.metadata)
        self.assertEqual(seq.metadata, {})
        assert_data_frame_almost_equal(seq.positional_metadata,
                                       pd.DataFrame(index=range(11)))
        self.assertEqual(seq.interval_metadata,
                         IntervalMetadata(len(seq)))

    def test_init_nondefault_parameters(self):
        s = '.ABC123xyz-'
        im = IntervalMetadata(len(s))
        im.add([(0, 1)], metadata={'gene': 'sagA'})
        seq = Sequence(s,
                       metadata={'id': 'foo', 'description': 'bar baz'},
                       positional_metadata={'quality': range(11)},
                       interval_metadata=im)

        self.assertEqual(seq.interval_metadata, im)

        npt.assert_equal(seq.values, np.array('.ABC123xyz-', dtype='c'))
        self.assertEqual(s, str(seq))

        self.assertTrue(seq.metadata)
        self.assertEqual(seq.metadata, {'id': 'foo', 'description': 'bar baz'})

        assert_data_frame_almost_equal(
            seq.positional_metadata,
            pd.DataFrame({'quality': range(11)}, index=range(11)))

    def test_init_empty_sequence(self):
        # Test constructing an empty sequence using each supported input type.
        for s in (b'',  # bytes
                  '',  # unicode
                  np.array('', dtype='c'),  # char vector
                  np.frombuffer(b'', dtype=np.uint8),  # byte vec
                  Sequence('')):  # another Sequence object
            seq = Sequence(s)

            self.assertIsInstance(seq.values, np.ndarray)
            self.assertEqual(seq.values.dtype, '|S1')
            self.assertEqual(seq.values.shape, (0, ))
            npt.assert_equal(seq.values, np.array('', dtype='c'))
            self.assertEqual(str(seq), '')
            self.assertEqual(len(seq), 0)

            self.assertFalse(seq.metadata)
            self.assertEqual(seq.metadata, {})

            self.assertEqual(seq.interval_metadata,
                             IntervalMetadata(0))

            assert_data_frame_almost_equal(seq.positional_metadata,
                                           pd.DataFrame(index=range(0)))

    def test_init_single_character_sequence(self):
        for s in (b'A',
                  'A',
                  np.array('A', dtype='c'),
                  np.frombuffer(b'A', dtype=np.uint8),
                  Sequence('A')):
            seq = Sequence(s)

            self.assertIsInstance(seq.values, np.ndarray)
            self.assertEqual(seq.values.dtype, '|S1')
            self.assertEqual(seq.values.shape, (1,))
            npt.assert_equal(seq.values, np.array('A', dtype='c'))
            self.assertEqual(str(seq), 'A')
            self.assertEqual(len(seq), 1)

            self.assertFalse(seq.metadata)
            self.assertEqual(seq.metadata, {})

            self.assertEqual(seq.interval_metadata,
                             IntervalMetadata(1))
            assert_data_frame_almost_equal(seq.positional_metadata,
                                           pd.DataFrame(index=range(1)))

    def test_init_multiple_character_sequence(self):
        for s in (b'.ABC\t123  xyz-',
                  '.ABC\t123  xyz-',
                  np.array('.ABC\t123  xyz-', dtype='c'),
                  np.frombuffer(b'.ABC\t123  xyz-', dtype=np.uint8),
                  Sequence('.ABC\t123  xyz-')):
            seq = Sequence(s)

            self.assertIsInstance(seq.values, np.ndarray)
            self.assertEqual(seq.values.dtype, '|S1')
            self.assertEqual(seq.values.shape, (14,))
            npt.assert_equal(seq.values,
                             np.array('.ABC\t123  xyz-', dtype='c'))
            self.assertEqual(str(seq), '.ABC\t123  xyz-')
            self.assertEqual(len(seq), 14)

            self.assertFalse(seq.metadata)
            self.assertEqual(seq.metadata, {})

            self.assertEqual(seq.interval_metadata,
                             IntervalMetadata(14))
            assert_data_frame_almost_equal(seq.positional_metadata,
                                           pd.DataFrame(index=range(14)))

    def test_init_from_sequence_object(self):
        # We're testing this in its simplest form in other tests. This test
        # exercises more complicated cases of building a sequence from another
        # sequence.

        # just the sequence, no other metadata
        seq = Sequence('ACGT')
        self.assertEqual(Sequence(seq), seq)

        # sequence with metadata should have everything propagated
        seq = Sequence('ACGT',
                       metadata={'id': 'foo', 'description': 'bar baz'},
                       positional_metadata={'quality': range(4)})
        seq.interval_metadata.add([(0, 1)], metadata={'gene': 'sagA'})
        self.assertEqual(Sequence(seq), seq)

        # should be able to override metadata
        im = IntervalMetadata(4)
        im.add([(0, 2)], metadata={'gene': 'sagB'})
        self.assertEqual(
            Sequence(seq, metadata={'id': 'abc', 'description': '123'},
                     positional_metadata={'quality': [42] * 4},
                     interval_metadata=im),
            Sequence('ACGT', metadata={'id': 'abc', 'description': '123'},
                     positional_metadata={'quality': [42] * 4},
                     interval_metadata=im))

        # subclasses work too
        im = IntervalMetadata(4)
        im.add([(0, 2)], metadata={'gene': 'sagB'})
        seq = SequenceSubclass('ACGT',
                               metadata={'id': 'foo',
                                         'description': 'bar baz'},
                               positional_metadata={'quality': range(4)},
                               interval_metadata=im)

        self.assertEqual(
            Sequence(seq),
            Sequence('ACGT', metadata={'id': 'foo', 'description': 'bar baz'},
                     positional_metadata={'quality': range(4)},
                     interval_metadata=im))

    def test_init_from_non_descendant_sequence_object(self):
        seq = SequenceSubclass('ACGT')
        with self.assertRaises(TypeError) as cm:
            SequenceSubclassTwo(seq)

        error = str(cm.exception)
        self.assertIn("SequenceSubclass", error)
        self.assertIn("SequenceSubclassTwo", error)
        self.assertIn("cast", error)

    def test_init_from_contiguous_sequence_bytes_view(self):
        bytes = np.array([65, 42, 66, 42, 65], dtype=np.uint8)
        view = bytes[:3]
        seq = Sequence(view)

        # sequence should be what we'd expect
        self.assertEqual(seq, Sequence('A*B'))

        # we shouldn't own the memory because no copy should have been made
        self.assertFalse(seq._owns_bytes)

        # can't mutate view because it isn't writeable anymore
        with self.assertRaises(ValueError):
            view[1] = 100

        # sequence shouldn't have changed
        self.assertEqual(seq, Sequence('A*B'))

        # mutate bytes (*not* the view)
        bytes[0] = 99

        # Sequence changed because we are only able to make the view read-only,
        # not its source (bytes). This is somewhat inconsistent behavior that
        # is (to the best of our knowledge) outside our control.
        self.assertEqual(seq, Sequence('c*B'))

    def test_init_from_noncontiguous_sequence_bytes_view(self):
        bytes = np.array([65, 42, 66, 42, 65], dtype=np.uint8)
        view = bytes[::2]
        seq = Sequence(view)

        # sequence should be what we'd expect
        self.assertEqual(seq, Sequence('ABA'))

        # we should own the memory because a copy should have been made
        self.assertTrue(seq._owns_bytes)

        # mutate bytes and its view
        bytes[0] = 99
        view[1] = 100

        # sequence shouldn't have changed
        self.assertEqual(seq, Sequence('ABA'))

    def test_init_no_copy_of_sequence(self):
        bytes = np.array([65, 66, 65], dtype=np.uint8)
        seq = Sequence(bytes)

        # should share the same memory
        self.assertIs(seq._bytes, bytes)

        # shouldn't be able to mutate the Sequence object's internals by
        # mutating the shared memory
        with self.assertRaises(ValueError):
            bytes[1] = 42

    def test_init_invalid_sequence(self):
        # invalid dtype (numpy.ndarray input)
        with self.assertRaises(TypeError):
            # int64
            Sequence(np.array([1, 2, 3]))
        with self.assertRaises(TypeError):
            # |S21
            Sequence(np.array([1, "23", 3]))
        with self.assertRaises(TypeError):
            # object
            Sequence(np.array([1, {}, ()]))

        # invalid input type (non-numpy.ndarray input)
        with self.assertRaisesRegex(TypeError, r'tuple'):
            Sequence(('a', 'b', 'c'))
        with self.assertRaisesRegex(TypeError, r'list'):
            Sequence(['a', 'b', 'c'])
        with self.assertRaisesRegex(TypeError, r'set'):
            Sequence({'a', 'b', 'c'})
        with self.assertRaisesRegex(TypeError, r'dict'):
            Sequence({'a': 42, 'b': 43, 'c': 44})
        with self.assertRaisesRegex(TypeError, r'int'):
            Sequence(42)
        with self.assertRaisesRegex(TypeError, r'float'):
            Sequence(4.2)
        with self.assertRaisesRegex(TypeError, r'int64'):
            Sequence(np.int_(50))
        with self.assertRaisesRegex(TypeError, r'float64'):
            Sequence(np.float_(50))
        with self.assertRaisesRegex(TypeError, r'Foo'):
            class Foo:
                pass
            Sequence(Foo())

        # out of ASCII range
        with self.assertRaises(UnicodeEncodeError):
            Sequence('abc\u1F30')

    def test_values_property(self):
        # Property tests are only concerned with testing the interface
        # provided by the property: that it can be accessed, can't be
        # reassigned or mutated in place, and that the correct type is
        # returned. More extensive testing of border cases (e.g., different
        # sequence lengths or input types, odd characters, etc.) are performed
        # in Sequence.__init__ tests.

        seq = Sequence('ACGT')

        # should get back a numpy.ndarray of '|S1' dtype
        self.assertIsInstance(seq.values, np.ndarray)
        self.assertEqual(seq.values.dtype, '|S1')
        npt.assert_equal(seq.values, np.array('ACGT', dtype='c'))

        # test that we can't mutate the property
        with self.assertRaises(ValueError):
            seq.values[1] = 'A'

        # test that we can't set the property
        with self.assertRaises(AttributeError):
            seq.values = np.array("GGGG", dtype='c')

    def test_sequence_numpy_compatibility(self):
        seq = Sequence('abc123')

        array = np.asarray(seq)

        self.assertIsInstance(array, np.ndarray)
        self.assertEqual(array.dtype, '|S1')
        npt.assert_equal(array, np.array('abc123', dtype='c'))
        npt.assert_equal(array, seq.values)

        with self.assertRaises(ValueError):
            array[1] = 'B'

    def test_observed_chars_property(self):
        self.assertEqual(Sequence('').observed_chars, set())
        self.assertEqual(Sequence('x').observed_chars, {'x'})
        self.assertEqual(Sequence('xYz').observed_chars, {'x', 'Y', 'z'})
        self.assertEqual(Sequence('zzz').observed_chars, {'z'})
        self.assertEqual(Sequence('xYzxxZz').observed_chars,
                         {'x', 'Y', 'z', 'Z'})
        self.assertEqual(Sequence('\t   ').observed_chars, {' ', '\t'})

        im = IntervalMetadata(6)
        im.add([(0, 2)], metadata={'gene': 'sagB'})
        self.assertEqual(
            Sequence('aabbcc', metadata={'foo': 'bar'},
                     positional_metadata={'foo': range(6)},
                     interval_metadata=im).observed_chars,
            {'a', 'b', 'c'})

        with self.assertRaises(AttributeError):
            Sequence('ACGT').observed_chars = {'a', 'b', 'c'}

    def test_eq_and_ne(self):
        seq_a = Sequence("A")
        seq_b = Sequence("B")

        im = IntervalMetadata(1)
        im.add([(0, 1)], metadata={'gene': 'sagA'})
        im2 = IntervalMetadata(1)
        im.add([(0, 1)], metadata={'gene': 'sagB'})

        self.assertTrue(seq_a == seq_a)
        self.assertTrue(Sequence("a") == Sequence("a"))
        self.assertTrue(Sequence("a", metadata={'id': 'b'}) ==
                        Sequence("a", metadata={'id': 'b'}))
        self.assertTrue(Sequence("a",
                                 metadata={'id': 'b', 'description': 'c'}) ==
                        Sequence("a",
                                 metadata={'id': 'b', 'description': 'c'}))
        self.assertTrue(Sequence("a", metadata={'id': 'b', 'description': 'c'},
                                 positional_metadata={'quality': [1]},
                                 interval_metadata=im) ==
                        Sequence("a", metadata={'id': 'b', 'description': 'c'},
                                 positional_metadata={'quality': [1]},
                                 interval_metadata=im))

        self.assertTrue(seq_a != seq_b)
        self.assertTrue(SequenceSubclass("a") != Sequence("a"))
        self.assertTrue(Sequence("a") != Sequence("b"))
        self.assertTrue(Sequence("a") != Sequence("a", metadata={'id': 'b'}))
        self.assertTrue(Sequence("a", metadata={'id': 'c'}) !=
                        Sequence("a",
                                 metadata={'id': 'c', 'description': 't'}))
        self.assertTrue(Sequence("a", positional_metadata={'quality': [1]}) !=
                        Sequence("a"))
        self.assertTrue(Sequence("a", interval_metadata=im) !=
                        Sequence("a"))

        self.assertTrue(Sequence("a", positional_metadata={'quality': [1]}) !=
                        Sequence("a", positional_metadata={'quality': [2]}))
        self.assertTrue(Sequence("a", interval_metadata=im) !=
                        Sequence("a", interval_metadata=im2))

        self.assertTrue(Sequence("c", positional_metadata={'quality': [3]}) !=
                        Sequence("b", positional_metadata={'quality': [3]}))
        self.assertTrue(Sequence("c", interval_metadata=im) !=
                        Sequence("b", interval_metadata=im))
        self.assertTrue(Sequence("a", metadata={'id': 'b'}) !=
                        Sequence("c", metadata={'id': 'b'}))

    def test_eq_sequences_without_metadata_compare_equal(self):
        self.assertTrue(Sequence('') == Sequence(''))
        self.assertTrue(Sequence('z') == Sequence('z'))
        self.assertTrue(
            Sequence('ACGT') == Sequence('ACGT'))

    def test_eq_sequences_with_metadata_compare_equal(self):
        seq1 = Sequence('ACGT', metadata={'id': 'foo', 'desc': 'abc'},
                        positional_metadata={'qual': [1, 2, 3, 4]})
        seq2 = Sequence('ACGT', metadata={'id': 'foo', 'desc': 'abc'},
                        positional_metadata={'qual': [1, 2, 3, 4]})
        self.assertTrue(seq1 == seq2)

        # order shouldn't matter
        self.assertTrue(seq2 == seq1)

    def test_eq_sequences_from_different_sources_compare_equal(self):
        # sequences that have the same data but are constructed from different
        # types of data should compare equal
        im = IntervalMetadata(4)
        im.add([(0, 2)], metadata={'gene': 'sagB'})
        seq1 = Sequence('ACGT', metadata={'id': 'foo', 'desc': 'abc'},
                        positional_metadata={'quality': (1, 2, 3, 4)},
                        interval_metadata=im)
        seq2 = Sequence(np.array([65, 67, 71, 84], dtype=np.uint8),
                        metadata={'id': 'foo', 'desc': 'abc'},
                        positional_metadata={'quality': np.array([1, 2, 3,
                                                                  4])},
                        interval_metadata=im)
        self.assertTrue(seq1 == seq2)

    def test_eq_type_mismatch(self):
        seq1 = Sequence('ACGT')
        seq2 = SequenceSubclass('ACGT')
        self.assertFalse(seq1 == seq2)

    def test_eq_positional_metadata_mismatch(self):
        # both provided
        seq1 = Sequence('ACGT', positional_metadata={'quality': [1, 2, 3, 4]})
        seq2 = Sequence('ACGT', positional_metadata={'quality': [1, 2, 3, 5]})
        self.assertFalse(seq1 == seq2)

        # one provided
        seq1 = Sequence('ACGT', positional_metadata={'quality': [1, 2, 3, 4]})
        seq2 = Sequence('ACGT')
        self.assertFalse(seq1 == seq2)

    def test_eq_interval_metadata_mismatch(self):
        im1 = IntervalMetadata(4)
        im1.add([(0, 3)], metadata={'gene': 'sagA'})
        im2 = IntervalMetadata(4)
        im2.add([(0, 2)], metadata={'gene': 'sagA'})
        # both provided
        seq1 = Sequence('ACGT', interval_metadata=im1)
        seq2 = Sequence('ACGT', interval_metadata=im2)
        self.assertFalse(seq1 == seq2)

        # one provided
        seq1 = Sequence('ACGT', interval_metadata=im1)
        seq2 = Sequence('ACGT')
        self.assertFalse(seq1 == seq2)

    def test_eq_sequence_mismatch(self):
        seq1 = Sequence('ACGT')
        seq2 = Sequence('TGCA')
        self.assertFalse(seq1 == seq2)

    def test_getitem_gives_new_sequence(self):
        seq = Sequence("Sequence string !1@2#3?.,")
        self.assertFalse(seq is seq[:])

    def test_getitem_drops_interval_metadata(self):
        s = "Sequence string !1@2#3?.,"
        seq = Sequence(s, metadata={'id': 'id', 'description': 'dsc'})
        seq.interval_metadata.add([(0, 3)], metadata={'gene': 'sagA'})

        eseq = Sequence('Se', metadata={'id': 'id', 'description': 'dsc'})
        self.assertEqual(seq[:2], eseq)

    def test_getitem_with_int_has_positional_metadata(self):
        s = "Sequence string !1@2#3?.,"
        length = len(s)
        seq = Sequence(s, metadata={'id': 'id', 'description': 'dsc'},
                       positional_metadata={'quality': np.arange(length)})

        eseq = Sequence("S", {'id': 'id', 'description': 'dsc'},
                        positional_metadata={'quality': np.array([0])})
        self.assertEqual(seq[0], eseq)

        eseq = Sequence(",", metadata={'id': 'id', 'description': 'dsc'},
                        positional_metadata={'quality':
                                             np.array([len(seq) - 1])})
        self.assertEqual(seq[len(seq) - 1], eseq)

        eseq = Sequence("t", metadata={'id': 'id', 'description': 'dsc'},
                        positional_metadata={'quality': [10]})
        self.assertEqual(seq[10], eseq)

    def test_single_index_to_slice(self):
        a = [1, 2, 3, 4]
        self.assertEqual(slice(0, 1), _single_index_to_slice(0))
        self.assertEqual([1], a[_single_index_to_slice(0)])
        self.assertEqual(slice(-1, None),
                         _single_index_to_slice(-1))
        self.assertEqual([4], a[_single_index_to_slice(-1)])

    def test_is_single_index(self):
        self.assertTrue(_is_single_index(0))
        self.assertFalse(_is_single_index(True))
        self.assertFalse(_is_single_index(bool()))
        self.assertFalse(_is_single_index('a'))

    def test_as_slice_if_single_index(self):
        self.assertEqual(slice(0, 1), _as_slice_if_single_index(0))
        slice_obj = slice(2, 3)
        self.assertIs(slice_obj,
                      _as_slice_if_single_index(slice_obj))

    def test_slice_positional_metadata(self):
        seq = Sequence('ABCDEFGHIJ',
                       positional_metadata={'foo': np.arange(10),
                                            'bar': np.arange(100, 110)})
        self.assertTrue(pd.DataFrame({'foo': [0], 'bar': [100]}).equals(
                        seq._slice_positional_metadata(0)))
        self.assertTrue(pd.DataFrame({'foo': [0], 'bar': [100]}).equals(
                        seq._slice_positional_metadata(slice(0, 1))))
        self.assertTrue(pd.DataFrame({'foo': [0, 1],
                                      'bar': [100, 101]}).equals(
                        seq._slice_positional_metadata(slice(0, 2))))
        self.assertTrue(pd.DataFrame(
            {'foo': [9], 'bar': [109]}, index=[9]).equals(
                seq._slice_positional_metadata(9)))

    def test_getitem_with_int_no_positional_metadata(self):
        seq = Sequence("Sequence string !1@2#3?.,",
                       metadata={'id': 'id2', 'description': 'no_qual'})

        eseq = Sequence("t", metadata={'id': 'id2', 'description': 'no_qual'})
        self.assertEqual(seq[10], eseq)

    def test_getitem_with_slice_has_positional_metadata(self):
        s = "0123456789abcdef"
        length = len(s)
        seq = Sequence(s, metadata={'id': 'id3', 'description': 'dsc3'},
                       positional_metadata={'quality': np.arange(length)})

        eseq = Sequence("012", metadata={'id': 'id3', 'description': 'dsc3'},
                        positional_metadata={'quality': np.arange(3)})
        self.assertEqual(seq[0:3], eseq)
        self.assertEqual(seq[:3], eseq)
        self.assertEqual(seq[:3:1], eseq)

        eseq = Sequence("def", metadata={'id': 'id3', 'description': 'dsc3'},
                        positional_metadata={'quality': [13, 14, 15]})
        self.assertEqual(seq[-3:], eseq)
        self.assertEqual(seq[-3::1], eseq)

        eseq = Sequence("02468ace",
                        metadata={'id': 'id3', 'description': 'dsc3'},
                        positional_metadata={'quality': [0, 2, 4, 6, 8, 10,
                                                         12, 14]})
        self.assertEqual(seq[0:length:2], eseq)
        self.assertEqual(seq[::2], eseq)

        eseq = Sequence(s[::-1], metadata={'id': 'id3', 'description': 'dsc3'},
                        positional_metadata={'quality':
                                             np.arange(length)[::-1]})
        self.assertEqual(seq[length::-1], eseq)
        self.assertEqual(seq[::-1], eseq)

        eseq = Sequence('fdb97531',
                        metadata={'id': 'id3', 'description': 'dsc3'},
                        positional_metadata={'quality': [15, 13, 11, 9, 7, 5,
                                                         3, 1]})
        self.assertEqual(seq[length::-2], eseq)
        self.assertEqual(seq[::-2], eseq)

        self.assertEqual(seq[0:500:], seq)

        eseq = Sequence('', metadata={'id': 'id3', 'description': 'dsc3'},
                        positional_metadata={'quality':
                                             np.array([], dtype=np.int64)})
        self.assertEqual(seq[length:0], eseq)
        self.assertEqual(seq[-length:0], eseq)
        self.assertEqual(seq[1:0], eseq)

        eseq = Sequence("0", metadata={'id': 'id3', 'description': 'dsc3'},
                        positional_metadata={'quality': [0]})
        self.assertEqual(seq[0:1], eseq)
        self.assertEqual(seq[0:1:1], eseq)
        self.assertEqual(seq[-length::-1], eseq)

    def test_getitem_with_slice_no_positional_metadata(self):
        s = "0123456789abcdef"
        length = len(s)
        seq = Sequence(s, metadata={'id': 'id4', 'description': 'no_qual4'})

        eseq = Sequence("02468ace",
                        metadata={'id': 'id4', 'description': 'no_qual4'})
        self.assertEqual(seq[0:length:2], eseq)
        self.assertEqual(seq[::2], eseq)

    def test_getitem_with_tuple_of_mixed_with_positional_metadata(self):
        s = "0123456789abcdef"
        length = len(s)
        seq = Sequence(s, metadata={'id': 'id5', 'description': 'dsc5'},
                       positional_metadata={'quality': np.arange(length)})

        eseq = Sequence("00000", metadata={'id': 'id5', 'description': 'dsc5'},
                        positional_metadata={'quality': [0, 0, 0, 0, 0]})
        self.assertEqual(seq[0, 0, 0, 0, 0], eseq)
        self.assertEqual(seq[0, 0:1, 0, 0, 0], eseq)
        self.assertEqual(seq[0, 0:1, 0, -length::-1, 0, 1:0], eseq)
        self.assertEqual(seq[0:1, 0:1, 0:1, 0:1, 0:1], eseq)
        self.assertEqual(seq[0:1, 0, 0, 0, 0], eseq)

        eseq = Sequence("0123fed9",
                        metadata={'id': 'id5', 'description': 'dsc5'},
                        positional_metadata={'quality': [0, 1, 2, 3, 15, 14,
                                                         13, 9]})
        self.assertEqual(seq[0, 1, 2, 3, 15, 14, 13, 9], eseq)
        self.assertEqual(seq[0, 1, 2, 3, :-4:-1, 9], eseq)
        self.assertEqual(seq[0:4, :-4:-1, 9, 1:0], eseq)
        self.assertEqual(seq[0:4, :-4:-1, 9:10], eseq)

    def test_getitem_with_tuple_of_mixed_no_positional_metadata(self):
        seq = Sequence("0123456789abcdef",
                       metadata={'id': 'id6', 'description': 'no_qual6'})
        eseq = Sequence("0123fed9",
                        metadata={'id': 'id6', 'description': 'no_qual6'})
        self.assertEqual(seq[0, 1, 2, 3, 15, 14, 13, 9], eseq)
        self.assertEqual(seq[0, 1, 2, 3, :-4:-1, 9], eseq)
        self.assertEqual(seq[0:4, :-4:-1, 9], eseq)
        self.assertEqual(seq[0:4, :-4:-1, 9:10], eseq)

    def test_getitem_with_tuple_of_mixed_no_metadata(self):
        seq = Sequence("0123456789abcdef")
        eseq = Sequence("0123fed9")
        self.assertEqual(seq[0, 1, 2, 3, 15, 14, 13, 9], eseq)
        self.assertEqual(seq[0, 1, 2, 3, :-4:-1, 9], eseq)
        self.assertEqual(seq[0:4, :-4:-1, 9], eseq)
        self.assertEqual(seq[0:4, :-4:-1, 9:10], eseq)

    def test_getitem_with_iterable_of_mixed_has_positional_metadata(self):
        s = "0123456789abcdef"
        length = len(s)
        seq = Sequence(s, metadata={'id': 'id7', 'description': 'dsc7'},
                       positional_metadata={'quality': np.arange(length)})

        def generator():
            yield slice(0, 4)
            yield slice(200, 400)
            yield -1
            yield slice(-2, -4, -1)
            yield 9

        eseq = Sequence("0123fed9",
                        metadata={'id': 'id7', 'description': 'dsc7'},
                        positional_metadata={'quality': [0, 1, 2, 3, 15, 14,
                                                         13, 9]})
        self.assertEqual(seq[[0, 1, 2, 3, 15, 14, 13, 9]], eseq)
        self.assertEqual(seq[generator()], eseq)
        self.assertEqual(seq[[slice(0, 4), slice(None, -4, -1), 9]], eseq)
        self.assertEqual(seq[
            [slice(0, 4), slice(None, -4, -1), slice(9, 10)]], eseq)

    def test_getitem_with_iterable_of_mixed_no_positional_metadata(self):
        s = "0123456789abcdef"
        seq = Sequence(s, metadata={'id': 'id7', 'description': 'dsc7'})

        def generator():
            yield slice(0, 4)
            yield slice(200, 400)
            yield slice(None, -4, -1)
            yield 9

        eseq = Sequence("0123fed9",
                        metadata={'id': 'id7', 'description': 'dsc7'})
        self.assertEqual(seq[[0, 1, 2, 3, 15, 14, 13, 9]], eseq)
        self.assertEqual(seq[generator()], eseq)
        self.assertEqual(seq[[slice(0, 4), slice(None, -4, -1), 9]], eseq)
        self.assertEqual(seq[
            [slice(0, 4), slice(None, -4, -1), slice(9, 10)]], eseq)

    def test_getitem_with_numpy_index_has_positional_metadata(self):
        s = "0123456789abcdef"
        length = len(s)
        seq = Sequence(s, metadata={'id': 'id9', 'description': 'dsc9'},
                       positional_metadata={'quality': np.arange(length)})

        eseq = Sequence("0123fed9",
                        metadata={'id': 'id9', 'description': 'dsc9'},
                        positional_metadata={'quality': [0, 1, 2, 3, 15, 14,
                                                         13, 9]})
        self.assertEqual(seq[np.array([0, 1, 2, 3, 15, 14, 13, 9])], eseq)

    def test_getitem_with_numpy_index_no_positional_metadata(self):
        s = "0123456789abcdef"
        seq = Sequence(s, metadata={'id': 'id10', 'description': 'dsc10'})

        eseq = Sequence("0123fed9",
                        metadata={'id': 'id10', 'description': 'dsc10'})
        self.assertEqual(seq[np.array([0, 1, 2, 3, 15, 14, 13, 9])], eseq)

    def test_getitem_with_empty_indices_empty_seq_no_pos_metadata(self):
        s = ""
        seq = Sequence(s, metadata={'id': 'id10', 'description': 'dsc10'})

        eseq = Sequence('', metadata={'id': 'id10', 'description': 'dsc10'})

        tested = 0
        for index in self.getitem_empty_indices:
            tested += 1
            self.assertEqual(seq[index], eseq)
        self.assertEqual(tested, 6)

    def test_getitem_with_empty_indices_non_empty_seq_no_pos_metadata(self):
        s = "0123456789abcdef"
        seq = Sequence(s, metadata={'id': 'id10', 'description': 'dsc10'})

        eseq = Sequence('', metadata={'id': 'id10', 'description': 'dsc10'})

        tested = 0
        for index in self.getitem_empty_indices:
            tested += 1
            self.assertEqual(seq[index], eseq)
        self.assertEqual(tested, 6)

    def test_getitem_with_boolean_vector_has_qual(self):
        s = "0123456789abcdef"
        length = len(s)
        seq = Sequence(s, metadata={'id': 'id11', 'description': 'dsc11'},
                       positional_metadata={'quality': np.arange(length)})

        eseq = Sequence("13579bdf",
                        metadata={'id': 'id11', 'description': 'dsc11'},
                        positional_metadata={'quality': [1, 3, 5, 7, 9, 11,
                                                         13, 15]})

        self.assertEqual(seq[np.array([False, True] * 8)], eseq)
        self.assertEqual(seq[[False, True] * 8], eseq)

    def test_getitem_with_boolean_vector_no_positional_metadata(self):
        s = "0123456789abcdef"
        seq = Sequence(s, metadata={'id': 'id11', 'description': 'dsc11'})

        eseq = Sequence("13579bdf",
                        metadata={'id': 'id11', 'description': 'dsc11'})

        self.assertEqual(seq[np.array([False, True] * 8)], eseq)

    def test_getitem_with_invalid(self):
        seq = Sequence("123456",
                       metadata={'id': 'idm', 'description': 'description'},
                       positional_metadata={'quality': [1, 2, 3, 4, 5, 6]})

        with self.assertRaises(IndexError):
            seq['not an index']

        with self.assertRaises(IndexError):
            seq[['1', '2']]

        with self.assertRaises(IndexError):
            seq[[1, slice(1, 2), 'a']]

        with self.assertRaises(IndexError):
            seq[[1, slice(1, 2), True]]

        with self.assertRaises(IndexError):
            seq[True]

        with self.assertRaises(IndexError):
            seq[np.array([True, False])]

        with self.assertRaises(IndexError):
            seq[999]

        with self.assertRaises(IndexError):
            seq[0, 0, 999]

        # numpy 1.8.1 and 1.9.2 raise different error types
        # (ValueError, IndexError).
        with self.assertRaises(Exception):
            seq[100 * [True, False, True]]

    def test_getitem_empty_positional_metadata(self):
        seq = Sequence('ACGT')
        seq.positional_metadata  # This will create empty positional_metadata
        self.assertEqual(Sequence('A'), seq[0])

    def test_len(self):
        self.assertEqual(len(Sequence("")), 0)
        self.assertEqual(len(Sequence("a")), 1)
        self.assertEqual(len(Sequence("abcdef")), 6)

    def test_nonzero(self):
        # blank
        self.assertFalse(Sequence(""))
        self.assertFalse(Sequence("",
                                  metadata={'id': 'foo'},
                                  positional_metadata={'quality': range(0)}))
        # single
        self.assertTrue(Sequence("A"))
        self.assertTrue(Sequence("A",
                                 metadata={'id': 'foo'},
                                 positional_metadata={'quality': range(1)}))
        # multi
        self.assertTrue(Sequence("ACGT"))
        self.assertTrue(Sequence("ACGT",
                                 metadata={'id': 'foo'},
                                 positional_metadata={'quality': range(4)}))

    def test_contains(self):
        seq = Sequence("#@ACGT,24.13**02")
        tested = 0
        for c in self.sequence_kinds:
            tested += 1
            self.assertTrue(c(',24') in seq)
            self.assertTrue(c('*') in seq)
            self.assertTrue(c('') in seq)

            self.assertFalse(c("$") in seq)
            self.assertFalse(c("AGT") in seq)

        self.assertEqual(tested, 4)

    def test_contains_sequence_subclass(self):
        with self.assertRaises(TypeError):
            SequenceSubclass("A") in Sequence("AAA")

        self.assertTrue(SequenceSubclass("A").values in Sequence("AAA"))

    def test_hash(self):
        with self.assertRaises(TypeError):
            hash(Sequence("ABCDEFG"))
        self.assertNotIsInstance(Sequence("ABCDEFG"), Hashable)

    def test_iter_has_positional_metadata(self):
        tested = False
        seq = Sequence("0123456789", metadata={'id': 'a', 'desc': 'b'},
                       positional_metadata={'qual': np.arange(10)})
        for i, s in enumerate(seq):
            tested = True
            self.assertEqual(s, Sequence(str(i),
                                         metadata={'id': 'a', 'desc': 'b'},
                                         positional_metadata={'qual': [i]}))
        self.assertTrue(tested)

    def test_iter_no_positional_metadata(self):
        tested = False
        seq = Sequence("0123456789", metadata={'id': 'a', 'desc': 'b'})
        for i, s in enumerate(seq):
            tested = True
            self.assertEqual(s, Sequence(str(i),
                                         metadata={'id': 'a', 'desc': 'b'}))
        self.assertTrue(tested)

    def test_reversed_has_positional_metadata(self):
        tested = False
        seq = Sequence("0123456789", metadata={'id': 'a', 'desc': 'b'},
                       positional_metadata={'qual': np.arange(10)})
        for i, s in enumerate(reversed(seq)):
            tested = True
            self.assertEqual(s, Sequence(str(9 - i),
                                         metadata={'id': 'a', 'desc': 'b'},
                                         positional_metadata={'qual':
                                                              [9 - i]}))
        self.assertTrue(tested)

    def test_reversed_no_positional_metadata(self):
        tested = False
        seq = Sequence("0123456789", metadata={'id': 'a', 'desc': 'b'})
        for i, s in enumerate(reversed(seq)):
            tested = True
            self.assertEqual(s, Sequence(str(9 - i),
                                         metadata={'id': 'a', 'desc': 'b'}))
        self.assertTrue(tested)

    def test_repr(self):
        # basic sanity checks -- more extensive testing of formatting and
        # special cases is performed in SequenceReprDoctests below. here we
        # only test that pieces of the repr are present. these tests also
        # exercise coverage in case doctests stop counting towards coverage in
        # the future

        # minimal
        obs = repr(Sequence(''))
        self.assertEqual(obs.count('\n'), 4)
        self.assertTrue(obs.startswith('Sequence'))
        self.assertIn('length: 0', obs)
        self.assertTrue(obs.endswith('-'))

        # no metadata
        obs = repr(Sequence('ACGT'))
        self.assertEqual(obs.count('\n'), 5)
        self.assertTrue(obs.startswith('Sequence'))
        self.assertIn('length: 4', obs)
        self.assertTrue(obs.endswith('0 ACGT'))

        # metadata and positional metadata of mixed types
        obs = repr(
            Sequence(
                'ACGT',
                metadata={'foo': 'bar', b'bar': 33.33, None: True, False: {},
                          (1, 2): 3, 'acb' * 100: "'", 10: 11},
                positional_metadata={'foo': range(4),
                                     42: ['a', 'b', [], 'c']}))
        self.assertEqual(obs.count('\n'), 16)
        self.assertTrue(obs.startswith('Sequence'))
        self.assertIn('None: True', obs)
        self.assertIn('\'foo\': \'bar\'', obs)
        self.assertIn('42: <dtype: object>', obs)
        self.assertIn('\'foo\': <dtype: int64>', obs)
        self.assertIn('length: 4', obs)
        self.assertTrue(obs.endswith('0 ACGT'))

        # sequence spanning > 5 lines
        obs = repr(Sequence('A' * 301))
        self.assertEqual(obs.count('\n'), 9)
        self.assertTrue(obs.startswith('Sequence'))
        self.assertIn('length: 301', obs)
        self.assertIn('...', obs)
        self.assertTrue(obs.endswith('300 A'))

    def test_str(self):
        self.assertEqual(str(Sequence("GATTACA")), "GATTACA")
        self.assertEqual(str(Sequence("ACCGGTACC")), "ACCGGTACC")
        self.assertEqual(str(Sequence("GREG")), "GREG")
        self.assertEqual(
            str(Sequence("ABC",
                         positional_metadata={'quality': [1, 2, 3]})),
            "ABC")
        self.assertIs(type(str(Sequence("A"))), str)

    def test_count(self):
        def construct_char_array(s):
            return np.frombuffer(s.encode('ascii'), dtype='|S1')

        def construct_uint8_array(s):
            return np.frombuffer(s.encode('ascii'), dtype=np.uint8)

        seq = Sequence("1234567899876555")
        tested = 0
        for c in self.sequence_kinds:
            tested += 1
            self.assertEqual(seq.count(c('4')), 1)
            self.assertEqual(seq.count(c('8')), 2)
            self.assertEqual(seq.count(c('5')), 4)
            self.assertEqual(seq.count(c('555')), 1)
            self.assertEqual(seq.count(c('555'), 0, 4), 0)
            self.assertEqual(seq.count(c('555'), start=0, end=4), 0)
            self.assertEqual(seq.count(c('5'), start=10), 3)
            self.assertEqual(seq.count(c('5'), end=10), 1)

            with self.assertRaises(ValueError):
                seq.count(c(''))

        self.assertEqual(tested, 4)

    def test_count_on_subclass(self):
        with self.assertRaises(TypeError) as cm:
            Sequence("abcd").count(SequenceSubclass("a"))

        self.assertIn("Sequence", str(cm.exception))
        self.assertIn("SequenceSubclass", str(cm.exception))

    def test_replace_sanity(self):
        seq = Sequence('AAGCATGCCCTTTACATTTG')
        index = self._make_index('10011011001111110111')
        obs = seq.replace(index, '_')
        exp = Sequence('_AG__T__CC______T___')
        self.assertEqual(obs, exp)

    def test_replace_index_array(self):
        seq = Sequence('TCGGGTGTTGTGCAACCACC')
        for _type in list, tuple, np.array, pd.Series:
            index = _type([0, 2, 5, 8, 9])
            obs = seq.replace(index, '-')
            exp = Sequence('-C-GG-GT--TGCAACCACC')
            self.assertEqual(obs, exp)

    def test_replace_iterable_slices(self):
        seq = Sequence('CATTATGGACCCAGCGTGCC')
        slices = (slice(0, 5), slice(8, 12), slice(15, 17))
        mixed_slices = (0, 1, 2, 3, 4, slice(8, 12), 15, 16)
        for _type in (lambda x: x, list, tuple, lambda x: np.array(tuple(x)),
                      lambda x: pd.Series(tuple(x))):
            index = (_type(slices), _type(mixed_slices))
            obs_slices = seq.replace(index[0], '-')
            obs_mixed = seq.replace(index[1], '-')
            exp = Sequence('-----TGG----AGC--GCC')
            self.assertEqual(obs_slices, exp)
            self.assertEqual(obs_mixed, exp)

    def test_replace_index_in_positional_metadata(self):
        positional_metadata = {'where': self._make_index('001110110'
                                                         '10001110000')}
        seq = Sequence('AAGATTGATACCACAGTTGT',
                       positional_metadata=positional_metadata)
        obs = seq.replace('where', '-')
        exp = Sequence('AA---T--T-CCA---TTGT',
                       positional_metadata=positional_metadata)
        self.assertEqual(obs, exp)

    def test_replace_does_not_mutate_original(self):
        seq = Sequence('ATCG')
        index = self._make_index('0011')
        seq.replace(index, '-')
        obs = seq
        exp = Sequence('ATCG')
        self.assertEqual(obs, exp)

    def test_replace_with_metadata(self):
        seq = Sequence('GCACGGCAAGAAGCGCCCCA',
                       metadata={'NM': 'Kestrel Gorlick'},
                       positional_metadata={'diff':
                                            list('01100001110010001100')})
        seq.interval_metadata.add([(0, 1)], metadata={'gene': 'sagA'})

        index = self._make_index('01100001110010001100')
        obs = seq.replace(index, '-')
        exp = Sequence('G--CGGC---AA-CGC--CA',
                       metadata={'NM': 'Kestrel Gorlick'},
                       positional_metadata={'diff':
                                            list('01100001110010001100')})
        exp.interval_metadata.add([(0, 1)], metadata={'gene': 'sagA'})

        self.assertEqual(obs, exp)

    def test_replace_with_subclass(self):
        seq = DNA('CGACAACCGATGTGCTGTAA')
        index = self._make_index('10101000111111110011')
        obs = seq.replace(index, '-')
        exp = DNA('-G-C-ACC--------GT--')
        self.assertEqual(obs, exp)

    def test_replace_with_bytes(self):
        seq = Sequence('ABC123')

        obs = seq.replace([1, 3, 5], b'*')

        self.assertEqual(obs, Sequence('A*C*2*'))

    def test_replace_invalid_char_for_type_error(self):
        seq = DNA('TAAACGGAACGCTACGTCTG')
        index = self._make_index('01000001101011001001')
        with self.assertRaisesRegex(ValueError, r"Invalid character.*'F'"):
            seq.replace(index, 'F')

    def test_replace_invalid_char_error(self):
        seq = Sequence('GGGAGCTAGA')
        index = self._make_index('1000101110')
        with self.assertRaisesRegex(UnicodeEncodeError,
                                    r"can't encode character.*not in "
                                    r"range\(128\)"):
            seq.replace(index, '\uFFFF')

    def test_replace_non_single_character_error(self):
        seq = Sequence('CCGAACTGTC')
        index = self._make_index('1100110011')
        with self.assertRaisesRegex(TypeError, r'string of length 2 found'):
            seq.replace(index, 'AB')

    def _make_index(self, bools):
        return [bool(int(char)) for char in bools]

    def test_lowercase_mungeable_key(self):
        # NOTE: This test relies on Sequence._munge_to_index_array working
        # properly. If the internal implementation of the lowercase method
        # changes to no longer use _munge_to_index_array, this test may need
        # to be updated to cover cases currently covered by
        # _munge_to_index_array
        self.assertEqual('AAAAaaaa', self.lowercase_seq.lowercase('key'))

    def test_lowercase_array_key(self):
        # NOTE: This test relies on Sequence._munge_to_index_array working
        # properly. If the internal implementation of the lowercase method
        # changes to no longer use _munge_to_index_array, this test may need
        # to be updated to cover cases currently covered by
        # _munge_to_index_array
        self.assertEqual('aaAAaaaa',
                         self.lowercase_seq.lowercase(
                             np.array([True, True, False, False, True, True,
                                       True, True])))
        self.assertEqual('AaAAaAAA',
                         self.lowercase_seq.lowercase([1, 4]))

    def test_matches(self):
        tested = 0
        for constructor in self.sequence_kinds:
            tested += 1
            seq1 = Sequence("AACCEEGG")
            seq2 = constructor("ABCDEFGH")
            expected = np.array([True, False] * 4)
            npt.assert_equal(seq1.matches(seq2), expected)

        self.assertEqual(tested, 4)

    def test_matches_on_subclass(self):
        seq1 = Sequence("AACCEEGG")
        seq2 = SequenceSubclass("ABCDEFGH")

        with self.assertRaises(TypeError):
            seq1.matches(seq2)

    def test_matches_unequal_length(self):
        seq1 = Sequence("AACCEEGG")
        seq2 = Sequence("TOOLONGTOCOMPARE")

        with self.assertRaises(ValueError):
            seq1.matches(seq2)

    def test_mismatches(self):
        tested = 0
        for constructor in self.sequence_kinds:
            tested += 1
            seq1 = Sequence("AACCEEGG")
            seq2 = constructor("ABCDEFGH")
            expected = np.array([False, True] * 4)
            npt.assert_equal(seq1.mismatches(seq2), expected)

        self.assertEqual(tested, 4)

    def test_mismatches_on_subclass(self):
        seq1 = Sequence("AACCEEGG")
        seq2 = SequenceSubclass("ABCDEFGH")

        with self.assertRaises(TypeError):
            seq1.mismatches(seq2)

    def test_mismatches_unequal_length(self):
        seq1 = Sequence("AACCEEGG")
        seq2 = Sequence("TOOLONGTOCOMPARE")

        with self.assertRaises(ValueError):
            seq1.mismatches(seq2)

    def test_mismatch_frequency(self):
        seq1 = Sequence("AACCEEGG")
        seq2 = Sequence("ABCDEFGH")
        seq3 = Sequence("TTTTTTTT")

        self.assertIs(type(seq1.mismatch_frequency(seq1)), int)
        self.assertEqual(seq1.mismatch_frequency(seq1), 0)
        self.assertEqual(seq1.mismatch_frequency(seq2), 4)
        self.assertEqual(seq1.mismatch_frequency(seq3), 8)

    def test_mismatch_frequency_relative(self):
        seq1 = Sequence("AACCEEGG")
        seq2 = Sequence("ABCDEFGH")
        seq3 = Sequence("TTTTTTTT")

        self.assertIs(type(seq1.mismatch_frequency(seq1, relative=True)),
                      float)
        self.assertEqual(seq1.mismatch_frequency(seq1, relative=True), 0.0)
        self.assertEqual(seq1.mismatch_frequency(seq2, relative=True), 0.5)
        self.assertEqual(seq1.mismatch_frequency(seq3, relative=True), 1.0)

    def test_mismatch_frequency_unequal_length(self):
        seq1 = Sequence("AACCEEGG")
        seq2 = Sequence("TOOLONGTOCOMPARE")

        with self.assertRaises(ValueError):
            seq1.mismatch_frequency(seq2)

    def test_mismatch_frequence_on_subclass(self):
        seq1 = Sequence("AACCEEGG")
        seq2 = SequenceSubclass("ABCDEFGH")

        with self.assertRaises(TypeError):
            seq1.mismatch_frequency(seq2)

    def test_match_frequency(self):
        seq1 = Sequence("AACCEEGG")
        seq2 = Sequence("ABCDEFGH")
        seq3 = Sequence("TTTTTTTT")

        self.assertIs(type(seq1.match_frequency(seq1)), int)
        self.assertEqual(seq1.match_frequency(seq1), 8)
        self.assertEqual(seq1.match_frequency(seq2), 4)
        self.assertEqual(seq1.match_frequency(seq3), 0)

    def test_match_frequency_relative(self):
        seq1 = Sequence("AACCEEGG")
        seq2 = Sequence("ABCDEFGH")
        seq3 = Sequence("TTTTTTTT")

        self.assertIs(type(seq1.match_frequency(seq1, relative=True)),
                      float)
        self.assertEqual(seq1.match_frequency(seq1, relative=True), 1.0)
        self.assertEqual(seq1.match_frequency(seq2, relative=True), 0.5)
        self.assertEqual(seq1.match_frequency(seq3, relative=True), 0.0)

    def test_match_frequency_unequal_length(self):
        seq1 = Sequence("AACCEEGG")
        seq2 = Sequence("TOOLONGTOCOMPARE")

        with self.assertRaises(ValueError):
            seq1.match_frequency(seq2)

    def test_match_frequency_on_subclass(self):
        seq1 = Sequence("AACCEEGG")
        seq2 = SequenceSubclass("ABCDEFGH")

        with self.assertRaises(TypeError):
            seq1.match_frequency(seq2)

    def test_index(self):
        tested = 0
        for c in self.sequence_kinds:
            tested += 1
            seq = Sequence("ABCDEFG@@ABCDFOO")
            self.assertEqual(seq.index(c("A")), 0)
            self.assertEqual(seq.index(c("@")), 7)
            self.assertEqual(seq.index(c("@@")), 7)

            with self.assertRaises(ValueError):
                seq.index("A", start=1, end=5)

        self.assertEqual(tested, 4)

    def test_index_on_subclass(self):
        with self.assertRaises(TypeError):
            Sequence("ABCDEFG").index(SequenceSubclass("A"))

        self.assertEqual(
            SequenceSubclass("ABCDEFG").index(SequenceSubclass("A")), 0)

    def test_frequencies_empty_sequence(self):
        seq = Sequence('')

        self.assertEqual(seq.frequencies(), {})
        self.assertEqual(seq.frequencies(relative=True), {})

        self.assertEqual(seq.frequencies(chars=set()), {})
        self.assertEqual(seq.frequencies(chars=set(), relative=True), {})

        self.assertEqual(seq.frequencies(chars={'a', 'b'}), {'a': 0, 'b': 0})

        # use npt.assert_equal to explicitly handle nan comparisons
        npt.assert_equal(seq.frequencies(chars={'a', 'b'}, relative=True),
                         {'a': np.nan, 'b': np.nan})

    def test_frequencies_observed_chars(self):
        seq = Sequence('x')
        self.assertEqual(seq.frequencies(), {'x': 1})
        self.assertEqual(seq.frequencies(relative=True), {'x': 1.0})

        seq = Sequence('xYz')
        self.assertEqual(seq.frequencies(), {'x': 1, 'Y': 1, 'z': 1})
        self.assertEqual(seq.frequencies(relative=True),
                         {'x': 1/3, 'Y': 1/3, 'z': 1/3})

        seq = Sequence('zzz')
        self.assertEqual(seq.frequencies(), {'z': 3})
        self.assertEqual(seq.frequencies(relative=True), {'z': 1.0})

        seq = Sequence('xYzxxZz')
        self.assertEqual(seq.frequencies(), {'x': 3, 'Y': 1, 'Z': 1, 'z': 2})
        self.assertEqual(seq.frequencies(relative=True),
                         {'x': 3/7, 'Y': 1/7, 'Z': 1/7, 'z': 2/7})

        seq = Sequence('\t   ')
        self.assertEqual(seq.frequencies(), {'\t': 1, ' ': 3})
        self.assertEqual(seq.frequencies(relative=True), {'\t': 1/4, ' ': 3/4})

        seq = Sequence('aabbcc', metadata={'foo': 'bar'},
                       positional_metadata={'foo': range(6)})
        self.assertEqual(seq.frequencies(), {'a': 2, 'b': 2, 'c': 2})
        self.assertEqual(seq.frequencies(relative=True),
                         {'a': 2/6, 'b': 2/6, 'c': 2/6})

    def test_frequencies_specified_chars(self):
        seq = Sequence('abcbca')

        self.assertEqual(seq.frequencies(chars=set()), {})
        self.assertEqual(seq.frequencies(chars=set(), relative=True), {})

        self.assertEqual(seq.frequencies(chars='a'), {'a': 2})
        self.assertEqual(seq.frequencies(chars='a', relative=True), {'a': 2/6})

        self.assertEqual(seq.frequencies(chars={'a'}), {'a': 2})
        self.assertEqual(seq.frequencies(chars={'a'}, relative=True),
                         {'a': 2/6})

        self.assertEqual(seq.frequencies(chars={'a', 'b'}), {'a': 2, 'b': 2})
        self.assertEqual(seq.frequencies(chars={'a', 'b'}, relative=True),
                         {'a': 2/6, 'b': 2/6})

        self.assertEqual(seq.frequencies(chars={'a', 'b', 'd'}),
                         {'a': 2, 'b': 2, 'd': 0})
        self.assertEqual(seq.frequencies(chars={'a', 'b', 'd'}, relative=True),
                         {'a': 2/6, 'b': 2/6, 'd': 0.0})

        self.assertEqual(seq.frequencies(chars={'x', 'y', 'z'}),
                         {'x': 0, 'y': 0, 'z': 0})
        self.assertEqual(seq.frequencies(chars={'x', 'y', 'z'}, relative=True),
                         {'x': 0.0, 'y': 0.0, 'z': 0.0})

    def test_frequencies_chars_varied_type(self):
        seq = Sequence('zabczzzabcz')

        # single character case (shortcut)
        chars = b'z'
        self.assertEqual(seq.frequencies(chars=chars), {b'z': 5})
        self.assertEqual(seq.frequencies(chars=chars, relative=True),
                         {b'z': 5/11})

        chars = 'z'
        self.assertEqual(seq.frequencies(chars=chars), {'z': 5})
        self.assertEqual(seq.frequencies(chars=chars, relative=True),
                         {'z': 5/11})

        chars = np.frombuffer('z'.encode('ascii'), dtype='|S1')[0]
        self.assertEqual(seq.frequencies(chars=chars), {b'z': 5})
        self.assertEqual(seq.frequencies(chars=chars, relative=True),
                         {b'z': 5/11})

        # set of characters, some present, some not
        chars = {b'x', b'z'}
        self.assertEqual(seq.frequencies(chars=chars), {b'x': 0, b'z': 5})
        self.assertEqual(seq.frequencies(chars=chars, relative=True),
                         {b'x': 0.0, b'z': 5/11})

        chars = {'x', 'z'}
        self.assertEqual(seq.frequencies(chars=chars), {'x': 0, 'z': 5})
        self.assertEqual(seq.frequencies(chars=chars, relative=True),
                         {'x': 0.0, 'z': 5/11})

        chars = {
            np.frombuffer('x'.encode('ascii'), dtype='|S1')[0],
            np.frombuffer('z'.encode('ascii'), dtype='|S1')[0]
        }
        self.assertEqual(seq.frequencies(chars=chars), {b'x': 0, b'z': 5})
        self.assertEqual(seq.frequencies(chars=chars, relative=True),
                         {b'x': 0.0, b'z': 5/11})

    def test_frequencies_equivalent_to_kmer_frequencies_k_of_1(self):
        seq = Sequence('abcabc')

        exp = {'a': 2, 'b': 2, 'c': 2}
        self.assertEqual(seq.frequencies(chars=None), exp)
        self.assertEqual(seq.kmer_frequencies(k=1), exp)

        exp = {'a': 2/6, 'b': 2/6, 'c': 2/6}
        self.assertEqual(seq.frequencies(chars=None, relative=True), exp)
        self.assertEqual(seq.kmer_frequencies(k=1, relative=True), exp)

    def test_frequencies_passing_observed_chars_equivalent_to_default(self):
        seq = Sequence('abcabc')

        exp = {'a': 2, 'b': 2, 'c': 2}
        self.assertEqual(seq.frequencies(chars=None), exp)
        self.assertEqual(seq.frequencies(chars=seq.observed_chars), exp)

        exp = {'a': 2/6, 'b': 2/6, 'c': 2/6}
        self.assertEqual(seq.frequencies(chars=None, relative=True), exp)
        self.assertEqual(
            seq.frequencies(chars=seq.observed_chars, relative=True),
            exp)

    def test_frequencies_invalid_chars(self):
        seq = Sequence('abcabc')

        with self.assertRaisesRegex(ValueError, r'0 characters'):
            seq.frequencies(chars='')

        with self.assertRaisesRegex(ValueError, r'0 characters'):
            seq.frequencies(chars={''})

        with self.assertRaisesRegex(ValueError, r'2 characters'):
            seq.frequencies(chars='ab')

        with self.assertRaisesRegex(ValueError, r'2 characters'):
            seq.frequencies(chars={'b', 'ab'})

        with self.assertRaisesRegex(TypeError, r'string.*NoneType'):
            seq.frequencies(chars={'a', None})

        with self.assertRaisesRegex(ValueError, r'outside the range'):
            seq.frequencies(chars='\u1F30')

        with self.assertRaisesRegex(ValueError, r'outside the range'):
            seq.frequencies(chars={'c', '\u1F30'})

        with self.assertRaisesRegex(TypeError, r'set.*int'):
            seq.frequencies(chars=42)

    def _compare_kmers_results(self, observed, expected):
        for obs, exp in itertools.zip_longest(observed, expected,
                                              fillvalue=None):
            self.assertEqual(obs, exp)

    def test_iter_kmers(self):
        seq = Sequence('GATTACA', positional_metadata={'quality': range(7)})

        expected = [
            Sequence('G', positional_metadata={'quality': [0]}),
            Sequence('A', positional_metadata={'quality': [1]}),
            Sequence('T', positional_metadata={'quality': [2]}),
            Sequence('T', positional_metadata={'quality': [3]}),
            Sequence('A', positional_metadata={'quality': [4]}),
            Sequence('C', positional_metadata={'quality': [5]}),
            Sequence('A', positional_metadata={'quality': [6]})
        ]
        self._compare_kmers_results(
            seq.iter_kmers(1, overlap=False), expected)

        expected = [
            Sequence('GA', positional_metadata={'quality': [0, 1]}),
            Sequence('TT', positional_metadata={'quality': [2, 3]}),
            Sequence('AC', positional_metadata={'quality': [4, 5]})
        ]
        self._compare_kmers_results(
            seq.iter_kmers(2, overlap=False), expected)

        expected = [
            Sequence('GAT', positional_metadata={'quality': [0, 1, 2]}),
            Sequence('TAC', positional_metadata={'quality': [3, 4, 5]})
        ]
        self._compare_kmers_results(
            seq.iter_kmers(3, overlap=False), expected)

        expected = [
            Sequence('GATTACA',
                     positional_metadata={'quality': [0, 1, 2, 3, 4, 5, 6]})
        ]
        self._compare_kmers_results(
            seq.iter_kmers(7, overlap=False), expected)

        expected = []
        self._compare_kmers_results(
            seq.iter_kmers(8, overlap=False), expected)

        self.assertIs(type(seq.iter_kmers(1)), GeneratorType)

    def test_iter_kmers_no_positional_metadata(self):
        seq = Sequence('GATTACA')

        expected = [
            Sequence('G'),
            Sequence('A'),
            Sequence('T'),
            Sequence('T'),
            Sequence('A'),
            Sequence('C'),
            Sequence('A')
        ]
        self._compare_kmers_results(
            seq.iter_kmers(1, overlap=False), expected)

        expected = [
            Sequence('GA'),
            Sequence('TT'),
            Sequence('AC')
        ]
        self._compare_kmers_results(
            seq.iter_kmers(2, overlap=False), expected)

        expected = [
            Sequence('GAT'),
            Sequence('TAC')
        ]
        self._compare_kmers_results(
            seq.iter_kmers(3, overlap=False), expected)

        expected = [
            Sequence('GATTACA')
        ]
        self._compare_kmers_results(
            seq.iter_kmers(7, overlap=False), expected)

        expected = []
        self._compare_kmers_results(
            seq.iter_kmers(8, overlap=False), expected)

        self.assertIs(type(seq.iter_kmers(1)), GeneratorType)

    def test_iter_kmers_with_overlap(self):
        seq = Sequence('GATTACA', positional_metadata={'quality': range(7)})
        expected = [
            Sequence('G', positional_metadata={'quality': [0]}),
            Sequence('A', positional_metadata={'quality': [1]}),
            Sequence('T', positional_metadata={'quality': [2]}),
            Sequence('T', positional_metadata={'quality': [3]}),
            Sequence('A', positional_metadata={'quality': [4]}),
            Sequence('C', positional_metadata={'quality': [5]}),
            Sequence('A', positional_metadata={'quality': [6]})
        ]
        self._compare_kmers_results(
            seq.iter_kmers(1, overlap=True), expected)

        expected = [
            Sequence('GA', positional_metadata={'quality': [0, 1]}),
            Sequence('AT', positional_metadata={'quality': [1, 2]}),
            Sequence('TT', positional_metadata={'quality': [2, 3]}),
            Sequence('TA', positional_metadata={'quality': [3, 4]}),
            Sequence('AC', positional_metadata={'quality': [4, 5]}),
            Sequence('CA', positional_metadata={'quality': [5, 6]})
        ]
        self._compare_kmers_results(
            seq.iter_kmers(2, overlap=True), expected)

        expected = [
            Sequence('GAT', positional_metadata={'quality': [0, 1, 2]}),
            Sequence('ATT', positional_metadata={'quality': [1, 2, 3]}),
            Sequence('TTA', positional_metadata={'quality': [2, 3, 4]}),
            Sequence('TAC', positional_metadata={'quality': [3, 4, 5]}),
            Sequence('ACA', positional_metadata={'quality': [4, 5, 6]})
        ]
        self._compare_kmers_results(
            seq.iter_kmers(3, overlap=True), expected)

        expected = [
            Sequence('GATTACA',
                     positional_metadata={'quality': [0, 1, 2, 3, 4, 5, 6]})
        ]
        self._compare_kmers_results(
            seq.iter_kmers(7, overlap=True), expected)

        expected = []
        self._compare_kmers_results(
            seq.iter_kmers(8, overlap=True), expected)

        self.assertIs(type(seq.iter_kmers(1)), GeneratorType)

    def test_iter_kmers_with_overlap_no_positional_metadata(self):
        seq = Sequence('GATTACA')
        expected = [
            Sequence('G'),
            Sequence('A'),
            Sequence('T'),
            Sequence('T'),
            Sequence('A'),
            Sequence('C'),
            Sequence('A')
        ]
        self._compare_kmers_results(
            seq.iter_kmers(1, overlap=True), expected)

        expected = [
            Sequence('GA'),
            Sequence('AT'),
            Sequence('TT'),
            Sequence('TA'),
            Sequence('AC'),
            Sequence('CA')
        ]
        self._compare_kmers_results(
            seq.iter_kmers(2, overlap=True), expected)

        expected = [
            Sequence('GAT'),
            Sequence('ATT'),
            Sequence('TTA'),
            Sequence('TAC'),
            Sequence('ACA')
        ]
        self._compare_kmers_results(
            seq.iter_kmers(3, overlap=True), expected)

        expected = [
            Sequence('GATTACA')
        ]
        self._compare_kmers_results(
            seq.iter_kmers(7, overlap=True), expected)

        expected = []
        self._compare_kmers_results(
            seq.iter_kmers(8, overlap=True), expected)

        self.assertIs(type(seq.iter_kmers(1)), GeneratorType)

    def test_iter_kmers_invalid_k(self):
        seq = Sequence('GATTACA', positional_metadata={'quality': range(7)})

        with self.assertRaises(ValueError):
            list(seq.iter_kmers(0))

        with self.assertRaises(ValueError):
            list(seq.iter_kmers(-42))

    def test_iter_kmers_invalid_k_no_positional_metadata(self):
        seq = Sequence('GATTACA')

        with self.assertRaises(ValueError):
            list(seq.iter_kmers(0))

        with self.assertRaises(ValueError):
            list(seq.iter_kmers(-42))

    def test_iter_kmers_different_sequences(self):
        seq = Sequence('HE..--..LLO',
                       metadata={'id': 'hello', 'desc': 'gapped hello'},
                       positional_metadata={'quality': range(11)})
        expected = [
            Sequence('HE.', positional_metadata={'quality': [0, 1, 2]},
                     metadata={'id': 'hello', 'desc': 'gapped hello'}),
            Sequence('.--', positional_metadata={'quality': [3, 4, 5]},
                     metadata={'id': 'hello', 'desc': 'gapped hello'}),
            Sequence('..L', positional_metadata={'quality': [6, 7, 8]},
                     metadata={'id': 'hello', 'desc': 'gapped hello'})
        ]
        self._compare_kmers_results(seq.iter_kmers(3, overlap=False), expected)

    def test_iter_kmers_different_sequences_no_positional_metadata(self):
        seq = Sequence('HE..--..LLO',
                       metadata={'id': 'hello', 'desc': 'gapped hello'})
        expected = [
            Sequence('HE.',
                     metadata={'id': 'hello', 'desc': 'gapped hello'}),
            Sequence('.--',
                     metadata={'id': 'hello', 'desc': 'gapped hello'}),
            Sequence('..L',
                     metadata={'id': 'hello', 'desc': 'gapped hello'})
        ]
        self._compare_kmers_results(seq.iter_kmers(3, overlap=False), expected)

    def test_iter_kmers_empty_sequence(self):
        seq = Sequence('')
        expected = []
        self._compare_kmers_results(seq.iter_kmers(3, overlap=False), expected)

    def test_iter_kmers_empty_sequence_with_positional_metadata(self):
        seq = Sequence('', positional_metadata={'quality': []})
        expected = []
        self._compare_kmers_results(seq.iter_kmers(3, overlap=False), expected)

    def test_kmer_frequencies_empty_sequence(self):
        seq = Sequence('')

        self.assertEqual(seq.kmer_frequencies(1), {})
        self.assertEqual(seq.kmer_frequencies(1, overlap=False), {})
        self.assertEqual(seq.kmer_frequencies(1, relative=True), {})
        self.assertEqual(seq.kmer_frequencies(1, relative=True, overlap=False),
                         {})

    def test_kmer_frequencies(self):
        seq = Sequence('GATTACA', positional_metadata={'quality': range(7)})

        # overlap = True
        expected = {'G': 1, 'A': 3, 'T': 2, 'C': 1}
        self.assertEqual(seq.kmer_frequencies(1, overlap=True), expected)

        expected = {'GAT': 1, 'ATT': 1, 'TTA': 1, 'TAC': 1, 'ACA': 1}
        self.assertEqual(seq.kmer_frequencies(3, overlap=True), expected)

        expected = {}
        self.assertEqual(seq.kmer_frequencies(8, overlap=True), expected)

        # overlap = False
        expected = {'GAT': 1, 'TAC': 1}
        self.assertEqual(seq.kmer_frequencies(3, overlap=False), expected)

        expected = {'GATTACA': 1}
        self.assertEqual(seq.kmer_frequencies(7, overlap=False), expected)

        expected = {}
        self.assertEqual(seq.kmer_frequencies(8, overlap=False), expected)

    def test_kmer_frequencies_relative(self):
        seq = Sequence('GATTACA', positional_metadata={'quality': range(7)})

        # overlap = True
        expected = {'A': 3/7, 'C': 1/7, 'G': 1/7, 'T': 2/7}
        self.assertEqual(seq.kmer_frequencies(1, overlap=True, relative=True),
                         expected)

        expected = {'GAT': 1/5, 'ATT': 1/5, 'TTA': 1/5, 'TAC': 1/5, 'ACA': 1/5}
        self.assertEqual(seq.kmer_frequencies(3, overlap=True, relative=True),
                         expected)

        expected = {}
        self.assertEqual(seq.kmer_frequencies(8, overlap=True, relative=True),
                         expected)

        # overlap = False
        expected = {'GAT': 1/2, 'TAC': 1/2}
        self.assertEqual(seq.kmer_frequencies(3, overlap=False, relative=True),
                         expected)

        expected = {'GATTACA': 1.0}
        self.assertEqual(seq.kmer_frequencies(7, overlap=False, relative=True),
                         expected)

        expected = {}
        self.assertEqual(seq.kmer_frequencies(8, overlap=False, relative=True),
                         expected)

    def test_kmer_frequencies_floating_point_precision(self):
        # Test that a sequence having no variation in k-words yields a
        # frequency of exactly 1.0. Note that it is important to use
        # self.assertEqual here instead of self.assertAlmostEqual because we
        # want to test for exactly 1.0. A previous implementation of
        # Sequence.kmer_frequencies(relative=True) added (1 / num_words) for
        # each occurrence of a k-word to compute the frequencies (see
        # https://github.com/biocore/scikit-bio/issues/801). In certain cases,
        # this yielded a frequency slightly less than 1.0 due to roundoff
        # error. The test case here uses a sequence with 10 characters that are
        # all identical and computes k-word frequencies with k=1. This test
        # case exposes the roundoff error present in the previous
        # implementation because there are 10 k-words (which are all
        # identical), so 1/10 added 10 times yields a number slightly less than
        # 1.0. This occurs because 1/10 cannot be represented exactly as a
        # floating point number.
        seq = Sequence('AAAAAAAAAA')
        self.assertEqual(seq.kmer_frequencies(1, relative=True), {'A': 1.0})

    def test_find_with_regex(self):
        seq = Sequence('GATTACA', positional_metadata={'quality': range(7)})
        pat = re.compile('(T+A)(CA)')

        obs = list(seq.find_with_regex(pat))
        exp = [slice(2, 5), slice(5, 7)]
        self.assertEqual(obs, exp)

        self.assertIs(type(seq.find_with_regex(pat)), GeneratorType)

    def test_find_with_regex_string_as_input(self):
        seq = Sequence('GATTACA', positional_metadata={'quality': range(7)})
        pat = '(T+A)(CA)'

        obs = list(seq.find_with_regex(pat))
        exp = [slice(2, 5), slice(5, 7)]
        self.assertEqual(obs, exp)

        self.assertIs(type(seq.find_with_regex(pat)), GeneratorType)

    def test_find_with_regex_no_groups(self):
        seq = Sequence('GATTACA', positional_metadata={'quality': range(7)})
        pat = re.compile('(FOO)')
        self.assertEqual(list(seq.find_with_regex(pat)), [])

    def test_find_with_regex_ignore_no_difference(self):
        seq = Sequence('..ABCDEFG..')
        pat = "([A-Z]+)"
        exp = [slice(2, 9)]
        self.assertEqual(list(seq.find_with_regex(pat)), exp)

        obs = seq.find_with_regex(
            pat, ignore=np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1],
                                 dtype=bool))
        self.assertEqual(list(obs), exp)

    def test_find_with_regex_ignore(self):
        obs = Sequence('A..A..BBAAB.A..AB..A.').find_with_regex(
            "(A+)", ignore=np.array([0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 1,
                                     1, 0, 0, 1, 1, 0, 1], dtype=bool))

        self.assertEqual(list(obs), [slice(0, 4), slice(8, 10), slice(12, 16),
                                     slice(19, 20)])

    def test_find_with_regex_ignore_index_array(self):
        obs = Sequence('A..A..BBAAB.A..AB..A.').find_with_regex(
            "(A+)", ignore=np.array([1, 2, 4, 5, 11, 13, 14, 17, 18, 20]))

        self.assertEqual(list(obs), [slice(0, 4), slice(8, 10), slice(12, 16),
                                     slice(19, 20)])

    def test_iter_contiguous_index_array(self):
        s = Sequence("0123456789abcdef")
        for c in list, tuple, np.array, pd.Series:
            exp = [Sequence("0123"), Sequence("89ab")]
            obs = s.iter_contiguous(c([0, 1, 2, 3, 8, 9, 10, 11]))
            self.assertEqual(list(obs), exp)

    def test_iter_contiguous_boolean_vector(self):
        s = Sequence("0123456789abcdef")
        for c in list, tuple, np.array, pd.Series:
            exp = [Sequence("0123"), Sequence("89ab")]
            obs = s.iter_contiguous(c(([True] * 4 + [False] * 4) * 2))
            self.assertEqual(list(obs), exp)

    def test_iter_contiguous_iterable_slices(self):
        def spaced_out():
            yield slice(0, 4)
            yield slice(8, 12)

        def contiguous():
            yield slice(0, 4)
            yield slice(4, 8)
            yield slice(12, 16)

        s = Sequence("0123456789abcdef")
        for c in (lambda x: x, list, tuple, lambda x: np.array(tuple(x)),
                  lambda x: pd.Series(tuple(x))):
            exp = [Sequence("0123"), Sequence("89ab")]
            obs = s.iter_contiguous(c(spaced_out()))
            self.assertEqual(list(obs), exp)

            exp = [Sequence("01234567"), Sequence("cdef")]
            obs = s.iter_contiguous(c(contiguous()))
            self.assertEqual(list(obs), exp)

    def test_iter_contiguous_with_max_length(self):
        s = Sequence("0123456789abcdef")
        for c in list, tuple, np.array, pd.Series:
            exp = [Sequence("234"), Sequence("678"), Sequence("abc")]
            obs = s.iter_contiguous(c([True, False, True, True] * 4),
                                    min_length=3)
            self.assertEqual(list(obs), exp)

            exp = [Sequence("0"), Sequence("234"), Sequence("678"),
                   Sequence("abc"), Sequence("ef")]
            obs1 = list(s.iter_contiguous(c([True, False, True, True] * 4),
                                          min_length=1))

            obs2 = list(s.iter_contiguous(c([True, False, True, True] * 4)))
            self.assertEqual(obs1, obs2)
            self.assertEqual(obs1, exp)

    def test_iter_contiguous_with_invert(self):
        def spaced_out():
            yield slice(0, 4)
            yield slice(8, 12)

        def contiguous():
            yield slice(0, 4)
            yield slice(4, 8)
            yield slice(12, 16)

        s = Sequence("0123456789abcdef")
        for c in (lambda x: x, list, tuple, lambda x: np.array(tuple(x)),
                  lambda x: pd.Series(tuple(x))):
            exp = [Sequence("4567"), Sequence("cdef")]
            obs = s.iter_contiguous(c(spaced_out()), invert=True)
            self.assertEqual(list(obs), exp)

            exp = [Sequence("89ab")]
            obs = s.iter_contiguous(c(contiguous()), invert=True)
            self.assertEqual(list(obs), exp)

    def test_copy_without_metadata(self):
        # shallow vs deep copy with sequence only should be equivalent
        for copy_method in copy.copy, copy.deepcopy:
            seq = Sequence('ACGT')
            seq_copy = copy_method(seq)

            self.assertEqual(seq_copy, seq)
            self.assertIsNot(seq_copy, seq)
            self.assertIsNot(seq_copy._bytes, seq._bytes)

    def test_copy_with_metadata_shallow(self):
        seq = Sequence('ACGT', metadata={'foo': [1]},
                       positional_metadata={'bar': [[], [], [], []],
                                            'baz': [42, 42, 42, 42]})
        seq.interval_metadata.add([(0, 3)], metadata={'gene': ['sagA']})

        seq_copy = copy.copy(seq)

        self.assertEqual(seq_copy, seq)
        self.assertIsNot(seq_copy, seq)
        self.assertIsNot(seq_copy._bytes, seq._bytes)
        self.assertIsNot(seq_copy._metadata, seq._metadata)
        self.assertIsNot(seq_copy._positional_metadata,
                         seq._positional_metadata)
        self.assertIsNot(seq_copy._positional_metadata.values,
                         seq._positional_metadata.values)
        self.assertIs(seq_copy._metadata['foo'], seq._metadata['foo'])
        self.assertIs(seq_copy._positional_metadata.loc[0, 'bar'],
                      seq._positional_metadata.loc[0, 'bar'])
        self.assertIsNot(seq_copy.interval_metadata, seq.interval_metadata)
        self.assertIsNot(seq_copy.interval_metadata._intervals[0],
                         seq.interval_metadata._intervals[0])
        self.assertIsNot(seq_copy.interval_metadata._intervals[0].metadata,
                         seq.interval_metadata._intervals[0].metadata)
        self.assertIs(
            seq_copy.interval_metadata._intervals[0].metadata['gene'],
            seq.interval_metadata._intervals[0].metadata['gene'])
        seq_copy.metadata['foo'].append(2)
        seq_copy.metadata['foo2'] = 42

        self.assertEqual(seq_copy.metadata, {'foo': [1, 2], 'foo2': 42})
        self.assertEqual(seq.metadata, {'foo': [1, 2]})

        seq_copy.positional_metadata.loc[0, 'bar'].append(1)
        seq_copy.positional_metadata.loc[0, 'baz'] = 43

        assert_data_frame_almost_equal(
            seq_copy.positional_metadata,
            pd.DataFrame({'bar': [[1], [], [], []],
                          'baz': [43, 42, 42, 42]}))
        assert_data_frame_almost_equal(
            seq.positional_metadata,
            pd.DataFrame({'bar': [[1], [], [], []],
                          'baz': [42, 42, 42, 42]}))

    def test_copy_with_metadata_deep(self):
        seq = Sequence('ACGT', metadata={'foo': [1]},
                       positional_metadata={'bar': [[], [], [], []],
                                            'baz': [42, 42, 42, 42]})
        seq.interval_metadata.add([(0, 3)], metadata={'gene': ['sagA']})
        seq_copy = copy.deepcopy(seq)

        self.assertEqual(seq_copy, seq)
        self.assertIsNot(seq_copy, seq)
        self.assertIsNot(seq_copy._bytes, seq._bytes)
        self.assertIsNot(seq_copy._metadata, seq._metadata)
        self.assertIsNot(seq_copy._positional_metadata,
                         seq._positional_metadata)
        self.assertIsNot(seq_copy._positional_metadata.values,
                         seq._positional_metadata.values)
        self.assertIsNot(seq_copy._metadata['foo'], seq._metadata['foo'])
        self.assertIsNot(seq_copy._positional_metadata.loc[0, 'bar'],
                         seq._positional_metadata.loc[0, 'bar'])
        self.assertIsNot(seq_copy.interval_metadata, seq.interval_metadata)
        self.assertIsNot(seq_copy.interval_metadata._intervals[0],
                         seq.interval_metadata._intervals[0])
        self.assertIsNot(seq_copy.interval_metadata._intervals[0].metadata,
                         seq.interval_metadata._intervals[0].metadata)
        self.assertIsNot(
            seq_copy.interval_metadata._intervals[0].metadata['gene'],
            seq.interval_metadata._intervals[0].metadata['gene'])
        seq_copy.metadata['foo'].append(2)
        seq_copy.metadata['foo2'] = 42

        self.assertEqual(seq_copy.metadata, {'foo': [1, 2], 'foo2': 42})
        self.assertEqual(seq.metadata, {'foo': [1]})

        seq_copy.positional_metadata.loc[0, 'bar'].append(1)
        seq_copy.positional_metadata.loc[0, 'baz'] = 43

        assert_data_frame_almost_equal(
            seq_copy.positional_metadata,
            pd.DataFrame({'bar': [[1], [], [], []],
                          'baz': [43, 42, 42, 42]}))
        assert_data_frame_almost_equal(
            seq.positional_metadata,
            pd.DataFrame({'bar': [[], [], [], []],
                          'baz': [42, 42, 42, 42]}))

    def test_copy_preserves_read_only_flag_on_bytes(self):
        seq = Sequence('ACGT')
        seq_copy = copy.copy(seq)

        with self.assertRaises(ValueError):
            seq_copy._bytes[0] = 'B'

    def test_deepcopy_memo_is_respected(self):
        # basic test to ensure deepcopy's memo is passed through to recursive
        # deepcopy calls
        seq = Sequence('ACGT', metadata={'foo': 'bar'})
        memo = {}
        copy.deepcopy(seq, memo)
        self.assertGreater(len(memo), 2)

    def test_munge_to_index_array_valid_index_array(self):
        s = Sequence('123456')

        for c in list, tuple, np.array, pd.Series:
            exp = np.array([1, 2, 3], dtype=int)
            obs = s._munge_to_index_array(c([1, 2, 3]))
            npt.assert_equal(obs, exp)

            exp = np.array([1, 3, 5], dtype=int)
            obs = s._munge_to_index_array(c([1, 3, 5]))
            npt.assert_equal(obs, exp)

    def test_munge_to_index_array_invalid_index_array(self):
        s = Sequence("12345678")
        for c in list, tuple, np.array, pd.Series:
            with self.assertRaises(ValueError):
                s._munge_to_index_array(c([3, 2, 1]))

            with self.assertRaises(ValueError):
                s._munge_to_index_array(c([5, 6, 7, 2]))

            with self.assertRaises(ValueError):
                s._munge_to_index_array(c([0, 1, 2, 1]))

    def test_munge_to_index_array_valid_bool_array(self):
        s = Sequence('123456')

        for c in list, tuple, np.array, pd.Series:
            exp = np.array([2, 3, 5], dtype=int)
            obs = s._munge_to_index_array(
                c([False, False, True, True, False, True]))
            npt.assert_equal(obs, exp)

            exp = np.array([], dtype=int)
            obs = s._munge_to_index_array(
                c([False] * 6))
            npt.assert_equal(obs, exp)

            exp = np.arange(6)
            obs = s._munge_to_index_array(
                c([True] * 6))
            npt.assert_equal(obs, exp)

    def test_munge_to_index_array_invalid_bool_array(self):
        s = Sequence('123456')

        for c in (list, tuple, lambda x: np.array(x, dtype=bool),
                  lambda x: pd.Series(x, dtype=bool)):

            with self.assertRaises(ValueError):
                s._munge_to_index_array(c([]))

            with self.assertRaises(ValueError):
                s._munge_to_index_array(c([True]))

            with self.assertRaises(ValueError):
                s._munge_to_index_array(c([True] * 10))

    def test_munge_to_index_array_valid_iterable(self):
        s = Sequence('')

        def slices_only():
            return (slice(i, i+1) for i in range(0, 10, 2))

        def mixed():
            return (slice(i, i+1) if i % 2 == 0 else i for i in range(10))

        def unthinkable():
            for i in range(10):
                if i % 3 == 0:
                    yield slice(i, i+1)
                elif i % 3 == 1:
                    yield i
                else:
                    yield np.array([i], dtype=int)
        for c in (lambda x: x, list, tuple, lambda x: np.array(tuple(x)),
                  lambda x: pd.Series(tuple(x))):
            exp = np.arange(10, dtype=int)
            obs = s._munge_to_index_array(c(mixed()))
            npt.assert_equal(obs, exp)

            exp = np.arange(10, dtype=int)
            obs = s._munge_to_index_array(c(unthinkable()))
            npt.assert_equal(obs, exp)

            exp = np.arange(10, step=2, dtype=int)
            obs = s._munge_to_index_array(c(slices_only()))
            npt.assert_equal(obs, exp)

    def test_munge_to_index_array_invalid_iterable(self):
        s = Sequence('')

        def bad1():
            yield "r"
            yield [1, 2, 3]

        def bad2():
            yield 1
            yield 'str'

        def bad3():
            yield False
            yield True
            yield 2

        def bad4():
            yield np.array([False, True])
            yield slice(2, 5)

        for c in (lambda x: x, list, tuple, lambda x: np.array(tuple(x)),
                  lambda x: pd.Series(tuple(x))):

            with self.assertRaises(TypeError):
                s._munge_to_index_array(bad1())

            with self.assertRaises(TypeError):
                s._munge_to_index_array(bad2())

            with self.assertRaises(TypeError):
                s._munge_to_index_array(bad3())

            with self.assertRaises(TypeError):
                s._munge_to_index_array(bad4())

    def test_munge_to_index_array_valid_string(self):
        seq = Sequence('ACGTACGT',
                       positional_metadata={'introns': [False, True, True,
                                                        False, False, True,
                                                        False, False]})
        npt.assert_equal(np.array([1, 2, 5]),
                         seq._munge_to_index_array('introns'))

        seq.positional_metadata['exons'] = ~seq.positional_metadata['introns']
        npt.assert_equal(np.array([0, 3, 4, 6, 7]),
                         seq._munge_to_index_array('exons'))

    def test_munge_to_index_array_invalid_string(self):
        seq_str = 'ACGT'
        seq = Sequence(seq_str,
                       positional_metadata={'quality': range(len(seq_str))})

        with self.assertRaisesRegex(ValueError,
                                    r"No positional metadata associated with "
                                    "key 'introns'"):
            seq._munge_to_index_array('introns')

        with self.assertRaisesRegex(TypeError,
                                    r"Column 'quality' in positional metadata "
                                    "does not correspond to a boolean "
                                    "vector"):
            seq._munge_to_index_array('quality')

    def test_munge_to_bytestring_return_bytes(self):
        seq = Sequence('')
        m = 'dummy_method'
        str_inputs = ('', 'a', 'acgt')
        unicode_inputs = ('', 'a', 'acgt')
        byte_inputs = (b'', b'a', b'acgt')
        seq_inputs = (Sequence(''), Sequence('a'), Sequence('acgt'))
        all_inputs = str_inputs + unicode_inputs + byte_inputs + seq_inputs
        all_expected = [b'', b'a', b'acgt'] * 4

        for input_, expected in zip(all_inputs, all_expected):
            observed = seq._munge_to_bytestring(input_, m)
            self.assertEqual(observed, expected)
            self.assertIs(type(observed), bytes)

    def test_munge_to_bytestring_unicode_out_of_ascii_range(self):
        seq = Sequence('')
        all_inputs = ('\x80', 'abc\x80', '\x80abc')
        for input_ in all_inputs:
            with self.assertRaisesRegex(UnicodeEncodeError,
                                        r"'ascii' codec can't encode character"
                                        r".*in position.*: ordinal not in"
                                        r" range\(128\)"):
                seq._munge_to_bytestring(input_, 'dummy_method')


class TestDistance(TestSequenceBase):
    def test_mungeable_inputs_to_sequence(self):
        def metric(a, b):
            self.assertEqual(a, Sequence("abcdef"))
            self.assertEqual(b, Sequence("12bcef"))
            return 42.0

        for constructor in self.sequence_kinds:
            seq1 = Sequence("abcdef")
            seq2 = constructor("12bcef")

            distance = seq1.distance(seq2, metric=metric)

            self.assertEqual(distance, 42.0)

    def test_mungeable_inputs_to_sequence_subclass(self):
        def metric(a, b):
            self.assertEqual(a, SequenceSubclass("abcdef"))
            self.assertEqual(b, SequenceSubclass("12bcef"))
            return -42.0

        sequence_kinds = frozenset([
            str, SequenceSubclass,
            lambda s: np.frombuffer(s.encode('ascii'), dtype='|S1'),
            lambda s: np.frombuffer(s.encode('ascii'), dtype=np.uint8)])

        for constructor in sequence_kinds:
            seq1 = SequenceSubclass("abcdef")
            seq2 = constructor("12bcef")

            distance = seq1.distance(seq2, metric=metric)

            self.assertEqual(distance, -42.0)

    def test_sequence_type_mismatch(self):
        seq1 = SequenceSubclass("abcdef")
        seq2 = Sequence("12bcef")

        with self.assertRaisesRegex(TypeError,
                                    r'SequenceSubclass.*Sequence.*`distance`'):
            seq1.distance(seq2)

        with self.assertRaisesRegex(TypeError,
                                    r'Sequence.*SequenceSubclass.*`distance`'):
            seq2.distance(seq1)

    def test_munging_invalid_characters_to_self_type(self):
        with self.assertRaisesRegex(ValueError, r'Invalid characters.*X'):
            DNA("ACGT").distance("WXYZ")

    def test_munging_invalid_type_to_self_type(self):
        with self.assertRaises(TypeError):
            Sequence("ACGT").distance(42)

    def test_return_type_coercion(self):
        def metric(a, b):
            return 42

        distance = Sequence('abc').distance('cba', metric=metric)

        self.assertIsInstance(distance, float)

    def test_invalid_return_type(self):
        def metric(a, b):
            return 'too far'

        with self.assertRaisesRegex(ValueError, r'string.*float'):
            Sequence('abc').distance('cba', metric=metric)

    def test_arbitrary_metric(self):
        def metric(x, y):
            return len(x) ** 2 + len(y) ** 2

        seq1 = Sequence("12345678")
        seq2 = Sequence("1234")

        distance = seq1.distance(seq2, metric=metric)

        self.assertEqual(distance, 80.0)

    def test_scipy_hamming_metric_with_metadata(self):
        # test for #1254
        seqs1 = [
            Sequence("ACGT"),
            Sequence("ACGT", metadata={'id': 'abc'}),
            Sequence("ACGT", positional_metadata={'qual': range(4)})
        ]
        seqs2 = [
            Sequence("AAAA"),
            Sequence("AAAA", metadata={'id': 'def'}),
            Sequence("AAAA", positional_metadata={'qual': range(4, 8)})
        ]

        for seqs in seqs1, seqs2:
            for seq1, seq2 in itertools.product(seqs, repeat=2):
                distance = seq1.distance(seq2,
                                         metric=scipy.spatial.distance.hamming)
                self.assertEqual(distance, 0.0)

        for seq1, seq2 in itertools.product(seqs1, seqs2):
            distance = seq1.distance(seq2,
                                     metric=scipy.spatial.distance.hamming)
            self.assertEqual(distance, 0.75)

    def test_default_metric_with_metadata(self):
        # test for #1254
        seqs1 = [
            Sequence("ACGT"),
            Sequence("ACGT", metadata={'id': 'abc'}),
            Sequence("ACGT", positional_metadata={'qual': range(4)})
        ]
        seqs2 = [
            Sequence("AAAA"),
            Sequence("AAAA", metadata={'id': 'def'}),
            Sequence("AAAA", positional_metadata={'qual': range(4, 8)})
        ]

        for seqs in seqs1, seqs2:
            for seq1, seq2 in itertools.product(seqs, repeat=2):
                distance = seq1.distance(seq2)
                self.assertEqual(distance, 0.0)

        for seq1, seq2 in itertools.product(seqs1, seqs2):
            distance = seq1.distance(seq2)
            self.assertEqual(distance, 0.75)

    def test_default_metric_matches_hamming(self):
        seq1 = Sequence("abcdef")
        seq2 = Sequence("12bcef")
        seq_wrong = Sequence("abcdefghijklmnop")

        distance1 = seq1.distance(seq2)
        distance2 = skbio.sequence.distance.hamming(seq1, seq2)

        self.assertEqual(distance1, distance2)

        with self.assertRaises(ValueError):
            seq1.distance(seq_wrong)

        with self.assertRaises(ValueError):
            seq_wrong.distance(seq1)


# NOTE: this must be a *separate* class for doctests only (no unit tests). nose
# will not run the unit tests otherwise
#
# these doctests exercise the correct formatting of Sequence's repr in a
# variety of situations. they are more extensive than the unit tests above
# (TestSequence.test_repr) but cannot be relied upon for coverage (the unit
# tests take care of this)
class SequenceReprDoctests:
    r""">>> import pandas as pd
    >>> from skbio import Sequence

    Empty (minimal) sequence:

    >>> Sequence('')
    Sequence
    -------------
    Stats:
        length: 0
    -------------

    Single character sequence:

    >>> Sequence('G')
    Sequence
    -------------
    Stats:
        length: 1
    -------------
    0 G

    Multicharacter sequence:

    >>> Sequence('ACGT')
    Sequence
    -------------
    Stats:
        length: 4
    -------------
    0 ACGT

    Full single line:

    >>> Sequence('A' * 60)
    Sequence
    -------------------------------------------------------------------
    Stats:
        length: 60
    -------------------------------------------------------------------
    0 AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA

    Full single line with 1 character overflow:

    >>> Sequence('A' * 61)
    Sequence
    --------------------------------------------------------------------
    Stats:
        length: 61
    --------------------------------------------------------------------
    0  AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    60 A

    Two full lines:

    >>> Sequence('T' * 120)
    Sequence
    --------------------------------------------------------------------
    Stats:
        length: 120
    --------------------------------------------------------------------
    0  TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT
    60 TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT

    Two full lines with 1 character overflow:

    >>> Sequence('T' * 121)
    Sequence
    ---------------------------------------------------------------------
    Stats:
        length: 121
    ---------------------------------------------------------------------
    0   TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT
    60  TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT
    120 T

    Five full lines (maximum amount of information):

    >>> Sequence('A' * 300)
    Sequence
    ---------------------------------------------------------------------
    Stats:
        length: 300
    ---------------------------------------------------------------------
    0   AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    60  AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    120 AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    180 AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    240 AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA

    Six lines starts "summarized" output:

    >>> Sequence('A' * 301)
    Sequence
    ---------------------------------------------------------------------
    Stats:
        length: 301
    ---------------------------------------------------------------------
    0   AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    60  AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    ...
    240 AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    300 A

    A naive algorithm would assume the width of the first column (noting
    position) based on the sequence's length alone. This can be off by one if
    the last position (in the last line) has a shorter width than the width
    calculated from the sequence's length. This test case ensures that only a
    single space is inserted between position 99960 and the first sequence
    chunk:

    >>> Sequence('A' * 100000)
    Sequence
    -----------------------------------------------------------------------
    Stats:
        length: 100000
    -----------------------------------------------------------------------
    0     AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    60    AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    ...
    99900 AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    99960 AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA

    The largest sequence that can be displayed using six chunks per line:

    >>> Sequence('A' * 100020)
    Sequence
    -----------------------------------------------------------------------
    Stats:
        length: 100020
    -----------------------------------------------------------------------
    0     AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    60    AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    ...
    99900 AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    99960 AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA

    A single character longer than the previous sequence causes the optimal
    number of chunks per line to be 5:

    >>> Sequence('A' * 100021)
    Sequence
    -------------------------------------------------------------
    Stats:
        length: 100021
    -------------------------------------------------------------
    0      AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    50     AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    ...
    99950  AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA AAAAAAAAAA
    100000 AAAAAAAAAA AAAAAAAAAA A

    Wide range of characters (locale-independent):

    >>> import string
    >>> Sequence((string.ascii_letters + string.punctuation + string.digits +
    ...          'a space') * 567)
    Sequence
    -----------------------------------------------------------------------
    Stats:
        length: 57267
    -----------------------------------------------------------------------
    0     abcdefghij klmnopqrst uvwxyzABCD EFGHIJKLMN OPQRSTUVWX YZ!"#$%&'(
    60    )*+,-./:;< =>?@[\]^_` {|}~012345 6789a spac eabcdefghi jklmnopqrs
    ...
    57180 opqrstuvwx yzABCDEFGH IJKLMNOPQR STUVWXYZ!" #$%&'()*+, -./:;<=>?@
    57240 [\]^_`{|}~ 0123456789 a space

    Supply horrendous metadata, positional, and interval metadata to
    exercise a variety of metadata formatting cases and rules. Sorting
    should be by type, then by value within each type (Python 3
    doesn't allow sorting of mixed types):

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
    >>> interval_metadata = IntervalMetadata(4)
    >>> _ = interval_metadata.add([(0, 2), (1, 3)],
    ...                           [(False, True), (False, False)],
    ...                           {'gene': 'p53'})
    >>> _ = interval_metadata.add([(1, 4)])
    >>> Sequence('ACGT', metadata=metadata,
    ...          positional_metadata=positional_metadata,
    ...          interval_metadata=interval_metadata)
    Sequence
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
    Interval metadata:
        2 interval features
    Stats:
        length: 4
    -----------------------------------------------------------------------
    0 ACGT

    """
    pass


if __name__ == "__main__":
    main()
