# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.standard_library import hooks
from six import string_types

import re
from types import GeneratorType
from collections import Counter, defaultdict, Hashable
from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
from scipy.spatial.distance import euclidean

from skbio import (
    Sequence, DNA, RNA)

with hooks():
    from itertools import zip_longest

class SequenceSubclass(Sequence):
    """Used for testing purposes."""
    pass

class SequenceTests(TestCase):
    def setUp(self):
        # DONT DELETE THIS!! JUST REMOVE THIS COMMENT
        self.sequence_kinds = frozenset([
            str, Sequence, lambda s: np.fromstring(s, dtype='|S1'),
            lambda s: np.fromstring(s, dtype=np.uint8)])
        # BELOW THIS CAN BE REMOVED
        self.b1 = Sequence('GATTACA', quality=range(7))
        self.b2 = Sequence(
            'ACCGGTACC', id="test-seq-2",
            description="A test sequence")
        self.b3 = Sequence(
            'GREG', id="test-seq-3", description="A protein sequence")
        self.b4 = Sequence(
            'PRTEIN', id="test-seq-4")
        self.b5 = Sequence(
            'LLPRTEIN', description="some description")
        self.b6 = Sequence('ACGTACGTACGT')
        self.b7 = Sequence('..--..', quality=range(6))
        self.b8 = Sequence('HE..--..LLO', id='hello',
                           description='gapped hello',
                           quality=range(11))


    def test_init_default_parameters(self):
        seq = Sequence('.ABC123xyz-')

        npt.assert_equal(seq.sequence, np.array('.ABC123xyz-', dtype='c'))
        self.assertEqual(seq.id, "")
        self.assertEqual(seq.description, "")
        self.assertIsNone(seq.quality)

    def test_init_nondefault_parameters(self):
        seq = Sequence('.ABC123xyz-', id='foo', description='bar baz',
                       quality=range(11))

        npt.assert_equal(seq.sequence, np.array('.ABC123xyz-', dtype='c'))
        self.assertEqual(seq.id, 'foo')
        self.assertEqual(seq.description, 'bar baz')
        npt.assert_equal(seq.quality, np.array(range(11), dtype='int'))

    def test_init_empty_sequence(self):
        # Test constructing an empty sequence using each supported input type.
        for s in (b'',  # bytes
                  u'',  # unicode
                  np.array('', dtype='c'),  # char vector
                  np.fromstring('', dtype=np.uint8),  # byte vec
                  Sequence('')):  # another Sequence object
            seq = Sequence(s)

            self.assertIsInstance(seq.sequence, np.ndarray)
            self.assertEqual(seq.sequence.dtype, '|S1')
            self.assertEqual(seq.sequence.shape, (0,))
            npt.assert_equal(seq.sequence, np.array('', dtype='c'))

    def test_init_single_character_sequence(self):
        for s in (b'A',
                  u'A',
                  np.array('A', dtype='c'),
                  np.fromstring('A', dtype=np.uint8),
                  Sequence('A')):
            seq = Sequence(s)

            self.assertIsInstance(seq.sequence, np.ndarray)
            self.assertEqual(seq.sequence.dtype, '|S1')
            self.assertEqual(seq.sequence.shape, (1,))
            npt.assert_equal(seq.sequence, np.array('A', dtype='c'))

    def test_init_multiple_character_sequence(self):
        for s in (b'.ABC\t123  xyz-',
                  u'.ABC\t123  xyz-',
                  np.array('.ABC\t123  xyz-', dtype='c'),
                  np.fromstring('.ABC\t123  xyz-', dtype=np.uint8),
                  Sequence('.ABC\t123  xyz-')):
            seq = Sequence(s)

            self.assertIsInstance(seq.sequence, np.ndarray)
            self.assertEqual(seq.sequence.dtype, '|S1')
            self.assertEqual(seq.sequence.shape, (14,))
            npt.assert_equal(seq.sequence,
                             np.array('.ABC\t123  xyz-', dtype='c'))

    def test_init_from_sequence_object(self):
        # We're testing this in its simplest form in other tests. This test
        # exercises more complicated cases of building a sequence from another
        # sequence.

        # just the sequence, no other metadata
        seq = Sequence('ACGT')
        self.assertEqual(Sequence(seq), seq)

        # sequence with metadata should have everything propagated
        seq = Sequence('ACGT', id='foo', description='bar baz',
                       quality=range(4))
        self.assertEqual(Sequence(seq), seq)

        # should be able to override metadata
        self.assertEqual(
            Sequence(seq, id='abc', description='123', quality=[42] * 4),
            Sequence('ACGT', id='abc', description='123', quality=[42] * 4))

        # subclasses work too
        seq = SequenceSubclass('ACGT', id='foo', description='bar baz',
                               quality=range(4))
        self.assertEqual(
            Sequence(seq),
            Sequence('ACGT', id='foo', description='bar baz',
                     quality=range(4)))

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

    def test_init_empty_id(self):
        seq = Sequence('', id='')

        self.assertIsInstance(seq.id, string_types)
        self.assertEqual(seq.id, '')

    def test_init_single_character_id(self):
        seq = Sequence('', id='z')

        self.assertIsInstance(seq.id, string_types)
        self.assertEqual(seq.id, 'z')

    def test_init_multiple_character_id(self):
        seq = Sequence('', id='\nabc\tdef  G123')

        self.assertIsInstance(seq.id, string_types)
        self.assertEqual(seq.id, '\nabc\tdef  G123')

    def test_init_empty_description(self):
        seq = Sequence('', description='')

        self.assertIsInstance(seq.description, string_types)
        self.assertEqual(seq.description, '')

    def test_init_single_character_description(self):
        seq = Sequence('', description='z')

        self.assertIsInstance(seq.description, string_types)
        self.assertEqual(seq.description, 'z')

    def test_init_multiple_character_description(self):
        seq = Sequence('', description='\nabc\tdef  G123')

        self.assertIsInstance(seq.description, string_types)
        self.assertEqual(seq.description, '\nabc\tdef  G123')

    def test_init_empty_quality(self):
        for q in ([], (), np.array([])):
            seq = Sequence('', quality=q)

            self.assertIsInstance(seq.quality, np.ndarray)
            self.assertEqual(seq.quality.dtype, np.int)
            self.assertEqual(seq.quality.shape, (0,))
            npt.assert_equal(seq.quality, np.array([]))

    def test_init_single_quality_score(self):
        for q in (2, [2], (2,), np.array([2])):
            seq = Sequence('G', quality=q)

            self.assertIsInstance(seq.quality, np.ndarray)
            self.assertEqual(seq.quality.dtype, np.int)
            self.assertEqual(seq.quality.shape, (1,))
            npt.assert_equal(seq.quality, np.array([2]))

    def test_init_multiple_quality_scores(self):
        for q in ([0, 42, 42, 1, 0, 8, 100, 0, 0],
                  (0, 42, 42, 1, 0, 8, 100, 0, 0),
                  np.array([0, 42, 42, 1, 0, 8, 100, 0, 0])):
            seq = Sequence('G' * 9, quality=q)

            self.assertIsInstance(seq.quality, np.ndarray)
            self.assertEqual(seq.quality.dtype, np.int)
            self.assertEqual(seq.quality.shape, (9,))
            npt.assert_equal(seq.quality,
                             np.array([0, 42, 42, 1, 0, 8, 100, 0, 0]))

    def test_init_no_copy_of_quality(self):
        qual = np.array([22, 22, 1])
        seq = Sequence('ACA', quality=qual)

        self.assertIs(seq.quality, qual)

        with self.assertRaises(ValueError):
            qual[1] = 42

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
        with self.assertRaisesRegexp(TypeError, 'tuple'):
            Sequence(('a', 'b', 'c'))
        with self.assertRaisesRegexp(TypeError, 'list'):
            Sequence(['a', 'b', 'c'])
        with self.assertRaisesRegexp(TypeError, 'set'):
            Sequence({'a', 'b', 'c'})
        with self.assertRaisesRegexp(TypeError, 'dict'):
            Sequence({'a': 42, 'b': 43, 'c': 44})
        with self.assertRaisesRegexp(TypeError, 'int'):
            Sequence(42)
        with self.assertRaisesRegexp(TypeError, 'float'):
            Sequence(4.2)
        with self.assertRaisesRegexp(TypeError, 'Foo'):
            class Foo(object):
                pass
            Sequence(Foo())

        # out of ASCII range
        with self.assertRaises(UnicodeEncodeError):
            Sequence(u'abc\u1F30')

    def test_init_invalid_id(self):
        with self.assertRaises(TypeError):
            Sequence('abc', id=('f', 'o', 'o'))

    def test_init_invalid_description(self):
        with self.assertRaises(TypeError):
            Sequence('abc', description=('f', 'o', 'o'))

    def test_init_invalid_quality(self):
        # invalid dtype
        with self.assertRaises(TypeError):
            Sequence('ACGT', quality=[2, 3, 4.1, 5])
        with self.assertRaises(TypeError):
            Sequence('ACGT', quality=[2, np.nan, 4, 5])

        # wrong number of dimensions
        with self.assertRaisesRegexp(ValueError, '2.*1-D'):
            Sequence('ACGT', quality=[[2, 3], [4, 5]])

        # wrong number of elements
        with self.assertRaisesRegexp(ValueError, '\(3\).*\(4\)'):
            Sequence('ACGT', quality=[2, 3, 4])

        # negatives
        with self.assertRaisesRegexp(ValueError,
                                     'Quality scores.*greater than.*zero'):
            Sequence('ACGT', quality=[2, 3, -1, 4])

    def test_sequence_property(self):
        # Property tests are only concerned with testing the interface
        # provided by the property: that it can be accessed, can't be
        # reassigned or mutated in place, and that the correct type is
        # returned. More extensive testing of border cases (e.g., different
        # sequence lengths or input types, odd characters, etc.) are performed
        # in Sequence.__init__ tests.

        seq = Sequence('ACGT')

        # should get back a numpy.ndarray of '|S1' dtype
        self.assertIsInstance(seq.sequence, np.ndarray)
        self.assertEqual(seq.sequence.dtype, '|S1')
        npt.assert_equal(seq.sequence, np.array('ACGT', dtype='c'))

        # test that we can't mutate the property
        with self.assertRaises(ValueError):
            seq.sequence[1] = 'A'

        # test that we can't set the property
        with self.assertRaises(AttributeError):
            seq.sequence = np.array("GGGG", dtype='c')

    def test_id_property(self):
        seq = Sequence('', id='foo')

        self.assertIsInstance(seq.id, string_types)
        self.assertEqual(seq.id, 'foo')

        with self.assertRaises(TypeError):
            seq.id[1] = 42

        with self.assertRaises(AttributeError):
            seq.id = 'bar'

    def test_description_property(self):
        seq = Sequence('', description='foo')

        self.assertIsInstance(seq.description, string_types)
        self.assertEqual(seq.description, 'foo')

        with self.assertRaises(TypeError):
            seq.description[1] = 42

        with self.assertRaises(AttributeError):
            seq.description = 'bar'

    def test_quality_property(self):
        seq = Sequence('ACA', quality=[22, 22, 0])

        self.assertIsInstance(seq.quality, np.ndarray)
        self.assertEqual(seq.quality.dtype, np.int)
        npt.assert_equal(seq.quality, np.array([22, 22, 0]))

        with self.assertRaises(ValueError):
            seq.quality[1] = 42

        with self.assertRaises(AttributeError):
            seq.quality = [22, 22, 42]

    def test_has_quality(self):
        seq = Sequence('')
        self.assertFalse(seq._has_quality())

        seq = Sequence('', quality=[])
        self.assertTrue(seq._has_quality())

        seq = Sequence('ACA', quality=(5, 4, 67))
        self.assertTrue(seq._has_quality())

        seq = Sequence('ACA')
        self.assertFalse(seq._has_quality())

    def test_eq_and_ne(self):
        seq_a = Sequence("A")
        seq_b = Sequence("B")

        self.assertTrue(seq_a == seq_a)
        self.assertTrue(Sequence("a") == Sequence("a"))
        self.assertTrue(Sequence("a", id='b') == Sequence("a", id='b'))
        self.assertTrue(Sequence("a", id='b', description='c') ==
                        Sequence("a", id='b', description='c'))
        self.assertTrue(Sequence("a", id='b', description='c', quality=[1]) ==
                        Sequence("a", id='b', description='c', quality=[1]))

        self.assertTrue(seq_a != seq_b)
        self.assertTrue(SequenceSubclass("a") != Sequence("a"))
        self.assertTrue(Sequence("a") != Sequence("b"))
        self.assertTrue(Sequence("a") != Sequence("a", id='b'))
        self.assertTrue(Sequence("a", id='c') !=
                        Sequence("a", id='c', description='t'))
        self.assertTrue(Sequence("a", quality=[1]) != Sequence("a"))
        self.assertTrue(Sequence("a", quality=[2]) !=
                        Sequence("a", quality=[1]))
        self.assertTrue(Sequence("c", quality=[3]) !=
                        Sequence("b", quality=[3]))
        self.assertTrue(Sequence("a", id='b') != Sequence("c", id='b'))

    def test_getitem_gives_new_sequence(self):
        seq = Sequence("Sequence string !1@2#3?.,")
        self.assertFalse(seq is seq[:])

    def test_getitem_with_int_has_qual(self):
        s = "Sequence string !1@2#3?.,"
        length = len(s)
        seq = Sequence(s, id='id', description='dsc',
                       quality=np.arange(length))

        eseq = Sequence("S", id='id', description='dsc', quality=np.array([0]))
        self.assertEqual(seq[0], eseq)

        eseq = Sequence(",", id='id', description='dsc',
                        quality=np.array([len(seq) - 1]))
        self.assertEqual(seq[len(seq) - 1], eseq)

        eseq = Sequence("t", id='id', description='dsc',
                        quality=[10])
        self.assertEqual(seq[10], eseq)

    def test_getitem_with_int_no_qual(self):
        seq = Sequence("Sequence string !1@2#3?.,", id='id2',
                       description='no_qual')

        eseq = Sequence("t", id='id2', description='no_qual')
        self.assertEqual(seq[10], eseq)

    def test_getitem_with_slice_has_qual(self):
        s = "0123456789abcdef"
        length = len(s)
        seq = Sequence(s, id='id3', description="dsc3",
                       quality=np.arange(length))

        eseq = Sequence("012", id='id3', description="dsc3",
                        quality=np.arange(3))
        self.assertEquals(seq[0:3], eseq)
        self.assertEquals(seq[:3], eseq)
        self.assertEquals(seq[:3:1], eseq)

        eseq = Sequence("def", id='id3', description="dsc3",
                        quality=[13, 14, 15])
        self.assertEquals(seq[-3:], eseq)
        self.assertEquals(seq[-3::1], eseq)

        eseq = Sequence("02468ace", id='id3', description='dsc3',
                        quality=[0, 2, 4, 6, 8, 10, 12, 14])
        self.assertEquals(seq[0:length:2], eseq)
        self.assertEquals(seq[::2], eseq)

        eseq = Sequence(s[::-1], id='id3', description='dsc3',
                        quality=np.arange(length)[::-1])
        self.assertEquals(seq[length::-1], eseq)
        self.assertEquals(seq[::-1], eseq)

        eseq = Sequence('fdb97531', id='id3', description='dsc3',
                        quality=[15, 13, 11, 9, 7, 5, 3, 1])
        self.assertEquals(seq[length::-2], eseq)
        self.assertEquals(seq[::-2], eseq)

        self.assertEquals(seq[0:500:], seq)

        eseq = Sequence('', id='id3', description='dsc3',
                        quality=[])
        self.assertEquals(seq[length:0], eseq)
        self.assertEquals(seq[-length:0], eseq)
        self.assertEquals(seq[1:0], eseq)

        eseq = Sequence("0", id='id3', description='dsc3',
                        quality=[0])
        self.assertEquals(seq[0:1], eseq)
        self.assertEquals(seq[0:1:1], eseq)
        self.assertEquals(seq[-length::-1], eseq)

    def test_getitem_with_slice_no_qual(self):
        s = "0123456789abcdef"
        length = len(s)
        seq = Sequence(s, id='id4', description="no_qual4")

        eseq = Sequence("02468ace", id='id4', description='no_qual4')
        self.assertEquals(seq[0:length:2], eseq)
        self.assertEquals(seq[::2], eseq)

    def test_getitem_with_tuple_of_mixed_with_qual(self):
        s = "0123456789abcdef"
        length = len(s)
        seq = Sequence(s, id='id5', description="dsc5",
                       quality=np.arange(length))

        eseq = Sequence("00000", id='id5', description='dsc5',
                        quality=[0, 0, 0, 0, 0])
        self.assertEquals(seq[0, 0, 0, 0, 0], eseq)
        self.assertEquals(seq[0, 0:1, 0, -length::-1, 0, 1:0], eseq)
        self.assertEquals(seq[0:1, 0:1, 0:1, 0:1, 0:1], eseq)
        self.assertEquals(seq[0:1, 0, 0, 0, 0], eseq)

        eseq = Sequence("0123fed9", id='id5', description='dsc5',
                        quality=[0, 1, 2, 3, 15, 14, 13, 9])
        self.assertEquals(seq[0, 1, 2, 3, 15, 14, 13, 9], eseq)
        self.assertEquals(seq[0, 1, 2, 3, :-4:-1, 9], eseq)
        self.assertEquals(seq[0:4, :-4:-1, 9, 1:0], eseq)
        self.assertEquals(seq[0:4, :-4:-1, 9:10], eseq)

    def test_getitem_with_tuple_of_mixed_no_qual(self):
        seq = Sequence("0123456789abcdef", id='id6', description="no_qual6")
        eseq = Sequence("0123fed9", id='id6', description='no_qual6')
        self.assertEquals(seq[0, 1, 2, 3, 15, 14, 13, 9], eseq)
        self.assertEquals(seq[0, 1, 2, 3, :-4:-1, 9], eseq)
        self.assertEquals(seq[0:4, :-4:-1, 9], eseq)
        self.assertEquals(seq[0:4, :-4:-1, 9:10], eseq)

    def test_getitem_with_iterable_of_mixed_has_qual(self):
        s = "0123456789abcdef"
        length = len(s)
        seq = Sequence(s, id='id7', description="dsc7",
                       quality=np.arange(length))

        def generator():
            yield slice(0, 4)
            yield slice(200, 400)
            yield -1
            yield slice(-2, -4, -1)
            yield 9

        eseq = Sequence("0123fed9", id='id7', description='dsc7',
                        quality=[0, 1, 2, 3, 15, 14, 13, 9])
        self.assertEquals(seq[[0, 1, 2, 3, 15, 14, 13, 9]], eseq)
        self.assertEquals(seq[generator()], eseq)
        self.assertEquals(seq[[slice(0, 4), slice(None, -4, -1), 9]], eseq)
        self.assertEquals(seq[
            [slice(0, 4), slice(None, -4, -1), slice(9, 10)]], eseq)

    def test_getitem_with_iterable_of_mixed_no_qual(self):
        s = "0123456789abcdef"
        seq = Sequence(s, id='id7', description="dsc7")

        def generator():
            yield slice(0, 4)
            yield slice(200, 400)
            yield slice(None, -4, -1)
            yield 9

        eseq = Sequence("0123fed9", id='id7', description='dsc7')
        self.assertEquals(seq[[0, 1, 2, 3, 15, 14, 13, 9]], eseq)
        self.assertEquals(seq[generator()], eseq)
        self.assertEquals(seq[[slice(0, 4), slice(None, -4, -1), 9]], eseq)
        self.assertEquals(seq[
            [slice(0, 4), slice(None, -4, -1), slice(9, 10)]], eseq)

    def test_getitem_with_numpy_index_has_qual(self):
        s = "0123456789abcdef"
        length = len(s)
        seq = Sequence(s, id='id9', description="dsc9",
                       quality=np.arange(length))

        eseq = Sequence("0123fed9", id='id9', description='dsc9',
                        quality=[0, 1, 2, 3, 15, 14, 13, 9])
        self.assertEquals(seq[np.array([0, 1, 2, 3, 15, 14, 13, 9])], eseq)

    def test_getitem_with_numpy_index_no_qual(self):
        s = "0123456789abcdef"
        seq = Sequence(s, id='id10', description="dsc10")

        eseq = Sequence("0123fed9", id='id10', description='dsc10')
        self.assertEquals(seq[np.array([0, 1, 2, 3, 15, 14, 13, 9])], eseq)

    def test_getitem_with_boolean_vector_has_qual(self):
        s = "0123456789abcdef"
        length = len(s)
        seq = Sequence(s, id='id11', description="dsc11",
                       quality=np.arange(length))

        eseq = Sequence("13579bdf", id='id11', description="dsc11",
                        quality=[1, 3, 5, 7, 9, 11, 13, 15])

        self.assertEqual(seq[np.array([False, True] * 8)], eseq)
        self.assertEqual(seq[[False, True] * 8], eseq)

    def test_getitem_with_boolean_vector_no_qual(self):
        s = "0123456789abcdef"
        seq = Sequence(s, id='id11', description="dsc11")

        eseq = Sequence("13579bdf", id='id11', description="dsc11")

        self.assertEqual(seq[np.array([False, True] * 8)], eseq)

    def test_getitem_with_invalid(self):
        seq = Sequence("123456", id='idm', description='description',
                       quality=[1, 2, 3, 4, 5, 6])

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
            seq[99999999999999999]

        with self.assertRaises(IndexError):
            seq[0, 0, 99999999999999999]

        # numpy 1.8.1 and 1.9.2 raise different error types
        # (ValueError, IndexError).
        with self.assertRaises(Exception):
            seq[100 * [True, False, True]]

    def test_len(self):
        self.assertEqual(len(Sequence("")), 0)
        self.assertEqual(len(Sequence("a")), 1)
        self.assertEqual(len(Sequence("abcdef")), 6)

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

        self.assertTrue(SequenceSubclass("A").sequence in Sequence("AAA"))

    def test_hash(self):
        with self.assertRaises(TypeError):
            hash(self.b1)
        self.assertNotIsInstance(self.b1, Hashable)

    def test_iter_has_quality(self):
        tested = False
        seq = Sequence("0123456789", id="a", description="b",
                       quality=np.arange(10))
        for i, s in enumerate(seq):
            tested = True
            self.assertEqual(s, Sequence(str(i), id='a', description='b',
                                         quality=[i]))
        self.assertTrue(tested)

    def test_iter_no_quality(self):
        tested = False
        seq = Sequence("0123456789", id="a", description="b")
        for i, s in enumerate(seq):
            tested = True
            self.assertEqual(s, Sequence(str(i), id='a', description='b'))
        self.assertTrue(tested)

    def test_reversed_has_quality(self):
        tested = False
        seq = Sequence("0123456789", id="a", description="b",
                       quality=np.arange(10))
        for i, s in enumerate(reversed(seq)):
            tested = True
            self.assertEqual(s, Sequence(str(9 - i), id='a', description='b',
                                         quality=[9 - i]))
        self.assertTrue(tested)

    def test_reversed_no_quality(self):
        tested = False
        seq = Sequence("0123456789", id="a", description="b")
        for i, s in enumerate(reversed(seq)):
            tested = True
            self.assertEqual(s, Sequence(str(9 - i), id='a', description='b'))
        self.assertTrue(tested)

    def test_repr(self):
        seq_simple = Sequence("ACGT")
        seq_med = Sequence("ACGT", id="id", description="desc",
                           quality=[1, 2, 3, 4])
        seq_complex = Sequence(("ASDKJHDJHFGUGF*&@KFHKHSDGKASDHGKDUYGKFHJ#&*YJ"
                                "FE&I@#JH@#ASJDHGF*&@#IG#*&IGUJKSADHAKSDJHI#*Y"
                                "LFUFLIU#RHL*Y#HHFLI#*FHL@#(*HJ"),
                               id="This is a long id", description="desc",
                               quality=([1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2] *
                                        10))
        self.assertEqual(repr(seq_simple), "Sequence('ACGT', length=4)")
        self.assertEqual(repr(seq_med),
                         ("Sequence('ACGT', length=4, id='id',"
                         " description='desc', quality=[1, 2, 3, 4])"))
        self.assertEqual(repr(seq_complex),
                         ("Sequence('ASDKJH ... @#(*HJ', length=120, id='This"
                          " is a long id', \n         description='desc', "
                          "quality=[1, 2, 3, 4, 5, 6, ..., 7, 8, 9, 0, 1, 2])")
                          )

    def test_str(self):
        self.assertEqual(str(Sequence("GATTACA")), "GATTACA")
        self.assertEqual(str(Sequence("ACCGGTACC")), "ACCGGTACC")
        self.assertEqual(str(Sequence("GREG")), "GREG")
        self.assertEqual(str(Sequence("ABC", quality=[1, 2, 3])), "ABC")
        self.assertIs(type(str(Sequence("A"))), str)

    def test_to_default_behavior(self):
        # minimal sequence, sequence with all optional attributes present, and
        # a subclass of Sequence
        for seq in self.b6, self.b8, RNA('ACGU', id='rna seq'):
            to = seq.to()
            self.assertTrue(seq.equals(to))
            self.assertFalse(seq is to)

    def test_to_update_single_attribute(self):
        to = self.b8.to(id='new id')
        self.assertFalse(self.b8 is to)

        # they don't compare equal when we compare all attributes...
        self.assertFalse(self.b8.equals(to))

        # ...but they *do* compare equal when we ignore id, as that was the
        # only attribute that changed
        self.assertTrue(self.b8.equals(to, ignore=['id']))

        # id should be what we specified in the to call...
        self.assertEqual(to.id, 'new id')

        # ..and shouldn't have changed on the original sequence
        self.assertEqual(self.b8.id, 'hello')

    def test_to_update_multiple_attributes(self):
        to = self.b8.to(id='new id', quality=range(20, 25),
                        sequence='ACGTA', description='new desc')
        self.assertFalse(self.b8 is to)
        self.assertFalse(self.b8.equals(to))

        # attributes should be what we specified in the to call...
        self.assertEqual(to.id, 'new id')
        npt.assert_array_equal(to.quality, np.array([20, 21, 22, 23, 24]))
        npt.assert_array_equal(to.sequence, np.array('ACGTA', dtype='c'))
        self.assertEqual(to.description, 'new desc')

        # ..and shouldn't have changed on the original sequence
        self.assertEqual(self.b8.id, 'hello')
        npt.assert_array_equal(self.b8.quality, range(11))
        npt.assert_array_equal(self.b8.sequence, np.array('HE..--..LLO',
                                                          dtype='c'))
        self.assertEqual(self.b8.description, 'gapped hello')

    def test_to_invalid_kwargs(self):
        with self.assertRaises(TypeError):
            self.b2.to(id='bar', unrecognized_kwarg='baz')

    def test_to_extra_non_attribute_kwargs(self):
        # test that we can pass through additional kwargs to the constructor
        # that aren't related to biological sequence attributes (i.e., they
        # aren't state that has to be copied)

        a = DNA('ACTG', description='foo')

        # should be able to `to` it b/c validate defaults to False
        b = a.to()
        self.assertTrue(a.equals(b))
        self.assertFalse(a is b)

        # specifying validate should raise an error when the copy is
        # instantiated
#        with self.assertRaises(SequenceError):
#            a.to(validate=True)

    def test_equals_sequences_without_metadata_compare_equal(self):
        self.assertTrue(Sequence('').equals(Sequence('')))
        self.assertTrue(Sequence('z').equals(Sequence('z')))
        self.assertTrue(
            Sequence('ACGT').equals(Sequence('ACGT')))

    def test_equals_sequences_with_metadata_compare_equal(self):
        seq1 = Sequence('ACGT', id='foo', description='abc',
                     quality=[1, 2, 3, 4])
        seq2 = Sequence('ACGT', id='foo', description='abc',
                     quality=[1, 2, 3, 4])
        self.assertTrue(seq1.equals(seq2))

        # order shouldn't matter
        self.assertTrue(seq2.equals(seq1))

    def test_equals_sequences_from_different_sources_compare_equal(self):
        # sequences that have the same data but are constructed from different
        # types of data should compare equal
        seq1 = Sequence('ACGT', id='foo', description='abc',
                        quality=(1, 2, 3, 4))
        seq2 = Sequence(np.array([65, 67, 71, 84], dtype=np.uint8),
                        id='foo', description='abc',
                        quality=np.array([1, 2, 3, 4]))
        self.assertTrue(seq1.equals(seq2))

    def test_equals_ignore_type(self):
        seq1 = Sequence('ACGT')
        seq2 = DNA('ACGT')
        self.assertTrue(seq1.equals(seq2, ignore=['type']))

    def test_equals_ignore_id(self):
        seq1 = Sequence('ACGT', id='foo')
        seq2 = Sequence('ACGT', id='bar')
        self.assertTrue(seq1.equals(seq2, ignore=['id']))

    def test_equals_ignore_description(self):
        seq1 = Sequence('ACGT', description='foo')
        seq2 = Sequence('ACGT', description='bar')
        self.assertTrue(seq1.equals(seq2, ignore=['description']))

    def test_equals_ignore_quality(self):
        seq1 = Sequence('ACGT', quality=[1, 2, 3, 4])
        seq2 = Sequence('ACGT', quality=[5, 6, 7, 8])
        self.assertTrue(seq1.equals(seq2, ignore=['quality']))

    def test_equals_ignore_sequence(self):
        seq1 = Sequence('ACGA')
        seq2 = Sequence('ACGT')
        self.assertTrue(seq1.equals(seq2, ignore=['sequence']))

    def test_equals_ignore_everything(self):
        seq1 = Sequence('ACGA', id='foo', description='abc',
                        quality=[1, 2, 3, 4])
        seq2 = DNA('ACGT', id='bar', description='def',
                   quality=[5, 6, 7, 8])
        self.assertTrue(seq1.equals(seq2,
                                    ignore=['quality', 'description', 'id',
                                            'sequence', 'type']))

    def test_equals_type_mismatch(self):
        seq1 = Sequence('ACGT', id='foo', description='abc',
                        quality=[1, 2, 3, 4])
        seq2 = DNA('ACGT', id='bar', description='def',
                   quality=[5, 6, 7, 8])
        self.assertFalse(seq1.equals(seq2,
                                     ignore=['quality', 'description', 'id']))

    def test_equals_id_mismatch(self):
        seq1 = Sequence('ACGT', id='foo')
        seq2 = Sequence('ACGT', id='bar')
        self.assertFalse(seq1.equals(seq2))

    def test_equals_description_mismatch(self):
        seq1 = Sequence('ACGT', description='foo')
        seq2 = Sequence('ACGT', description='bar')
        self.assertFalse(seq1.equals(seq2))

    def test_equals_quality_mismatch(self):
        # both provided
        seq1 = Sequence('ACGT', quality=[1, 2, 3, 4])
        seq2 = Sequence('ACGT', quality=[1, 2, 3, 5])
        self.assertFalse(seq1.equals(seq2))

        # one provided
        seq1 = Sequence('ACGT', quality=[1, 2, 3, 4])
        seq2 = Sequence('ACGT')
        self.assertFalse(seq1.equals(seq2))

    def test_equals_sequence_mismatch(self):
        seq1 = Sequence('ACGT')
        seq2 = Sequence('TGCA')
        self.assertFalse(seq1.equals(seq2))

    def test_count(self):
        def construct_char_array(s):
            return np.fromstring(s, dtype='|S1')

        def construct_uint8_array(s):
            return np.fromstring(s, dtype=np.uint8)

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

        self.assertEquals(tested, 4)

    def test_count_on_subclass(self):
        with self.assertRaises(TypeError) as cm:
            Sequence("abcd").count(SequenceSubclass("a"))

        self.assertIn("Sequence", str(cm.exception))
        self.assertIn("SequenceSubclass", str(cm.exception))

    def test_distance(self):
        tested = 0
        for constructor in self.sequence_kinds:
            tested += 1
            seq1 = Sequence("abcdef")
            seq2 = constructor("12bcef")

            self.assertIsInstance(seq1.distance(seq1), float)
            self.assertEqual(seq1.distance(seq2), 2.0/3.0)

        self.assertEqual(tested, 4)

    def test_distance_arbitrary_function(self):
        def metric(x, y):
            return len(x) ** 2 + len(y) ** 2

        seq1 = Sequence("12345678")
        seq2 = Sequence("1234")
        result = seq1.distance(seq2, metric=metric)
        self.assertIsInstance(result, float)
        self.assertEqual(result, 80.0)

    def test_distance_default_metric(self):
        seq1 = Sequence("abcdef")
        seq2 = Sequence("12bcef")
        seq_wrong = Sequence("abcdefghijklmnop")

        self.assertIsInstance(seq1.distance(seq1), float)
        self.assertEqual(seq1.distance(seq1), 0.0)
        self.assertEqual(seq1.distance(seq2), 2.0/3.0)

        with self.assertRaises(ValueError):
            seq1.distance(seq_wrong)

        with self.assertRaises(ValueError):
            seq_wrong.distance(seq1)

    def test_distance_on_subclass(self):
        seq1 = Sequence("abcdef")
        seq2 = SequenceSubclass("12bcef")

        with self.assertRaises(TypeError):
            seq1.distance(seq2)

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
                seq.index("A", begin=1, end=5)

        self.assertEqual(tested, 4)

    def test_index_on_subclass(self):
        with self.assertRaises(TypeError):
            Sequence("ABCDEFG").index(SequenceSubclass("A"))

        self.assertEqual(
            SequenceSubclass("ABCDEFG").index(SequenceSubclass("A")), 0)

    def _compare_kmers_results(self, observed, expected):
        for obs, exp in zip_longest(observed, expected, fillvalue=None):
            self.assertEqual(obs, exp)

    def test_kmers(self):
        seq = Sequence('GATTACA', quality=range(7))

        expected = [
            Sequence('G', quality=[0]),
            Sequence('A', quality=[1]),
            Sequence('T', quality=[2]),
            Sequence('T', quality=[3]),
            Sequence('A', quality=[4]),
            Sequence('C', quality=[5]),
            Sequence('A', quality=[6])
        ]
        self._compare_kmers_results(
            seq.kmers(1, overlap=False), expected)

        expected = [
            Sequence('GA', quality=[0, 1]),
            Sequence('TT', quality=[2, 3]),
            Sequence('AC', quality=[4, 5])
        ]
        self._compare_kmers_results(
            seq.kmers(2, overlap=False), expected)

        expected = [
            Sequence('GAT', quality=[0, 1, 2]),
            Sequence('TAC', quality=[3, 4, 5])
        ]
        self._compare_kmers_results(
            seq.kmers(3, overlap=False), expected)

        expected = [
            Sequence('GATTACA', quality=[0, 1, 2, 3, 4, 5, 6])
        ]
        self._compare_kmers_results(
            seq.kmers(7, overlap=False), expected)


        self.assertIs(type(seq.kmers(1)), GeneratorType)

    def test_kmers_with_overlap(self):
        seq = Sequence('GATTACA', quality=range(7))
        expected = [
            Sequence('G', quality=[0]),
            Sequence('A', quality=[1]),
            Sequence('T', quality=[2]),
            Sequence('T', quality=[3]),
            Sequence('A', quality=[4]),
            Sequence('C', quality=[5]),
            Sequence('A', quality=[6])
        ]
        self._compare_kmers_results(
            seq.kmers(1, overlap=True), expected)

        expected = [
            Sequence('GA', quality=[0, 1]),
            Sequence('AT', quality=[1, 2]),
            Sequence('TT', quality=[2, 3]),
            Sequence('TA', quality=[3, 4]),
            Sequence('AC', quality=[4, 5]),
            Sequence('CA', quality=[5, 6])
        ]
        self._compare_kmers_results(
            seq.kmers(2, overlap=True), expected)

        expected = [
            Sequence('GAT', quality=[0, 1, 2]),
            Sequence('ATT', quality=[1, 2, 3]),
            Sequence('TTA', quality=[2, 3, 4]),
            Sequence('TAC', quality=[3, 4, 5]),
            Sequence('ACA', quality=[4, 5, 6])
        ]
        self._compare_kmers_results(
            seq.kmers(3, overlap=True), expected)

        expected = [
            Sequence('GATTACA', quality=[0, 1, 2, 3, 4, 5, 6])
        ]
        self._compare_kmers_results(
            seq.kmers(7, overlap=True), expected)

    def test_kmers_invalid_k(self):
        seq = Sequence('GATTACA', quality=range(7))

        with self.assertRaises(ValueError):
            list(self.b1.kmers(0))

        with self.assertRaises(ValueError):
            list(self.b1.kmers(-42))

        with self.assertRaises(ValueError):
            list(self.b1.kmers(8))

    def test_kmers_different_sequences(self):
        seq = Sequence('HE..--..LLO', id='hello', description='gapped hello',
                       quality=range(11))
        expected = [
            Sequence('HE.', quality=[0, 1, 2], id='hello',
                     description='gapped hello'),
            Sequence('.--', quality=[3, 4, 5], id='hello',
                     description='gapped hello'),
            Sequence('..L', quality=[6, 7, 8], id='hello',
                     description='gapped hello')
        ]
        self._compare_kmers_results(seq.kmers(3, overlap=False), expected)

    def test_kmer_frequencies(self):
        seq = Sequence('GATTACA', quality=range(7))
        # overlap = True
        expected = Counter('GATTACA')
        self.assertEqual(seq.kmer_frequencies(1, overlap=True), expected)
        expected = Counter(['GAT', 'ATT', 'TTA', 'TAC', 'ACA'])
        self.assertEqual(seq.kmer_frequencies(3, overlap=True), expected)

        # overlap = False
        expected = Counter(['GAT', 'TAC'])
        self.assertEqual(seq.kmer_frequencies(3, overlap=False), expected)
        expected = Counter(['GATTACA'])
        self.assertEqual(seq.kmer_frequencies(7, overlap=False), expected)

    def test_kmer_frequencies_relative(self):
        seq = Sequence('GATTACA', quality=range(7))
        # overlap = True
        expected = defaultdict(float)
        expected['A'] = 3/7.
        expected['C'] = 1/7.
        expected['G'] = 1/7.
        expected['T'] = 2/7.
        self.assertEqual(seq.kmer_frequencies(1, overlap=True, relative=True),
                         expected)
        expected = defaultdict(float)
        expected['GAT'] = 1/5.
        expected['ATT'] = 1/5.
        expected['TTA'] = 1/5.
        expected['TAC'] = 1/5.
        expected['ACA'] = 1/5.
        self.assertEqual(seq.kmer_frequencies(3, overlap=True, relative=True),
                         expected)

        # overlap = False
        expected = defaultdict(float)
        expected['GAT'] = 1/2.
        expected['TAC'] = 1/2.
        self.assertEqual(seq.kmer_frequencies(3, overlap=False, relative=True),
                         expected)
        expected = defaultdict(float)
        expected['GATTACA'] = 1.0
        self.assertEqual(seq.kmer_frequencies(7, overlap=False, relative=True),
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
        self.assertEqual(seq.kmer_frequencies(1, relative=True),
                         defaultdict(float, {'A': 1.0}))

    def test_slices_from_regex(self):
        seq = Sequence('GATTACA', quality=range(7))
        pat = re.compile('(T+A)(CA)')

        obs = list(seq.slices_from_regex(pat))
        exp = [slice(2, 5), slice(5, 7)]
        self.assertEqual(obs, exp)

        self.assertIs(type(seq.slices_from_regex(pat)), GeneratorType)

    def test_slices_from_regex_no_groups(self):
        seq = Sequence('GATTACA', quality=range(7))
        pat = re.compile('(FOO)')
        self.assertEqual(list(seq.slices_from_regex(pat)), [])


if __name__ == "__main__":
    main()
