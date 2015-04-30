# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.standard_library import hooks

from re import compile as re_compile
from collections import Counter, defaultdict, Hashable
from unittest import TestCase, main
from itertools import product, chain

import numpy as np
import numpy.testing as npt
from scipy.spatial.distance import euclidean

from skbio import (
    Sequence, DNA, RNA,
    Protein)
from skbio.sequence import SequenceError
from skbio.util._testing import IDValidationTests

with hooks():
    from itertools import zip_longest

class SequenceTests(TestCase):
    def setUp(self):
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

        # byte overflow
        #Sequence()

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

    def test_init_varied_input(self):
        # init as string
        b = Sequence('ACCGGXZY')
        self.assertEqual(str(b), 'ACCGGXZY')
        self.assertEqual(b.id, "")
        self.assertEqual(b.description, "")

        # init as string with optional values
        b = Sequence(
            'ACCGGXZY', 'test-seq-1', 'The first test sequence')
        self.assertEqual(str(b), 'ACCGGXZY')
        self.assertEqual(b.id, "test-seq-1")
        self.assertEqual(b.description, "The first test sequence")

        # test init as a different string
        b = Sequence('WRRTY')
        self.assertEqual(str(b), 'WRRTY')

    def test_sequence(self):
        npt.assert_array_equal(self.b1.sequence, np.array("GATTACA", dtype='c'))
        npt.assert_array_equal(self.b2.sequence,  np.array("ACCGGTACC", dtype='c'))
        npt.assert_array_equal(self.b3.sequence,  np.array("GREG", dtype='c'))

    def test_id(self):
        self.assertEqual(self.b1.id, "")
        self.assertEqual(self.b2.id, "test-seq-2")
        self.assertEqual(self.b3.id, "test-seq-3")

    def test_description(self):
        self.assertEqual(self.b1.description, "")
        self.assertEqual(self.b2.description, "A test sequence")
        self.assertEqual(self.b3.description, "A protein sequence")

    def test_quality(self):
        a = Sequence('ACA', quality=(22, 22, 1))

        # should get back a read-only numpy array of int dtype
        self.assertIsInstance(a.quality, np.ndarray)
        self.assertEqual(a.quality.dtype, np.int)
        npt.assert_equal(a.quality, np.array((22, 22, 1)))

        # test that we can't mutate the quality scores
        with self.assertRaises(ValueError):
            a.quality[1] = 42

        # test that we can't set the property
        with self.assertRaises(AttributeError):
            a.quality = (22, 22, 42)

    def test_quality_not_provided(self):
        b = Sequence('ACA')
        self.assertIs(b.quality, None)

    def test_quality_scalar(self):
        b = Sequence('G', quality=2)

        self.assertIsInstance(b.quality, np.ndarray)
        self.assertEqual(b.quality.dtype, np.int)
        self.assertEqual(b.quality.shape, (1,))
        npt.assert_equal(b.quality, np.array([2]))

    def test_quality_no_copy(self):
        qual = np.array([22, 22, 1])
        a = Sequence('ACA', quality=qual)
        self.assertIs(a.quality, qual)

        with self.assertRaises(ValueError):
            a.quality[1] = 42

        with self.assertRaises(ValueError):
            qual[1] = 42

    def test_has_quality(self):
        a = Sequence('ACA', quality=(5, 4, 67))
        self.assertTrue(a._has_quality())

        b = Sequence('ACA')
        self.assertFalse(b._has_quality())

    def test_eq_and_ne(self):
        seq_a = Sequence("A")
        seq_b = Sequence("B")

        class SequenceSubclass(Sequence):
            pass

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
            yield slice(None, -4, -1)
            yield 9

        eseq = Sequence("0123fed9", id='id7', description='dsc7',
                        quality=[0, 1, 2, 3, 15, 14, 13, 9])
        self.assertEquals(seq[[0, 1, 2, 3, 15, 14, 13, 9]], eseq)
        self.assertEquals(seq[generator()], eseq)
        self.assertEquals(seq[[slice(0, 4), slice(None, -4,-1), 9]], eseq)
        self.assertEquals(seq[[slice(0, 4), slice(None, -4,-1), slice(9, 10)]],
                          eseq)

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
        self.assertEquals(seq[[slice(0, 4), slice(None, -4,-1), 9]], eseq)
        self.assertEquals(seq[[slice(0, 4), slice(None, -4,-1), slice(9, 10)]],
                          eseq)

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
            seq[99999999999999999]

        with self.assertRaises(IndexError):
            seq[0, 0, 99999999999999999]

        # numpy 1.8.1 and 1.9.2 raise different error types
        # (ValueError, IndexError).
        with self.assertRaises(Exception):
            seq[100 * [True, False, True]]


    def test_contains(self):
        self.assertTrue('G' in self.b1)
        self.assertFalse('g' in self.b1)

    def test_hash(self):
        with self.assertRaises(TypeError):
            hash(self.b1)
        self.assertNotIsInstance(self.b1, Hashable)

    def test_iter(self):
        b1_iter = iter(self.b1)
        for actual, exp_c, exp_q in zip(b1_iter, "GATTACA", range(7)):
            expected = Sequence(exp_c, quality=exp_q)
            self.assertTrue(actual.equals(expected, descriptive=True))

        self.assertRaises(StopIteration, lambda: next(b1_iter))

    def _compare_kmers_results(self, observed, expected):
        for obs, exp in zip_longest(observed, expected, fillvalue=None):
            # use equals to compare quality, id, description, sequence, and
            # type
            self.assertTrue(obs.equals(exp))

    def test_kmers_overlap_true(self):
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
            self.b1.kmers(1, overlap=True), expected)

        expected = [
            Sequence('GA', quality=[0, 1]),
            Sequence('AT', quality=[1, 2]),
            Sequence('TT', quality=[2, 3]),
            Sequence('TA', quality=[3, 4]),
            Sequence('AC', quality=[4, 5]),
            Sequence('CA', quality=[5, 6])
        ]
        self._compare_kmers_results(
            self.b1.kmers(2, overlap=True), expected)

        expected = [
            Sequence('GAT', quality=[0, 1, 2]),
            Sequence('ATT', quality=[1, 2, 3]),
            Sequence('TTA', quality=[2, 3, 4]),
            Sequence('TAC', quality=[3, 4, 5]),
            Sequence('ACA', quality=[4, 5, 6])
        ]
        self._compare_kmers_results(
            self.b1.kmers(3, overlap=True), expected)

        expected = [
            Sequence('GATTACA', quality=[0, 1, 2, 3, 4, 5, 6])
        ]
        self._compare_kmers_results(
            self.b1.kmers(7, overlap=True), expected)

        self.assertEqual(list(self.b1.kmers(8, overlap=True)), [])

    def test_kmers_overlap_false(self):
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
            self.b1.kmers(1, overlap=False), expected)

        expected = [
            Sequence('GA', quality=[0, 1]),
            Sequence('TT', quality=[2, 3]),
            Sequence('AC', quality=[4, 5])
        ]
        self._compare_kmers_results(
            self.b1.kmers(2, overlap=False), expected)

        expected = [
            Sequence('GAT', quality=[0, 1, 2]),
            Sequence('TAC', quality=[3, 4, 5])
        ]
        self._compare_kmers_results(
            self.b1.kmers(3, overlap=False), expected)

        expected = [
            Sequence('GATTACA', quality=[0, 1, 2, 3, 4, 5, 6])
        ]
        self._compare_kmers_results(
            self.b1.kmers(7, overlap=False), expected)

        self.assertEqual(list(self.b1.kmers(8, overlap=False)), [])

    def test_kmers_invalid_k(self):
        with self.assertRaises(ValueError):
            list(self.b1.kmers(0))

        with self.assertRaises(ValueError):
            list(self.b1.kmers(-42))

    def test_kmers_different_sequences(self):
        expected = [
            Sequence('HE.', quality=[0, 1, 2], id='hello',
                               description='gapped hello'),
            Sequence('.--', quality=[3, 4, 5], id='hello',
                               description='gapped hello'),
            Sequence('..L', quality=[6, 7, 8], id='hello',
                               description='gapped hello')
        ]
        self._compare_kmers_results(
            self.b8.kmers(3, overlap=False), expected)


    def test_kmer_frequencies(self):
        # overlap = True
        expected = Counter('GATTACA')
        self.assertEqual(self.b1.kmer_frequencies(1, overlap=True),
                         expected)
        expected = Counter(['GAT', 'ATT', 'TTA', 'TAC', 'ACA'])
        self.assertEqual(self.b1.kmer_frequencies(3, overlap=True),
                         expected)

        # overlap = False
        expected = Counter(['GAT', 'TAC'])
        self.assertEqual(self.b1.kmer_frequencies(3, overlap=False),
                         expected)
        expected = Counter(['GATTACA'])
        self.assertEqual(self.b1.kmer_frequencies(7, overlap=False),
                         expected)

    def test_kmer_frequencies_relative(self):
        # overlap = True
        expected = defaultdict(float)
        expected['A'] = 3/7.
        expected['C'] = 1/7.
        expected['G'] = 1/7.
        expected['T'] = 2/7.
        self.assertEqual(self.b1.kmer_frequencies(1, overlap=True,
                                                  relative=True),
                         expected)
        expected = defaultdict(float)
        expected['GAT'] = 1/5.
        expected['ATT'] = 1/5.
        expected['TTA'] = 1/5.
        expected['TAC'] = 1/5.
        expected['ACA'] = 1/5.
        self.assertEqual(self.b1.kmer_frequencies(3, overlap=True,
                                                  relative=True),
                         expected)

        # overlap = False
        expected = defaultdict(float)
        expected['GAT'] = 1/2.
        expected['TAC'] = 1/2.
        self.assertEqual(self.b1.kmer_frequencies(3, overlap=False,
                                                  relative=True),
                         expected)
        expected = defaultdict(float)
        expected['GATTACA'] = 1.0
        self.assertEqual(self.b1.kmer_frequencies(7, overlap=False,
                                                  relative=True),
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

    def test_len(self):
        self.assertEqual(len(self.b1), 7)
        self.assertEqual(len(self.b2), 9)
        self.assertEqual(len(self.b3), 4)

    def test_repr(self):
        pass
        # self.assertEqual(repr(self.b1),
        #                  "<Sequence: GATTACA (length: 7)>")
        # self.assertEqual(repr(self.b6),
        #                  "<Sequence: ACGTACGTAC... (length: 12)>")

    def test_reversed(self):
        b1_reversed = reversed(self.b1)
        loop_ran = False
        for actual, exp_c, exp_q in zip(b1_reversed, "ACATTAG", reversed(range(7))):
            loop_ran = True
            expected = Sequence(exp_c, quality=exp_q)
            self.assertTrue(actual.equals(expected, descriptive=True))

        self.assertTrue(loop_ran)
        self.assertRaises(StopIteration, lambda: next(b1_reversed))

    def test_str(self):
        self.assertEqual(str(self.b1), "GATTACA")
        self.assertEqual(str(self.b2), "ACCGGTACC")
        self.assertEqual(str(self.b3), "GREG")

#    def test_gap_chars(self):
#        self.assertEqual(self.b1.gap_chars, set('-.'))

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
        npt.assert_array_equal(self.b8.sequence, np.array('HE..--..LLO', dtype='c'))
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

    def test_equals_true(self):
        # sequences match, all other attributes are not provided
        self.assertTrue(
            Sequence('ACGT').equals(Sequence('ACGT')))

        # all attributes are provided and match
        a = Sequence('ACGT', id='foo', description='abc',
                               quality=[1, 2, 3, 4])
        b = Sequence('ACGT', id='foo', description='abc',
                               quality=[1, 2, 3, 4])
        self.assertTrue(a.equals(b))

        # ignore type
        a = Sequence('ACGT')
        b = DNA('ACGT')
        self.assertTrue(a.equals(b, ignore=['type']))

        # ignore id
        a = Sequence('ACGT', id='foo')
        b = Sequence('ACGT', id='bar')
        self.assertTrue(a.equals(b, ignore=['id']))

        # ignore description
        a = Sequence('ACGT', description='foo')
        b = Sequence('ACGT', description='bar')
        self.assertTrue(a.equals(b, ignore=['description']))

        # ignore quality
        a = Sequence('ACGT', quality=[1, 2, 3, 4])
        b = Sequence('ACGT', quality=[5, 6, 7, 8])
        self.assertTrue(a.equals(b, ignore=['quality']))

        # ignore sequence
        a = Sequence('ACGA')
        b = Sequence('ACGT')
        self.assertTrue(a.equals(b, ignore=['sequence']))

        # ignore everything
        a = Sequence('ACGA', id='foo', description='abc',
                               quality=[1, 2, 3, 4])
        b = DNA('ACGT', id='bar', description='def',
                        quality=[5, 6, 7, 8])
        self.assertTrue(a.equals(b, ignore=['quality', 'description', 'id',
                                            'sequence', 'type']))

    def test_equals_false(self):
        # type mismatch
        a = Sequence('ACGT', id='foo', description='abc',
                               quality=[1, 2, 3, 4])
        b = DNA('ACGT', id='bar', description='def',
                        quality=[5, 6, 7, 8])
        self.assertFalse(a.equals(b, ignore=['quality', 'description', 'id']))

        # id mismatch
        a = Sequence('ACGT', id='foo')
        b = Sequence('ACGT', id='bar')
        self.assertFalse(a.equals(b))

        # description mismatch
        a = Sequence('ACGT', description='foo')
        b = Sequence('ACGT', description='bar')
        self.assertFalse(a.equals(b))

        # quality mismatch (both provided)
        a = Sequence('ACGT', quality=[1, 2, 3, 4])
        b = Sequence('ACGT', quality=[1, 2, 3, 5])
        self.assertFalse(a.equals(b))

        # quality mismatch (one provided)
        a = Sequence('ACGT', quality=[1, 2, 3, 4])
        b = Sequence('ACGT')
        self.assertFalse(a.equals(b))

        # sequence mismatch
        a = Sequence('ACGT')
        b = Sequence('TGCA')
        self.assertFalse(a.equals(b))

    def test_count(self):
        self.assertEqual(self.b1.count('A'), 3)
        self.assertEqual(self.b1.count('T'), 2)
        self.assertEqual(self.b1.count('TT'), 1)

#    def test_degap(self):
#        # use equals method to ensure that id, description, and filtered
#        # quality are correctly propagated to the resulting sequence
#
#        # no filtering, has quality
#        self.assertTrue(self.b1.degap().equals(self.b1))
#
#        # no filtering, doesn't have quality
#        self.assertTrue(self.b2.degap().equals(self.b2))
#
#        # everything is filtered, has quality
#        self.assertTrue(self.b7.degap().equals(
#            Sequence('', quality=[])))
#
#        # some filtering, has quality
#        self.assertTrue(self.b8.degap().equals(
#            Sequence('HELLO', id='hello', description='gapped hello',
#                               quality=[0, 1, 8, 9, 10])))

    def test_distance(self):
        # note that test_hamming_distance covers default behavior more
        # extensively
        self.assertEqual(self.b1.distance(self.b1), 0.0)
        self.assertEqual(self.b1.distance(Sequence('GATTACC')), 1./7)

        def dumb_distance(x, y):
            return 42

        self.assertEqual(
            self.b1.distance(self.b1, distance_fn=dumb_distance), 42)

    def test_distance_unequal_length(self):
        # distance requires sequences to be of equal length
        # While some functions passed to distance may throw an error not all
        # will. Therefore an error will be raised for sequences of unequal
        # length regardless of the function being passed.
        # With default hamming distance function
        with self.assertRaises(ValueError):
            self.b1.distance(self.b2)

        # Alternate functions should also raise an error
        # Another distance function from scipy:
        with self.assertRaises(ValueError):
            self.b1.distance(self.b2, distance_fn=euclidean)

        # Any other function should raise an error as well
        def dumb_distance(x, y):
            return 42

        with self.assertRaises(ValueError):
            self.b1.distance(self.b2, distance_fn=dumb_distance)

    def test_mismatch_frequency(self):
        # relative = False (default)
        self.assertEqual(self.b1.mismatch_frequency(self.b1), 0)
        self.assertEqual(self.b1.mismatch_frequency(Sequence('GATTACC')), 1)
        # relative = True
        self.assertEqual(self.b1.mismatch_frequency(self.b1, relative=True),
                         0., 5)
        self.assertEqual(self.b1.mismatch_frequency(Sequence('GATTACC'),
                                                    relative=True),
                         1. / 7., 5)

    def test_match_frequency(self):
        # relative = False (default)
        self.assertAlmostEqual(self.b1.match_frequency(self.b1), 7)
        self.assertAlmostEqual(
            self.b1.match_frequency(Sequence('GATTACC')), 6)
        # relative = True
        self.assertAlmostEqual(self.b1.match_frequency(self.b1, relative=True),
                               1., 5)
        self.assertAlmostEqual(self.b1.match_frequency(Sequence('GATTACC'),
                                                       relative=True),
                               6. / 7., 5)

    def test_index(self):
        self.assertEqual(self.b1.index('G'), 0)
        self.assertEqual(self.b1.index('A'), 1)
        self.assertEqual(self.b1.index('AC'), 4)
        self.assertRaises(ValueError, self.b1.index, 'x')

    # def test_has_gaps(self):
    #     self.assertFalse(self.b1.has_gaps())
    #     self.assertFalse(self.b2.has_gaps())
    #     self.assertTrue(self.b7.has_gaps())
    #     self.assertTrue(self.b8.has_gaps())

    def test_regex_iter(self):
        pat = re_compile('(T+A)(CA)')

        obs = list(self.b1.regex_iter(pat))
        exp = [slice(2, 5), slice(5, 7)]
        self.assertEqual(obs, exp)

        obs = list(self.b1.regex_iter(pat, retrieve_group_0=True))
        exp = [slice(2, 7), slice(2, 5), slice(5, 7)]
        self.assertEqual(obs, exp)


if __name__ == "__main__":
    main()
