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
from collections import Counter, defaultdict
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


# Interface Tests
class SequenceInterfaceTests(object):
    def test_property_id(self):
        pass

    def test_property_sequence(self):
        pass

    def test_property_description(self):
        pass

    def test_property_quality(self):
        pass

    def test___contains__(self):
        pass

    def test___eq__(self):
        pass

    def test___ne__(self):
        pass

    def _generate_valid_slices(self, length):
        """Helper function used for getitem tests."""
        full_range = range(-length, length)
        # Cartesian product gives us all possible combinations of full_range.
        for start, stop, step in product(full_range, repeat=3):
            # Slice syntax is hard, so take the descriptivist approach.
            if step != 0 and ('a'*length)[start:stop:step] != '':
                yield slice(start, stop, step)
            else:
                # Would have yielded a slice that would be empty.
                # An empty sequence will raise an error from the
                # constructor.
                pass

    def _getitem_compenents(self):
        """Return the Class, a string, a quality array, and the length"""
        # Get our class and sample string
        Seq = self.cls
        string = self.sample_sequence
        # This will make an arbitrary numpy int array of correct length
        qual = np.fromstring(string, dtype=np.uint8)
        max_range = len(string)
        return Seq, string, qual, max_range


    def test___getitem___gives_new_sequence(self):
        Seq, string, _, _ = self._getitem_compenents()

        seq = Seq(string)
        self.assertFalse(seq is seq[:])

    def test___getitem___with_int_and_slice(self):
        Seq, string, qual, max_range = self._getitem_compenents()
        # These are the objects to test.
        seq = Seq(string, id='id', description='dsc', quality=qual)
        no_qual = Seq(string, id='idq', description='no_qual')

        # We should be able to index and slice into the entire range.
        for index in (list(range(-max_range, max_range)) +
                      list(self._generate_valid_slices(max_range))):
            # Expected string and quality data
            e_string = string[index]
            e_qual = qual[index]

            # Expected objects
            exp_seq = Seq(e_string, id='id', description='dsc', quality=e_qual)
            exp_no_qual = Seq(e_string, id='idq', description='no_qual')

            # Descriptive equality will raise an error if they are not equal.
            # Correctness of `equals(..., descriptive=True)` tested elsewhere.
            # self.assertTrue used as a fail-safe.
            self.assertTrue(seq[index].equals(exp_seq, descriptive=True))
            self.assertTrue(no_qual[index].equals(exp_no_qual,
                                                  descriptive=True))

    def test___getitem___with_tuple_and_list_of_int_and_slice(self):
        Seq, string, qual, max_range = self._getitem_compenents()
        # These are the objects to test.
        seq = Seq(string, id='id', description='dsc', quality=qual)
        no_qual = Seq(string, id='idq', description='no_qual')

        # Different cases to try
        slices = list(self._generate_valid_slices(max_range))
        indices = list(range(-max_range, max_range))
        mixed = (list(chain.from_iterable(zip(indices, slices))) +
                 list(chain.from_iterable(zip(indices, slices[::-1]))))

        for indexable in [slices, indices, mixed]:
            # Build the expected string and quality data
            e_string = ''
            e_qual = np.array([], dtype=np.int)
            for i in indexable:
                e_string += string[i]
                e_qual = np.append(e_qual, qual[i])

            # Expected objects
            exp_seq = Seq(e_string, id='id', description='dsc', quality=e_qual)
            exp_no_qual = Seq(e_string, id='idq', description='no_qual')

            # Descriptive equality will raise an error if they are not equal.
            # Correctness of `equals(..., descriptive=True)` tested elsewhere.
            # self.assertTrue used as a fail-safe.
            self.assertTrue(seq[tuple(indexable)].equals(exp_seq,
                                                         descriptive=True))
            self.assertTrue(seq[list(indexable)].equals(exp_seq,
                                                        descriptive=True))
            self.assertTrue(no_qual[tuple(indexable)].equals(exp_no_qual,
                                                             descriptive=True))
            self.assertTrue(no_qual[list(indexable)].equals(exp_no_qual,
                                                            descriptive=True))

    def test___getitem___with_mask(self):
        Seq, string, qual, max_range = self._getitem_compenents()

        # These are the objects to test.
        seq = Seq(string, id='id', description='dsc', quality=qual)
        no_qual = Seq(string, id='idq', description='no_qual')

        # Create an arbitrary mask, will be True, False, False repeated.
        np_mask = np.bincount(list(range(0, len(self.sample_sequence),
                                         3))).astype(np.bool)
        py_mask = list(np_mask)


        # Build the expected string and quality data
        e_string = ""
        e_qual = np.array([], dtype=np.int)
        for i in np.where(np_mask)[0]: # np.where is the inverse of bincount
            e_string += string[i]
            e_qual = np.append(e_qual, qual[i])

        # Expected objects
        exp_seq = Seq(e_string, id='id', description='dsc', quality=e_qual)
        exp_no_qual = Seq(e_string, id='idq', description='no_qual')

        # Descriptive equality will raise an error if they are not equal.
        # Correctness of `equals(..., descriptive=True)` tested elsewhere.
        # self.assertTrue used as a fail-safe.
        self.assertTrue(seq[np_mask].equals(exp_seq, descriptive=True))
        self.assertTrue(no_qual[np_mask].equals(exp_no_qual, descriptive=True))
        self.assertTrue(seq[py_mask].equals(exp_seq, descriptive=True))
        self.assertTrue(no_qual[py_mask].equals(exp_no_qual, descriptive=True))

    def test___getitem___with_short_mask(self):
        Seq, string, qual, max_range = self._getitem_compenents()

        # These are the objects to test.
        seq = Seq(string, id='id', description='dsc', quality=qual)
        no_qual = Seq(string, id='idq', description='no_qual')

        py_mask = [True, True, False, True]
        np_mask = np.asarray(py_mask)

        e_string = string[0] + string[1] + string[3]
        e_qual = np.array([], dtype=np.int)
        for a in [qual[0], qual[1], qual[3]]:
            e_qual = np.append(e_qual, a)

        # Expected objects
        exp_seq = Seq(e_string, id='id', description='dsc', quality=e_qual)
        exp_no_qual = Seq(e_string, id='idq', description='no_qual')

        # Descriptive equality will raise an error if they are not equal.
        # Correctness of `equals(..., descriptive=True)` tested elsewhere.
        # self.assertTrue used as a fail-safe.
        self.assertTrue(seq[np_mask].equals(exp_seq, descriptive=True))
        self.assertTrue(no_qual[np_mask].equals(exp_no_qual, descriptive=True))
        self.assertTrue(seq[py_mask].equals(exp_seq, descriptive=True))
        self.assertTrue(no_qual[py_mask].equals(exp_no_qual, descriptive=True))

    def test___getitem___with_invalid(self):
        Seq, string, qual, max_range = self._getitem_compenents()

        seq = Seq(string, id='idm', description='description', quality=qual)

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

        with self.assertRaises(ValueError):
            seq[99999999999999999:2]

        # numpy 1.8.1 and 1.9.2 raise different error types
        # (ValueError, IndexError).
        with self.assertRaises(Exception):
            seq[100 * [True, False, True]]

    def test___hash__(self):
        pass

    def test___iter__(self):
        pass

    def test___len__(self):
        pass

    def test___repr__(self):
        pass

    def test___reversed__(self):
        pass

    def test___str__(self):
        pass

    def test_to(self):
        pass

    def test_count(self):
        pass

    def test_distance(self):
        pass

    def test_equals(self):
        pass

    def test_hamming(self):
        pass

    def test_index(self):
        pass

    def test_rindex(self):
        pass

    def test_kmers(self):
        pass

    def test_kmers_counts(self):
        pass

    def test_iregex(self):
        pass

    def test_startswith(self):
        pass

    def test_endswith(self):
        pass

    def test_find(self):
        pass

    def test_rfind(self):
        pass

    def test_replace(self):
        pass

    def test_split(self):
        pass

    def test_rsplit(self):
        pass

    def test_reverse(self):
        pass


class IUPACSequenceInterfaceTests(SequenceInterfaceTests):
    def test___str__(self):
        seq = ''.join(self.cls.alphabet).lower() \
            + ''.join(self.cls.alphabet).upper()
        # IUPAC sequences are uppercase
        self.assertEqual(seq.upper(), str(self.cls(seq, case_insensitive=True)))

    def test_degap(self):
        pass

    def test_gaps(self):
        pass

    def test_has_gaps(self):
        pass

    def test_expand_degenerates(self):
        pass

    def test_split_gaps(self):
        pass

    def test_has_degenerates(self):
        pass

class NucleotideInterfaceTests(IUPACSequenceInterfaceTests):
    def test_property_rc(self):
        pass

    def test_complement(self):
        pass

    def test_find_features(self):
        pass

    def test_gc_content(self):
        pass

    def test_translate(self):
        pass


# Concrete Tests
class TestSequence(SequenceInterfaceTests, IDValidationTests, TestCase):
    def setUp(self):
        self.cls = Sequence
        self.sample_sequence = "This is a sample sequence."

        # These are defined for IDValidationTests
        self.id_cls = Sequence
        self.id_kwargs = {'sequence':'A'}


    def test___init__(self):
        pass


class TestProtein(IUPACSequenceInterfaceTests, IDValidationTests, TestCase):
    def setUp(self):
        self.cls = Protein
        self.sample_sequence = ''.join(Protein.alphabet)

        # These are defined for IDValidationTests
        self.id_cls = Protein
        self.id_kwargs = {'sequence':'A'}


    def test_property_alphabet(self):
        pass

    def test_property_gap_chars(self):
        pass

    def test_property_nondegenerate_chars(self):
        pass

    def test_property_degenerate_chars(self):
        pass

    def test_property_degenerate_map(self):
        pass

    def test___init__(self):
        pass


class TestDNA(NucleotideInterfaceTests, IDValidationTests, TestCase):
    def setUp(self):
        self.cls = DNA
        self.sample_sequence = ''.join(DNA.alphabet)

        # These are defined for IDValidationTests
        self.id_cls = DNA
        self.id_kwargs = {'sequence':'A'}


    def test_property_alphabet(self):
        pass

    def test_property_gap_chars(self):
        pass

    def test_property_nondegenerate_chars(self):
        pass

    def test_property_degenerate_chars(self):
        pass

    def test_property_degenerate_map(self):
        pass

    def test_property_complement_map(self):
        pass

    def test_transcribe(self):
        pass

    def test___init__(self):
        pass


class TestRNA(NucleotideInterfaceTests, IDValidationTests, TestCase):
    def setUp(self):
        self.cls = RNA
        self.sample_sequence = ''.join(RNA.alphabet)

        # These are defined for IDValidationTests
        self.id_cls = RNA
        self.id_kwargs = {'sequence':'A'}

    def test_property_alphabet(self):
        pass

    def test_property_gap_chars(self):
        pass

    def test_property_nondegenerate_chars(self):
        pass

    def test_property_degenerate_chars(self):
        pass

    def test_property_degenerate_map(self):
        pass

    def test_property_complement_map(self):
        pass

    def test_reverse_transcribe(self):
        pass

    def test___init__(self):
        pass

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

    def test_init_with_invalid_quality(self):
        # invalid dtype
        with self.assertRaises(TypeError):
            Sequence('ACGT', quality=[2, 3, 4.1, 5])

        # wrong number of dimensions (2-D)
        with self.assertRaisesRegexp(SequenceError, '1-D'):
            Sequence('ACGT', quality=[[2, 3], [4, 5]])

        # wrong number of elements
        with self.assertRaisesRegexp(SequenceError, '\(3\).*\(4\)'):
            Sequence('ACGT', quality=[2, 3, 4])

        # negatives
        with self.assertRaisesRegexp(SequenceError,
                                     'quality scores.*greater than.*zero'):
            Sequence('ACGT', quality=[2, 3, -1, 4])

    def test_contains(self):
        self.assertTrue('G' in self.b1)
        self.assertFalse('g' in self.b1)

    def test_eq_and_ne(self):
        self.assertTrue(self.b1 == self.b1)
        self.assertTrue(self.b2 == self.b2)
        self.assertTrue(self.b3 == self.b3)

        self.assertTrue(self.b1 != self.b3)
        self.assertTrue(self.b1 != self.b2)
        self.assertTrue(self.b2 != self.b3)

        # identicial sequences of the same type are equal, even if they have
        # different ids, descriptions, and/or quality
        self.assertTrue(
            Sequence('ACGT') == Sequence('ACGT'))
        self.assertTrue(
            Sequence('ACGT', id='a') !=
            Sequence('ACGT', id='b'))
        self.assertTrue(
            Sequence('ACGT', description='c') !=
            Sequence('ACGT', description='d'))
        self.assertTrue(
            Sequence('ACGT', id='a', description='c') !=
            Sequence('ACGT', id='b', description='d'))
        self.assertTrue(
            Sequence('ACGT', id='a', description='c',
                               quality=[1, 2, 3, 4]) !=
            Sequence('ACGT', id='b', description='d',
                               quality=[5, 6, 7, 8]))

        # different type causes sequences to not be equal
        self.assertFalse(
            Sequence('ACGT') == DNA('ACGT'))

    def test_getitem(self):
        # use equals method to ensure that id, description, and sliced
        # quality are correctly propagated to the resulting sequence
        self.assertTrue(self.b1[0].equals(
            Sequence('G', quality=(0,))))

        self.assertTrue(self.b1[:].equals(
            Sequence('GATTACA', quality=range(7))))

        self.assertTrue(self.b1[::-1].equals(
            Sequence('ACATTAG', quality=range(7)[::-1])))

        # test a sequence without quality scores
        b = Sequence('ACGT', id='foo', description='bar')
        self.assertTrue(b[2:].equals(
            Sequence('GT', id='foo', description='bar')))
        self.assertTrue(b[2].equals(
            Sequence('G', id='foo', description='bar')))

    def test_getitem_indices(self):
        # no ordering, repeated items
        self.assertTrue(self.b1[[3, 5, 4, 0, 5, 0]].equals(
            Sequence('TCAGCG', quality=(3, 5, 4, 0, 5, 0))))

        # single item
        self.assertTrue(
            self.b1[[2]].equals(Sequence('T', quality=(2,))))

        # negatives
        self.assertTrue(self.b1[[2, -2, 4]].equals(
            Sequence('TCA', quality=(2, 5, 4))))

        # tuple
        self.assertTrue(self.b1[1, 2, 3].equals(
            Sequence('ATT', quality=(1, 2, 3))))
        self.assertTrue(self.b1[(1, 2, 3)].equals(
            Sequence('ATT', quality=(1, 2, 3))))

        # test a sequence without quality scores
        self.assertTrue(self.b2[5, 4, 1].equals(
            Sequence('TGC', id='test-seq-2',
                               description='A test sequence')))

    def test_getitem_wrong_type(self):
        with self.assertRaises(IndexError):
            self.b1['1']

    def test_getitem_out_of_range(self):
        # seq with quality
        with self.assertRaises(IndexError):
            self.b1[42]
        with self.assertRaises(IndexError):
            self.b1[[1, 0, 23, 3]]

        # seq without quality
        with self.assertRaises(IndexError):
            self.b2[43]
        with self.assertRaises(IndexError):
            self.b2[[2, 3, 22, 1]]

    def test_hash(self):
        self.assertTrue(isinstance(hash(self.b1), int))

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


    def test_kmer_counts(self):
        # overlap = True
        expected = Counter('GATTACA')
        self.assertEqual(self.b1.kmer_counts(1, overlap=True),
                         expected)
        expected = Counter(['GAT', 'ATT', 'TTA', 'TAC', 'ACA'])
        self.assertEqual(self.b1.kmer_counts(3, overlap=True),
                         expected)

        # overlap = False
        expected = Counter(['GAT', 'TAC'])
        self.assertEqual(self.b1.kmer_counts(3, overlap=False),
                         expected)
        expected = Counter(['GATTACA'])
        self.assertEqual(self.b1.kmer_counts(7, overlap=False),
                         expected)

    def test_k_word_frequencies(self):
        # overlap = True
        expected = defaultdict(float)
        expected['A'] = 3/7.
        expected['C'] = 1/7.
        expected['G'] = 1/7.
        expected['T'] = 2/7.
        self.assertEqual(self.b1.k_word_frequencies(1, overlap=True),
                         expected)
        expected = defaultdict(float)
        expected['GAT'] = 1/5.
        expected['ATT'] = 1/5.
        expected['TTA'] = 1/5.
        expected['TAC'] = 1/5.
        expected['ACA'] = 1/5.
        self.assertEqual(self.b1.k_word_frequencies(3, overlap=True),
                         expected)

        # overlap = False
        expected = defaultdict(float)
        expected['GAT'] = 1/2.
        expected['TAC'] = 1/2.
        self.assertEqual(self.b1.k_word_frequencies(3, overlap=False),
                         expected)
        expected = defaultdict(float)
        expected['GATTACA'] = 1.0
        self.assertEqual(self.b1.k_word_frequencies(7, overlap=False),
                         expected)

    def test_k_word_frequencies_floating_point_precision(self):
        # Test that a sequence having no variation in k-words yields a
        # frequency of exactly 1.0. Note that it is important to use
        # self.assertEqual here instead of self.assertAlmostEqual because we
        # want to test for exactly 1.0. A previous implementation of
        # Sequence.k_word_frequencies added (1 / num_words) for each
        # occurrence of a k-word to compute the frequencies (see
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
        self.assertEqual(seq.k_word_frequencies(1),
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
        self.assertTrue(a.has_quality())

        b = Sequence('ACA')
        self.assertFalse(b.has_quality())

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


#class NucelotideSequenceTests(TestCase):
#
#    def setUp(self):
#        self.empty = NucleotideSequence('')
#        self.b1 = NucleotideSequence('GATTACA')
#        self.b2 = NucleotideSequence(
#            'ACCGGUACC', id="test-seq-2",
#            description="A test sequence")
#        self.b3 = NucleotideSequence('G-AT-TG.AT.T')
#
#    def test_alphabet(self):
#        exp = {
#            'A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'U', 'T',
#            'W', 'V', 'Y', 'a', 'c', 'b', 'd', 'g', 'h', 'k', 'm', 'n', 's',
#            'r', 'u', 't', 'w', 'v', 'y'
#        }
#
#        # Test calling from an instance and purely static context.
#        self.assertEqual(self.b1.alphabet, exp)
#        self.assertEqual(NucleotideSequence.alphabet, exp)
#
#    def test_gap_chars(self):
#        self.assertEqual(self.b1.gap_chars, set('-.'))
#
#    def test_complement_map(self):
#        exp = {}
#        self.assertEqual(self.b1.complement_map, exp)
#        self.assertEqual(NucleotideSequence.complement_map, exp)
#
#    def test_nondegenerate_chars(self):
#        exp = set("ACGTUacgtu")
#        self.assertEqual(self.b1.nondegenerate_chars, exp)
#        self.assertEqual(NucleotideSequence.nondegenerate_chars, exp)
#
#    def test_degenerate_map(self):
#        exp = {
#            # upper
#            'B': set(['C', 'U', 'T', 'G']), 'D': set(['A', 'U', 'T', 'G']),
#            'H': set(['A', 'C', 'U', 'T']), 'K': set(['U', 'T', 'G']),
#            'M': set(['A', 'C']), 'N': set(['A', 'C', 'U', 'T', 'G']),
#            'S': set(['C', 'G']), 'R': set(['A', 'G']),
#            'W': set(['A', 'U', 'T']), 'V': set(['A', 'C', 'G']),
#            'Y': set(['C', 'U', 'T']),
#            # lower
#            'b': set(['c', 'u', 't', 'g']), 'd': set(['a', 'u', 't', 'g']),
#            'h': set(['a', 'c', 'u', 't']), 'k': set(['u', 't', 'g']),
#            'm': set(['a', 'c']), 'n': set(['a', 'c', 'u', 't', 'g']),
#            's': set(['c', 'g']), 'r': set(['a', 'g']),
#            'w': set(['a', 'u', 't']), 'v': set(['a', 'c', 'g']),
#            'y': set(['c', 'u', 't'])
#        }
#        self.assertEqual(self.b1.degenerate_map, exp)
#        self.assertEqual(NucleotideSequence.degenerate_map, exp)
#
#        # Test that we can modify a copy of the mapping without altering the
#        # canonical representation.
#        degen = NucleotideSequence.degenerate_map
#        degen.update({'V': set("BRO"), 'Z': set("ZORRO")})
#        self.assertNotEqual(degen, exp)
#        self.assertEqual(NucleotideSequence.degenerate_map, exp)
#
#    def test_degenerate_chars(self):
#        exp = set(['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y',
#                   'b', 'd', 'h', 'k', 'm', 'n', 's', 'r', 'w', 'v', 'y'])
#        self.assertEqual(self.b1.degenerate_chars, exp)
#        self.assertEqual(NucleotideSequence.degenerate_chars, exp)
#
#    def test_complement(self):
#        self.assertRaises(SequenceError,
#                          self.b1.complement)
#
#    def test_reverse_complement(self):
#        self.assertRaises(SequenceError,
#                          self.b1.reverse_complement)
#
#    def test_is_reverse_complement(self):
#        self.assertRaises(SequenceError,
#                          self.b1.is_reverse_complement, self.b1)
#
#    def test_nondegenerates_invalid(self):
#        with self.assertRaises(SequenceError):
#            list(NucleotideSequence('AZA').nondegenerates())
#
#    def test_nondegenerates_empty(self):
#        self.assertEqual(list(self.empty.nondegenerates()), [self.empty])
#
#    def test_nondegenerates_no_degens(self):
#        self.assertEqual(list(self.b1.nondegenerates()), [self.b1])
#
#    def test_nondegenerates_all_degens(self):
#        # Same chars.
#        exp = [NucleotideSequence('CC'), NucleotideSequence('CG'),
#               NucleotideSequence('GC'), NucleotideSequence('GG')]
#        # Sort based on sequence string, as order is not guaranteed.
#        obs = sorted(NucleotideSequence('SS').nondegenerates(), key=str)
#        self.assertEqual(obs, exp)
#
#        # Different chars.
#        exp = [NucleotideSequence('AC'), NucleotideSequence('AG'),
#               NucleotideSequence('GC'), NucleotideSequence('GG')]
#        obs = sorted(NucleotideSequence('RS').nondegenerates(), key=str)
#        self.assertEqual(obs, exp)
#
#        # Odd number of chars.
#        obs = list(NucleotideSequence('NNN').nondegenerates())
#        self.assertEqual(len(obs), 5**3)
#
#    def test_nondegenerates_mixed_degens(self):
#        exp = [NucleotideSequence('AGC'), NucleotideSequence('AGT'),
#               NucleotideSequence('AGU'), NucleotideSequence('GGC'),
#               NucleotideSequence('GGT'), NucleotideSequence('GGU')]
#        obs = sorted(NucleotideSequence('RGY').nondegenerates(), key=str)
#        self.assertEqual(obs, exp)
#
#    def test_nondegenerates_gap_mixed_case(self):
#        exp = [NucleotideSequence('-A.a'), NucleotideSequence('-A.c'),
#               NucleotideSequence('-C.a'), NucleotideSequence('-C.c')]
#        obs = sorted(NucleotideSequence('-M.m').nondegenerates(), key=str)
#        self.assertEqual(obs, exp)
#
#    def test_find_features(self):
#        exp = [(0, 2, 'GA'), (4, 5, 'A'), (6, 7, 'A')]
#        obs = list(self.b1.find_features('purine_run'))
#        self.assertEqual(obs, exp)
#
#        exp = [(2, 4, 'TT'), (5, 6, 'C')]
#        obs = list(self.b1.find_features('pyrimidine_run'))
#        self.assertEqual(obs, exp)
#
#        exp = [(0, 1, 'A'), (3, 5, 'GG'), (6, 7, 'A')]
#        obs = list(self.b2.find_features('purine_run'))
#        self.assertEqual(obs, exp)
#
#        exp = [(1, 3, 'CC'), (5, 6, 'U'), (7, 9, 'CC')]
#        obs = list(self.b2.find_features('pyrimidine_run'))
#        self.assertEqual(obs, exp)
#
#    def test_find_features_min_length(self):
#        exp = [(0, 2, 'GA')]
#        obs = list(self.b1.find_features('purine_run', 2))
#        self.assertEqual(obs, exp)
#
#        exp = [(2, 4, 'TT')]
#        obs = list(self.b1.find_features('pyrimidine_run', 2))
#        self.assertEqual(obs, exp)
#
#        exp = [(3, 5, 'GG')]
#        obs = list(self.b2.find_features('purine_run', 2))
#        self.assertEqual(obs, exp)
#
#        exp = [(1, 3, 'CC'), (7, 9, 'CC')]
#        obs = list(self.b2.find_features('pyrimidine_run', 2))
#        self.assertEqual(obs, exp)
#
#    def test_find_features_no_feature_type(self):
#        with self.assertRaises(ValueError):
#            list(self.b1.find_features('nonexistent_feature_type'))
#
#    def test_find_features_allow_gaps(self):
#        exp = [(0, 3, 'G-A'), (6, 9, 'G.A')]
#        obs = list(self.b3.find_features('purine_run', 2, True))
#        self.assertEqual(obs, exp)
#
#        exp = [(3, 6, 'T-T'), (9, 12, 'T.T')]
#        obs = list(self.b3.find_features('pyrimidine_run', 2, True))
#        self.assertEqual(obs, exp)
#
#    def test_nondegenerates_propagate_optional_properties(self):
#        seq = NucleotideSequence('RS', id='foo', description='bar',
#                                 quality=[42, 999])
#
#        exp = [
#            NucleotideSequence('AC', id='foo', description='bar',
#                               quality=[42, 999]),
#            NucleotideSequence('AG', id='foo', description='bar',
#                               quality=[42, 999]),
#            NucleotideSequence('GC', id='foo', description='bar',
#                               quality=[42, 999]),
#            NucleotideSequence('GG', id='foo', description='bar',
#                               quality=[42, 999])
#        ]
#
#        obs = sorted(seq.nondegenerates(), key=str)
#
#        for o, e in zip(obs, exp):
#            # use equals method to ensure that id, description, and quality are
#            # correctly propagated to the resulting sequence
#            self.assertTrue(o.equals(e))


class DNATests(TestCase):

    def setUp(self):
        self.b1 = DNA('GATTACA')
        self.b2 = DNA('ACCGGTACC', id="test-seq-2",
                              description="A test sequence", quality=range(9))
        self.b4 = DNA(
            'MRWSYKVHDBN', id="degen",
            description="All of the degenerate bases")
        self.b5 = DNA('.G--ATTAC-A...')

    def test_alphabet(self):
        exp = {
            'A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'T', 'W',
            'V', 'Y', '-', '.'
        }

        self.assertEqual(self.b1.alphabet, exp)
        self.assertEqual(DNA.alphabet, exp)

    def test_gap_chars(self):
        self.assertEqual(self.b1.gap_chars, set('-.'))

    def test_complement_map(self):
        exp = {
            '-': '-', '.': '.', 'A': 'T', 'C': 'G', 'B': 'V', 'D': 'H',
            'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', 'S': 'S',
            'R': 'Y', 'T': 'A', 'W': 'W', 'V': 'B', 'Y': 'R'
        }
        self.assertEqual(self.b1.complement_map, exp)
        self.assertEqual(DNA.complement_map, exp)

    def test_nondegenerate_chars(self):
        exp = set("ACGT")
        self.assertEqual(self.b1.nondegenerate_chars, exp)
        self.assertEqual(DNA.nondegenerate_chars, exp)

    def test_degenerate_map(self):
        exp = {
            'B': set(['C', 'T', 'G']), 'D': set(['A', 'T', 'G']),
            'H': set(['A', 'C', 'T']), 'K': set(['T', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'T', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']), 'W': set(['A', 'T']),
            'V': set(['A', 'C', 'G']), 'Y': set(['C', 'T'])
        }
        self.assertEqual(self.b1.degenerate_map, exp)
        self.assertEqual(DNA.degenerate_map, exp)

    def test_degenerate_chars(self):
        exp = set(['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y'])
        self.assertEqual(self.b1.degenerate_chars, exp)
        self.assertEqual(DNA.degenerate_chars, exp)

    def test_complement(self):
        # use equals method to ensure that id, description, and quality are
        # correctly propagated to the resulting sequence
        self.assertTrue(self.b1.complement().equals(DNA("CTAATGT")))

        self.assertTrue(self.b2.complement().equals(
            DNA("TGGCCATGG", id="test-seq-2",
                        description="A test sequence", quality=range(9))))

        self.assertTrue(self.b4.complement().equals(
            DNA("KYWSRMBDHVN", id="degen",
                        description="All of the degenerate bases")))

        self.assertTrue(self.b5.complement().equals(
            DNA(".C--TAATG-T...")))

    def test_reverse_complement(self):
        # use equals method to ensure that id, description, and (reversed)
        # quality scores are correctly propagated to the resulting sequence
        self.assertTrue(self.b1.reverse_complement().equals(
            DNA("TGTAATC")))

        self.assertTrue(self.b2.reverse_complement().equals(
            DNA("GGTACCGGT", id="test-seq-2",
                        description="A test sequence",
                        quality=range(9)[::-1])))

        self.assertTrue(self.b4.reverse_complement().equals(
            DNA("NVHDBMRSWYK", id="degen",
                        description="All of the degenerate bases")))

    def test_is_reverse_complement(self):
        self.assertFalse(self.b1.is_reverse_complement(self.b1))

        # id, description, and quality scores should be ignored (only sequence
        # data and type should be compared)
        self.assertTrue(self.b1.is_reverse_complement(
            DNA('TGTAATC', quality=range(7))))

        self.assertTrue(
            self.b4.is_reverse_complement(DNA('NVHDBMRSWYK')))



    def test_nondegenerates_no_degens(self):
        self.assertEqual(list(self.b1.nondegenerates()), [self.b1])

    def test_nondegenerates_all_degens(self):
        # Same chars.
        exp = [DNA('CC'), DNA('CG'), DNA('GC'),
               DNA('GG')]
        # Sort based on sequence string, as order is not guaranteed.
        obs = sorted(DNA('SS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Different chars.
        exp = [DNA('AC'), DNA('AG'), DNA('GC'),
               DNA('GG')]
        obs = sorted(DNA('RS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Odd number of chars.
        obs = list(DNA('NNN').nondegenerates())
        self.assertEqual(len(obs), 4**3)

    def test_nondegenerates_mixed_degens(self):
        exp = [DNA('AGC'), DNA('AGT'), DNA('GGC'),
               DNA('GGT')]
        obs = sorted(DNA('RGY').nondegenerates(), key=str)
        self.assertEqual(obs, exp)


class RNATests(TestCase):

    def setUp(self):
        self.b1 = RNA('GAUUACA')
        self.b2 = RNA('ACCGGUACC', id="test-seq-2",
                              description="A test sequence", quality=range(9))
        self.b4 = RNA(
            'MRWSYKVHDBN', id="degen",
            description="All of the degenerate bases")
        self.b5 = RNA('.G--AUUAC-A...')

    def test_alphabet(self):
        exp = {
            'A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'U', 'W',
            'V', 'Y', '-', '.'
        }

        self.assertEqual(self.b1.alphabet, exp)
        self.assertEqual(RNA.alphabet, exp)

    def test_gap_chars(self):
        self.assertEqual(self.b1.gap_chars, set('-.'))

    def test_complement_map(self):
        exp = {
            '-': '-', '.': '.', 'A': 'U', 'C': 'G', 'B': 'V', 'D': 'H',
            'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', 'S': 'S',
            'R': 'Y', 'U': 'A', 'W': 'W', 'V': 'B', 'Y': 'R'
        }
        self.assertEqual(self.b1.complement_map, exp)
        self.assertEqual(RNA.complement_map, exp)

    def test_nondegenerate_chars(self):
        exp = set("ACGU")
        self.assertEqual(self.b1.nondegenerate_chars, exp)
        self.assertEqual(RNA.nondegenerate_chars, exp)

    def test_degenerate_map(self):
        exp = {
            'B': set(['C', 'U', 'G']), 'D': set(['A', 'U', 'G']),
            'H': set(['A', 'C', 'U']), 'K': set(['U', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'U', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']), 'W': set(['A', 'U']),
            'V': set(['A', 'C', 'G']), 'Y': set(['C', 'U'])
        }
        self.assertEqual(self.b1.degenerate_map, exp)
        self.assertEqual(RNA.degenerate_map, exp)

    def test_degenerate_chars(self):
        exp = set(['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y'])
        self.assertEqual(self.b1.degenerate_chars, exp)
        self.assertEqual(RNA.degenerate_chars, exp)

    def test_complement(self):
        # use equals method to ensure that id, description, and quality are
        # correctly propagated to the resulting sequence
        self.assertTrue(self.b1.complement().equals(RNA("CUAAUGU")))

        self.assertTrue(self.b2.complement().equals(
            RNA("UGGCCAUGG", id="test-seq-2",
                        description="A test sequence", quality=range(9))))


        self.assertTrue(self.b4.complement().equals(
            RNA("KYWSRMBDHVN", id="degen",
                        description="All of the degenerate bases")))

        self.assertTrue(self.b5.complement().equals(
            RNA(".C--UAAUG-U...")))

    def test_reverse_complement(self):
        # use equals method to ensure that id, description, and (reversed)
        # quality scores are correctly propagated to the resulting sequence
        self.assertTrue(self.b1.reverse_complement().equals(
            RNA("UGUAAUC")))

        self.assertTrue(self.b2.reverse_complement().equals(
            RNA("GGUACCGGU", id="test-seq-2",
                        description="A test sequence",
                        quality=range(9)[::-1])))


        self.assertTrue(self.b4.reverse_complement().equals(
            RNA("NVHDBMRSWYK", id="degen",
                        description="All of the degenerate bases")))

    def test_is_reverse_complement(self):
        self.assertFalse(self.b1.is_reverse_complement(self.b1))

        # id, description, and quality scores should be ignored (only sequence
        # data and type should be compared)
        self.assertTrue(self.b1.is_reverse_complement(
            RNA('UGUAAUC', quality=range(7))))

        self.assertTrue(
            self.b4.is_reverse_complement(RNA('NVHDBMRSWYK')))

    def test_nondegenerates_no_degens(self):
        self.assertEqual(list(self.b1.nondegenerates()), [self.b1])

    def test_nondegenerates_all_degens(self):
        # Same chars.
        exp = [RNA('CC'), RNA('CG'), RNA('GC'),
               RNA('GG')]
        # Sort based on sequence string, as order is not guaranteed.
        obs = sorted(RNA('SS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Different chars.
        exp = [RNA('AC'), RNA('AG'), RNA('GC'),
               RNA('GG')]
        obs = sorted(RNA('RS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Odd number of chars.
        obs = list(RNA('NNN').nondegenerates())
        self.assertEqual(len(obs), 4**3)

    def test_nondegenerates_mixed_degens(self):
        exp = [RNA('AGC'), RNA('AGU'), RNA('GGC'),
               RNA('GGU')]
        obs = sorted(RNA('RGY').nondegenerates(), key=str)
        self.assertEqual(obs, exp)


class ProteinTests(TestCase):

    def setUp(self):
        self.p1 = Protein('GREG')
        self.p2 = Protein(
            'PRTEINSEQNCE', id="test-seq-2",
            description="A test sequence")

    def test_alphabet(self):
        exp = {
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
            'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', '-', '.'
        }

        self.assertEqual(self.p1.alphabet, exp)
        self.assertEqual(Protein.alphabet, exp)

    def test_gap_chars(self):
        self.assertEqual(self.p1.gap_chars, set('-.'))

    def test_nondegenerate_chars(self):
        exp = set("ACDEFGHIKLMNPQRSTVWY")
        self.assertEqual(self.p1.nondegenerate_chars, exp)
        self.assertEqual(Protein.nondegenerate_chars, exp)

    def test_degenerate_map(self):
        exp = {
            'B': set(['D', 'N']), 'Z': set(['E', 'Q']),
            'X': set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                      'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
        }
        self.assertEqual(self.p1.degenerate_map, exp)
        self.assertEqual(Protein.degenerate_map, exp)

    def test_degenerate_chars(self):
        exp = set(['B', 'X', 'Z'])
        self.assertEqual(self.p1.degenerate_chars, exp)
        self.assertEqual(Protein.degenerate_chars, exp)

    def test_nondegenerates(self):
        exp = [Protein('AD'), Protein('AN')]
        # Sort based on sequence string, as order is not guaranteed.
        obs = sorted(Protein('AB').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

if __name__ == "__main__":
    main()
