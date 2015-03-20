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

import numpy as np
import numpy.testing as npt
from scipy.spatial.distance import euclidean

from skbio import (
    Sequence, DNA, RNA,
    Protein)
from skbio.sequence import SequenceError

with hooks():
    from itertools import zip_longest


# Interface Tests
class SequenceInterfaceTest(object):
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

    def test___getitem__(self):
        pass

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

    def test_kmer(self):
        pass

    def test_kmer_counts(self):
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


class IUPACSequenceInterfaceTest(SequenceInterfaceTest):
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

class NucleotideInterfaceTest(IUPACSequenceInterfaceTest):
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
class TestSequence(SequenceInterfaceTest, TestCase):
    def setUp(self):
        pass

    def test___init__(self):
        pass


class TestProtein(IUPACSequenceInterfaceTest, TestCase):
    def setUp(self):
        pass

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


class TestDNA(NucleotideInterfaceTest, TestCase):
    def setUp(self):
        pass

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


class TestRNA(NucleotideInterfaceTest, TestCase):
    def setUp(self):
        pass

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

        # init as list
        b = Sequence(list('ACCGGXZY'))
        self.assertEqual(str(b), 'ACCGGXZY')
        self.assertEqual(b.id, "")
        self.assertEqual(b.description, "")

        # init as tuple
        b = Sequence(tuple('ACCGGXZY'))
        self.assertEqual(str(b), 'ACCGGXZY')
        self.assertEqual(b.id, "")
        self.assertEqual(b.description, "")

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
            Sequence('ACGT', id='a') ==
            Sequence('ACGT', id='b'))
        self.assertTrue(
            Sequence('ACGT', description='c') ==
            Sequence('ACGT', description='d'))
        self.assertTrue(
            Sequence('ACGT', id='a', description='c') ==
            Sequence('ACGT', id='b', description='d'))
        self.assertTrue(
            Sequence('ACGT', id='a', description='c',
                               quality=[1, 2, 3, 4]) ==
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
        with self.assertRaises(TypeError):
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
        for actual, expected in zip(b1_iter, "GATTACA"):
            self.assertEqual(actual, expected)

        self.assertRaises(StopIteration, lambda: next(b1_iter))

    def _compare_k_words_results(self, observed, expected):
        for obs, exp in zip_longest(observed, expected, fillvalue=None):
            # use equals to compare quality, id, description, sequence, and
            # type
            self.assertTrue(obs.equals(exp))

    def test_k_words_overlapping_true(self):
        expected = [
            Sequence('G', quality=[0]),
            Sequence('A', quality=[1]),
            Sequence('T', quality=[2]),
            Sequence('T', quality=[3]),
            Sequence('A', quality=[4]),
            Sequence('C', quality=[5]),
            Sequence('A', quality=[6])
        ]
        self._compare_k_words_results(
            self.b1.k_words(1, overlapping=True), expected)

        expected = [
            Sequence('GA', quality=[0, 1]),
            Sequence('AT', quality=[1, 2]),
            Sequence('TT', quality=[2, 3]),
            Sequence('TA', quality=[3, 4]),
            Sequence('AC', quality=[4, 5]),
            Sequence('CA', quality=[5, 6])
        ]
        self._compare_k_words_results(
            self.b1.k_words(2, overlapping=True), expected)

        expected = [
            Sequence('GAT', quality=[0, 1, 2]),
            Sequence('ATT', quality=[1, 2, 3]),
            Sequence('TTA', quality=[2, 3, 4]),
            Sequence('TAC', quality=[3, 4, 5]),
            Sequence('ACA', quality=[4, 5, 6])
        ]
        self._compare_k_words_results(
            self.b1.k_words(3, overlapping=True), expected)

        expected = [
            Sequence('GATTACA', quality=[0, 1, 2, 3, 4, 5, 6])
        ]
        self._compare_k_words_results(
            self.b1.k_words(7, overlapping=True), expected)

        self.assertEqual(list(self.b1.k_words(8, overlapping=True)), [])

    def test_k_words_overlapping_false(self):
        expected = [
            Sequence('G', quality=[0]),
            Sequence('A', quality=[1]),
            Sequence('T', quality=[2]),
            Sequence('T', quality=[3]),
            Sequence('A', quality=[4]),
            Sequence('C', quality=[5]),
            Sequence('A', quality=[6])
        ]
        self._compare_k_words_results(
            self.b1.k_words(1, overlapping=False), expected)

        expected = [
            Sequence('GA', quality=[0, 1]),
            Sequence('TT', quality=[2, 3]),
            Sequence('AC', quality=[4, 5])
        ]
        self._compare_k_words_results(
            self.b1.k_words(2, overlapping=False), expected)

        expected = [
            Sequence('GAT', quality=[0, 1, 2]),
            Sequence('TAC', quality=[3, 4, 5])
        ]
        self._compare_k_words_results(
            self.b1.k_words(3, overlapping=False), expected)

        expected = [
            Sequence('GATTACA', quality=[0, 1, 2, 3, 4, 5, 6])
        ]
        self._compare_k_words_results(
            self.b1.k_words(7, overlapping=False), expected)

        self.assertEqual(list(self.b1.k_words(8, overlapping=False)), [])

    def test_k_words_invalid_k(self):
        with self.assertRaises(ValueError):
            list(self.b1.k_words(0))

        with self.assertRaises(ValueError):
            list(self.b1.k_words(-42))

    def test_k_words_different_sequences(self):
        expected = [
            Sequence('HE.', quality=[0, 1, 2], id='hello',
                               description='gapped hello'),
            Sequence('.--', quality=[3, 4, 5], id='hello',
                               description='gapped hello'),
            Sequence('..L', quality=[6, 7, 8], id='hello',
                               description='gapped hello')
        ]
        self._compare_k_words_results(
            self.b8.k_words(3, overlapping=False), expected)


    def test_k_word_counts(self):
        # overlapping = True
        expected = Counter('GATTACA')
        self.assertEqual(self.b1.k_word_counts(1, overlapping=True),
                         expected)
        expected = Counter(['GAT', 'ATT', 'TTA', 'TAC', 'ACA'])
        self.assertEqual(self.b1.k_word_counts(3, overlapping=True),
                         expected)

        # overlapping = False
        expected = Counter(['GAT', 'TAC'])
        self.assertEqual(self.b1.k_word_counts(3, overlapping=False),
                         expected)
        expected = Counter(['GATTACA'])
        self.assertEqual(self.b1.k_word_counts(7, overlapping=False),
                         expected)

    def test_k_word_frequencies(self):
        # overlapping = True
        expected = defaultdict(float)
        expected['A'] = 3/7.
        expected['C'] = 1/7.
        expected['G'] = 1/7.
        expected['T'] = 2/7.
        self.assertEqual(self.b1.k_word_frequencies(1, overlapping=True),
                         expected)
        expected = defaultdict(float)
        expected['GAT'] = 1/5.
        expected['ATT'] = 1/5.
        expected['TTA'] = 1/5.
        expected['TAC'] = 1/5.
        expected['ACA'] = 1/5.
        self.assertEqual(self.b1.k_word_frequencies(3, overlapping=True),
                         expected)

        # overlapping = False
        expected = defaultdict(float)
        expected['GAT'] = 1/2.
        expected['TAC'] = 1/2.
        self.assertEqual(self.b1.k_word_frequencies(3, overlapping=False),
                         expected)
        expected = defaultdict(float)
        expected['GATTACA'] = 1.0
        self.assertEqual(self.b1.k_word_frequencies(7, overlapping=False),
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
        self.assertEqual(repr(self.b1),
                         "<Sequence: GATTACA (length: 7)>")
        self.assertEqual(repr(self.b6),
                         "<Sequence: ACGTACGTAC... (length: 12)>")

    def test_reversed(self):
        b1_reversed = reversed(self.b1)
        for actual, expected in zip(b1_reversed, "ACATTAG"):
            self.assertEqual(actual, expected)

        self.assertRaises(StopIteration, lambda: next(b1_reversed))

    def test_str(self):
        self.assertEqual(str(self.b1), "GATTACA")
        self.assertEqual(str(self.b2), "ACCGGTACC")
        self.assertEqual(str(self.b3), "GREG")

#    def test_gap_chars(self):
#        self.assertEqual(self.b1.gap_chars, set('-.'))

    def test_sequence(self):
        self.assertEqual(self.b1.sequence, "GATTACA")
        self.assertEqual(self.b2.sequence, "ACCGGTACC")
        self.assertEqual(self.b3.sequence, "GREG")

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

    def test_copy_default_behavior(self):
        # minimal sequence, sequence with all optional attributes present, and
        # a subclass of Sequence
        for seq in self.b6, self.b8, RNA('ACGU', id='rna seq'):
            copy = seq.copy()
            self.assertTrue(seq.equals(copy))
            self.assertFalse(seq is copy)

    def test_copy_update_single_attribute(self):
        copy = self.b8.copy(id='new id')
        self.assertFalse(self.b8 is copy)

        # they don't compare equal when we compare all attributes...
        self.assertFalse(self.b8.equals(copy))

        # ...but they *do* compare equal when we ignore id, as that was the
        # only attribute that changed
        self.assertTrue(self.b8.equals(copy, ignore=['id']))

        # id should be what we specified in the copy call...
        self.assertEqual(copy.id, 'new id')

        # ..and shouldn't have changed on the original sequence
        self.assertEqual(self.b8.id, 'hello')

    def test_copy_update_multiple_attributes(self):
        copy = self.b8.copy(id='new id', quality=range(20, 25),
                            sequence='ACGTA', description='new desc')
        self.assertFalse(self.b8 is copy)
        self.assertFalse(self.b8.equals(copy))

        # attributes should be what we specified in the copy call...
        self.assertEqual(copy.id, 'new id')
        npt.assert_equal(copy.quality, np.array([20, 21, 22, 23, 24]))
        self.assertEqual(copy.sequence, 'ACGTA')
        self.assertEqual(copy.description, 'new desc')

        # ..and shouldn't have changed on the original sequence
        self.assertEqual(self.b8.id, 'hello')
        npt.assert_equal(self.b8.quality, range(11))
        self.assertEqual(self.b8.sequence, 'HE..--..LLO')
        self.assertEqual(self.b8.description, 'gapped hello')

    def test_copy_invalid_kwargs(self):
        with self.assertRaises(TypeError):
            self.b2.copy(id='bar', unrecognized_kwarg='baz')

    def test_copy_extra_non_attribute_kwargs(self):
        # test that we can pass through additional kwargs to the constructor
        # that aren't related to biological sequence attributes (i.e., they
        # aren't state that has to be copied)

        a = DNA('ACTG', description='foo')

        # should be able to copy it b/c validate defaults to False
        b = a.copy()
        self.assertTrue(a.equals(b))
        self.assertFalse(a is b)

        # specifying validate should raise an error when the copy is
        # instantiated
#        with self.assertRaises(SequenceError):
#            a.copy(validate=True)

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
        with self.assertRaises(SequenceError):
            self.b1.distance(self.b2)

        # Alternate functions should also raise an error
        # Another distance function from scipy:
        with self.assertRaises(SequenceError):
            self.b1.distance(self.b2, distance_fn=euclidean)

        # Any other function should raise an error as well
        def dumb_distance(x, y):
            return 42

        with self.assertRaises(SequenceError):
            self.b1.distance(self.b2, distance_fn=dumb_distance)

    def test_fraction_diff(self):
        self.assertEqual(self.b1.fraction_diff(self.b1), 0., 5)
        self.assertEqual(
            self.b1.fraction_diff(Sequence('GATTACC')), 1. / 7., 5)

    def test_fraction_same(self):
        self.assertAlmostEqual(self.b1.fraction_same(self.b1), 1., 5)
        self.assertAlmostEqual(
            self.b1.fraction_same(Sequence('GATTACC')), 6. / 7., 5)

#    def test_gap_maps(self):
#        # in sequence with no gaps, the gap_maps are identical
#        self.assertEqual(self.b1.gap_maps(),
#                         ([0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 6]))
#        # in sequence with all gaps, the map of degapped to gapped is the empty
#        # list (bc its length is 0), and the map of gapped to degapped is all
#        # None
#        self.assertEqual(self.b7.gap_maps(),
#                         ([], [None, None, None, None, None, None]))
#
#        self.assertEqual(self.b8.gap_maps(),
#                         ([0, 1, 8, 9, 10],
#                          [0, 1, None, None, None, None, None, None, 2, 3, 4]))
#
#        # example from the gap_maps doc string
#        self.assertEqual(Sequence('-ACCGA-TA-').gap_maps(),
#                         ([1, 2, 3, 4, 5, 7, 8],
#                          [None, 0, 1, 2, 3, 4, None, 5, 6, None]))

#    def test_gap_vector(self):
#        self.assertEqual(self.b1.gap_vector(),
#                         [False] * len(self.b1))
#        self.assertEqual(self.b7.gap_vector(),
#                         [True] * len(self.b7))
#        self.assertEqual(self.b8.gap_vector(),
#                         [False, False, True, True, True, True,
#                          True, True, False, False, False])

    def test_index(self):
        self.assertEqual(self.b1.index('G'), 0)
        self.assertEqual(self.b1.index('A'), 1)
        self.assertEqual(self.b1.index('AC'), 4)
        self.assertRaises(ValueError, self.b1.index, 'x')

#    def test_is_gap(self):
#        self.assertTrue(self.b1.is_gap('.'))
#        self.assertTrue(self.b1.is_gap('-'))
#        self.assertFalse(self.b1.is_gap('A'))
#        self.assertFalse(self.b1.is_gap('x'))
#        self.assertFalse(self.b1.is_gap(' '))
#        self.assertFalse(self.b1.is_gap(''))

#    def test_is_gapped(self):
#        self.assertFalse(self.b1.is_gapped())
#        self.assertFalse(self.b2.is_gapped())
#        self.assertTrue(self.b7.is_gapped())
#        self.assertTrue(self.b8.is_gapped())

    def test_regex_iter(self):
        pat = re_compile('(T+A)(CA)')

        obs = list(self.b1.regex_iter(pat))
        exp = [(2, 5, 'TTA'), (5, 7, 'CA')]
        self.assertEqual(obs, exp)

        obs = list(self.b1.regex_iter(pat, retrieve_group_0=True))
        exp = [(2, 7, 'TTACA'), (2, 5, 'TTA'), (5, 7, 'CA')]
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
            'V', 'Y', 'a', 'c', 'b', 'd', 'g', 'h', 'k', 'm', 'n', 's', 'r',
            't', 'w', 'v', 'y', '-', '.'
        }

        self.assertEqual(self.b1.alphabet, exp)
        self.assertEqual(DNA.alphabet, exp)

    def test_gap_chars(self):
        self.assertEqual(self.b1.gap_chars, set('-.'))

    def test_complement_map(self):
        exp = {
            '-': '-', '.': '.', 'A': 'T', 'C': 'G', 'B': 'V', 'D': 'H',
            'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', 'S': 'S',
            'R': 'Y', 'T': 'A', 'W': 'W', 'V': 'B', 'Y': 'R', 'a': 't',
            'c': 'g', 'b': 'v', 'd': 'h', 'g': 'c', 'h': 'd', 'k': 'm',
            'm': 'k', 'n': 'n', 's': 's', 'r': 'y', 't': 'a', 'w': 'w',
            'v': 'b', 'y': 'r'
        }
        self.assertEqual(self.b1.complement_map, exp)
        self.assertEqual(DNA.complement_map, exp)

    def test_nondegenerate_chars(self):
        exp = set("ACGTacgt")
        self.assertEqual(self.b1.nondegenerate_chars, exp)
        self.assertEqual(DNA.nondegenerate_chars, exp)

    def test_degenerate_map(self):
        exp = {
            'B': set(['C', 'T', 'G']), 'D': set(['A', 'T', 'G']),
            'H': set(['A', 'C', 'T']), 'K': set(['T', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'T', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']), 'W': set(['A', 'T']),
            'V': set(['A', 'C', 'G']), 'Y': set(['C', 'T']),
            'b': set(['c', 't', 'g']), 'd': set(['a', 't', 'g']),
            'h': set(['a', 'c', 't']), 'k': set(['t', 'g']),
            'm': set(['a', 'c']), 'n': set(['a', 'c', 't', 'g']),
            's': set(['c', 'g']), 'r': set(['a', 'g']), 'w': set(['a', 't']),
            'v': set(['a', 'c', 'g']), 'y': set(['c', 't'])
        }
        self.assertEqual(self.b1.degenerate_map, exp)
        self.assertEqual(DNA.degenerate_map, exp)

    def test_degenerate_chars(self):
        exp = set(['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y',
                   'b', 'd', 'h', 'k', 'm', 'n', 's', 'r', 'w', 'v', 'y'])
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

    def test_nondegenerates_gap_mixed_case(self):
        exp = [DNA('-A.a'), DNA('-A.c'),
               DNA('-C.a'), DNA('-C.c')]
        obs = sorted(DNA('-M.m').nondegenerates(), key=str)
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
            'V', 'Y', 'a', 'c', 'b', 'd', 'g', 'h', 'k', 'm', 'n', 's', 'r',
            'u', 'w', 'v', 'y', '-', '.'
        }

        self.assertEqual(self.b1.alphabet, exp)
        self.assertEqual(RNA.alphabet, exp)

    def test_gap_chars(self):
        self.assertEqual(self.b1.gap_chars, set('-.'))

    def test_complement_map(self):
        exp = {
            '-': '-', '.': '.', 'A': 'U', 'C': 'G', 'B': 'V', 'D': 'H',
            'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', 'S': 'S',
            'R': 'Y', 'U': 'A', 'W': 'W', 'V': 'B', 'Y': 'R', 'a': 'u',
            'c': 'g', 'b': 'v', 'd': 'h', 'g': 'c', 'h': 'd', 'k': 'm',
            'm': 'k', 'n': 'n', 's': 's', 'r': 'y', 'u': 'a', 'w': 'w',
            'v': 'b', 'y': 'r'
        }
        self.assertEqual(self.b1.complement_map, exp)
        self.assertEqual(RNA.complement_map, exp)

    def test_nondegenerate_chars(self):
        exp = set("ACGUacgu")
        self.assertEqual(self.b1.nondegenerate_chars, exp)
        self.assertEqual(RNA.nondegenerate_chars, exp)

    def test_degenerate_map(self):
        exp = {
            'B': set(['C', 'U', 'G']), 'D': set(['A', 'U', 'G']),
            'H': set(['A', 'C', 'U']), 'K': set(['U', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'U', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']), 'W': set(['A', 'U']),
            'V': set(['A', 'C', 'G']), 'Y': set(['C', 'U']),
            'b': set(['c', 'u', 'g']), 'd': set(['a', 'u', 'g']),
            'h': set(['a', 'c', 'u']), 'k': set(['u', 'g']),
            'm': set(['a', 'c']), 'n': set(['a', 'c', 'u', 'g']),
            's': set(['c', 'g']), 'r': set(['a', 'g']), 'w': set(['a', 'u']),
            'v': set(['a', 'c', 'g']), 'y': set(['c', 'u'])
        }
        self.assertEqual(self.b1.degenerate_map, exp)
        self.assertEqual(RNA.degenerate_map, exp)

    def test_degenerate_chars(self):
        exp = set(['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y',
                   'b', 'd', 'h', 'k', 'm', 'n', 's', 'r', 'w', 'v', 'y'])
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

    def test_nondegenerates_gap_mixed_case(self):
        exp = [RNA('-A.a'), RNA('-A.c'),
               RNA('-C.a'), RNA('-C.c')]
        obs = sorted(RNA('-M.m').nondegenerates(), key=str)
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
            'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c',
            'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm', 'n', 'p', 'q', 'r',
            's', 't', 'v', 'w', 'x', 'y', 'z', '-', '.'
        }

        self.assertEqual(self.p1.alphabet, exp)
        self.assertEqual(Protein.alphabet, exp)

    def test_gap_chars(self):
        self.assertEqual(self.p1.gap_chars, set('-.'))

    def test_nondegenerate_chars(self):
        exp = set("ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy")
        self.assertEqual(self.p1.nondegenerate_chars, exp)
        self.assertEqual(Protein.nondegenerate_chars, exp)

    def test_degenerate_map(self):
        exp = {
            'B': set(['D', 'N']), 'Z': set(['E', 'Q']),
            'X': set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                      'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']),
            'b': set(['d', 'n']), 'z': set(['e', 'q']),
            'x': set(['a', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'k', 'l', 'm',
                      'n', 'p', 'q', 'r', 's', 't', 'v', 'w', 'y']),
        }
        self.assertEqual(self.p1.degenerate_map, exp)
        self.assertEqual(Protein.degenerate_map, exp)

    def test_degenerate_chars(self):
        exp = set(['B', 'X', 'Z', 'b', 'x', 'z'])
        self.assertEqual(self.p1.degenerate_chars, exp)
        self.assertEqual(Protein.degenerate_chars, exp)

    def test_nondegenerates(self):
        exp = [Protein('AD'), Protein('AN')]
        # Sort based on sequence string, as order is not guaranteed.
        obs = sorted(Protein('AB').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

if __name__ == "__main__":
    main()
