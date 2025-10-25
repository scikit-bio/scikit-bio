# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from itertools import product

import numpy as np
import numpy.testing as npt

from skbio.sequence import Sequence, GrammaredSequence, DNA, RNA, Protein
from skbio.util import classproperty
from skbio.util._decorator import overrides

from skbio.sequence.distance import (
    _metric_specs, hamming, p_dist, kmer_distance, jc69_correct
)


class TestMetricSpecs(TestCase):
    def test_is_metric(self):
        """Test if the wrapped function is marked as a metric."""

        def metric1(seq1, seq2):
            return 1

        self.assertFalse(hasattr(metric1, "_is_metric"))

        @_metric_specs()
        def metric2(seq1, seq2):
            return 1

        self.assertTrue(hasattr(metric2, "_is_metric"))
        self.assertTrue(metric2._is_metric)

    def test_instance(self):
        """Test if input sequences are of matching Sequence objects."""

        @_metric_specs()
        def metric1(seq1, seq2):
            return 1
        
        seq1, seq2 = DNA("ACGT"), DNA("AGTC")
        self.assertEqual(metric1(seq1, seq2), 1)

        seq1, seq2 = Sequence("Hello"), Sequence("There")
        self.assertEqual(metric1(seq1, seq2), 1)

        seq1, seq2 = "Hello", "There"
        with self.assertRaises(TypeError) as cm:
            metric1(seq1, seq2)
        msg = "Sequences must be skbio.sequence.Sequence instances, not 'str'."
        self.assertEqual(str(cm.exception), msg)

        seq1, seq2 = DNA("ACGT"), Protein("MKVS")
        with self.assertRaises(TypeError) as cm:
            metric1(seq1, seq2)
        msg = "Sequences must have matching type. 'DNA' does not match 'Protein'."
        self.assertEqual(str(cm.exception), msg)

    def test_seqtype(self):
        """Test if input sequences are of expected sequence types."""

        @_metric_specs(seqtype=Protein)
        def metric1(seq1, seq2):
            return 1

        self.assertIs(metric1._seqtype, Protein)
        self.assertEqual(metric1(Protein("MVR"), Protein("TPD")), 1)

        with self.assertRaises(TypeError) as cm:
            metric1(DNA("GGC"), DNA("CAT"))
        msg = "Sequences must be 'Protein' instances, not 'DNA'."
        self.assertEqual(str(cm.exception), msg)

        @_metric_specs(seqtype=(DNA, RNA))
        def metric2(seq1, seq2):
            return 1

        self.assertTupleEqual(metric2._seqtype, (DNA, RNA))
        self.assertEqual(metric2(RNA("AUCG"), RNA("UAAC")), 1)

        with self.assertRaises(TypeError) as cm:
            metric2(Protein("MVR"), Protein("TPD"))
        msg = "Sequences must be ('DNA', 'RNA') instances, not 'Protein'."
        self.assertEqual(str(cm.exception), msg)

        @_metric_specs(seqtype=GrammaredSequence)
        def metric3(seq1, seq2):
            return 1

        class CustomSequence(GrammaredSequence):
            @classproperty
            @overrides(GrammaredSequence)
            def gap_chars(cls):
                return set('^$')

            @classproperty
            @overrides(GrammaredSequence)
            def default_gap_char(cls):
                return '^'

            @classproperty
            @overrides(GrammaredSequence)
            def definite_chars(cls):
                return set('WXYZ')

            @classproperty
            @overrides(GrammaredSequence)
            def degenerate_map(cls):
                return {}

        self.assertEqual(metric3(DNA("GGC"), DNA("CAT")), 1)
        self.assertEqual(metric3(RNA("AUCG"), RNA("UAAC")), 1)
        self.assertEqual(metric3(Protein("MVR"), Protein("TPD")), 1)
        self.assertEqual(metric3(CustomSequence("XXY"), CustomSequence("WWZ")), 1)

        with self.assertRaises(TypeError) as cm:
            metric3(Sequence("hello"), Sequence("there"))
        msg = "Sequences must be 'GrammaredSequence' instances, not 'Sequence'."
        self.assertEqual(str(cm.exception), msg)

    def test_equal(self):
        """Test if input sequences are of equal length."""
        seq1, seq2 = DNA("ACGT"), DNA("GATGC")

        @_metric_specs()
        def metric1(seq1, seq2):
            return 1

        self.assertFalse(metric1._equal)
        self.assertEqual(metric1(seq1, seq2), 1)

        @_metric_specs(equal=True)
        def metric2(seq1, seq2):
            return 1

        self.assertTrue(metric2._equal)
        with self.assertRaises(ValueError) as cm:
            metric2(seq1, seq2)
        msg = ("'metric2' can only be calculated between equal-length sequences. "
               "4 != 5.")
        self.assertEqual(str(cm.exception), msg)

    def test_alphabet(self):
        """Filter input sequences by a given alphabet."""
        seq1, seq2 = Sequence("X.,123Y"), Sequence("?'ZmYY0")

        @_metric_specs()
        def metric1(seq1, seq2):
            return len(seq1) * len(seq2)

        self.assertIsNone(metric1._alphabet)
        self.assertEqual(metric1(seq1, seq2), 49)

        @_metric_specs(alphabet="XYZ")
        def metric2(seq1, seq2):
            return len(seq1) * len(seq2)

        self.assertEqual(metric2._alphabet, "XYZ")
        self.assertEqual(metric2(seq1, seq2), 6)

        # DNA sequences with gaps and degenerate characters ("N" and "R")
        seq1, seq2 = DNA("ACCR--GT"), DNA("A-CGCANT")

        @_metric_specs(seqtype=(DNA, RNA), alphabet="nongap")
        def metric3(seq1, seq2):
            return len(seq1) * len(seq2)

        self.assertEqual(metric3(seq1, seq2), 42)

        @_metric_specs(seqtype=(DNA, RNA), equal=True, alphabet="nongap")
        def metric4(seq1, seq2):
            return len(seq1) * len(seq2)

        self.assertEqual(metric4(seq1, seq2), 25)

        @_metric_specs(seqtype=(DNA, RNA), equal=True, alphabet="definite")
        def metric5(seq1, seq2):
            return len(seq1) * len(seq2)

        self.assertEqual(metric5(seq1, seq2), 9)

        # protein sequences with degenerate ("X") and non-canonical characters ("O")
        seq1, seq2 = Protein("MNXSQ"), Protein("MKPWO")

        @_metric_specs(seqtype=Protein, equal=True, alphabet="definite")
        def metric6(seq1, seq2):
            return len(seq1) * len(seq2)

        self.assertEqual(metric6(seq1, seq2), 16)

        @_metric_specs(seqtype=Protein, equal=True, alphabet="canonical")
        def metric7(seq1, seq2):
            return len(seq1) * len(seq2)

        self.assertEqual(metric7(seq1, seq2), 9)

        # alphabet only works with grammared sequence
        @_metric_specs(alphabet="nongap")
        def metric8(seq1, seq2):
            return len(seq1) * len(seq2)

        self.assertRaises(AttributeError, metric8, Sequence("X"), Sequence("Y"))


class TestHamming(TestCase):
    def test_non_sequence(self):
        seq1 = Sequence('abc')
        seq2 = 'abc'
        with self.assertRaisesRegex(TypeError, "Sequence instances, not 'str'."):
            hamming(seq1, seq2)
        with self.assertRaisesRegex(TypeError, "Sequence instances, not 'str'."):
            hamming(seq2, seq1)

    def test_type_mismatch(self):
        seq1 = Sequence('ABC')
        seq2 = DNA('ACG')
        with self.assertRaisesRegex(TypeError, "'Sequence' does not match 'DNA'."):
            hamming(seq1, seq2)

    def test_length_mismatch(self):
        seq1 = Sequence('ABC')
        seq2 = Sequence('ABCD')
        with self.assertRaisesRegex(ValueError, "equal-length sequences. 3 != 4."):
            hamming(seq1, seq2)

    def test_return_type(self):
        seq1 = Sequence('ABC')
        seq2 = Sequence('ABC')
        obs = hamming(seq1, seq2)
        self.assertIsInstance(obs, float)
        self.assertEqual(obs, 0.0)

    def test_minimum_distance(self):
        seq1 = Sequence('ABC')
        seq2 = Sequence('ABC')
        obs = hamming(seq1, seq2)
        self.assertEqual(obs, 0.0)

    def test_mid_range_distance(self):
        seq1 = Sequence("abcdefgh")
        seq2 = Sequence("1b23ef45")
        obs = hamming(seq1, seq2)
        self.assertEqual(obs, 5.0/8.0)

    def test_maximum_distance(self):
        seq1 = Sequence('ABC')
        seq2 = Sequence('CAB')
        obs = hamming(seq1, seq2)
        self.assertEqual(obs, 1.0)

    def test_empty_sequences(self):
        seq1 = Sequence('')
        seq2 = Sequence('')
        obs = hamming(seq1, seq2)
        self.assertTrue(np.isnan(obs))

    def test_single_characters(self):
        seq1 = Sequence('a')
        seq2 = Sequence('b')
        self.assertEqual(hamming(seq1, seq1), 0.0)
        self.assertEqual(hamming(seq1, seq2), 1.0)

    def test_sequence_subclass(self):
        seq1 = DNA('ACG-T')
        seq2 = DNA('ACCTT')
        obs = hamming(seq1, seq2)
        self.assertEqual(obs, 2.0/5.0)

    def test_sequences_with_metadata(self):
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
            for seq1, seq2 in product(seqs, repeat=2):
                obs = hamming(seq1, seq2)
                self.assertEqual(obs, 0.0)

        for seq1, seq2 in product(seqs1, seqs2):
            obs = hamming(seq1, seq2)
            self.assertEqual(obs, 0.75)

    def test_raw_count(self):
        seq1 = DNA("AAGTC")
        seq2 = DNA("ACGAC")
        obs = hamming(seq1, seq2)
        self.assertEqual(obs, 0.4)
        obs = hamming(seq1, seq2, proportion=False)
        self.assertIsInstance(obs, float)
        self.assertEqual(obs, 2.0)

    def test_ignore_gap_degen(self):
        seq1 = DNA("ACGT")
        seq2 = DNA("TCGA")
        obs = hamming(seq1, seq2)
        self.assertEqual(obs, 0.5)
        obs = hamming(seq1, seq2, gaps=False)
        self.assertEqual(obs, 0.5)

        seq1 = DNA("AGCGT")
        seq2 = DNA("CG-AT")
        obs = hamming(seq1, seq2)
        self.assertEqual(obs, 0.6)
        obs = hamming(seq1, seq2, gaps=False)
        self.assertEqual(obs, 0.5)

        seq1 = DNA("AGCGT")
        seq2 = DNA("CGNAT")
        obs = hamming(seq1, seq2)
        self.assertEqual(obs, 0.6)
        obs = hamming(seq1, seq2, degenerates=False)
        self.assertEqual(obs, 0.5)

        seq1 = DNA("CARGT")
        seq2 = DNA("BAAGD")
        obs = hamming(seq1, seq2)
        self.assertEqual(obs, 0.6)
        obs = hamming(seq1, seq2, degenerates=False)
        self.assertEqual(obs, 0.0)

        seq1 = DNA("TARSTG-G")
        seq2 = DNA("C--ATNAG")
        obs = hamming(seq1, seq2)
        self.assertEqual(obs, 0.75)
        obs = hamming(seq1, seq2, gaps=False, degenerates=False)
        self.assertAlmostEqual(obs, 1 / 3)

    def test_non_left(self):
        seq1 = DNA("AAA---")
        seq2 = DNA("---TTT")
        obs = hamming(seq1, seq2)
        self.assertEqual(obs, 1.0)
        obs = hamming(seq1, seq2, gaps=False)
        self.assertTrue(np.isnan(obs))

    def test_gap_degen_undefined(self):
        seq1 = Sequence("AGCNT")
        seq2 = Sequence("CG-AT")
        with self.assertRaisesRegex(AttributeError, r"has no attribute 'gaps'"):
            hamming(seq1, seq2, gaps=False)
        with self.assertRaisesRegex(AttributeError, r"has no attribute 'degenerates'"):
            hamming(seq1, seq2, degenerates=False)


class TestPDist(TestCase):
    def test_p_dist(self):
        # sequences of canonical characters
        seq1 = DNA("AGATC")
        seq2 = DNA("AGATG")
        self.assertEqual(p_dist(seq1, seq2), 0.2)

        seq1 = RNA("AUCG")
        seq2 = RNA("UACG")
        self.assertEqual(p_dist(seq1, seq2), 0.5)

        seq1 = Protein("RCKMAF")
        seq2 = Protein("SCPTAA")
        self.assertAlmostEqual(p_dist(seq1, seq2), 2 / 3)

        # sequences with gaps
        seq1 = DNA("A--ACGG")
        seq2 = DNA("AGAAT-G")
        self.assertEqual(p_dist(seq1, seq2), 0.25)

        seq1 = Protein("-PYCRNG")
        seq2 = Protein("MPYAKC-")
        self.assertEqual(p_dist(seq1, seq2), 0.6)

        # sequences with degenerate characters
        seq1 = DNA("ANGCRT")
        seq2 = DNA("CCSMTT")
        self.assertEqual(p_dist(seq1, seq2), 0.5)

        seq1 = Protein("NBMKK")
        seq2 = Protein("HEMYX")
        self.assertAlmostEqual(p_dist(seq1, seq2), 2 / 3)

        # sequences with non-canonical characters
        seq1 = Protein("NKOC")
        seq2 = Protein("UKPA")
        self.assertAlmostEqual(p_dist(seq1, seq2), 0.5)

        # identical sequences
        seq1 = DNA("ACGT")
        seq2 = DNA("ACGT")
        self.assertEqual(p_dist(seq1, seq2), 0.0)

        # distinct sequences
        seq1 = DNA("ACGT")
        seq2 = DNA("TGCA")
        self.assertEqual(p_dist(seq1, seq2), 1.0)

        # single-character sequences
        seq1 = RNA("U")
        seq2 = RNA("G")
        self.assertEqual(p_dist(seq1, seq2), 1.0)

        # empty sequences
        seq1 = RNA("")
        seq2 = RNA("")
        self.assertTrue(np.isnan(p_dist(seq1, seq2)))

        # empty sequences after trimming
        seq1 = DNA("AAA---")
        seq2 = DNA("---TTT")
        self.assertTrue(np.isnan(p_dist(seq1, seq2)))

        seq1 = Protein("MGCPS")
        seq2 = Protein("XXXXX")
        self.assertTrue(np.isnan(p_dist(seq1, seq2)))

        # non-grammared sequences
        seq1 = Sequence("AGCNT")
        seq2 = Sequence("CG-AT")
        with self.assertRaises(TypeError):
            p_dist(seq1, seq2)


class TestKmerDistance(TestCase):
    def test_default_kwargs(self):
        seq1 = Sequence('AACCTAGCAATGGAT')
        seq2 = Sequence('CAGGCAGTTCTCACC')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 0.9130434782608695
        self.assertAlmostEqual(obs, exp)

    def test_nondefault_k(self):
        seq1 = Sequence('GCTTATGGAGAGAGA')
        seq2 = Sequence('CTCGAACTCCAGCCA')
        obs = kmer_distance(seq1, seq2, 2)
        exp = 0.7333333333333333
        self.assertAlmostEqual(obs, exp)
        seq1 = Sequence('EADDECAEECDEACD')
        seq2 = Sequence('DCBCBADADABCCDA')
        obs = kmer_distance(seq1, seq2, 1)
        exp = 0.4
        self.assertAlmostEqual(obs, exp)

    def test_overlap_false(self):
        seq1 = Sequence('CGTTATGTCTGTGAT')
        seq2 = Sequence('CTGAATCGGTAGTGT')
        obs = kmer_distance(seq1, seq2, 3, overlap=False)
        exp = 0.8888888888888888
        self.assertAlmostEqual(obs, exp)

    def test_entirely_different_sequences(self):
        seq1 = Sequence('CCGTGGTCGTATAAG')
        seq2 = Sequence('CGCCTTCCACATCAG')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 1.0
        self.assertEqual(obs, exp)

    def test_same_sequence(self):
        seq1 = Sequence('CTGCGACAGTTGGTA')
        seq2 = Sequence('CTGCGACAGTTGGTA')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 0.0
        self.assertEqual(obs, exp)

    def test_differing_length_seqs(self):
        seq1 = Sequence('AGAAATCTGAGCAAGGATCA')
        seq2 = Sequence('TTAGTGCGTAATCCG')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 0.9285714285714286
        self.assertAlmostEqual(obs, exp)

    def test_with_sequence_subclass(self):
        seq1 = DNA('GATGGTACTGTAGGT')
        seq2 = DNA('AGGGTGAAGGTATCA')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 0.8421052631578947
        self.assertAlmostEqual(obs, exp)

    def test_with_metadata_sanity(self):
        seq1 = Sequence('AACCTAGCAATGGAT',
                        metadata={'Name': 'Kestrel Gorlick'},
                        positional_metadata={'seq': list('ACTCAAGCTACGAAG')})
        seq2 = Sequence('CAGGCAGTTCTCACC')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 0.9130434782608695
        self.assertAlmostEqual(obs, exp)

    def test_return_type(self):
        seq1 = Sequence('ATCG')
        seq2 = Sequence('ATCG')
        obs = kmer_distance(seq1, seq2, 3)
        self.assertIsInstance(obs, float)
        self.assertEqual(obs, 0.0)

    def test_empty_sequences(self):
        seq1 = Sequence('')
        seq2 = Sequence('')
        obs = kmer_distance(seq1, seq2, 3)
        self.assertTrue(np.isnan(obs))

    def test_one_empty_sequence(self):
        seq1 = Sequence('')
        seq2 = Sequence('CGGGCAGCTCCTACCTGCTA')
        obs = kmer_distance(seq1, seq2, 3)
        exp = 1.0
        self.assertAlmostEqual(obs, exp)

    def test_no_kmers_found(self):
        seq1 = Sequence('ATCG')
        seq2 = Sequence('ACGT')
        obs = kmer_distance(seq1, seq2, 5)
        self.assertTrue(np.isnan(obs))

    def test_k_less_than_one_error(self):
        seq1 = Sequence('ATCG')
        seq2 = Sequence('ACTG')
        with self.assertRaisesRegex(ValueError, r'k must be greater than 0.'):
            kmer_distance(seq1, seq2, 0)

    def test_type_mismatch_error(self):
        seq1 = Sequence('ABC')
        seq2 = DNA('ATC')
        with self.assertRaisesRegex(TypeError, r"'Sequence' does not match 'DNA'"):
            kmer_distance(seq1, seq2, 3)

    def test_non_sequence_error(self):
        seq1 = Sequence('ATCG')
        seq2 = 'ATCG'
        with self.assertRaisesRegex(TypeError, r"not 'str'"):
            kmer_distance(seq1, seq2, 3)


class TestJC69(TestCase):
    def test_jc69_correct_scalar(self):
        obs = jc69_correct(0.1)
        self.assertIsInstance(obs, float)
        self.assertEqual(round(obs, 3), 0.107)
        self.assertEqual(round(jc69_correct(0.2), 3), 0.233)
        self.assertEqual(round(jc69_correct(0.5), 3), 0.824)
        self.assertEqual(round(jc69_correct(0.7), 3), 2.031)
        self.assertEqual(jc69_correct(0.0), 0.0)
        self.assertTrue(np.isnan(jc69_correct(0.8)))

    def test_jc69_correct_array(self):
        lst = [0.0, 0.1, 0.2, 0.5, 0.7, 1.0]
        obs = jc69_correct(lst)
        self.assertIsInstance(obs, np.ndarray)
        self.assertTupleEqual(obs.shape, (6,))
        exp = np.array([0.0, 0.107, 0.233, 0.824, 2.031, np.nan])
        npt.assert_array_equal(obs.round(3), exp)

        shape = (2, 3)
        arr = np.reshape(lst, shape)
        obs = jc69_correct(arr)
        self.assertTupleEqual(obs.shape, shape)
        exp = np.reshape(exp, shape)
        npt.assert_array_equal(obs.round(3), exp)

    def test_jc69_correct_alt_chars(self):
        self.assertEqual(round(jc69_correct(0.5, chars=5), 3), 0.785)
        self.assertEqual(round(jc69_correct(0.8, chars=9), 3), 2.047)

    def test_jc69_correct_error(self):
        with self.assertRaisesRegex(ValueError, r"`chars` must be at least 2."):
            jc69_correct(0.5, chars=1)


if __name__ == "__main__":
    main()
