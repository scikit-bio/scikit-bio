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
    _metric_specs, _char_hash, _char_freqs, hamming, pdist, kmer_distance, jc69,
    jc69_correct, f81, k2p, f84, tn93
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
            metric1("GGC", "CAT")
        msg = "Sequences must be skbio.sequence.Sequence instances, not 'str'."
        self.assertEqual(str(cm.exception), msg)

        with self.assertRaises(TypeError) as cm:
            metric1(Protein("MVR"), DNA("CAT"))
        msg = "Sequences must have matching type. 'Protein' does not match 'DNA'."
        self.assertEqual(str(cm.exception), msg)

        with self.assertRaises(TypeError) as cm:
            metric1(DNA("GGC"), DNA("CAT"))
        msg = "'metric1' is compatible with 'Protein' sequences, not 'DNA'."
        self.assertEqual(str(cm.exception), msg)

        @_metric_specs(seqtype=(DNA, RNA))
        def metric2(seq1, seq2):
            return 1

        self.assertTupleEqual(metric2._seqtype, (DNA, RNA))
        self.assertEqual(metric2(RNA("AUCG"), RNA("UAAC")), 1)

        with self.assertRaises(TypeError) as cm:
            metric2(Protein("MVR"), Protein("TPD"))
        msg = "'metric2' is compatible with ('DNA', 'RNA') sequences, not 'Protein'."
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
        msg = ("'metric3' is compatible with 'GrammaredSequence' sequences, not "
               "'Sequence'.")
        self.assertEqual(str(cm.exception), msg)

    def test_char_hash(self):
        """Hash table of valid characters."""
        # Note: A more intuitive test of this function is in skbio.alignment.distance.
        # DNA sequence
        obs = _char_hash(None, DNA)
        self.assertIsNone(obs)

        obs = _char_hash("ACGT", DNA)
        exp = np.zeros(128, dtype=bool)
        exp[[65, 67, 71, 84]] = True
        npt.assert_array_equal(obs, exp)

        obs = _char_hash("definite", DNA)
        npt.assert_array_equal(obs, exp)

        obs = _char_hash("canonical", DNA)
        npt.assert_array_equal(obs, exp)

        obs = _char_hash("ACGTN", DNA)
        exp[78] = True
        npt.assert_array_equal(obs, exp)

        obs = _char_hash("nongap", DNA)
        exp[sorted(map(ord, DNA.degenerate_chars))] = True
        npt.assert_array_equal(obs, exp)

        # RNA sequence
        obs = _char_hash("definite", RNA)
        exp = np.zeros(128, dtype=bool)
        exp[[65, 67, 71, 85]] = True
        npt.assert_array_equal(obs, exp)

        # protein sequence
        obs = _char_hash("definite", Protein)
        exp = np.zeros(128, dtype=bool)
        exp[sorted(map(ord, Protein.definite_chars))] = True
        npt.assert_array_equal(obs, exp)

        obs = _char_hash("canonical", Protein)
        exp[sorted(map(ord, Protein.noncanonical_chars))] = False
        npt.assert_array_equal(obs, exp)

        obs = _char_hash("ABCD", Protein)
        exp = np.zeros(128, dtype=bool)
        exp[65:69] = True
        npt.assert_array_equal(obs, exp)

    def test_char_freqs(self):
        # DNA characters
        valid = _char_hash("definite", DNA)
        
        # one 1D array
        seq1 = DNA("CGATCATCTA")
        obs = _char_freqs(seq1._bytes, valid)
        exp = np.array([.3, .3, .1, .3])
        npt.assert_array_equal(obs, exp)

        # two 1D arrays
        seq2 = DNA("CTGGCACCGA")
        obs = _char_freqs((seq1._bytes, seq2._bytes), valid)
        exp = np.array([.25, .35, .2, .2])
        npt.assert_array_equal(obs, exp)

        # 2D array
        seqs = np.vstack((seq1._bytes, seq2._bytes))
        obs = _char_freqs(seqs, valid)
        exp = np.array([.25, .35, .2, .2])
        npt.assert_array_equal(obs, exp)

        # all ASCII codes
        obs = _char_freqs(seqs)
        self.assertTupleEqual(obs.shape, (128,))
        self.assertEqual(obs[65], .25)
        self.assertEqual(obs[67], .35)
        self.assertEqual(obs[71], .2)
        self.assertEqual(obs[84], .2)

        # with gaps
        seq3 = DNA("CGATC---ATCTA")
        obs = _char_freqs(seq3._bytes, valid)
        exp = np.array([.3, .3, .1, .3])
        npt.assert_array_equal(obs, exp)

        # RNA characters
        valid = _char_hash("definite", RNA)
        seq4 = RNA("CGAUCAUCUA")
        obs = _char_freqs(seq4._bytes, valid)
        exp = np.array([.3, .3, .1, .3])
        npt.assert_array_equal(obs, exp)

        # empty sequence
        seq5 = RNA("")
        obs = _char_freqs(seq5._bytes, valid)
        exp = np.full(4, np.nan)
        npt.assert_array_equal(obs, exp)

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

    # def test_ignore_gap_degen(self):
    #     seq1 = DNA("ACGT")
    #     seq2 = DNA("TCGA")
    #     obs = hamming(seq1, seq2)
    #     self.assertEqual(obs, 0.5)
    #     obs = hamming(seq1, seq2, gaps=False)
    #     self.assertEqual(obs, 0.5)

    #     seq1 = DNA("AGCGT")
    #     seq2 = DNA("CG-AT")
    #     obs = hamming(seq1, seq2)
    #     self.assertEqual(obs, 0.6)
    #     obs = hamming(seq1, seq2, gaps=False)
    #     self.assertEqual(obs, 0.5)

    #     seq1 = DNA("AGCGT")
    #     seq2 = DNA("CGNAT")
    #     obs = hamming(seq1, seq2)
    #     self.assertEqual(obs, 0.6)
    #     obs = hamming(seq1, seq2, degenerates=False)
    #     self.assertEqual(obs, 0.5)

    #     seq1 = DNA("CARGT")
    #     seq2 = DNA("BAAGD")
    #     obs = hamming(seq1, seq2)
    #     self.assertEqual(obs, 0.6)
    #     obs = hamming(seq1, seq2, degenerates=False)
    #     self.assertEqual(obs, 0.0)

    #     seq1 = DNA("TARSTG-G")
    #     seq2 = DNA("C--ATNAG")
    #     obs = hamming(seq1, seq2)
    #     self.assertEqual(obs, 0.75)
    #     obs = hamming(seq1, seq2, gaps=False, degenerates=False)
    #     self.assertAlmostEqual(obs, 1 / 3)

    # def test_non_left(self):
    #     seq1 = DNA("AAA---")
    #     seq2 = DNA("---TTT")
    #     obs = hamming(seq1, seq2)
    #     self.assertEqual(obs, 1.0)
    #     obs = hamming(seq1, seq2, gaps=False)
    #     self.assertTrue(np.isnan(obs))

    # def test_gap_degen_undefined(self):
    #     seq1 = Sequence("AGCNT")
    #     seq2 = Sequence("CG-AT")
    #     with self.assertRaisesRegex(AttributeError, r"has no attribute 'gaps'"):
    #         hamming(seq1, seq2, gaps=False)
    #     with self.assertRaisesRegex(AttributeError, r"has no attribute 'degenerates'"):
    #         hamming(seq1, seq2, degenerates=False)


class TestPDist(TestCase):
    def test_pdist(self):
        # sequences of canonical characters
        seq1 = DNA("AGATC")
        seq2 = DNA("AGATG")
        self.assertEqual(pdist(seq1, seq2), 0.2)

        seq1 = RNA("AUCG")
        seq2 = RNA("UACG")
        self.assertEqual(pdist(seq1, seq2), 0.5)

        seq1 = Protein("RCKMAF")
        seq2 = Protein("SCPTAA")
        self.assertAlmostEqual(pdist(seq1, seq2), 2 / 3)

        # sequences with gaps
        seq1 = DNA("A--ACGG")
        seq2 = DNA("AGAAT-G")
        self.assertEqual(pdist(seq1, seq2), 0.25)

        seq1 = Protein("-PYCRNG")
        seq2 = Protein("MPYAKC-")
        self.assertEqual(pdist(seq1, seq2), 0.6)

        # sequences with degenerate characters
        seq1 = DNA("ANGCRT")
        seq2 = DNA("CCSMTT")
        self.assertEqual(pdist(seq1, seq2), 0.5)

        seq1 = Protein("NBMKK")
        seq2 = Protein("HEMYX")
        self.assertAlmostEqual(pdist(seq1, seq2), 2 / 3)

        # sequences with non-canonical characters
        # seq1 = Protein("NKOC")
        # seq2 = Protein("UKPA")
        # self.assertAlmostEqual(pdist(seq1, seq2), 0.5)

        # identical sequences
        seq1 = DNA("ACGT")
        seq2 = DNA("ACGT")
        self.assertEqual(pdist(seq1, seq2), 0.0)

        # distinct sequences
        seq1 = DNA("ACGT")
        seq2 = DNA("TGCA")
        self.assertEqual(pdist(seq1, seq2), 1.0)

        # single-character sequences
        seq1 = RNA("U")
        seq2 = RNA("G")
        self.assertEqual(pdist(seq1, seq2), 1.0)

        # empty sequences
        seq1 = RNA("")
        seq2 = RNA("")
        self.assertTrue(np.isnan(pdist(seq1, seq2)))

        # empty sequences after trimming
        seq1 = DNA("AAA---")
        seq2 = DNA("---TTT")
        self.assertTrue(np.isnan(pdist(seq1, seq2)))

        seq1 = Protein("MGCPS")
        seq2 = Protein("XXXXX")
        self.assertTrue(np.isnan(pdist(seq1, seq2)))

        # non-grammared sequences
        seq1 = Sequence("AGCNT")
        seq2 = Sequence("CG-AT")
        with self.assertRaises(TypeError):
            pdist(seq1, seq2)


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


class TestJC69(TestCase):
    def test_jc69(self):
        # normal case
        seq1 = DNA("AGATC")
        seq2 = DNA("AGATG")
        obs = jc69(seq1, seq2)
        self.assertIsInstance(obs, float)
        self.assertEqual(round(obs, 5), 0.23262)

        # RNA sequences
        seq1, seq2 = RNA("AUCG"), RNA("UACG")
        self.assertEqual(round(jc69(seq1, seq2), 5), 0.82396)

        # sequences with gaps
        seq1, seq2 = DNA("A--ACGG"), DNA("AGAAT-G")
        self.assertEqual(round(jc69(seq1, seq2), 5), 0.30410)

        # sequences with degenerate characters
        seq1, seq2 = DNA("ANGCRT"), DNA("CCSMTT")
        self.assertEqual(round(jc69(seq1, seq2), 5), 0.82396)

        # identical sequences
        seq1, seq2 = DNA("ACGT"), DNA("ACGT")
        self.assertEqual(jc69(seq1, seq2), 0.0)

        # distinct sequences
        seq1, seq2 = DNA("ACGT"), DNA("TGCA")
        self.assertTrue(np.isnan(jc69(seq1, seq2)))

        # highly divergent sequences (p = 0.7)
        seq1, seq2 = DNA("ACGAGCTCCT"), DNA("GCTTGAGTCA")
        self.assertEqual(round(jc69(seq1, seq2), 5), 2.03104)

        # overly divergent sequences (p = 0.8)
        seq1, seq2 = DNA("GACTA"), DNA("CTCAG")
        self.assertTrue(np.isnan(jc69(seq1, seq2)))

        # empty sequences
        seq1, seq2 = DNA(""), DNA("")
        self.assertTrue(np.isnan(jc69(seq1, seq2)))

        # protein sequences
        seq1, seq2 = Protein("-PYCRNG"), Protein("MPYAKC-")
        with self.assertRaises(TypeError):
            jc69(seq1, seq2)

        # non-grammared sequences
        seq1, seq2 = Sequence("AGCNT"), Sequence("CG-AT")
        with self.assertRaises(TypeError):
            jc69(seq1, seq2)

    def test_jc69_correct(self):
        # scalar input
        obs = jc69_correct(0.1)
        self.assertIsInstance(obs, float)
        self.assertEqual(round(obs, 5), 0.10733)
        self.assertEqual(round(jc69_correct(0.2), 5), 0.23262)
        self.assertEqual(round(jc69_correct(0.5), 5), 0.82396)
        self.assertEqual(round(jc69_correct(0.7), 5), 2.03104)
        self.assertEqual(jc69_correct(0.0), 0.0)
        self.assertTrue(np.isnan(jc69_correct(0.8)))

        # list input
        lst = [0.0, 0.1, 0.2, 0.5, 0.7, 1.0]
        obs = jc69_correct(lst)
        self.assertIsInstance(obs, np.ndarray)
        self.assertTupleEqual(obs.shape, (6,))
        exp = np.array([0.0, 0.107, 0.233, 0.824, 2.031, np.nan])
        npt.assert_array_equal(obs.round(3), exp)

        # inplace has no effect on non-array
        obs = jc69_correct(lst, inplace=True)
        npt.assert_array_equal(obs.round(3), exp)
        self.assertIsNot(obs, lst)

        # 1D array input
        arr = np.array(lst)
        obs = jc69_correct(arr, inplace=False)
        npt.assert_array_equal(obs.round(3), exp)
        self.assertIsNot(obs, arr)

        # modify array in-place
        obs = jc69_correct(arr, inplace=True)
        npt.assert_array_equal(obs.round(3), exp)
        self.assertIs(obs, arr)
        npt.assert_array_equal(arr.round(3), exp)

        # 2D array input
        shape = (2, 3)
        arr = np.reshape(lst, shape)
        obs = jc69_correct(arr)
        self.assertTupleEqual(obs.shape, shape)
        exp = np.reshape(exp, shape)
        npt.assert_array_equal(obs.round(3), exp)

        # alternative character count
        self.assertEqual(round(jc69_correct(0.5, chars=5), 3), 0.785)
        self.assertEqual(round(jc69_correct(0.8, chars=9), 3), 2.047)

        with self.assertRaisesRegex(ValueError, r"`chars` must be at least 2."):
            jc69_correct(0.5, chars=1)


class TestF81(TestCase):
    def test_f81(self):
        # normal case
        seq1 = DNA("AT-ACGGCGA-C")
        seq2 = DNA("AGAAT--CAACC")
        obs = f81(seq1, seq2)
        self.assertIsInstance(obs, float)
        self.assertEqual(round(obs, 5), 0.53708)

        # even base frequencies (equivalent to JC69)
        obs = f81(seq1, seq2, freqs=(.25, .25, .25, .25))
        self.assertEqual(round(obs, 5), 0.51986)
        exp = jc69(seq1, seq2)
        self.assertAlmostEqual(obs, exp)

        # identical sequences after trimming
        self.assertEqual(f81(DNA("AACGTY"), DNA("WACGTT")), 0.0)

        # empty sequences after trimming
        self.assertTrue(np.isnan(f81(DNA("AAA---"), DNA("---TTT"))))

        # highly divergent sequences
        seq1, seq2 = DNA("ACGAGCTCCT"), DNA("GCTTGAGTCA")
        self.assertEqual(round(f81(seq1, seq2), 5), 2.09101)

        # overly divergent sequences
        seq1, seq2 = DNA("GACTA"), DNA("CTCAG")
        self.assertTrue(np.isnan(f81(seq1, seq2)))

        # RNA sequences
        seq1, seq2 = RNA("AUCU-CGGU"), RNA("AGGUUCA--")
        self.assertEqual(round(f81(seq1, seq2), 5), 0.83539)

        # non-nucleotide sequences
        with self.assertRaises(TypeError):
            f81(Protein("-PYCRNG"), Protein("MPYAKC-"))
        with self.assertRaises(TypeError):
            f81(Sequence("AGCNT"), Sequence("CG-AT"))


class TestK2P(TestCase):
    def test_k2p(self):
        # normal case (8 sites, 2 transitions, 1 transversion)
        seq1 = DNA("AT-ACGGCGA-C")
        seq2 = DNA("AGAAT--CAACC")
        obs = k2p(seq1, seq2)
        self.assertIsInstance(obs, float)
        self.assertEqual(round(obs, 5), 0.56234)

        # identical sequences after trimming
        self.assertEqual(k2p(DNA("AACGTY"), DNA("WACGTT")), 0.0)

        # empty sequences after trimming
        seq1, seq2 = DNA("AAA---"), DNA("---TTT")
        self.assertTrue(np.isnan(k2p(seq1, seq2)))

        # too many transversions (2Q > 1)
        seq1 = DNA("ACGTACGT")
        seq2 = DNA("AGCATATT")
        self.assertTrue(np.isnan(k2p(seq1, seq2)))

        # too many transitions (2P + Q > 1)
        seq1 = DNA("ACGTATGT")
        seq2 = DNA("GTCTACAT")
        self.assertTrue(np.isnan(k2p(seq1, seq2)))

        # RNA sequences
        seq1, seq2 = RNA("AUCU-CGGU"), RNA("AGGUUCA--")
        self.assertEqual(round(k2p(seq1, seq2), 5), 0.82396)

        # non-nucleotide sequences
        with self.assertRaises(TypeError):
            k2p(Protein("-PYCRNG"), Protein("MPYAKC-"))
        with self.assertRaises(TypeError):
            k2p(Sequence("AGCNT"), Sequence("CG-AT"))


class TestF84(TestCase):
    def test_f84(self):
        # normal case
        seq1 = DNA("AT-ACGGCGA-C")
        seq2 = DNA("AGAAT--CAACC")
        obs = f84(seq1, seq2)
        self.assertIsInstance(obs, float)
        self.assertEqual(round(obs, 5), 0.62024)

        # even base frequencies (equivalent to K2P)
        obs = f84(seq1, seq2, freqs=(.25, .25, .25, .25))
        self.assertEqual(round(obs, 5), 0.56234)
        exp = k2p(seq1, seq2)
        self.assertAlmostEqual(obs, exp)

        # identical sequences after trimming
        self.assertEqual(f84(DNA("AACGTY"), DNA("WACGTT")), 0.0)

        # empty sequences after trimming
        seq1, seq2 = DNA("AAA---"), DNA("---TTT")
        self.assertTrue(np.isnan(f84(seq1, seq2)))

        # too many transversions
        seq1 = DNA("ACGTACGT")
        seq2 = DNA("AGCATATT")
        self.assertTrue(np.isnan(f84(seq1, seq2)))

        # too many transitions
        seq1 = DNA("ACGTATGT")
        seq2 = DNA("GTCTACAT")
        self.assertTrue(np.isnan(f84(seq1, seq2)))

        # RNA sequences
        seq1, seq2 = RNA("AUCU-CGGU"), RNA("AGGUUCA--")
        self.assertEqual(round(f84(seq1, seq2), 5), 0.83551)

        # non-nucleotide sequences
        with self.assertRaises(TypeError):
            f84(Protein("-PYCRNG"), Protein("MPYAKC-"))
        with self.assertRaises(TypeError):
            f84(Sequence("AGCNT"), Sequence("CG-AT"))


class TestTN93(TestCase):
    def test_tn93(self):
        # normal case: 8 sites, 1 purine transition, 1 pyrimidine transition,
        # 1 transversion; use observed base frequencies
        seq1 = DNA("AT-ACGGCGA-C")
        seq2 = DNA("AGAAT--CAACC")
        obs = tn93(seq1, seq2)
        self.assertIsInstance(obs, float)
        self.assertEqual(round(obs, 5), 0.99700)

        # specify base frequencies
        obs = tn93(seq1, seq2, freqs=(.2, .2, .3, .3))
        self.assertEqual(round(obs, 5), 0.57303)

        # even base frequencies
        obs = tn93(seq1, seq2, freqs=(.25, .25, .25, .25))
        self.assertEqual(round(obs, 5), 0.56234)

        # zero base frequency
        obs = tn93(seq1, seq2, freqs=(.0, .2, .3, .5))
        self.assertTrue(np.isnan(obs))

        # identical sequences after trimming
        self.assertEqual(tn93(DNA("AACGTY"), DNA("WACGTT")), 0.0)

        # empty sequences after trimming
        seq1, seq2 = DNA("AAA---"), DNA("---TTT")
        self.assertTrue(np.isnan(tn93(seq1, seq2)))

        # too many purine transitions (a1 < 0)
        seq1 = DNA("ACGTATGT")
        seq2 = DNA("GCATGCAT")
        self.assertTrue(np.isnan(tn93(seq1, seq2)))

        # too many pyrimidine transitions (a2 < 0)
        seq1 = DNA("ACCTTGCC")
        seq2 = DNA("ATTTCGCT")
        self.assertTrue(np.isnan(tn93(seq1, seq2)))

        # too many transversions (a3 < 0)
        seq1 = DNA("ACGTACGT")
        seq2 = DNA("AGCATATT")
        self.assertTrue(np.isnan(tn93(seq1, seq2)))

        # RNA sequences
        seq1 = RNA("AUCU-CGCAGU")
        seq2 = RNA("AGGUUCAUA--")
        self.assertEqual(round(tn93(seq1, seq2), 5), 0.88543)

        # non-nucleotide sequences
        with self.assertRaises(TypeError):
            tn93(Protein("-PYCRNG"), Protein("MPYAKC-"))
        with self.assertRaises(TypeError):
            tn93(Sequence("AGCNT"), Sequence("CG-AT"))


if __name__ == "__main__":
    main()
