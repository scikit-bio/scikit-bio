# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main

import numpy as np
import numpy.testing as npt

from skbio.alignment import TabularMSA
from skbio.sequence import DNA, RNA, Protein
from skbio.sequence.distance import (
    hamming, pdist, jc69, f81, k2p, f84, tn93, logdet, paralin, jc69_correct
)
from skbio.stats.distance import DistanceMatrix
from skbio.util import get_data_path

from skbio.alignment._distance import align_dists, _char_hash


nan = np.nan

class AlignDistsTests(TestCase):
    def setUp(self):
        # a DNA sequence alignment
        self.msa_nucl = TabularMSA([
            DNA("CTGGTG---CCTGCTGCCTTCCCTGCCCCAGTACCCGGAGAAG--TCCAAAGATGTA"),
            DNA("CTGGTG---ACCACGGCC------ACTTCACAAGT----GGAGACTTCACAGAGGAT"),
            DNA("CTG-CGCCCGCCG-CGCCGTCC--CTGCCCGCCGCCGCGGACTCGTCCG----GGTT"),
        ])

        # a protein sequence alignment
        self.msa_prot = TabularMSA([
            Protein("MLGLLLVLPAAF-APVPPGEDMF"),
            Protein("MLPLLLVMATAFDTPGRL-EDML"),
            Protein("--GLMLQCTTAFPTSQVRRGDIL"),
        ])

        # a DNA sequence alignment with ambiguous codes
        self.msa1 = TabularMSA([
            DNA("ATC-GTATCGY"),  # degenerate char: "Y"
            DNA("ATGCG--CCGC"),
            DNA("GT-CGTACG-T"),
            DNA("GT-NTTACAGT"),  # degenerate char: "N"
        ], index=list("abcd"))

        # a protein sequence alignment with ambiguous codes
        self.msa2 = TabularMSA([
            Protein("MQVEGATSI"),
            Protein("MKNJ-PTSL"),  # degenerate char: "J"
            Protein("MKDE-PYRX"),  # degenerate char: "X"
            Protein("MQUENPT--"),  # non-canonical char: "U"
        ], index=list("abcd"))

        # an RNA sequence alignment where no site is shared across all sequences.
        self.msa3 = TabularMSA([
            RNA("AG-UGCU-CA"),
            RNA("--AUC--GC-"),
            RNA("ACU--AUG-A"),
        ], index=list("abc"))

        # a DNA sequence alignment where one pair of sequences shares no site.
        self.msa4 = TabularMSA([
            DNA("ACGT-CT"),
            DNA("AG--A-T"),
            DNA("--CT-C-"),
        ], index=list("abc"))

        # Identical sequences:
        # 0 and 1 are identical
        # 1 and 2 are identical after complete/pairwise deletion
        # 2 and 3 are identical after complete deletion
        self.msa_ident = TabularMSA([
            DNA("ACG-TCTGCC"),
            DNA("ACG-TCTGCC"),
            DNA("ACGGTCT-CC"),
            DNA("ACGATCT-CC"),
        ], index=list("abcd"))

    def test_align_dists_error(self):
        # wrong sequence type
        with self.assertRaises(TypeError) as cm:
            align_dists(self.msa_prot, "jc69")
        msg = "'jc69' is compatible with ('DNA', 'RNA') sequences, not 'Protein'."
        self.assertEqual(str(cm.exception), msg)

        # wrong metric type
        with self.assertRaises(TypeError) as cm:
            align_dists(self.msa_prot, 123)
        msg = "`metric` must be a function or a string."
        self.assertEqual(str(cm.exception), msg)

        # wrong metric name
        with self.assertRaises(ValueError) as cm:
            align_dists(self.msa_prot, "hello")
        msg = ("'hello' is not an available sequence distance metric name. " 
               "Refer to `skbio.sequence.distance` for a list of available metrics.")
        self.assertEqual(str(cm.exception), msg)

        # is a preset function but not a distance metric
        with self.assertRaises(ValueError) as cm:
            align_dists(self.msa_prot, jc69_correct)
        msg = ("`jc69_correct` is not a sequence distance metric. Refer to "
               "`skbio.sequence.distance` for a list of available metrics.")
        self.assertEqual(str(cm.exception), msg)

        # no sequence
        with self.assertRaises(ValueError) as cm:
            align_dists(TabularMSA(()), "pdist")
        msg = "Alignment contains no sequence."
        self.assertEqual(str(cm.exception), msg)

    def test_align_dists_edge(self):
        # one sequence
        msa = TabularMSA((DNA("ACGT"),), index=["a"])
        obs = align_dists(msa, "pdist")
        exp = np.array([])
        npt.assert_array_equal(obs.condensed_form(), exp)
        self.assertTupleEqual(obs.ids, ("a",))

        # zero-length sequences
        msa = TabularMSA((DNA(""), DNA(""), DNA("")), index=list("abc"))
        obs = align_dists(msa, "pdist")
        exp = np.array([nan, nan, nan])
        npt.assert_array_equal(obs.condensed_form(), exp)
        self.assertTupleEqual(obs.ids, ("a", "b", "c"))

    def test_align_dists_custom(self):

        # custom function
        def metric1(seq1, seq2):
            return seq1.count("A") + seq2.count("A")

        obs = align_dists(self.msa1, metric1)
        exp = np.array([2, 1, 2, 1, 2, 1], dtype=float)
        npt.assert_array_equal(obs.condensed_form(), exp)

        obs = align_dists(self.msa1, metric1, shared_by_all=False)
        exp = np.array([2, 3, 4, 1, 2, 3], dtype=float)
        npt.assert_array_equal(obs.condensed_form(), exp)

    def test_align_dists_hamming(self):
        """Hamming distance."""
        obs = align_dists(self.msa_nucl, "hamming")
        exp = np.array([0.5087719, 0.5614035, 0.6666667])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "hamming", shared_by_all=False)
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, hamming)
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, hamming, shared_by_all=False)
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "hamming", proportion=False)
        exp = np.array([29, 32, 38], dtype=float)
        npt.assert_array_equal(obs.condensed_form(), exp)

    def test_align_dists_pdist(self):
        """p-distance."""
        obs = align_dists(self.msa_nucl, "pdist")
        exp = np.array([0.4444444, 0.4444444, 0.5277778])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "pdist", shared_by_all=False)
        exp = np.array([0.4047619, 0.4318182, 0.5526316])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        # complete deletion (5 sites)
        obs = align_dists(self.msa1, "pdist")
        self.assertIsInstance(obs, DistanceMatrix)
        self.assertTupleEqual(obs.ids, tuple("abcd"))
        exp = np.array([0.2, 0.6, 0.8, 0.4, 0.6, 0.4])
        npt.assert_array_equal(obs.condensed_form(), exp)

        # supply function
        obs = align_dists(self.msa1, pdist)
        npt.assert_array_equal(obs.condensed_form(), exp)

        # pairwise deletion
        obs = align_dists(self.msa1, "pdist", shared_by_all=False)
        exp = np.array([2/7, 3/7, 4/8, 3/7, 4/7, 2/8])
        npt.assert_allclose(obs.condensed_form(), exp)

        # supply function
        obs = align_dists(self.msa1, pdist, shared_by_all=False)
        npt.assert_array_equal(obs.condensed_form(), exp)

        # edge case: no site left after complete deletion
        obs = align_dists(self.msa3, "pdist")
        self.assertTrue(np.isnan(obs.condensed_form()).all())

        obs = align_dists(self.msa3, "pdist", shared_by_all=False)
        exp = np.array([0.333, 0.4, 0.5])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        # edge case: no site left after pairwise deletion
        obs = align_dists(self.msa4, "pdist")
        self.assertTrue(np.isnan(obs.condensed_form()).all())

        obs = align_dists(self.msa4, "pdist", shared_by_all=False)
        exp = np.array([0.333, 0.333, nan])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

    def test_align_dists_jc69(self):
        """JC69 distance."""
        obs = align_dists(self.msa_nucl, "jc69")
        exp = np.array([0.6734562, 0.6734562, 0.9122965])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "jc69", shared_by_all=False)
        exp = np.array([0.5818792, 0.6430877, 1.0012508])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa1, "jc69")
        exp = np.array([0.233, 1.207, nan, 0.572, 1.207, 0.572])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(self.msa1, jc69)
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(self.msa1, "jc69", shared_by_all=False)
        exp = np.array([0.360, 0.635, 0.824, 0.635, 1.076, 0.304])
        npt.assert_allclose(obs.condensed_form().round(3), exp)

        obs = align_dists(self.msa1, jc69, shared_by_all=False)
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(self.msa3, "jc69")
        self.assertTrue(np.isnan(obs.condensed_form()).all())

        obs = align_dists(self.msa3, "jc69", shared_by_all=False)
        exp = np.array([0.441, 0.572, 0.824])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        msg = ("'jc69' is compatible with ('DNA', 'RNA') sequences, not 'Protein'.")
        with self.assertRaises(TypeError) as cm:
            align_dists(self.msa2, "jc69")
        self.assertEqual(str(cm.exception), msg)
        with self.assertRaises(TypeError) as cm:
            align_dists(self.msa2, jc69)
        self.assertEqual(str(cm.exception), msg)

    def test_align_dists_f81(self):
        """F81 distance."""
        # normal DNA sequences
        obs = align_dists(self.msa_nucl, "f81")
        exp = np.array([0.6922797, 0.6922797, 0.9524649])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, f81)
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "f81", shared_by_all=False)
        exp = np.array([0.5951889, 0.6599422, 1.0526322])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "f81", freqs=np.full(4, .25))
        exp = np.array([0.6734562, 0.6734562, 0.9122965])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        # F81 is equivalent to JC69 when base frequencies are equal.
        exp = align_dists(self.msa_nucl, "jc69")
        npt.assert_allclose(obs.condensed_form(), exp.condensed_form())

        # short DNA sequences with ambiguity
        obs = align_dists(self.msa1, "f81")
        exp = np.array([0.233, 1.234, nan, 0.576, 1.234, 0.576])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(self.msa1, "f81", shared_by_all=False)
        exp = np.array([0.361, 0.641, 0.834, 0.641, 1.096, 0.305])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        # short, hollow RNA sequences
        obs = align_dists(self.msa3, "f81")
        self.assertTrue(np.isnan(obs.condensed_form()).all())

        obs = align_dists(self.msa3, "f81", shared_by_all=False)
        exp = np.array([0.442, 0.574, 0.829])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        with self.assertRaises(TypeError):
            align_dists(self.msa2, "f81")

    def test_align_dists_k2p(self):
        """K2P distance."""
        # normal DNA sequences
        obs = align_dists(self.msa_nucl, "k2p")
        exp = np.array([0.6962528, 0.6749634, 0.9547713])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, k2p)
        exp = np.array([0.6962528, 0.6749634, 0.9547713])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "k2p", shared_by_all=False)
        exp = np.array([0.5921321, 0.6440233, 1.0287045])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        # short DNA sequences with ambiguity
        obs = align_dists(self.msa1, "k2p")
        exp = np.array([0.255, nan, nan, 0.586, 1.207, 0.586])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(self.msa1, "k2p", shared_by_all=False)
        exp = np.array([0.364, 0.710, 0.866, 0.710, 1.185, 0.307])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        # short, hollow RNA sequences
        obs = align_dists(self.msa3, "k2p")
        self.assertTrue(np.isnan(obs.condensed_form()).all())

        obs = align_dists(self.msa3, "k2p", shared_by_all=False)
        exp = np.array([0.477, 0.658, nan])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        with self.assertRaises(TypeError):
            align_dists(self.msa2, "k2p")

    def test_align_dists_f84(self):
        """F84 distance."""
        # normal DNA sequences
        obs = align_dists(self.msa_nucl, "f84")
        exp = np.array([0.7322711, 0.6966819, 1.0393658])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, f84)
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "f84", shared_by_all=False)
        exp = np.array([0.6132254, 0.6629780, 1.1168612])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "f84", freqs=np.full(4, .25))
        exp = np.array([0.6962528, 0.6749634, 0.9547713])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        # F84 is equivalent to K2P when base frequencies are equal.
        exp = align_dists(self.msa_nucl, "k2p")
        npt.assert_allclose(obs.condensed_form(), exp.condensed_form())

        # short DNA sequences with ambiguity
        obs = align_dists(self.msa1, "f84")
        exp = np.array([0.258, nan, nan, 0.592, 1.234, 0.592])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(self.msa1, "f84", shared_by_all=False)
        exp = np.array([0.366, 0.727, 0.884, 0.727, 1.230, 0.308])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        # short, hollow RNA sequences
        obs = align_dists(self.msa3, "f84")
        self.assertTrue(np.isnan(obs.condensed_form()).all())

        obs = align_dists(self.msa3, "f84", shared_by_all=False)
        exp = np.array([0.479, 0.661, nan])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        with self.assertRaises(TypeError):
            align_dists(self.msa2, "f84")

    def test_align_dists_tn93(self):
        """TN93 distance."""
        obs = align_dists(self.msa_nucl, "tn93")
        exp = np.array([0.7404173, 0.7636621, 1.0432672])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, tn93)
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "tn93", shared_by_all=False)
        exp = np.array([0.6164222, 0.6757994, 1.1236656])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "tn93", freqs=np.full(4, .25))
        exp = np.array([0.7256986, 0.7737915, 0.9709059])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa1, "tn93")
        exp = np.array([0.357, nan, nan, nan, nan, nan])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(self.msa3, "tn93")
        self.assertTrue(np.isnan(obs.condensed_form()).all())

        obs = align_dists(self.msa3, "tn93", shared_by_all=False)
        exp = np.array([0.479, 0.661, nan])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

    def test_align_dists_logdet(self):
        """LogDet distance."""
        obs = align_dists(self.msa_nucl, "logdet")
        exp = np.array([0.7422038, 0.9347606, nan])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, logdet)
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "logdet", shared_by_all=False)
        exp = np.array([0.6549783, 0.9558150, nan])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        aln = TabularMSA.read(get_data_path("tp53.nucl.aln"), constructor=DNA)
        obs = align_dists(aln, "logdet")
        exp = np.array([
            0.049, 0.242, 0.195, 0.171, 0.253, 0.204, 0.175, 0.289, 0.248, 0.219])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(aln, "logdet", shared_by_all=False)
        exp = np.array([
            0.051, 0.250, 0.199, 0.177, 0.262, 0.208, 0.181, 0.291, 0.247, 0.219])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        aln = TabularMSA.read(get_data_path("tp53.prot.aln"), constructor=Protein)
        obs = align_dists(aln, "logdet")
        exp = np.array([
            0.189, 0.393, 0.342, 0.277, 0.406, 0.348, 0.267, 0.450, 0.374, 0.347])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(aln, "logdet", pseudocount=0.5)
        exp = np.array([
            0.585, 0.787, 0.737, 0.672, 0.800, 0.743, 0.663, 0.844, 0.768, 0.741])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(self.msa4, "logdet", shared_by_all=False)
        self.assertTrue(np.isnan(obs.condensed_form()).all())

    def test_align_dists_paralin(self):
        """Paralinear distance."""
        obs = align_dists(self.msa_nucl, "paralin")
        exp = np.array([0.7212106, 0.7054428, nan])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, paralin)
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        obs = align_dists(self.msa_nucl, "paralin", shared_by_all=False)
        exp = np.array([0.6378228, 0.6976025, nan])
        npt.assert_array_equal(obs.condensed_form().round(7), exp)

        aln = TabularMSA.read(get_data_path("tp53.nucl.aln"), constructor=DNA)
        obs = align_dists(aln, "paralin")
        exp = np.array([
            0.035, 0.231, 0.181, 0.151, 0.242, 0.190, 0.154, 0.278, 0.230, 0.199])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(aln, "paralin", shared_by_all=False)
        exp = np.array([
            0.037, 0.238, 0.186, 0.158, 0.25 , 0.195, 0.161, 0.279, 0.229, 0.199])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        aln = TabularMSA.read(get_data_path("tp53.prot.aln"), constructor=Protein)
        obs = align_dists(aln, "paralin")
        exp = np.array([
            0.048, 0.244, 0.19 , 0.127, 0.258, 0.197, 0.119, 0.291, 0.218, 0.188])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(aln, "paralin", pseudocount=0.5)
        exp = np.array([
            0.531, 0.732, 0.678, 0.615, 0.745, 0.685, 0.606, 0.784, 0.710, 0.681])
        npt.assert_array_equal(obs.condensed_form().round(3), exp)

        obs = align_dists(self.msa4, "paralin", shared_by_all=False)
        self.assertTrue(np.isnan(obs.condensed_form()).all())

    def test_site_mask(self):
        # NOTE: This test complements the original test of `_char_hash` in skbio.
        # sequence.distance, as it demonstrates the outcome of site masking.
        # DNA sequences
        seqs = np.vstack([seq._bytes for seq in self.msa1])

        obs = _char_hash("ACGT", DNA)[seqs]
        exp = np.array([[1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0],
                        [1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1],
                        [1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                        [1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        obs = _char_hash("ACGTN", DNA)[seqs]
        exp = np.array([[1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0],
                        [1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1],
                        [1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                        [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        obs = _char_hash("ABCD", DNA)[seqs]
        exp = np.array([[1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0],
                        [1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1],
                        [0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        obs = _char_hash("nongap", DNA)[seqs]
        exp = np.array([[1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1],
                        [1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                        [1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        obs = _char_hash("definite", DNA)[seqs]
        exp = np.array([[1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0],
                        [1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1],
                        [1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                        [1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        # protein sequences
        seqs = np.vstack([seq._bytes for seq in self.msa2])

        obs = _char_hash("nongap", Protein)[seqs]
        exp = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 0, 1, 1, 1, 1],
                        [1, 1, 1, 1, 0, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 0, 0]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        obs = _char_hash("definite", Protein)[seqs]
        exp = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 0, 0, 1, 1, 1, 1],
                        [1, 1, 1, 1, 0, 1, 1, 1, 0],
                        [1, 1, 1, 1, 1, 1, 1, 0, 0]], dtype=bool)
        npt.assert_array_equal(obs, exp)

        obs = _char_hash("canonical", Protein)[seqs]
        exp = np.array([[1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 0, 0, 1, 1, 1, 1],
                        [1, 1, 1, 1, 0, 1, 1, 1, 0],
                        [1, 1, 0, 1, 1, 1, 1, 0, 0]], dtype=bool)
        npt.assert_array_equal(obs, exp)


if __name__ == "__main__":
    main()
