# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy as np
import numpy.testing as npt

from skbio.sequence import DNA, Protein, SubstitutionMatrix
from skbio.alignment import TabularMSA
from skbio.util import get_data_path

from skbio.alignment._score import trim_terminal_gaps, align_score


class TestScore(unittest.TestCase):

    def test_trim_terminal_gaps(self):
        aln = [
            "GAGTCCCA",  # no gap
            "--GTTCCA",  # leading gap
            "CAATT---",  # trailing gap
            "-AGGTCC-",  # both gaps
            "-AG--CC-",  # with internal gap
            "--------",  # gap-only
        ]
        obs = trim_terminal_gaps(aln)
        exp = ([0, 2, 0, 1, 1, 0], [8, 8, 5, 7, 7, 0])
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

    def test_align_score_basic(self):
        """Test on a basic case."""
        # This is a complex example, with terminal gap, continuous internal gap, dual
        # gaps, and terminal gap swapping sequences.
        aln = ["CGGTCGTAACGCGTA---CA",
               "CAG--GTAAG-CATACCTCA",
               "CGGTCGTCAC-TGTACAC--"]
        
        # manual calculation:

        # 0: CGGTCGTAACGCGTA---CA
        # 1: CAG--GTAAG-CATACCTCA
        # 5 * 11 - 4 * 3 - 5 * 3 - 2 * 3 = 22

        # 0: CGGTCGTAACGCGTA---|CA
        # 2: CGGTCGTCAC-TGTACAC|--
        # 5 * 12 - 4 * 2 - 5 * 2 - 2 * 2 = 38

        # 1: CAG--GTAAG-CATACCT|CA
        # 2: CGGTCGTCAC-TGTACAC|--
        # 5 * 8 - 4 * 7 - 5 * 1 - 2 * 1 = 5

        # total: 22 + 38 + 5 = 65

        # typical setting: match/mismatch with affine gap penalty
        obs = align_score(aln, (5, -4), (5, 2))
        self.assertEqual(obs, 65)

        # linear gap penalty
        obs = align_score(aln, (5, -4), 10)
        self.assertEqual(obs, -13)

        # substitution matrix
        submat = SubstitutionMatrix.by_name("NUC.4.4")
        obs = align_score(aln, submat, (5, 2))
        self.assertEqual(obs, 65)

        # substitution matrix by name
        obs = align_score(aln, "NUC.4.4", (5, 2))
        self.assertEqual(obs, 65)

        # penalize terminal gaps
        obs = align_score(aln, (5, -4), (5, 2), True)
        self.assertEqual(obs, 51)

    def test_align_score_input(self):
        """Test on various input formats."""
        aln = ["CGGTCGTAACGCGTA---CA",
               "CAG--GTAAG-CATACCTCA",
               "CGGTCGTCAC-TGTACAC--"]

        # list of raw strings
        obs = align_score(aln, (5, -4), (5, 2))
        self.assertEqual(obs, 65)

        # list of Sequence objects
        seqs = list(map(DNA, aln))
        obs = align_score(seqs, (5, -4), (5, 2))
        self.assertEqual(obs, 65)

        # custom gap_chars
        obs = align_score(seqs, (5, -4), (5, 2), gap_chars="A")
        self.assertEqual(obs, 8)

        # which is equivalent to (- becomes X and A becomes -):
        aln2 = [
            "CGGTCGT--CGCGT-XXXC-",
            "C-GXXGT--GXC-T-CCTC-",
            "CGGTCGTC-CXTGT-C-CXX",
        ]
        obs = align_score(aln2, (5, -4), (5, 2))
        self.assertEqual(obs, 8)

        # TabularMSA object
        msa = TabularMSA(map(DNA, aln))
        obs = align_score(msa, (5, -4), (5, 2))
        self.assertEqual(obs, 65)

        # gap_chars doesn't impact TabularMSA
        obs = align_score(msa, (5, -4), (5, 2), gap_chars="A")
        self.assertEqual(obs, 65)

    def test_align_score_pair(self):
        """Test on pairwise alignments."""
        # single internal gap, no terminal gap
        aln = ["CCTCAT-C",
               "CGTCGTGC"]
        obs = align_score(aln, (5, -4), (5, 2))
        self.assertEqual(obs, 12)
        obs = align_score(aln, "NUC.4.4", (5, 2))
        self.assertEqual(obs, 12)
        obs = align_score(aln, (5, -4), (5, 2), terminal_gaps=True)
        self.assertEqual(obs, 12)
        obs = align_score(aln, (5, -4), 3)
        self.assertEqual(obs, 14)

        # continuous internal gap, no terminal gap
        aln = ["ACGT---ACGT",
               "ACGTCCCACGT"]
        obs = align_score(aln, (5, -4), (5, 2))
        self.assertEqual(obs, 31)

        # with terminal gaps
        aln = ["CGTCAT--",
               "--TCATGC"]
        obs = align_score(aln, (5, -4), (5, 2))
        self.assertEqual(obs, 20)
        obs = align_score(aln, (5, -4), (5, 2), terminal_gaps=True)
        self.assertEqual(obs, 6)

        # very complex gaps
        aln = ["--AA-AAA-----CC--TCAT--G------",
               "----C-GG--G--C-TC--GT---CCA---"]
        obs = align_score(aln, (5, -4), (5, 2))
        self.assertEqual(obs, -41)
        obs = align_score(aln, (5, -4), (5, 2), terminal_gaps=True)
        self.assertEqual(obs, -57)

    def test_align_score_real(self):
        """Test on real-world datasets."""
        # protein
        aln = TabularMSA.read(get_data_path("il6.prot.aln"), constructor=Protein)
        exp = 11144
        obs = align_score(aln, "BLOSUM62", (5, 2))
        self.assertEqual(obs, exp)
        exp = 11030
        obs = align_score(aln, "BLOSUM62", (5, 2), terminal_gaps=True)
        self.assertEqual(obs, exp)

        # DNA
        aln = TabularMSA.read(get_data_path("il6.nucl.aln"), constructor=DNA)
        exp = 26658
        obs = align_score(aln, (5, -4), (5, 2))
        self.assertEqual(obs, exp)
        obs = align_score(aln, "NUC.4.4", (5, 2))
        self.assertEqual(obs, exp)

    def test_align_score_edge(self):
        """Test on edge cases."""
        # single character, same across sequences
        obs = align_score(["A"] * 3, (5, -4), (5, 2))
        self.assertEqual(obs, 15)

        # single character, different between sequences
        obs = align_score(["A", "B", "C"], (5, -4), (5, 2))
        self.assertEqual(obs, -12)

        # all gaps are terminal
        obs = align_score(["A--", "-A-", "--A"], (5, -4), (5, 2))
        self.assertEqual(obs, 0)

    def test_align_score_error(self):
        """Test on inputs that cause errors."""
        msg = "There is no sequence in the alignment."
        with self.assertRaises(ValueError) as cm:
            align_score([], (5, -4), (5, 2))
        self.assertEqual(str(cm.exception), msg)

        msg = "There is only one sequence in the alignment."
        with self.assertRaises(ValueError) as cm:
            align_score(["ABC"], (5, -4), (5, 2))
        self.assertEqual(str(cm.exception), msg)

        msg = "The alignment contains gap-only sequence(s)."
        with self.assertRaises(ValueError) as cm:
            align_score(["ABC", "---"], (5, -4), (5, 2))
        self.assertEqual(str(cm.exception), msg)

        msg = "Sequence lengths do not match."
        with self.assertRaises(ValueError) as cm:
            align_score(["A", "AA", "AAA"], (5, -4), (5, 2))
        self.assertEqual(str(cm.exception), msg)

        msg = "The alignment has a length of 0."
        with self.assertRaises(ValueError) as cm:
            align_score([""] * 3, (5, -4), (5, 2))
        self.assertEqual(str(cm.exception), msg)

        aln = ["MKQ-PSV", "MKIDTS-"]
        obs = align_score(aln, "BLOSUM62", (5, 2))
        self.assertEqual(obs, 5)
        msg = ("Sequences contain characters that are not present in the provided "
               "substitution matrix.")
        with self.assertRaises(ValueError) as cm:
            align_score(aln, "NUC.4.4", (5, 2))
        self.assertEqual(str(cm.exception), msg)


if __name__ == "__main__":
    unittest.main()
