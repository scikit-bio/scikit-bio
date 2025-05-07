# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

import numpy.testing as npt

from skbio.sequence import DNA, Protein, SubstitutionMatrix
from skbio.alignment import TabularMSA, AlignPath
from skbio.util import get_data_path

from skbio.alignment._score import trim_end_gaps, align_score


class TestScore(unittest.TestCase):

    def test_trim_end_gaps(self):
        aln = [
            "GAGTCCCA",  # no gap
            "--GTTCCA",  # leading gap
            "CAATT---",  # trailing gap
            "-AGGTCC-",  # both gaps
            "-AG--CC-",  # with internal gap
            "--------",  # gap-only
        ]
        obs = trim_end_gaps(aln)
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

        # Typical setting: match/mismatch with affine gap penalty.
        # Parameters were adopted from NCBI BLASTN:
        # match = 2, mismatch = -3, gap_open = 5, gap_extend = 2
        obs = align_score(aln, (2, -3), (5, 2))
        self.assertEqual(obs, -28.0)

        # Manual calculation:
        # 0: CGGTCGTAACGCGTA---CA
        # 1: CAG--GTAAG-CATACCTCA
        # => 2 * 11 - 3 * 3 - 5 * 3 - 2 * 6 = -14
        # 0: CGGTCGTAACGCGTA---|CA
        # 2: CGGTCGTCAC-TGTACAC|xx
        # => 2 * 12 - 3 * 2 - 5 * 2 - 2 * 4 = 0
        # 1: CAG--GTAAG|x|CATACCT|CA
        # 2: CGGTCGTCAC|x|TGTACAC|xx
        # => 2 *  8 - 3 * 7 - 5 * 1 - 2 * 2 = -14
        # total => -14 + 0 -14 = -28

        # penalize terminal gaps
        obs = align_score(aln, (2, -3), (5, 2), free_ends=False)
        self.assertEqual(obs, -46.0)

        # Alternative parameters: adopted from EMBOSS Needle:
        # match = 5, mismatch = -4, gap_open = 10, gap_extend = 0.5
        # Because EMBOSS doesn't apply gap_extend to the first gap position, we
        # subtract gap_extend from gap_open here to mimic its behavior.
        # 5 * 11 - 4 * 3 - 9.5 * 3 - 0.5 * 6 = 11.5
        # 5 * 12 - 4 * 2 - 9.5 * 2 - 0.5 * 4 = 31
        # 5 *  8 - 4 * 7 - 9.5 * 1 - 0.5 * 2 = 1.5
        # total: 11.5 + 31 + 1.5 = 44
        obs = align_score(aln, (5, -4), (9.5, 0.5))
        self.assertEqual(obs, 44.0)

        # The match/mismatch scores are from substitution matrix "NUC.4.4".
        submat = SubstitutionMatrix.by_name("NUC.4.4")
        obs = align_score(aln, submat, (9.5, 0.5))
        self.assertEqual(obs, 44.0)

        # substitution matrix by name
        obs = align_score(aln, "NUC.4.4", (9.5, 0.5))
        self.assertEqual(obs, 44.0)

        # linear gap penalty
        # 1 * (11 + 12 + 8) - 2 * (6 + 4 + 2) = 7
        obs = align_score(aln, (1, 0), 2)
        self.assertEqual(obs, 7.0)

    def test_align_score_input(self):
        """Test on various input formats."""
        aln = ["CGGTCGTAACGCGTA---CA",
               "CAG--GTAAG-CATACCTCA",
               "CGGTCGTCAC-TGTACAC--"]

        # list of raw strings
        obs = align_score(aln, (5, -4), (5, 2))
        self.assertEqual(obs, 53)

        # list of Sequence objects
        seqs = list(map(DNA, aln))
        obs = align_score(seqs, (5, -4), (5, 2))
        self.assertEqual(obs, 53)

        # custom gap_chars
        obs = align_score(seqs, (5, -4), (5, 2), gap_chars="A")
        self.assertEqual(obs, -8)

        # which is equivalent to (- becomes X and A becomes -):
        aln2 = [
            "CGGTCGT--CGCGT-XXXC-",
            "C-GXXGT--GXC-T-CCTC-",
            "CGGTCGTC-CXTGT-C-CXX",
        ]
        obs = align_score(aln2, (5, -4), (5, 2))
        self.assertEqual(obs, -8)

        # TabularMSA object
        msa = TabularMSA(map(DNA, aln))
        obs = align_score(msa, (5, -4), (5, 2))
        self.assertEqual(obs, 53)

        # gap_chars doesn't impact TabularMSA
        obs = align_score(msa, (5, -4), (5, 2), gap_chars="A")
        self.assertEqual(obs, 53)

        # alignment path
        path = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                         states=[0, 2, 0, 6, 0, 1, 4],
                         starts=[0, 3, 0])
        seqs = ["CGGTCGTAACGCGTACA",
                "GGGCAGGTAAGCATACCTCA",
                "CGGTCGTCACTGTACACAAA"]
        obs = align_score((path, seqs), (5, -4), (5, 2))
        self.assertEqual(obs, 53)

        # gap_chars doesn't impact AlignPath
        obs = align_score((path, seqs), (5, -4), (5, 2), gap_chars="A")
        self.assertEqual(obs, 53)

    def test_align_score_pair(self):
        """Test on pairwise alignments."""
        # single internal gap, no terminal gap
        aln = ["CCTCAT-C",
               "CGTCGTGC"]
        obs = align_score(aln, (5, -4), (5, 2))
        self.assertEqual(obs, 10)
        obs = align_score(aln, "NUC.4.4", (5, 2))
        self.assertEqual(obs, 10)
        obs = align_score(aln, (5, -4), (5, 2), free_ends=False)
        self.assertEqual(obs, 10)
        obs = align_score(aln, (5, -4), 3)
        self.assertEqual(obs, 14)

        # continuous internal gap, no terminal gap
        aln = ["ACGT---ACGT",
               "ACGTCCCACGT"]
        obs = align_score(aln, (5, -4), (5, 2))
        self.assertEqual(obs, 29)

        # with terminal gaps
        aln = ["CGTCAT--",
               "--TCATGC"]
        obs = align_score(aln, (5, -4), (5, 2))
        self.assertEqual(obs, 20)
        obs = align_score(aln, (5, -4), (5, 2), free_ends=False)
        self.assertEqual(obs, 2)

        # very complex gaps
        aln = ["--AA-AAA-----CC--TCAT--G------",
               "----C-GG--G--C-TC--GT---CCA---"]
        obs = align_score(aln, (5, -4), (5, 2))
        self.assertEqual(obs, -55)
        obs = align_score(aln, (5, -4), (5, 2), free_ends=False)
        self.assertEqual(obs, -75)

    def test_align_score_real(self):
        """Test on real-world datasets."""
        # protein
        aln = TabularMSA.read(get_data_path("il6.prot.aln"), constructor=Protein)
        exp = 11044
        obs = align_score(aln, "BLOSUM62", (5, 2))
        self.assertEqual(obs, exp)
        exp = 10918
        obs = align_score(aln, "BLOSUM62", (5, 2), free_ends=False)
        self.assertEqual(obs, exp)

        # DNA
        aln = TabularMSA.read(get_data_path("il6.nucl.aln"), constructor=DNA)
        exp = 26550
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

        # one pair doesn't align
        aln = ["AAAAA", "AA---", "---AA"]
        obs = align_score(aln, (5, -4), (5, 2), free_ends=False)
        self.assertEqual(obs, -20)
        obs = align_score(aln, (5, -4), (5, 2), free_ends=True)
        self.assertEqual(obs, 20)

    def test_align_score_error(self):
        """Test on inputs that cause errors."""
        msg = "Sequences must be strings or Sequence objects."
        with self.assertRaises(ValueError) as cm:
            align_score([1, None, 0.5], (5, -4), (5, 2))
        self.assertEqual(str(cm.exception), msg)

        msg = "Alignment must contain at least two sequences."
        with self.assertRaises(ValueError) as cm:
            align_score([], (5, -4), (5, 2))
        self.assertEqual(str(cm.exception), msg)

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

        msg = "The alignment has a length of zero."
        with self.assertRaises(ValueError) as cm:
            align_score([""] * 3, (5, -4), (5, 2))
        self.assertEqual(str(cm.exception), msg)

        aln = ["MKQ-PSV", "MKIDTS-"]
        obs = align_score(aln, "BLOSUM62", (5, 2))
        self.assertEqual(obs, 3)
        msg = ("Sequences contain characters that are not present in the provided "
               "substitution matrix.")
        with self.assertRaises(ValueError) as cm:
            align_score(aln, "NUC.4.4", (5, 2))
        self.assertEqual(str(cm.exception), msg)

        aln = [chr(0) + chr(12) + chr(34) + chr(56) + chr(78) + chr(90),
               chr(123) + chr(98) + chr(76) + "-" + chr(54) + chr(32)]
        with self.assertRaises(ValueError) as cm:
            align_score(aln, "BLOSUM62", (5, 2))
        self.assertEqual(str(cm.exception), msg)


if __name__ == "__main__":
    unittest.main()
