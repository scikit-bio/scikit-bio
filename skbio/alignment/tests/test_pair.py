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

import skbio.io
from skbio.util import get_data_path
from skbio.alignment import TabularMSA
from skbio.sequence import Sequence, GrammaredSequence, DNA, Protein
from skbio.sequence import SubstitutionMatrix
from skbio.util import classproperty
from skbio.util._decorator import overrides
from skbio.alignment._utils import encode_sequences

from skbio.alignment._pair import (
    pair_align,
    pair_align_nucl,
    pair_align_prot,
    _prep_free_ends,
    _alloc_matrices,
    _init_matrices,
    _one_stop,
    _all_stops,
    _traceback_one,
    _traceback_all,
)
from skbio.alignment._cutils import (
    _fill_linear_matrix,
    _fill_affine_matrices,
)


class PairAlignTests(unittest.TestCase):
    def test_pair_align_nucl(self):
        """Walk-through of a nucleotide alignment case."""
        # Recognition sites of restriction enzymes EcoRI (GAATTC) and BglII (AGATCT).
        seq1 = DNA("GAATTC")
        seq2 = DNA("AGATCT")

        # Default behavior is semi-global alignment (free ends) with match score = 1,
        # mismatch score = -1 and linear gap penalty = 2.
        res = pair_align(seq1, seq2)

        # Up to one alignment path is returned.
        self.assertEqual(len(res.paths), 1)

        # Print the aligned sequences.
        path = res.paths[0]
        obs = path.to_aligned((seq1, seq2))
        exp = ["-GAATTC",
               "AGATCT-"]
        self.assertListEqual(obs, exp)

        # A path can be represented by a CIGAR string.
        self.assertEqual(path.to_cigar(), "1I5M1D")

        # There are 3 matches and 2 mismatches, so score = 3 - 2 = 1.
        self.assertEqual(res.score, 1)

        # By default, the alignment matrix is discarded.
        self.assertIsNone(res.matrices)

        # Return all optimal alignments. There are two.
        res = pair_align(seq1, seq2, max_paths=None)
        self.assertEqual(len(res.paths), 2)

        # The 1st path is identical to the previous one. But the 2nd is different.
        obs = res.paths[0].to_aligned((seq1, seq2))
        self.assertListEqual(obs, exp)
        obs = res.paths[1].to_aligned((seq1, seq2))
        exp = ["GAATTC-",
               "-AGATCT"]
        self.assertListEqual(obs, exp)
        self.assertEqual(res.score, 1)

        # Disable traceback.
        res = pair_align(seq1, seq2, max_paths=0)
        self.assertIsNone(res.paths)

        # Keep the alignment matrix.
        res = pair_align(seq1, seq2, keep_matrices=True)
        self.assertEqual(len(res.matrices), 1)
        obs = res.matrices[0]
        exp = np.array([[ 0,  0,  0,  0,  0,  0,  0],
                        [ 0, -1,  1, -1, -1, -1, -1],
                        [ 0,  1, -1,  2,  0, -2, -2],
                        [ 0,  1,  0,  0,  1, -1, -3],
                        [ 0, -1,  0, -1,  1,  0,  0],
                        [ 0, -1, -2, -1,  0,  0,  1],
                        [ 0, -1, -2, -3, -2,  1, -1]])
        npt.assert_array_equal(obs, exp)

        # Disable free ends. This is the true global alignment (i.e, the Needleman-
        # Wunsch algorithm in textbooks).
        res = pair_align(seq1, seq2, mode="global", free_ends=False)
        obs = res.paths[0].to_aligned((seq1, seq2))
        exp = ["GAATTC",
               "AGATCT"]
        self.assertListEqual(obs, exp)
        self.assertEqual(res.score, -2)

        # Local alignment (i.e., the Smith-Waterman algorithm in textbooks). It gets
        # subsequences that are best aligned.
        res = pair_align(seq1, seq2, mode="local")
        obs = res.paths[0].to_aligned((seq1, seq2))
        exp = ["GA", "GA"]
        self.assertListEqual(obs, exp)
        self.assertEqual(res.score, 2)

        # Return all alignments.
        res = pair_align(seq1, seq2, mode="local", max_paths=None)
        self.assertEqual(len(res.paths), 3)
        obs = res.paths[1].to_aligned((seq1, seq2))
        self.assertListEqual(obs, ["AT", "AT"])
        obs = res.paths[2].to_aligned((seq1, seq2))
        self.assertListEqual(obs, ["TC", "TC"])

        # Tweak the parameters: match score = 1, mismatch score = -2, linear gap
        # penalty = 2.5. These are the default MegaBLAST parameters.
        res = pair_align(seq1, seq2, sub_score=(1, -2), gap_cost=2.5)

        # With greater mismatch and gap penalties, the two sequences are not aligned
        # at all. (Semi-)global alignment returns the unaligned "alignment" anyway.
        obs = res.paths[0].to_aligned((seq1, seq2))
        exp = ["------GAATTC",
               "AGATCT------"]
        self.assertListEqual(obs, exp)
        self.assertEqual(res.score, 0)

        # Try a different scoring scheme with affine gap penalty: match score = 5,
        # mismatch score = -4, gap opening penalty = 5, gap extension penalty = 2.
        # This triggers the Gotoh algorithm.
        res = pair_align(seq1, seq2, sub_score=(5, -4), gap_cost=(5, 2))
        obs = res.paths[0].to_aligned((seq1, seq2))
        exp = ["-GAATTC-",
               "AGA--TCT"]
        self.assertListEqual(obs, exp)
        self.assertEqual(res.score, 11)

        # The scores match = 5 and mismatch = -4 correspond to the standard
        # nucleotide substitution matrix "NUC.4.4" (a.k.a., DNAfull). We can
        # specify it and get the same result.
        res = pair_align(seq1, seq2, sub_score="NUC.4.4", gap_cost=(5, 2))
        obs = res.paths[0].to_aligned((seq1, seq2))
        self.assertListEqual(obs, exp)

        # Or load an instance of SubstitutionMatrix.
        submat = SubstitutionMatrix.by_name("NUC.4.4")
        res = pair_align(seq1, seq2, sub_score=submat, gap_cost=(5, 2))
        obs = res.paths[0].to_aligned((seq1, seq2))
        self.assertListEqual(obs, exp)

    def test_pair_align_demo(self):
        """Demonstration of very simple cases."""
        # perfect match
        obs = pair_align("ACG", "ACG")
        self.assertEqual(obs.score, 3)
        self.assertEqual(obs.paths[0].to_cigar(), "3M")

        # one mismatch
        obs = pair_align("ACG", "ATG")
        self.assertEqual(obs.score, 1)
        self.assertEqual(obs.paths[0].to_cigar(), "3M")

        # leading insertion
        obs = pair_align("CGT", "ACGT")
        self.assertEqual(obs.score, 3)
        self.assertEqual(obs.paths[0].to_cigar(), "1I3M")

        # trailing deletion
        obs = pair_align("ACGT", "ACG")
        self.assertEqual(obs.score, 3)
        self.assertEqual(obs.paths[0].to_cigar(), "3M1D")

        # overlap
        obs = pair_align("CGT", "ACGTA")
        self.assertEqual(obs.score, 3)
        self.assertEqual(obs.paths[0].to_cigar(), "1I3M1I")

        # one gap
        obs = pair_align("CGT", "ACGTA")
        self.assertEqual(obs.score, 3)
        self.assertEqual(obs.paths[0].to_cigar(), "1I3M1I")

    def test_pair_align_prot(self):
        """Align protein sequences."""
        # Nuclear localization signals (NLS') of SV40 (PKKKRKV) (Kalderon et al.,
        # Cell, 1984) and c-Myc (PAAKRVKLD) (Dang & Lee, Mol Cell Biol, 1988).
        seq1 = Protein("PKKKRKV")
        seq2 = Protein("PAAKRVKLD")

        # default BLASTP parameters
        obs = pair_align(seq1, seq2, sub_score="BLOSUM62", gap_cost=(11, 1))
        self.assertEqual(obs.score, 11)
        self.assertEqual(obs.paths[0].to_cigar(), "7M2I")
        exp = ["PKKKRKV--",
               "PAAKRVKLD"]
        self.assertListEqual(obs.paths[0].to_aligned((seq1, seq2)), exp)

        # Defensins (a category of antimicrobial peptides) of human (2PM4) and rat
        # (NP_775421.1).
        seq1 = Protein("DCYCRIPACIAGEKKYGTCIYQGKLWAFCC")
        seq2 = Protein("CYCRIGACVSGERLTGACGLNGRIYRLCC")
        obs = pair_align(seq1, seq2, sub_score="BLOSUM62", gap_cost=(11, 1))
        self.assertEqual(obs.score, 97)
        self.assertEqual(obs.paths[0].to_cigar(), "1D29M")

        # alternative substitution matrix
        obs = pair_align(seq1, seq2, sub_score="PAM70", gap_cost=(11, 1))
        self.assertEqual(obs.score, 86)
        self.assertEqual(obs.paths[0].to_cigar(), "1D29M")

    def test_pair_align_custom(self):
        """Align custom sequence types."""
        # This example is taken from the previous scikit-bio code.

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

        seq1 = CustomSequence("WXYZ")
        seq2 = CustomSequence("WXYYZZ")
        obs = pair_align(seq1, seq2)
        self.assertEqual(obs.score, 2)
        self.assertListEqual(obs.paths[0].to_aligned((seq1, seq2), gap_char="^"), [
            "WXYZ^^",
            "WXYYZZ",
        ])
        obs = pair_align(seq1, seq2, free_ends=False)
        self.assertEqual(obs.score, 0)
        self.assertListEqual(obs.paths[0].to_aligned((seq1, seq2), gap_char="^"), [
            "WXY^Z^",
            "WXYYZZ",
        ])
        obs = pair_align(seq1, seq2, mode="local")
        self.assertEqual(obs.score, 3)
        self.assertListEqual(obs.paths[0].to_aligned((seq1, seq2), gap_char="^"), [
            "WXY",
            "WXY",
        ])

        # also from the previous scikit-bio code
        seq1 = CustomSequence("YWXXZZYWXXWYYZWXX")
        seq2 = CustomSequence("YWWXZZZYWXYZWWX")
        obs = pair_align(seq1, seq2, mode="local", sub_score=(5, -4),
                         gap_cost=(4.5, 0.5), max_paths=None)
        self.assertEqual(obs.score, 41)
        self.assertEqual(len(obs.paths), 5)
        self.assertEqual(
            TabularMSA.from_path_seqs(obs.paths[4], [seq1, seq2]),
            TabularMSA([CustomSequence("WXXZZYWXXWYYZWXX"),
                        CustomSequence("WXZZZYWX^^^YZWWX")]))

    def test_pair_align_no_grammar(self):
        """Align Sequence objects that are not GrammaredSequence."""
        seq1 = Sequence("GAATTC")
        seq2 = Sequence("AGATCT")
        obs = pair_align(seq1, seq2)
        self.assertEqual(obs.score, 1)
        self.assertListEqual(obs.paths[0].to_aligned((seq1, seq2)), [
            "-GAATTC",
            "AGATCT-",
        ])

    def test_pair_align_strings(self):
        """Align raw strings."""
        seq1 = "accttgcaaaa"
        seq2 = "cctgca"
        obs = pair_align(seq1, seq2, mode="local")
        self.assertEqual(obs.score, 4)
        self.assertEqual(obs.paths[0].to_cigar(), "3M1D3M")
        self.assertListEqual(obs.paths[0].to_aligned((seq1, seq2)), [
            "ccttgca",
            "cct-gca",
        ])

    def test_pair_align_unicode(self):
        """Align Unicode characters."""
        seq1 = "äëïïöëüëö"
        seq2 = "äïïööüëëö"
        obs = pair_align(seq1, seq2)
        self.assertEqual(obs.score, 2)
        self.assertListEqual(obs.paths[0].to_aligned((seq1, seq2)), [
            "äëïïöëüëö-",
            "ä-ïïööüëëö",
        ])

    def test_pair_align_words(self):
        """Align sentences of words."""
        seq1 = "lorem ipsum sit amet tempor".split()
        seq2 = "ipsum sit dolor sed eiusmod".split()
        obs = pair_align(seq1, seq2, mode="local")
        self.assertEqual(obs.score, 2)
        self.assertEqual(obs.paths[0].to_cigar(), "2M")
        npt.assert_array_equal(obs.paths[0].starts, [1, 0])

    def test_pair_align_numbers(self):
        """Align arrays of numbers."""
        seq1 = np.array([2, 13, 52, 11, 27, 8, 33, 77, 25])
        seq2 = np.array([95, 11, 27, 8, 62, 33, 77, 25, 5])
        obs = pair_align(seq1, seq2, mode="local")
        self.assertEqual(obs.score, 4)
        self.assertEqual(obs.paths[0].to_cigar(), "3M1I3M")
        npt.assert_array_equal(obs.paths[0].starts, [3, 1])

    def test_pair_align_degen(self):
        """Align sequences with degenerate codes."""
        # This works because "N" has been defined as a wildcard by the DNA class, and
        # it is also included in the substitution matrix.
        submat = SubstitutionMatrix.identity("ACGTN", 1, -1)
        seq1 = DNA("CGGRATCCA")
        seq2 = DNA("GGAATYCT")
        obs = pair_align(seq1, seq2, sub_score=submat)
        self.assertEqual(obs.score, 2)
        self.assertListEqual(obs.paths[0].to_aligned((seq1, seq2)), [
            "CGGRATCCA",
            "-GGAATYCT",
        ])

        # Won't work if the substitution matrix doesn't contain "N".
        submat = SubstitutionMatrix.identity("ACGT", 1, -1)
        msg = "Sequence 1 contain character(s) absent from the substitution matrix."
        with self.assertRaises(ValueError) as cm:
            _ = pair_align(seq1, seq2, sub_score=submat)
        self.assertEqual(str(cm.exception), msg)

    def test_pair_align_trick_1(self):
        """Tricky case 1."""
        # This example is from Altschul & Erickson, Bull Math Biol, 1986 (Fig. 2),
        # which showed that the original Gotoh algorithm (1982) could fail in some
        # situations. Specifically, if the traceback process does not jump between
        # matrices, it won't be able to distinguish ties between creating a new gap
        # vs. extending an existing gap. In this example, such a tie happens in cell
        # (2, 5). This test could be useful during future attempts to improve the
        # algorithm (e.g., by introducing a traceback matrix).
        seq1 = DNA("AGT")
        seq2 = DNA("TGAGTT")
        obs = pair_align(seq1, seq2, sub_score=(0, -1), gap_cost=(1, 1),
                         free_ends=False, max_paths=None, keep_matrices=True)
        self.assertEqual(obs.score, -5)
        self.assertEqual(len(obs.paths), 3)
        self.assertEqual(obs.paths[0].to_cigar(), "2I3M1I")
        self.assertEqual(obs.paths[1].to_cigar(), "2M3I1M")
        self.assertEqual(obs.paths[2].to_cigar(), "2I2M1I1M")
        npt.assert_array_equal(obs.matrices[0], np.array([
            [ 0, -2, -3, -4, -5, -6, -7],
            [-2, -1, -3, -3, -5, -6, -7],
            [-3, -3, -1, -3, -3, -5, -6],
            [-4, -3, -3, -2, -4, -3, -5]]))
        npt.assert_array_equal(obs.matrices[1][1:, 1:], np.array([
            [-4, -3, -4, -5, -6, -7],
            [-5, -5, -3, -4, -5, -6],
            [-6, -5, -5, -4, -5, -5]]))
        npt.assert_array_equal(obs.matrices[2][1:, 1:], np.array([
            [-4, -5, -6, -7, -8, -9],
            [-3, -5, -5, -7, -8, -9],
            [-4, -3, -5, -5, -7, -8]]))

    def test_pair_align_trick_2(self):
        """Tricky case 2."""
        # This example is from Flouri et al., BioRxiv, 2015, which stated that the
        # original Gotoh algorithm (1982) had a mistake in the initialization of the
        # matrices. Specifically, the insertion and deletion matrices should be
        # initialized with at least 2 gap_open + k gap_extend. If one does this with
        # gap_open + k * gap_extend, as the Gotoh paper suggested, this example will
        # break.
        obs = pair_align("A", "G", sub_score=(0, -5), gap_cost=(2, 1), free_ends=False)
        self.assertEqual(obs.score, -5)
        self.assertEqual(obs.paths[0].to_cigar(), "1M")

    def test_pair_align_semi(self):
        """Customization of end gap policy."""
        # case 1: search longer seq1 for shorter seq1
        seq1 = "TCGA"
        seq2 = "ATCGAG"

        # true global: end gaps penalized
        obs = pair_align(seq1, seq2, free_ends=False)
        self.assertEqual(obs.score, 0)
        self.assertEqual(obs.paths[0].to_cigar(), "1I4M1I")

        # ends completely free
        obs = pair_align(seq1, seq2, free_ends=True)
        self.assertEqual(obs.score, 4)
        self.assertEqual(obs.paths[0].to_cigar(), "1I4M1I")

        # seq1 free (typical semi-global)
        obs = pair_align(seq1, seq2, free_ends=(True, False))
        self.assertEqual(obs.score, 4)

        # seq2 free (not helpful here)
        obs = pair_align(seq1, seq2, free_ends=(False, True))
        self.assertEqual(obs.score, 0)

        # case 2: seq1 can align to either side of seq2
        seq1 = "ACGT"
        seq2 = "ACGTACGT"
        obs = pair_align(seq1, seq2, free_ends=True, max_paths=None)
        self.assertEqual(obs.score, 4)
        self.assertEqual(len(obs.paths), 2)
        self.assertEqual(obs.paths[0].to_cigar(), "4M4I")
        self.assertEqual(obs.paths[1].to_cigar(), "4I4M")

        # force align to left
        obs = pair_align(seq1, seq2, free_ends=(False, True, False, False), max_paths=None)
        self.assertEqual(obs.score, 4)
        path, = obs.paths
        self.assertEqual(path.to_cigar(), "4M4I")

        # force align to right
        obs = pair_align(seq1, seq2, free_ends=(True, False, False, False), max_paths=None)
        self.assertEqual(obs.score, 4)
        path, = obs.paths
        self.assertEqual(path.to_cigar(), "4I4M")

        # case 3: seqs overlap on either side
        seq1 = "GAATTC"
        seq2 = "TTCGAA"
        obs = pair_align(seq1, seq2, free_ends=True, max_paths=None)
        path1, path2 = obs.paths
        self.assertEqual(path1.to_cigar(), "3I3M3D")
        self.assertEqual(path2.to_cigar(), "3D3M3I")

        # seq1 left connects seq2 right
        obs = pair_align(seq1, seq2, free_ends=(True, False, False, True), max_paths=None)
        path, = obs.paths
        self.assertEqual(path.to_cigar(), "3I3M3D")

        # seq1 right connects seq2 left
        obs = pair_align(seq1, seq2, free_ends=(False, True, True, False), max_paths=None)
        path, = obs.paths
        self.assertEqual(path.to_cigar(), "3D3M3I")

    def test_pair_align_text(self):
        """Text analysis."""
        # edit (Levenshtein) distance
        txts = ["here", "hear"]
        obs = pair_align(*txts, sub_score=(0, -1), gap_cost=1, free_ends=False)
        self.assertEqual(obs.score, -2)
        self.assertListEqual(obs.paths[0].to_aligned(txts), ["he-re", "hear-"])

        # longest common subsequence (LCS)
        txts = ["strawberry", "stranger"]
        obs = pair_align(*txts, mode="local", sub_score=(1, -np.inf), gap_cost=0)
        exp = ["stra--wber", "strang--er"]  # LCS is "straer"; gap allowed
        self.assertListEqual(obs.paths[0].to_aligned(txts), exp)

        # longest common substring
        obs = pair_align(*txts, mode="local", sub_score=(1, -np.inf), gap_cost=np.inf)
        exp = ["stra", "stra"]  # "stra"; gap not allowed
        self.assertListEqual(obs.paths[0].to_aligned(txts), exp)

    def test_pair_align_trim(self):
        """Trim free gaps."""
        # make match score large such that "AAA" in both sequences alway align
        scoring = dict(sub_score=(9, -1))
        # baseline (-AAA- vs TAAAT)
        obs = pair_align("AAA", "TAAAT", **scoring)
        self.assertEqual(obs.paths[0].to_cigar(), "1I3M1I")
        # trim both ends
        obs = pair_align("AAA", "TAAAT", **scoring, trim_ends=True)
        self.assertEqual(obs.paths[0].to_cigar(), "3M")
        # no free, no trim
        obs = pair_align("AAA", "TAAAT", **scoring, free_ends=False, trim_ends=True)
        self.assertEqual(obs.paths[0].to_cigar(), "1I3M1I")
        # trim seq1 (can trim)
        obs = pair_align("AAA", "TAAAT", **scoring, free_ends=(True, False),
                         trim_ends=True)
        self.assertEqual(obs.paths[0].to_cigar(), "3M")
        # trim seq2 (nothing to trim)
        obs = pair_align("AAA", "TAAAT", **scoring, free_ends=(False, True),
                         trim_ends=True)
        self.assertEqual(obs.paths[0].to_cigar(), "1I3M1I")
        # trim leading gap
        obs = pair_align("AAA", "TAAAT", **scoring, free_ends=(True, False, False, False),
                         trim_ends=True)
        self.assertEqual(obs.paths[0].to_cigar(), "3M1I")
        # no impact
        obs = pair_align("AAA", "TAAAT", **scoring, free_ends=(True, False, True, False),
                         trim_ends=True)
        self.assertEqual(obs.paths[0].to_cigar(), "3M1I")
        # trim trailing gap
        obs = pair_align("AAA", "TAAAT", **scoring, free_ends=(False, True, False, False),
                         trim_ends=True)
        self.assertEqual(obs.paths[0].to_cigar(), "1I3M")
        # no impact
        obs = pair_align("AAA", "TAAAT", **scoring, free_ends=(False, True, False, True),
                         trim_ends=True)
        self.assertEqual(obs.paths[0].to_cigar(), "1I3M")

    def test_pair_align_search(self):
        """Search for short query in long target."""
        obs = pair_align("r", "strawberry", trim_ends=True, max_paths=None).paths
        self.assertEqual(len(obs), 3)
        self.assertEqual(obs[0].starts[1], 2)
        self.assertEqual(obs[1].starts[1], 7)
        self.assertEqual(obs[2].starts[1], 8)

    def test_pair_align_64bit(self):
        """Float64 scoring scheme."""
        obs = pair_align("GAATTC", "AGATCT", keep_matrices=True)
        self.assertEqual(obs.paths[0].to_cigar(), "1I5M1D")
        self.assertEqual(obs.matrices[0].dtype, np.float32)

        submat = SubstitutionMatrix.identity("ACGTN", 1, -1, dtype="float64")
        obs = pair_align("GAATTC", "AGATCT", sub_score=submat, keep_matrices=True)
        self.assertEqual(obs.paths[0].to_cigar(), "1I5M1D")
        self.assertEqual(obs.matrices[0].dtype, np.float64)

    def test_pair_align_edge(self):
        """Edge cases."""
        # distinct sequences: global alignment completely misaligned
        obs = pair_align("AAA", "TTT", sub_score=(1, -1), gap_cost=2)
        self.assertEqual(obs.paths[0].to_cigar(), "3I3D")

        # local alignment returns no path
        obs = pair_align("AAA", "TTT", mode="local", sub_score=(1, -1), gap_cost=2)
        self.assertEqual(len(obs.paths), 0)

        # single character, identical
        obs = pair_align("A", "A")
        self.assertEqual(obs.paths[0].to_cigar(), "1M")
        obs = pair_align("A", "A", free_ends=False)
        self.assertEqual(obs.paths[0].to_cigar(), "1M")
        obs = pair_align("A", "A", mode="local")
        self.assertEqual(obs.paths[0].to_cigar(), "1M")

        # single character, different
        obs = pair_align("A", "T")
        self.assertEqual(obs.paths[0].to_cigar(), "1I1D")
        obs = pair_align("A", "T", free_ends=False)
        self.assertEqual(obs.paths[0].to_cigar(), "1M")
        obs = pair_align("A", "T", mode="local")
        self.assertEqual(len(obs.paths), 0)  # no path returned

        # gap cost is negative (more favorable than match)
        obs = pair_align("A", "A", sub_score=(1, -1), gap_cost=-2).paths[0]
        self.assertEqual(obs.to_cigar(), "1I1D")
        self.assertListEqual(obs.to_aligned(("A", "A")), ["-A", "A-"])

        # mismatch score greater than match score
        obs = pair_align("ATG", "ATG", sub_score=(1, 2)).paths[0]
        self.assertEqual(obs.to_cigar(), "1I2M1D")
        self.assertListEqual(obs.to_aligned(("ATG", "ATG")), ["-ATG", "ATG-"])

        # Make sure path length is always <= sum of sequence lengths. The maximum
        # happens when all characters are misaligned. This check is to make sure
        # the pre-allocation of path in the trackback functions does not lead to
        # overflow.
        seq1, seq2 = "AAAAA", "TTTTT"
        obs = pair_align(seq1, seq2).paths[0]
        self.assertEqual(obs.to_cigar(), "5I5D")
        self.assertLessEqual(obs.lengths.sum(), len(seq1) + len(seq2))

        # complete misalignment becomes empty after trimming ends
        obs = pair_align(seq1, seq2, trim_ends=True).paths[0]
        self.assertEqual(obs.to_cigar(), "")
        self.assertListEqual(obs.ranges.tolist(), [[0, 0], [5, 5]])

    def test_pair_align_inf(self):
        """Infinite scores."""
        # infinite match score (returned score is also infinite)
        obs = pair_align("ATG", "ATG", sub_score=(np.inf, 0))
        self.assertTrue(np.isposinf(obs.score))

        # infinite match score not working (there is no match)
        obs = pair_align("AT", "GC", sub_score=(np.inf, 0))
        self.assertEqual(obs.score, 0)

        # infinite gap cost (which disables gaps)
        obs = pair_align("ATT", "GGC", sub_score=(1, -100), gap_cost=np.inf,
                         free_ends=False)
        self.assertEqual(obs.score, -300)
        self.assertEqual(obs.paths[0].to_cigar(), "3M")

        # make sure things don't break
        obs = pair_align("AT", "A", sub_score=(np.inf, np.inf),
                         gap_cost=(np.inf, np.inf))
        self.assertTrue(np.isposinf(obs.score))

    def test_pair_align_nan(self):
        """NaN scores."""
        # this is nonsense; just want to make sure things don't break
        obs = pair_align("AT", "A", sub_score=(np.nan, np.nan),
                         gap_cost=(np.nan, np.nan))
        self.assertTrue(np.isnan(obs.score))

    def test_pair_align_error(self):
        """Errors."""
        with self.assertRaises(ValueError) as cm:
            _ = pair_align("", "ATG")
        self.assertEqual(str(cm.exception), "Sequence 1 has a length of zero.")
        with self.assertRaises(ValueError) as cm:
            _ = pair_align("ATG", "")
        self.assertEqual(str(cm.exception), "Sequence 2 has a length of zero.")
        with self.assertRaises(ValueError) as cm:
            _ = pair_align("ACGT", "TCGA", mode="xyz")
        self.assertEqual(str(cm.exception), "Invalid mode: xyz.")

    def test_pair_align_real_1(self):
        # 16S rRNA gene sequences of Escherichia coli and Streptococcus pneumoniae
        with open(get_data_path("16s.frn"), "r") as fh:
            seq1, seq2 = skbio.io.read(fh, format="fasta", constructor=DNA)

        # Perform semi-global alignment using BLASTN's default parameters.
        obs = pair_align(seq1, seq2, sub_score=(2, -3), gap_cost=(5, 2))
        self.assertEqual(obs.score, 1199)

        # Perform global alignment using EMBOSS needle's default parameters.
        # Note: the score should be exact because it only involves half-integers.
        # `assertAlmostEqual` is used here just in case.
        obs = pair_align(seq1, seq2, sub_score="NUC.4.4", gap_cost=(9.5, 0.5),
                         free_ends=False)
        self.assertAlmostEqual(obs.score, 4355.5)

        # Map the 16S rRNA V4 primers (515F and 806R) to the E. coli sequence. Note
        # that there are degenerate characters in the primer sequences.
        # Source: https://earthmicrobiome.org/protocols-and-standards/16s/
        seq515F = DNA("GTGYCAGCMGCCGCGGTAA")
        seq806R = DNA("GGACTACNVGGGTWTCTAAT").reverse_complement()

        # Perform semi-global alignment using BLASTN's default parameters.
        # There is a unique alignment.
        obs = pair_align(seq515F, seq1, sub_score=(2, -3), gap_cost=(5, 2), max_paths=10)
        self.assertEqual(obs.score, 28)
        self.assertEqual(len(obs.paths), 1)
        self.assertEqual(obs.paths[0].to_cigar(), "505I19M926I")
        msa = TabularMSA.from_path_seqs(obs.paths[0], [seq515F, seq1])[:, 500:530]
        self.assertEqual(str(msa[0]), "-----GTGYCAGCMGCCGCGGTAA------")
        self.assertEqual(str(msa[1]), "ACTCCGTGCCAGCAGCCGCGGTAATACGGA")

        obs = pair_align(seq806R, seq1, sub_score=(2, -3), gap_cost=(5, 2), max_paths=10)
        self.assertEqual(obs.score, 25)
        self.assertEqual(len(obs.paths), 1)
        self.assertEqual(obs.paths[0].to_cigar(), "777I20M653I")
        msa = TabularMSA.from_path_seqs(obs.paths[0], [seq806R, seq1])[:, 770:805]
        self.assertEqual(str(msa[0]), "-------ATTAGAWACCCBNGTAGTCC--------")
        self.assertEqual(str(msa[1]), "AAACAGGATTAGATACCCTGGTAGTCCACGCCGTA")

        # Perform local alignment instead.
        obs = pair_align(seq515F, seq1, mode="local", sub_score=(2, -3),
                         gap_cost=(5, 2), max_paths=10)
        self.assertEqual(obs.score, 28)
        self.assertEqual(len(obs.paths), 1)
        self.assertEqual(obs.paths[0].to_cigar(), "19M")

    def test_pair_align_real_2(self):
        # Insulin sequences of human and mouse
        with open(get_data_path("insulin.faa"), "r") as fh:
            seq1, seq2 = skbio.io.read(fh, format="fasta", constructor=Protein)

        # Perform semi-global alignment using BLASTP's default parameters.
        obs = pair_align(seq1, seq2, sub_score="BLOSUM62", gap_cost=(11, 1))
        self.assertEqual(obs.score, 440)

    def test_pair_align_nucl_wrap(self):
        """Nucleotide sequence alignment."""
        seq1 = DNA("GATCGTC")
        seq2 = DNA("ATCGCTC")
        obs = pair_align_nucl(seq1, seq2)
        self.assertEqual(obs.score, 5)
        self.assertEqual(obs.paths[0].to_cigar(), "1D4M1I2M")
        self.assertListEqual(obs.paths[0].to_aligned((seq1, seq2)), [
            "GATCG-TC",
            "-ATCGCTC"])

        kwargs = dict(sub_score="NUC.4.4", gap_cost=7)
        obs = pair_align_nucl(seq1, seq2, **kwargs)
        self.assertEqual(obs.score, 23)

    def test_pair_align_prot_wrap(self):
        """Protein sequence alignment."""
        seq1 = Protein("PKKKRKV")
        seq2 = Protein("PAAKRVKLD")
        obs = pair_align_prot(seq1, seq2)
        self.assertEqual(obs.score, 11)
        self.assertEqual(obs.paths[0].to_cigar(), "7M2I")
        self.assertListEqual(obs.paths[0].to_aligned((seq1, seq2)), [
            "PKKKRKV--",
            "PAAKRVKLD"])

        kwargs = dict(sub_score="PAM70", gap_cost=6)
        obs = pair_align_prot(seq1, seq2, **kwargs)
        self.assertEqual(obs.score, 13)

    def test_prep_free_ends(self):
        obs = _prep_free_ends(local=True, free_ends=True, trim_ends=False)
        self.assertTupleEqual(obs, (2, 2, 2, 2))
        obs = _prep_free_ends(local=False, free_ends=False, trim_ends=False)
        self.assertTupleEqual(obs, (0, 0, 0, 0))
        obs = _prep_free_ends(local=False, free_ends=True, trim_ends=False)
        self.assertTupleEqual(obs, (1, 1, 1, 1))
        obs = _prep_free_ends(local=False, free_ends=True, trim_ends=True)
        self.assertTupleEqual(obs, (2, 2, 2, 2))
        obs = _prep_free_ends(local=False, free_ends=(True, False), trim_ends=False)
        self.assertTupleEqual(obs, (1, 1, 0, 0))
        obs = _prep_free_ends(local=False, free_ends=(True, False), trim_ends=True)
        self.assertTupleEqual(obs, (2, 2, 0, 0))
        obs = _prep_free_ends(
            local=False, free_ends=(True, False, False, True), trim_ends=True)
        self.assertTupleEqual(obs, (2, 0, 0, 2))
        obs = _prep_free_ends(
            local=False, free_ends=(True, False, True, False), trim_ends=False)
        self.assertTupleEqual(obs, (1, 0, 1, 0))

        # non-Boolean inputs
        obs = _prep_free_ends(local=False, free_ends=(0, 1, 2, 3), trim_ends="yes")
        self.assertTupleEqual(obs, (0, 2, 2, 2))

        msg = "`free_ends` must be one, two or four Booleans."
        with self.assertRaises(ValueError) as cm:
            _ = _prep_free_ends(local=False, free_ends=(1, 2, 3), trim_ends=False)
        self.assertEqual(str(cm.exception), msg)

    def test_alloc_matrices(self):
        # linear
        obs = _alloc_matrices(3, 4, affine=False)
        self.assertEqual(len(obs), 1)
        self.assertTrue(isinstance(obs[0], np.ndarray))
        self.assertTupleEqual(obs[0].shape, (4, 5))
        self.assertEqual(obs[0].dtype, np.float32)

        # ensure matrix is C-contiguous
        self.assertTrue(obs[0].flags.c_contiguous)

        # affine
        obs = _alloc_matrices(3, 4, affine=True, dtype=np.float32)
        self.assertEqual(len(obs), 3)
        for mat in obs:
            self.assertTupleEqual(mat.shape, (4, 5))
            self.assertEqual(mat.dtype, np.float32)

        # float64
        obs = _alloc_matrices(3, 4, affine=False, dtype=np.float64)
        self.assertEqual(obs[0].dtype, np.float64)

    def test_init_matrices(self):
        # classic global alignment
        obs = _alloc_matrices(3, 4, affine=False)
        _init_matrices(obs, 0, 2, local=False, lead1=False, lead2=False)
        npt.assert_array_equal(obs[0][:, 0], [0, -2, -4, -6])
        npt.assert_array_equal(obs[0][0, :], [0, -2, -4, -6, -8])

        # local alignment
        obs = _alloc_matrices(3, 4, affine=False)
        _init_matrices(obs, 0, 2, local=True, lead1=False, lead2=False)
        self.assertTrue((obs[0][:, 0] == 0).all())
        self.assertTrue((obs[0][0, :] == 0).all())

        # semi-global alignment (free ends)
        obs = _alloc_matrices(3, 4, affine=False)
        _init_matrices(obs, 0, 2, local=False, lead1=True, lead2=True)
        self.assertTrue((obs[0][:, 0] == 0).all())
        self.assertTrue((obs[0][0, :] == 0).all())

        # free end for seq1
        obs = _alloc_matrices(3, 4, affine=False)
        _init_matrices(obs, 0, 2, local=False, lead1=True, lead2=False)
        npt.assert_array_equal(obs[0][:, 0], [0, -2, -4, -6])
        self.assertTrue((obs[0][0, :] == 0).all())

        # free end for seq2
        obs = _alloc_matrices(3, 4, affine=False)
        _init_matrices(obs, 0, 2, local=False, lead1=False, lead2=True)
        self.assertTrue((obs[0][:, 0] == 0).all())
        npt.assert_array_equal(obs[0][0, :], [0, -2, -4, -6, -8])

        # affine and global
        obs = _alloc_matrices(3, 4, affine=True)
        _init_matrices(obs, 5, 2, local=False, lead1=False, lead2=False)
        npt.assert_array_equal(obs[0][:, 0], [0, -7, -9, -11])
        npt.assert_array_equal(obs[0][0, :], [0, -7, -9, -11, -13])
        self.assertTrue(np.isneginf(obs[1][1:, 0]).all())
        self.assertTrue(np.isneginf(obs[2][0, 1:]).all())

        # affine and local
        obs = _alloc_matrices(3, 4, affine=True)
        _init_matrices(obs, 5, 2, local=True, lead1=False, lead2=False)
        self.assertTrue((obs[0][:, 0] == 0).all())
        self.assertTrue((obs[0][0, :] == 0).all())
        self.assertTrue(np.isneginf(obs[1][1:, 0]).all())
        self.assertTrue(np.isneginf(obs[2][0, 1:]).all())

    def test_fill_linear_matrix(self):
        seq1 = DNA("ACGT")
        seq2 = DNA("AACTG")
        (seq1, seq2), submat, _ = encode_sequences((seq1, seq2), (1, -1))
        query, target = np.ascontiguousarray(submat[seq1]), seq2
        dtype=query.dtype.type
        m, n = seq1.size, seq2.size

        # global
        obs = _alloc_matrices(m, n, False, dtype=dtype)
        _init_matrices(obs, 0, 2, local=False, lead1=False, lead2=False)
        _fill_linear_matrix(obs[0], query, target, 2, local=False)
        exp = np.array([[  0,  -2,  -4,  -6,  -8, -10],
                        [ -2,   1,  -1,  -3,  -5,  -7],
                        [ -4,  -1,   0,   0,  -2,  -4],
                        [ -6,  -3,  -2,  -1,  -1,  -1],
                        [ -8,  -5,  -4,  -3,   0,  -2]])
        npt.assert_array_equal(obs[0], exp)

        # local
        obs = _alloc_matrices(m, n, False, dtype=dtype)
        _init_matrices(obs, 0, 2, local=True, lead1=False, lead2=False)
        _fill_linear_matrix(obs[0], query, target, 2, local=True)
        exp = np.array([[0, 0, 0, 0, 0, 0],
                        [0, 1, 1, 0, 0, 0],
                        [0, 0, 0, 2, 0, 0],
                        [0, 0, 0, 0, 1, 1],
                        [0, 0, 0, 0, 1, 0]])
        npt.assert_array_equal(obs[0], exp)

        # semi-global
        obs = _alloc_matrices(m, n, False, dtype=dtype)
        _init_matrices(obs, 0, 2, local=False, lead1=True, lead2=True)
        _fill_linear_matrix(obs[0], query, target, 2, local=False)
        exp = np.array([[ 0,  0,  0,  0,  0,  0],
                        [ 0,  1,  1, -1, -1, -1],
                        [ 0, -1,  0,  2,  0, -2],
                        [ 0, -1, -2,  0,  1,  1],
                        [ 0, -1, -2, -2,  1,  0]])
        npt.assert_array_equal(obs[0], exp)

    def test_fill_affine_matrices(self):
        seq1 = DNA("ACGT")
        seq2 = DNA("AACTG")
        (seq1, seq2), submat, _ = encode_sequences((seq1, seq2), (1, -1))
        query, target = np.ascontiguousarray(submat[seq1]), seq2
        dtype=query.dtype.type
        m, n = seq1.size, seq2.size

        # global
        obs = _alloc_matrices(m, n, True, dtype=dtype)
        _init_matrices(obs, 3, 1, local=False, lead1=False, lead2=False)
        _fill_affine_matrices(*obs, query, target, 3, 1, local=False)
        exp0 = np.array([[ 0, -4, -5, -6, -7, -8],
                         [-4,  1, -3, -4, -5, -6],
                         [-5, -3,  0, -2, -5, -6],
                         [-6, -4, -4, -1, -3, -4],
                         [-7, -5, -5, -5,  0, -4]])
        npt.assert_array_equal(obs[0], exp0)
        exp1 = np.array([[ -8,  -3,  -4,  -5,  -6],
                         [ -9,  -7,  -4,  -5,  -6],
                         [-10,  -8,  -8,  -5,  -6],
                         [-11,  -9,  -9,  -9,  -4]])
        npt.assert_array_equal(obs[1][1:, 1:], exp1)
        exp2 = np.array([[ -8,  -9, -10, -11, -12],
                         [ -3,  -7,  -8,  -9, -10],
                         [ -4,  -4,  -6,  -9, -10],
                         [ -5,  -5,  -5,  -7,  -8]])
        npt.assert_array_equal(obs[2][1:, 1:], exp2)

        # local
        obs = _alloc_matrices(m, n, True, dtype=dtype)
        _init_matrices(obs, 3, 1, local=True, lead1=False, lead2=False)
        _fill_affine_matrices(*obs, query, target, 3, 1, local=True)
        exp0 = np.array([[0, 0, 0, 0, 0, 0],
                         [0, 1, 1, 0, 0, 0],
                         [0, 0, 0, 2, 0, 0],
                         [0, 0, 0, 0, 1, 1],
                         [0, 0, 0, 0, 1, 0]])
        npt.assert_array_equal(obs[0], exp0)
        exp1 = np.array([[-4, -3, -3, -4, -4],
                         [-4, -4, -4, -2, -3],
                         [-4, -4, -4, -4, -3],
                         [-4, -4, -4, -4, -3]])
        npt.assert_array_equal(obs[1][1:, 1:], exp1)
        exp2 = np.array([[-4, -4, -4, -4, -4],
                         [-3, -3, -4, -4, -4],
                         [-4, -4, -2, -4, -4],
                         [-4, -4, -3, -3, -3]])
        npt.assert_array_equal(obs[2][1:, 1:], exp2)

        # semi-global
        obs = _alloc_matrices(m, n, True, dtype=dtype)
        _init_matrices(obs, 3, 1, local=False, lead1=True, lead2=True)
        _fill_affine_matrices(*obs, query, target, 3, 1, local=False)
        exp0 = np.array([[ 0,  0,  0,  0,  0,  0],
                         [ 0,  1,  1, -1, -1, -1],
                         [ 0, -1,  0,  2, -2, -2],
                         [ 0, -1, -2, -1,  1, -1],
                         [ 0, -1, -2, -3,  0,  0]])
        npt.assert_array_equal(obs[0], exp0)
        exp1 = np.array([[-4, -3, -3, -4, -5],
                         [-4, -5, -4, -2, -3],
                         [-4, -5, -6, -5, -3],
                         [-4, -5, -6, -7, -4]])
        npt.assert_array_equal(obs[1][1:, 1:], exp1)
        exp2 = np.array([[-4, -4, -4, -4, -4],
                         [-3, -3, -5, -5, -5],
                         [-4, -4, -2, -6, -6],
                         [-5, -5, -3, -3, -5]])
        npt.assert_array_equal(obs[2][1:, 1:], exp2)

    def test_one_stop(self):
        # ACGT vs AACTG, sub=(1, -1), gap=2, same below, unless otherwise stated
        # global (bottom-right corner)
        scomat = np.array([[  0,  -2,  -4,  -6,  -8, -10],
                           [ -2,   1,  -1,  -3,  -5,  -7],
                           [ -4,  -1,   0,   0,  -2,  -4],
                           [ -6,  -3,  -2,  -1,  -1,  -1],
                           [ -8,  -5,  -4,  -3,   0,  -2]])
        obs = _one_stop(scomat, local=False, trail1=False, trail2=False)
        self.assertEqual(obs[0], -2)
        npt.assert_array_equal(obs[1], [[4, 5]])

        # local (maximum in the matrix)
        scomat = np.array([[0, 0, 0, 0, 0, 0],
                           [0, 1, 1, 0, 0, 0],
                           [0, 0, 0, 2, 0, 0],
                           [0, 0, 0, 0, 1, 1],
                           [0, 0, 0, 0, 1, 0]])
        obs = _one_stop(scomat, local=True, trail1=False, trail2=False)
        self.assertEqual(obs[0], 2)
        npt.assert_array_equal(obs[1], [[2, 3]])

        # local, tie (will pick the smaller index)
        scomat = np.array([[0, 0, 0, 0, 0, 0],
                           [0, 0, 1, 0, 0, 0],
                           [0, 0, 0, 2, 0, 1],
                           [0, 1, 0, 0, 1, 0],
                           [0, 0, 2, 0, 0, 0]])
        obs = _one_stop(scomat, local=True, trail1=False, trail2=False)
        self.assertEqual(obs[0], 2)
        npt.assert_array_equal(obs[1], [[2, 3]])

        # semi-global (maximum in bottom row)
        scomat = np.array([[ 0,  0,  0,  0,  0,  0,  0],
                           [ 0, -1,  1, -1, -1, -1, -1],
                           [ 0, -1, -2,  2, -2, -2, -2],
                           [ 0, -1, -2, -2,  3, -1, -1],
                           [ 0,  1, -2, -3, -1,  4,  0]])
        obs = _one_stop(scomat, local=False, trail1=True, trail2=True)
        self.assertEqual(obs[0], 4)
        npt.assert_array_equal(obs[1], [[4, 5]])

        # bottom row only
        obs = _one_stop(scomat, local=False, trail1=True, trail2=False)
        self.assertEqual(obs[0], 4)
        npt.assert_array_equal(obs[1], [[4, 5]])

        # right-most column only
        obs = _one_stop(scomat, local=False, trail1=False, trail2=True)
        self.assertEqual(obs[0], 0)
        npt.assert_array_equal(obs[1], [[0, 6]])

        # transpose and maximum in right-most column
        scomat = np.ascontiguousarray(scomat.T)
        obs = _one_stop(scomat, local=False, trail1=True, trail2=True)
        self.assertEqual(obs[0], 4)
        npt.assert_array_equal(obs[1], [[5, 4]])

        # semi-global, tie
        scomat = np.array([[ 0,  0,  0,  0,  0,  0],
                           [ 0,  1,  1, -1, -1, -1],
                           [ 0, -1,  0,  2,  0, -2],
                           [ 0, -1, -2,  0,  1,  1],
                           [ 0, -1, -2, -2,  1,  0]])
        obs = _one_stop(scomat, local=False, trail1=True, trail2=True)
        self.assertEqual(obs[0], 1)
        npt.assert_array_equal(obs[1], [[3, 5]])

    def test_all_stops(self):
        # ACGT vs AACTG, sub=(1, -1), gap=2, same below, unless otherwise stated
        # global (always unique)
        scomat = np.array([[  0,  -2,  -4,  -6,  -8, -10],
                           [ -2,   1,  -1,  -3,  -5,  -7],
                           [ -4,  -1,   0,   0,  -2,  -4],
                           [ -6,  -3,  -2,  -1,  -1,  -1],
                           [ -8,  -5,  -4,  -3,   0,  -2]])
        obs = _all_stops(scomat, local=False, trail1=False, trail2=False)
        self.assertEqual(obs[0], -2)
        npt.assert_array_equal(obs[1], [[4, 5]])

        # local, unique
        scomat = np.array([[0, 0, 0, 0, 0, 0],
                           [0, 1, 1, 0, 0, 0],
                           [0, 0, 0, 2, 0, 0],
                           [0, 0, 0, 0, 1, 1],
                           [0, 0, 0, 0, 1, 0]])
        obs = _all_stops(scomat, local=True, trail1=False, trail2=False)
        self.assertEqual(obs[0], 2)
        npt.assert_array_equal(obs[1], [[2, 3]])

        # local, tie (CGTC vs. TCGAG)
        scomat = np.array([[0, 0, 0, 0, 0, 0],
                           [0, 0, 1, 0, 0, 0],
                           [0, 0, 0, 2, 0, 1],
                           [0, 1, 0, 0, 1, 0],
                           [0, 0, 2, 0, 0, 0]])
        obs = _all_stops(scomat, local=True, trail1=False, trail2=False)
        self.assertEqual(obs[0], 2)
        npt.assert_array_equal(obs[1], [[2, 3], [4, 2]])

        # exact comparison
        obs = _all_stops(scomat, local=True, trail1=False, trail2=False, eps=None)
        self.assertEqual(obs[0], 2)
        npt.assert_array_equal(obs[1], [[2, 3], [4, 2]])

        # semi-global, unique (TCGA vs. ATCGAG)
        scomat = np.array([[ 0,  0,  0,  0,  0,  0,  0],
                           [ 0, -1,  1, -1, -1, -1, -1],
                           [ 0, -1, -2,  2, -2, -2, -2],
                           [ 0, -1, -2, -2,  3, -1, -1],
                           [ 0,  1, -2, -3, -1,  4,  0]])
        obs = _all_stops(scomat, local=False, trail1=True, trail2=True)
        self.assertEqual(obs[0], 4)
        npt.assert_array_equal(obs[1], [[4, 5]])

        # semi-global, tie
        scomat = np.array([[ 0,  0,  0,  0,  0,  0],
                           [ 0,  1,  1, -1, -1, -1],
                           [ 0, -1,  0,  2,  0, -2],
                           [ 0, -1, -2,  0,  1,  1],
                           [ 0, -1, -2, -2,  1,  0]])
        obs = _all_stops(scomat, local=False, trail1=True, trail2=True)
        self.assertEqual(obs[0], 1)
        npt.assert_array_equal(obs[1], [[3, 5], [4, 4]])

        # bottom row only
        obs = _all_stops(scomat, local=False, trail1=True, trail2=False)
        self.assertEqual(obs[0], 1)
        npt.assert_array_equal(obs[1], [[4, 4]])

        # right-most column only
        obs = _all_stops(scomat, local=False, trail1=False, trail2=True)
        self.assertEqual(obs[0], 1)
        npt.assert_array_equal(obs[1], [[3, 5]])

        # exact comparison
        obs = _all_stops(scomat, local=False, trail1=True, trail2=True, eps=None)
        self.assertEqual(obs[0], 1)
        npt.assert_array_equal(obs[1], [[3, 5], [4, 4]])

        # ensure ascending order
        scomat = np.array([[0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 1],
                           [0, 0, 0, 0, 1],
                           [0, 0, 1, 1, 1]])
        obs = _all_stops(scomat, local=False, trail1=True, trail2=True)
        self.assertEqual(obs[0], 1)
        npt.assert_array_equal(obs[1], [[1, 4], [2, 4], [3, 2], [3, 3], [3, 4]])

    def test_traceback_one(self):
        # ACTCA vs CAGAG, sub=(1, -1), gap=2 (linear)
        # global: unique path, bottom-right corner to left-most column
        scomat = np.array([[  0,  -2,  -4,  -6,  -8, -10],
                           [ -2,  -1,  -1,  -3,  -5,  -7],
                           [ -4,  -1,  -2,  -2,  -4,  -6],
                           [ -6,  -3,  -2,  -3,  -3,  -5],
                           [ -8,  -5,  -4,  -3,  -4,  -4],
                           [-10,  -7,  -4,  -5,  -2,  -4]], dtype=float)
        obs = _traceback_one(5, 5, (scomat,), 0., 2.)
        self.assertEqual(obs.to_cigar(), "1D4M1I")

        # transpose: bottom-right corner to top row
        scomat = np.ascontiguousarray(scomat.T)
        obs = _traceback_one(5, 5, (scomat,), 0., 2.)
        self.assertEqual(obs.to_cigar(), "1I4M1D")

        # semi-global: unique path, bottom row to left-most column
        scomat = np.array([[ 0,  0,  0,  0,  0,  0],
                           [ 0, -1,  1, -1,  1, -1],
                           [ 0,  1, -1,  0, -1,  0],
                           [ 0, -1,  0, -2, -1, -2],
                           [ 0,  1, -1, -1, -3, -2],
                           [ 0, -1,  2,  0,  0, -2]], dtype=float)
        obs = _traceback_one(5, 2, (scomat,), 0., 2.)
        self.assertEqual(obs.to_cigar(), "3D2M3I")

        # transpose: right-most column to top row
        scomat = np.ascontiguousarray(scomat.T)
        obs = _traceback_one(2, 5, (scomat,), 0., 2.)
        self.assertEqual(obs.to_cigar(), "3I2M3D")

        # local: unique path, bottom row to left-most column
        scomat = np.array([[0, 0, 0, 0, 0, 0],
                           [0, 0, 1, 0, 1, 0],
                           [0, 1, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0, 0],
                           [0, 1, 0, 0, 0, 0],
                           [0, 0, 2, 0, 1, 0]], dtype=float)
        obs = _traceback_one(5, 2, (scomat,), 0., 2., True, 2, 2, 2, 2)
        self.assertEqual(obs.to_cigar(), "2M")
        npt.assert_array_equal(obs.starts, [3, 0])

        # transpose: right-most column to top row
        scomat = np.ascontiguousarray(scomat.T)
        obs = _traceback_one(2, 5, (scomat,), 0., 2., True, 2, 2, 2, 2)
        self.assertEqual(obs.to_cigar(), "2M")
        npt.assert_array_equal(obs.starts, [0, 3])

        # GATCA vs CTTC, sub=(1, -1), gap=2 (linear)
        # local: unique path, right-most column to middle of the matrix
        scomat = np.array([[0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 0],
                           [0, 0, 1, 1, 0],
                           [0, 1, 0, 0, 2],
                           [0, 0, 0, 0, 0]], dtype=float)
        obs = _traceback_one(4, 4, (scomat,), 0., 2., True, 2, 2, 2, 2)
        self.assertEqual(obs.to_cigar(), "2M")
        npt.assert_array_equal(obs.starts, [2, 2])

        # ACGT vs AACTG, sub=(1, -1), gap=2, global
        # multiple paths: algorithm favors D > I > M in reverse order
        scomat = np.array([[  0,  -2,  -4,  -6,  -8, -10],
                           [ -2,   1,  -1,  -3,  -5,  -7],
                           [ -4,  -1,   0,   0,  -2,  -4],
                           [ -6,  -3,  -2,  -1,  -1,  -1],
                           [ -8,  -5,  -4,  -3,   0,  -2]], dtype=float)
        obs = _traceback_one(4, 5, (scomat,), 0., 2., local=False)
        self.assertEqual(obs.to_cigar(), "4M1I")

        # TAGCATC vs TCAGTC, sub=(2, -1), gap=(3, 1) (affine), global
        # involves tie, gap opening, gap extension
        NAN, INF = np.nan, -np.inf
        scomat = np.array([[  0,  -4,  -5,  -6,  -7,  -8,  -9],
                           [ -4,   2,  -2,  -3,  -4,  -5,  -6],
                           [ -5,  -2,   1,   0,  -4,  -5,  -6],
                           [ -6,  -3,  -3,   0,   2,  -2,  -3],
                           [ -7,  -4,  -1,  -4,  -1,   1,   0],
                           [ -8,  -5,  -5,   1,  -3,  -2,   0],
                           [ -9,  -6,  -6,  -3,   0,  -1,  -3],
                           [-10,  -7,  -4,  -4,  -4,  -1,   1]], dtype=float)
        insmat = np.array([[NAN, NAN, NAN, NAN, NAN, NAN, NAN],
                           [INF,  -8,  -2,  -3,  -4,  -5,  -6],
                           [INF,  -9,  -6,  -3,  -4,  -5,  -6],
                           [INF, -10,  -7,  -7,  -4,  -2,  -3],
                           [INF, -11,  -8,  -5,  -6,  -5,  -3],
                           [INF, -12,  -9,  -9,  -3,  -4,  -5],
                           [INF, -13, -10, -10,  -7,  -4,  -5],
                           [INF, -14, -11,  -8,  -8,  -8,  -5]], dtype=float)
        delmat = np.array([[NAN, INF, INF, INF, INF, INF, INF],
                           [NAN,  -8,  -9, -10, -11, -12, -13],
                           [NAN,  -2,  -6,  -7,  -8,  -9, -10],
                           [NAN,  -3,  -3,  -4,  -8,  -9, -10],
                           [NAN,  -4,  -4,  -4,  -2,  -6,  -7],
                           [NAN,  -5,  -5,  -5,  -3,  -3,  -4],
                           [NAN,  -6,  -6,  -3,  -4,  -4,  -4],
                           [NAN,  -7,  -7,  -4,  -4,  -5,  -5]], dtype=float)
        matrices = (scomat, insmat, delmat)
        obs = _traceback_one(7, 6, matrices, 3., 1.)
        self.assertEqual(obs.to_cigar(), "1M1I2M2D2M")

        # transpose (needs to transpose all 3 matrices and shuffle the latter 2)
        scomat = np.ascontiguousarray(scomat.T)
        insmat = np.ascontiguousarray(insmat.T)
        delmat = np.ascontiguousarray(delmat.T)
        matrices = (scomat, delmat, insmat)
        obs = _traceback_one(6, 7, matrices, 3., 1.)
        self.assertEqual(obs.to_cigar(), "1M2I2M1D2M")

        # GTCGG vs ATCGA, sub=(2, -1), gap=(3, 1), local
        # starts from middle, ends in middle of the matrix
        scomat = np.array([[0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 2, 0],
                           [0, 0, 2, 0, 0, 1],
                           [0, 0, 0, 4, 0, 0],
                           [0, 0, 0, 0, 6, 2],
                           [0, 0, 0, 0, 2, 5]], dtype=float)
        insmat = np.array([[NAN, NAN, NAN, NAN, NAN, NAN],
                           [INF,  -4,  -4,  -4,  -4,  -2],
                           [INF,  -4,  -4,  -2,  -3,  -4],
                           [INF,  -4,  -4,  -4,   0,  -1],
                           [INF,  -4,  -4,  -4,  -4,   2],
                           [INF,  -4,  -4,  -4,  -4,  -2]], dtype=float)
        delmat = np.array([[NAN, INF, INF, INF, INF, INF],
                           [NAN,  -4,  -4,  -4,  -4,  -4],
                           [NAN,  -4,  -4,  -4,  -2,  -4],
                           [NAN,  -4,  -2,  -4,  -3,  -3],
                           [NAN,  -4,  -3,   0,  -4,  -4],
                           [NAN,  -4,  -4,  -1,   2,  -2]], dtype=float)
        matrices = (scomat, insmat, delmat)
        obs = _traceback_one(4, 4, matrices, 3., 1., True, 2, 2, 2, 2)
        self.assertEqual(obs.to_cigar(), "3M")

        # floating-point tolerance
        scomat = np.array([[0, -2], [-2, -4]], dtype=float)
        args = (1, 1, (scomat,), 0., 2.,)
        self.assertEqual(_traceback_one(*args).to_cigar(), "1I1D")
        scomat[0, 1] = -2.0000000001
        self.assertEqual(_traceback_one(*args).to_cigar(), "1I1D")
        scomat[0, 1] = -2.0001
        self.assertEqual(_traceback_one(*args).to_cigar(), "1D1I")
        self.assertEqual(_traceback_one(*args, eps=1.0e-2).to_cigar(), "1I1D")

    def test_traceback_all(self):
        # linear, global, one stop, unique path
        seq1 = DNA("ACTCA")
        seq2 = DNA("CAGAG")
        (seq1, seq2), submat, _ = encode_sequences((seq1, seq2), (1, -1))
        query, target = np.ascontiguousarray(submat[seq1]), seq2

        scomat = np.array([[  0,  -2,  -4,  -6,  -8, -10],
                           [ -2,  -1,  -1,  -3,  -5,  -7],
                           [ -4,  -1,  -2,  -2,  -4,  -6],
                           [ -6,  -3,  -2,  -3,  -3,  -5],
                           [ -8,  -5,  -4,  -3,  -4,  -4],
                           [-10,  -7,  -4,  -5,  -2,  -4]])
        obs = _traceback_all([[5, 5]], None, (scomat,), query, target, 0, 2)
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0].to_cigar(), "1D4M1I")

        # one stop, three paths (ordering: D > I > M in reverse order)
        seq1 = DNA("ACGT")
        seq2 = DNA("AACTG")
        (seq1, seq2), submat, _ = encode_sequences((seq1, seq2), (1, -1))
        scomat = np.array([[  0,  -2,  -4,  -6,  -8, -10],
                           [ -2,   1,  -1,  -3,  -5,  -7],
                           [ -4,  -1,   0,   0,  -2,  -4],
                           [ -6,  -3,  -2,  -1,  -1,  -1],
                           [ -8,  -5,  -4,  -3,   0,  -2]])
        obs = _traceback_all([[4, 5]], None, (scomat,), submat[seq1], seq2, 0, 2)
        self.assertEqual(len(obs), 3)
        self.assertEqual(obs[0].to_cigar(), "4M1I")
        self.assertEqual(obs[1].to_cigar(), "1M1I3M")
        self.assertEqual(obs[2].to_cigar(), "1I4M")

        # limit number of paths to return
        obs = _traceback_all([[4, 5]], 2, (scomat,), submat[seq1], seq2, 0, 2)
        self.assertEqual(len(obs), 2)
        self.assertEqual(obs[0].to_cigar(), "4M1I")
        self.assertEqual(obs[1].to_cigar(), "1M1I3M")
        obs = _traceback_all([[4, 5]], 1, (scomat,), submat[seq1], seq2, 0, 2)
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0].to_cigar(), "4M1I")

        # semi-global, two stops, one path each
        scomat = np.array([[ 0,  0,  0,  0,  0,  0],
                           [ 0,  1,  1, -1, -1, -1],
                           [ 0, -1,  0,  2,  0, -2],
                           [ 0, -1, -2,  0,  1,  1],
                           [ 0, -1, -2, -2,  1,  0]])
        obs = _traceback_all(
            [[3, 5], [4, 4]], None, (scomat,), submat[seq1], seq2, 0, 2)
        self.assertEqual(len(obs), 2)
        self.assertEqual(obs[0].to_cigar(), "1I2M1I1M1D")
        self.assertEqual(obs[1].to_cigar(), "1I2M1D1M1I")

        # local, two stops, one path each
        seq1 = DNA("CGTC")
        seq2 = DNA("TCGAG")
        (seq1, seq2), submat, _ = encode_sequences((seq1, seq2), (1, -1))
        scomat = np.array([[0, 0, 0, 0, 0, 0],
                           [0, 0, 1, 0, 0, 0],
                           [0, 0, 0, 2, 0, 1],
                           [0, 1, 0, 0, 1, 0],
                           [0, 0, 2, 0, 0, 0]])
        obs = _traceback_all([[2, 3], [4, 2]], None, (scomat,), submat[seq1], seq2,
                             0, 2, True, 2, 2, 2, 2)
        self.assertEqual(len(obs), 2)
        self.assertEqual(obs[0].to_cigar(), "2M")
        self.assertEqual(obs[1].to_cigar(), "2M")

        # affine, one stop, two paths (involving both opening and extension)
        seq1 = DNA("TAGCATC")
        seq2 = DNA("TCAGTC")
        (seq1, seq2), submat, _ = encode_sequences((seq1, seq2), (2, -1))
        NAN, INF = np.nan, -np.inf
        scomat = np.array([[  0,  -4,  -5,  -6,  -7,  -8,  -9],
                           [ -4,   2,  -2,  -3,  -4,  -5,  -6],
                           [ -5,  -2,   1,   0,  -4,  -5,  -6],
                           [ -6,  -3,  -3,   0,   2,  -2,  -3],
                           [ -7,  -4,  -1,  -4,  -1,   1,   0],
                           [ -8,  -5,  -5,   1,  -3,  -2,   0],
                           [ -9,  -6,  -6,  -3,   0,  -1,  -3],
                           [-10,  -7,  -4,  -4,  -4,  -1,   1]])
        insmat = np.array([[NAN, NAN, NAN, NAN, NAN, NAN, NAN],
                           [INF,  -8,  -2,  -3,  -4,  -5,  -6],
                           [INF,  -9,  -6,  -3,  -4,  -5,  -6],
                           [INF, -10,  -7,  -7,  -4,  -2,  -3],
                           [INF, -11,  -8,  -5,  -6,  -5,  -3],
                           [INF, -12,  -9,  -9,  -3,  -4,  -5],
                           [INF, -13, -10, -10,  -7,  -4,  -5],
                           [INF, -14, -11,  -8,  -8,  -8,  -5]])
        delmat = np.array([[NAN, INF, INF, INF, INF, INF, INF],
                           [NAN,  -8,  -9, -10, -11, -12, -13],
                           [NAN,  -2,  -6,  -7,  -8,  -9, -10],
                           [NAN,  -3,  -3,  -4,  -8,  -9, -10],
                           [NAN,  -4,  -4,  -4,  -2,  -6,  -7],
                           [NAN,  -5,  -5,  -5,  -3,  -3,  -4],
                           [NAN,  -6,  -6,  -3,  -4,  -4,  -4],
                           [NAN,  -7,  -7,  -4,  -4,  -5,  -5]])
        matrices = [scomat, insmat, delmat]
        obs = _traceback_all([[7, 6]], None, matrices, submat[seq1], seq2, 3, 1)
        self.assertEqual(len(obs), 2)
        self.assertEqual(obs[0].to_cigar(), "1M1I2M2D2M")
        self.assertEqual(obs[1].to_cigar(), "1M2D2M1I2M")

        # transpose
        submat = np.ascontiguousarray(submat.T)
        scomat = np.ascontiguousarray(scomat.T)
        insmat = np.ascontiguousarray(insmat.T)
        delmat = np.ascontiguousarray(delmat.T)
        matrices = [scomat, delmat, insmat]
        obs = _traceback_all([[6, 7]], None, matrices, submat[seq2], seq1, 3, 1)
        self.assertEqual(len(obs), 2)
        self.assertEqual(obs[0].to_cigar(), "1M2I2M1D2M")
        self.assertEqual(obs[1].to_cigar(), "1M1D2M2I2M")

        # local, unique path, middle to middle
        seq1 = DNA("GTCGG")
        seq2 = DNA("ATCGA")
        (seq1, seq2), submat, _ = encode_sequences((seq1, seq2), (2, -1))
        scomat = np.array([[0, 0, 0, 0, 0, 0],
                           [0, 0, 0, 0, 2, 0],
                           [0, 0, 2, 0, 0, 1],
                           [0, 0, 0, 4, 0, 0],
                           [0, 0, 0, 0, 6, 2],
                           [0, 0, 0, 0, 2, 5]])
        insmat = np.array([[NAN, NAN, NAN, NAN, NAN, NAN],
                           [INF,  -4,  -4,  -4,  -4,  -2],
                           [INF,  -4,  -4,  -2,  -3,  -4],
                           [INF,  -4,  -4,  -4,   0,  -1],
                           [INF,  -4,  -4,  -4,  -4,   2],
                           [INF,  -4,  -4,  -4,  -4,  -2]])
        delmat = np.array([[NAN, INF, INF, INF, INF, INF],
                           [NAN,  -4,  -4,  -4,  -4,  -4],
                           [NAN,  -4,  -4,  -4,  -2,  -4],
                           [NAN,  -4,  -2,  -4,  -3,  -3],
                           [NAN,  -4,  -3,   0,  -4,  -4],
                           [NAN,  -4,  -4,  -1,   2,  -2]])
        matrices = [scomat, insmat, delmat]
        obs = _traceback_all(
            [[4, 4]], None, matrices, submat[seq1], seq2, 3, 1, True, 2, 2, 2, 2)
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0].to_cigar(), "3M")

        # floating-point tolerance
        seq1 = DNA("A")
        seq2 = DNA("A")
        (seq1, seq2), submat, _ = encode_sequences((seq1, seq2), (2, -1))
        scomat = np.array([[1, 5], [5, 3]], dtype=float)
        args = ([[1, 1]], None, (scomat,), submat[seq1], seq2, 0., 2.)
        obs = _traceback_all(*args)
        self.assertEqual(len(obs), 3)
        for i, exp in enumerate(["1I1D", "1D1I", "1M"]):
            self.assertEqual(obs[i].to_cigar(), exp)

        scomat[0, 1] -= 0.00000001
        scomat[1, 0] += 0.00000001
        obs = _traceback_all(*args)
        self.assertEqual(len(obs), 3)
        for i, exp in enumerate(["1I1D", "1D1I", "1M"]):
            self.assertEqual(obs[i].to_cigar(), exp)

        scomat[0, 1] -= 0.0001
        scomat[1, 0] += 0.0001
        obs = _traceback_all(*args)
        self.assertEqual(len(obs), 1)
        self.assertEqual(obs[0].to_cigar(), "1M")

        obs = _traceback_all(*args, eps=1.0e-2)
        self.assertEqual(len(obs), 3)
        for i, exp in enumerate(["1I1D", "1D1I", "1M"]):
            self.assertEqual(obs[i].to_cigar(), exp)


if __name__ == "__main__":
    unittest.main()
