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

from skbio.sequence import Sequence, DNA, Protein, SubstitutionMatrix
from skbio.alignment import TabularMSA
from skbio.alignment._utils import (
    encode_sequences,
    encode_alignment,
    prep_gapcost,
    prep_identity_matrix,
    _check_seqtype,
)


class UtilsTests(unittest.TestCase):

    def test_encode_sequences_submat(self):
        """Encode sequences with a substitution matrix."""
        # nucleotide substitution matrix
        sm = SubstitutionMatrix.by_name("NUC.4.4")

        # DNA sequences
        seqs = (DNA("GAATCC"),
                DNA("AGATCT"),
                DNA("CCGG"))
        obs = encode_sequences(seqs, sm)
        exp = [[2, 0, 0, 1, 3, 3],
               [0, 2, 0, 1, 3, 1],
               [3, 3, 2, 2]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertIs(obs[1], sm._data)

        # non-grammared sequences
        seqs = (Sequence("GAATCC"),
                Sequence("AGATCT"),
                Sequence("CCGG"))
        obs = encode_sequences(seqs, sm)
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)

        # raw strings
        seqs = ("GAATCC",
                "AGATCT",
                "CCGG")
        obs = encode_sequences(seqs, sm)
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)

        # degenerate code becomes wildcard ("N")
        sm = SubstitutionMatrix.identity("ACGTN", 1, -1)
        seqs = (DNA("CGGRATCCA"),
                DNA("GGAATYCT"))
        obs = encode_sequences(seqs, sm)
        exp = [[1, 2, 2, 4, 0, 3, 1, 1, 0],
               [2, 2, 0, 0, 3, 4, 1, 3]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertIs(obs[1], sm._data)

        # won't work if wildcard is not defined
        msg = "Sequence 1 contain character(s) absent from the substitution matrix."
        seqs = (Sequence("CGGRATCCA"),
                Sequence("GGAATYCT"))
        with self.assertRaises(ValueError) as cm:
            _ = encode_sequences(seqs, sm)
        self.assertEqual(str(cm.exception), msg)

        # won't work if the substitution matrix doesn't contain "N".
        sm = SubstitutionMatrix.identity("ACGT", 1, -1)
        seqs = (DNA("CGGRATCCA"),
                DNA("GGAATYCT"))
        with self.assertRaises(ValueError) as cm:
            _ = encode_sequences(seqs, sm)
        self.assertEqual(str(cm.exception), msg)

        # Unicode sequences against ASCII matrix
        self.assertTrue(sm.is_ascii)
        seqs = ("①②③", "①③④⑤")
        msg = ("Substitution matrix has an ASCII alphabet, but sequences cannot be "
               "fully encoded into ASCII codes.")
        with self.assertRaises(ValueError) as cm:
            _ = encode_sequences(seqs, sm)
        self.assertEqual(str(cm.exception), msg)

        # Unicode sequences against Unicode matrix
        sm = SubstitutionMatrix.identity("①②③④⑤", 1, -1)
        self.assertFalse(sm.is_ascii)
        obs = encode_sequences(seqs, sm)
        exp = [[0, 1, 2], [0, 2, 3, 4]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertIs(obs[1], sm._data)

        # with missing characters
        msg = "Sequence 2 contain character(s) absent from the substitution matrix."
        seqs = ("①③⑤", "④ⓐ②")
        with self.assertRaises(ValueError) as cm:
            _ = encode_sequences(seqs, sm)
        self.assertEqual(str(cm.exception), msg)

        # ASCII sequences against Unicode matrix
        alphabet = "ACGT" + chr(100) + chr(200) + chr(300)
        sm = SubstitutionMatrix.identity(alphabet, 1, -1)
        seqs = (DNA("GAATCC"),
                DNA("AGATCT"))
        obs = encode_sequences(seqs, sm)
        exp = [[2, 0, 0, 3, 1, 1],
               [0, 2, 0, 3, 1, 3]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertIs(obs[1], sm._data)

    def test_encode_sequences_edit(self):
        """Encode sequences with match/mismatch scores."""
        match, mismatch = 1, -1

        # grammared sequences: submat will cover their alphabet
        seqs = (DNA("GAATCC"),
                DNA("AGATCT"),
                DNA("CCGG"))
        obs = encode_sequences(seqs, (match, mismatch))
        exp = [[6, 2, 2, 13, 4, 4],
               [2, 6, 2, 13, 4, 13],
               [4, 4, 6, 6]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1].shape[0], 17)
        self.assertEqual(obs[1][0, 0], match)
        self.assertEqual(obs[1][0, 1], mismatch)

        # non-grammared sequences: submat will cover all ASCII codes
        seqs = (Sequence("GAATCC"),
                Sequence("AGATCT"),
                Sequence("CCGG"))
        obs = encode_sequences(seqs, (match, mismatch))
        exp = [[71, 65, 65, 84, 67, 67],
               [65, 71, 65, 84, 67, 84],
               [67, 67, 71, 71]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1].shape[0], 128)
        submat = obs[1]

        # raw strings: same as above, plus the same submat will be reloaded from cache
        seqs = ("GAATCC",
                "AGATCT",
                "CCGG")
        obs = encode_sequences(seqs, (match, mismatch))
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertIs(submat, obs[1])

        # Unicode strings: submat will cover unique elements
        seqs = ("①②③", "①③④⑤")
        exp = [[0, 1, 2], [0, 2, 3, 4]]
        obs = encode_sequences(seqs, (match, mismatch))
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1].shape[0], 5)
        self.assertEqual(obs[1][0, 0], match)
        self.assertEqual(obs[1][0, 1], mismatch)
        self.assertIsNot(submat, obs[1])

        # lists of words
        seq1 = "lorem ipsum sit amet tempor".split()
        seq2 = "ipsum sit dolor sed eiusmod".split()
        obs = encode_sequences([seq1, seq2], (match, mismatch))
        exp = [[4, 3, 6, 0, 7], [3, 6, 1, 5, 2]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1].shape[0], 8)

        # arrays of numbers that are within ASCII range
        seq1 = np.array([9, 4, 4, 2, 3])
        seq2 = np.array([5, 2, 4, 3, 3, 3, 7])
        obs = encode_sequences([seq1, seq2], (match, mismatch))
        npt.assert_array_equal(obs[0][0], seq1)
        npt.assert_array_equal(obs[0][1], seq2)
        self.assertEqual(obs[1].shape[0], 128)

        # arrays of numbers that are outside ASCII range
        seq1 = np.array([100, 150, 200, 200, 300])
        seq2 = np.array([150, 300, 300, 200, 100, 350])
        obs = encode_sequences([seq1, seq2], (match, mismatch))
        exp = [[0, 1, 2, 2, 3], [1, 3, 3, 2, 0, 4]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1].shape[0], 5)

    def test_encode_sequences_error(self):
        msg = "Sequences are of different types."
        with self.assertRaises(TypeError) as cm:
            _ = encode_sequences([DNA("CGT"), Protein("MKY")], (1, -1))
        self.assertEqual(str(cm.exception), msg)

        msg = "No sequence is provided."
        with self.assertRaises(ValueError) as cm:
            _ = encode_sequences([], (1, -1))
        self.assertEqual(str(cm.exception), msg)

        msg = "Sequence {} has a length of zero."
        with self.assertRaises(ValueError) as cm:
            _ = encode_sequences(["", "ATG"], (1, -1))
        self.assertEqual(str(cm.exception), msg.format(1))
        with self.assertRaises(ValueError) as cm:
            _ = encode_sequences(["MKY", ""], (1, -1))
        self.assertEqual(str(cm.exception), msg.format(2))

    def test_encode_alignment_aligned(self):
        # Note: the unit-tests of this function aren't adequate. It still has 100%
        # coverage because all conditions have been tested in upstream and downstream
        # functions. But it would be good to enrich tests in the future.
        match, mismatch = 1, -1

        # raw sequences
        seqs = ["GCG-TCAGT",
                "AGG-TCATT",
                "ACGTTC--T"]
        obs = encode_alignment(seqs, (match, mismatch))
        exp0 = np.array([[71, 67, 71, 45, 84, 67, 65, 71, 84],
                         [65, 71, 71, 45, 84, 67, 65, 84, 84],
                         [65, 67, 71, 84, 84, 67, 45, 45, 84]])
        npt.assert_array_equal(obs[0], exp0)
        self.assertEqual(obs[1].shape[0], 128)
        exp2 = np.array([[0, 1, 0, 0, 0],
                         [0, 1, 0, 0, 0],
                         [0, 0, 0, 1, 0]]).astype(bool)
        npt.assert_array_equal(obs[2], exp2)
        exp3 = np.array([3, 1, 2, 2, 1])
        npt.assert_array_equal(obs[3], exp3)

        # grammared sequences
        seqs = list(map(DNA, seqs))
        obs = encode_alignment(seqs, (match, mismatch))
        exp0 = np.array([[6,  4,  6,  0, 13,  4,  2,  6, 13],
                         [2,  6,  6,  0, 13,  4,  2, 13, 13],
                         [2,  4,  6, 13, 13,  4,  0,  0, 13]])
        npt.assert_array_equal(obs[0], exp0)
        self.assertEqual(obs[1].shape[0], 17)

        # tabular alignment
        msa = TabularMSA(map(DNA, seqs))
        obs = encode_alignment(msa, (match, mismatch))
        npt.assert_array_equal(obs[0], exp0)
        self.assertEqual(obs[1].shape[0], 17)

        # Unicode characters
        seqs = ["ïëï-öëäïö",
                "äïï-öëäöö",
                "äëïööë--ö"]
        obs = encode_alignment(seqs, (match, mismatch))
        exp0 = np.array([[3, 2, 3, 0, 4, 2, 1, 3, 4],
                         [1, 3, 3, 0, 4, 2, 1, 4, 4],
                         [1, 2, 3, 4, 4, 2, 0, 0, 4]])
        npt.assert_array_equal(obs[0], exp0)

        # substitution matrix
        alphabet = "ACGT"
        sm = SubstitutionMatrix.identity(alphabet, 1, -1)
        seqs = (DNA("GAA-CC"),
                DNA("AGATCT"))
        obs = encode_alignment(seqs, sm)
        exp = [[2, 0, 0, -1, 1, 1],
               [0, 2, 0, 3, 1, 3]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertIs(obs[1], sm._data)

        # ASCII sequences against Unicode matrix
        alphabet += chr(100) + chr(200) + chr(300)
        sm = SubstitutionMatrix.identity(alphabet, 1, -1)
        obs = encode_alignment(seqs, sm)
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)

    def test_prep_gapcost(self):
        obs = prep_gapcost(3)
        self.assertTupleEqual(obs, (0.0, 3.0))
        self.assertIsInstance(obs[0], np.float32)
        self.assertIsInstance(obs[1], np.float32)

        obs = prep_gapcost((2, 5))
        self.assertTupleEqual(obs, (2.0, 5.0))
        self.assertIsInstance(obs[0], np.float32)
        self.assertIsInstance(obs[1], np.float32)

        obs = prep_gapcost((0.5, 2.5))
        self.assertTupleEqual(obs, (0.5, 2.5))
        self.assertIsInstance(obs[0], np.float32)
        self.assertIsInstance(obs[1], np.float32)

        obs = prep_gapcost((0.5, 2.5), dtype=np.float64)
        self.assertTupleEqual(obs, (0.5, 2.5))
        self.assertIsInstance(obs[0], np.float64)
        self.assertIsInstance(obs[1], np.float64)

        obs = prep_gapcost((np.float64(0.5), np.float64(2.5)))
        self.assertTupleEqual(obs, (0.5, 2.5))
        self.assertIsInstance(obs[0], np.float32)
        self.assertIsInstance(obs[1], np.float32)

        obs = prep_gapcost((np.float32(0.5), np.float32(2.5)), dtype=np.float64)
        self.assertTupleEqual(obs, (0.5, 2.5))
        self.assertIsInstance(obs[0], np.float64)
        self.assertIsInstance(obs[1], np.float64)

    def testprep_identity_matrix(self):
        # Grammared (DNA) matrix
        seqs = [np.array([71, 65, 65, 84, 67, 67]),
                np.array([65, 71, 65, 84, 67, 84]),
                np.array([67, 67, 71, 71])]
        obs = prep_identity_matrix(seqs, DNA, 1.0, -1.0)
        exp = [np.array([6, 2, 2, 13, 4, 4]),
                np.array([2, 6, 2, 13, 4, 13]),
                np.array([4, 4, 6, 6])]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1].shape[0], 17)
        self.assertEqual(obs[1][2, 2], 1.0)
        self.assertEqual(obs[1][3, 1], -1.0)
        submat = obs[1]

        # update cached matrix without recreating it
        obs = prep_identity_matrix(seqs, DNA, 2.0, -3.0)
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1].shape[0], 17)
        self.assertEqual(obs[1][2, 2], 2.0)
        self.assertEqual(obs[1][3, 1], -3.0)
        self.assertIs(obs[1], submat)

        # ASCII matrix
        seqs = [np.array([9, 4, 4, 2, 3], dtype=np.uint8),
                np.array([5, 2, 4, 3, 3, 3, 7], dtype=np.uint8)]
        obs = prep_identity_matrix(seqs, "ascii", 2.0, -3.0)
        for o, e in zip(obs[0], seqs):
            npt.assert_array_equal(o, e)
            self.assertIs(o.dtype.type, np.intp)
        self.assertEqual(obs[1].shape[0], 128)
        self.assertEqual(obs[1][5, 5], 2.0)
        self.assertEqual(obs[1][7, 8], -3.0)
        submat = obs[1]

        # update cached matrix
        obs = prep_identity_matrix(seqs, "ascii", 1.0, -2.5)
        for o, e in zip(obs[0], seqs):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1][5, 5], 1.0)
        self.assertEqual(obs[1][7, 8], -2.5)
        self.assertIs(obs[1], submat)

        # arbitrary matrix
        seqs = [np.array([0, 1, 2, 2, 3]),
                np.array([1, 3, 3, 2, 0, 4])]
        obs = prep_identity_matrix(seqs, 5, 1.0, -2.5)
        for o, e in zip(obs[0], seqs):
            npt.assert_array_equal(o, e)
            self.assertIs(o.dtype.type, np.intp)
        self.assertEqual(obs[1].shape[0], 5)
        self.assertEqual(obs[1][2, 2], 1.0)
        self.assertEqual(obs[1][3, 1], -2.5)
        submat = obs[1]

        # create every time (no update)
        obs = prep_identity_matrix(seqs, 5, 3.5, -4.0)
        self.assertEqual(obs[1][2, 2], 3.5)
        self.assertEqual(obs[1][3, 1], -4.0)
        self.assertIsNot(obs[1], submat)

    def test_check_seqtype(self):
        obs = _check_seqtype([1, 2, 3, 4])
        self.assertIs(obs, int)
        obs = _check_seqtype(["abc", "def"])
        self.assertIs(obs, str)
        obs = _check_seqtype([DNA("CGT"), DNA("AAT"), DNA("GC")])
        self.assertIs(obs, DNA)
        obs = _check_seqtype(TabularMSA(list(map(DNA, ["CGT", "AAT", "GCA"]))))
        self.assertIs(obs, DNA)

        msg = "Sequences are of different types."
        with self.assertRaises(TypeError) as cm:
            _ = _check_seqtype(["a", 1])
        self.assertEqual(str(cm.exception), msg)
        with self.assertRaises(TypeError) as cm:
            _ = _check_seqtype([DNA("CGT"), Protein("MKY")])
        self.assertEqual(str(cm.exception), msg)

        msg = "No sequence is provided."
        with self.assertRaises(ValueError) as cm:
            _ = _check_seqtype([])
        self.assertEqual(str(cm.exception), msg)
        with self.assertRaises(ValueError) as cm:
            _ = _check_seqtype(TabularMSA([]))
        self.assertEqual(str(cm.exception), msg)


if __name__ == "__main__":
    unittest.main()
