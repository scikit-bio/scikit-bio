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
from skbio.alignment._utils import (
    _prep_seqs_submat,
    _prep_gap_cost,
    _parse_seqs_with_submat,
    _parse_seqs_alone,
    _prep_idmat,
    _check_same_type,
)


class UtilsTests(unittest.TestCase):
    def test_prep_seqs_submat(self):
        # Note: most situations have been covered by other tests.
        # substitution matrix instance
        sm = SubstitutionMatrix.by_name("NUC.4.4")
        seqs = (DNA("GAATCC"),
                DNA("AGATCT"),
                DNA("CCGG"))
        obs = _prep_seqs_submat(seqs, sm)
        exp = [[2, 0, 0, 1, 3, 3],
               [0, 2, 0, 1, 3, 1],
               [3, 3, 2, 2]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertIs(obs[1], sm._data)
        self.assertEqual(obs[1].shape[0], 15)

        # substitution matrix name
        obs = _prep_seqs_submat(seqs, "NUC.4.4")
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertIs(obs[1], sm._data)
        self.assertEqual(obs[1].shape[0], 15)

        # match/mismatch scores
        obs = _prep_seqs_submat(seqs, (1, -1))
        exp = [[6, 2, 2, 13, 4, 4],
               [2, 6, 2, 13, 4, 13],
               [4, 4, 6, 6]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1].shape[0], 17)
        self.assertEqual(obs[1][2, 2], 1.0)
        self.assertEqual(obs[1][3, 1], -1.0)

    def test_prep_gap_cost(self):
        obs = _prep_gap_cost(3)
        self.assertTupleEqual(obs, (0.0, 3.0))
        self.assertIsInstance(obs[0], float)
        self.assertIsInstance(obs[1], float)

        obs = _prep_gap_cost((2, 5))
        self.assertTupleEqual(obs, (2.0, 5.0))
        self.assertIsInstance(obs[0], float)
        self.assertIsInstance(obs[1], float)

        obs = _prep_gap_cost((0.5, -2.5))
        self.assertTupleEqual(obs, (0.5, -2.5))
        self.assertIsInstance(obs[0], float)
        self.assertIsInstance(obs[1], float)

    def test_parse_seqs_with_submat(self):
        sm = SubstitutionMatrix.by_name("NUC.4.4")

        # DNA sequences
        seqs = (DNA("GAATCC"),
                DNA("AGATCT"),
                DNA("CCGG"))
        obs = _parse_seqs_with_submat(seqs, sm)
        exp = [[2, 0, 0, 1, 3, 3],
               [0, 2, 0, 1, 3, 1],
               [3, 3, 2, 2]]
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # non-grammared sequences
        seqs = (Sequence("GAATCC"),
                Sequence("AGATCT"),
                Sequence("CCGG"))
        obs = _parse_seqs_with_submat(seqs, sm)
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # raw strings
        seqs = ("GAATCC",
                "AGATCT",
                "CCGG")
        obs = _parse_seqs_with_submat(seqs, sm)
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # degenerate code becomes wildcard ("N")
        sm = SubstitutionMatrix.identity("ACGTN", 1, -1)
        seqs = (DNA("CGGRATCCA"),
                DNA("GGAATYCT"))
        obs = _parse_seqs_with_submat(seqs, sm)
        exp = [[1, 2, 2, 4, 0, 3, 1, 1, 0],
               [2, 2, 0, 0, 3, 4, 1, 3]]
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # won't work if wildcard is not defined
        msg = "One or more characters in the sequence are absent from the alphabet."
        seqs = (Sequence("CGGRATCCA"),
                Sequence("GGAATYCT"))
        with self.assertRaises(ValueError) as cm:
            _ = _parse_seqs_with_submat(seqs, sm)
        self.assertEqual(str(cm.exception), msg)

        # won't work if the substitution matrix doesn't contain "N".
        msg = 'Wildcard character "N" is not in the alphabet.'
        sm = SubstitutionMatrix.identity("ACGT", 1, -1)
        seqs = (DNA("CGGRATCCA"),
                DNA("GGAATYCT"))
        with self.assertRaises(ValueError) as cm:
            _ = _parse_seqs_with_submat(seqs, sm)
        self.assertEqual(str(cm.exception), msg)

        # Unicode sequences against ASCII matrix
        self.assertTrue(sm.is_ascii)
        seqs = ("あいう", "あうえお")
        msg = "Cannot encode ASCII characters."
        with self.assertRaises(ValueError) as cm:
            _ = _parse_seqs_with_submat(seqs, sm)
        self.assertEqual(str(cm.exception), msg)

        # Unicode sequences against Unicode matrix
        sm = SubstitutionMatrix.identity("あいうえお", 1, -1)
        self.assertFalse(sm.is_ascii)
        obs = _parse_seqs_with_submat(seqs, sm)
        exp = [[0, 1, 2], [0, 2, 3, 4]]
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # with missing characters
        msg = "One or more characters in the sequence are absent from the alphabet."
        seqs = ("かさた", "あうお")
        with self.assertRaises(ValueError) as cm:
            _ = _parse_seqs_with_submat(seqs, sm)
        self.assertEqual(str(cm.exception), msg)

    def test_parse_seqs_alone(self):
        # grammared sequences
        seqs = (DNA("GAATCC"),
                DNA("AGATCT"),
                DNA("CCGG"))
        obs = _parse_seqs_alone(seqs)
        exp = [[71, 65, 65, 84, 67, 67],
               [65, 71, 65, 84, 67, 84],
               [67, 67, 71, 71]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertIs(obs[1], DNA)

        # non-grammared sequences
        seqs = (Sequence("GAATCC"),
                Sequence("AGATCT"),
                Sequence("CCGG"))
        obs = _parse_seqs_alone(seqs)
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1], "ascii")

        # raw strings
        seqs = ("GAATCC",
                "AGATCT",
                "CCGG")
        obs = _parse_seqs_alone(seqs)
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1], "ascii")

        # Unicode strings
        seqs = ("あいう", "あうえお")
        exp = [[0, 1, 2], [0, 2, 3, 4]]
        obs = _parse_seqs_alone(seqs)
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1], 5)

        # lists of words
        seq1 = "lorem ipsum sit amet tempor".split()
        seq2 = "ipsum sit dolor sed eiusmod".split()
        obs = _parse_seqs_alone([seq1, seq2])
        exp = [[4, 3, 6, 0, 7], [3, 6, 1, 5, 2]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1], 8)

        # arrays of numbers that are within ASCII range
        seq1 = np.array([9, 4, 4, 2, 3])
        seq2 = np.array([5, 2, 4, 3, 3, 3, 7])
        obs = _parse_seqs_alone([seq1, seq2])
        npt.assert_array_equal(obs[0][0], seq1)
        npt.assert_array_equal(obs[0][1], seq2)
        self.assertEqual(obs[1], "ascii")

        # arrays of numbers that are outside ASCII range
        seq1 = np.array([100, 150, 200, 200, 300])
        seq2 = np.array([150, 300, 300, 200, 100, 350])
        obs = _parse_seqs_alone([seq1, seq2])
        exp = [[0, 1, 2, 2, 3], [1, 3, 3, 2, 0, 4]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1], 5)

    def test_prep_idmat(self):
        # Grammared (DNA) matrix
        seqs = [np.array([71, 65, 65, 84, 67, 67]),
                np.array([65, 71, 65, 84, 67, 84]),
                np.array([67, 67, 71, 71])]
        obs = _prep_idmat(seqs, DNA, 1.0, -1.0)
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
        obs = _prep_idmat(seqs, DNA, 2.0, -3.0)
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1].shape[0], 17)
        self.assertEqual(obs[1][2, 2], 2.0)
        self.assertEqual(obs[1][3, 1], -3.0)
        self.assertIs(obs[1], submat)

        # ASCII matrix
        seqs = [np.array([9, 4, 4, 2, 3], dtype=np.uint8),
                np.array([5, 2, 4, 3, 3, 3, 7], dtype=np.uint8)]
        obs = _prep_idmat(seqs, "ascii", 2.0, -3.0)
        for o, e in zip(obs[0], seqs):
            npt.assert_array_equal(o, e)
            self.assertIs(o.dtype.type, np.intp)
        self.assertEqual(obs[1].shape[0], 128)
        self.assertEqual(obs[1][5, 5], 2.0)
        self.assertEqual(obs[1][7, 8], -3.0)
        submat = obs[1]

        # update cached matrix
        obs = _prep_idmat(seqs, "ascii", 1.0, -2.5)
        for o, e in zip(obs[0], seqs):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1][5, 5], 1.0)
        self.assertEqual(obs[1][7, 8], -2.5)
        self.assertIs(obs[1], submat)

        # arbitrary matrix
        seqs = [np.array([0, 1, 2, 2, 3]),
                np.array([1, 3, 3, 2, 0, 4])]
        obs = _prep_idmat(seqs, 5, 1.0, -2.5)
        for o, e in zip(obs[0], seqs):
            npt.assert_array_equal(o, e)
            self.assertIs(o.dtype.type, np.intp)
        self.assertEqual(obs[1].shape[0], 5)
        self.assertEqual(obs[1][2, 2], 1.0)
        self.assertEqual(obs[1][3, 1], -2.5)
        submat = obs[1]

        # create every time (no update)
        obs = _prep_idmat(seqs, 5, 3.5, -4.0)
        self.assertEqual(obs[1][2, 2], 3.5)
        self.assertEqual(obs[1][3, 1], -4.0)
        self.assertIsNot(obs[1], submat)

    def test_check_same_type(self):
        obs = _check_same_type([1, 2, 3, 4])
        self.assertIs(obs, int)
        obs = _check_same_type(["abc", "def"])
        self.assertIs(obs, str)
        obs = _check_same_type([DNA("CGT"), DNA("AAT"), DNA("GC")])
        self.assertIs(obs, DNA)
        msg = "Variables are of different types."
        with self.assertRaises(TypeError) as cm:
            _ = _check_same_type(["a", 1])
        self.assertEqual(str(cm.exception), msg)
        with self.assertRaises(TypeError) as cm:
            _ = _check_same_type([DNA("CGT"), Protein("MKY")])
        self.assertEqual(str(cm.exception), msg)


if __name__ == "__main__":
    unittest.main()
