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
from skbio.sequence import Sequence, DNA, Protein, SubstitutionMatrix
from skbio.alignment._utils import (
    _check_same_type,
    _parse_seqs,
    _parse_seqs_submat,
    _prep_seqs_submat,
    _prep_gap_cost,
)


class UtilsTests(unittest.TestCase):
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

    def test_parse_seqs(self):
        # grammared sequences
        seqs = (DNA("GAATCC"),
                DNA("AGATCT"),
                DNA("CCGG"))
        obs = _parse_seqs(seqs)
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
        obs = _parse_seqs(seqs)
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1], "ascii")

        # raw strings
        seqs = ("GAATCC",
                "AGATCT",
                "CCGG")
        obs = _parse_seqs(seqs)
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1], "ascii")

        # Unicode strings
        seqs = ("あいう", "あうえお")
        exp = [[0, 1, 2], [0, 2, 3, 4]]
        obs = _parse_seqs(seqs)
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1], 5)

        # lists of words
        seq1 = "lorem ipsum sit amet tempor".split()
        seq2 = "ipsum sit dolor sed eiusmod".split()
        obs = _parse_seqs([seq1, seq2])
        exp = [[4, 3, 6, 0, 7], [3, 6, 1, 5, 2]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1], 8)

        # arrays of numbers that are within ASCII range
        seq1 = np.array([9, 4, 4, 2, 3])
        seq2 = np.array([5, 2, 4, 3, 3, 3, 7])
        obs = _parse_seqs([seq1, seq2])
        npt.assert_array_equal(obs[0][0], seq1)
        npt.assert_array_equal(obs[0][1], seq2)
        self.assertEqual(obs[1], "ascii")

        # arrays of numbers that are outside ASCII range
        seq1 = np.array([100, 150, 200, 200, 300])
        seq2 = np.array([150, 300, 300, 200, 100, 350])
        obs = _parse_seqs([seq1, seq2])
        exp = [[0, 1, 2, 2, 3], [1, 3, 3, 2, 0, 4]]
        for o, e in zip(obs[0], exp):
            npt.assert_array_equal(o, e)
        self.assertEqual(obs[1], 5)

    def test_parse_seqs_submat(self):
        sm = SubstitutionMatrix.by_name("NUC.4.4")

        # DNA sequences
        seqs = (DNA("GAATCC"),
                DNA("AGATCT"),
                DNA("CCGG"))
        obs = _parse_seqs_submat(seqs, sm)
        exp = [[2, 0, 0, 1, 3, 3],
               [0, 2, 0, 1, 3, 1],
               [3, 3, 2, 2]]
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # non-grammared sequences
        seqs = (Sequence("GAATCC"),
                Sequence("AGATCT"),
                Sequence("CCGG"))
        obs = _parse_seqs_submat(seqs, sm)
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # raw strings
        seqs = ("GAATCC",
                "AGATCT",
                "CCGG")
        obs = _parse_seqs_submat(seqs, sm)
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # degenerate code becomes wildcard ("N")
        sm = SubstitutionMatrix.identity("ACGTN", 1, -1)
        seqs = (DNA("CGGRATCCA"),
                DNA("GGAATYCT"))
        obs = _parse_seqs_submat(seqs, sm)
        exp = [[1, 2, 2, 4, 0, 3, 1, 1, 0],
               [2, 2, 0, 0, 3, 4, 1, 3]]
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # won't work if wildcard is not defined
        msg = "One or more characters in the sequence are absent from the alphabet."
        seqs = (Sequence("CGGRATCCA"),
                Sequence("GGAATYCT"))
        with self.assertRaises(ValueError) as cm:
            _ = _parse_seqs_submat(seqs, sm)
        self.assertEqual(str(cm.exception), msg)

        # won't work if the substitution matrix doesn't contain "N".
        msg = 'Wildcard character "N" is not in the alphabet.'
        sm = SubstitutionMatrix.identity("ACGT", 1, -1)
        seqs = (DNA("CGGRATCCA"),
                DNA("GGAATYCT"))
        with self.assertRaises(ValueError) as cm:
            _ = _parse_seqs_submat(seqs, sm)
        self.assertEqual(str(cm.exception), msg)

        # Unicode sequences against ASCII matrix
        self.assertTrue(sm.is_ascii)
        seqs = ("あいう", "あうえお")
        msg = "Cannot encode ASCII characters."
        with self.assertRaises(ValueError) as cm:
            _ = _parse_seqs_submat(seqs, sm)
        self.assertEqual(str(cm.exception), msg)

        # Unicode sequences against Unicode matrix
        sm = SubstitutionMatrix.identity("あいうえお", 1, -1)
        self.assertFalse(sm.is_ascii)
        obs = _parse_seqs_submat(seqs, sm)
        exp = [[0, 1, 2], [0, 2, 3, 4]]
        for o, e in zip(obs, exp):
            npt.assert_array_equal(o, e)

        # with missing characters
        msg = "One or more characters in the sequence are absent from the alphabet."
        seqs = ("かさた", "あうお")
        with self.assertRaises(ValueError) as cm:
            _ = _parse_seqs_submat(seqs, sm)
        self.assertEqual(str(cm.exception), msg)


    # >>> tokens = 'lorem ipsum dolor sit amet'.split()
    # >>> sm = SubstitutionMatrix(tokens, np.array([
    # ...     [3, 1, 2, 0, 0],
    # ...     [1, 3, 1, 2, 0],
    # ...     [2, 1, 3, 1, 2],
    # ...     [0, 2, 1, 3, 1],
    # ...     [0, 0, 2, 1, 3]]))
    # >>> sm.alphabet
    # ('lorem', 'ipsum', 'dolor', 'sit', 'amet')
    # >>> sm['lorem', 'ipsum']
    # 1.0



if __name__ == "__main__":
    unittest.main()
