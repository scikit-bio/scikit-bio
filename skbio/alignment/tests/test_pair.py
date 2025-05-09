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

from skbio.alignment import PairAlignPath
from skbio.alignment import TabularMSA
from skbio.sequence import Sequence, DNA, Protein

from skbio.alignment._pair import PairAligner

from skbio.alignment._pair import (
    _alloc_matrices,
    _init_matrices,
)



class TestPairAlign(unittest.TestCase):
    def setUp(self):
        pass

    def test_alloc_matrices(self):
        # linear
        obs = _alloc_matrices(3, 4, affine=False)
        self.assertTrue(isinstance(obs[0], np.ndarray))
        self.assertTupleEqual(obs[0].shape, (4, 5))
        self.assertEqual(obs[0].dtype, np.float64)
        self.assertIsNone(obs[1])
        self.assertIsNone(obs[2])

        # ensure matrix is C-contiguous
        self.assertTrue(obs[0].flags.c_contiguous)

        # affine
        obs = _alloc_matrices(3, 4, affine=True)
        for mat in obs:
            self.assertTupleEqual(mat.shape, (4, 5))
            self.assertEqual(mat.dtype, np.float64)

        # float32
        obs = _alloc_matrices(3, 4, affine=False, dtype=np.float32)
        self.assertEqual(obs[0].dtype, np.float32)

    def test_init_matrices(self):
        # classical global alignment
        obs = _alloc_matrices(3, 4, affine=False)
        _init_matrices(*obs, 0, 2, local=False, free_ends=False)
        npt.assert_array_equal(obs[0][:, 0], [0, -2, -4, -6])
        npt.assert_array_equal(obs[0][0, :], [0, -2, -4, -6, -8])

        # local alignment
        obs = _alloc_matrices(3, 4, affine=False)
        _init_matrices(*obs, 0, 2, local=True, free_ends=False)
        self.assertTrue((obs[0][:, 0] == 0).all())
        self.assertTrue((obs[0][0, :] == 0).all())

        # semi-global alignment (free ends)
        obs = _alloc_matrices(3, 4, affine=False)
        _init_matrices(*obs, 0, 2, local=False, free_ends=True)
        self.assertTrue((obs[0][:, 0] == 0).all())
        self.assertTrue((obs[0][0, :] == 0).all())

        # affine and global
        obs = _alloc_matrices(3, 4, affine=True)
        _init_matrices(*obs, 5, 2, local=False, free_ends=False)
        npt.assert_array_equal(obs[0][:, 0], [0, -7, -9, -11])
        npt.assert_array_equal(obs[0][0, :], [0, -7, -9, -11, -13])
        self.assertTrue(np.isneginf(obs[1][1:, 0]).all())
        self.assertTrue(np.isneginf(obs[2][0, 1:]).all())

        # affine and local
        obs = _alloc_matrices(3, 4, affine=True)
        _init_matrices(*obs, 5, 2, local=True, free_ends=False)
        self.assertTrue((obs[0][:, 0] == 0).all())
        self.assertTrue((obs[0][0, :] == 0).all())
        self.assertTrue((obs[1][1:, 0] == 0).all())
        self.assertTrue((obs[2][0, 1:] == 0).all())


class TestPairAligner(unittest.TestCase):
    def setUp(self):
        pass

    def test_init(self):
        # default parameters
        obs = PairAligner()
        self.assertFalse(obs._local)
        npt.assert_array_equal(obs._submat, np.empty((0, 0)))
        self.assertIsNone(obs._charmap)
        npt.assert_array_equal(obs._mthmis, np.array([1, -1]))
        self.assertFalse(obs._affine)
        self.assertEqual(obs._gap_open, 0)
        self.assertEqual(obs._gap_extend, -2)
        self.assertFalse(obs._terminal)

        # custom parameters
        obs = PairAligner(
            mode="local",
            sub_score="BLOSUM62",
            gap_score=(-11, -1),
            terminal_gaps=True,
        )
        self.assertTrue(obs._local)
        npt.assert_array_equal(obs._submat[0, :5], [4, -1, -2, -2, 0])
        npt.assert_array_equal(obs._charmap[65:70], [0, 20, 4, 3, 6])
        npt.assert_array_equal(obs._mthmis, [])
        self.assertTrue(obs._affine)
        self.assertEqual(obs._gap_open, -11)
        self.assertEqual(obs._gap_extend, -1)
        self.assertTrue(obs._terminal)

    def test_mode(self):
        self.assertEqual(PairAligner().mode, "global")
        self.assertEqual(PairAligner(mode="local").mode, "local")

    def test_sub_score(self):
        self.assertTupleEqual(PairAligner().sub_score, (1, -1))

    def test_gap_score(self):
        self.assertEqual(PairAligner().gap_score, -2)
        self.assertTupleEqual(PairAligner(gap_score=(-5, -2)).gap_score, (-5, -2))

    def test_terminal(self):
        self.assertFalse(PairAligner().terminal)
        self.assertTrue(PairAligner(terminal_gaps=True).terminal)

    def test_align(self):
        obs = PairAligner()
        obs.align("AATCG", "AACG")
        npt.assert_array_equal(obs._seq1, [65, 65, 84, 67, 71])
        npt.assert_array_equal(obs._seq2, [65, 65, 67, 71])
        self.assertEqual(obs._len1, 5)
        self.assertEqual(obs._len2, 4)
        self.assertEqual(obs._score, 2)
        exp = np.array([
            [  0,  -2,  -4,  -6,  -8],
            [ -2,   1,  -1,  -3,  -5],
            [ -4,  -1,   2,   0,  -2],
            [ -6,  -3,   0,   1,  -1],
            [ -8,  -5,  -2,   1,   0],
            [-10,  -7,  -4,  -1,   2],
        ])
        npt.assert_array_equal(obs._alnmat, exp)

        obs = PairAligner()
        obs.align("AAGCGTC", "AGCGGTA")
        self.assertEqual(obs._score, 0)
        exp = np.array([
            [  0,  -2,  -4,  -6,  -8, -10, -12, -14],
            [ -2,   1,  -1,  -3,  -5,  -7,  -9, -11],
            [ -4,  -1,   0,  -2,  -4,  -6,  -8,  -8],
            [ -6,  -3,   0,  -1,  -1,  -3,  -5,  -7],
            [ -8,  -5,  -2,   1,  -1,  -2,  -4,  -6],
            [-10,  -7,  -4,  -1,   2,   0,  -2,  -4],
            [-12,  -9,  -6,  -3,   0,   1,   1,  -1],
            [-14, -11,  -8,  -5,  -2,  -1,   0,   0],
        ])
        npt.assert_array_equal(obs._alnmat, exp)

    def test_alloc_linear_matrix(self):
        obs = PairAligner()
        obs._len1 = 5
        obs._len2 = 4
        obs._alloc_linear_matrix()
        self.assertTupleEqual(obs._alnmat.shape, (6, 5))
        self.assertTrue(obs._alnmat.flags.c_contiguous)
        self.assertIsNone(obs._insmat)
        self.assertIsNone(obs._delmat)

    def test_alloc_affine_matrix(self):
        obs = PairAligner()
        obs._len1 = 5
        obs._len2 = 4
        obs._alloc_affine_matrix()
        for key in ("aln", "ins", "del"):
            mat = getattr(obs, f"_{key}mat")
            self.assertTupleEqual(mat.shape, (6, 5))
            self.assertTrue(mat.flags.c_contiguous)

    def test_init_global_linear_matrix(self):
        obs = PairAligner()
        obs._len1 = 3
        obs._len2 = 4
        obs._alloc_linear_matrix()
        obs._gap_extend = -2
        obs._init_global_linear_matrix()
        npt.assert_array_equal(obs._alnmat[:, 0], [0, -2, -4, -6])
        npt.assert_array_equal(obs._alnmat[0, :], [0, -2, -4, -6, -8])

    def test_init_local_linear_matrix(self):
        obs = PairAligner()
        obs._len1 = 3
        obs._len2 = 4
        obs._alloc_linear_matrix()
        obs._gap_extend = -2
        obs._init_local_linear_matrix()
        npt.assert_array_equal(obs._alnmat[:, 0], [0, 0, 0, 0])
        npt.assert_array_equal(obs._alnmat[0, :], [0, 0, 0, 0, 0])

    def test_trace_global_linear_one(self):
        pass


if __name__ == "__main__":
    unittest.main()
