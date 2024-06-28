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

from skbio.alignment._path import PairAlignPath, AlignPath, _run_length_encode
from skbio.alignment import TabularMSA
from skbio.sequence import DNA


class TestAlignPath(unittest.TestCase):
    def test_init(self):
        # test 1-D starts vector
        with self.assertRaises(TypeError, msg="`starts` must be a 1-D vector."):
            path = AlignPath(lengths=[1, 2, 3], states=[1, 2, 3], starts=[[0], [0]])

        # test states and starts matching
        with self.assertRaises(ValueError, msg="Sizes of `starts` and `states` do not "
                               "match."):
            path = AlignPath(lengths=[1, 2, 3], states=[1, 2, 3],
                             starts=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    def test_lengths(self):
        obs = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                        states=[0, 2, 0, 6, 0, 1, 0],
                        starts=[0, 0, 0]).lengths
        npt.assert_array_equal(obs, np.array([3, 2, 5, 1, 4, 3, 2], dtype=np.int64))

    def test_states(self):
        obs = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                        states=[0, 2, 0, 6, 0, 1, 0],
                        starts=[0, 0, 0]).states
        npt.assert_array_equal(obs, np.array([[0, 2, 0, 6, 0, 1, 0]], dtype=np.uint8))

    def test_starts(self):
        obs = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                        states=[0, 2, 0, 6, 0, 1, 0],
                        starts=[0, 0, 0]).starts
        npt.assert_array_equal(obs, np.array([0, 0, 0], dtype=np.int64))

    def test_shape(self):
        obs = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                        states=[0, 2, 0, 6, 0, 1, 0],
                        starts=[0, 0, 0]).shape
        self.assertEqual(obs.sequence, 3)
        self.assertEqual(obs.position, 20)

    def test_to_bits(self):
        path = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                         states=[0, 2, 0, 6, 0, 1, 0],
                         starts=[0, 0, 0])
        exp = np.array(([0, 0, 0, 0, 0, 1, 0],
                        [0, 1, 0, 1, 0, 0, 0],
                        [0, 0, 0, 1, 0, 0, 0]))
        obs = path.to_bits()
        npt.assert_array_equal(obs, exp)

    def test_from_bits(self):
        # test 1D base case, less than 8 sequences
        bits = np.array(([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0],
                         [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
        exp = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                        states=[0, 2, 0, 6, 0, 1, 0],
                        starts=[0, 0, 0])
        obs = AlignPath.from_bits(bits)
        npt.assert_array_equal(obs.lengths, exp.lengths)
        npt.assert_array_equal(obs.states, exp.states)

        # test starts parameter
        starts = [1, 2, 3]
        obs = AlignPath.from_bits(bits, starts)
        npt.assert_array_equal(obs.lengths, exp.lengths)
        npt.assert_array_equal(obs.states, exp.states)
        npt.assert_array_equal(obs.starts, starts)

        # test 2D base case, more than 8 sequences
        rng = np.random.default_rng(seed=42)
        bits = rng.choice([0, 1], size=(10, 10), p=[0.85, 0.15])
        exp = AlignPath(lengths=[1, 1, 1, 1, 1, 1, 3, 1],
                        states=[[0, 10, 133, 4, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 2]],
                        starts=[0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        obs = AlignPath.from_bits(bits)
        npt.assert_array_equal(obs.lengths, exp.lengths)
        npt.assert_array_equal(obs.states, exp.states)

    def test_from_tabular(self):
        msa = ('CGGTCGTAACGCGTA---CA',
               'CAG--GTAAG-CATACCTCA',
               'CGGTCGTCAC-TGTACACTA')
        tabular = TabularMSA([DNA(x) for x in msa])
        path = AlignPath.from_tabular(tabular)
        lengths = [3, 2, 5, 1, 4, 3, 2]
        states = [0, 2, 0, 6, 0, 1, 0]
        npt.assert_array_equal(lengths, path.lengths)
        npt.assert_array_equal(states, np.squeeze(path.states))

    def test_to_indices(self):
        # test gap = -1
        path = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                         states=[0, 2, 0, 6, 0, 1, 0],
                         starts=[0, 0, 0])
        exp = np.array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, -1, -1, -1,
                         15, 16],
                        [0, 1, 2, -1, -1, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, 12, 13, 14,
                         15, 16],
                        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, 10, 11, 12, 13, 14, 15, 16,
                         17, 18]])
        obs = path.to_indices()
        npt.assert_array_equal(obs, exp)

        # test gap = 'del'
        exp = np.array([[0, 1, 2, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15, 16],
                        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 16],
                        [0, 1, 2, 5, 6, 7, 8, 9, 10, 11, 12, 13, 17, 18]])
        obs = path.to_indices(gap='del')
        npt.assert_array_equal(obs, exp)

        # test gap = 'mask'
        exp = np.ma.array(data=[[0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
                                 13, 14, 14, 14, 14, 15, 16],
                                [0,  1,  2,  2,  2,  3,  4,  5,  6,  7,  7,  8,  9,
                                 10, 11, 12, 13, 14, 15, 16],
                                [0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  9, 10, 11,
                                 12, 13, 14, 15, 16, 17, 18]],
                          mask=[[False, False, False, False, False, False, False,
                                 False, False, False, False, False, False, False,
                                 False,  True,  True,  True, False, False],
                                [False, False, False,  True,  True, False, False,
                                 False, False, False,  True, False, False, False,
                                 False, False, False, False, False, False],
                                [False, False, False, False, False, False, False,
                                 False, False, False,  True, False, False, False,
                                 False, False, False, False, False, False]],
                          fill_value=999999)
        obs = path.to_indices(gap='mask')
        npt.assert_array_equal(obs, exp)

        # test with starts as non-zero
        path = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                         states=[0, 2, 0, 6, 0, 1, 0],
                         starts=[1, 35, 28])
        exp = np.array([[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, -1, -1, -1,
                         16, 17],
                        [35, 36, 37, -1, -1, 38, 39, 40, 41, 42, -1, 43, 44, 45, 46,
                         47, 48, 49, 50, 51],
                        [28, 29, 30, 31, 32, 33, 34, 35, 36, 37, -1, 38, 39, 40, 41,
                         42, 43, 44, 45, 46]])
        obs = path.to_indices()
        npt.assert_array_equal(obs, exp)

        # test 'del' with non-zero starts
        exp = np.array([[1,  2,  3,  6,  7,  8,  9, 10, 12, 13, 14, 15, 16, 17],
                        [35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 50, 51],
                        [28, 29, 30, 33, 34, 35, 36, 37, 38, 39, 40, 41, 45, 46]])
        obs = path.to_indices(gap='del')
        npt.assert_array_equal(obs, exp)

        # test 'mask' with non-zero starts
        exp = np.ma.array(data=[[1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,
                                 14, 15, 15, 15, 15, 16, 17],
                                [35, 36, 37, 37, 37, 38, 39, 40, 41, 42, 42, 43, 44,
                                 45, 46, 47, 48, 49, 50, 51],
                                [28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 37, 38, 39,
                                 40, 41, 42, 43, 44, 45, 46]],
                          mask=[[False, False, False, False, False, False, False,
                                 False, False, False, False, False, False, False,
                                 False,  True,  True,  True, False, False],
                                [False, False, False,  True,  True, False, False,
                                 False, False, False,  True, False, False, False,
                                 False, False, False, False, False, False],
                                [False, False, False, False, False, False, False,
                                 False, False, False,  True, False, False, False,
                                 False, False, False, False, False, False]],
                          fill_value=999999)
        obs = path.to_indices(gap='mask')
        npt.assert_array_equal(obs, exp)

        # test invalid gap
        with self.assertRaises(TypeError,
                               msg="Gap must be an integer, np.nan, np.inf, 'del', "
                               "or 'mask'."):
            path.to_indices(gap="no")

    def test_from_indices(self):
        # test no mask
        indices = np.array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, -1, -1,
                             -1, 15, 16],
                            [0, 1, 2, -1, -1, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, 12, 13,
                             14, 15, 16],
                            [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, -1, 10, 11, 12, 13, 14, 15,
                             16, 17, 18]])
        path = AlignPath.from_indices(indices)
        lengths = [3, 2, 5, 1, 4, 3, 2]
        states = [0, 2, 0, 6, 0, 1, 0]
        npt.assert_array_equal(lengths, path.lengths)
        npt.assert_array_equal(states, np.squeeze(path.states))

        # test masked array
        masked = np.ma.array(indices, mask=(indices == -1))
        path = AlignPath.from_indices(masked, gap="mask")
        npt.assert_array_equal(lengths, path.lengths)
        npt.assert_array_equal(states, np.squeeze(path.states))

        # test non-zero indices
        indices = np.array([[25, 26, -1, -1, 27, 28, 29, 30],
                            [-1, 79, 80, 81, 82, 83, 84, -1]])
        path = AlignPath.from_indices(indices)
        lengths = [1, 1, 2, 3, 1]
        states = [2, 0, 1, 0, 2]
        npt.assert_array_equal(lengths, path.lengths)
        npt.assert_array_equal(states, np.squeeze(path.states))
        npt.assert_array_equal(path.starts, [25, 79])

        # test masked array and non-zero indices
        # TODO

        # test indices all gaps
        indices = np.array([[25, 26, -1, -1, 27, 28, 29, 30],
                            [-1, -1, -1, -1, -1, -1, -1, -1]])
        path = AlignPath.from_indices(indices)
        lengths = [2, 2, 4]
        states = [2, 3, 2]
        starts = [25, -1]
        npt.assert_array_equal(path.lengths, lengths)
        npt.assert_array_equal(np.squeeze(path.states), states)
        npt.assert_array_equal(path.starts, starts)

    def test_to_coordinates(self):
        # test base case
        exp = np.array([[0, 3, 5, 10, 11, 15, 15, 17],
                        [0, 3, 3,  8,  8, 12, 15, 17],
                        [0, 3, 5, 10, 10, 14, 17, 19]])
        path = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                         states=[0, 2, 0, 6, 0, 1, 0],
                         starts=[0, 0, 0])
        obs = path.to_coordinates()
        npt.assert_array_equal(obs, exp)

        # test non-zero starts
        exp = np.array([[2,  5,  7, 12, 13, 17, 17, 19],
                        [512, 515, 515, 520, 520, 524, 527, 529],
                        [28, 31, 33, 38, 38, 42, 45, 47]])
        path = AlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                         states=[0, 2, 0, 6, 0, 1, 0],
                         starts=[2, 512, 28])
        obs = path.to_coordinates()
        npt.assert_array_equal(obs, exp)

    def test_from_coordinates(self):
        # test base case
        coords = np.array([[0, 3, 5, 10, 11, 15, 15, 17],
                           [0, 3, 3,  8,  8, 12, 15, 17],
                           [0, 3, 5, 10, 10, 14, 17, 19]])
        path = AlignPath.from_coordinates(coords)
        lengths = [3, 2, 5, 1, 4, 3, 2]
        states = [0, 2, 0, 6, 0, 1, 0]
        starts = [0, 0, 0]
        npt.assert_array_equal(lengths, path.lengths)
        npt.assert_array_equal(states, np.squeeze(path.states))
        npt.assert_array_equal(starts, path.starts)

        # test non-zero starts
        coords = np.array([[2,  5,  7, 12, 13, 17, 17, 19],
                           [512, 515, 515, 520, 520, 524, 527, 529],
                           [28, 31, 33, 38, 38, 42, 45, 47]])
        path = AlignPath.from_coordinates(coords)
        lengths = [3, 2, 5, 1, 4, 3, 2]
        states = [0, 2, 0, 6, 0, 1, 0]
        starts = [2, 512, 28]
        npt.assert_array_equal(lengths, path.lengths)
        npt.assert_array_equal(states, np.squeeze(path.states))
        npt.assert_array_equal(starts, path.starts)


class TestPairAlignPath(unittest.TestCase):
    def test_from_cigar(self):
        # test valid cigar with no = or X
        cigar = "3M42I270M32D"
        path = PairAlignPath.from_cigar(cigar)
        lengths = [3, 42, 270, 32]
        states = [0, 1, 0, 2]
        npt.assert_array_equal(lengths, path.lengths)
        npt.assert_array_equal(states, np.squeeze(path.states))

        # test valid cigar with = or X
        cigar = "3M42I270M23X663=32D24X43="
        path = PairAlignPath.from_cigar(cigar)
        lengths = [3, 42, 956, 32, 67]
        states = [0, 1, 0, 2, 0]
        npt.assert_array_equal(lengths, path.lengths)
        npt.assert_array_equal(states, np.squeeze(path.states))

        # test empty cigar string
        with self.assertRaises(ValueError, msg="CIGAR string must not be empty."):
            PairAlignPath.from_cigar("")

        # test invalid cigar string
        with self.assertRaises(ValueError, msg="Invalid characters in CIGAR string."):
            PairAlignPath.from_cigar("23M45B13X")

        # test valid cigar with no 1's
        cigar = "MID12MI"
        path = PairAlignPath.from_cigar(cigar)
        lengths = [1, 1, 1, 12, 1]
        states = [0, 1, 2, 0, 1]
        npt.assert_array_equal(lengths, path.lengths)
        npt.assert_array_equal(states, np.squeeze(path.states))

        # test cigar with all possible valid codes
        cigar = "1M2I3D4P5=6X7N8S9H"
        path = PairAlignPath.from_cigar(cigar)
        lengths = [1, 2, 3, 4, 11, 7, 8, 9]
        states = [0, 1, 2, 3, 0, 2, 1, 3]
        npt.assert_array_equal(lengths, path.lengths)
        npt.assert_array_equal(states, np.squeeze(path.states))

    def test_to_cigar(self):
        # test base case
        lengths = [1, 2, 3, 4, 50, 234]
        gaps = [1, 0, 2, 1, 0, 1]
        ExamplePairAlignPath = PairAlignPath(lengths=lengths,
                                             states=gaps,
                                             starts=[0, 0])
        obs = ExamplePairAlignPath.to_cigar()
        exp = "1I2M3D4I50M234I"
        self.assertEqual(obs, exp)

        # test if seqs are provided
        seq1 = "-AATCT----" + "C"*50 + "-"*234
        seq2 = "TAC---GGCC" + "C"*20 + "A"*264
        seqs = [DNA(seq1), DNA(seq2)]
        obs = ExamplePairAlignPath.to_cigar(seqs=seqs)
        exp = "1I2X3D4I5X18=27X234I"
        self.assertEqual(obs, exp)

        # test if alignment has two gaps in same position
        lengths = [1, 2, 3, 4, 1]
        gaps = [1, 0, 2, 1, 3]
        path = PairAlignPath(lengths=lengths, states=gaps, starts=[0, 0])
        obs = path.to_cigar()
        exp = "1I2M3D4I1P"
        self.assertEqual(obs, exp)

        # two gaps with seqs provided
        seq1 = '-ATCGC-----'
        seq2 = 'GTA---ATTA-'
        seqs = [DNA(seq1), DNA(seq2)]
        obs = path.to_cigar(seqs=seqs)
        exp = "1I1X1=3D4I1P"
        self.assertEqual(obs, exp)

        # test sequences as strings
        seq1 = '-ATCGC-----'
        seq2 = 'GTA---ATTA-'
        seqs = [seq1, seq2]
        obs = path.to_cigar(seqs=seqs)
        exp = "1I1X1=3D4I1P"
        self.assertEqual(obs, exp)

        # test invalid sequence input
        seq1 = 1
        seq2 = 'GTA---ATTA-'
        seqs = [seq1, seq2]
        with self.assertRaises(TypeError,
                               msg="`seqs` must be of type string or Sequence object."):
            obs = path.to_cigar(seqs=seqs)

    def test_from_bits(self):
        # test base case
        bits = np.array(([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0],
                         [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
        exp = PairAlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                            states=[0, 2, 0, 2, 0, 1, 0],
                            starts=[0, 0])
        obs = PairAlignPath.from_bits(bits)
        npt.assert_array_equal(obs.lengths, exp.lengths)
        npt.assert_array_equal(obs.states, exp.states)

        # test empty bit array
        bits = np.array(([], []))
        with self.assertRaises(TypeError,
                               msg="Input 'bits' must be a non-empty 2D numpy array."):
            PairAlignPath.from_bits(bits)

        # test 1D bit array
        bits = np.array([0, 0, 1])
        with self.assertRaises(TypeError,
                               msg="Input 'bits' must be a non-empty 2D numpy array."):
            PairAlignPath.from_bits(bits)

        # test array with invalid values
        bits = np.array(([1, 2, 3], [0, 5, 1]))
        with self.assertRaises(ValueError,
                               msg="Input 'bits' must contain only zeros and ones."):
            PairAlignPath.from_bits(bits)

        # test non numpy array input
        bits = [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0],
                [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        exp = PairAlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                            states=[0, 2, 0, 2, 0, 1, 0],
                            starts=[0, 0, 0])
        obs = PairAlignPath.from_bits(bits)
        npt.assert_array_equal(obs.lengths, exp.lengths)
        npt.assert_array_equal(obs.states, exp.states)

    def test_to_bits(self):
        # test input with invalid values
        with self.assertRaises(ValueError,
                               msg="For pairwise alignment, `states` must only "
                               "contain zeros, ones, twos, or threes."):
            PairAlignPath(lengths=[1, 2, 3], states=[1, 2, 4], starts=[0, 0]).to_bits()

        # test base case
        path = PairAlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                             states=[0, 2, 0, 2, 0, 1, 0],
                             starts=[0, 0])
        exp = np.array(([0, 0, 0, 0, 0, 1, 0], [0, 1, 0, 1, 0, 0, 0]))
        obs = path.to_bits()
        npt.assert_array_equal(np.squeeze(obs), exp)

    def test_run_length_encode(self):
        obs = _run_length_encode("ABBCCCDDDD")
        exp = "1A2B3C4D"
        self.assertEqual(obs, exp)


if __name__ == "__main__":
    unittest.main()
