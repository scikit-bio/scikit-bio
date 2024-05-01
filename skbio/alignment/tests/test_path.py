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

from skbio.alignment._path import PairAlignPath, AlignPath

class TestPairAlignPath(unittest.TestCase):
    def test_from_cigar(self):
        # test valid cigar with no = or X
        cigar = "3M42I270M32D"
        obs = PairAlignPath.from_cigar(cigar)
        exp = ([3, 42, 270, 32], [0, 1, 0, 2])
        npt.assert_array_equal(obs, exp)

        # test valid cigar with = or X
        cigar = "3M42I270M23X663=32D24X43="
        obs = PairAlignPath.from_cigar(cigar)
        exp = ([3, 42, 956, 32, 67], [0, 1, 0, 2, 0])
        npt.assert_array_equal(obs, exp)

        # test empty cigar string
        with self.assertRaises(ValueError, msg="CIGAR string must not be empty."):
            PairAlignPath.from_cigar("")

        # test invalid cigar string
        # Do we want from_cigar to handle other characters besides
        # 'M', 'I', 'D', '=', and 'X'?
        with self.assertRaises(ValueError, msg="Invalid characters in CIGAR string."):
            PairAlignPath.from_cigar("23M45B13X")

        # test valid cigar with no 1's
        cigar = "MID12MI"
        obs = PairAlignPath.from_cigar(cigar)
        exp = ([1, 1, 1, 12, 1], [0, 1, 2, 0, 1])
        npt.assert_array_equal(obs, exp)

    def test_to_cigar(self):
        # test base case
        lengths = [1, 2, 3, 4, 50, 234]
        gaps = [1, 0, 2, 1, 0, 1]
        ExamplePairAlignPath = PairAlignPath(lengths=lengths, states=gaps, n_seqs=2)
        obs = ExamplePairAlignPath.to_cigar()
        exp = "1I2M3D4I50M234I"
        # do we want to include or omit 1?
        # include 1 by default, maybe have option to not include 1
        self.assertEqual(obs, exp)

        # test if seqs are provided
        seq1 = "-AATCT----" + "C"*50 + "-"*234
        seq2 = "TAC---GGCC" + "C"*20 + "A"*264
        seqs = [seq1, seq2]
        obs = ExamplePairAlignPath.to_cigar(seqs=seqs)
        exp = "1I1=1X3D4I20=30X234I"
        self.assertEqual(obs, exp)

        # test invalid input

    def test_from_bits(self):
        # test base case
        bits = np.array(([0, 1, 0, 0, 0, 1], [1, 0, 0, 0, 1, 0]))
        exp = PairAlignPath(lengths=[1, 1, 2, 1, 1], states=[2, 1, 0, 2, 1], n_seqs=2)
        obs = PairAlignPath.from_bits(bits)
        npt.assert_array_equal(obs.lengths, exp.lengths)
        npt.assert_array_equal(obs.states, exp.states)

        # test empty bit array
        bits = np.array(([],[]))
        with self.assertRaises(ValueError,
                               msg="Input 'bits' must be a non-empty 2D numpy array."):
            PairAlignPath.from_bits(bits)

        # test 1D bit array
        bits = np.array([0, 0, 1])
        with self.assertRaises(ValueError,
                               msg="Input 'bits' must be a non-empty 2D numpy array."):
            PairAlignPath.from_bits(bits)

        # test array with invalid values
        bits = np.array(([1, 2, 3], [0, 5, 1]))
        with self.assertRaises(ValueError,
                               msg="Input 'bits' must contain only zeros and ones."):
            PairAlignPath.from_bits(bits)

    def test_to_bits(self):
        # test input with invalid values
        with self.assertRaises(ValueError,
                               msg="For pairwise alignment, `states` must only "
                               "contain zeros, ones, or twos."):
            PairAlignPath(lengths=[1, 2, 3], states=[1, 2, 3], n_seqs=2).to_bits()

        # test base case
        obj = PairAlignPath(lengths=[3, 2, 5, 1, 4, 3, 2],
                            states=[0, 2, 0, 2, 0, 1, 0],
                            n_seqs=2)
        exp = np.array(([0, 0, 0, 0, 0, 1, 0], [0, 1, 0, 1, 0, 0, 0]))
        obs = obj.to_bits()
        npt.assert_array_equal(obs, exp)

if __name__ == "__main__":
    unittest.main()
