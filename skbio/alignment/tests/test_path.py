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
        with self.assertRaises(ValueError):
            PairAlignPath.from_cigar("")

        # test invalid cigar string
        # Do we want from_cigar to handle other characters besides
        # 'M', 'I', 'D', '=', and 'X'?
        with self.assertRaises(ValueError):
            PairAlignPath.from_cigar("23M45B13X")

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

    def from_bits(self):
        # test base case
        bits = np.array(([0, 1, 0, 0, 0, 1], [1, 0, 0, 0, 1, 0]))
        exp = PairAlignPath(lengths=[1, 1, 2, 1, 1], states=[2, 1, 0, 2, 1], n_seqs=2)
        obs = PairAlignPath.from_bits(bits)
        self.assertEqual(obs, exp)

        # test empty bit array
        bits = np.array(([],[]))


if __name__ == "__main__":
    unittest.main()
