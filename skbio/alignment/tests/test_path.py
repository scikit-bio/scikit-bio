# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from skbio.alignment._path import PairAlignPath, AlignPath

class TestPairAlignPath(unittest.TestCase):
    def test_from_cigar(self):
        # test valid cigar with no = or X
        cigar = "3M42I270M32D"
        obs = PairAlignPath.from_cigar(cigar)
        exp = ([3, 42, 270, 32], [0, 1, 0, 2])
        self.assertEqual(obs, exp)

        # test valid cigar with = or X
        cigar = cigar + "24X43="
        obs = PairAlignPath.from_cigar(cigar)
        exp = ([3, 42, 270, 32, 24, 43], [0, 1, 0, 2, 0, 0])
        self.assertEqual(obs, exp)

        # test empty cigar string

        # test invalid cigar string


if __name__ == "__main__":
    unittest.main()
