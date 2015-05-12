# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

from skbio import Protein


class TestProtein(unittest.TestCase):
    def test_nondegenerate_chars(self):
        exp = set("ACDEFGHIKLMNPQRSTVWY*")
        self.assertEqual(Protein("").nondegenerate_chars, exp)
        self.assertEqual(Protein.nondegenerate_chars, exp)

    def test_degenerate_map(self):
        exp = {
            'B': set(['D', 'N']), 'Z': set(['E', 'Q']),
            'X': set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                      'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
        }
        self.assertEqual(Protein("").degenerate_map, exp)
        self.assertEqual(Protein.degenerate_map, exp)

    def test_motif_n_glycosylation(self):
        seq = Protein("ACDFFACGNPSL")
        self.assertEqual(list(seq.find_motifs("N-glycosylation")), [])

        seq = Protein("ACDFNFTACGNPSL")
        self.assertEqual(list(seq.find_motifs("N-glycosylation")),
                         [slice(4, 8)])

        seq = Protein("AC-DFN-FTACGNPSL")
        self.assertEqual(list(seq.find_motifs("N-glycosylation",
                                              ignore=seq.gaps())),
                         [slice(5, 10)])


if __name__ == "__main__":
    unittest.main()
