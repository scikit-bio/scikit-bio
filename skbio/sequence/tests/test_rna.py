# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

from skbio import RNA


class TestRNA(unittest.TestCase):
    def test_nondegenerate_chars(self):
        exp = set("ACGU")
        self.assertEqual(RNA('').nondegenerate_chars, exp)
        self.assertEqual(RNA.nondegenerate_chars, exp)

    def test_degenerate_map(self):
        exp = {
            'B': set(['C', 'U', 'G']), 'D': set(['A', 'U', 'G']),
            'H': set(['A', 'C', 'U']), 'K': set(['U', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'U', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']), 'W': set(['A', 'U']),
            'V': set(['A', 'C', 'G']), 'Y': set(['C', 'U'])
        }
        self.assertEqual(RNA('').degenerate_map, exp)
        self.assertEqual(RNA.degenerate_map, exp)

    def test_complement_map(self):
        exp = {
            '-': '-', '.': '.', 'A': 'U', 'C': 'G', 'B': 'V', 'D': 'H',
            'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', 'S': 'S',
            'R': 'Y', 'U': 'A', 'W': 'W', 'V': 'B', 'Y': 'R'
        }
        self.assertEqual(RNA('').complement_map, exp)
        self.assertEqual(RNA.complement_map, exp)

    def test_motif_purine_run(self):
        seq = RNA("")
        self.assertEqual(list(seq.find_motifs("purine-run")), [])

        seq = RNA("AARC--UCRG")
        self.assertEqual(list(seq.find_motifs("purine-run")),
                         [slice(0, 3), slice(8, 10)])

        seq = RNA("AA-RC--UCR-G")
        self.assertEqual(list(seq.find_motifs("purine-run", min_length=3,
                                              exclude=seq.gaps())),
                         [slice(0, 4)])

    def test_motif_pyrimidine_run(self):
        seq = RNA("")
        self.assertEqual(list(seq.find_motifs("pyrimidine-run")), [])

        seq = RNA("AARC--UCRG")
        self.assertEqual(list(seq.find_motifs("pyrimidine-run")),
                         [slice(3, 4), slice(6, 8)])

        seq = RNA("AA-RC--UCR-G")
        self.assertEqual(list(seq.find_motifs("pyrimidine-run", min_length=3,
                                              exclude=seq.gaps())),
                         [slice(4, 9)])

if __name__ == "__main__":
    unittest.main()
