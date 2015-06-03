# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

from skbio import DNA


class TestDNA(unittest.TestCase):
    def test_nondegenerate_chars(self):
        exp = set("ACGT")
        self.assertEqual(DNA('').nondegenerate_chars, exp)
        self.assertEqual(DNA.nondegenerate_chars, exp)

    def test_degenerate_map(self):
        exp = {
            'B': set(['C', 'T', 'G']), 'D': set(['A', 'T', 'G']),
            'H': set(['A', 'C', 'T']), 'K': set(['T', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'T', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']), 'W': set(['A', 'T']),
            'V': set(['A', 'C', 'G']), 'Y': set(['C', 'T'])
        }
        self.assertEqual(DNA('').degenerate_map, exp)
        self.assertEqual(DNA.degenerate_map, exp)

    def test_complement_map(self):
        exp = {
            '-': '-', '.': '.', 'A': 'T', 'C': 'G', 'B': 'V', 'D': 'H',
            'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', 'S': 'S',
            'R': 'Y', 'T': 'A', 'W': 'W', 'V': 'B', 'Y': 'R'
        }
        self.assertEqual(DNA('').complement_map, exp)
        self.assertEqual(DNA.complement_map, exp)

    def test_motif_purine_run(self):
        seq = DNA("")
        self.assertEqual(list(seq.find_motifs("purine-run")), [])

        seq = DNA("AARC--TCRG")
        self.assertEqual(list(seq.find_motifs("purine-run")),
                         [slice(0, 3), slice(8, 10)])

        seq = DNA("AA-RC--TCR-G")
        self.assertEqual(list(seq.find_motifs("purine-run", min_length=3,
                                              ignore=seq.gaps())),
                         [slice(0, 4)])

    def test_motif_pyrimidine_run(self):
        seq = DNA("")
        self.assertEqual(list(seq.find_motifs("pyrimidine-run")), [])

        seq = DNA("AARC--TCRA")
        self.assertEqual(list(seq.find_motifs("pyrimidine-run")),
                         [slice(3, 4), slice(6, 8)])

        seq = DNA("AA-RC--TCR-A")
        self.assertEqual(list(seq.find_motifs("pyrimidine-run", min_length=3,
                                              ignore=seq.gaps())),
                         [slice(4, 9)])


if __name__ == "__main__":
    unittest.main()
