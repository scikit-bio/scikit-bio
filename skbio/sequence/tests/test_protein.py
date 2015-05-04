# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

from unittest import TestCase, main

from skbio import Protein


class ProteinTests(TestCase):

    def setUp(self):
        self.p1 = Protein('GREG')
        self.p2 = Protein(
            'PRTEINSEQNCE', id="test-seq-2",
            description="A test sequence")

    def test_alphabet(self):
        exp = {
            'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N',
            'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z', '-', '.', '*'
        }

        self.assertEqual(self.p1.alphabet, exp)
        self.assertEqual(Protein.alphabet, exp)

    def test_gap_chars(self):
        self.assertEqual(self.p1.gap_chars, set('-.'))

    def test_nondegenerate_chars(self):
        exp = set("ACDEFGHIKLMNPQRSTVWY*")
        self.assertEqual(self.p1.nondegenerate_chars, exp)
        self.assertEqual(Protein.nondegenerate_chars, exp)

    def test_degenerate_map(self):
        exp = {
            'B': set(['D', 'N']), 'Z': set(['E', 'Q']),
            'X': set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                      'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
        }
        self.assertEqual(self.p1.degenerate_map, exp)
        self.assertEqual(Protein.degenerate_map, exp)

    def test_degenerate_chars(self):
        exp = set(['B', 'X', 'Z'])
        self.assertEqual(self.p1.degenerate_chars, exp)
        self.assertEqual(Protein.degenerate_chars, exp)

    def test_nondegenerates(self):
        exp = [Protein('AD'), Protein('AN')]
        # Sort based on sequence string, as order is not guaranteed.
        obs = sorted(Protein('AB').nondegenerates(), key=str)
        self.assertEqual(obs, exp)


if __name__ == "__main__":
    main()
