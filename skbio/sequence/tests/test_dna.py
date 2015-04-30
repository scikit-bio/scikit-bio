# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function
from future.standard_library import hooks

from re import compile as re_compile
from collections import Counter, defaultdict, Hashable
from unittest import TestCase, main
from itertools import product, chain

import numpy as np
import numpy.testing as npt
from scipy.spatial.distance import euclidean

from skbio import DNA
from skbio.sequence import SequenceError
from skbio.util._testing import IDValidationTests

with hooks():
    from itertools import zip_longest


class DNATests(TestCase):

    def setUp(self):
        self.b1 = DNA('GATTACA')
        self.b2 = DNA('ACCGGTACC', id="test-seq-2",
                              description="A test sequence", quality=range(9))
        self.b4 = DNA(
            'MRWSYKVHDBN', id="degen",
            description="All of the degenerate bases")
        self.b5 = DNA('.G--ATTAC-A...')

    def test_alphabet(self):
        exp = {
            'A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'T', 'W',
            'V', 'Y', '-', '.'
        }

        self.assertEqual(self.b1.alphabet, exp)
        self.assertEqual(DNA.alphabet, exp)

    def test_gap_chars(self):
        self.assertEqual(self.b1.gap_chars, set('-.'))

    def test_complement_map(self):
        exp = {
            '-': '-', '.': '.', 'A': 'T', 'C': 'G', 'B': 'V', 'D': 'H',
            'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', 'S': 'S',
            'R': 'Y', 'T': 'A', 'W': 'W', 'V': 'B', 'Y': 'R'
        }
        self.assertEqual(self.b1.complement_map, exp)
        self.assertEqual(DNA.complement_map, exp)

    def test_nondegenerate_chars(self):
        exp = set("ACGT")
        self.assertEqual(self.b1.nondegenerate_chars, exp)
        self.assertEqual(DNA.nondegenerate_chars, exp)

    def test_degenerate_map(self):
        exp = {
            'B': set(['C', 'T', 'G']), 'D': set(['A', 'T', 'G']),
            'H': set(['A', 'C', 'T']), 'K': set(['T', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'T', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']), 'W': set(['A', 'T']),
            'V': set(['A', 'C', 'G']), 'Y': set(['C', 'T'])
        }
        self.assertEqual(self.b1.degenerate_map, exp)
        self.assertEqual(DNA.degenerate_map, exp)

    def test_degenerate_chars(self):
        exp = set(['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y'])
        self.assertEqual(self.b1.degenerate_chars, exp)
        self.assertEqual(DNA.degenerate_chars, exp)

    def test_complement(self):
        # use equals method to ensure that id, description, and quality are
        # correctly propagated to the resulting sequence
        self.assertTrue(self.b1.complement().equals(DNA("CTAATGT")))

        self.assertTrue(self.b2.complement().equals(
            DNA("TGGCCATGG", id="test-seq-2",
                        description="A test sequence", quality=range(9))))

        self.assertTrue(self.b4.complement().equals(
            DNA("KYWSRMBDHVN", id="degen",
                        description="All of the degenerate bases")))

        self.assertTrue(self.b5.complement().equals(
            DNA(".C--TAATG-T...")))

    def test_reverse_complement(self):
        # use equals method to ensure that id, description, and (reversed)
        # quality scores are correctly propagated to the resulting sequence
        self.assertTrue(self.b1.reverse_complement().equals(
            DNA("TGTAATC")))

        self.assertTrue(self.b2.reverse_complement().equals(
            DNA("GGTACCGGT", id="test-seq-2",
                        description="A test sequence",
                        quality=range(9)[::-1])))

        self.assertTrue(self.b4.reverse_complement().equals(
            DNA("NVHDBMRSWYK", id="degen",
                        description="All of the degenerate bases")))

    def test_is_reverse_complement(self):
        self.assertFalse(self.b1.is_reverse_complement(self.b1))

        # id, description, and quality scores should be ignored (only sequence
        # data and type should be compared)
        self.assertTrue(self.b1.is_reverse_complement(
            DNA('TGTAATC', quality=range(7))))

        self.assertTrue(
            self.b4.is_reverse_complement(DNA('NVHDBMRSWYK')))



    def test_nondegenerates_no_degens(self):
        self.assertEqual(list(self.b1.nondegenerates()), [self.b1])

    def test_nondegenerates_all_degens(self):
        # Same chars.
        exp = [DNA('CC'), DNA('CG'), DNA('GC'),
               DNA('GG')]
        # Sort based on sequence string, as order is not guaranteed.
        obs = sorted(DNA('SS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Different chars.
        exp = [DNA('AC'), DNA('AG'), DNA('GC'),
               DNA('GG')]
        obs = sorted(DNA('RS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Odd number of chars.
        obs = list(DNA('NNN').nondegenerates())
        self.assertEqual(len(obs), 4**3)

    def test_nondegenerates_mixed_degens(self):
        exp = [DNA('AGC'), DNA('AGT'), DNA('GGC'),
               DNA('GGT')]
        obs = sorted(DNA('RGY').nondegenerates(), key=str)
        self.assertEqual(obs, exp)


if __name__ == "__main__":
    main()
