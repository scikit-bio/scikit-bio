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

from skbio import RNA
from skbio.sequence import SequenceError
from skbio.util._testing import IDValidationTests

with hooks():
    from itertools import zip_longest


class RNATests(TestCase):

    def setUp(self):
        self.b1 = RNA('GAUUACA')
        self.b2 = RNA('ACCGGUACC', id="test-seq-2",
                              description="A test sequence", quality=range(9))
        self.b4 = RNA(
            'MRWSYKVHDBN', id="degen",
            description="All of the degenerate bases")
        self.b5 = RNA('.G--AUUAC-A...')

    def test_alphabet(self):
        exp = {
            'A', 'C', 'B', 'D', 'G', 'H', 'K', 'M', 'N', 'S', 'R', 'U', 'W',
            'V', 'Y', '-', '.'
        }

        self.assertEqual(self.b1.alphabet, exp)
        self.assertEqual(RNA.alphabet, exp)

    def test_gap_chars(self):
        self.assertEqual(self.b1.gap_chars, set('-.'))

    def test_complement_map(self):
        exp = {
            '-': '-', '.': '.', 'A': 'U', 'C': 'G', 'B': 'V', 'D': 'H',
            'G': 'C', 'H': 'D', 'K': 'M', 'M': 'K', 'N': 'N', 'S': 'S',
            'R': 'Y', 'U': 'A', 'W': 'W', 'V': 'B', 'Y': 'R'
        }
        self.assertEqual(self.b1.complement_map, exp)
        self.assertEqual(RNA.complement_map, exp)

    def test_nondegenerate_chars(self):
        exp = set("ACGU")
        self.assertEqual(self.b1.nondegenerate_chars, exp)
        self.assertEqual(RNA.nondegenerate_chars, exp)

    def test_degenerate_map(self):
        exp = {
            'B': set(['C', 'U', 'G']), 'D': set(['A', 'U', 'G']),
            'H': set(['A', 'C', 'U']), 'K': set(['U', 'G']),
            'M': set(['A', 'C']), 'N': set(['A', 'C', 'U', 'G']),
            'S': set(['C', 'G']), 'R': set(['A', 'G']), 'W': set(['A', 'U']),
            'V': set(['A', 'C', 'G']), 'Y': set(['C', 'U'])
        }
        self.assertEqual(self.b1.degenerate_map, exp)
        self.assertEqual(RNA.degenerate_map, exp)

    def test_degenerate_chars(self):
        exp = set(['B', 'D', 'H', 'K', 'M', 'N', 'S', 'R', 'W', 'V', 'Y'])
        self.assertEqual(self.b1.degenerate_chars, exp)
        self.assertEqual(RNA.degenerate_chars, exp)

    def test_complement(self):
        # use equals method to ensure that id, description, and quality are
        # correctly propagated to the resulting sequence
        self.assertTrue(self.b1.complement().equals(RNA("CUAAUGU")))

        self.assertTrue(self.b2.complement().equals(
            RNA("UGGCCAUGG", id="test-seq-2",
                        description="A test sequence", quality=range(9))))


        self.assertTrue(self.b4.complement().equals(
            RNA("KYWSRMBDHVN", id="degen",
                        description="All of the degenerate bases")))

        self.assertTrue(self.b5.complement().equals(
            RNA(".C--UAAUG-U...")))

    def test_reverse_complement(self):
        # use equals method to ensure that id, description, and (reversed)
        # quality scores are correctly propagated to the resulting sequence
        self.assertTrue(self.b1.reverse_complement().equals(
            RNA("UGUAAUC")))

        self.assertTrue(self.b2.reverse_complement().equals(
            RNA("GGUACCGGU", id="test-seq-2",
                        description="A test sequence",
                        quality=range(9)[::-1])))


        self.assertTrue(self.b4.reverse_complement().equals(
            RNA("NVHDBMRSWYK", id="degen",
                        description="All of the degenerate bases")))

    def test_is_reverse_complement(self):
        self.assertFalse(self.b1.is_reverse_complement(self.b1))

        # id, description, and quality scores should be ignored (only sequence
        # data and type should be compared)
        self.assertTrue(self.b1.is_reverse_complement(
            RNA('UGUAAUC', quality=range(7))))

        self.assertTrue(
            self.b4.is_reverse_complement(RNA('NVHDBMRSWYK')))

    def test_nondegenerates_no_degens(self):
        self.assertEqual(list(self.b1.nondegenerates()), [self.b1])

    def test_nondegenerates_all_degens(self):
        # Same chars.
        exp = [RNA('CC'), RNA('CG'), RNA('GC'),
               RNA('GG')]
        # Sort based on sequence string, as order is not guaranteed.
        obs = sorted(RNA('SS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Different chars.
        exp = [RNA('AC'), RNA('AG'), RNA('GC'),
               RNA('GG')]
        obs = sorted(RNA('RS').nondegenerates(), key=str)
        self.assertEqual(obs, exp)

        # Odd number of chars.
        obs = list(RNA('NNN').nondegenerates())
        self.assertEqual(len(obs), 4**3)

    def test_nondegenerates_mixed_degens(self):
        exp = [RNA('AGC'), RNA('AGU'), RNA('GGC'),
               RNA('GGU')]
        obs = sorted(RNA('RGY').nondegenerates(), key=str)
        self.assertEqual(obs, exp)


if __name__ == "__main__":
    main()
