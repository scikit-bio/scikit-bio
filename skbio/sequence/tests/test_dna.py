# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

from __future__ import absolute_import, division, print_function

import unittest

from skbio import DNA, RNA


# tests specific to DNA go here. tests for functionality shared by DNA and RNA
# go in test_nucleotide_sequences.py
class TestDNA(unittest.TestCase):
    def test_transcribe(self):
        # without changes
        self.assertEqual(DNA('').transcribe(), RNA(''))
        self.assertEqual(DNA('A').transcribe(), RNA('A'))
        self.assertEqual(DNA('.ACGW-').transcribe(), RNA('.ACGW-'))

        # with changes
        self.assertEqual(DNA('T').transcribe(), RNA('U'))
        self.assertEqual(DNA('TT').transcribe(), RNA('UU'))
        self.assertEqual(DNA('ATCTG').transcribe(), RNA('AUCUG'))
        self.assertEqual(DNA('TTTG').transcribe(), RNA('UUUG'))

    def test_transcribe_preserves_all_metadata(self):
        exp = RNA('AGUU', metadata={'foo': 'bar'},
                  positional_metadata={'foo': range(4)})
        seq = DNA('AGTT', metadata={'foo': 'bar'},
                  positional_metadata={'foo': range(4)})
        self.assertEqual(seq.transcribe(), exp)

    def test_transcribe_does_not_modify_input(self):
        seq = DNA('ATAT')
        self.assertEqual(seq.transcribe(), RNA('AUAU'))
        self.assertEqual(seq, DNA('ATAT'))


if __name__ == '__main__':
    unittest.main()
