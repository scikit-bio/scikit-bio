# ----------------------------------------------------------------------------
# Copyright (c) 2013--, scikit-bio development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE.txt, distributed with this software.
# ----------------------------------------------------------------------------

import unittest

from skbio import DNA, RNA
from skbio.metadata import IntervalMetadata


# tests specific to RNA go here. tests for functionality shared by DNA and RNA
# go in test_nucleotide_sequences.py
class TestRNA(unittest.TestCase):
    def test_reverse_transcribe(self):
        # without changes
        self.assertEqual(RNA('').reverse_transcribe(), DNA(''))
        self.assertEqual(RNA('A').reverse_transcribe(), DNA('A'))
        self.assertEqual(RNA('.ACGW-').reverse_transcribe(), DNA('.ACGW-'))

        # with changes
        self.assertEqual(DNA('T'), RNA('U').reverse_transcribe())
        self.assertEqual(DNA('TT'), RNA('UU').reverse_transcribe())
        self.assertEqual(DNA('ATCTG'), RNA('AUCUG').reverse_transcribe())
        self.assertEqual(DNA('TTTG'), RNA('UUUG').reverse_transcribe())

    def test_reverse_transcribe_preserves_all_metadata(self):
        im = IntervalMetadata(4)
        im.add([(0, 2)], metadata={'gene': 'p53'})

        seq = RNA('AGUU', metadata={'foo': 'bar'},
                  positional_metadata={'foo': range(4)},
                  interval_metadata=im)
        exp = DNA('AGTT', metadata={'foo': 'bar'},
                  positional_metadata={'foo': range(4)},
                  interval_metadata=im)
        self.assertEqual(seq.reverse_transcribe(), exp)

    def test_reverse_transcribe_does_not_modify_input(self):
        seq = RNA('AUAU')
        self.assertEqual(seq.reverse_transcribe(), DNA('ATAT'))
        self.assertEqual(seq, RNA('AUAU'))


if __name__ == '__main__':
    unittest.main()
